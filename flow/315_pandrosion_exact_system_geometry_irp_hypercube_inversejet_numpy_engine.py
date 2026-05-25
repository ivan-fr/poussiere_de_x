#!/usr/bin/env python3
"""
315_pandrosion_exact_system_geometry_irp_hypercube_inversejet_numpy_engine.py

But du fichier
--------------
315 est une couche d'orchestration compacte autour du moteur numerique 314.
Le point important est la separation stricte entre:

  1. le systeme exact a resoudre;
  2. une geometrie optionnelle qui sert seulement a choisir de meilleurs
     points de depart;
  3. le correcteur Pandrosion/IRP, qui valide toujours avec le systeme exact.

Pourquoi ce fichier existe
--------------------------
Dans 314, `--system-mode geometry-kernel` construit un systeme geometrique
Kostlan. C'est utile pour les tres grands `KS(n,d)`, mais ce n'est pas le meme
probleme qu'un polynome donne explicitement par l'utilisateur.

Dans 315, un polynome utilisateur reste le polynome exact:

    --system-source polynomial --polys "x^2 - 3*x - 10"

Si on active la geometrie:

    --geometry-mode kernel

ou, par compatibilite:

    --system-mode geometry-kernel

la geometrie est ajustee autour des valeurs exactes F(anchor). Elle peut guider
les starts, mais elle ne remplace jamais F. Les residus exportes, les racines
acceptees et les doublons sont calcules avec le systeme exact.

Ce qui est reutilise
--------------------
Le moteur lourd 304/306/312/314 reste dans:

    flow/314_pandrosion_geometry_kostlan_irp_hypercube_inversejet_numpy_engine.py

315 importe ce moteur localement par chemin de fichier et ne recopie pas les
milliers de lignes de correcteurs:

  - atlas universel 304;
  - hypercube inverse-jet 306;
  - IRP direct/lazy 312;
  - backends KS dense/lazy/geometrie de 314;
  - helpers NumPy, JSON, clustering et parametrage.

Ce fichier ajoute seulement:

  - un parseur d'expressions polynomiales volontairement petit et sans eval;
  - un systeme polynomial exact `ExpressionPolynomialSystem`;
  - une geometrie ajustee `FittedGeometryKernelSystem`;
  - un wrapper `GeometryWrappedSystem` qui expose F exact au correcteur;
  - un polissage local fini-difference, optionnel et valide sur F exact;
  - une CLI 315 qui combine ces pieces.

Exemple de verification simple
------------------------------
La commande suivante active la geometrie, mais valide sur le vrai polynome:

    python3 flow/315_pandrosion_exact_system_geometry_irp_hypercube_inversejet_numpy_engine.py \
      --system-source polynomial \
      --system-mode geometry-kernel \
      --cases 1,2 \
      --polys "x^2 - 3*x - 10" \
      --starts=-8,4 \
      --count 2 \
      --startopt-steps 0

Les racines attendues sont `-2` et `5`. Une partie imaginaire de l'ordre de
1e-18 est seulement du bruit flottant.

Dependances
-----------
Python standard library + NumPy, plus le fichier 314 local. Aucune dependance
SciPy, sympy, homotopy, Newton analytique ou solveur externe n'est ajoutee.
"""
from __future__ import annotations

import argparse
import ast
import dataclasses
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


def _load_core314() -> Any:
    path = Path(__file__).with_name(
        "314_pandrosion_geometry_kostlan_irp_hypercube_inversejet_numpy_engine.py"
    )
    spec = importlib.util.spec_from_file_location("_pandrosion_314_core", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load 314 core from {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


core = _load_core314()
core.ensure_numpy()
np = core.np


class SafeExpression:
    """Tiny AST evaluator for polynomial-style expressions."""

    _allowed = (
        ast.Expression,
        ast.BinOp,
        ast.UnaryOp,
        ast.Constant,
        ast.Name,
        ast.Add,
        ast.Sub,
        ast.Mult,
        ast.Div,
        ast.Pow,
        ast.USub,
        ast.UAdd,
        ast.Load,
    )

    def __init__(self, raw: str) -> None:
        self.raw = str(raw).strip()
        tree = ast.parse(self.raw.replace("^", "**"), mode="eval")
        self._validate(tree)
        self.tree = tree

    def _validate(self, node: ast.AST) -> None:
        if not isinstance(node, self._allowed):
            raise ValueError(f"unsupported expression node: {type(node).__name__}")
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            if not isinstance(node.right, ast.Constant):
                raise ValueError("polynomial powers must be numeric constants")
            if isinstance(node.right.value, complex):
                raise ValueError("polynomial powers must be real constants")
            if float(node.right.value) != int(node.right.value):
                raise ValueError("polynomial powers must be integer constants")
        for child in ast.iter_child_nodes(node):
            self._validate(child)

    def eval(self, env: dict[str, Any]) -> Any:
        return self._eval_node(self.tree.body, env)

    def _eval_node(self, node: ast.AST, env: dict[str, Any]) -> Any:
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id not in env:
                raise ValueError(f"unknown variable {node.id!r}")
            return env[node.id]
        if isinstance(node, ast.UnaryOp):
            val = self._eval_node(node.operand, env)
            if isinstance(node.op, ast.USub):
                return -val
            if isinstance(node.op, ast.UAdd):
                return val
        if isinstance(node, ast.BinOp):
            left = self._eval_node(node.left, env)
            right = self._eval_node(node.right, env)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                return left ** int(right)
        raise ValueError(f"unsupported expression node: {type(node).__name__}")


@dataclasses.dataclass
class ExpressionPolynomialSystem:
    n: int
    d: int
    seed: int
    expressions_raw: list[str]
    expressions: list[SafeExpression]
    variable_names: list[str]
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(
        cls,
        n: int,
        d: int,
        raw: str,
        variable_names: Optional[Sequence[str]] = None,
        seed_index: int = 0,
    ) -> "ExpressionPolynomialSystem":
        parts = [p.strip() for p in str(raw).replace("|", ";").split(";") if p.strip()]
        if not parts:
            raise ValueError("--polys/--poly is required for --system-source polynomial")
        names = [str(p).strip() for p in (variable_names or []) if str(p).strip()]
        if not names:
            names = ["x"] if int(n) == 1 else [f"x{i + 1}" for i in range(int(n))]
        if len(names) != int(n):
            raise ValueError(f"expected {n} variable names, got {len(names)}")
        if len(parts) != int(n):
            raise ValueError(f"expected {n} polynomial expression(s), got {len(parts)}")
        exprs = [SafeExpression(p) for p in parts]
        seed = core.stable_seed(n, d, seed_index, salt=0x315EAC7)
        return cls(int(n), int(d), seed, parts, exprs, names)

    @property
    def terms_per_poly(self) -> int:
        return core.exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d**self.n)

    @property
    def generation_seconds(self) -> float:
        return 0.0

    def _env(self, zz: Any) -> dict[str, Any]:
        env: dict[str, Any] = {"I": 1j, "j": 1j}
        for idx, name in enumerate(self.variable_names):
            env[name] = zz[:, idx]
        if self.n == 1:
            env.setdefault("x", zz[:, 0])
            env.setdefault("z", zz[:, 0])
        for idx in range(self.n):
            env.setdefault(f"x{idx + 1}", zz[:, idx])
            env.setdefault(f"z{idx + 1}", zz[:, idx])
        return env

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = time.time()
        out = self.eval_batch(np.asarray(z, dtype=np.complex128)[None, :])[0]
        self.seconds_eval += time.time() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = time.time()
        zz = np.asarray(Z, dtype=np.complex128)
        if zz.ndim == 1:
            zz = zz[None, :]
        env = self._env(zz)
        cols = []
        for expr in self.expressions:
            val = expr.eval(env)
            arr = np.asarray(val, dtype=np.complex128)
            if arr.ndim == 0:
                arr = np.full(int(zz.shape[0]), complex(arr), dtype=np.complex128)
            cols.append(arr.reshape(int(zz.shape[0])))
        out = np.stack(cols, axis=1)
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += int(zz.shape[0])
        self.seconds_eval += time.time() - t0
        return out

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        t0 = time.time()
        aa = np.asarray(a, dtype=np.complex128)
        bb = np.asarray(b, dtype=np.complex128)
        cur = aa.copy()
        f_prev = self.eval(cur)
        Q = np.zeros((self.n, self.n), dtype=np.complex128)
        for j in range(self.n):
            old = cur[j]
            cur[j] = bb[j]
            f_next = self.eval(cur)
            dz = bb[j] - old
            if abs(dz) > 1e-300:
                Q[:, j] = (f_next - f_prev) / dz
            else:
                h = 1e-6 * max(1.0, abs(old))
                plus = cur.copy()
                minus = cur.copy()
                plus[j] = old + h
                minus[j] = old - h
                Q[:, j] = (self.eval(plus) - self.eval(minus)) / (2.0 * h)
            f_prev = f_next
        self.slope_count += 1
        self.seconds_slope += time.time() - t0
        return Q

    def stats(self) -> dict[str, Any]:
        return {
            "eval_count": int(self.eval_count),
            "slope_count": int(self.slope_count),
            "seconds_eval": float(self.seconds_eval),
            "seconds_slope": float(self.seconds_slope),
            "terms_per_poly": self.terms_per_poly,
            "total_terms": self.total_terms,
            "expressions": list(self.expressions_raw),
            "variables": list(self.variable_names),
        }


@dataclasses.dataclass
class FittedGeometryKernelSystem:
    base: Any
    kernel: Any
    weights: Any
    ridge: float
    fit_seconds: float
    eval_count: int = 0
    seconds_eval: float = 0.0

    @classmethod
    def make(cls, base: Any, n: int, d: int, args: argparse.Namespace) -> "FittedGeometryKernelSystem":
        t0 = time.time()
        kernel = core.GeometryKernelKostlanSystem.make(
            n,
            d,
            seed_index=int(getattr(args, "seed_index", 0)),
            equation_normalize=False,
            geometry_anchors=int(getattr(args, "geometry_anchors", 0)),
            geometry_anchor_cap=int(getattr(args, "geometry_anchor_cap", 4096)),
            geometry_anchor_scales=core.parse_float_list(
                getattr(args, "geometry_anchor_scales", None),
                [0.25, 0.5, 1.0, 2.0, 4.0],
                positive=True,
            ),
            dynamic_normalize=bool(getattr(args, "geometry_dynamic_normalize", True)),
            self_normalize=bool(getattr(args, "geometry_self_normalize", True)),
            eval_block=int(getattr(args, "geometry_eval_block", 128)),
        )
        A = kernel._kernel_block(kernel.anchors)
        F = base.eval_batch(kernel.anchors)
        ridge = max(0.0, float(getattr(args, "geometry_fit_ridge", 1e-8)))
        eye = np.eye(int(A.shape[0]), dtype=np.complex128)
        try:
            W = np.linalg.solve(A + ridge * eye, F)
            method = "solve"
        except Exception:
            W, _, _, _ = np.linalg.lstsq(A + ridge * eye, F, rcond=None)
            method = "lstsq"
        obj = cls(base, kernel, np.asarray(W, dtype=np.complex128), ridge, time.time() - t0)
        obj._fit_method = method
        return obj

    @property
    def n(self) -> int:
        return int(self.base.n)

    @property
    def d(self) -> int:
        return int(self.base.d)

    @property
    def seed(self) -> int:
        return int(self.kernel.seed)

    @property
    def terms_per_poly(self) -> int:
        return int(self.base.terms_per_poly)

    @property
    def total_terms(self) -> int:
        return int(self.base.total_terms)

    @property
    def bezout(self) -> int:
        return int(self.base.bezout)

    @property
    def generation_seconds(self) -> float:
        return float(self.fit_seconds)

    def eval(self, z: Sequence[complex]) -> Any:
        return self.eval_batch(np.asarray(z, dtype=np.complex128)[None, :])[0]

    def eval_batch(self, Z: Any) -> Any:
        t0 = time.time()
        K = self.kernel._kernel_block(np.asarray(Z, dtype=np.complex128))
        out = K @ self.weights
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += int(out.shape[0])
        self.seconds_eval += time.time() - t0
        return out

    def residual(self, z: Sequence[complex]) -> float:
        return float(np.linalg.norm(self.eval(z)))

    def residuals_batch(self, Z: Any) -> Any:
        F = self.eval_batch(Z)
        R = np.linalg.norm(F, axis=1)
        R = np.asarray(R, dtype=float)
        R[~np.isfinite(R)] = np.inf
        return R

    def stats(self) -> dict[str, Any]:
        return {
            "geometry_fit_eval_count": int(self.eval_count),
            "geometry_fit_seconds_eval": float(self.seconds_eval),
            "geometry_fit_seconds": float(self.fit_seconds),
            "geometry_fit_ridge": float(self.ridge),
            "geometry_fit_method": str(getattr(self, "_fit_method", "unknown")),
            "geometry_anchors": int(self.kernel.geometry_anchors),
        }


@dataclasses.dataclass
class GeometryWrappedSystem:
    exact: Any
    geometry_system: FittedGeometryKernelSystem
    geometry_mode: str = "kernel"

    @property
    def n(self) -> int:
        return int(self.exact.n)

    @property
    def d(self) -> int:
        return int(self.exact.d)

    @property
    def seed(self) -> int:
        return int(self.exact.seed)

    @property
    def terms_per_poly(self) -> int:
        return int(self.exact.terms_per_poly)

    @property
    def total_terms(self) -> int:
        return int(self.exact.total_terms)

    @property
    def bezout(self) -> int:
        return int(self.exact.bezout)

    @property
    def generation_seconds(self) -> float:
        return float(self.exact.generation_seconds) + float(self.geometry_system.fit_seconds)

    def eval(self, z: Sequence[complex]) -> Any:
        return self.exact.eval(z)

    def eval_batch(self, Z: Any) -> Any:
        return self.exact.eval_batch(Z)

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return self.exact.slope_matrix(a, b)

    def stats(self) -> dict[str, Any]:
        out = dict(self.exact.stats())
        out.update({"geometry_mode": self.geometry_mode, **self.geometry_system.stats()})
        return out


def parse_complex_token(raw: str) -> complex:
    s = str(raw).strip().replace("i", "j")
    if not s:
        raise ValueError("empty complex token")
    return complex(s)


def parse_start_points(raw: Optional[str], n: int) -> list[Any]:
    if raw is None or str(raw).strip() == "":
        return []
    text = str(raw).strip()
    if int(n) == 1 and ";" not in text and "|" not in text:
        return [np.asarray([parse_complex_token(p)], dtype=np.complex128) for p in text.split(",") if p.strip()]
    points = []
    for part in text.replace("|", ";").split(";"):
        if not part.strip():
            continue
        coords = [p.strip() for p in part.split(",") if p.strip()]
        if len(coords) != int(n):
            raise ValueError(f"start point {part!r} has {len(coords)} coord(s), expected {n}")
        points.append(np.asarray([parse_complex_token(p) for p in coords], dtype=np.complex128))
    return points


def make_system_315(args: argparse.Namespace, n: int, d: int) -> tuple[Any, str, str]:
    source = str(getattr(args, "system_source", "ks")).strip().lower().replace("_", "-")
    geometry_mode = str(getattr(args, "geometry_mode", "none")).strip().lower().replace("_", "-")
    if source in {"poly", "polynomial", "expr", "expression"}:
        var_raw = str(getattr(args, "variables", "") or "").strip()
        variables = [p.strip() for p in var_raw.replace(";", ",").split(",") if p.strip()] if var_raw else None
        exact = ExpressionPolynomialSystem.make(n, d, str(getattr(args, "polys", "") or ""), variables, int(args.seed_index))
        mode = str(getattr(args, "system_mode", "auto")).strip().lower().replace("_", "-")
        if geometry_mode in {"", "none"} and mode in {"geometry", "geometry-kernel", "kernel", "projective-kernel"}:
            geometry_mode = "kernel"
        if geometry_mode in {"kernel", "geometry", "geometry-kernel", "projective-kernel"}:
            guide = FittedGeometryKernelSystem.make(exact, n, d, args)
            return GeometryWrappedSystem(exact, guide, "kernel"), "polynomial", "exact+geometry-kernel"
        return exact, "polynomial", "exact"
    if source not in {"ks", "kostlan"}:
        raise ValueError(f"unknown --system-source {source!r}; use ks or polynomial")
    system = core.make_kostlan_system(args, n, d)
    if isinstance(system, core.GeometryKernelKostlanSystem):
        backend = "geometry-kernel"
    elif isinstance(system, core.LazyFeatureKostlanSystem):
        backend = "lazy-feature"
    else:
        backend = "dense"
    return system, "ks", backend


def polish_exact_root(
    system: Any,
    chart: Any,
    y: Any,
    accept: float,
    max_steps: int,
) -> tuple[Any, Any, float, dict[str, Any]]:
    z = chart.z_from_y(y)
    residual = float(np.linalg.norm(system.eval(z)))
    meta = {
        "exact_polish_steps": 0,
        "exact_polish_r0": residual,
        "exact_polish_r1": residual,
        "exact_polish_enabled": bool(max_steps > 0),
    }
    if max_steps <= 0 or not math.isfinite(residual):
        return z, np.asarray(y, dtype=np.complex128), residual, meta

    best_z = np.asarray(z, dtype=np.complex128).copy()
    best_f = system.eval(best_z)
    best_r = float(np.linalg.norm(best_f))
    gains = (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
    for step in range(int(max_steps)):
        if best_r <= float(accept):
            break
        try:
            J = system.slope_matrix(best_z, best_z)
            delta, _, _, _ = np.linalg.lstsq(np.asarray(J, dtype=np.complex128), -best_f, rcond=None)
        except Exception:
            meta["exact_polish_status"] = "linear-solve-failed"
            break
        delta = np.asarray(delta, dtype=np.complex128)
        if not np.all(np.isfinite(delta)):
            meta["exact_polish_status"] = "non-finite-step"
            break
        improved = False
        for gain in gains:
            cand_z = best_z + complex(gain) * delta
            cand_f = system.eval(cand_z)
            cand_r = float(np.linalg.norm(cand_f))
            if math.isfinite(cand_r) and cand_r < best_r:
                best_z, best_f, best_r = cand_z, cand_f, cand_r
                meta["exact_polish_steps"] = step + 1
                improved = True
                break
        if not improved:
            meta["exact_polish_status"] = "no-decrease"
            break
    meta["exact_polish_r1"] = best_r
    return best_z, chart.y_from_z(best_z), best_r, meta


def _source_kind(args: argparse.Namespace) -> str:
    return str(getattr(args, "system_source", "ks")).strip().lower().replace("_", "-")


def _patch_ks_result(result: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    result = dict(result)
    result["script"] = Path(__file__).name
    result["mode"] = "315-ks-via-314-core"
    result["system_source"] = "ks"
    params = dict(result.get("parameters", {}))
    params.update({
        "system_source": str(args.system_source),
        "geometry_mode": str(args.geometry_mode),
        "geometry_starts": bool(args.geometry_starts),
        "starts": str(args.starts or ""),
    })
    result["parameters"] = params
    return result


def _base_meta(args: argparse.Namespace, system: Any, source: str, backend: str, n: int, d: int, chart: Any) -> dict[str, Any]:
    geom_kernel = getattr(getattr(system, "geometry_system", None), "kernel", system)
    return {
        "script": Path(__file__).name,
        "autonomous": True,
        "dependencies": {"python_scripts": ["314"], "numpy": bool(np is not None)},
        "mode": "315-exact-system-plus-geometry-guided-312-lazy-irp-hypercube-inversejet",
        "flow_formula": "exact F -> optional fitted geometry guide for starts -> 312 lazy IRP / 306 hypercube inverse-jet on exact F -> exact residual validation",
        "case": f"{n},{d}",
        "family": source,
        "system_source": source,
        "system_backend": backend,
        "seed_index": int(args.seed_index),
        "seed": int(system.seed),
        "n": int(n),
        "degree": int(d),
        "terms_per_poly": int(system.terms_per_poly),
        "terms": int(system.total_terms),
        "bezout": int(system.bezout),
        "equation_normalize": bool(args.equation_normalize),
        "linear_A": [[core.cjson(chart.A[i, j]) for j in range(n)] for i in range(n)],
        "parameters": {
            "system_source": str(args.system_source),
            "polys": str(args.polys or ""),
            "variables": str(args.variables or ""),
            "geometry_mode": str(args.geometry_mode),
            "geometry_starts": bool(args.geometry_starts),
            "geometry_fit_ridge": float(args.geometry_fit_ridge),
            "exact_polish_steps": int(args.exact_polish_steps),
            "starts": str(args.starts or ""),
            "system_mode": str(args.system_mode),
            "count": int(args.count),
            "pool": int(args.pool),
            "accept": float(args.accept),
            "tol": float(args.tol),
            "epochs": int(args.epochs),
            "cluster_sep": float(args.cluster_sep),
            "line_search": int(args.line_search),
            "hypercube_nodes": int(args.hypercube_nodes),
            "startopt_steps": int(args.startopt_steps),
            "startopt_candidates": int(args.startopt_candidates),
            "geometry_anchors": int(getattr(geom_kernel, "geometry_anchors", 0) or 0),
            "geometry_anchor_cap": int(args.geometry_anchor_cap),
            "geometry_anchor_scales": core.parse_float_list(args.geometry_anchor_scales, [0.25, 0.5, 1.0, 2.0, 4.0], positive=True),
        },
    }


def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    if _source_kind(args) in {"ks", "kostlan"}:
        return _patch_ks_result(core.run_case(args, case_raw), args)

    t_case = time.time()
    n, d = core.parse_case(case_raw)
    system, source, backend = make_system_315(args, n, d)
    chart = core.LinearChart.identity(n, scale=float(args.linear_scale))
    target = core.TargetTrack(system, chart)
    guide_system = getattr(system, "geometry_system", None)
    guide_target = core.TargetTrack(guide_system, chart) if bool(args.geometry_starts) and guide_system is not None else target
    starts = parse_start_points(getattr(args, "starts", None), n)
    powers = sorted(set(round(float(x), 16) for x in core.parse_float_list(args.powers, core.DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in core.parse_float_list(args.angles, core.DEFAULT_ANGLES_DEG)]
    radii = core.parse_float_list(args.rays, core.DEFAULT_RADII, positive=True)
    gains = core.parse_float_list(args.startopt_gains, core.DEFAULT_GAINS, positive=True)
    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    failures = 0
    duplicates = 0
    t_extract = time.time()

    for trial in range(int(args.pool)):
        if len(roots) >= int(args.count):
            break
        explicit = trial < len(starts)
        if explicit:
            y_raw = np.asarray(starts[trial], dtype=np.complex128).copy()
            geom = {
                "chart": "315-explicit-start",
                "atlas_mode": "315-explicit-start",
                "atlas_selected_geometry_cell": "explicit-start",
                "homothety": 1.0,
                "theta_deg": None,
            }
        else:
            y_raw, geom = core.universal_atlas_start(
                guide_target,
                n,
                trial,
                system.seed + 0x113000,
                powers,
                angles,
                radii,
                float(args.power_cap),
                roots_found=len(roots),
                duplicates=duplicates,
                failures=failures,
                target_count=int(args.count),
                universal_cells=int(args.universal_cells),
                universal_shells=int(args.universal_shells),
                cell_probe_radius=float(args.cell_probe_radius),
                cell_descent_min=float(args.cell_descent_min),
                cell_equal_gap_min=float(args.cell_equal_gap_min),
                cell_log_max=float(args.cell_log_max),
                universal_cycle=bool(args.universal_cycle),
            )
        start_target = target if explicit else guide_target
        y0, smeta = core.startopt(
            start_target,
            y_raw,
            trial,
            system.seed + 0x112555,
            int(args.startopt_steps),
            int(args.startopt_candidates),
            gains,
            int(args.startopt_micro_epochs),
        )
        loc = core.lazy_irp_hypercube_inversejet_corrector_312(
            target,
            y0,
            max_epochs=int(args.epochs),
            tol=float(args.tol),
            accept=float(args.accept),
            trial_timeout=float(args.trial_timeout),
            line_search=int(args.line_search),
            line_grid=core.parse_float_list(args.line_grid, []),
            direction_seed=system.seed + 7919 * trial,
            cloud_nodes=int(args.hypercube_nodes),
            irp_layers=int(args.irp_layers),
            irp_inner_epochs=int(args.irp_inner_epochs),
            irp_scales=core._311_complex_scale_palette(
                core.parse_float_list(args.irp_gains, [1.0, 0.5, 2.0, 0.25, 4.0, 0.125, 8.0], positive=True),
                core.parse_float_list(args.irp_phases, [0.0, 0.08, -0.08, 0.19, -0.19]),
                int(args.irp_top),
            ),
            irp_chart_top=int(args.irp_chart_top),
            irp_inversion=bool(args.irp_inversion),
            collapse=bool(args.collapse),
            collapse_residual=float(args.collapse_residual),
            collapse_drop=float(args.collapse_drop),
            collapse_rel_step=float(args.collapse_rel_step),
            collapse_after=int(args.collapse_after),
            local_inner_epochs=int(args.local_inner_epochs),
            lazy_direct_epochs=int(args.lazy_direct_epochs),
            lazy_trigger_drop=float(args.lazy_trigger_drop),
            lazy_trigger_after=int(args.lazy_trigger_after),
            lazy_bad_cond=float(args.lazy_bad_cond),
            lazy_log_energy=float(args.lazy_log_energy),
            eager_irp=bool(args.eager_irp),
            rescue_collapsed=bool(args.rescue_collapsed),
        )
        z, y_polished, residual, polish_meta = polish_exact_root(
            system,
            chart,
            loc["y"],
            float(args.accept),
            int(args.exact_polish_steps),
        )
        accepted = bool(math.isfinite(residual) and residual < float(args.accept))
        rec = {
            "trial": int(trial),
            "accepted": accepted,
            "status": loc.get("status"),
            "r1": residual,
            "epochs": int(loc.get("epochs", 0)),
            "seconds": float(loc.get("seconds", 0.0)),
            **geom,
            **smeta,
            **polish_meta,
        }
        if bool(args.verbose_trials):
            rec["z"] = core.root_to_json(z)
            rec["y0"] = core.root_to_json(y0)
        if not accepted:
            failures += 1
            trials.append(rec)
            continue
        dup = core.cluster_index(roots, z, float(args.cluster_sep))
        if dup is not None:
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = int(dup)
            trials.append(rec)
            continue
        root = {
            "id": len(roots),
            "source": "315-exact-system-geometry-guided",
            "trial": int(trial),
            "z_complex": np.asarray(z, dtype=np.complex128).copy(),
            "y_complex": np.asarray(y_polished, dtype=np.complex128).copy(),
            "residual": residual,
            "realness": core.realness(z),
            "cond": core.slope_condition_from_corrector(loc),
            "epochs": int(loc.get("epochs", 0)),
            "seconds": float(loc.get("seconds", 0.0)),
            **geom,
            **smeta,
            **polish_meta,
        }
        root["score"] = core.score_root(root["residual"], root["realness"], root["cond"])
        roots.append(root)
        rec["status"] = "new-root"
        rec["root_id"] = root["id"]
        trials.append(rec)

    encoded_roots = []
    for root in sorted(roots, key=lambda q: (float(q.get("score", float("inf"))), int(q.get("id", 0)))):
        rr = dict(root)
        rr["z"] = core.root_to_json(rr.pop("z_complex"))
        rr["y"] = core.root_to_json(rr.pop("y_complex"))
        encoded_roots.append(rr)
    result = _base_meta(args, system, source, backend, n, d, chart)
    result.update({
        "roots": encoded_roots,
        "trials": trials if bool(args.verbose_trials) else trials[: min(len(trials), int(args.keep_trials))],
        "summary": {
            "requested_roots": int(args.count),
            "unique_roots": len(roots),
            "success": bool(len(roots) >= int(args.count)),
            "trials_used": len(trials),
            "duplicates": int(duplicates),
            "failures": int(failures),
            "generation_seconds": float(system.generation_seconds),
            "extract_seconds": float(time.time() - t_extract),
            "total_seconds": float(time.time() - t_case),
            "eval_stats": system.stats(),
        },
    })
    return result


def build_parser() -> argparse.ArgumentParser:
    p = core.build_parser()
    p.description = (
        "315 exact-system Pandrosion engine: exact KS or user polynomials, "
        "optional fitted geometry-kernel starts, exact residual validation."
    )
    p.add_argument("--system-source", choices=["ks", "kostlan", "polynomial", "poly", "expr", "expression"], default="ks")
    p.add_argument("--polys", "--poly", default=None)
    p.add_argument("--variables", default=None)
    p.add_argument("--geometry-mode", choices=["none", "kernel", "geometry", "geometry-kernel", "projective-kernel"], default="none")
    p.add_argument("--geometry-starts", action="store_true", default=True)
    p.add_argument("--no-geometry-starts", dest="geometry_starts", action="store_false")
    p.add_argument("--geometry-fit-ridge", type=float, default=1e-8)
    p.add_argument("--exact-polish-steps", type=int, default=6)
    p.add_argument("--starts", default=None)
    p.set_defaults(outdir="/mnt/data/315_exact_system_geometry_out")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(getattr(args, "self_test", False)):
        args.cases = "2,2"
        args.count = 4
        args.pool = min(int(args.pool), 512)
        args.epochs = min(int(args.epochs), 16)
        args.accept = min(float(args.accept), 1e-8)
        args.keep_trials = min(int(args.keep_trials), 20)
        args.out = args.out or "/mnt/data/315_exact_system_geometry_out/self_test_315.json"
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    outputs = [run_case(args, c) for c in cases]
    final = outputs[0] if len(outputs) == 1 else {"script": Path(__file__).name, "autonomous": True, "cases": outputs}
    out = Path(args.out) if args.out else Path(args.outdir) / f"315_exact_geometry_{cases[0].replace(',', 'x') if cases else 'case'}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2), encoding="utf-8")

    print("=" * 120, flush=True)
    print("315 EXACT SYSTEM + GEOMETRY-GUIDED LAZY IRP / HYPERCUBE INVERSE-JET", flush=True)
    print("Exact eval/validation; fitted geometry is only a guide for starts.", flush=True)
    print("=" * 120, flush=True)
    for r in outputs:
        s = r["summary"]
        backend = r.get("system_backend", "dense")
        family = r.get("family", "ks")
        geom = r.get("parameters", {}).get("geometry_anchors", 0)
        suffix = f", geometry_anchors={geom}" if geom else ""
        print(f"case={family}({r['n']},{r['degree']}), backend={backend}{suffix}, seed={r['seed']}", flush=True)
        print(f"roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} trials={s['trials_used']} duplicates={s['duplicates']} failures={s['failures']}", flush=True)
        print(f"seconds: generation={s['generation_seconds']:.2f}, extract={s['extract_seconds']:.2f}, total={s['total_seconds']:.2f}", flush=True)
        if r.get("roots"):
            best = r["roots"][0]
            print(f"best_root: residual={float(best.get('residual', float('inf'))):.3e}, trial={best.get('trial')}", flush=True)
    print(f"out={out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
