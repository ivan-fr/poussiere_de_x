#!/usr/bin/env python3
"""
116_pandrosion_multifamily_vectorized_pure_pandrosion.py

Thin multifamily test harness for the 115 vectorized PURE Pandrosion engine.

This file deliberately reuses the 115 flow/corrector:
  - Mobius/Riemann + Thales homothety starts
  - StartOpt
  - vectorized exact telescopic slope Q(a,b)
  - PURE Pandrosion local corrector
  - clustering and JSON encoding conventions

The only new layer is the polynomial-system generator.  It can run the same
case (n,d) across several multivariate families: dense Kostlan, sparse Kostlan,
IID dense/sparse, fewnomial, degree-shell, mixed-degree, real-coefficient,
phase-only, and structured diagonal/chain/cyclic systems.
"""
from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any, Optional, Sequence


def _load_engine_115() -> Any:
    path = Path(__file__).resolve().with_name("115_pandrosion_vectorized_pure_pandrosion.py")
    spec = importlib.util.spec_from_file_location("pandrosion_115_vectorized_pure", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load 115 engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ENGINE = _load_engine_115()
np = ENGINE.np


FAMILY_DESCRIPTIONS: dict[str, str] = {
    "ks": "Dense total-degree Kostlan-Shub-Smale coefficients, identical to 115 for the same seed.",
    "ks_sparse": "Common sparse support sampled from total-degree Kostlan terms, Kostlan weights kept.",
    "dense_iid": "Dense total-degree complex IID Gaussian coefficients without Kostlan weights.",
    "sparse_iid": "Common sparse total-degree support with complex IID Gaussian coefficients.",
    "fewnomial": "Very sparse common support: constants, pure powers, and a few random monomials.",
    "degree_shell_ks": "Kostlan-weighted constants/linear terms plus the total-degree d shell.",
    "mixed_degree": "Equation i has its own degree <= d, with Kostlan weights for that row degree.",
    "real_ks": "Dense total-degree Kostlan support with real Gaussian coefficients.",
    "phase_ks": "Dense total-degree Kostlan support with random unit-modulus phases.",
    "diagonal": "Structured system x_i^d = c_i with known Bezout shape.",
    "chain": "Sparse structured chain x_i^d + eps*x_{i+1} = c_i.",
    "cyclic": "Sparse structured cyclic chain x_i^d + eps*x_{i+1 mod n} = c_i.",
    "ill_scaled": "Dense Kostlan support with deterministic term and row scale spread.",
}

DEFAULT_FAMILIES = [
    "ks",
    "ks_sparse",
    "dense_iid",
    "sparse_iid",
    "fewnomial",
    "degree_shell_ks",
    "mixed_degree",
    "real_ks",
    "phase_ks",
    "diagonal",
    "chain",
    "cyclic",
    "ill_scaled",
]

FAMILY_GROUPS: dict[str, list[str]] = {
    "all": DEFAULT_FAMILIES,
    "random": [
        "ks",
        "ks_sparse",
        "dense_iid",
        "sparse_iid",
        "fewnomial",
        "degree_shell_ks",
        "mixed_degree",
        "real_ks",
        "phase_ks",
        "ill_scaled",
    ],
    "sparse": ["ks_sparse", "sparse_iid", "fewnomial", "diagonal", "chain", "cyclic"],
    "structured": ["diagonal", "chain", "cyclic"],
    "dense": ["ks", "dense_iid", "real_ks", "phase_ks", "ill_scaled"],
    "fast": ["ks", "ks_sparse", "fewnomial", "diagonal", "chain"],
}

FAMILY_ALIASES = {
    "kostlan": "ks",
    "kss": "ks",
    "ks_dense": "ks",
    "kostlan_sparse": "ks_sparse",
    "sparse_ks": "ks_sparse",
    "iid": "dense_iid",
    "iid_dense": "dense_iid",
    "dense_gaussian": "dense_iid",
    "sparse_gaussian": "sparse_iid",
    "few": "fewnomial",
    "fewnomials": "fewnomial",
    "shell": "degree_shell_ks",
    "degree_shell": "degree_shell_ks",
    "near_homogeneous_ks": "degree_shell_ks",
    "mixed": "mixed_degree",
    "real": "real_ks",
    "phase": "phase_ks",
    "binomial": "diagonal",
    "diag": "diagonal",
    "ill": "ill_scaled",
    "illscaled": "ill_scaled",
}


def ensure_numpy() -> None:
    ENGINE.ensure_numpy()


def _canonical_family(raw: str) -> str:
    key = str(raw).strip().lower().replace("-", "_")
    key = FAMILY_ALIASES.get(key, key)
    if key not in FAMILY_DESCRIPTIONS:
        valid = ", ".join(sorted(FAMILY_DESCRIPTIONS))
        groups = ", ".join(sorted(FAMILY_GROUPS))
        raise ValueError(f"unknown family {raw!r}; valid families: {valid}; groups: {groups}")
    return key


def parse_families(raw: Optional[str]) -> list[str]:
    if raw is None or str(raw).strip() == "":
        return list(DEFAULT_FAMILIES)
    out: list[str] = []
    text = str(raw).replace("|", ",").replace(";", ",")
    for part in text.split(","):
        token = part.strip().lower().replace("-", "_")
        if not token:
            continue
        group = FAMILY_GROUPS.get(token)
        names = group if group is not None else [_canonical_family(token)]
        for name in names:
            if name not in out:
                out.append(name)
    return out or list(DEFAULT_FAMILIES)


def _family_salt(name: str) -> int:
    if name == "ks":
        return 0
    h = 0x116D00D5EED
    for b in name.encode("utf-8"):
        h = ENGINE.splitmix64(h ^ b)
    return int(h & 0x7FFFFFFF)


def _family_seed(n: int, d: int, seed_index: int, family: str) -> int:
    return int(ENGINE.stable_seed(n, d, seed_index, salt=_family_salt(family)))


def _complex_gaussian(rng: Any, shape: tuple[int, ...]) -> Any:
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / math.sqrt(2.0)


def _row_normalize(coeff: Any) -> Any:
    norms = np.linalg.norm(coeff, axis=1)
    norms = np.where(norms > 0, norms, 1.0)
    return coeff / norms[:, None]


def _active_counts(coeff: Any) -> list[int]:
    return [int(x) for x in np.count_nonzero(np.abs(coeff) > 0, axis=1).tolist()]


def _exponent_index(exps: Any) -> dict[tuple[int, ...], int]:
    return {tuple(int(v) for v in row): i for i, row in enumerate(exps.tolist())}


def _pure_exp(n: int, j: int, degree: int) -> tuple[int, ...]:
    row = [0] * n
    row[j] = int(degree)
    return tuple(row)


def _linear_exp(n: int, j: int) -> tuple[int, ...]:
    return _pure_exp(n, j, 1)


def _unique_exps(rows: Sequence[Sequence[int]]) -> tuple[Any, dict[tuple[int, ...], int]]:
    seen: dict[tuple[int, ...], int] = {}
    out: list[tuple[int, ...]] = []
    for row in rows:
        key = tuple(int(x) for x in row)
        if key not in seen:
            seen[key] = len(out)
            out.append(key)
    max_degree = max((max(r) for r in out), default=0)
    dtype = np.int16 if max_degree < 32767 else np.int32
    return np.asarray(out, dtype=dtype), seen


def _required_total_degree_indices(exps: Any, n: int, d: int, include_linear: bool = False) -> list[int]:
    index = _exponent_index(exps)
    required: list[int] = []
    zero = tuple([0] * n)
    if zero in index:
        required.append(index[zero])
    for j in range(n):
        pure = _pure_exp(n, j, d)
        if pure in index:
            required.append(index[pure])
    if include_linear:
        for j in range(n):
            lin = _linear_exp(n, j)
            if lin in index:
                required.append(index[lin])
    return sorted(set(required))


def _auto_sparse_terms(family: str, n: int, full_terms: int, required_count: int, sparse_terms: int, sparse_frac: float) -> int:
    if sparse_terms and sparse_terms > 0:
        return max(required_count, min(full_terms, int(sparse_terms)))
    frac_terms = int(math.ceil(max(0.0, float(sparse_frac)) * full_terms))
    if family == "fewnomial":
        wanted = max(required_count, 2 * n + 5)
    else:
        wanted = max(required_count + 4, 4 * n + 6, frac_terms)
    return max(1, min(full_terms, wanted))


def _choose_support(exps: Any, rng: Any, count: int, required: Sequence[int]) -> Any:
    full_terms = int(exps.shape[0])
    required_sorted = sorted(set(int(i) for i in required if 0 <= int(i) < full_terms))
    count = max(len(required_sorted), min(full_terms, int(count)))
    pool = np.asarray([i for i in range(full_terms) if i not in set(required_sorted)], dtype=np.int64)
    if count > len(required_sorted) and len(pool) > 0:
        extra = rng.choice(pool, size=min(count - len(required_sorted), len(pool)), replace=False)
        idx = required_sorted + [int(i) for i in extra.tolist()]
    else:
        idx = required_sorted
    return np.asarray(sorted(set(idx)), dtype=np.int64)


def _structural_constants(rng: Any, n: int) -> list[complex]:
    constants: list[complex] = []
    for _ in range(n):
        radius = 0.72 + 0.72 * float(rng.random())
        theta = 2.0 * math.pi * float(rng.random())
        constants.append(radius * ENGINE.phase(theta))
    return constants


def _make_structured_family(family: str, n: int, d: int, rng: Any) -> tuple[Any, Any, Any, dict[str, Any]]:
    rows: list[tuple[int, ...]] = [tuple([0] * n)]
    for j in range(n):
        rows.append(_pure_exp(n, j, d))
    if family in {"chain", "cyclic"}:
        for j in range(n):
            rows.append(_linear_exp(n, j))
    exps, index = _unique_exps(rows)
    coeff = np.zeros((n, int(exps.shape[0])), dtype=np.complex128)
    constants = _structural_constants(rng, n)
    coupling = 0.0 if family == "diagonal" else 0.18
    zero_idx = index[tuple([0] * n)]
    for i in range(n):
        coeff[i, zero_idx] = -constants[i]
        coeff[i, index[_pure_exp(n, i, d)]] += 1.0 + 0.0j
        if family == "chain" and i + 1 < n:
            coeff[i, index[_linear_exp(n, i + 1)]] += coupling * ENGINE.phase(0.41 + i)
        elif family == "cyclic":
            coeff[i, index[_linear_exp(n, (i + 1) % n)]] += coupling * ENGINE.phase(0.41 + i)
    weights = np.ones(int(exps.shape[0]), dtype=np.float64)
    metadata = {
        "structured": True,
        "constants": [ENGINE.cjson(c) for c in constants],
        "coupling": float(coupling),
        "expected_bezout": int(d ** n),
    }
    return exps, weights, coeff, metadata


@dataclasses.dataclass
class MultiFamilySystem(ENGINE.DenseKostlanSystem):
    family: str = "ks"
    family_description: str = ""
    active_terms: int = 0
    active_terms_per_equation: list[int] = dataclasses.field(default_factory=list)
    generator_metadata: dict[str, Any] = dataclasses.field(default_factory=dict)

    @classmethod
    def make(
        cls,
        n: int,
        d: int,
        seed_index: int = 0,
        equation_normalize: bool = True,
        family: str = "ks",
        sparse_terms: int = 0,
        sparse_frac: float = 0.18,
    ) -> "MultiFamilySystem":
        ensure_numpy()
        family = _canonical_family(family)
        t0 = ENGINE.now()
        seed = _family_seed(n, d, seed_index, family)
        rng = np.random.default_rng(seed)
        metadata: dict[str, Any] = {"family": family}

        if family in {"ks", "ks_sparse", "dense_iid", "sparse_iid", "fewnomial", "real_ks", "phase_ks", "ill_scaled"}:
            full_exps = ENGINE.compositions_leq(d, n)
            full_weights = ENGINE.multinomial_kostlan_weights(full_exps, d)
            if family in {"ks_sparse", "sparse_iid", "fewnomial"}:
                required = _required_total_degree_indices(full_exps, n, d, include_linear=(family == "fewnomial"))
                count = _auto_sparse_terms(family, n, int(full_exps.shape[0]), len(required), sparse_terms, sparse_frac)
                support = _choose_support(full_exps, rng, count, required)
                exps = full_exps[support]
                weights = full_weights[support] if family == "ks_sparse" else np.ones(len(support), dtype=np.float64)
                metadata.update({
                    "support": "common-random",
                    "support_size": int(len(support)),
                    "full_terms_per_poly": int(full_exps.shape[0]),
                    "sparse_frac": float(sparse_frac),
                })
            else:
                exps = full_exps
                weights = full_weights if family in {"ks", "real_ks", "phase_ks", "ill_scaled"} else np.ones(int(full_exps.shape[0]), dtype=np.float64)
                metadata.update({"support": "dense-total-degree"})

            if family == "real_ks":
                coeff = rng.standard_normal((n, int(exps.shape[0]))) * weights[None, :]
            elif family == "phase_ks":
                angles = rng.uniform(0.0, 2.0 * math.pi, size=(n, int(exps.shape[0])))
                coeff = np.exp(1j * angles) * weights[None, :]
            else:
                coeff = _complex_gaussian(rng, (n, int(exps.shape[0]))) * weights[None, :]

            if family == "ill_scaled":
                totals = np.sum(exps, axis=1).astype(np.float64)
                degree_pos = totals / max(1.0, float(d))
                coord_pos = exps[:, 0].astype(np.float64) / max(1.0, float(d)) if n else 0.0
                term_scales = np.power(10.0, 6.0 * (degree_pos - 0.5) + 2.0 * (coord_pos - 0.5))
                row_scales = np.power(10.0, np.linspace(-2.0, 2.0, n))
                coeff = coeff * term_scales[None, :] * row_scales[:, None]
                metadata.update({
                    "term_scale_min": float(np.min(term_scales)),
                    "term_scale_max": float(np.max(term_scales)),
                    "row_scale_min": float(np.min(row_scales)),
                    "row_scale_max": float(np.max(row_scales)),
                })

        elif family == "degree_shell_ks":
            full_exps = ENGINE.compositions_leq(d, n)
            totals = np.sum(full_exps, axis=1)
            mask = (totals == 0) | (totals == 1) | (totals == d)
            exps = full_exps[mask]
            weights = ENGINE.multinomial_kostlan_weights(exps, d)
            coeff = _complex_gaussian(rng, (n, int(exps.shape[0]))) * weights[None, :]
            metadata.update({
                "support": "degree-shell",
                "degrees": sorted(set(int(x) for x in np.sum(exps, axis=1).tolist())),
                "full_terms_per_poly": int(full_exps.shape[0]),
            })

        elif family == "mixed_degree":
            exps = ENGINE.compositions_leq(d, n)
            weights = np.ones(int(exps.shape[0]), dtype=np.float64)
            coeff = np.zeros((n, int(exps.shape[0])), dtype=np.complex128)
            totals = np.sum(exps, axis=1)
            span = max(1, min(d, n))
            degree_vector = [max(1, d - (i % span)) for i in range(n)]
            for i, row_degree in enumerate(degree_vector):
                mask = totals <= row_degree
                row_weights = ENGINE.multinomial_kostlan_weights(exps[mask], row_degree)
                coeff[i, mask] = _complex_gaussian(rng, (int(np.count_nonzero(mask)),)) * row_weights
            metadata.update({
                "support": "row-mixed-total-degree",
                "degree_vector": [int(x) for x in degree_vector],
                "expected_bezout": int(math.prod(degree_vector)),
            })

        elif family in {"diagonal", "chain", "cyclic"}:
            exps, weights, coeff, metadata = _make_structured_family(family, n, d, rng)

        else:  # pragma: no cover - guarded by _canonical_family
            raise ValueError(f"unsupported family {family!r}")

        coeff = coeff.astype(np.complex128)
        if equation_normalize:
            coeff = _row_normalize(coeff)
        active_counts = _active_counts(coeff)
        obj = cls(
            n=n,
            d=d,
            seed=seed,
            exps=exps,
            coeff=coeff,
            weights=weights,
            equation_normalize=equation_normalize,
            family=family,
            family_description=FAMILY_DESCRIPTIONS[family],
            active_terms=int(sum(active_counts)),
            active_terms_per_equation=active_counts,
            generator_metadata=metadata,
        )
        obj._generation_seconds = ENGINE.now() - t0
        return obj

    @property
    def bezout(self) -> int:
        value = self.generator_metadata.get("expected_bezout")
        if value is not None:
            return int(value)
        return int(self.d ** self.n)

    def stats(self) -> dict[str, Any]:
        stats = super().stats()
        stats.update({
            "family": self.family,
            "active_terms": int(self.active_terms),
            "active_terms_per_equation": list(self.active_terms_per_equation),
            "generator_metadata": dict(self.generator_metadata),
        })
        return stats


def run_case(args: argparse.Namespace, case_raw: str, family: str) -> dict[str, Any]:
    ensure_numpy()
    family = _canonical_family(family)
    t_case = ENGINE.now()
    n, d = ENGINE.parse_case(case_raw)
    system = MultiFamilySystem.make(
        n,
        d,
        seed_index=int(args.seed_index),
        equation_normalize=bool(args.equation_normalize),
        family=family,
        sparse_terms=int(args.sparse_terms),
        sparse_frac=float(args.sparse_frac),
    )
    chart = ENGINE.LinearChart.identity(n, scale=float(args.linear_scale))
    target = ENGINE.TargetTrack(system, chart)

    powers = sorted(set(round(float(x), 16) for x in ENGINE.parse_float_list(args.powers, ENGINE.DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in ENGINE.parse_float_list(args.angles, ENGINE.DEFAULT_ANGLES_DEG)]
    angles_deg = [math.degrees(x) for x in angles]
    radii = ENGINE.parse_float_list(args.rays, ENGINE.DEFAULT_RADII, positive=True)
    gains = ENGINE.parse_float_list(args.startopt_gains, ENGINE.DEFAULT_GAINS, positive=True)

    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    duplicates = 0
    failures = 0
    t_extract = ENGINE.now()

    for trial in range(int(args.pool)):
        if len(roots) >= int(args.count):
            break
        y_raw, geom = ENGINE.mobius_homothety_start(
            n,
            trial,
            system.seed + 0x116000,
            powers,
            angles,
            radii,
            float(args.power_cap),
            roots_found=len(roots),
            duplicates=duplicates,
            failures=failures,
            target_count=int(args.count),
        )
        y0, smeta = ENGINE.startopt(
            target,
            y_raw,
            trial,
            system.seed + 0x116555,
            int(args.startopt_steps),
            int(args.startopt_candidates),
            gains,
            int(args.startopt_micro_epochs),
        )
        loc = ENGINE.pandrosion_corrector(
            target,
            y0,
            max_epochs=int(args.epochs),
            tol=float(args.tol),
            accept=float(args.accept),
            trial_timeout=float(args.trial_timeout),
            line_search=int(args.line_search),
            probe_scale=float(getattr(args, "probe_scale", 0.035)),
            direction_seed=system.seed + 7919 * trial,
        )
        z = chart.z_from_y(loc["y"])
        r_orig = float(np.linalg.norm(system.eval(z)))
        accepted = bool(math.isfinite(r_orig) and r_orig < float(args.accept))
        rec = {
            "trial": int(trial),
            "accepted": accepted,
            "status": loc.get("status"),
            "r1": r_orig,
            "epochs": int(loc.get("epochs", 0)),
            "seconds": float(loc.get("seconds", 0.0)),
            "corrector": loc.get("corrector", "pure-pandrosion-exact-telescopic-slope"),
            "slope_cond": loc.get("slope_cond"),
            **geom,
            **smeta,
        }
        if bool(args.verbose_trials):
            rec["z"] = ENGINE.root_to_json(z)
            rec["y0"] = ENGINE.root_to_json(y0)
        if not accepted:
            failures += 1
            trials.append(rec)
            continue
        dup = ENGINE.cluster_index(roots, z, float(args.cluster_sep))
        if dup is not None:
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = int(dup)
            trials.append(rec)
            continue
        rid = len(roots)
        cond = ENGINE.slope_condition_from_corrector(loc)
        realv = ENGINE.realness(z)
        roots.append({
            "id": rid,
            "source": f"116-multifamily-vectorized-pure-pandrosion:{family}",
            "trial": int(trial),
            "z_complex": np.asarray(z, dtype=np.complex128).copy(),
            "y_complex": np.asarray(loc["y"], dtype=np.complex128).copy(),
            "residual": float(r_orig),
            "realness": float(realv),
            "cond": cond,
            "score": ENGINE.score_root(float(r_orig), realv, cond),
            "epochs": int(loc.get("epochs", 0)),
            "seconds": float(loc.get("seconds", 0.0)),
            "corrector": loc.get("corrector", "pure-pandrosion-exact-telescopic-slope"),
            "slope_cond": loc.get("slope_cond"),
            **geom,
            **smeta,
        })
        rec["status"] = "new-root"
        rec["root_id"] = rid
        trials.append(rec)

    encoded_roots = []
    for r in sorted(roots, key=lambda q: (float(q.get("score", float("inf"))), int(q.get("id", 0)))):
        rr = dict(r)
        zc = rr.pop("z_complex")
        yc = rr.pop("y_complex")
        rr["z"] = ENGINE.root_to_json(zc)
        rr["y"] = ENGINE.root_to_json(yc)
        encoded_roots.append(rr)

    result = {
        "script": "116_pandrosion_multifamily_vectorized_pure_pandrosion.py",
        "autonomous": False,
        "dependencies": {"python_scripts": ["flow/115_pandrosion_vectorized_pure_pandrosion.py"], "numpy": bool(np is not None)},
        "mode": "multifamily-harness/115-vectorized-pure-pandrosion-single-flow",
        "flow_formula": "115 engine: u -> Mobius_Riemann(theta,pole)(u) -> y=Lambda*Mobius(u) -> StartOpt(y) -> z=A*y -> PURE Pandrosion Q_G(a,b) on F(Ay), no Jacobian",
        "case": f"{n},{d}",
        "family": family,
        "family_description": system.family_description,
        "generator_metadata": dict(system.generator_metadata),
        "seed_index": int(args.seed_index),
        "seed": int(system.seed),
        "n": int(n),
        "degree": int(d),
        "terms_per_poly": system.terms_per_poly,
        "active_terms": int(system.active_terms),
        "active_terms_per_equation": list(system.active_terms_per_equation),
        "terms": system.total_terms,
        "bezout": system.bezout,
        "equation_normalize": bool(args.equation_normalize),
        "linear_A": [[ENGINE.cjson(chart.A[i, j]) for j in range(n)] for i in range(n)],
        "parameters": {
            "count": int(args.count),
            "pool": int(args.pool),
            "accept": float(args.accept),
            "tol": float(args.tol),
            "cluster_sep": float(args.cluster_sep),
            "epochs": int(args.epochs),
            "trial_timeout": float(args.trial_timeout),
            "line_search": int(args.line_search),
            "probe_scale": float(args.probe_scale),
            "powers": powers,
            "power_cap": float(args.power_cap),
            "angles_deg": angles_deg,
            "base_rays": radii,
            "startopt_steps": int(args.startopt_steps),
            "startopt_candidates": int(args.startopt_candidates),
            "startopt_gains": gains,
            "startopt_micro_epochs": int(args.startopt_micro_epochs),
            "families": parse_families(args.families),
            "sparse_terms": int(args.sparse_terms),
            "sparse_frac": float(args.sparse_frac),
        },
        "roots": encoded_roots,
        "trials": trials if bool(args.verbose_trials) else trials[: min(len(trials), int(args.keep_trials))],
        "summary": {
            "requested_roots": int(args.count),
            "unique_roots": len(roots),
            "success": bool(len(roots) >= int(args.count)),
            "trials_used": len(trials),
            "duplicates": int(duplicates),
            "failures": int(failures),
            "generation_seconds": system.generation_seconds,
            "extract_seconds": float(ENGINE.now() - t_extract),
            "total_seconds": float(ENGINE.now() - t_case),
            "eval_stats": system.stats(),
        },
    }
    return result


def build_parser() -> argparse.ArgumentParser:
    p = ENGINE.build_parser()
    p.description = "Multifamily harness for the 115 vectorized pure Thales Pandrosion engine."
    p.set_defaults(outdir="verification/116_multifamily_out")
    p.add_argument("--families", default=",".join(DEFAULT_FAMILIES), help="family list, or group: all/random/sparse/structured/dense/fast")
    p.add_argument("--list-families", action="store_true", help="print available generator families and exit")
    p.add_argument("--sparse-terms", type=int, default=0, help="common sparse support size; 0 means automatic")
    p.add_argument("--sparse-frac", type=float, default=0.18, help="automatic sparse support fraction for sparse random families")
    return p


def _output_path(args: argparse.Namespace, cases: Sequence[str], families: Sequence[str]) -> Path:
    if args.out:
        return Path(args.out)
    first_case = cases[0].replace(",", "x") if cases else "case"
    fam = families[0] if len(families) == 1 else "multi"
    return Path(args.outdir) / f"116_multifamily_{first_case}_{fam}.json"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(args.list_families):
        print("116 multifamily generators")
        for name in DEFAULT_FAMILIES:
            print(f"{name}: {FAMILY_DESCRIPTIONS[name]}")
        print("groups: " + ", ".join(f"{k}={','.join(v)}" for k, v in FAMILY_GROUPS.items()))
        return 0

    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    families = parse_families(args.families)
    outputs = [run_case(args, case, family) for case in cases for family in families]
    if len(outputs) == 1:
        final: dict[str, Any] = outputs[0]
    else:
        final = {
            "script": "116_pandrosion_multifamily_vectorized_pure_pandrosion.py",
            "dependencies": {"python_scripts": ["flow/115_pandrosion_vectorized_pure_pandrosion.py"], "numpy": bool(np is not None)},
            "families": list(families),
            "cases": list(cases),
            "runs": outputs,
            "summary": {
                "runs": len(outputs),
                "successes": int(sum(1 for r in outputs if r.get("summary", {}).get("success"))),
                "total_roots": int(sum(int(r.get("summary", {}).get("unique_roots", 0)) for r in outputs)),
            },
        }

    out = _output_path(args, cases, families)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2), encoding="utf-8")

    print("=" * 120, flush=True)
    print("116 multifamily VECTORIZED PURE Pandrosion harness using 115 engine", flush=True)
    print("115 flow/corrector reused; only the multivariate system family generator changes.", flush=True)
    print("=" * 120, flush=True)
    for r in outputs:
        s = r["summary"]
        print(
            f"family={r['family']} case=ks({r['n']},{r['degree']}), seed={r['seed']}, "
            f"active_terms={r['active_terms']}, support_terms={r['terms']}, Bezout={r['bezout']}",
            flush=True,
        )
        print(
            f"roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} "
            f"trials={s['trials_used']} duplicates={s['duplicates']} failures={s['failures']}",
            flush=True,
        )
        print(
            f"seconds: generation={s['generation_seconds']:.2f}, extract={s['extract_seconds']:.2f}, "
            f"total={s['total_seconds']:.2f}",
            flush=True,
        )
        if r.get("roots"):
            best = r["roots"][0]
            print(
                f"best_root: residual={float(best.get('residual', float('inf'))):.3e}, "
                f"trial={best.get('trial')}, Lambda={best.get('homothety')}, "
                f"theta={best.get('theta_deg')}, startopt_ratio={best.get('startopt_ratio')}",
                flush=True,
            )
    print(f"out={out}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    finally:
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except Exception:
            pass
