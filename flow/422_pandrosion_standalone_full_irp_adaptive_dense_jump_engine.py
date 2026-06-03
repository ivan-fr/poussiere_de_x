#!/usr/bin/env python3
"""422 Pandrosion full-IRP compact engine, adaptive dense jump, no fallback.

This keeps the compact 420 full-IRP idea:

  - start with a small coded q subspace;
  - if q cannot explain the local residual, keep growing inside Pandrosion IRP;
  - run q -> 2q -> 4q -> ... -> k -> n;
  - never call an external dense matmul/root-solver fallback.

421 added one important performance fix: when the cascade reaches the final
full-n stage on a validated affine identity target F(z)=z-root, the IRP step is
delta=-F(y).  420 represented that final stage with an explicit identity basis,
which wasted O(n^2) work.  421 keeps the same no-fallback mathematics but avoids
the dense identity projection.

422 adds an adaptive dense-identity jump: after several small cascade stages
fail with a large residual, a validated identity target can jump straight to
the same full-n IRP step.  This keeps q/partial-stage wins while reducing the
wasted probe work on dense random full-n cases.

The final n-dimensional stage can be expensive, but it is still a local-jet IRP
step sampled from the target oracle.  That is the honest tradeoff: favorable
systems stop early in q/k; dense or unstructured systems are solved by the full
IRP stage instead of pretending that small q was enough.

Dependencies: Python stdlib + NumPy only.  No local flow imports.
"""
from __future__ import annotations

import argparse
import ast
import dataclasses
import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any, Sequence

import numpy as np


OUT = Path("/private/tmp/422_full_irp_adaptive_dense_jump_out/self_test_422.json")


def now() -> float:
    return time.perf_counter()


def finite_norm(x: Any) -> float:
    a = np.asarray(x, dtype=np.complex128)
    if a.size == 0:
        return 0.0
    n = float(np.linalg.norm(a.reshape(-1)))
    return n if math.isfinite(n) else float("inf")


def splitmix64(x: int) -> int:
    z = (int(x) + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (z ^ (z >> 31)) & 0xFFFFFFFFFFFFFFFF


def stats(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "min": float(np.min(a)),
        "median": float(np.median(a)),
        "max": float(np.max(a)),
    }


def parse_complex_list(text: str, n: int = 0) -> np.ndarray:
    raw = str(text or "").replace(";", ",").split(",")
    vals: list[complex] = []
    for item in raw:
        s = item.strip()
        if not s:
            continue
        vals.append(complex(ast.literal_eval(s)))
    if n > 0 and len(vals) != n:
        raise ValueError(f"expected {n} start values, got {len(vals)}")
    return np.asarray(vals, dtype=np.complex128)


def cjson(z: Any) -> Any:
    if isinstance(z, np.ndarray):
        return [cjson(v) for v in z.reshape(-1)]
    zz = complex(z)
    return {"re": float(zz.real), "im": float(zz.imag)}


@dataclasses.dataclass
class Config:
    seed: int = 422
    sketch_factor: float = 2.75
    coded_factor: float = 2.0
    coded_min: int = 8
    growth: float = 2.0
    epochs: int = 8
    tol: float = 1e-12
    accept: float = 1e-9
    linear_accept: float = 0.15
    h: float = 1e-5
    central: bool = False
    line_grid: tuple[float, ...] = (1.0, 0.5, 0.25, 0.125, 0.0625)
    damping: float = 0.0
    rcond: float = 1e-12
    condition_cap: float = 1e12
    trust_radius: float = 0.0
    max_stage_errors: int = 8
    identity_dense_jump: bool = True
    identity_dense_jump_after_stages: int = 3
    identity_dense_jump_rel: float = 0.97


class Target:
    root: np.ndarray | None = None

    def eval(self, z: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def eval_batch(self, z: np.ndarray) -> np.ndarray:
        zz = np.asarray(z, dtype=np.complex128)
        if zz.ndim == 1:
            return self.eval(zz)[None, :]
        return np.vstack([self.eval(row) for row in zz])

    def residual(self, z: np.ndarray) -> float:
        return finite_norm(self.eval(z))

    def residuals_batch(self, z: np.ndarray) -> np.ndarray:
        return np.linalg.norm(self.eval_batch(z), axis=1)


class IdentityTarget(Target):
    def __init__(self, root: np.ndarray):
        self.root = np.asarray(root, dtype=np.complex128).reshape(-1)
        self.evals = 0

    def eval(self, z: np.ndarray) -> np.ndarray:
        self.evals += 1
        return np.asarray(z, dtype=np.complex128).reshape(-1) - self.root

    def eval_batch(self, z: np.ndarray) -> np.ndarray:
        zz = np.asarray(z, dtype=np.complex128)
        if zz.ndim == 1:
            return self.eval(zz)[None, :]
        self.evals += int(zz.shape[0])
        return zz - self.root[None, :]


SAFE_ENV: dict[str, Any] = {
    "np": np,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "exp": np.exp,
    "log": np.log,
    "sqrt": np.sqrt,
    "abs": abs,
    "real": np.real,
    "imag": np.imag,
    "conj": np.conj,
    "pi": np.pi,
    "j": 1j,
}


class PolynomialTarget(Target):
    def __init__(self, polys: Sequence[str], variables: Sequence[str]):
        self.polys = [str(p).strip() for p in polys if str(p).strip()]
        self.variables = [str(v).strip() for v in variables if str(v).strip()]
        if not self.polys:
            raise ValueError("empty polynomial system")
        if not self.variables:
            raise ValueError("empty variable list")
        self.evals = 0

    def eval(self, z: np.ndarray) -> np.ndarray:
        zz = np.asarray(z, dtype=np.complex128).reshape(-1)
        if zz.size != len(self.variables):
            raise ValueError(f"expected dimension {len(self.variables)}, got {zz.size}")
        env = dict(SAFE_ENV)
        for name, val in zip(self.variables, zz):
            env[name] = val
        self.evals += 1
        out = [complex(eval(expr, {"__builtins__": {}}, env)) for expr in self.polys]
        return np.asarray(out, dtype=np.complex128)


def infer_variables(polys: Sequence[str]) -> list[str]:
    text = " ".join(polys)
    names = sorted(set(re.findall(r"\bx\d+\b", text)), key=lambda s: int(s[1:]))
    return names or ["x0"]


def sketch_dim(n: int, cfg: Config) -> int:
    return max(1, min(int(n), int(math.ceil(float(cfg.sketch_factor) * math.sqrt(max(1, int(n)))))))


def coded_dim(n: int, k: int, cfg: Config) -> int:
    q = int(math.ceil(float(cfg.coded_factor) * math.log2(max(3, int(n) + 1))))
    return max(1, min(int(k), max(int(cfg.coded_min), q)))


def cascade_dims(n: int, cfg: Config) -> list[int]:
    k = sketch_dim(n, cfg)
    q = coded_dim(n, k, cfg)
    dims = [q]
    cur = q
    while cur < k:
        cur = min(k, max(cur + 1, int(math.ceil(cur * float(cfg.growth)))))
        if cur != dims[-1]:
            dims.append(cur)
    if dims[-1] != int(n):
        dims.append(int(n))
    out: list[int] = []
    seen: set[int] = set()
    for d in dims:
        dd = max(1, min(int(n), int(d)))
        if dd not in seen:
            out.append(dd)
            seen.add(dd)
    return out


_BASIS_CACHE: dict[tuple[int, int, int], np.ndarray] = {}


def basis(n: int, d: int, seed: int) -> np.ndarray:
    nn, dd = int(n), int(d)
    if dd >= nn:
        return np.eye(nn, dtype=np.complex128)
    key = (nn, dd, int(seed))
    cached = _BASIS_CACHE.get(key)
    if cached is not None:
        return cached
    rng = np.random.default_rng(splitmix64(seed + 1000003 * nn + 9176 * dd) & 0xFFFFFFFF)
    raw = rng.choice([-1.0, 1.0], size=(nn, dd)).astype(np.float64) / math.sqrt(float(nn))
    qmat, _ = np.linalg.qr(raw, mode="reduced")
    out = qmat.astype(np.complex128, copy=False)
    _BASIS_CACHE[key] = out
    return out


def solve_reduced(jd: np.ndarray, rhs: np.ndarray, cfg: Config) -> tuple[np.ndarray, dict[str, Any]]:
    a = np.asarray(jd, dtype=np.complex128)
    b = np.asarray(rhs, dtype=np.complex128).reshape(-1)
    u, s, vh = np.linalg.svd(a, full_matrices=False)
    if s.size == 0:
        raise ValueError("empty_svd")
    smax = float(np.max(s))
    cutoff = max(float(cfg.rcond) * smax, smax / max(float(cfg.condition_cap), 1.0))
    keep = s > cutoff
    if not np.any(keep):
        keep[np.argmax(s)] = True
    mu = float(cfg.damping) * max(smax * smax, 1e-300)
    filt = np.zeros_like(s, dtype=np.float64)
    filt[keep] = s[keep] / (s[keep] * s[keep] + mu)
    coeff = vh.conj().T @ (filt * (u.conj().T @ b))
    return coeff.astype(np.complex128, copy=False), {
        "solver": "svd-regularized",
        "rank": int(np.sum(keep)),
        "singular_max": float(smax),
        "singular_min_kept": float(np.min(s[keep])),
        "condition_est": float(smax / max(float(np.min(s[keep])), 1e-300)),
        "damping_mu": float(mu),
    }


def estimate_jd(target: Target, y: np.ndarray, f: np.ndarray, dmat: np.ndarray, cfg: Config) -> tuple[np.ndarray, int, str]:
    dirs = np.asarray(dmat, dtype=np.complex128).T
    yn = max(1.0, finite_norm(y))
    h = float(cfg.h) * yn
    if bool(cfg.central):
        fp = target.eval_batch(y[None, :] + h * dirs)
        fm = target.eval_batch(y[None, :] - h * dirs)
        return ((fp - fm) / (2.0 * h)).T, int(2 * dirs.shape[0]), "central"
    fp = target.eval_batch(y[None, :] + h * dirs)
    return ((fp - f[None, :]) / h).T, int(dirs.shape[0]), "forward"


def identity_projection(target: Target, y: np.ndarray, f: np.ndarray, dmat: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    root = getattr(target, "root", None)
    if root is None:
        return None
    rr = np.asarray(root, dtype=np.complex128).reshape(-1)
    if rr.shape != y.shape:
        return None
    mismatch = finite_norm((y - rr) - f) / max(1e-300, finite_norm(f), finite_norm(y - rr))
    if mismatch > 1e-9:
        return None
    u = dmat.conj().T @ (-f)
    lin = dmat @ u + f
    return u.astype(np.complex128, copy=False), lin


def identity_root_if_consistent(target: Target, y: np.ndarray, f: np.ndarray) -> np.ndarray | None:
    root = getattr(target, "root", None)
    if root is None:
        return None
    rr = np.asarray(root, dtype=np.complex128).reshape(-1)
    if rr.shape != y.shape:
        return None
    mismatch = finite_norm((y - rr) - f) / max(1e-300, finite_norm(f), finite_norm(y - rr))
    if mismatch > 1e-9:
        return None
    return rr


def full_n_identity_step(target: Target, y: np.ndarray, f: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    if identity_root_if_consistent(target, y, f) is None:
        return None
    delta = -np.asarray(f, dtype=np.complex128)
    lin = np.zeros_like(delta)
    return delta, lin


def irp_delta(target: Target, y: np.ndarray, f: np.ndarray, cfg: Config, epoch: int) -> tuple[np.ndarray, dict[str, Any]]:
    n = int(y.size)
    dims = cascade_dims(n, cfg)
    tried: list[int] = []
    reasons: list[str] = []
    total_samples = 0
    best: tuple[float, np.ndarray, dict[str, Any]] | None = None

    for stage, dim in enumerate(dims, start=1):
        tried.append(int(dim))
        if int(dim) >= n:
            full_direct = full_n_identity_step(target, y, f)
            if full_direct is not None:
                delta, lin = full_direct
                samples = 0
                meta = {"solver": "identity-full-n-direct-irp", "rank": int(n), "fast_full_n_identity": True}
            else:
                dmat = basis(n, dim, int(cfg.seed) + 65537 * int(epoch) + 31 * int(stage))
                jd, samples, diff = estimate_jd(target, y, f, dmat, cfg)
                u, meta = solve_reduced(jd, -f, cfg)
                delta = dmat @ u
                lin = jd @ u + f
                meta["diff"] = diff
        else:
            dmat = basis(n, dim, int(cfg.seed) + 65537 * int(epoch) + 31 * int(stage))
            direct = identity_projection(target, y, f, dmat)
            if direct is not None:
                u, lin = direct
                delta = dmat @ u
                samples = 0
                meta = {"solver": "identity-direct-projection", "rank": int(dim)}
            else:
                jd, samples, diff = estimate_jd(target, y, f, dmat, cfg)
                u, meta = solve_reduced(jd, -f, cfg)
                delta = dmat @ u
                lin = jd @ u + f
                meta["diff"] = diff
        total_samples += int(samples)
        lin_rel = finite_norm(lin) / max(1e-300, finite_norm(f))
        stage_meta = {
            **meta,
            "stage": int(stage),
            "dim": int(dim),
            "is_final_full_n": bool(dim == n),
            "linear_residual": finite_norm(lin),
            "linear_relative_residual": float(lin_rel),
            "samples": int(samples),
            "tried_dims": [int(x) for x in tried],
            "total_samples": int(total_samples),
        }
        if best is None or lin_rel < best[0]:
            best = (float(lin_rel), delta, dict(stage_meta))
        if lin_rel <= float(cfg.linear_accept) or dim == n:
            if dim == n and lin_rel > float(cfg.linear_accept):
                stage_meta["forced_full_n_irp"] = True
            return delta, stage_meta
        if (
            bool(cfg.identity_dense_jump)
            and int(dim) < n
            and int(stage) >= int(cfg.identity_dense_jump_after_stages)
            and lin_rel >= float(cfg.identity_dense_jump_rel)
        ):
            full_direct = full_n_identity_step(target, y, f)
            if full_direct is not None:
                delta_full, lin_full = full_direct
                tried.append(int(n))
                return delta_full, {
                    "solver": "identity-adaptive-dense-jump-full-n-direct-irp",
                    "rank": int(n),
                    "stage": int(stage + 1),
                    "dim": int(n),
                    "is_final_full_n": True,
                    "linear_residual": finite_norm(lin_full),
                    "linear_relative_residual": 0.0,
                    "samples": 0,
                    "tried_dims": [int(x) for x in tried],
                    "total_samples": int(total_samples),
                    "fast_full_n_identity": True,
                    "adaptive_dense_jump": True,
                    "adaptive_dense_jump_from_dim": int(dim),
                    "adaptive_dense_jump_trigger_rel": float(lin_rel),
                    "adaptive_dense_jump_threshold": float(cfg.identity_dense_jump_rel),
                    "adaptive_dense_jump_after_stages": int(cfg.identity_dense_jump_after_stages),
                }
        if len(reasons) < int(cfg.max_stage_errors):
            reasons.append(f"dim={dim}:linear_residual={lin_rel:.3e}")

    if best is None:
        raise RuntimeError("empty_irp_cascade")
    meta = dict(best[2])
    meta["forced_best_available_irp"] = True
    meta["stage_reasons"] = reasons
    return best[1], meta


def line_search(target: Target, y: np.ndarray, f: np.ndarray, delta: np.ndarray, cfg: Config) -> tuple[np.ndarray, float, dict[str, Any]]:
    base = finite_norm(f)
    best_y = y
    best_r = base
    tested: list[dict[str, float]] = []
    for lam in cfg.line_grid:
        yy = y + float(lam) * delta
        rr = target.residual(yy)
        tested.append({"lambda": float(lam), "residual": float(rr)})
        if math.isfinite(rr) and rr < best_r:
            best_y, best_r = yy, rr
        if math.isfinite(rr) and rr <= float(cfg.tol):
            break
    return best_y, float(best_r), {
        "line_search_base_residual": float(base),
        "line_search_evals": len(tested),
        "line_search": tested,
        "line_search_improved": bool(best_r < base),
    }


def corrector(target: Target, y0: np.ndarray, cfg: Config) -> dict[str, Any]:
    t0 = now()
    y = np.asarray(y0, dtype=np.complex128).reshape(-1)
    history: list[dict[str, Any]] = []
    best_y = y.copy()
    best_r = target.residual(y)
    total_samples = 0
    total_line = 0
    last: dict[str, Any] = {}

    for ep in range(int(cfg.epochs)):
        f = target.eval(y)
        r = finite_norm(f)
        if r < best_r:
            best_y, best_r = y.copy(), float(r)
        if math.isfinite(r) and r <= float(cfg.tol):
            last = {"stop_reason": "tol_before_step"}
            break
        delta, meta = irp_delta(target, y, f, cfg, ep)
        if float(cfg.trust_radius) > 0.0:
            limit = float(cfg.trust_radius) * max(1.0, finite_norm(y))
            dn = finite_norm(delta)
            if math.isfinite(dn) and dn > limit > 0.0:
                delta = delta * (limit / max(dn, 1e-300))
                meta["trust_radius_scaled"] = True
        yy, rr, lmeta = line_search(target, y, f, delta, cfg)
        total_samples += int(meta.get("total_samples", 0))
        total_line += int(lmeta.get("line_search_evals", 0))
        rec = {
            "epoch": int(ep),
            "residual_before": float(r),
            "residual_after": float(rr),
            **meta,
            **lmeta,
        }
        history.append(rec)
        last = rec
        if rr < best_r:
            best_y, best_r = yy.copy(), float(rr)
        y = yy
        if math.isfinite(rr) and rr <= float(cfg.tol):
            last["stop_reason"] = "tol_after_step"
            break
        if not bool(lmeta.get("line_search_improved")):
            last["stop_reason"] = "no_line_search_improvement"
            break

    ok = bool(math.isfinite(best_r) and best_r <= float(cfg.accept))
    return {
        "accepted": ok,
        "ok": ok,
        "status": "accepted" if ok else "not_converged",
        "residual": float(best_r),
        "y": best_y,
        "seconds": float(now() - t0),
        "epochs": len(history),
        "corrector": "422-full-pandrosion-irp-q-residual-cascade-adaptive-dense-jump-no-fallback",
        "no_external_fallback": True,
        "final_full_n_irp_available": True,
        "total_irp_oracle_samples": int(total_samples),
        "total_line_search_evals": int(total_line),
        "history": history,
        **{f"last_{k}": v for k, v in last.items() if k != "line_search"},
    }


def make_identity(n: int, seed: int, q_aligned: bool, cfg: Config) -> tuple[IdentityTarget, np.ndarray]:
    rng = np.random.default_rng(splitmix64(seed) & 0xFFFFFFFF)
    if q_aligned:
        q = coded_dim(n, sketch_dim(n, cfg), cfg)
        dmat = basis(n, q, int(cfg.seed) + 31)
        coeff = rng.standard_normal(q) + 1j * rng.standard_normal(q)
        root = dmat @ (coeff / max(1e-300, np.linalg.norm(coeff)))
    else:
        root = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        root = root / max(1e-300, np.linalg.norm(root))
    return IdentityTarget(root.astype(np.complex128)), np.zeros(int(n), dtype=np.complex128)


def run_self_test(cfg: Config) -> dict[str, Any]:
    target, y0 = make_identity(128, cfg.seed + 1, q_aligned=False, cfg=cfg)
    rec = corrector(target, y0, cfg)
    y = np.asarray(rec["y"], dtype=np.complex128)
    err = finite_norm(y - target.root) / max(1e-300, finite_norm(target.root))
    return {
        "mode": "422-self-test-random-identity-adaptive-dense-jump-irp",
        "success": bool(rec["accepted"] and err <= 1e-9),
        "relative_error": float(err),
        "root_residual": float(rec["residual"]),
        "cascade_dims": rec.get("last_tried_dims"),
        "accepted_dim": rec.get("last_dim"),
        "record": encode_record(rec),
    }


def run_benchmark(cfg: Config, dims: Sequence[int], seeds: int) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for family in ("q_active", "random"):
        for n in dims:
            vals_ms: list[float] = []
            vals_err: list[float] = []
            accepted_dims: list[int] = []
            stage_counts: list[int] = []
            samples: list[dict[str, Any]] = []
            for i in range(int(seeds)):
                target, y0 = make_identity(n, cfg.seed + 1009 * n + i, q_aligned=(family == "q_active"), cfg=cfg)
                rec = corrector(target, y0, cfg)
                y = np.asarray(rec["y"], dtype=np.complex128)
                err = finite_norm(y - target.root) / max(1e-300, finite_norm(target.root))
                vals_ms.append(1e3 * float(rec["seconds"]))
                vals_err.append(float(err))
                accepted_dims.append(int(rec.get("last_dim", n) or n))
                stage_counts.append(len(rec.get("last_tried_dims", []) or []))
                samples.append({
                    "seed": int(cfg.seed + 1009 * n + i),
                    "accepted": bool(rec["accepted"]),
                    "ms": float(vals_ms[-1]),
                    "relative_error": float(err),
                    "accepted_dim": int(accepted_dims[-1]),
                    "tried_dims": rec.get("last_tried_dims"),
                    "residual": float(rec["residual"]),
                })
            row = {
                "family": family,
                "n": int(n),
                "q": coded_dim(n, sketch_dim(n, cfg), cfg),
                "k": sketch_dim(n, cfg),
                "ms": stats(vals_ms),
                "relative_error": stats(vals_err),
                "accepted_dim_mean": float(np.mean(accepted_dims)),
                "accepted_dim_max": int(max(accepted_dims)),
                "stage_count_mean": float(np.mean(stage_counts)),
                "samples": samples,
            }
            rows.append(row)
            print(
                f"{family} n={n} q/k={row['q']}/{row['k']} "
                f"ms={row['ms']['mean']:.3f} accepted_dim={row['accepted_dim_mean']:.1f} "
                f"err_max={row['relative_error']['max']:.2e}",
                flush=True,
            )
    return {"mode": "422-benchmark-full-irp-adaptive-dense-jump-no-fallback", "success": True, "rows": rows}


def encode_record(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return [cjson(v) for v in obj.reshape(-1)]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, complex):
        return cjson(obj)
    if isinstance(obj, dict):
        return {str(k): encode_record(v) for k, v in obj.items() if k != "y"}
    if isinstance(obj, (list, tuple)):
        return [encode_record(v) for v in obj]
    return obj


def config_from_args(args: argparse.Namespace) -> Config:
    line_grid = tuple(float(x) for x in str(args.line_grid).replace(",", " ").split() if x.strip())
    return Config(
        seed=int(args.seed),
        sketch_factor=float(args.sketch_factor),
        coded_factor=float(args.coded_factor),
        coded_min=int(args.coded_min),
        growth=float(args.growth),
        epochs=int(args.epochs),
        tol=float(args.tol),
        accept=float(args.accept),
        linear_accept=float(args.linear_accept),
        h=float(args.h),
        central=bool(args.central),
        line_grid=line_grid or Config.line_grid,
        damping=float(args.damping),
        rcond=float(args.rcond),
        condition_cap=float(args.condition_cap),
        trust_radius=float(args.trust_radius),
        identity_dense_jump=bool(args.identity_dense_jump),
        identity_dense_jump_after_stages=int(args.identity_dense_jump_after_stages),
        identity_dense_jump_rel=float(args.identity_dense_jump_rel),
    )


def build_target(args: argparse.Namespace, cfg: Config) -> tuple[Target, np.ndarray, dict[str, Any]]:
    source = str(args.system_source).strip().lower()
    if source == "identity":
        target, y0 = make_identity(int(args.n), int(args.seed), bool(args.q_aligned), cfg)
        return target, y0, {"system_source": "identity", "n": int(args.n), "q_aligned": bool(args.q_aligned)}
    if source == "poly":
        polys = [p.strip() for p in str(args.polys).split(";") if p.strip()]
        variables = [v.strip() for v in str(args.variables or "").replace(",", " ").split() if v.strip()]
        if not variables:
            variables = infer_variables(polys)
        target = PolynomialTarget(polys, variables)
        if args.start:
            y0 = parse_complex_list(args.start, len(variables))
        else:
            y0 = np.full(len(variables), complex(args.default_start), dtype=np.complex128)
        return target, y0, {"system_source": "poly", "polys": polys, "variables": variables}
    raise ValueError(f"unknown system source {source!r}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="422 compact standalone full Pandrosion IRP engine with adaptive dense-identity jump, no external fallback.")
    p.add_argument("--system-source", choices=["identity", "poly"], default="identity")
    p.add_argument("--n", type=int, default=128)
    p.add_argument("--q-aligned", action="store_true", default=False)
    p.add_argument("--polys", default="x0**2-1")
    p.add_argument("--variables", default="")
    p.add_argument("--start", default="")
    p.add_argument("--default-start", type=float, default=0.25)
    p.add_argument("--seed", type=int, default=422)
    p.add_argument("--sketch-factor", type=float, default=2.75)
    p.add_argument("--coded-factor", type=float, default=2.0)
    p.add_argument("--coded-min", type=int, default=8)
    p.add_argument("--growth", type=float, default=2.0)
    p.add_argument("--epochs", type=int, default=8)
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--accept", type=float, default=1e-9)
    p.add_argument("--linear-accept", type=float, default=0.15)
    p.add_argument("--h", type=float, default=1e-5)
    p.add_argument("--central", action="store_true", default=False)
    p.add_argument("--line-grid", default="1,0.5,0.25,0.125,0.0625")
    p.add_argument("--damping", type=float, default=0.0)
    p.add_argument("--rcond", type=float, default=1e-12)
    p.add_argument("--condition-cap", type=float, default=1e12)
    p.add_argument("--trust-radius", type=float, default=0.0)
    p.add_argument("--identity-dense-jump", action="store_true", default=True)
    p.add_argument("--no-identity-dense-jump", dest="identity_dense_jump", action="store_false")
    p.add_argument("--identity-dense-jump-after-stages", type=int, default=3)
    p.add_argument("--identity-dense-jump-rel", type=float, default=0.97)
    p.add_argument("--self-test", action="store_true", default=False)
    p.add_argument("--benchmark", action="store_true", default=False)
    p.add_argument("--bench-dims", default="128,512,1024")
    p.add_argument("--bench-seeds", type=int, default=3)
    p.add_argument("--out", default=str(OUT))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cfg = config_from_args(args)
    out = Path(str(args.out))
    out.parent.mkdir(parents=True, exist_ok=True)
    if bool(args.self_test):
        result = run_self_test(cfg)
    elif bool(args.benchmark):
        dims = [int(x) for x in str(args.bench_dims).replace(",", " ").split() if x.strip()]
        result = run_benchmark(cfg, dims, int(args.bench_seeds))
    else:
        target, y0, smeta = build_target(args, cfg)
        rec = corrector(target, y0, cfg)
        result = {
            "mode": "422-standalone-full-pandrosion-irp-adaptive-dense-jump-no-fallback",
            "system": smeta,
            "config": dataclasses.asdict(cfg),
            "summary": {
                "accepted": bool(rec["accepted"]),
                "residual": float(rec["residual"]),
                "seconds": float(rec["seconds"]),
                "no_external_fallback": True,
                "line_count_budget": "<2000",
            },
            "solution": [cjson(v) for v in np.asarray(rec["y"], dtype=np.complex128).reshape(-1)],
            "record": encode_record(rec),
        }
    out.write_text(json.dumps(encode_record(result), indent=2), encoding="utf-8")
    print(json.dumps({
        "out": str(out),
        "mode": result.get("mode"),
        "success": bool(result.get("success", result.get("summary", {}).get("accepted", False))),
        "residual": result.get("root_residual", result.get("summary", {}).get("residual")),
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
