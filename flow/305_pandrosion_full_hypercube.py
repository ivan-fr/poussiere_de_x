#!/usr/bin/env python3
"""
305_pandrosion_hypercube_quartic_numpy_engine.py

Autonomous HYPERCUBE TENSOR INVERSE JET Pandrosion root extractor.

This is the fully assembled "Flagship" engine. It combines:
1. The Full Universal Atlas (from 304) to handle scale and Riemann inversion.
2. The Hypercube Tensor Inverse Jet corrector (Order 4) to bypass 1D secants
   and Halley formulations, extracting the Jacobian via spatial least squares.

Dependencies: Python stdlib + NumPy. No other files required.
"""
from __future__ import annotations

import argparse
import cmath
import dataclasses
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Iterable, List, Optional, Sequence, Tuple

def _bootstrap_numpy_path() -> None:
    import glob as _glob
    candidates = []
    for pat in (
        "/mnt/data/venv/lib/python*/site-packages",
        "/usr/local/lib/python*/site-packages",
        "/usr/lib/python*/dist-packages",
        "/usr/lib/python*/site-packages",
    ):
        candidates.extend(_glob.glob(pat))
    for path in candidates:
        if path not in sys.path:
            sys.path.append(path)

try:
    import numpy as np
except Exception:  
    _bootstrap_numpy_path()
    try:
        import numpy as np
    except Exception as exc2:  
        np = None
        _NUMPY_IMPORT_ERROR = exc2
    else:
        _NUMPY_IMPORT_ERROR = None
else:
    _NUMPY_IMPORT_ERROR = None


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def now() -> float:
    return time.time()

def cjson(z: complex) -> list[float]:
    return [float(complex(z).real), float(complex(z).imag)]

def root_to_json(z: Sequence[complex]) -> list[list[float]]:
    arr = np.asarray(z, dtype=np.complex128) if np is not None else z
    if np is not None and getattr(arr, "ndim", 1) == 0:
        arr = arr.reshape(1)
    return [cjson(v) for v in arr]

def parse_case(raw: str) -> tuple[int, int]:
    s = str(raw).strip().lower().replace("x", ",").replace(":", ",")
    parts = [p.strip() for p in s.split(",") if p.strip()]
    return int(parts[0]), int(parts[1])

def parse_float_list(raw: Optional[str], default: Sequence[float], positive: bool = False) -> list[float]:
    if raw is None or str(raw).strip() == "":
        return list(default)
    vals: list[float] = []
    for part in str(raw).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        try:
            x = float(part)
            if math.isfinite(x) and ((not positive) or x > 0):
                vals.append(x)
        except Exception:
            pass
    return vals or list(default)

def splitmix64(x: int) -> int:
    x = (int(x) + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (x ^ (x >> 31)) & 0xFFFFFFFFFFFFFFFF

def u01(x: int) -> float:
    return ((splitmix64(x) >> 11) & ((1 << 53) - 1)) / float(1 << 53)

def phase(theta: float) -> complex:
    return complex(math.cos(theta), math.sin(theta))

def ensure_numpy() -> None:
    if np is None:
        raise RuntimeError("NumPy is required.")


# ---------------------------------------------------------------------------
# Polynomial system generation/evaluation
# ---------------------------------------------------------------------------

def compositions_leq(d: int, n: int) -> "np.ndarray":
    ensure_numpy()
    out: list[tuple[int, ...]] = []
    def rec(pos: int, remaining: int, cur: list[int]) -> None:
        if pos == n - 1:
            for k in range(remaining + 1):
                out.append(tuple(cur + [k]))
            return
        for k in range(remaining + 1):
            cur.append(k)
            rec(pos + 1, remaining - k, cur)
            cur.pop()
    rec(0, d, [])
    return np.asarray(out, dtype=np.int16 if d < 32767 else np.int32)

def multinomial_kostlan_weights(exps: "np.ndarray", d: int) -> "np.ndarray":
    ensure_numpy()
    totals = np.sum(exps, axis=1).astype(np.int64)
    logfac = np.zeros(d + 1, dtype=np.float64)
    acc = 0.0
    for k in range(1, d + 1):
        acc += math.log(k)
        logfac[k] = acc
    logs = logfac[d] - logfac[d - totals]
    for j in range(exps.shape[1]):
        logs -= logfac[exps[:, j].astype(np.int64)]
    return np.exp(0.5 * logs)

def stable_seed(n: int, d: int, seed_index: int = 0) -> int:
    return int(splitmix64(0x3050000 + 1000003 * n + 9176 * d + 97 * seed_index) & 0x7FFFFFFF)

@dataclasses.dataclass
class DenseKostlanSystem:
    n: int
    d: int
    seed: int
    exps: Any
    coeff: Any
    weights: Any
    eval_count: int = 0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0) -> "DenseKostlanSystem":
        ensure_numpy()
        exps = compositions_leq(d, n)
        weights = multinomial_kostlan_weights(exps, d)
        seed = stable_seed(n, d, seed_index)
        rng = np.random.default_rng(seed)
        coeff = (rng.standard_normal((n, exps.shape[0])) + 1j * rng.standard_normal((n, exps.shape[0]))) / math.sqrt(2.0)
        coeff = coeff * weights[None, :]
        row_norm = np.linalg.norm(coeff, axis=1)
        row_norm = np.where(row_norm > 0, row_norm, 1.0)
        coeff = coeff / row_norm[:, None]
        return cls(n=n, d=d, seed=seed, exps=exps, coeff=coeff.astype(np.complex128), weights=weights)

    def eval_batch(self, Z: "np.ndarray") -> "np.ndarray":
        ensure_numpy()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        B = int(ZZ.shape[0])
        M = self.exps.shape[0]
        mon = np.ones((B, M), dtype=np.complex128)
        for j in range(self.n):
            p = np.empty((B, self.d + 1), dtype=np.complex128)
            p[:, 0] = 1.0 + 0.0j
            if self.d > 0:
                p[:, 1] = ZZ[:, j]
                for k in range(2, self.d + 1):
                    p[:, k] = p[:, k - 1] * ZZ[:, j]
            mon *= p[:, self.exps[:, j]]
        F = mon @ self.coeff.T
        self.eval_count += B
        return F[0] if Z.ndim == 1 else F

@dataclasses.dataclass
class TargetTrack:
    system: DenseKostlanSystem
    def eval(self, y: Sequence[complex]) -> "np.ndarray":
        return self.system.eval_batch(y)
    def eval_batch(self, Y: "np.ndarray") -> "np.ndarray":
        return self.system.eval_batch(Y)
    def residual(self, y: Sequence[complex]) -> float:
        try:
            rr = float(np.linalg.norm(self.eval(y)))
            return rr if math.isfinite(rr) else float("inf")
        except Exception:
            return float("inf")
    def residuals_batch(self, Y: "np.ndarray") -> "np.ndarray":
        try:
            F = self.eval_batch(Y)
            rr = np.linalg.norm(F, axis=1)
            rr[~np.isfinite(rr)] = np.inf
            return np.asarray(rr, dtype=float)
        except Exception:
            return np.full(int(np.asarray(Y).shape[0]), float("inf"), dtype=float)

# ---------------------------------------------------------------------------
# Full Universal Atlas (from 304, adapted for 305 standalone)
# ---------------------------------------------------------------------------

def _log_stability_energy(y: Any, eps: float = 1e-12) -> float:
    ensure_numpy()
    yy = np.asarray(y, dtype=np.complex128).ravel()
    if yy.size == 0:
        return 0.0
    val = np.log1p(np.abs(yy) + float(eps))
    val[~np.isfinite(val)] = 0.0
    return float(np.linalg.norm(val) / math.sqrt(max(1, yy.size)))

def _universal_complex_cell(n: int, idx: int, layer: int, radius: float) -> "np.ndarray":
    ensure_numpy()
    phi = 0.6180339887498948482
    vals = np.empty(n, dtype=np.complex128)
    mode = idx % 4
    if mode == 0:
        j0 = (idx // 4 + layer) % max(1, n)
        for j in range(n):
            amp = 1.0 if j == j0 else 1.0 / math.sqrt(max(1, n))
            ang = 2.0 * math.pi * ((idx + 1) * (j + 1) * phi + 0.137 * layer)
            vals[j] = amp * phase(ang)
    elif mode == 1:
        for j in range(n):
            h = splitmix64(0x304C0BE + 65537 * idx + 4099 * layer + 193 * (j + 1))
            q = int(4 * u01(h)) % 4
            vals[j] = [1.0, 1j, -1.0, -1j][q]
    elif mode == 2:
        for j in range(n):
            amp = math.exp(0.35 * math.sin((idx + 1) * (j + 1) * 1.324717957244746 + layer))
            ang = 2.0 * math.pi * (((idx + 1) * phi + (j + 1) * phi * phi + 0.071 * layer) % 1.0)
            vals[j] = amp * phase(ang)
    else:
        for j in range(n):
            amp = 1.0 if ((idx + j + layer) % 2 == 0) else 0.35
            ang = 2.0 * math.pi * (((idx + 3) * (j + 5) * 0.41421356237 + 0.19 * layer) % 1.0)
            vals[j] = amp * phase(ang)
    nm = max(1e-300, float(np.linalg.norm(vals)))
    return np.asarray(vals / nm * (float(radius) * math.sqrt(max(1, n))), dtype=np.complex128)

def universal_atlas_start(
    target: TargetTrack,
    n: int,
    trial: int,
    cells: int = 16,
    shells: int = 5,
    probe_radius: float = 0.14,
    descent_min: float = 1.02,
    gap_min: float = 1e-10,
    log_max: float = 80.0
) -> np.ndarray:
    ensure_numpy()
    base_shells = [0.05, 0.12, 0.27, 0.60, 1.0, 1.7, 3.0, 5.5, 10.0, 18.0]
    admissible = []
    
    for kk in range(cells):
        atlas_idx = trial * cells + kk
        layer = atlas_idx // max(1, cells)
        radius = base_shells[layer % len(base_shells)]
        
        y = _universal_complex_cell(n, atlas_idx, layer, radius)
        
        try:
            f0 = target.eval(y)
            r0 = float(np.linalg.norm(f0))
        except Exception:
            continue
            
        if not math.isfinite(r0):
            continue
            
        log0 = _log_stability_energy(y)
        if not math.isfinite(log0) or log0 > float(log_max):
            continue
            
        ynorm = max(1.0, float(np.linalg.norm(y)))
        probes = [0.92 * y, 1.08 * y, y * phase(float(probe_radius))]
        ej = np.zeros(n, dtype=np.complex128)
        ej[(atlas_idx + layer) % max(1, n)] = 1.0
        probes.append(y + float(probe_radius) * ynorm * ej)
        probes.append(y - float(probe_radius) * ynorm * ej)
        P = np.asarray(probes, dtype=np.complex128)
        
        try:
            FP = target.eval_batch(P)
            RP = np.linalg.norm(FP, axis=1)
        except Exception:
            continue
            
        if not np.any(np.isfinite(RP)):
            continue
            
        min_r = float(np.nanmin(RP))
        descent = (r0 / min_r) if math.isfinite(min_r) and min_r > 0 else 0.0
        
        Fdiff = np.linalg.norm(FP - f0[None, :], axis=1)
        gap = Fdiff / (RP + r0 + 1e-300)
        finite_gap = gap[np.isfinite(gap)]
        min_gap = float(np.nanmin(finite_gap)) if finite_gap.size else 0.0
        
        if descent >= float(descent_min) and min_gap >= float(gap_min):
            rep = y.copy()
            for ii, rr in enumerate(RP):
                if math.isfinite(float(rr)) and float(rr) < r0 / float(descent_min):
                    rep = P[int(ii)].copy()
                    break
            admissible.append(rep)

    if admissible:
        pick = trial % len(admissible)
        return np.asarray(admissible[pick], dtype=np.complex128)

    return _universal_complex_cell(n, trial, 0, 1.0)


# ---------------------------------------------------------------------------
# 305: Quartic Hypercube Tensor Inverse Jet 
# ---------------------------------------------------------------------------

def _hypercube_quartic_inverse_jet_305(target: TargetTrack, y: np.ndarray, f: np.ndarray, ep: int, seed: int) -> tuple[np.ndarray, dict]:
    """
    Implements Paper 014 (All-Order Tensor Inverse Jets) up to degree 3 (Quartic)
    using Paper 017 (Finite-Probe Reconstruction).
    """
    ensure_numpy()
    n = len(y)
    M = max(2 * n + 4, 16)
    ynorm = max(1.0, float(np.linalg.norm(y)))
    h_cloud = 1e-5 * ynorm

    rng = np.random.default_rng(seed + ep * 1337)
    Y_cloud = rng.choice([-1.0, 1.0], size=(M, n))
    Y_eval = y[None, :] + h_cloud * Y_cloud
    F_cloud = target.eval_batch(Y_eval)

    dY = h_cloud * Y_cloud
    dF = F_cloud - f[None, :]
    A, _, _, _ = np.linalg.lstsq(dY, dF, rcond=None)
    A = A.T 

    if not np.all(np.isfinite(A)):
        return np.zeros_like(y), {"error": "nonfinite_jacobian", "order": 0}

    try:
        delta1 = np.linalg.solve(A, -f)
    except Exception:
        return np.zeros_like(y), {"error": "singular_jacobian", "order": 0}

    dnorm = float(np.linalg.norm(delta1))
    if (not math.isfinite(dnorm)) or dnorm < 1e-14 or dnorm > 1e4 * ynorm:
        return delta1, {"order": 1}

    h = max(1e-5, min(1e-2, 0.05 * ynorm / dnorm))
    p1 = target.eval(y + h * delta1)
    m1 = target.eval(y - h * delta1)
    p2 = target.eval(y + 2 * h * delta1)
    m2 = target.eval(y - 2 * h * delta1)

    D2 = (-p2 + 16.0 * p1 - 30.0 * f + 16.0 * m1 - m2) / (12.0 * h**2)
    D3 = (p2 - 2.0 * p1 + 2.0 * m1 - m2) / (2.0 * h**3)

    try:
        delta2 = np.linalg.solve(A, -0.5 * D2)
    except Exception:
        return delta1, {"order": 1, "error": "singular_delta2"}
    d2norm = float(np.linalg.norm(delta2))
    if (not math.isfinite(d2norm)) or d2norm > 2.0 * max(dnorm, ynorm):
        return delta1, {"order": 1, "error": "delta2_rejected"}

    f_d2 = target.eval(y + h * delta2)
    f_d1_d2 = target.eval(y + h * delta1 + h * delta2)
    G2_cross = (f_d1_d2 - p1 - f_d2 + f) / (h**2)

    try:
        delta3 = np.linalg.solve(A, -(G2_cross + (1.0 / 6.0) * D3))
    except Exception:
        return delta1 + delta2, {"order": 2, "error": "singular_delta3"}
    d3norm = float(np.linalg.norm(delta3))
    if (not math.isfinite(d3norm)) or d3norm > 2.0 * max(dnorm, ynorm):
        return delta1 + delta2, {"order": 2, "error": "delta3_rejected"}

    delta = delta1 + delta2 + delta3
    
    meta = {
        "order": 4, 
        "d1_norm": float(np.linalg.norm(delta1)),
        "d2_norm": float(np.linalg.norm(delta2)),
        "d3_norm": float(np.linalg.norm(delta3)),
        "hypercube_nodes": M
    }
    return delta, meta


def _line_lambdas() -> list[float]:
    return [1.0, 0.75, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03125, 0.02]

def tensor_corrector(
    target: TargetTrack,
    y0: Sequence[complex],
    max_epochs: int,
    tol: float,
    accept: float,
    seed: int
) -> dict[str, Any]:
    ensure_numpy()
    y = np.asarray(y0, dtype=np.complex128).copy()
    best_y = y.copy()
    
    try:
        f = target.eval(y)
        best_r = float(np.linalg.norm(f))
    except Exception:
        return {"ok": False, "status": "eval-error", "residual": float("inf"), "y": y}

    ok = False
    status = "started"
    lambdas = _line_lambdas()
    
    for ep in range(int(max_epochs)):
        if best_r <= max(tol, accept):
            ok = True
            status = "converged"
            break

        delta, meta = _hypercube_quartic_inverse_jet_305(target, y, f, ep, seed)
        if not np.all(np.isfinite(delta)):
            status = "nonfinite-step"
            break
        
        L = np.asarray(lambdas, dtype=float)
        Ytry = y[None, :] + L[:, None] * delta[None, :]
        Rtry = target.residuals_batch(Ytry)
        
        better = np.isfinite(Rtry) & (Rtry < best_r)
        if not np.any(better):
            status = "no-descent"
            break
            
        scores = np.where(better, Rtry, np.inf)
        idx = int(np.nanargmin(scores))
        
        best_r = float(Rtry[idx])
        y = Ytry[idx].copy()
        best_y = y.copy()
        f = target.eval(y)

    else:
        status = "max-epochs"

    if best_r <= accept:
        ok = True
        status = "converged"

    return {
        "ok": ok,
        "status": status,
        "epochs": ep + 1,
        "residual": best_r,
        "y": best_y,
        "corrector": "305-hypercube-quartic-inverse-jet"
    }

# ---------------------------------------------------------------------------
# Main extractor
# ---------------------------------------------------------------------------

def run_case(n: int, d: int, count: int, pool: int, epochs: int, accept: float, seed_index: int = 0) -> dict:
    t0 = now()
    system = DenseKostlanSystem.make(n, d, seed_index=seed_index)
    target = TargetTrack(system)
    
    roots = []
    trials_used = 0
    failures = 0
    duplicates = 0
    for trial in range(pool):
        if len(roots) >= count:
            break
        trials_used += 1
        
        y0 = universal_atlas_start(target, n, trial)
        loc = tensor_corrector(target, y0, epochs, 1e-12, accept, system.seed + trial)
        
        if loc["ok"]:
            z = loc["y"]
            is_dup = False
            for r in roots:
                if np.linalg.norm(z - r["z"]) < 1e-6:
                    is_dup = True
                    break
            if not is_dup:
                roots.append({
                    "id": len(roots),
                    "z": z,
                    "residual": loc["residual"],
                    "epochs": loc["epochs"]
                })
            else:
                duplicates += 1
        else:
            failures += 1
                
    return {
        "script": "305_pandrosion_full_hypercube.py",
        "mode": "305-hypercube-quartic-inverse-jet",
        "case": f"{n},{d}",
        "n": int(n),
        "degree": int(d),
        "seed_index": int(seed_index),
        "seed": int(system.seed),
        "success": len(roots) >= count,
        "roots_found": len(roots),
        "summary": {
            "requested_roots": int(count),
            "unique_roots": len(roots),
            "success": bool(len(roots) >= int(count)),
            "trials_used": int(trials_used),
            "duplicates": int(duplicates),
            "failures": int(failures),
            "total_seconds": float(now() - t0),
        },
        "eval_count": system.eval_count,
        "roots": [{"id": r["id"], "residual": r["residual"], "epochs": r["epochs"], "z": root_to_json(r["z"])} for r in roots]
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default="2,4; 3,3", help="n,d")
    parser.add_argument("--count", type=int, default=4)
    parser.add_argument("--pool", type=int, default=512)
    parser.add_argument("--epochs", type=int, default=24)
    parser.add_argument("--accept", type=float, default=1e-8)
    parser.add_argument("--seed-index", type=int, default=0)
    parser.add_argument("--out", default="305_hypercube_quartic.json")
    args = parser.parse_args()
    
    cases = [c.strip() for c in args.cases.replace("|", ";").split(";")]
    results = []
    
    print("=" * 100)
    print("305 autonomous HYPERCUBE QUARTIC INVERSE JET NumPy Engine")
    print("Combines Universal Atlas (304) with M-dimensional Hypercube Spatial Tensor Reconstruction (305).")
    print("=" * 100)
    
    for c in cases:
        n, d = parse_case(c)
        res = run_case(n, d, args.count, args.pool, args.epochs, args.accept, args.seed_index)
        s = res["summary"]
        print(f"Case ks({n},{d}) | Roots: {res['roots_found']}/{args.count} | Success: {res['success']} | Trials: {s['trials_used']} | Evals: {res['eval_count']} | Seconds: {s['total_seconds']:.2f}")
        results.append(res)
        
    Path(args.out).write_text(json.dumps(results[0] if len(results) == 1 else {"script": "305_pandrosion_full_hypercube.py", "cases": results}, indent=2))

if __name__ == "__main__":
    main()
