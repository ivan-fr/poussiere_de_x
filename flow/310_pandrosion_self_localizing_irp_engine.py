#!/usr/bin/env python3
"""
310_pandrosion_self_localizing_irp_engine.py

SELF-LOCALIZING STANDALONE NUMPY ENGINE FOR ITERATED RENORMALIZED PANDROSION (IRP).

This file intentionally rebuilds the engine instead of adapting the 304/308
Halley stack.  It keeps only the minimum geometry needed for a universal dense
polynomial extractor:

  - autonomous dense Kostlan polynomial system generator for ks(n,d),
  - NumPy polynomial evaluation and batched evaluation,
  - exact derivative-free telescopic Pandrosion slope Q(a,b),
  - deterministic complex-shell atlas starts,
  - pure IRP local cascade

        geometric renormalization -> raw Pandrosion ->
        geometric renormalization -> raw Pandrosion -> ...

No analytic Jacobian, Hessian, Newton step, Halley/Householder correction,
Steffensen quotient, Broyden update, SciPy, homotopy continuation, or imports
from previous Pandrosion scripts are used.

The only local direction is the raw anchored Pandrosion solve

        Q_G(a,b) delta = -G(a),

where Q_G(a,b) is constructed exactly from finite telescopic monomial sums so
that G(b)-G(a)=Q_G(a,b)(b-a).  Higher-order behavior is sought only through
repeated geometric renormalization of this first-order geometry.

Dependencies: Python stdlib + NumPy only.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


def _bootstrap_numpy_path() -> None:
    import glob as _glob
    for pat in (
        "/mnt/data/venv/lib/python*/site-packages",
        "/usr/local/lib/python*/site-packages",
        "/usr/lib/python*/dist-packages",
        "/usr/lib/python*/site-packages",
    ):
        for path in _glob.glob(pat):
            if path not in sys.path:
                sys.path.append(path)


try:
    import numpy as np
except Exception as exc:  # pragma: no cover
    _bootstrap_numpy_path()
    try:
        import numpy as np
    except Exception as exc2:  # pragma: no cover
        np = None
        _NUMPY_IMPORT_ERROR = exc2
    else:
        _NUMPY_IMPORT_ERROR = None
else:
    _NUMPY_IMPORT_ERROR = None


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def ensure_numpy() -> None:
    if np is None:
        raise RuntimeError(f"NumPy is required; import error={_NUMPY_IMPORT_ERROR!r}")


def now() -> float:
    return time.time()


def cjson(z: complex) -> list[float]:
    zz = complex(z)
    return [float(zz.real), float(zz.imag)]


def root_to_json(z: Sequence[complex]) -> list[list[float]]:
    return [cjson(v) for v in z]


def parse_case(raw: str) -> tuple[int, int]:
    s = str(raw).strip().lower().replace("x", ",").replace(":", ",")
    parts = [p.strip() for p in s.split(",") if p.strip()]
    if len(parts) != 2:
        raise ValueError(f"case must be n,d, got {raw!r}")
    n, d = int(parts[0]), int(parts[1])
    if n <= 0 or d <= 0:
        raise ValueError("n,d must be positive")
    return n, d


def parse_float_list(raw: Optional[str], default: Sequence[float], positive: bool = False) -> list[float]:
    if raw is None or str(raw).strip() == "":
        return list(default)
    out: list[float] = []
    for part in str(raw).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        try:
            x = float(part)
            if math.isfinite(x) and ((not positive) or x > 0):
                out.append(x)
        except Exception:
            pass
    return out or list(default)


def splitmix64(x: int) -> int:
    x = (int(x) + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (x ^ (x >> 31)) & 0xFFFFFFFFFFFFFFFF


def u01(x: int) -> float:
    return ((splitmix64(x) >> 11) & ((1 << 53) - 1)) / float(1 << 53)


def phase(theta: float) -> complex:
    return complex(math.cos(theta), math.sin(theta))


def norm2(v: Sequence[complex]) -> float:
    return float(np.linalg.norm(np.asarray(v, dtype=np.complex128)))


def stable_seed(n: int, d: int, seed_index: int = 0, salt: int = 0) -> int:
    base = 0x50414E44524F5349  # "PANDROSI" fragment
    return int(splitmix64(base + 1000003 * n + 9176 * d + 97 * seed_index + salt) & 0x7FFFFFFF)


def finite_norm(x: Any) -> float:
    try:
        r = float(np.linalg.norm(np.asarray(x, dtype=np.complex128)))
        return r if math.isfinite(r) else float("inf")
    except Exception:
        return float("inf")


# ---------------------------------------------------------------------------
# Dense Kostlan system
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


def kostlan_weights(exps: "np.ndarray", d: int) -> "np.ndarray":
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


@dataclasses.dataclass
class DenseKostlanSystem:
    n: int
    d: int
    seed: int
    exps: Any
    coeff: Any
    weights: Any
    equation_normalize: bool = False
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = False) -> "DenseKostlanSystem":
        ensure_numpy()
        t0 = now()
        exps = compositions_leq(d, n)
        weights = kostlan_weights(exps, d)
        seed = stable_seed(n, d, seed_index)
        rng = np.random.default_rng(seed)
        coeff = (rng.standard_normal((n, exps.shape[0])) + 1j * rng.standard_normal((n, exps.shape[0]))) / math.sqrt(2.0)
        coeff = coeff * weights[None, :]
        if equation_normalize:
            rn = np.linalg.norm(coeff, axis=1)
            rn = np.where(rn > 0, rn, 1.0)
            coeff = coeff / rn[:, None]
        obj = cls(n=n, d=d, seed=seed, exps=exps, coeff=coeff.astype(np.complex128), weights=weights, equation_normalize=equation_normalize)
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return int(self.exps.shape[0])

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def _powers(self, z: "np.ndarray") -> list[Any]:
        tables = []
        for zj in z:
            p = np.empty(self.d + 1, dtype=np.complex128)
            p[0] = 1.0 + 0.0j
            if self.d >= 1:
                p[1] = complex(zj)
                for k in range(2, self.d + 1):
                    p[k] = p[k - 1] * zj
            tables.append(p)
        return tables

    def monomials(self, z: Sequence[complex]) -> "np.ndarray":
        ensure_numpy()
        zz = np.asarray(z, dtype=np.complex128)
        pows = self._powers(zz)
        mon = np.ones(self.terms_per_poly, dtype=np.complex128)
        for j in range(self.n):
            mon *= pows[j][self.exps[:, j]]
        return mon

    def eval(self, z: Sequence[complex]) -> "np.ndarray":
        ensure_numpy()
        t0 = now()
        mon = self.monomials(z)
        f = self.coeff @ mon
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return f

    def monomials_batch(self, Z: "np.ndarray") -> "np.ndarray":
        ensure_numpy()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.monomials(ZZ)[None, :]
        B = int(ZZ.shape[0])
        M = self.terms_per_poly
        mon = np.ones((B, M), dtype=np.complex128)
        for j in range(self.n):
            p = np.empty((B, self.d + 1), dtype=np.complex128)
            p[:, 0] = 1.0 + 0.0j
            if self.d >= 1:
                p[:, 1] = ZZ[:, j]
                for k in range(2, self.d + 1):
                    p[:, k] = p[:, k - 1] * ZZ[:, j]
            mon *= p[:, self.exps[:, j]]
        return mon

    def eval_batch(self, Z: "np.ndarray") -> "np.ndarray":
        ensure_numpy()
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        mon = self.monomials_batch(ZZ)
        F = mon @ self.coeff.T
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return F

    def _slope_power_table(self, aj: complex, bj: complex, pows_b_j: "np.ndarray") -> "np.ndarray":
        """slope[m] = sum_{r=0}^{m-1} bj^(m-1-r) aj^r, slope[0]=0."""
        slope = np.empty(self.d + 1, dtype=np.complex128)
        slope[0] = 0.0 + 0.0j
        acc = 0.0 + 0.0j
        ajc = complex(aj)
        for m in range(1, self.d + 1):
            acc = pows_b_j[m - 1] + ajc * acc
            slope[m] = acc
        return slope

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> "np.ndarray":
        """Exact derivative-free telescopic Pandrosion slope matrix Q(a,b).

        For all dense polynomial rows G_i, Q satisfies

            G(b)-G(a) = Q(a,b) @ (b-a).

        The construction uses only finite monomial identities and not an analytic
        Jacobian formula.  The removable case b=a is filled by the same finite
        recurrence, so no symbolic derivative is called.
        """
        ensure_numpy()
        t0 = now()
        aa = np.asarray(a, dtype=np.complex128)
        bb = np.asarray(b, dtype=np.complex128)
        pows_a = self._powers(aa)
        pows_b = self._powers(bb)
        n = self.n
        M = self.terms_per_poly

        pa_cols = [pows_a[j][self.exps[:, j]] for j in range(n)]
        pb_cols = [pows_b[j][self.exps[:, j]] for j in range(n)]

        prefix_b = [None] * (n + 1)
        suffix_a = [None] * (n + 1)
        prefix_b[0] = np.ones(M, dtype=np.complex128)
        for j in range(n):
            prefix_b[j + 1] = prefix_b[j] * pb_cols[j]
        suffix_a[n] = np.ones(M, dtype=np.complex128)
        for j in range(n - 1, -1, -1):
            suffix_a[j] = suffix_a[j + 1] * pa_cols[j]

        terms = np.empty((M, n), dtype=np.complex128)
        for j in range(n):
            slope_table = self._slope_power_table(aa[j], bb[j], pows_b[j])
            terms[:, j] = prefix_b[j] * suffix_a[j + 1] * slope_table[self.exps[:, j]]
        Q = self.coeff @ terms

        self.slope_count += 1
        self.seconds_slope += now() - t0
        return Q

    def residual(self, z: Sequence[complex]) -> float:
        try:
            return finite_norm(self.eval(z))
        except Exception:
            return float("inf")

    def residuals_batch(self, Z: "np.ndarray") -> "np.ndarray":
        try:
            F = self.eval_batch(Z)
            r = np.linalg.norm(F, axis=1).astype(float)
            r[~np.isfinite(r)] = np.inf
            return r
        except Exception:
            ZZ = np.asarray(Z, dtype=np.complex128)
            return np.full(int(ZZ.shape[0]) if ZZ.ndim > 1 else 1, np.inf, dtype=float)

    def stats(self) -> dict[str, Any]:
        return {
            "eval_count": int(self.eval_count),
            "slope_count": int(self.slope_count),
            "seconds_eval": float(self.seconds_eval),
            "seconds_slope": float(self.seconds_slope),
            "terms_per_poly": int(self.terms_per_poly),
            "total_terms": int(self.total_terms),
        }


# ---------------------------------------------------------------------------
# Pure IRP geometry
# ---------------------------------------------------------------------------


def raw_direction(n: int, seed: int) -> "np.ndarray":
    ensure_numpy()
    vals = []
    for j in range(n):
        r = math.sqrt(max(1e-300, -math.log(max(1e-300, u01(seed + 1009 * j + 17)))))
        th = 2.0 * math.pi * u01(seed + 1009 * j + 4099)
        vals.append(r * phase(th))
    v = np.asarray(vals, dtype=np.complex128)
    nv = max(1e-300, float(np.linalg.norm(v)))
    return v / nv


def atlas_start(n: int, trial: int, seed: int, radii: Sequence[float]) -> "np.ndarray":
    ensure_numpy()
    rr = [float(r) for r in radii if math.isfinite(float(r)) and float(r) > 0]
    if not rr:
        rr = [0.25, 0.5, 0.85, 1.25, 1.8, 2.6, 3.8]
    shell = rr[trial % len(rr)]
    direction = raw_direction(n, seed + 104729 * (trial + 1))
    # Mild anisotropic Thales/Riemann twist: coordinates get deterministic
    # phase and magnitude offsets, but the whole start remains shell-local.
    z = np.empty(n, dtype=np.complex128)
    for j in range(n):
        mag = shell * (0.72 + 0.56 * u01(seed + 7919 * (trial + 1) + 193 * (j + 1)))
        ph = phase(2.0 * math.pi * u01(seed + 15485863 * (trial + 1) + 271 * (j + 1)))
        z[j] = mag * ph * direction[j] * math.sqrt(max(1, n))
    return z


def complex_scale_palette(gains: Sequence[float], phases: Sequence[float], top: int) -> list[complex]:
    out: list[complex] = [1.0 + 0.0j]
    for g in gains:
        try:
            gg = float(g)
        except Exception:
            continue
        if not (math.isfinite(gg) and gg > 0):
            continue
        for ph in phases:
            try:
                pp = float(ph)
            except Exception:
                continue
            if not math.isfinite(pp):
                continue
            z = complex(gg * phase(pp))
            if all(abs(z - w) > 1e-14 for w in out):
                out.append(z)
            if len(out) >= max(1, int(top)):
                return out
    return out[: max(1, int(top))]


def log_energy(z: Sequence[complex]) -> float:
    zz = np.asarray(z, dtype=np.complex128)
    if zz.size == 0:
        return 0.0
    r = np.abs(zz)
    # scale-invariant-ish projective log energy with clipping for numerical safety
    vals = np.log1p(r) + np.log1p(1.0 / np.maximum(r, 1e-300))
    e = float(np.mean(vals))
    return e if math.isfinite(e) else float("inf")


def renormalize(system: DenseKostlanSystem, z: "np.ndarray", current_r: float, scales: Sequence[complex]) -> tuple["np.ndarray", float, dict[str, Any]]:
    ensure_numpy()
    zz = np.asarray(z, dtype=np.complex128)
    if len(scales) <= 1:
        return zz.copy(), float(current_r), {"renorm_changed": False, "renorm_evals": 0, "renorm_scale": cjson(1.0)}
    Z = np.asarray([complex(s) * zz for s in scales], dtype=np.complex128)
    R = system.residuals_batch(Z)
    # Residual dominates; log energy breaks ties toward less distorted charts.
    E = np.asarray([log_energy(row) for row in Z], dtype=float)
    scores = np.log1p(np.maximum(R, 0.0)) + 1e-12 * np.log1p(np.maximum(E, 0.0))
    scores[~np.isfinite(scores)] = np.inf
    if not np.any(np.isfinite(scores)):
        return zz.copy(), float(current_r), {"renorm_changed": False, "renorm_evals": int(len(scales)), "renorm_scale": cjson(1.0), "renorm_status": "no-finite-scale"}
    idx = int(np.nanargmin(scores))
    best_r = float(R[idx]) if math.isfinite(float(R[idx])) else float(current_r)
    changed = idx != 0 and best_r <= float(current_r) * (1.0 + 1e-14)
    if changed:
        return Z[idx].copy(), best_r, {"renorm_changed": True, "renorm_evals": int(len(scales)), "renorm_scale": cjson(scales[idx]), "renorm_residual": best_r, "renorm_log_energy": float(E[idx])}
    return zz.copy(), float(current_r), {"renorm_changed": False, "renorm_evals": int(len(scales)), "renorm_scale": cjson(1.0), "renorm_residual": float(current_r), "renorm_log_energy": float(E[0]) if len(E) else None}


def line_lambdas(raw: Optional[str], default_count: int = 6) -> list[float]:
    vals = parse_float_list(raw, [], positive=True)
    if not vals:
        vals = [1.0, 0.65, 0.42, 0.25, 0.15, 0.09, 0.055, 0.033]
    out: list[float] = []
    for v in vals:
        vv = min(max(float(v), 1e-15), 4.0)
        if all(abs(vv - u) > 1e-15 for u in out):
            out.append(vv)
        if len(out) >= max(1, int(default_count)):
            break
    return out or [1.0]


def endpoint_for_layer(z: "np.ndarray", prev_delta: Optional["np.ndarray"], seed: int, probe_scale: float, attempt: int) -> "np.ndarray":
    ensure_numpy()
    n = int(z.size)
    zn = max(1.0, float(np.linalg.norm(z)))
    if prev_delta is not None and np.all(np.isfinite(prev_delta)) and float(np.linalg.norm(prev_delta)) > 1e-14:
        direction = np.asarray(prev_delta, dtype=np.complex128)
        direction = direction / max(1e-300, float(np.linalg.norm(direction)))
    else:
        direction = raw_direction(n, seed + 1299709 * (attempt + 1))
    # Attempts rotate the same local scale; no residual-scored probe search.
    rot = phase(2.399963229728653 * attempt + 0.1732050807)
    jitter = raw_direction(n, seed + 104729 * (attempt + 11))
    direction = direction * rot + 0.18 * jitter
    direction = direction / max(1e-300, float(np.linalg.norm(direction)))
    h = float(probe_scale) * zn * (1.0 + 0.45 * attempt)
    b = z + h * direction
    # Avoid exact coordinate equality that can produce zero columns in unlucky cases.
    tiny = 1e-13 * zn
    for j in range(n):
        if abs(b[j] - z[j]) < tiny:
            b[j] += tiny * phase(0.37 + j + attempt)
    return b.astype(np.complex128)


@dataclasses.dataclass
class IRPResult:
    z: Any
    residual: float
    best_z: Any
    best_residual: float
    ok: bool
    status: str
    epochs: int
    layers_completed: int
    raw_steps: int
    renorm_steps: int
    line_evals: int
    seconds: float
    last_lambda: Optional[float]
    last_cond: Optional[float]
    last_delta_norm: Optional[float]
    collapsed: bool
    collapse_epoch: Optional[int]
    collapse_reason: Optional[str]
    capture_steps: int
    local_steps: int
    recoveries: int
    trace: list[float]


def _raw_pandrosion_layer(
    system: DenseKostlanSystem,
    z: "np.ndarray",
    r: float,
    prev_delta: Optional["np.ndarray"],
    lambdas: Sequence[float],
    probe_scale: float,
    probe_tries: int,
    trust: float,
    seed: int,
    ep: int,
    layer: int,
    local_mode: bool,
) -> tuple[Optional["np.ndarray"], float, Optional["np.ndarray"], dict[str, Any]]:
    """One raw derivative-free Pandrosion layer.

    This is intentionally still the 000/309 raw geometry: a finite telescopic
    slope matrix Q_G(a,b) and a solve Q_G delta = -G(a).  No Newton, Halley,
    Steffensen, Hessian, or analytic Jacobian is used.  The only difference in
    local_mode is that the probe/line grids are collapsed to a cheaper local
    variant once the orbit appears self-localized.
    """
    ensure_numpy()
    f = system.eval(z)
    r0 = finite_norm(f)
    if not math.isfinite(r0):
        return None, float("inf"), None, {"status": "nonfinite-residual"}

    best_candidate = None
    best_candidate_r = float("inf")
    best_delta = None
    best_meta: dict[str, Any] = {"status": "no-candidate"}

    attempts = max(1, int(probe_tries))
    if local_mode:
        attempts = 1
    for attempt in range(attempts):
        b = endpoint_for_layer(z, prev_delta, seed + 104729 * (ep + 1) + 7919 * (layer + 1), probe_scale, attempt)
        try:
            Q = system.slope_matrix(z, b)
            cond = float(np.linalg.cond(Q))
            delta = np.linalg.solve(Q, -f)
        except Exception as exc:
            best_meta = {"status": "solve-failed", "error": str(exc)[:120]}
            continue
        if not np.all(np.isfinite(delta)):
            best_meta = {"status": "nonfinite-delta"}
            continue
        zn = max(1.0, float(np.linalg.norm(z)))
        dn = float(np.linalg.norm(delta))
        if not math.isfinite(dn) or dn <= 0:
            best_meta = {"status": "bad-delta-norm"}
            continue

        # Geometric trust cap.  In local collapsed mode we make the cap tighter:
        # the orbit is assumed to be in a self-localized cell, so huge jumps are
        # usually harmful and expensive to line-search.
        cond_eff = cond if math.isfinite(cond) else 1e12
        cap = max(1e-14, float(trust) * zn / (1.0 + 0.02 * math.log1p(max(0.0, cond_eff))))
        if local_mode:
            cap = min(cap, 0.75 * zn)
        if dn > cap:
            delta = delta * (cap / dn)
            dn = cap

        # Local mode uses the leading part of the line grid.  This is the actual
        # collapse: after capture, the solver no longer pays for a large family of
        # candidate charts/line steps every layer.
        if local_mode:
            line = list(lambdas[: max(1, min(3, len(lambdas)))])
            if 1.0 not in line:
                line = [1.0] + line
        else:
            line = list(lambdas)
        Zcand = np.asarray([z + float(lam) * delta for lam in line], dtype=np.complex128)
        Rcand = system.residuals_batch(Zcand)
        if not np.any(np.isfinite(Rcand)):
            best_meta = {"status": "no-finite-line"}
            continue
        idx = int(np.nanargmin(Rcand))
        rr = float(Rcand[idx])
        if rr < best_candidate_r:
            best_candidate_r = rr
            best_candidate = Zcand[idx].copy()
            best_delta = best_candidate - z
            best_meta = {
                "status": "candidate",
                "lambda": float(line[idx]),
                "cond": cond,
                "delta_norm": float(np.linalg.norm(best_delta)),
                "raw_delta_norm": dn,
                "attempt": int(attempt),
                "line_evals": int(len(line)),
                "local_mode": bool(local_mode),
                "residual_before": float(r0),
                "residual_after": float(rr),
            }
    return best_candidate, best_candidate_r, best_delta, best_meta


def self_localizing_irp_corrector(
    system: DenseKostlanSystem,
    z0: Sequence[complex],
    epochs: int,
    layers: int,
    accept: float,
    tol: float,
    trial_timeout: float,
    scales: Sequence[complex],
    lambdas: Sequence[float],
    probe_scale: float,
    probe_tries: int,
    trust: float,
    seed: int,
    collapse: bool = True,
    collapse_residual: float = 1e-4,
    collapse_drop: float = 0.42,
    collapse_rel_step: float = 0.35,
    collapse_after: int = 2,
    local_probe_scale: float = 0.022,
    local_trust: float = 1.15,
    recoveries_max: int = 2,
    trace: bool = False,
) -> IRPResult:
    """Self-localizing IRP.

    Phase A (capture): full IRP cascade

        R -> P -> R -> P -> ...

    where R is a scalar geometric renormalization over a small complex scale
    palette and P is the raw finite-telescopic Pandrosion layer.

    Phase B (collapse): when the orbit looks local, R is frozen to the identity
    and P is run with a smaller probe and shorter line grid.  Failed local steps
    can temporarily recover to capture mode.  This implements the 310 idea:

        use geometric renormalization only until the fixed-point geometry becomes
        self-localizing.
    """
    ensure_numpy()
    t0 = now()
    deadline = t0 + float(trial_timeout) if trial_timeout and trial_timeout > 0 else None
    z = np.asarray(z0, dtype=np.complex128).copy()
    r = system.residual(z)
    best_z = z.copy()
    best_r = float(r)
    prev_delta: Optional[np.ndarray] = None
    ok = False
    status = "started"
    line_evals = 0
    renorm_steps = 0
    raw_steps = 0
    layers_completed = 0
    capture_steps = 0
    local_steps = 0
    recoveries = 0
    collapsed = False
    collapse_epoch: Optional[int] = None
    collapse_reason: Optional[str] = None
    good_locality_hits = 0
    last_lambda: Optional[float] = None
    last_cond: Optional[float] = None
    last_delta_norm: Optional[float] = None
    hist: list[float] = [float(r)] if trace else []

    if math.isfinite(r) and r <= max(float(accept), float(tol)):
        return IRPResult(z=z, residual=r, best_z=best_z, best_residual=best_r, ok=True, status="converged-initial", epochs=0, layers_completed=0, raw_steps=0, renorm_steps=0, line_evals=0, seconds=now()-t0, last_lambda=None, last_cond=None, last_delta_norm=None, collapsed=False, collapse_epoch=None, collapse_reason=None, capture_steps=0, local_steps=0, recoveries=0, trace=hist)

    max_epochs = max(1, int(epochs))
    for ep in range(max_epochs):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break
        improved_epoch = False
        layers_this_epoch = 1 if collapsed else max(1, int(layers))
        for layer in range(layers_this_epoch):
            if deadline is not None and now() > deadline:
                status = "timeout"
                break

            # R: active only during capture.  In collapsed mode the chart is
            # frozen; this is the cost-saving locality collapse.
            if not collapsed:
                z, r, rmeta = renormalize(system, z, r, scales)
                renorm_steps += int(bool(rmeta.get("renorm_changed", False)))
                if math.isfinite(r) and r < best_r:
                    best_r = r
                    best_z = z.copy()
                    improved_epoch = True
                if trace:
                    hist.append(float(r))
                if math.isfinite(r) and r <= max(float(accept), float(tol)):
                    ok = True
                    status = "converged"
                    break

            r_before = float(r)
            z_before = z.copy()
            local_mode = bool(collapsed)
            raw_probe = float(local_probe_scale if local_mode else probe_scale)
            raw_trust = float(local_trust if local_mode else trust)
            cand, cand_r, delta, meta = _raw_pandrosion_layer(
                system=system,
                z=z,
                r=r,
                prev_delta=prev_delta,
                lambdas=lambdas,
                probe_scale=raw_probe,
                probe_tries=probe_tries,
                trust=raw_trust,
                seed=seed,
                ep=ep,
                layer=layer,
                local_mode=local_mode,
            )
            line_evals += int(meta.get("line_evals", 0))
            if cand is None:
                if collapsed and recoveries < max(0, int(recoveries_max)):
                    # Locality detector was too optimistic.  Reopen the chart.
                    recoveries += 1
                    collapsed = False
                    collapse_reason = (collapse_reason or "") + "; recovered-from-local-failure"
                    status = "recovered-to-capture"
                    break
                status = str(meta.get("status", "raw-solve-failed"))
                break

            # Accept strict descent, or any improvement over the global best.
            if cand_r <= r_before * (1.0 - 1e-14) or cand_r < best_r:
                z = cand
                r = float(cand_r)
                raw_steps += 1
                layers_completed += 1
                if collapsed:
                    local_steps += 1
                else:
                    capture_steps += 1
                prev_delta = np.asarray(delta, dtype=np.complex128) if delta is not None else (z - z_before)
                last_lambda = meta.get("lambda")
                last_cond = meta.get("cond")
                last_delta_norm = meta.get("delta_norm")
                if r < best_r:
                    best_r = r
                    best_z = z.copy()
                improved_epoch = True
                if trace:
                    hist.append(float(r))
                if r <= max(float(accept), float(tol)):
                    ok = True
                    status = "converged"
                    break

                # Self-locality detector.  It uses only observed geometric data:
                # residual drop and relative step size.  No derivative order test
                # and no Halley/Steffensen condition is involved.
                step_norm = finite_norm(prev_delta)
                rel_step = step_norm / max(1.0, finite_norm(z_before))
                drop = r / max(r_before, 1e-300)
                hit = False
                reason = None
                if bool(collapse):
                    if r <= float(collapse_residual):
                        hit = True
                        reason = "residual-threshold"
                    elif drop <= float(collapse_drop) and rel_step <= float(collapse_rel_step):
                        good_locality_hits += 1
                        if good_locality_hits >= max(1, int(collapse_after)):
                            hit = True
                            reason = "observed-local-contraction"
                    else:
                        good_locality_hits = 0
                if hit and not collapsed:
                    collapsed = True
                    collapse_epoch = ep
                    collapse_reason = reason
                    # Shrink the next probe to the recent geometry.
                    prev_delta = prev_delta / max(1.0, float(np.linalg.norm(prev_delta))) * min(finite_norm(prev_delta), float(local_probe_scale) * max(1.0, finite_norm(z)))
            else:
                if collapsed and recoveries < max(0, int(recoveries_max)):
                    recoveries += 1
                    collapsed = False
                    good_locality_hits = 0
                    status = "recovered-to-capture-no-descent"
                    break
                status = "no-raw-descent"
                break

        if ok or status in {"timeout", "nonfinite-residual", "raw-solve-failed", "solve-failed", "nonfinite-delta", "bad-delta-norm", "no-finite-line", "no-raw-descent"}:
            break
        if not improved_epoch and status not in {"recovered-to-capture", "recovered-to-capture-no-descent"}:
            status = "stagnated"
            break
        status = "local-collapsed" if collapsed else "capture-running"

    final_r = system.residual(z)
    if math.isfinite(final_r) and final_r < best_r:
        best_r = final_r
        best_z = z.copy()
    if math.isfinite(best_r) and best_r <= max(float(accept), float(tol)):
        ok = True
        status = "converged"
        z = best_z.copy()
        final_r = best_r
    return IRPResult(
        z=z, residual=float(final_r), best_z=best_z, best_residual=float(best_r), ok=ok,
        status=status, epochs=ep + 1 if 'ep' in locals() else 0,
        layers_completed=int(layers_completed), raw_steps=int(raw_steps), renorm_steps=int(renorm_steps),
        line_evals=int(line_evals), seconds=now() - t0, last_lambda=last_lambda,
        last_cond=last_cond, last_delta_norm=last_delta_norm, collapsed=bool(collapsed),
        collapse_epoch=collapse_epoch, collapse_reason=collapse_reason,
        capture_steps=int(capture_steps), local_steps=int(local_steps), recoveries=int(recoveries), trace=hist,
    )

# ---------------------------------------------------------------------------
# Root handling and case runner
# ---------------------------------------------------------------------------


def cluster_index(roots: list[dict[str, Any]], z: Sequence[complex], sep: float) -> Optional[int]:
    zz = np.asarray(z, dtype=np.complex128)
    zn = max(1.0, float(np.linalg.norm(zz)))
    for k, root in enumerate(roots):
        rr = np.asarray(root["z_complex"], dtype=np.complex128)
        rn = max(1.0, float(np.linalg.norm(rr)))
        if float(np.linalg.norm(zz - rr)) <= float(sep) * max(zn, rn):
            return k
    return None


def realness(z: Sequence[complex]) -> float:
    zz = np.asarray(z, dtype=np.complex128)
    return float(np.linalg.norm(zz.imag) / max(1.0, np.linalg.norm(zz)))


def score_root(residual: float, realv: float, cond: Optional[float]) -> float:
    c = 0.0 if cond is None or not math.isfinite(float(cond)) else 1e-10 * math.log1p(float(cond))
    return float(math.log10(max(residual, 1e-300)) + 0.01 * realv + c)


def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    ensure_numpy()
    n, d = parse_case(case_raw)
    t_case = now()
    system = DenseKostlanSystem.make(n=n, d=d, seed_index=int(args.seed_index), equation_normalize=bool(args.equation_normalize))

    radii = parse_float_list(args.rays, [0.18, 0.32, 0.55, 0.85, 1.25, 1.8, 2.6, 3.8, 5.4], positive=True)
    gains = parse_float_list(args.renorm_gains, [1.0, 0.76, 1.32, 0.58, 1.72, 0.43, 2.30], positive=True)
    phases = parse_float_list(args.renorm_phases, [0.0, 0.035, -0.035, 0.085, -0.085])
    scales = complex_scale_palette(gains, phases, int(args.renorm_top))
    start_gains = parse_float_list(args.startopt_gains, [1.0, 0.62, 1.62, 0.38, 2.62], positive=True)
    start_phases = parse_float_list(args.startopt_phases, [0.0, 0.11, -0.11, 0.23, -0.23])
    start_scales = complex_scale_palette(start_gains, start_phases, int(args.startopt_top))
    lambdas = line_lambdas(args.line_grid, int(args.line_search))

    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    failures = 0
    duplicates = 0
    t_extract = now()

    for trial in range(max(1, int(args.pool))):
        if len(roots) >= int(args.count):
            break
        z0 = atlas_start(n, trial, system.seed + 17 * int(args.seed_index), radii)
        r0 = system.residual(z0)
        # Start optimization is a one-shot geometric renormalization of the atlas point.
        z0, r0, smeta = renormalize(system, z0, r0, start_scales)
        loc = self_localizing_irp_corrector(
            system=system,
            z0=z0,
            epochs=int(args.epochs),
            layers=int(args.irp_layers),
            accept=float(args.accept),
            tol=float(args.tol),
            trial_timeout=float(args.trial_timeout),
            scales=scales,
            lambdas=lambdas,
            probe_scale=float(args.probe_scale),
            probe_tries=int(args.probe_tries),
            trust=float(args.trust),
            seed=system.seed + 131071 * (trial + 1),
            collapse=bool(args.collapse),
            collapse_residual=float(args.collapse_residual),
            collapse_drop=float(args.collapse_drop),
            collapse_rel_step=float(args.collapse_rel_step),
            collapse_after=int(args.collapse_after),
            local_probe_scale=float(args.local_probe_scale),
            local_trust=float(args.local_trust),
            recoveries_max=int(args.recoveries),
            trace=bool(args.trace_trials),
        )
        z = np.asarray(loc.best_z if loc.best_residual <= loc.residual else loc.z, dtype=np.complex128)
        r = float(min(loc.best_residual, loc.residual))
        accepted = bool(loc.ok and math.isfinite(r) and r <= float(args.accept))
        rec = {
            "trial": int(trial),
            "status": loc.status,
            "accepted": accepted,
            "initial_residual": float(r0),
            "residual": float(r),
            "epochs": int(loc.epochs),
            "layers_completed": int(loc.layers_completed),
            "raw_steps": int(loc.raw_steps),
            "renorm_steps": int(loc.renorm_steps),
            "line_evals": int(loc.line_evals),
            "seconds": float(loc.seconds),
            "last_lambda": loc.last_lambda,
            "last_cond": loc.last_cond,
            "last_delta_norm": loc.last_delta_norm,
            "collapsed": bool(loc.collapsed),
            "collapse_epoch": loc.collapse_epoch,
            "collapse_reason": loc.collapse_reason,
            "capture_steps": int(loc.capture_steps),
            "local_steps": int(loc.local_steps),
            "recoveries": int(loc.recoveries),
            "start_renorm_changed": bool(smeta.get("renorm_changed", False)),
        }
        if bool(args.trace_trials):
            rec["z"] = root_to_json(z)
            rec["trace"] = loc.trace[:200]
        if not accepted:
            failures += 1
            if len(trials) < int(args.keep_trials):
                trials.append(rec)
            continue
        dup = cluster_index(roots, z, float(args.cluster_sep))
        if dup is not None:
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = int(dup)
            if len(trials) < int(args.keep_trials):
                trials.append(rec)
            continue
        rid = len(roots)
        cond = loc.last_cond
        rv = realness(z)
        roots.append({
            "id": int(rid),
            "source": "310-self-localizing-iterated-renormalized-pandrosion",
            "trial": int(trial),
            "z_complex": z.copy(),
            "residual": float(r),
            "realness": float(rv),
            "cond": cond,
            "score": score_root(float(r), rv, cond),
            "epochs": int(loc.epochs),
            "layers_completed": int(loc.layers_completed),
            "raw_steps": int(loc.raw_steps),
            "renorm_steps": int(loc.renorm_steps),
            "line_evals": int(loc.line_evals),
            "seconds": float(loc.seconds),
            "last_lambda": loc.last_lambda,
            "last_delta_norm": loc.last_delta_norm,
            "collapsed": bool(loc.collapsed),
            "collapse_epoch": loc.collapse_epoch,
            "collapse_reason": loc.collapse_reason,
            "capture_steps": int(loc.capture_steps),
            "local_steps": int(loc.local_steps),
            "recoveries": int(loc.recoveries),
            "halley_enabled": False,
            "steffensen_enabled": False,
            "newton_enabled": False,
        })
        rec["status"] = "new-root"
        rec["root_id"] = rid
        if len(trials) < int(args.keep_trials):
            trials.append(rec)

    encoded_roots = []
    for root in sorted(roots, key=lambda q: (float(q.get("score", float("inf"))), int(q.get("id", 0)))):
        rr = dict(root)
        zc = rr.pop("z_complex")
        rr["z"] = root_to_json(zc)
        encoded_roots.append(rr)

    result = {
        "script": "310_pandrosion_self_localizing_irp_engine.py",
        "autonomous": True,
        "dependencies": {"python_scripts": [], "numpy": bool(np is not None)},
        "mode": "310-self-localizing-iterated-renormalized-pandrosion-irp",
        "flow_formula": "atlas start -> geometric start renormalization -> IRP cascade (R∘P)^k: scalar homothety renormalization R, exact finite telescopic Pandrosion slope Q_G(a,b), raw solve Q delta=-G(a), small residual-gated line step; no Halley, no Steffensen, no Newton/Jacobian/Hessian fallback",
        "case": f"{n},{d}",
        "family": "ks",
        "seed_index": int(args.seed_index),
        "seed": int(system.seed),
        "n": int(n),
        "degree": int(d),
        "terms_per_poly": int(system.terms_per_poly),
        "terms": int(system.total_terms),
        "bezout": int(system.bezout),
        "equation_normalize": bool(args.equation_normalize),
        "parameters": {
            "count": int(args.count),
            "pool": int(args.pool),
            "epochs": int(args.epochs),
            "irp_layers": int(args.irp_layers),
            "accept": float(args.accept),
            "tol": float(args.tol),
            "cluster_sep": float(args.cluster_sep),
            "probe_scale": float(args.probe_scale),
            "probe_tries": int(args.probe_tries),
            "trust": float(args.trust),
            "collapse": bool(args.collapse),
            "collapse_residual": float(args.collapse_residual),
            "collapse_drop": float(args.collapse_drop),
            "collapse_rel_step": float(args.collapse_rel_step),
            "collapse_after": int(args.collapse_after),
            "local_probe_scale": float(args.local_probe_scale),
            "local_trust": float(args.local_trust),
            "recoveries": int(args.recoveries),
            "line_search": int(args.line_search),
            "line_grid": lambdas,
            "renorm_gains": gains,
            "renorm_phases": phases,
            "renorm_top": int(args.renorm_top),
            "startopt_top": int(args.startopt_top),
            "rays": radii,
        },
        "roots": encoded_roots,
        "trials": trials,
        "summary": {
            "requested_roots": int(args.count),
            "unique_roots": int(len(roots)),
            "success": bool(len(roots) >= int(args.count)),
            "trials_used": int(min(int(args.pool), len(trials) + failures + duplicates + len(roots))),
            "duplicates": int(duplicates),
            "failures": int(failures),
            "generation_seconds": float(system.generation_seconds),
            "extract_seconds": float(now() - t_extract),
            "total_seconds": float(now() - t_case),
            "eval_stats": system.stats(),
        },
    }
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="310 self-localizing NumPy engine for Iterated Renormalized Pandrosion (IRP).")
    p.add_argument("--cases", default="2,4", help="case n,d; multiple cases separated by ';'")
    p.add_argument("--seed-index", type=int, default=0)
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--no-equation-normalize", dest="equation_normalize", action="store_false")
    p.add_argument("--count", type=int, default=8)
    p.add_argument("--pool", type=int, default=2048)
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--irp-layers", type=int, default=2)
    p.add_argument("--accept", type=float, default=1e-8)
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--cluster-sep", type=float, default=1e-8)
    p.add_argument("--trial-timeout", type=float, default=0.0)
    p.add_argument("--probe-scale", type=float, default=0.045)
    p.add_argument("--probe-tries", type=int, default=2)
    p.add_argument("--trust", type=float, default=2.8)
    p.add_argument("--collapse", action="store_true", default=True, help="enable self-localizing collapse after capture")
    p.add_argument("--no-collapse", dest="collapse", action="store_false", help="disable collapse and run 309-style IRP")
    p.add_argument("--collapse-residual", type=float, default=1e-4, help="absolute residual threshold for switching to local mode")
    p.add_argument("--collapse-drop", type=float, default=0.42, help="observed residual drop ratio needed for locality hit")
    p.add_argument("--collapse-rel-step", type=float, default=0.35, help="relative step size threshold for locality hit")
    p.add_argument("--collapse-after", type=int, default=2, help="number of consecutive locality hits before collapse")
    p.add_argument("--local-probe-scale", type=float, default=0.022, help="probe scale after IRP collapse")
    p.add_argument("--local-trust", type=float, default=1.15, help="trust cap after IRP collapse")
    p.add_argument("--recoveries", type=int, default=2, help="number of allowed returns from local mode to capture mode")
    p.add_argument("--line-search", type=int, default=6)
    p.add_argument("--line-grid", default="1,0.65,0.42,0.25,0.15,0.09,0.055")
    p.add_argument("--renorm-gains", default="1,0.76,1.32,0.58,1.72,0.43,2.3")
    p.add_argument("--renorm-phases", default="0,0.035,-0.035,0.085,-0.085")
    p.add_argument("--renorm-top", type=int, default=12)
    p.add_argument("--startopt-gains", default="1,0.62,1.62,0.38,2.62")
    p.add_argument("--startopt-phases", default="0,0.11,-0.11,0.23,-0.23")
    p.add_argument("--startopt-top", type=int, default=10)
    p.add_argument("--rays", default="0.18,0.32,0.55,0.85,1.25,1.8,2.6,3.8,5.4")
    p.add_argument("--out", default=None)
    p.add_argument("--outdir", default="/mnt/data/310_self_localizing_irp_out")
    p.add_argument("--keep-trials", type=int, default=160)
    p.add_argument("--trace-trials", action="store_true", default=False)
    p.add_argument("--self-test", action="store_true", help="run a small ks(2,2) smoke test")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(args.self_test):
        args.cases = "2,2"
        args.count = 2
        args.pool = min(int(args.pool), 256)
        args.epochs = min(int(args.epochs), 16)
        args.accept = min(float(args.accept), 1e-8)
        args.keep_trials = min(int(args.keep_trials), 40)
        args.out = args.out or "/mnt/data/310_self_localizing_irp_out/self_test_310.json"
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    outputs = [run_case(args, c) for c in cases]
    if len(outputs) == 1:
        final = outputs[0]
    else:
        final = {"script": "310_pandrosion_self_localizing_irp_engine.py", "autonomous": True, "mode": "310-self-localizing-irp-multicase", "cases": outputs}
    if args.out:
        out = Path(args.out)
    else:
        first = cases[0].replace(",", "x") if cases else "case"
        out = Path(args.outdir) / f"310_self_localizing_irp_{first}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2), encoding="utf-8")

    print("=" * 116, flush=True)
    print("310 SELF-LOCALIZING ITERATED RENORMALIZED PANDROSION NumPy engine", flush=True)
    print("IRP cascade only: geometric renormalization -> raw exact telescopic Pandrosion slope solve; no Halley/Steffensen/Newton fallback.", flush=True)
    print("=" * 116, flush=True)
    for r in outputs:
        s = r["summary"]
        print(f"case=ks({r['n']},{r['degree']}), seed={r['seed']}, terms={r['terms']}, Bezout={r['bezout']}", flush=True)
        print(f"roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} failures={s['failures']} duplicates={s['duplicates']}", flush=True)
        print(f"seconds: generation={s['generation_seconds']:.3f}, extract={s['extract_seconds']:.3f}, total={s['total_seconds']:.3f}", flush=True)
        if r.get("roots"):
            best = r["roots"][0]
            print(f"best_root: residual={float(best.get('residual', float('inf'))):.3e}, trial={best.get('trial')}, raw_steps={best.get('raw_steps')}, layers={best.get('layers_completed')}, collapsed={best.get('collapsed')}", flush=True)
    print(f"out={out}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        raise SystemExit(130)
