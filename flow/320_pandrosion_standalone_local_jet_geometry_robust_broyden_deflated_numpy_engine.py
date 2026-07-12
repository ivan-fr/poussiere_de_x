"""
320_pandrosion_standalone_local_jet_geometry_robust_broyden_deflated_numpy_engine.py

Standalone NumPy engine that realises pandrosion as a FULLY (jet-)GEOMETRIC method.

NEW IN 320 (over 319) -- same Pandrosion architecture, hardened invariants:

  [320-0] ACCEPTANCE-SAFE DEFLATION.  The deflated line search selects only
      candidates that passed its admissibility mask, and accumulates repulsion
      in log-space so harvesting many roots cannot overflow the merit function.

  [320-1] ANCHORED BROYDEN STATE.  A recycled Jacobian carries the point/value
      pair at which it is anchored, separately from the current residual reuse
      state.  Cross-call secant updates are therefore real Broyden updates, not
      silent reuse of a Jacobian sampled at the previous iterate.  Adaptive LM
      damping is persisted through the same state.

  [320-2] DEADLINE PROPAGATION.  Trial deadlines cover direct correction,
      chart scoring and every IRP chart correction rather than only the direct
      phase.

  [320-3] HONEST CONDITIONING.  Every finite accepted candidate gets a final
      finite-difference SVD estimate even when polishing is disabled.  A zero
      singular value is reported as singular/undefined conditioning, never as
      cond=0.  The result is explicitly called an estimate, not a certificate.

  [320-4] AUDITABLE ORACLE ACCOUNTING.  Jet and inverse-jet metadata report
      actual oracle evaluations (2n+1 for a paired local jet; 2n for a paired
      correction when F(y) is already known), distinct from stencil nodes.

  [320-5] RADIUS-SAFE LRU JETS.  The local-jet cache key includes the requested
      radius, preventing a jet sampled at one scale from being returned for
      another.

The 319 improvements retained below are Broyden-recycled jets, deflated merit,
adaptive LM, parabolic refinement, residual-vector reuse, vectorised swarm,
lazy IRP chart rescue, ghost escapes, singular-direction seeds and final SVD
conditioning estimates.

Historical 319 feature notes:

  [320-A] BROYDEN-RECYCLED JETS.  Rebuilding the 2n-sample paired jet at EVERY
      correction epoch is the dominant oracle cost.  320 keeps the last Jacobian
      alive and, on epochs where the previous step descended, replaces the full
      rebuild by a rank-1 Broyden secant update A += ((df - A dy) dy^H)/(dy^H dy)
      -- ZERO oracle samples for the Jacobian on those epochs.  The true paired
      jet is rebuilt every --jet-refresh epochs, after any stall, and after any
      chart/ghost jump, so superlinear local convergence is preserved while the
      per-trial oracle budget drops sharply (often 2-3x on the direct phase).

  [320-B] DEFLATED LINE SEARCH (root repulsion).  318 could only detect a
      duplicate after drifting into a known root's neighbourhood ([318-E]).
      320 additionally DEFLATES the merit function of the line search on the
      base chart: candidate residuals are multiplied by
      prod_k (1 + alpha * s_k / ||y - y_k||), alpha = --deflation-alpha, where
      y_k are the accepted roots.  Convergence/acceptance still uses the
      faithful raw residual, but the descent dynamics is now REPELLED from
      basins already harvested, so trials are spent on new roots instead of
      rediscovering old ones.  (Classic deflation a la Brown-Gearhart /
      Farrell-Birkisson-Funke, applied here purely on the merit function so the
      oracle and the Newton direction stay untouched.)

  [320-C] ADAPTIVE LEVENBERG-MARQUARDT.  The fixed --lm-damping of 318 is now a
      live trust parameter: mu shrinks (x0.4) after every strongly descending
      step and grows (x5, then retry with a FRESH full jet) after a rejected
      one, in the classic LM nu-schedule.  Far from a root the step is safely
      damped; near a root mu -> mu_min and the pure inverse-jet order is
      recovered.  Enabled by default (--no-adaptive-lm restores 318 behaviour).

  [320-D] PARABOLIC LINE REFINEMENT.  After the geometric lambda grid the line
      search fits a parabola through the best bracket (in lambda, r^2) and
      spends exactly ONE extra oracle sample on the predicted minimiser,
      routinely gaining an extra digit per epoch for free.

  [320-E] ORACLE-SAMPLE REUSE.  The line search now returns the full residual
      VECTORS, and the winning candidate's F(y) is carried into the next epoch
      instead of being re-evaluated; the stall-retry path reuses it as well.
      One sample saved per epoch per corrector, exactly zero change in maths.

  [320-G] SWARM ORACLE REUSE + FAIR IRP TRIGGER.  The batched swarm prefilter
      now carries the residual VECTORS of its line search into the next
      iteration as the ready-made F0 of the Jacobian build, saving one full
      batch evaluation per iteration.  (A batched secant variant of the swarm
      Jacobians was tried and REJECTED: it measurably degraded basin
      separation.)  The 312 lazy-IRP trigger judges the direct-corrector drop
      PER inner epoch (geometric compounding), and --lazy-direct-epochs
      defaults to 2, so multi-epoch secant runs face the same per-step bar as
      318's single full-jet steps instead of spuriously triggering the
      expensive chart rescue.

  [320-F] FINITE-DIFFERENCE SVD CONDITIONING ESTIMATE.  Every accepted root is finished with an
      SVD of its final local jet: sigma_min and cond(J) are reported per root
      (318 reported cond=None), feed score_root, and flag near-multiple roots.

Inherited from 318: paired central-difference jets with per-coordinate radii and
truncated-SVD fallback [318-A]; vectorised swarm prefilter [318-B]; ghost-trap
escapes [318-C]; guaranteed chart-kind diversity in the IRP shortlist [318-D];
early duplicate abort [318-E]; singular-direction seeds + bounded LRU jet cache
[318-F].

History of this file: an earlier attempt built a single fixed-size GLOBAL projective
(Fubini-Study) sketch of the system and ran everything on it.  An empirical study
(see 317_probe_*.py) showed that representation does NOT converge to the system's
zero set at any anchor budget -- O(1), O(log n) or O(k) -- because the projective
kernel compactifies C and conflates the regions where the roots live.  A Euclidean
*local* approximation, by contrast, reaches machine precision with ~16 samples.

Conclusion adopted here: the geometry that faithfully *represents* a system is not a
global projective sketch but its FIELD OF LOCAL JETS.  The base system is treated
as a black-box geometric oracle and, around the current point, sampled on a small
paired stencil (2n+1 points) to build a local affine/curvature jet.  The whole
pandrosion stack (universal Mobius atlas starts, StartOpt, 306 hypercube inverse-jet,
312 lazy IRP chart rescue, jet-Newton polish) navigates on these jets.  Per correction
step the oracle is queried O(n) times on rebuild epochs and O(1) times on Broyden
epochs -- independent of the global monomial enumeration -- and the local-jet
layer adds no global precompute (the chosen base oracle may still have its own
construction cost).

Because the jets are sampled from the TRUE system, the residual ||F|| is faithful, so
geometric acceptance and exact acceptance coincide: this is "100% geometric" in the
jet sense while still locating genuine roots.  A per-step oracle-sample budget is
reported so you can verify the cost.

Systems supported as the geometric oracle:
  - exact user polynomial systems (--system-source poly --polys "...");
  - KS/Kostlan systems (dense / lazy-feature / geometry-kernel backends).

Dependencies: Python stdlib + NumPy only.  No local flow imports.
"""
from __future__ import annotations

import argparse
import ast
import dataclasses
import json
import math
import time
from collections import OrderedDict
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def now() -> float:
    return time.time()


def cjson(z: complex) -> list[float]:
    w = complex(z)
    return [float(w.real), float(w.imag)]


def root_to_json(z: Sequence[complex]) -> list[list[float]]:
    return [cjson(v) for v in z]


def json_safe(value: Any) -> Any:
    """Convert diagnostics to strict RFC-compatible JSON values."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, complex):
        return cjson(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def parse_case(raw: str) -> tuple[int, int]:
    s = str(raw).strip().lower().replace("x", ",").replace(":", ",")
    p = [q.strip() for q in s.split(",") if q.strip()]
    if len(p) != 2:
        raise ValueError(f"case must be n,d, got {raw!r}")
    n, d = int(p[0]), int(p[1])
    if n <= 0 or d <= 0:
        raise ValueError("n,d must be positive")
    return n, d


def parse_float_list(raw: Optional[str], default: Sequence[float], positive: bool = False) -> list[float]:
    if raw is None or str(raw).strip() == "":
        return list(default)
    out: list[float] = []
    for part in str(raw).replace(";", ",").split(","):
        try:
            x = float(part.strip())
        except Exception:
            continue
        if math.isfinite(x) and ((not positive) or x > 0):
            out.append(x)
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


def stable_seed(n: int, d: int, seed_index: int = 0, salt: int = 0) -> int:
    return int(splitmix64(0x50414E44524F5349 + 1000003 * n + 9176 * d + 97 * seed_index + salt) & 0x7FFFFFFF)


def exact_kostlan_terms(n: int, d: int) -> int:
    return int(math.comb(int(n) + int(d), int(d)))


def auto_lazy_feature_count(n: int, d: int, cap: int) -> int:
    base = int(32 * math.sqrt(max(1, (int(n) + 1) * (int(d) + 1))))
    return int(max(int(n) + 8, min(max(1, int(cap)), max(512, base))))


def auto_geometry_anchor_count(n: int, d: int, cap: int) -> int:
    base = int(24 * math.sqrt(max(1, (int(n) + 1) * (int(d) + 1))))
    return int(max(int(n) + 16, min(max(1, int(cap)), max(512, base))))


def finite_norm(v: Any) -> float:
    try:
        r = float(np.linalg.norm(v))
        return r if math.isfinite(r) else float("inf")
    except Exception:
        return float("inf")


def log_energy(y: Sequence[complex]) -> float:
    yy = np.asarray(y, dtype=np.complex128)
    if yy.size == 0:
        return 0.0
    return float(np.mean(np.abs(np.log(np.maximum(np.abs(yy), 1e-300)))))


def realness(z: Sequence[complex]) -> float:
    zz = np.asarray(z, dtype=np.complex128)
    den = max(1e-300, float(np.linalg.norm(zz)))
    return float(np.linalg.norm(zz.imag) / den)


def score_root(residual: float, realness_value: float, cond: Optional[float]) -> float:
    c = float(cond) if cond is not None and math.isfinite(float(cond)) else 1e300
    return float(math.log1p(max(0.0, residual)) + 0.01 * math.log1p(max(0.0, realness_value)) + 0.001 * math.log1p(max(0.0, c)))


def cluster_index(roots: list[dict[str, Any]], z: Any, sep: float) -> Optional[int]:
    zz = np.asarray(z, dtype=np.complex128)
    for i, r in enumerate(roots):
        if float(np.linalg.norm(zz - np.asarray(r["z_complex"], dtype=np.complex128))) <= float(sep):
            return i
    return None


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
    pts = []
    for part in text.replace("|", ";").split(";"):
        if not part.strip():
            continue
        coords = [p.strip() for p in part.split(",") if p.strip()]
        if len(coords) != int(n):
            raise ValueError(f"start point {part!r} has {len(coords)} coord(s), expected {n}")
        pts.append(np.asarray([parse_complex_token(p) for p in coords], dtype=np.complex128))
    return pts


def safe_norms(F: Any) -> Any:
    R = np.linalg.norm(np.asarray(F, dtype=np.complex128), axis=1)
    R = np.asarray(R, dtype=float)
    R[~np.isfinite(R)] = np.inf
    return R


def finite_slope_matrix(system: Any, a: Sequence[complex], b: Sequence[complex]) -> Any:
    t0 = now()
    aa = np.asarray(a, dtype=np.complex128)
    bb = np.asarray(b, dtype=np.complex128)
    cur = aa.copy()
    f_prev = system.eval(cur)
    Q = np.zeros((int(system.n), int(system.n)), dtype=np.complex128)
    for j in range(int(system.n)):
        old = cur[j]
        cur[j] = bb[j]
        f_next = system.eval(cur)
        dz = bb[j] - old
        if abs(dz) > 1e-300:
            Q[:, j] = (f_next - f_prev) / dz
        else:
            h = 1e-6 * max(1.0, abs(old))
            plus, minus = cur.copy(), cur.copy()
            plus[j] = old + h
            minus[j] = old - h
            Q[:, j] = (system.eval(plus) - system.eval(minus)) / (2.0 * h)
        f_prev = f_next
    system.slope_count = int(getattr(system, "slope_count", 0)) + 1
    system.seconds_slope = float(getattr(system, "seconds_slope", 0.0)) + (now() - t0)
    return Q


# ---------------------------------------------------------------------------
# KS systems (used only as a black-box geometric oracle)
# ---------------------------------------------------------------------------

def compositions_leq(d: int, n: int) -> Any:
    out: list[tuple[int, ...]] = []

    def rec(pos: int, rem: int, cur: list[int]) -> None:
        if pos == n - 1:
            for k in range(rem + 1):
                out.append(tuple(cur + [k]))
            return
        for k in range(rem + 1):
            cur.append(k)
            rec(pos + 1, rem - k, cur)
            cur.pop()

    rec(0, int(d), [])
    return np.asarray(out, dtype=np.int16 if d < 32767 else np.int32)


def multinomial_kostlan_weights(exps: Any, d: int) -> Any:
    totals = np.sum(exps, axis=1).astype(np.int64)
    logfac = np.zeros(int(d) + 1, dtype=np.float64)
    for k in range(1, int(d) + 1):
        logfac[k] = logfac[k - 1] + math.log(k)
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
    equation_normalize: bool = True
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True) -> "DenseKostlanSystem":
        t0 = now()
        exps = compositions_leq(d, n)
        weights = multinomial_kostlan_weights(exps, d)
        seed = stable_seed(n, d, seed_index)
        rng = np.random.default_rng(seed)
        coeff = (rng.standard_normal((n, exps.shape[0])) + 1j * rng.standard_normal((n, exps.shape[0]))) / math.sqrt(2.0)
        coeff = coeff * weights[None, :]
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None]
        obj = cls(int(n), int(d), seed, exps, coeff.astype(np.complex128), weights, bool(equation_normalize))
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

    def monomials_batch(self, Z: Any) -> Any:
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        B, M = int(ZZ.shape[0]), self.terms_per_poly
        mon = np.ones((B, M), dtype=np.complex128)
        for j in range(self.n):
            p = np.empty((B, self.d + 1), dtype=np.complex128)
            p[:, 0] = 1.0
            if self.d > 0:
                p[:, 1] = ZZ[:, j]
                for k in range(2, self.d + 1):
                    p[:, k] = p[:, k - 1] * ZZ[:, j]
            mon *= p[:, self.exps[:, j]]
        return mon

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = (self.monomials_batch(np.asarray(z, dtype=np.complex128))[0] @ self.coeff.T)
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        out = self.monomials_batch(ZZ) @ self.coeff.T
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms}


@dataclasses.dataclass
class LazyFeatureKostlanSystem:
    n: int
    d: int
    seed: int
    feature_exps: Any
    feature_exps_t: Any
    feature_log_scales: Any
    coeff: Any
    lazy_features: int
    projective_normalize: bool = True
    dynamic_normalize: bool = True
    equation_normalize: bool = True
    eval_block: int = 128
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True, lazy_features: int = 0, lazy_feature_cap: int = 8192, projective_normalize: bool = True, dynamic_normalize: bool = True, eval_block: int = 128) -> "LazyFeatureKostlanSystem":
        t0 = now()
        n, d = int(n), int(d)
        m = max(n + 1, int(lazy_features) if int(lazy_features) > 0 else auto_lazy_feature_count(n, d, int(lazy_feature_cap)))
        seed = stable_seed(n, d, seed_index, salt=0x314314)
        rng = np.random.default_rng(seed)
        exps = np.zeros((m, n), dtype=np.int16 if d < 32767 else np.int32)
        degrees = np.zeros(m, dtype=np.int64)
        idx = 1
        for j in range(n):
            if idx >= m:
                break
            exps[idx, j] = 1 if d >= 1 else 0
            degrees[idx] = 1 if d >= 1 else 0
            idx += 1
        probs = np.full(n, 1.0 / max(1, n))
        while idx < m:
            q = (idx - (n + 1) + 0.5) / max(1, m - (n + 1))
            k = int(min(d, max(0, math.floor(q * (d + 1)))))
            if d > 0 and idx % 7 == 0:
                k = int(rng.integers(0, d + 1))
            if k:
                exps[idx, :] = rng.multinomial(k, probs).astype(exps.dtype)
            degrees[idx] = k
            idx += 1
        log_scales = np.empty(m, dtype=np.float64)
        log_m, log_deg, log_n = math.log(max(1, m)), math.log(max(1, d + 1)), math.log(max(1, n))
        for i, k in enumerate(degrees):
            lc = math.lgamma(d + 1) - math.lgamma(int(k) + 1) - math.lgamma(d - int(k) + 1)
            log_scales[i] = 0.5 * (lc + log_deg + int(k) * log_n - log_m)
        coeff = (rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m))) / math.sqrt(2.0)
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None] * math.sqrt(float(m))
        obj = cls(n, d, seed, exps, exps.astype(np.float64).T, log_scales, coeff.astype(np.complex128), m, bool(projective_normalize), bool(dynamic_normalize), bool(equation_normalize), max(1, int(eval_block)))
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def _eval_block(self, ZZ: Any) -> Any:
        ZZ = np.asarray(ZZ, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        az = np.abs(ZZ)
        log_amp = np.log(np.maximum(az, 1e-300)) @ self.feature_exps_t + self.feature_log_scales[None, :]
        phase_arg = np.angle(ZZ) @ self.feature_exps_t
        if self.projective_normalize:
            log_amp -= (0.5 * float(self.d) * np.log1p(np.sum(az * az, axis=1)))[:, None]
        if self.dynamic_normalize:
            shift = np.max(np.where(np.isfinite(log_amp), log_amp, -np.inf), axis=1)
            log_amp = log_amp - np.where(np.isfinite(shift), shift, 0.0)[:, None]
            log_amp = np.clip(log_amp, -745.0, 0.0)
        else:
            log_amp = np.clip(log_amp, -745.0, 700.0)
        Phi = np.exp(log_amp) * np.exp(1j * phase_arg)
        out = Phi @ self.coeff.T
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        return out

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = self._eval_block(z)[0]
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        chunks = [self._eval_block(ZZ[i:i + self.eval_block]) for i in range(0, int(ZZ.shape[0]), self.eval_block)]
        out = np.vstack(chunks) if chunks else np.empty((0, self.n), dtype=np.complex128)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "lazy_features": int(self.lazy_features)}


@dataclasses.dataclass
class GeometryKernelKostlanSystem:
    n: int
    d: int
    seed: int
    anchors: Any
    anchor_conj_t: Any
    anchor_den: Any
    coeff: Any
    geometry_anchors: int
    anchor_scales: list[float]
    dynamic_normalize: bool = True
    self_normalize: bool = True
    equation_normalize: bool = True
    eval_block: int = 128
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True, geometry_anchors: int = 0, geometry_anchor_cap: int = 4096, geometry_anchor_scales: Optional[Sequence[float]] = None, dynamic_normalize: bool = True, self_normalize: bool = True, eval_block: int = 128) -> "GeometryKernelKostlanSystem":
        t0 = now()
        n, d = int(n), int(d)
        m = max(n + 2, int(geometry_anchors) if int(geometry_anchors) > 0 else auto_geometry_anchor_count(n, d, int(geometry_anchor_cap)))
        seed = stable_seed(n, d, seed_index, salt=0x314C0DE)
        rng = np.random.default_rng(seed)
        scales = [float(x) for x in (geometry_anchor_scales or []) if math.isfinite(float(x)) and float(x) > 0] or [0.25, 0.5, 1.0, 2.0, 4.0]
        anchors = np.zeros((m, n), dtype=np.complex128)
        idx = 1
        axis_scale = scales[min(2, len(scales) - 1)]
        for j in range(n):
            if idx >= m:
                break
            anchors[idx, j] = complex(axis_scale)
            idx += 1
        sqrt_n = math.sqrt(max(1, n))
        while idx < m:
            shell = scales[(idx - 1 - n) % len(scales)]
            v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
            v = v / max(1e-300, float(np.linalg.norm(v))) * (shell * sqrt_n)
            anchors[idx, :] = v
            idx += 1
        coeff = (rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m))) / math.sqrt(2.0 * max(1, m))
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None]
        anchor_norm2 = np.sum(np.abs(anchors) ** 2, axis=1)
        obj = cls(n, d, seed, anchors, np.conjugate(anchors).T, np.sqrt(1.0 + anchor_norm2), coeff.astype(np.complex128), m, [float(x) for x in scales], bool(dynamic_normalize), bool(self_normalize), bool(equation_normalize), max(1, int(eval_block)))
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def _kernel_block(self, ZZ: Any) -> Any:
        ZZ = np.asarray(ZZ, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        dot = ZZ @ self.anchor_conj_t
        zn = np.sqrt(1.0 + np.sum(np.abs(ZZ) ** 2, axis=1))
        base = (1.0 + dot) / np.maximum(1e-300, zn[:, None] * self.anchor_den[None, :])
        mag = np.abs(base)
        base = np.where(mag > 1.0, base / np.maximum(mag, 1e-300), base)
        mag = np.minimum(mag, 1.0)
        log_amp = float(self.d) * np.log(np.maximum(mag, 1e-300))
        phase_arg = float(self.d) * np.angle(base)
        if self.dynamic_normalize:
            shift = np.max(np.where(np.isfinite(log_amp), log_amp, -np.inf), axis=1)
            log_amp = log_amp - np.where(np.isfinite(shift), shift, 0.0)[:, None]
        K = np.exp(np.clip(log_amp, -745.0, 0.0)) * np.exp(1j * phase_arg)
        if self.self_normalize:
            row = np.sqrt(np.mean(np.abs(K) ** 2, axis=1))
            K = K / np.maximum(row[:, None], 1e-300)
        return K

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = (self._kernel_block(z) @ self.coeff.T)[0]
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        chunks = []
        for i in range(0, int(ZZ.shape[0]), self.eval_block):
            F = self._kernel_block(ZZ[i:i + self.eval_block]) @ self.coeff.T
            F[~np.isfinite(F)] = complex(1e300, 0.0)
            chunks.append(F)
        out = np.vstack(chunks) if chunks else np.empty((0, self.n), dtype=np.complex128)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "geometry_anchors": int(self.geometry_anchors), "geometry_anchor_scales": list(self.anchor_scales)}


def make_kostlan_base(args: argparse.Namespace, n: int, d: int) -> Any:
    mode = str(getattr(args, "system_mode", "auto")).strip().lower().replace("_", "-")
    dense_terms = exact_kostlan_terms(n, d)
    if mode == "dense" or (mode == "auto" and dense_terms <= int(args.dense_max_terms)):
        return DenseKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize))
    if mode in {"auto", "geometry", "geometry-kernel", "kernel", "projective-kernel"}:
        return GeometryKernelKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize), int(args.geometry_anchors), int(args.geometry_anchor_cap), parse_float_list(args.geometry_anchor_scales, [0.25, 0.5, 1.0, 2.0, 4.0], positive=True), bool(args.geometry_dynamic_normalize), bool(args.geometry_self_normalize), int(args.geometry_eval_block))
    if mode not in {"lazy", "lazy-feature", "feature", "stream"}:
        raise ValueError(f"unknown system mode {mode!r}")
    return LazyFeatureKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize), int(args.lazy_features), int(args.lazy_feature_cap), bool(args.lazy_projective_normalize), bool(args.lazy_dynamic_normalize), int(args.lazy_eval_block))


# ---------------------------------------------------------------------------
# Exact expression polynomial systems (used only as a black-box geometric oracle)
# ---------------------------------------------------------------------------

class SafeExpression:
    _allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Add, ast.Sub, ast.Mult, ast.Pow, ast.USub, ast.UAdd, ast.Load)

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
            if isinstance(node.right.value, complex) or float(node.right.value) != int(node.right.value):
                raise ValueError("polynomial powers must be integer constants")
            if int(node.right.value) < 0:
                raise ValueError("polynomial powers must be non-negative")
        for child in ast.iter_child_nodes(node):
            self._validate(child)

    def eval(self, env: dict[str, Any]) -> Any:
        return self._eval(self.tree.body, env)

    def _eval(self, node: ast.AST, env: dict[str, Any]) -> Any:
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id not in env:
                raise ValueError(f"unknown variable {node.id!r}")
            return env[node.id]
        if isinstance(node, ast.UnaryOp):
            v = self._eval(node.operand, env)
            return -v if isinstance(node.op, ast.USub) else v
        if isinstance(node, ast.BinOp):
            a, b = self._eval(node.left, env), self._eval(node.right, env)
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                return a * b
            if isinstance(node.op, ast.Div):
                return a / b
            if isinstance(node.op, ast.Pow):
                return a ** int(b)
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
    def make(cls, n: int, d: int, raw: str, variable_names: Optional[Sequence[str]] = None, seed_index: int = 0) -> "ExpressionPolynomialSystem":
        parts = [p.strip() for p in str(raw).replace("|", ";").split(";") if p.strip()]
        if not parts:
            raise ValueError("--polys/--poly is required for polynomial systems")
        names = [str(p).strip() for p in (variable_names or []) if str(p).strip()] or (["x"] if int(n) == 1 else [f"x{i + 1}" for i in range(int(n))])
        if len(names) != int(n) or len(parts) != int(n):
            raise ValueError(f"expected {n} variables and {n} polynomials")
        return cls(int(n), int(d), stable_seed(n, d, seed_index, salt=0x315EAC7), parts, [SafeExpression(p) for p in parts], names)

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return 0.0

    def _env(self, ZZ: Any) -> dict[str, Any]:
        env: dict[str, Any] = {"I": 1j, "j": 1j}
        for i, name in enumerate(self.variable_names):
            env[name] = ZZ[:, i]
        if self.n == 1:
            env.setdefault("x", ZZ[:, 0])
            env.setdefault("z", ZZ[:, 0])
        for i in range(self.n):
            env.setdefault(f"x{i + 1}", ZZ[:, i])
            env.setdefault(f"z{i + 1}", ZZ[:, i])
        return env

    def eval(self, z: Sequence[complex]) -> Any:
        return self.eval_batch(np.asarray(z, dtype=np.complex128)[None, :])[0]

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        env = self._env(ZZ)
        cols = []
        for expr in self.expressions:
            val = expr.eval(env)
            arr = np.asarray(val, dtype=np.complex128)
            if arr.ndim == 0:
                arr = np.full(int(ZZ.shape[0]), complex(arr), dtype=np.complex128)
            cols.append(arr.reshape(int(ZZ.shape[0])))
        out = np.stack(cols, axis=1)
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "expressions": list(self.expressions_raw), "variables": list(self.variable_names)}


# ---------------------------------------------------------------------------
# The local-jet geometry -- the heart of 317
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class LocalJetGeometry:
    """The geometry that faithfully represents a system: its field of local jets.

    The local-jet layer adds no global precompute.  The base system is treated
    as a black-box geometric oracle; around any point the geometry is described by a
    locally sampled jet (value + affine Jacobian + optional curvature) built from a
    small paired stencil of exactly 2n+1 oracle queries.  Per correction step the oracle is hit
    O(n) times, without enumerating the global comb(n+d, d) monomial budget.

    Because jets are sampled from the TRUE system, eval()/residual() are faithful;
    geometric acceptance therefore coincides with exact acceptance.  Optionally, a
    small jet cache reuses jets across nearby queries.
    """

    base: Any
    n: int
    d: int
    seed: int
    use_quadratic: bool = True
    cache_enabled: bool = True
    cache_decimals: int = 9
    cache_cap: int = 4096
    jet_radius: float = 1e-5
    eval_count: int = 0
    slope_count: int = 0
    oracle_samples: int = 0
    jets_built: int = 0
    cache_hits: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, base: Any, n: int, d: int, args: argparse.Namespace) -> "LocalJetGeometry":
        obj = cls(
            base, int(n), int(d), int(getattr(base, "seed", stable_seed(n, d, int(args.seed_index)))),
            bool(args.jet_quadratic), bool(args.jet_cache), max(1, int(args.jet_cache_decimals)),
            max(16, int(getattr(args, "jet_cache_cap", 4096))),
            max(1e-12, float(args.jet_radius)),
        )
        obj._cache = OrderedDict()
        return obj

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
        return float(getattr(self.base, "generation_seconds", 0.0))

    @property
    def samples_per_jet(self) -> int:
        return int(2 * self.n + 1)

    # The oracle is queried only through eval/eval_batch -- never via the algebra.
    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = self.base.eval(np.asarray(z, dtype=np.complex128))
        self.eval_count += 1
        self.oracle_samples += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        out = self.base.eval_batch(ZZ)
        self.eval_count += int(ZZ.shape[0])
        self.oracle_samples += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def local_jet(self, z: Sequence[complex], radius: Optional[float] = None) -> dict[str, Any]:
        """Sample the oracle on a small hypercube around z and fit a local jet.

        Returns {f0, J, (Q,) z}, the explicit local geometry at z.  This is the same
        object the hypercube inverse-jet corrector builds; exposing it makes the
        jet-field nature of the geometry concrete and supports a jet-Newton polish.
        """
        zz = np.asarray(z, dtype=np.complex128)
        key = None
        base_h = float(radius) if radius is not None else self.jet_radius
        if self.cache_enabled:
            rounded_z = tuple(np.round(np.concatenate([zz.real, zz.imag]), self.cache_decimals))
            key = (float(base_h),) + rounded_z
            hit = self._cache.get(key)
            if hit is not None:
                self._cache.move_to_end(key)
                self.cache_hits += 1
                return hit
        # [318-A] per-coordinate radius: better conditioned when coordinates have
        # disparate magnitudes; the column-wise division keeps the O(h^2) accuracy.
        hvec = base_h * np.maximum(1.0, np.abs(zz))
        f0 = self.eval(zz)
        E = np.diag(hvec).astype(np.complex128)
        Pp = zz[None, :] + E
        Pm = zz[None, :] - E
        Fp = self.eval_batch(Pp)
        Fm = self.eval_batch(Pm)
        J = ((Fp - Fm) / (2.0 * hvec[:, None])).T
        jet: dict[str, Any] = {"z": zz, "f0": f0, "J": J, "h": float(np.max(hvec)), "hvec": hvec}
        if self.use_quadratic:
            jet["Qdiag"] = ((Fp - 2.0 * f0[None, :] + Fm) / (hvec[:, None] ** 2)).T
        self.jets_built += 1
        if key is not None:
            self._cache[key] = jet
            self._cache.move_to_end(key)
            while len(self._cache) > int(self.cache_cap):
                self._cache.popitem(last=False)
        return jet

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        t0 = now()
        jet = self.local_jet(np.asarray(a, dtype=np.complex128))
        self.slope_count += 1
        self.seconds_slope += now() - t0
        return jet["J"]

    def stats(self) -> dict[str, Any]:
        return {
            "eval_count": int(self.eval_count), "slope_count": int(self.slope_count),
            "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope),
            "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms,
            "geometry": "local-jet-field", "oracle_samples": int(self.oracle_samples),
            "jets_built": int(self.jets_built), "jet_cache_hits": int(self.cache_hits),
            "samples_per_jet": int(self.samples_per_jet), "use_quadratic": bool(self.use_quadratic),
            "base_stats": self.base.stats(),
        }


def describe_base_backend(base: Any) -> str:
    if isinstance(base, ExpressionPolynomialSystem):
        return "polynomial-exact"
    if isinstance(base, DenseKostlanSystem):
        return "ks-dense"
    if isinstance(base, LazyFeatureKostlanSystem):
        return "ks-lazy-feature"
    if isinstance(base, GeometryKernelKostlanSystem):
        return "ks-geometry-kernel"
    return "unknown"


def make_system_317(args: argparse.Namespace, n: int, d: int) -> tuple[LocalJetGeometry, Any, str, str]:
    source = str(args.system_source).strip().lower().replace("_", "-")
    if source in {"poly", "polynomial", "expr", "expression"}:
        vars_raw = str(args.variables or "").strip()
        variables = [p.strip() for p in vars_raw.replace(";", ",").split(",") if p.strip()] if vars_raw else None
        base: Any = ExpressionPolynomialSystem.make(n, d, str(args.polys or ""), variables, int(args.seed_index))
        source = "polynomial"
    elif source in {"ks", "kostlan"}:
        base = make_kostlan_base(args, n, d)
        source = "ks"
    else:
        raise ValueError(f"unknown --system-source {source!r}")
    geometry = LocalJetGeometry.make(base, n, d, args)
    return geometry, base, source, describe_base_backend(base)


# ---------------------------------------------------------------------------
# Charts, starts, and correctors (operate on the geometry through a generic target)
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class LinearChart:
    A: Any
    Ainv: Any

    @classmethod
    def identity(cls, n: int, scale: float = 1.0) -> "LinearChart":
        A = np.eye(int(n), dtype=np.complex128) * complex(scale)
        Ainv = np.eye(int(n), dtype=np.complex128) / complex(scale)
        return cls(A, Ainv)

    def z_from_y(self, y: Sequence[complex]) -> Any:
        return self.A @ np.asarray(y, dtype=np.complex128)

    def y_from_z(self, z: Sequence[complex]) -> Any:
        return self.Ainv @ np.asarray(z, dtype=np.complex128)


@dataclasses.dataclass
class TargetTrack:
    system: Any
    chart: LinearChart

    def eval(self, y: Sequence[complex]) -> Any:
        return self.system.eval(self.chart.z_from_y(y))

    def eval_batch(self, Y: Any) -> Any:
        YY = np.asarray(Y, dtype=np.complex128)
        if YY.ndim == 1:
            return self.eval(YY)[None, :]
        return self.system.eval_batch(YY @ np.asarray(self.chart.A, dtype=np.complex128).T)

    def residual(self, y: Sequence[complex]) -> float:
        return finite_norm(self.eval(y))

    def residuals_batch(self, Y: Any) -> Any:
        try:
            return safe_norms(self.eval_batch(Y))
        except Exception:
            YY = np.asarray(Y, dtype=np.complex128)
            return np.full(1 if YY.ndim == 1 else int(YY.shape[0]), np.inf, dtype=float)


DEFAULT_POWERS = [2.0 ** k for k in range(-20, 21)] + [3.0 * (2.0 ** k) for k in range(-12, 18)] + [10.0 ** k for k in range(-6, 7)]
DEFAULT_ANGLES_DEG = [0, 6, 12, 18, 24, 32, 40, 48, 56, 64, 72, 80, 86, 89, 90, 91, 94, 100, 108, 116, 128, 140, 152, 164, 172]


def monomial_scale_ladder(subdivisions: int, octaves: int, base: float = 2.0) -> list[float]:
    """The 'complete logarithmic ladder': equally spaced logarithmic scale offsets.

    Between consecutive powers ``base**m`` the ladder inserts the ``subdivisions - 1`` geometric
    (proportional) means ``base**(m + k/subdivisions)``, k = 1..subdivisions-1.  For base 2 and
    subdivisions 3 those means are the cube-root proportional means 2**(1/3), 2**(2/3) -- exactly
    the x**(1/p), x**(2/p), ... ladder of the Pandrosion construction (the diagram that yields
    cbrt(2) yields cbrt(4) = cbrt(2)**2 in the same figure).

    The monomial scale-palette theorem proves that equally spaced logarithmic offsets minimise the
    worst-case raw residual multiplier 1 - p/S_p(q**(1/(2K))).  Using the ladder as the multi-start
    homothety palette (and as the IRP chart palette) therefore makes |log y| uniformly small before
    iteration, so a start is more likely to land inside a convergence basin.
    """
    p = max(1, int(subdivisions))
    m = max(1, int(octaves))
    b = float(base) if float(base) > 1.0 else 2.0
    return [b ** (j / p) for j in range(-m * p, m * p + 1)]
DEFAULT_RADII = [0.025, 0.04, 0.06, 0.08, 0.12, 0.18, 0.27, 0.4, 0.6, 0.85, 1.15, 1.55, 2.05, 2.75, 3.6, 4.8, 6.4]
DEFAULT_GAINS = [0.035, 0.055, 0.085, 0.12, 0.18, 0.27, 0.4, 0.58, 0.78, 1.0, 1.28, 1.65, 2.2, 3.0, 4.2, 6.0, 8.5, 12.0]


def raw_direction(n: int, trial: int, seed: int) -> Any:
    vals = []
    for j in range(int(n)):
        h1 = splitmix64(seed + 0xD1A5E + 1000003 * trial + 4099 * (j + 1))
        h2 = splitmix64(seed + 0xBADC0DE + 1000033 * trial + 9176 * (j + 1))
        vals.append(math.exp(0.45 * (2.0 * u01(h2) - 1.0)) * phase(2.0 * math.pi * u01(h1)))
    v = np.asarray(vals, dtype=np.complex128)
    return v / max(1e-300, float(np.linalg.norm(v))) * math.sqrt(max(1, int(n)))


def universal_atlas_start(target: TargetTrack, n: int, trial: int, seed: int, powers: Sequence[float], angles: Sequence[float], radii: Sequence[float], cap: float, roots_found: int, duplicates: int, failures: int, target_count: int, universal_cells: int = 16, universal_shells: int = 5, atlas_selection: str = "diverse-shell", **_: Any) -> tuple[Any, dict[str, Any]]:
    cells = max(1, int(universal_cells))
    shells = max(1, int(universal_shells))
    pows = [min(max(float(x), 1e-300), float(cap)) for x in powers if float(x) > 0] or [1.0]
    mode = str(atlas_selection or "diverse-shell").strip().lower().replace("_", "-")
    if mode in {"diverse", "diverse-shell", "shell", "raw-shell"}:
        rr = [float(x) for x in radii if math.isfinite(float(x)) and float(x) > 0] or DEFAULT_RADII
        idx = int(trial)
        rad = rr[(idx + 3 * roots_found + failures) % len(rr)]
        y = float(rad) * raw_direction(n, idx, seed + 65537 * idx + 104729 * duplicates)
        r0 = target.residual(y)
        meta = {
            "chart": "320-standalone-universal-diverse-shell-atlas",
            "atlas_mode": "320-diverse-shell",
            "atlas_selection": mode,
            "atlas_startopt_bypass_recommended": True,
            "atlas_cells_tested": 1,
            "atlas_admissible_cells": int(math.isfinite(r0)),
            "atlas_cell_residual": float(r0),
            "atlas_selected_index": int(idx),
            "atlas_selected_layer": int(idx % shells),
            "homothety": float(rad),
            "base_homothety": 1.0,
            "thales_thrust": 1.0,
            "theta_deg": None,
            "base_radius": float(rad),
            "dup_pressure": float((duplicates + 1.0) / (roots_found + 1.0)),
            "fail_pressure": float((failures + 1.0) / (trial + 1.0)),
            "progress": float(min(1.0, roots_found / max(1.0, float(target_count)))),
        }
        return np.asarray(y, dtype=np.complex128).copy(), meta

    candidates = []
    metas = []
    for k in range(cells):
        idx = trial * cells + k
        layer = (idx // cells) % shells
        pwr = pows[(idx * 37 + 11 * layer + roots_found) % len(pows)]
        rad = radii[(idx * 13 + failures) % len(radii)] if radii else 1.0
        theta = angles[(idx * 19 + duplicates) % len(angles)] if angles else 0.0
        thrust = [1.0, 1.6, 2.5, 4.0, 6.5, 10.0, 16.0][idx % 7]
        lam = min(float(cap), pwr * (thrust ** (0.2 + 0.8 * min(1.0, roots_found / max(1.0, float(target_count))))))
        d = raw_direction(n, idx, seed + 65537 * k)
        y = rad * d
        c, s = math.cos(theta), math.sin(theta)
        pole = phase(0.37 * (idx + 1))
        w = y / pole
        den = -s * w + c
        den = np.where(np.abs(den) < 1e-12, den + 1e-12, den)
        y = lam * pole * ((c * w + s) / den)
        candidates.append(np.asarray(y, dtype=np.complex128))
        metas.append({"homothety": float(lam), "base_homothety": float(pwr), "thales_thrust": float(thrust), "theta_deg": float(math.degrees(theta)), "base_radius": float(rad), "atlas_selected_index": int(idx), "atlas_selected_layer": int(layer)})
    Y = np.asarray(candidates, dtype=np.complex128)
    R = target.residuals_batch(Y)
    idx_best = int(np.nanargmin(R)) if np.any(np.isfinite(R)) else 0
    meta = dict(metas[idx_best])
    meta.update({"chart": "320-standalone-universal-mobius-kernel-atlas", "atlas_mode": "320-compact-batched-score", "atlas_selection": mode, "atlas_startopt_bypass_recommended": False, "atlas_cells_tested": int(cells), "atlas_admissible_cells": int(np.sum(np.isfinite(R))), "atlas_cell_residual": float(R[idx_best]) if idx_best < len(R) else float("inf"), "dup_pressure": float((duplicates + 1.0) / (roots_found + 1.0)), "fail_pressure": float((failures + 1.0) / (trial + 1.0)), "progress": float(min(1.0, roots_found / max(1.0, float(target_count))))})
    return Y[idx_best].copy(), meta


def origin_affine_start(target: TargetTrack, n: int, h: float, max_norm: float) -> tuple[Optional[Any], dict[str, Any]]:
    h0 = max(1e-12, float(h))
    y0 = np.zeros(int(n), dtype=np.complex128)
    try:
        f0 = target.eval(y0)
        P = np.eye(int(n), dtype=np.complex128) * h0
        J = ((target.eval_batch(P) - f0[None, :]) / h0).T
        try:
            y = np.linalg.solve(J, -f0)
            method = "solve"
        except Exception:
            y, _, _, _ = np.linalg.lstsq(J, -f0, rcond=None)
            method = "lstsq"
        if not np.all(np.isfinite(y)):
            return None, {"origin_seed_enabled": True, "origin_seed_status": "nonfinite"}
        cap = float(max_norm)
        yn = float(np.linalg.norm(y))
        if math.isfinite(cap) and cap > 0 and yn > cap:
            y = y * (cap / max(yn, 1e-300))
        return y, {"origin_seed_enabled": True, "origin_seed_status": "ok", "origin_seed_method": method, "origin_seed_h": float(h0), "origin_seed_r0": finite_norm(f0), "origin_seed_r1": target.residual(y), "origin_seed_norm": float(np.linalg.norm(y))}
    except Exception as exc:
        return None, {"origin_seed_enabled": True, "origin_seed_status": f"error:{type(exc).__name__}", "origin_seed_h": float(h0)}


def startopt(target: TargetTrack, y0: Any, trial: int, seed: int, steps: int, candidates: int, gains: Sequence[float], micro_epochs: int) -> tuple[Any, dict[str, Any]]:
    best = np.asarray(y0, dtype=np.complex128).copy()
    best_r = target.residual(best)
    initial = best_r
    evals = 1
    chosen_gain = 1.0
    for step in range(max(0, int(steps))):
        pool = [(1.0, best)]
        for c in range(max(0, int(candidates) - 1)):
            gain = float(gains[(trial + 3 * step + c) % len(gains)])
            cand = gain * best
            if c % 3 == 1:
                noise = raw_direction(len(best), trial + 31 * step + c, seed)
                cand = cand + 0.08 * max(1.0, float(np.linalg.norm(best))) * noise
            elif c % 3 == 2:
                cand = cand * phase(0.11 * (c + 1))
            pool.append((gain, np.asarray(cand, dtype=np.complex128)))
        Y = np.asarray([p[1] for p in pool], dtype=np.complex128)
        R = target.residuals_batch(Y)
        evals += len(pool)
        idx = int(np.nanargmin(R)) if np.any(np.isfinite(R)) else 0
        if float(R[idx]) < best_r:
            best, best_r, chosen_gain = Y[idx].copy(), float(R[idx]), float(pool[idx][0])
    micro_done = 0
    for _ in range(max(0, int(micro_epochs))):
        f0 = target.eval(best)
        evals += 1
        _, J = batched_paired_jacobians(target, best[None, :], F0=f0[None, :])
        evals += 2 * len(best)
        try:
            delta = np.linalg.solve(J[0], -f0)
        except Exception:
            delta = np.linalg.pinv(J[0]) @ (-f0)
        if not np.all(np.isfinite(delta)):
            break
        Y = np.asarray([best + lam * delta for lam in (1.0, 0.5, 0.25)], dtype=np.complex128)
        R = target.residuals_batch(Y)
        evals += len(Y)
        idx = int(np.argmin(np.where(np.isfinite(R), R, np.inf)))
        if not math.isfinite(float(R[idx])) or float(R[idx]) >= best_r:
            break
        best, best_r = Y[idx].copy(), float(R[idx])
        micro_done += 1
    return best, {"startopt_enabled": bool(steps > 0 or micro_epochs > 0), "startopt_r0": float(initial), "startopt_r1": float(best_r), "startopt_ratio": (float(best_r / initial) if math.isfinite(best_r) and math.isfinite(initial) and initial > 0 else None), "startopt_steps": int(max(0, steps)), "startopt_evals": int(evals), "startopt_micro_epochs": int(max(0, micro_epochs)), "startopt_micro_epochs_completed": int(micro_done), "startopt_gain": float(chosen_gain), "startopt_batch_numpy": True}



def batched_paired_jacobians(target: TargetTrack, Y: Any, h_rel: float = 1e-5, F0: Optional[Any] = None) -> tuple[Any, Any]:
    """[318-B] Central-difference Jacobians for a whole batch of points at once.

    One stacked eval_batch for all +/- probes of all points: the oracle cost per
    point is identical to the sequential paired jet (2n+1) but the Python/dispatch
    overhead is amortised across the batch.  [320-G] a known F0 (residual vectors
    already evaluated at Y) can be handed in, saving B oracle samples.
    """
    YY = np.asarray(Y, dtype=np.complex128)
    B, n = int(YY.shape[0]), int(YY.shape[1])
    F0 = target.eval_batch(YY) if F0 is None else np.asarray(F0, dtype=np.complex128)
    H = float(h_rel) * np.maximum(1.0, np.abs(YY))
    idx = np.arange(n)
    Pp = np.repeat(YY[:, None, :], n, axis=1)
    Pp[:, idx, idx] += H
    Pm = np.repeat(YY[:, None, :], n, axis=1)
    Pm[:, idx, idx] -= H
    F = target.eval_batch(np.vstack([Pp.reshape(B * n, n), Pm.reshape(B * n, n)]))
    Fp = F[: B * n].reshape(B, n, n)
    Fm = F[B * n :].reshape(B, n, n)
    J = np.transpose(Fp - Fm, (0, 2, 1)) / (2.0 * H[:, None, :])
    return F0, J

def swarm_starts(n: int, size: int, seed: int, radii: Sequence[float]) -> Any:
    rr = [float(x) for x in radii if math.isfinite(float(x)) and float(x) > 0] or DEFAULT_RADII
    return np.stack([rr[i % len(rr)] * raw_direction(n, i, seed + 65537 * i) for i in range(max(1, int(size)))])



def swarm_prefilter(target: TargetTrack, n: int, seed: int, radii: Sequence[float], size: int, iters: int, keep: int, line_lams: Sequence[float], sep: float) -> tuple[list[Any], dict[str, Any]]:
    """[318-B] Vectorised damped-Newton prefilter over a batch of atlas starts.

    Runs a few batched Newton iterations on `size` diverse-shell starts at once
    (stacked Jacobians, batched solve with pinv fallback, batched line search),
    culls the points that stop descending, and returns up to `keep` survivors
    ordered by residual under a pairwise diversity separation, so the survivors
    sample DISTINCT basins rather than `keep` copies of the easiest one.

    [320-G] residual VECTORS from the line search are reused as the next F0.
    Jacobians deliberately remain exact paired jets on every swarm iteration;
    batched Broyden was rejected because it degraded basin separation.
    """
    t0 = now()
    Y = swarm_starts(n, size, seed, radii)
    B = int(Y.shape[0])
    F = target.eval_batch(Y)
    R = safe_norms(F)
    alive = np.isfinite(R)
    L = np.asarray([float(x) for x in line_lams if math.isfinite(float(x)) and float(x) > 0] or [1.0, 0.5, 0.25, 0.1], dtype=float)
    iters_done = 0
    jets_full = broyden_upd = 0
    for it in range(max(0, int(iters))):
        idx = np.where(alive)[0]
        if idx.size == 0:
            break
        # [320-G] the residual vectors carried from the line search serve as F0, so
        # each iteration saves one full-batch oracle evaluation; the Jacobians stay
        # EXACT paired jets (a batched secant variant was tried and measurably
        # degraded basin separation, inflating duplicates).
        _, J = batched_paired_jacobians(target, Y[idx], F0=F[idx])
        jets_full += int(idx.size)
        F0 = F[idx]
        with np.errstate(all="ignore"):
            try:
                delta = np.linalg.solve(J, -F0[..., None])[..., 0]
            except np.linalg.LinAlgError:
                delta = np.matmul(np.linalg.pinv(J), -F0[..., None])[..., 0]
        bad = ~np.all(np.isfinite(delta), axis=1)
        if np.any(bad):
            delta[bad] = 0.0
        nd = np.linalg.norm(delta, axis=1)
        yn = np.maximum(1.0, np.linalg.norm(Y[idx], axis=1))
        lim = 10.0 * yn
        scl = np.where(nd > lim, lim / np.maximum(nd, 1e-300), 1.0)
        delta = delta * scl[:, None]
        cand = Y[idx][:, None, :] + L[None, :, None] * delta[:, None, :]
        Fc = target.eval_batch(cand.reshape(-1, n)).reshape(idx.size, int(L.size), n)
        Rc = safe_norms(Fc.reshape(-1, n)).reshape(idx.size, int(L.size))
        bestl = np.argmin(np.where(np.isfinite(Rc), Rc, np.inf), axis=1)
        rows = np.arange(idx.size)
        Rbest = Rc[rows, bestl]
        improved = np.isfinite(Rbest) & (Rbest < R[idx])
        if np.any(improved):
            sel = idx[improved]
            Ynew = cand[rows, bestl][improved]
            Fnew = Fc[rows, bestl][improved]
            Y[sel] = Ynew
            F[sel] = Fnew
            R[sel] = Rbest[improved]
        if it >= 1:
            alive[idx[~improved]] = False
        iters_done = it + 1
    order = np.argsort(np.where(np.isfinite(R), R, np.inf))
    picked: list[int] = []
    for i in order:
        if not np.isfinite(R[int(i)]):
            break
        if all(float(np.linalg.norm(Y[int(i)] - Y[j])) > float(sep) for j in picked):
            picked.append(int(i))
        if len(picked) >= max(1, int(keep)):
            break
    meta = {
        "swarm_enabled": True, "swarm_size": B, "swarm_iters": iters_done,
        "swarm_survivors": int(np.sum(alive)), "swarm_kept": len(picked),
        "swarm_sep": float(sep), "swarm_best_residual": (float(R[picked[0]]) if picked else None),
        "swarm_full_jets": int(jets_full), "swarm_broyden_updates": int(broyden_upd),
        "swarm_seconds": float(now() - t0),
    }
    return [Y[i].copy() for i in picked], meta

def singular_direction_starts(geometry: Any, n: int, count: int, max_norm: float) -> list[Any]:
    """[318-F] Starts along the smallest right singular directions of the origin jet.

    The smallest singular directions of J(0) are the locally flattest directions of
    the system -- where the affine model predicts the zero set extends farthest --
    so seeding the Newton point of the origin jet plus excursions along them probes
    basins that radial atlas shells reach only by luck.
    """
    if int(count) <= 0:
        return []
    try:
        jet0 = geometry.local_jet(np.zeros(int(n), dtype=np.complex128))
        J0, f00 = np.asarray(jet0["J"], dtype=np.complex128), np.asarray(jet0["f0"], dtype=np.complex128)
        U, S, Vh = np.linalg.svd(J0)
        try:
            ybase = np.linalg.solve(J0, -f00)
        except Exception:
            cut = float(S[0]) * 1e-12 if S.size else 0.0
            Sinv = np.where(S > cut, 1.0 / np.maximum(S, 1e-300), 0.0)
            ybase = Vh.conj().T @ (Sinv * (U.conj().T @ (-f00)))
        if not np.all(np.isfinite(ybase)):
            return []
        cap = float(max_norm) if float(max_norm) > 0 else 2.0 * math.sqrt(max(1, n))
        bn = float(np.linalg.norm(ybase))
        if bn > cap:
            ybase = ybase * (cap / max(bn, 1e-300))
        out: list[Any] = []
        k = min(int(count), 2 * int(n))
        amp = 0.8 * math.sqrt(max(1, n))
        for i in range(k):
            v = Vh[int(n) - 1 - (i // 2)].conj()
            out.append(np.asarray(ybase + (amp if i % 2 == 0 else -amp) * v, dtype=np.complex128))
        return out
    except Exception:
        return []


def _line_lambdas(line_search: int, line_grid: Sequence[float]) -> list[float]:
    vals = [float(x) for x in line_grid if math.isfinite(float(x)) and float(x) > 0]
    if vals:
        return vals
    return [1.0, 0.75, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03125, 0.02][:max(1, int(line_search))]


def deflation_log_weights(Y: Any, known: Sequence[Any], alpha: float) -> Any:
    YY = np.asarray(Y, dtype=np.complex128)
    if YY.ndim == 1:
        YY = YY[None, :]
    log_w = np.zeros(int(YY.shape[0]), dtype=float)
    for ky in known:
        dist = np.linalg.norm(YY - np.asarray(ky, dtype=np.complex128)[None, :], axis=1)
        log_w += np.log1p(float(alpha) / np.maximum(dist, 1e-12))
    return log_w


def deflated_line_choice(R: Any, Y: Any, y: Any, current_r: float, best_r: float, known: Sequence[Any], alpha: float) -> tuple[Optional[int], Any, Any]:
    """Return the best admissible line candidate using overflow-safe merit."""
    RR = np.asarray(R, dtype=float)
    log_w = deflation_log_weights(Y, known, alpha) if known else np.zeros(int(RR.size), dtype=float)
    log_w0 = float(deflation_log_weights(np.asarray(y)[None, :], known, alpha)[0]) if known else 0.0
    log_merit = np.where(np.isfinite(RR), np.log(np.maximum(RR, 1e-300)) + log_w, np.inf)
    log_merit0 = math.log(max(float(current_r), 1e-300)) + log_w0
    admissible = np.isfinite(RR) & ((log_merit < log_merit0) | (RR < float(best_r)))
    if not np.any(admissible):
        return None, log_w, log_merit
    return int(np.argmin(np.where(admissible, log_merit, np.inf))), log_w, log_merit



def hypercube_delta(target: TargetTrack, y: Any, f: Any, ep: int, seed: int, cloud_nodes: int = 0, lm_damping: float = 0.0, trust_radius: float = 0.0, sampling: str = "paired", A_pre: Optional[Any] = None, lm_mu: Optional[float] = None, max_order: int = 4, return_jac: bool = False) -> tuple[Any, dict[str, Any]]:
    n = len(y)
    yn = max(1.0, float(np.linalg.norm(y)))
    mode = str(sampling or "paired").strip().lower()
    if A_pre is not None:
        # [320-A] recycled (Broyden-updated) Jacobian: ZERO oracle samples spent here.
        A = np.asarray(A_pre, dtype=np.complex128)
        M = 0
        oracle_evals = 0
        mode = "broyden-reuse"
    elif mode == "cloud":
        # Legacy 317 random-sign cloud (kept for A/B comparison).
        M = max(int(cloud_nodes), 2 * n + 4, 16)
        h = 1e-5 * yn
        rng = np.random.default_rng(int(seed) + int(ep) * 1337)
        signs = rng.choice([-1.0, 1.0], size=(M, n))
        dY = h * signs
        dF = target.eval_batch(y[None, :] + dY) - f[None, :]
        oracle_evals = M
        A, _, _, _ = np.linalg.lstsq(dY, dF, rcond=None)
        A = A.T
    else:
        # [318-A] paired central differences with per-coordinate radii: exact same
        # oracle budget class (2n samples) but O(h^2) truncation and deterministic.
        M = 2 * n + 1
        hvec = 1e-5 * np.maximum(1.0, np.abs(y))
        E = np.diag(hvec).astype(np.complex128)
        FpFm = target.eval_batch(np.vstack([y[None, :] + E, y[None, :] - E]))
        oracle_evals = 2 * n
        Fp, Fm = FpFm[:n], FpFm[n:]
        A = ((Fp - Fm) / (2.0 * hvec[:, None])).T
    if not np.all(np.isfinite(A)):
        return np.zeros_like(y), {"hypercube_error": "nonfinite_jacobian", "hypercube_order": 0, "hypercube_nodes": M, "hypercube_sampling": mode, "hypercube_oracle_evals": int(oracle_evals)}
    try:
        slope_cond = float(np.linalg.cond(A))
    except Exception:
        slope_cond = 1e300
    if not math.isfinite(slope_cond):
        slope_cond = 1e300
    jac_extra = {"jac": A} if return_jac else {}
    # Inverse-jet inversion of A.  Order of defences: (adaptive) LM damping ->
    # plain solve -> [318-A] truncated-SVD pseudo-inverse when the plain solve is
    # singular or yields an unusable step.  The trust radius caps whatever solver won.
    # [320-C] lm_mu, when provided, is the LIVE adaptive trust parameter and
    # overrides the static --lm-damping.
    lm_val = float(lm_mu) if (lm_mu is not None and float(lm_mu) > 0.0) else float(lm_damping)
    lm = lm_val > 0.0
    tr = float(trust_radius) > 0.0
    svd_used = False
    if lm:
        Ah = A.conj().T
        AhA = Ah @ A
        mu = lm_val * (float(np.trace(AhA).real) / max(1, n) + 1e-300)
        damped = AhA + mu * np.eye(n, dtype=np.complex128)

        def solve_step(rhs: Any) -> Any:
            return np.linalg.solve(damped, Ah @ rhs)
    else:
        def _svd_solve(rhs: Any) -> Any:
            U, S, Vh = np.linalg.svd(A)
            cut = float(S[0]) * 1e-12 if S.size else 0.0
            Sinv = np.where(S > cut, 1.0 / np.maximum(S, 1e-300), 0.0)
            return Vh.conj().T @ (Sinv * (U.conj().T @ rhs))

        def solve_step(rhs: Any) -> Any:
            nonlocal svd_used
            try:
                v = np.linalg.solve(A, rhs)
            except Exception:
                svd_used = True
                return _svd_solve(rhs)
            nv = float(np.linalg.norm(v))
            if not math.isfinite(nv) or nv > 1e8 * yn:
                svd_used = True
                return _svd_solve(rhs)
            return v

    lim = float(trust_radius) * yn if tr else float("inf")

    def cap(v: Any) -> Any:
        if tr:
            nv = float(np.linalg.norm(v))
            if math.isfinite(nv) and nv > lim > 0.0:
                return v * (lim / nv)
        return v

    extra = {"hypercube_sampling": mode, "hypercube_oracle_evals": int(oracle_evals), "slope_cond": slope_cond}
    if lm or tr:
        extra.update({"hypercube_lm_damping": float(lm_val), "hypercube_trust_radius": float(trust_radius)})
    try:
        d1 = cap(solve_step(-f))
    except Exception:
        return np.zeros_like(y), {"hypercube_error": "singular_jacobian", "hypercube_order": 0, "hypercube_nodes": M, **extra, **jac_extra}
    if svd_used:
        extra["hypercube_svd_fallback"] = True
    nd1 = float(np.linalg.norm(d1))
    if int(max_order) <= 1 or not math.isfinite(nd1) or nd1 < 1e-14 or (not tr and nd1 > 1e4 * yn):
        return d1, {"hypercube_order": 1, "hypercube_delta1_norm": nd1, "hypercube_nodes": M, **extra, **jac_extra}
    hh = max(1e-5, min(1e-2, 0.05 * yn / nd1))
    p1, m1 = target.eval(y + hh * d1), target.eval(y - hh * d1)
    p2, m2 = target.eval(y + 2.0 * hh * d1), target.eval(y - 2.0 * hh * d1)
    oracle_evals += 4
    extra["hypercube_oracle_evals"] = int(oracle_evals)
    D2 = (-p2 + 16.0 * p1 - 30.0 * f + 16.0 * m1 - m2) / (12.0 * hh**2)
    D3 = (p2 - 2.0 * p1 + 2.0 * m1 - m2) / (2.0 * hh**3)
    try:
        d2 = cap(solve_step(-0.5 * D2))
    except Exception:
        return d1, {"hypercube_order": 1, "hypercube_error": "singular_delta2", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, **extra, **jac_extra}
    nd2 = float(np.linalg.norm(d2))
    if not math.isfinite(nd2) or (not tr and nd2 > 2.0 * max(nd1, yn)):
        return d1, {"hypercube_order": 1, "hypercube_error": "delta2_rejected", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, **extra, **jac_extra}
    try:
        f2 = target.eval(y + hh * d2)
        f12 = target.eval(y + hh * d1 + hh * d2)
        oracle_evals += 2
        extra["hypercube_oracle_evals"] = int(oracle_evals)
        cross = (f12 - p1 - f2 + f) / (hh**2)
        d3 = cap(solve_step(-(cross + (1.0 / 6.0) * D3)))
    except Exception:
        return d1 + d2, {"hypercube_order": 2, "hypercube_error": "singular_delta3", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, **extra, **jac_extra}
    nd3 = float(np.linalg.norm(d3))
    if not math.isfinite(nd3) or (not tr and nd3 > 2.0 * max(nd1, yn)):
        return d1 + d2, {"hypercube_order": 2, "hypercube_error": "delta3_rejected", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, "hypercube_delta3_norm": nd3, **extra, **jac_extra}
    return cap(d1 + d2 + d3), {"hypercube_order": 4, "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, "hypercube_delta3_norm": nd3, **extra, **jac_extra}


def hypercube_inversejet_corrector(target: TargetTrack, y0: Sequence[complex], max_epochs: int, tol: float, accept: float, trial_timeout: float, line_search: int = 12, line_grid: Sequence[float] = (), direction_seed: int = 0, cloud_nodes: int = 0, lm_damping: float = 0.0, trust_radius: float = 0.0, sampling: str = "paired", known_y: Optional[Sequence[Any]] = None, defl_alpha: float = 0.0, broyden: bool = True, jet_refresh: int = 4, adaptive_lm: bool = True, parabolic: bool = True, jet_state: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    """306 hypercube inverse-jet corrector, 320 edition.

    [320-A] the local Jacobian is kept alive across epochs; after a descending
    step it is refreshed by a rank-1 Broyden secant update (0 oracle samples)
    and only rebuilt from a full paired jet every `jet_refresh` epochs or after
    a stall.  [320-B] on the base chart, line-search candidates are scored by a
    DEFLATED merit r * prod_k (1 + alpha/||y-y_k||) so descent is repelled from
    already-accepted roots -- convergence acceptance stays on the raw faithful
    residual.  [320-C] the LM damping is a live trust parameter (shrinks after
    strong descent, grows before a fresh-jet retry after a rejection).
    [320-D] a single-sample parabolic refinement follows the lambda grid.
    [320-E] the winning candidate's residual VECTOR is carried into the next
    epoch, so F(y) is never evaluated twice at the same point.
    """
    t0 = now()
    deadline = t0 + float(trial_timeout) if trial_timeout and trial_timeout > 0 else None
    y = np.asarray(y0, dtype=np.complex128).copy()
    # [320-A/E] the caller may hand over the live jet state of a previous call on
    # the SAME chart (312 calls this corrector once per outer epoch with tiny
    # inner budgets; without hand-over the Broyden recycling could never engage).
    state = jet_state if isinstance(jet_state, dict) else None
    f = None
    if state is not None and state.get("current_y") is not None and state.get("current_f") is not None and np.array_equal(np.asarray(state["current_y"]), y):
        f = np.asarray(state["current_f"], dtype=np.complex128).copy()
    if f is None:
        f = target.eval(y)
    r = finite_norm(f)
    best_y, best_f, best_r = y.copy(), f.copy(), r
    status, ok, epochs = "started", False, 0
    total_line, total_hyper, used, last = 0, 0, 0, {}
    lambdas = _line_lambdas(line_search, line_grid)
    L_full = np.asarray(lambdas, dtype=float)
    # [320-A/D] secant (Broyden) epochs use a SHORT lambda grid: a good secant step
    # wants lambda ~ 1, and the parabolic vertex probe recovers the fine-tuning the
    # long grid used to buy; a stalled short grid falls back to a full rebuild.
    L_short = L_full[: min(4, len(L_full))] if len(L_full) > 4 else L_full
    alpha = max(0.0, float(defl_alpha))
    known = [np.asarray(k, dtype=np.complex128) for k in (known_y or [])] if alpha > 0.0 else []

    A: Optional[Any] = None
    A_age = 0
    prev_y: Optional[Any] = None
    prev_f: Optional[Any] = None
    if state is not None and state.get("A") is not None:
        A = np.asarray(state["A"], dtype=np.complex128)
        A_age = int(state.get("age", 0))
        prev_y = np.asarray(state["A_y"], dtype=np.complex128) if state.get("A_y") is not None else None
        prev_f = np.asarray(state["A_f"], dtype=np.complex128) if state.get("A_f") is not None else None
    # [320-C] the adaptive trust parameter starts UNDAMPED (pure inverse-jet order);
    # damping only switches on after a rejected step.
    mu = float(state.get("mu", lm_damping)) if (adaptive_lm and state is not None) else (float(lm_damping) if float(lm_damping) > 0.0 else 0.0)
    mu = min(1.0, max(0.0, mu))
    mu_max = 1.0
    rebuilds = broyden_updates = parabolic_evals = 0
    retries_left = 2  # whole-call budget for stall retries (fresh jet / LM growth).
    refresh = max(1, int(jet_refresh))
    for ep in range(max(1, int(max_epochs))):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break
        if r < best_r:
            best_y, best_f, best_r = y.copy(), f.copy(), r
        if r <= max(float(tol), float(accept)) and (accept <= 0 or r < accept):
            ok, status = True, "converged"
            break
        rebuild = (A is None) or (prev_y is None) or (prev_f is None) or (not bool(broyden)) or (A_age >= refresh)
        if not rebuild and prev_y is not None and prev_f is not None:
            dy = y - prev_y
            df = f - prev_f
            denom = complex(np.vdot(dy, dy))
            if abs(denom) > 1e-300 and np.all(np.isfinite(df)) and np.all(np.isfinite(dy)):
                A = A + np.outer(df - A @ dy, dy.conj()) / denom
                A_age += 1
                broyden_updates += 1
                prev_y, prev_f = y.copy(), f.copy()
            elif float(np.linalg.norm(dy)) > 0.0:
                rebuild = True  # non-trivial but degenerate move: distrust the secant.
        moved = False
        for attempt in range(2):
            if rebuild:
                delta, meta = hypercube_delta(target, y, f, ep, int(direction_seed), int(cloud_nodes), float(lm_damping), float(trust_radius), str(sampling), lm_mu=(mu if (adaptive_lm and mu > 0.0) else None), return_jac=True)
                jac = meta.pop("jac", None)
                if jac is not None:
                    A, A_age = jac, 0
                    prev_y, prev_f = y.copy(), f.copy()
                rebuilds += 1
            else:
                delta, meta = hypercube_delta(target, y, f, ep, int(direction_seed), int(cloud_nodes), float(lm_damping), float(trust_radius), str(sampling), A_pre=A, lm_mu=(mu if (adaptive_lm and mu > 0.0) else None), max_order=1)
            meta.pop("jac", None)
            last = dict(meta)
            total_hyper += int(meta.get("hypercube_oracle_evals", 0))
            used += 1
            if not np.all(np.isfinite(delta)) or float(np.linalg.norm(delta)) <= 0.0:
                if not rebuild:
                    rebuild = True
                    continue
                status = "nonfinite-hypercube-step"
                break
            L = L_full if rebuild else L_short
            Y = y[None, :] + L[:, None] * delta[None, :]
            F = target.eval_batch(Y)
            R = safe_norms(F)
            total_line += int(Y.shape[0])
            idx, log_w, log_merit = deflated_line_choice(R, Y, y, r, best_r, known, alpha)
            if idx is None:
                if not rebuild and retries_left > 0:
                    # [320-A] a stale Broyden jet earns a fresh-jet retry.
                    rebuild = True
                    retries_left -= 1
                    continue
                if adaptive_lm and attempt == 0 and mu < mu_max and retries_left > 0:
                    # [320-C] grow the trust damping and retry once on a fresh jet.
                    mu = min(mu_max, max(mu * 5.0, 1e-3))
                    retries_left -= 1
                    continue
                status = "no-hypercube-decrease"
                break
            lam, cand_y, cand_f, cand_r = float(L[idx]), Y[idx], F[idx], float(R[idx])
            if parabolic and math.isfinite(cand_r):
                # [320-D] one-sample parabolic refinement: fit r^2(lambda) through the
                # current point and the best bracket, probe the vertex.
                j = idx - 1 if (idx > 0 and math.isfinite(float(R[idx - 1]))) else (idx + 1 if (idx + 1 < len(L) and math.isfinite(float(R[idx + 1]))) else None)
                if j is not None:
                    x1, q1 = 0.0, r * r
                    x2, q2 = float(L[idx]), float(R[idx]) ** 2
                    x3, q3 = float(L[j]), float(R[j]) ** 2
                    den = (x1 - x2) * (x1 - x3) * (x2 - x3)
                    if abs(den) > 1e-300:
                        a2 = (x3 * (q2 - q1) + x2 * (q1 - q3) + x1 * (q3 - q2)) / den
                        b1 = (x3 * x3 * (q1 - q2) + x2 * x2 * (q3 - q1) + x1 * x1 * (q2 - q3)) / den
                        if a2 > 0.0:
                            lam_star = -b1 / (2.0 * a2)
                            if math.isfinite(lam_star) and 0.0 < lam_star <= 1.5 * float(np.max(L)) and abs(lam_star - lam) > 0.05 * max(lam, 1e-12):
                                yc = y + lam_star * delta
                                fc = target.eval(yc)
                                rc = finite_norm(fc)
                                total_line += 1
                                parabolic_evals += 1
                                log_wc = float(deflation_log_weights(yc[None, :], known, alpha)[0]) if known else 0.0
                                log_cand = math.log(max(cand_r, 1e-300)) + float(log_w[idx])
                                if math.isfinite(rc) and math.log(max(rc, 1e-300)) + log_wc < log_cand:
                                    lam, cand_y, cand_f, cand_r = float(lam_star), yc, fc, rc
            r_before = r
            prev_y, prev_f = y, f
            y, f, r = np.asarray(cand_y, dtype=np.complex128).copy(), np.asarray(cand_f, dtype=np.complex128).copy(), cand_r
            if r < best_r:
                best_y, best_f, best_r = y.copy(), f.copy(), r
            if adaptive_lm and math.isfinite(r_before) and r <= 0.7 * r_before:
                mu = 0.0 if mu < 1e-8 else mu * 0.25
            epochs = ep + 1
            last["hypercube_chosen_lambda"] = lam
            moved = True
            break
        if not moved:
            break
    else:
        status = "max-epochs"
    if state is not None:
        state["current_y"], state["current_f"], state["mu"] = best_y.copy(), best_f.copy(), float(mu)
        if np.array_equal(best_y, y) and A is not None and prev_y is not None and prev_f is not None:
            state["A"], state["age"] = np.asarray(A, dtype=np.complex128).copy(), int(A_age)
            state["A_y"], state["A_f"] = np.asarray(prev_y, dtype=np.complex128).copy(), np.asarray(prev_f, dtype=np.complex128).copy()
        else:
            state["A"], state["A_y"], state["A_f"], state["age"] = None, None, None, 0
    final_r = finite_norm(best_f)
    if final_r <= max(float(tol), float(accept)) and (accept <= 0 or final_r < accept):
        ok, status, best_r = True, "converged", final_r
    return {"accepted": bool(ok if accept <= 0 else (math.isfinite(best_r) and best_r < accept)), "ok": bool(ok), "status": status, "epochs": int(epochs), "residual": float(best_r), "y": best_y, "seconds": float(now() - t0), "slope_cond": None, "corrector": "320-306-hypercube-inversejet-broyden", "line_search_evals": int(total_line), "line_lambdas": [float(x) for x in lambdas], "hypercube_total_evals": int(total_hyper), "hypercube_used_count": int(used), "jet_rebuilds": int(rebuilds), "broyden_updates": int(broyden_updates), "parabolic_evals": int(parabolic_evals), "deflation_active": bool(known), "lm_mu_final": float(mu), "halley_enabled": False, "halley_total_evals": 0, "halley_used_count": 0, **{k: v for k, v in last.items() if k not in ("y", "jac")}}

class NonlinearChartTarget:
    def __init__(self, base: TargetTrack, scale: complex, reciprocal: bool) -> None:
        self.base, self.scale, self.reciprocal = base, complex(scale), bool(reciprocal)

    def to_base(self, u: Any) -> Any:
        uu = np.asarray(u, dtype=np.complex128)
        if not self.reciprocal:
            return self.scale * uu
        den = self.scale * uu
        den = np.where(np.abs(den) < 1e-14, den + 1e-14, den)
        return 1.0 / den

    def from_base(self, y: Any) -> Any:
        yy = np.asarray(y, dtype=np.complex128)
        if not self.reciprocal:
            return yy / self.scale
        den = self.scale * yy
        den = np.where(np.abs(den) < 1e-14, den + 1e-14, den)
        return 1.0 / den

    def eval(self, u: Sequence[complex]) -> Any:
        return self.base.eval(self.to_base(u))

    def eval_batch(self, U: Any) -> Any:
        UU = np.asarray(U, dtype=np.complex128)
        return self.base.eval_batch(self.to_base(UU))

    def residual(self, u: Sequence[complex]) -> float:
        return finite_norm(self.eval(u))

    def residuals_batch(self, U: Any) -> Any:
        return safe_norms(self.eval_batch(U))


def complex_scale_palette(gains: Sequence[float], phases: Sequence[float], top: int) -> list[complex]:
    out: list[complex] = []
    for g in gains:
        for ph in phases:
            out.append(float(g) * phase(float(ph)))
            if len(out) >= int(top):
                return out
    return out or [1.0 + 0.0j]


def chart_candidates(base: TargetTrack, y: Any, scales: Sequence[complex], inversion: bool, top: int) -> list[tuple[float, NonlinearChartTarget, Any, dict[str, Any]]]:
    recs = []
    for s in list(scales)[:max(1, int(top))]:
        for recip in ([False, True] if inversion else [False]):
            ct = NonlinearChartTarget(base, complex(s), recip)
            u0 = ct.from_base(y)
            r = ct.residual(u0)
            e = log_energy(u0)
            if math.isfinite(r) and math.isfinite(e):
                recs.append((r * (1.0 + 0.002 * e), ct, u0, {"kind": "reciprocal" if recip else "direct", "scale": s, "score": r, "log_energy": e}))
    recs.sort(key=lambda x: x[0])
    return recs


def select_chart_shortlist(recs: list[tuple[float, Any, Any, dict[str, Any]]], top: int) -> list[tuple[float, Any, Any, dict[str, Any]]]:
    """[318-D] Keep the `top` best-scored charts but guarantee kind diversity.

    The scalar threshold theorem (direct chart attractive iff the anchor multiplier
    stays below 1, with a hard instability region below x* = (1+sqrt(3))^-3 for p=3)
    shows residual score alone can shortlist two charts of the SAME unstable kind.
    If both kinds exist among the candidates, the shortlist always contains at least
    one direct and one reciprocal chart.
    """
    if not recs:
        return []
    top = max(1, int(top))
    picked = list(recs[:top])
    kinds_all = {r[3].get("kind") for r in recs}
    kinds_picked = {r[3].get("kind") for r in picked}
    if len(kinds_all) > 1 and len(kinds_picked) == 1:
        missing = next(r for r in recs if r[3].get("kind") not in kinds_picked)
        if len(picked) >= 2:
            picked[-1] = missing
        else:
            picked.append(missing)
    return picked


def lazy_irp_hypercube_inversejet_corrector_312(base: TargetTrack, y0: Sequence[complex], max_epochs: int, tol: float, accept: float, trial_timeout: float, line_search: int = 12, line_grid: Sequence[float] = (), direction_seed: int = 0, cloud_nodes: int = 0, irp_layers: int = 2, irp_inner_epochs: int = 2, irp_scales: Optional[Sequence[complex]] = None, irp_chart_top: int = 2, irp_inversion: bool = True, collapse: bool = True, collapse_residual: float = 1e-4, collapse_drop: float = 0.42, collapse_rel_step: float = 0.35, collapse_after: int = 2, local_inner_epochs: int = 3, lazy_direct_epochs: int = 1, lazy_trigger_drop: float = 0.82, lazy_trigger_after: int = 1, lazy_bad_cond: float = 1e10, lazy_log_energy: float = 8.0, eager_irp: bool = False, rescue_collapsed: bool = False, lm_damping: float = 0.0, trust_radius: float = 0.0, sampling: str = "paired", known_y: Optional[Sequence[Any]] = None, dup_sep: float = 0.0, ghost_escapes: int = 2, ghost_kick: float = 0.45, ghost_gate: float = 1e4, defl_alpha: float = 0.0, broyden: bool = True, jet_refresh: int = 4, adaptive_lm: bool = True, parabolic: bool = True) -> dict[str, Any]:
    t0 = now()
    deadline = t0 + float(trial_timeout) if trial_timeout and trial_timeout > 0 else None
    y = np.asarray(y0, dtype=np.complex128).copy()
    best_y, best_r = y.copy(), base.residual(y)
    status, ok, epochs_done = "started", False, 0
    scales = list(irp_scales or [1.0 + 0.0j])
    total_line = total_hyper = hyper_used = triggers = rescues = direct_steps = skipped = chart_switches = recip_uses = direct_uses = 0
    collapsed, locality_hits, stagnation_hits = False, 0, 0
    collapse_epoch = collapse_reason = None
    ghost_left, ghost_kicks, dup_early = max(0, int(ghost_escapes)), 0, False
    ghost_gate_r = max(float(accept) * float(ghost_gate), 0.0)
    r_init = best_r
    known = [np.asarray(k, dtype=np.complex128) for k in (known_y or [])]
    jet_state: dict[str, Any] = {}  # [320-A] Broyden jet handed over between inner-corrector calls.
    last: dict[str, Any] = {}
    last_chart: dict[str, Any] = {"kind": "direct", "scale": 1.0 + 0j, "score": best_r, "log_energy": log_energy(y)}

    def absorb(loc: dict[str, Any]) -> None:
        nonlocal total_line, total_hyper, hyper_used
        total_line += int(loc.get("line_search_evals", 0) or 0)
        total_hyper += int(loc.get("hypercube_total_evals", 0) or 0)
        hyper_used += int(loc.get("hypercube_used_count", 0) or 0)

    for ep in range(max(1, int(max_epochs))):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break
        r0 = base.residual(y)
        if r0 < best_r:
            best_y, best_r = y.copy(), r0
        if known and float(dup_sep) > 0.0:
            dmin = min(float(np.linalg.norm(y - ky)) for ky in known)
            if dmin <= float(dup_sep):
                status, dup_early = "duplicate-early", True
                best_y, best_r = y.copy(), r0
                break
        if r0 <= max(float(tol), float(accept)) and (accept <= 0 or r0 < accept):
            ok, status = True, "converged"
            break
        remaining = 0.0 if deadline is None else deadline - now()
        if deadline is not None and remaining <= 0.0:
            status = "timeout"
            break
        loc = hypercube_inversejet_corrector(base, y, max(1, int(local_inner_epochs if collapsed else lazy_direct_epochs)), tol, accept, remaining, line_search, line_grid, int(direction_seed) + 1000003 * ep + 17, cloud_nodes, lm_damping, trust_radius, sampling, known_y=(known if float(defl_alpha) > 0.0 else None), defl_alpha=float(defl_alpha), broyden=bool(broyden), jet_refresh=int(jet_refresh), adaptive_lm=bool(adaptive_lm), parabolic=bool(parabolic), jet_state=jet_state)  # [320-A/B/C/D/E]
        absorb(loc)
        last = dict(loc)
        yd = np.asarray(loc.get("y", y), dtype=np.complex128)
        rd = base.residual(yd)
        if deadline is not None and now() >= deadline:
            status = "timeout"
            if rd < best_r:
                best_y, best_r = yd.copy(), rd
            break
        if known and float(dup_sep) > 0.0 and (accept <= 0 or not (math.isfinite(rd) and rd < accept)):
            dmin = min(float(np.linalg.norm(yd - ky)) for ky in known)
            if dmin <= float(dup_sep):
                status, dup_early = "duplicate-early", True
                best_y, best_r = yd.copy(), rd
                break
        improved = math.isfinite(rd) and rd <= r0 * (1.0 - 1e-14)
        if improved:
            direct_steps += 1
            stagnation_hits = 0
            if rd < best_r:
                best_y, best_r = yd.copy(), rd
        else:
            stagnation_hits += 1
        # [320-A] the drop threshold is judged PER inner epoch (geometric compounding),
        # so multi-epoch secant runs face the same per-step bar as 318's single steps.
        ep_used = min(4, max(1, int(loc.get("epochs", 1) or 1)))
        strong = bool(improved and rd / max(r0, 1e-300) <= float(lazy_trigger_drop) ** ep_used and log_energy(yd) <= float(lazy_log_energy))
        try:
            loc_cond = float(loc.get("slope_cond"))
            cond_trigger = bool(float(lazy_bad_cond) > 0.0 and math.isfinite(loc_cond) and loc_cond >= float(lazy_bad_cond))
        except (TypeError, ValueError):
            cond_trigger = False
        trigger = bool(eager_irp and not collapsed) or cond_trigger or (not improved and stagnation_hits >= max(1, int(lazy_trigger_after))) or (improved and not strong)
        if collapsed and not rescue_collapsed:
            trigger = False
        if not trigger:
            if improved:
                step = float(np.linalg.norm(yd - y)) / max(1.0, float(np.linalg.norm(y)))
                drop = rd / max(r0, 1e-300)
                y = yd.copy()
                direct_uses += 1
                epochs_done = ep + 1
                if collapse and not collapsed:
                    if rd <= float(collapse_residual):
                        collapsed, collapse_epoch, collapse_reason = True, ep, "residual-threshold"
                    elif drop <= float(collapse_drop) and step <= float(collapse_rel_step):
                        locality_hits += 1
                        if locality_hits >= max(1, int(collapse_after)):
                            collapsed, collapse_epoch, collapse_reason = True, ep, "observed-local-contraction"
                status = "local-collapsed" if collapsed else "direct-306-running"
                continue
            if ghost_left > 0 and (accept <= 0 or not (math.isfinite(best_r) and best_r < accept)) and math.isfinite(best_r) and (best_r <= ghost_gate_r or best_r <= 1e-3 * max(r_init, 1e-300)):
                # [318-C] stalled while unconverged: likely trapped near a ghost
                # attractor (attracting non-root cycle).  Kick the best iterate by a
                # phase rotation plus a small raw direction and resume.
                ghost_left -= 1
                ghost_kicks += 1
                kick = raw_direction(len(best_y), 9001 + ghost_kicks, int(direction_seed))
                y = best_y * (1.0 + float(ghost_kick) * phase(2.4 * ghost_kicks)) + 0.5 * float(ghost_kick) * max(1e-3, float(np.linalg.norm(best_y))) * kick
                stagnation_hits, locality_hits = 0, 0
                jet_state.clear()  # [320-A] the kick is a discontinuous jump: the secant is void.
                epochs_done = ep + 1
                status = "ghost-escape"
                continue
            status = "no-lazy-direct-decrease"
            break
        triggers += 1
        rescues += 1 if not improved else 0
        base_y, base_r = (yd.copy(), rd) if improved else (y.copy(), r0)
        irp_y, irp_r = base_y.copy(), base_r
        any_irp = False
        for layer in range(max(1, int(irp_layers))):
            if deadline is not None and now() >= deadline:
                status = "timeout"
                break
            # [318-D] score a slightly wider pool cheaply (one residual each), then
            # actually correct only a kind-diverse shortlist of size irp_chart_top
            # (317 corrected every direct AND reciprocal variant of the top scales).
            cands = chart_candidates(base, irp_y, scales, bool(irp_inversion), int(irp_chart_top) + 3)
            if not cands:
                status = "no-admissible-lazy-irp-chart"
                break
            shortlist = select_chart_shortlist(cands, int(irp_chart_top))
            best_cy = None
            best_cr = irp_r
            best_cloc: dict[str, Any] = {}
            best_chart: dict[str, Any] = {}
            for _, ct, u0, cm in shortlist:
                remaining = 0.0 if deadline is None else deadline - now()
                if deadline is not None and remaining <= 0.0:
                    status = "timeout"
                    break
                loc2 = hypercube_inversejet_corrector(ct, u0, max(1, int(irp_inner_epochs)), tol, accept, remaining, line_search, line_grid, int(direction_seed) + 1000003 * ep + 65537 * layer, cloud_nodes, lm_damping, trust_radius, sampling, broyden=bool(broyden), jet_refresh=int(jet_refresh), adaptive_lm=bool(adaptive_lm), parabolic=bool(parabolic))
                absorb(loc2)
                yb = ct.to_base(loc2.get("y", u0))
                rb = base.residual(yb)
                if math.isfinite(rb) and rb < best_cr:
                    best_cy, best_cr, best_cloc, best_chart = yb.copy(), rb, dict(loc2), dict(cm)
            if status == "timeout":
                break
            if best_cy is None:
                status = "lazy-irp-no-extra-decrease"
                break
            irp_y, irp_r, last, last_chart = best_cy, best_cr, best_cloc, best_chart
            any_irp = True
            chart_switches += 1
            recip_uses += 1 if best_chart.get("kind") == "reciprocal" else 0
            direct_uses += 0 if best_chart.get("kind") == "reciprocal" else 1
            if irp_r < best_r:
                best_y, best_r = irp_y.copy(), irp_r
            if irp_r <= max(float(tol), float(accept)) and (accept <= 0 or irp_r < accept):
                ok, status = True, "converged"
                break
        if status == "timeout":
            break
        if ok:
            break
        if any_irp and irp_r < r0:
            y = irp_y.copy()
            jet_state.clear()  # [320-A] chart round-trip: the base-chart secant is void.
            epochs_done = ep + 1
            status = "lazy-irp-running"
            continue
        if improved:
            y = yd.copy()
            epochs_done = ep + 1
            skipped += 1
            status = "direct-after-weak-irp"
            continue
        if ghost_left > 0 and (accept <= 0 or not (math.isfinite(best_r) and best_r < accept)) and math.isfinite(best_r) and (best_r <= ghost_gate_r or best_r <= 1e-3 * max(r_init, 1e-300)):
            ghost_left -= 1
            ghost_kicks += 1
            kick = raw_direction(len(best_y), 9001 + ghost_kicks, int(direction_seed))
            y = best_y * (1.0 + float(ghost_kick) * phase(2.4 * ghost_kicks)) + 0.5 * float(ghost_kick) * max(1e-3, float(np.linalg.norm(best_y))) * kick
            stagnation_hits, locality_hits = 0, 0
            jet_state.clear()  # [320-A] the kick is a discontinuous jump: the secant is void.
            epochs_done = ep + 1
            status = "ghost-escape"
            continue
        status = status if status != "started" else "no-lazy-irp-or-direct-decrease"
        break
    final_r = float(best_r)
    if final_r <= max(float(tol), float(accept)) and (accept <= 0 or final_r < accept):
        ok, status, best_r = True, "converged", final_r
    return {"accepted": bool(ok if accept <= 0 else (math.isfinite(best_r) and best_r < accept)), "ok": bool(ok), "status": status, "epochs": int(epochs_done), "residual": float(best_r), "y": best_y, "seconds": float(now() - t0), "slope_cond": None, "corrector": "320-312-lazy-irp-plus-306-hypercube-inversejet-broyden-deflated", "deflation_alpha": float(defl_alpha), "duplicate_early": bool(dup_early), "ghost_kicks": int(ghost_kicks), "ghost_escapes_left": int(ghost_left), "line_search_evals": int(total_line), "hypercube_total_evals": int(total_hyper), "hypercube_used_count": int(hyper_used), "irp_layers_completed": int(chart_switches), "irp_chart_switches": int(chart_switches), "irp_collapsed": bool(collapsed), "irp_collapse_epoch": collapse_epoch, "irp_collapse_reason": collapse_reason, "irp_reciprocal_uses": int(recip_uses), "irp_direct_uses": int(direct_uses), "irp_lazy_triggers": int(triggers), "irp_lazy_rescues": int(rescues), "irp_lazy_direct_steps": int(direct_steps), "irp_lazy_direct_good": int(direct_steps - triggers if direct_steps >= triggers else 0), "irp_lazy_direct_weak": int(triggers), "irp_lazy_skipped": int(skipped), "irp_last_chart_kind": last_chart.get("kind"), "irp_last_chart_scale": cjson(last_chart.get("scale", 1.0 + 0j)), "irp_last_chart_log_energy": last_chart.get("log_energy"), "halley_enabled": False, "halley_total_evals": 0, "halley_used_count": 0, **{k: v for k, v in last.items() if k != "y"}}


def polish_geometric_root(geometry: Any, chart: LinearChart, y: Any, accept: float, max_steps: int) -> tuple[Any, Any, float, dict[str, Any]]:
    """Newton polish on the GEOMETRY (finite-difference slope of F_geo). Residual is geometric."""
    z = chart.z_from_y(y)
    initial_f = geometry.eval(z)
    r = finite_norm(initial_f)
    meta = {"geometric_polish_enabled": bool(max_steps > 0), "geometric_polish_steps": 0, "geometric_polish_r0": float(r), "geometric_polish_r1": float(r)}
    if not math.isfinite(r):
        return z, np.asarray(y, dtype=np.complex128), r, meta
    best_z, best_f, best_r = np.asarray(z, dtype=np.complex128).copy(), np.asarray(initial_f, dtype=np.complex128).copy(), r
    for step in range(int(max_steps)):
        if best_r <= float(accept):
            break
        try:
            J = geometry.slope_matrix(best_z, best_z)
            delta, _, _, _ = np.linalg.lstsq(np.asarray(J, dtype=np.complex128), -best_f, rcond=None)
        except Exception:
            meta["geometric_polish_status"] = "linear-solve-failed"
            break
        improved = False
        for gain in (1.0, 0.5, 0.25, 0.125, 0.0625):
            cand_z = best_z + gain * delta
            cand_f = geometry.eval(cand_z)
            cand_r = finite_norm(cand_f)
            if cand_r < best_r:
                best_z, best_f, best_r = cand_z, cand_f, cand_r
                meta["geometric_polish_steps"] = step + 1
                improved = True
                break
        if not improved:
            meta["geometric_polish_status"] = "no-decrease"
            break
    meta["geometric_polish_r1"] = float(best_r)
    # [320-F] Honest conditioning estimate: SVD of a finite-difference local jet.
    # This is diagnostic, not a mathematical certificate.
    try:
        Jf = np.asarray(geometry.slope_matrix(best_z, best_z), dtype=np.complex128)
        S = np.linalg.svd(Jf, compute_uv=False)
        if S.size:
            smin, smax = float(S[-1]), float(S[0])
            meta["root_smin"] = smin
            meta["root_smax"] = smax
            meta["root_condition_method"] = "finite-difference-jet-svd-estimate"
            if smin <= 0.0:
                meta["root_cond"] = None
                meta["root_condition_status"] = "singular"
            else:
                cond = float(smax / smin)
                meta["root_cond"] = cond if math.isfinite(cond) else None
                meta["root_condition_status"] = "finite" if math.isfinite(cond) else "overflow"
            meta["root_near_multiple"] = bool(smin <= 1e-8 * max(smax, 1e-300))
            meta["root_singular"] = bool(smin <= np.finfo(float).eps * max(smax, 1e-300))
    except Exception:
        meta["root_cond"] = None
        meta["root_condition_status"] = "unavailable"
    return best_z, chart.y_from_z(best_z), best_r, meta


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    t_case = now()
    n, d = parse_case(case_raw)
    geometry, base, source, backend = make_system_317(args, n, d)
    chart = LinearChart.identity(n, float(args.linear_scale))
    target = TargetTrack(geometry, chart)  # the WHOLE pandrosion stack runs on the geometry
    user_starts = parse_start_points(args.starts, n)
    # [318-B] vectorised swarm prefilter -> priority starts (one batched Newton
    # sweep over many atlas starts before any sequential trial runs).
    swarm_meta: dict[str, Any] = {"swarm_enabled": False}
    swarm_starts_list: list[Any] = []
    if bool(args.swarm) and int(args.count) > 0 and not user_starts:
        ssize = int(args.swarm_size) if int(args.swarm_size) > 0 else min(int(args.pool), max(64, 16 * int(args.count)))
        skeep = int(args.swarm_keep) if int(args.swarm_keep) > 0 else max(4 * int(args.count), 8)
        ssep = float(args.swarm_sep) if float(args.swarm_sep) > 0 else 0.2 * math.sqrt(max(1, n))
        swarm_starts_list, swarm_meta = swarm_prefilter(
            target, n, geometry.seed + 0x320000, parse_float_list(args.rays, DEFAULT_RADII, positive=True),
            ssize, int(args.swarm_iters), skeep,
            parse_float_list(args.swarm_line, [1.0, 0.5, 0.25, 0.1], positive=True), ssep)
    # [318-F] singular-direction seeds from the origin jet.
    eigen_starts = singular_direction_starts(geometry, n, int(args.eigen_starts), float(args.origin_seed_max_norm))
    priority_starts: list[tuple[str, Any]] = (
        [("320-explicit-start", s) for s in user_starts]
        + [("320-swarm-start", s) for s in swarm_starts_list]
        + [("320-singular-direction-seed", s) for s in eigen_starts]
    )
    # The 'complete logarithmic ladder' (proportional-means / x^{k/p} ladder) is the default
    # scale geometry for 317. It gives even log coverage so a start lands in an easy basin.
    ladder_on = bool(args.scale_ladder)
    default_powers = monomial_scale_ladder(int(args.ladder_subdiv), int(args.ladder_octaves), float(args.ladder_base)) if ladder_on else DEFAULT_POWERS
    powers = sorted(set(round(float(x), 16) for x in parse_float_list(args.powers, default_powers, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in parse_float_list(args.angles, DEFAULT_ANGLES_DEG)]
    radii = parse_float_list(args.rays, DEFAULT_RADII, positive=True)
    default_startopt_gains = sorted(monomial_scale_ladder(int(args.ladder_subdiv), min(3, max(1, int(args.ladder_octaves))), float(args.ladder_base)), key=lambda g: abs(math.log(g))) if ladder_on else DEFAULT_GAINS
    gains = parse_float_list(args.startopt_gains, default_startopt_gains, positive=True)
    startopt_gains_source = "user" if args.startopt_gains not in (None, "") else ("pandrosion-ladder" if ladder_on else "legacy")
    # IRP chart palette: a tighter proportional-means ladder around 1, ordered by |log s| so the
    # near-basin scales are tried first within --irp-top, else the legacy dyadic palette.
    default_irp_gains = sorted(monomial_scale_ladder(int(args.ladder_subdiv), min(3, max(1, int(args.ladder_octaves))), float(args.ladder_base)), key=lambda g: abs(math.log(g))) if ladder_on else [1.0, 0.5, 2.0, 0.25, 4.0, 0.125, 8.0]
    irp_gain_values = parse_float_list(args.irp_gains, default_irp_gains, positive=True)
    irp_gains_source = "user" if args.irp_gains not in (None, "") else ("pandrosion-ladder" if ladder_on else "legacy")
    irp_scales = complex_scale_palette(irp_gain_values, parse_float_list(args.irp_phases, [0.0, 0.08, -0.08, 0.19, -0.19]), int(args.irp_top))
    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    failures = duplicates = 0
    atlas_calls = 0
    t_extract = now()
    for trial in range(int(args.pool)):
        if len(roots) >= int(args.count):
            break
        explicit = trial < len(priority_starts)
        if explicit:
            tag, sp = priority_starts[trial]
            y_raw = np.asarray(sp, dtype=np.complex128)
            if tag == "320-swarm-start" and roots:
                # [320-B/E] swarm survivors have already descended several Newton
                # steps, so a survivor sitting near an accepted root will almost
                # surely rediscover it: skip it for free instead of correcting it.
                skip_sep = float(args.swarm_sep) if float(args.swarm_sep) > 0 else 0.2 * math.sqrt(max(1, n))
                skip_sep = max(skip_sep, float(args.early_dup_sep))
                dmin = min(float(np.linalg.norm(y_raw - np.asarray(r["y_complex"], dtype=np.complex128))) for r in roots)
                if dmin <= skip_sep:
                    duplicates += 1
                    trials.append({"trial": int(trial), "accepted": False, "status": "swarm-start-skipped-near-root", "r1": None, "epochs": 0, "seconds": 0.0, "chart": tag, "skip_distance": dmin})
                    continue
            geom = {"chart": tag, "atlas_mode": tag, "homothety": 1.0, "theta_deg": None}
            if tag in {"320-swarm-start", "320-singular-direction-seed"}:
                # Swarm survivors are already locally optimised and basin-diverse;
                # StartOpt's gain rescaling could pull them out of their basin.
                geom["atlas_startopt_bypass_recommended"] = True
        else:
            y_raw = None
            geom: dict[str, Any] = {}
            if bool(args.origin_seed) and (trial == 0 if int(args.origin_seed_period) <= 0 else trial % max(1, int(args.origin_seed_period)) == 0):
                cap = float(args.origin_seed_max_norm) if float(args.origin_seed_max_norm) > 0 else 2.0 * math.sqrt(max(1, n))
                y_origin, ometa = origin_affine_start(target, n, float(args.origin_seed_h), cap)
                if y_origin is not None:
                    y_raw = y_origin
                    geom = {"chart": "320-origin-affine-finite-difference-seed", "atlas_mode": "origin-affine-seed", "homothety": 1.0, "theta_deg": None, **ometa}
            if y_raw is None:
                y_raw, geom = universal_atlas_start(target, n, atlas_calls, geometry.seed + 0x113000, powers, angles, radii, float(args.power_cap), len(roots), duplicates, failures, int(args.count), int(args.universal_cells), int(args.universal_shells), atlas_selection=str(args.atlas_selection))
                atlas_calls += 1
        bypass_startopt = bool(args.atlas_bypass_startopt) and bool(geom.get("atlas_startopt_bypass_recommended")) and args.startopt_gains in (None, "")
        if bypass_startopt:
            y0 = np.asarray(y_raw, dtype=np.complex128).copy()
            r0 = target.residual(y0)
            smeta = {"startopt_enabled": False, "startopt_skipped": "atlas-diversity", "startopt_r0": float(r0), "startopt_r1": float(r0), "startopt_ratio": 1.0 if math.isfinite(r0) and r0 > 0 else None, "startopt_steps": 0, "startopt_evals": 1, "startopt_micro_epochs": 0, "startopt_gain": 1.0, "startopt_batch_numpy": False}
        else:
            y0, smeta = startopt(target, y_raw, trial, geometry.seed + 0x112555, int(args.startopt_steps), int(args.startopt_candidates), gains, int(args.startopt_micro_epochs))
        known_y = [np.asarray(r["y_complex"], dtype=np.complex128) for r in roots]
        loc = lazy_irp_hypercube_inversejet_corrector_312(target, y0, int(args.epochs), float(args.tol), float(args.accept), float(args.trial_timeout), int(args.line_search), parse_float_list(args.line_grid, []), geometry.seed + 7919 * trial, int(args.hypercube_nodes), int(args.irp_layers), int(args.irp_inner_epochs), irp_scales, int(args.irp_chart_top), bool(args.irp_inversion), bool(args.collapse), float(args.collapse_residual), float(args.collapse_drop), float(args.collapse_rel_step), int(args.collapse_after), int(args.local_inner_epochs), int(args.lazy_direct_epochs), float(args.lazy_trigger_drop), int(args.lazy_trigger_after), float(args.lazy_bad_cond), float(args.lazy_log_energy), bool(args.eager_irp), bool(args.rescue_collapsed), float(args.lm_damping), float(args.trust_radius), sampling=str(args.jet_sampling), known_y=known_y, dup_sep=max(0.0, float(args.early_dup_sep)), ghost_escapes=int(args.ghost_escapes), ghost_kick=float(args.ghost_kick), ghost_gate=float(args.ghost_gate), defl_alpha=float(args.deflation_alpha), broyden=bool(args.broyden), jet_refresh=int(args.jet_refresh), adaptive_lm=bool(args.adaptive_lm), parabolic=bool(args.parabolic_line))
        if bool(loc.get("duplicate_early")):
            # [318-E] the iterate entered a known root's neighbourhood: record the
            # duplicate immediately, skip the polish, save the remaining budget.
            duplicates += 1
            z_dup = chart.z_from_y(np.asarray(loc["y"], dtype=np.complex128))
            rec = {"trial": int(trial), "accepted": False, "status": "duplicate-early", "r1": float(loc.get("residual", float("inf"))), "geometric_residual": float(loc.get("residual", float("inf"))), "epochs": int(loc.get("epochs", 0)), "seconds": float(loc.get("seconds", 0.0)), **{k: v for k, v in loc.items() if k != "y"}, **geom, **smeta}
            dup = cluster_index(roots, z_dup, max(float(args.cluster_sep), float(args.early_dup_sep)))
            if dup is not None:
                rec["cluster"] = int(dup)
            trials.append(rec)
            continue
        # Optional jet-Newton polish, still entirely on the local-jet geometry.
        z, y_final, geo_residual, polish_meta = polish_geometric_root(geometry, chart, loc["y"], float(args.accept), int(args.geometric_polish_steps))
        # The jet residual is sampled from the TRUE oracle, so it IS the faithful residual.
        accepted = bool(math.isfinite(geo_residual) and geo_residual < float(args.accept))
        rec = {"trial": int(trial), "accepted": accepted, "status": loc.get("status"), "r1": float(geo_residual), "geometric_residual": float(geo_residual), "epochs": int(loc.get("epochs", 0)), "seconds": float(loc.get("seconds", 0.0)), **{k: v for k, v in loc.items() if k != "y"}, **geom, **smeta, **polish_meta}
        if bool(args.verbose_trials):
            rec["z"] = root_to_json(z)
            rec["y0"] = root_to_json(y0)
        if not accepted:
            failures += 1
            trials.append(rec)
            continue
        dup = cluster_index(roots, z, float(args.cluster_sep))
        if dup is not None:
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = int(dup)
            trials.append(rec)
            continue
        root_cond = polish_meta.get("root_cond")
        root = {"id": len(roots), "source": "320-standalone-local-jet-geometry-fully-geometric", "trial": int(trial), "z_complex": np.asarray(z, dtype=np.complex128).copy(), "y_complex": np.asarray(y_final, dtype=np.complex128).copy(), "residual": float(geo_residual), "geometric_residual": float(geo_residual), "realness": realness(z), "cond": (float(root_cond) if root_cond is not None and math.isfinite(float(root_cond)) else None), "epochs": int(loc.get("epochs", 0)), "seconds": float(loc.get("seconds", 0.0)), **{k: v for k, v in loc.items() if k != "y"}, **geom, **smeta, **polish_meta}
        root["score"] = score_root(float(root["residual"]), float(root["realness"]), root["cond"])
        roots.append(root)
        rec["status"] = "new-root"
        rec["root_id"] = int(root["id"])
        trials.append(rec)
    encoded_roots = []
    for root in sorted(roots, key=lambda q: (float(q.get("score", float("inf"))), int(q.get("id", 0)))):
        rr = dict(root)
        rr["z"] = root_to_json(rr.pop("z_complex"))
        rr["y"] = root_to_json(rr.pop("y_complex"))
        encoded_roots.append(rr)
    gstats = geometry.stats()
    extract_seconds = float(now() - t_extract)
    samples_per_trial = (float(geometry.oracle_samples) / max(1, len(trials))) if trials else 0.0
    result = {
        "script": Path(__file__).name, "autonomous": True, "dependencies": {"python_scripts": [], "numpy": True},
        "mode": "320-standalone-local-jet-geometry-fully-geometric-swarm-lazy-irp-hypercube-inversejet",
        "flow_formula": "any system as black-box geometric oracle -> lazily-sampled local-jet field (paired central differences + [320-A] Broyden-recycled jets, O(n) samples on rebuild epochs / O(1) on secant epochs, no additional global jet precompute) -> [318-B] batched swarm prefilter -> universal Mobius atlas + singular-direction starts -> StartOpt -> 312 lazy IRP (kind-diverse charts, ghost escapes) / 306 hypercube inverse-jet ([320-B] deflated line search, [320-C] adaptive LM, [320-D] parabolic refinement, [320-E] oracle-sample reuse, SVD fallback) -> early duplicate abort -> jet-Newton polish + [320-F] finite-difference SVD conditioning estimate, ALL on the jet geometry -> faithful (jet) residual acceptance",
        "case": f"{n},{d}", "family": source, "system_source": source, "base_backend": backend, "fully_geometric": True,
        "geometry_kind": "local-jet-field", "residual_is_faithful": True,
        "seed_index": int(args.seed_index), "seed": int(geometry.seed), "n": int(n), "degree": int(d),
        "terms_per_poly": int(geometry.terms_per_poly), "terms": int(geometry.total_terms), "bezout": int(geometry.bezout),
        "equation_normalize": bool(args.equation_normalize),
        "linear_A": [[cjson(chart.A[i, j]) for j in range(n)] for i in range(n)],
        "geometry": {"kind": "local-jet-field", "use_quadratic": bool(args.jet_quadratic), "jet_radius": float(args.jet_radius), "jet_cache": bool(args.jet_cache), "samples_per_jet": int(geometry.samples_per_jet), "jets_built": int(gstats.get("jets_built", 0)), "jet_cache_hits": int(gstats.get("jet_cache_hits", 0)), "oracle_samples_total": int(geometry.oracle_samples), "construction_complexity": "no additional global jet precompute; O(n) oracle samples on paired-jet rebuilds, without monomial enumeration"},
        "parameters": {"system_source": str(args.system_source), "polys": str(args.polys or ""), "variables": str(args.variables or ""), "jet_quadratic": bool(args.jet_quadratic), "jet_radius": float(args.jet_radius), "jet_cache": bool(args.jet_cache), "jet_cache_cap": int(args.jet_cache_cap), "jet_sampling": str(args.jet_sampling), "swarm": bool(args.swarm), "swarm_iters": int(args.swarm_iters), "ghost_escapes": int(args.ghost_escapes), "ghost_kick": float(args.ghost_kick), "early_dup_sep": float(args.early_dup_sep), "eigen_starts": int(args.eigen_starts), "deflation_alpha": float(args.deflation_alpha), "broyden": bool(args.broyden), "jet_refresh": int(args.jet_refresh), "adaptive_lm": bool(args.adaptive_lm), "parabolic_line": bool(args.parabolic_line), "geometric_polish_steps": int(args.geometric_polish_steps), "starts": str(args.starts or ""), "system_mode": str(args.system_mode), "count": int(args.count), "pool": int(args.pool), "accept": float(args.accept), "tol": float(args.tol), "epochs": int(args.epochs), "cluster_sep": float(args.cluster_sep), "line_search": int(args.line_search), "hypercube_nodes": int(args.hypercube_nodes), "atlas_selection": str(args.atlas_selection), "atlas_bypass_startopt": bool(args.atlas_bypass_startopt), "atlas_calls": int(atlas_calls), "startopt_steps": int(args.startopt_steps), "startopt_candidates": int(args.startopt_candidates), "startopt_gains_source": startopt_gains_source, "startopt_gains_count": len(gains), "scale_ladder": bool(args.scale_ladder), "ladder_subdiv": int(args.ladder_subdiv), "ladder_octaves": int(args.ladder_octaves), "ladder_base": float(args.ladder_base), "powers_count": len(powers), "irp_gains_source": irp_gains_source, "irp_gains_count": len(irp_gain_values), "irp_scales_count": len(irp_scales)},
        "swarm": swarm_meta,
        "roots": encoded_roots,
        "trials": trials if bool(args.verbose_trials) else trials[: min(len(trials), int(args.keep_trials))],
        "summary": {"requested_roots": int(args.count), "unique_roots": len(roots), "success": bool(len(roots) >= int(args.count)), "trials_used": len(trials), "duplicates": int(duplicates), "duplicates_early": int(sum(1 for t in trials if t.get("status") == "duplicate-early")), "ghost_kicks_total": int(sum(int(t.get("ghost_kicks", 0) or 0) for t in trials)), "failures": int(failures), "priority_starts": len(priority_starts), "swarm_kept": int(swarm_meta.get("swarm_kept", 0) or 0), "eigen_starts": len(eigen_starts), "generation_seconds": float(geometry.generation_seconds), "extract_seconds": extract_seconds, "total_seconds": float(now() - t_case), "oracle_samples_total": int(geometry.oracle_samples), "oracle_samples_per_trial": samples_per_trial, "eval_stats": gstats},
    }
    return result


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="320 standalone NumPy local-jet-geometry fully-geometric swarm + lazy IRP / hypercube inverse-jet engine.")
    p.add_argument("--cases", default="2,4")
    p.add_argument("--seed-index", type=int, default=0)
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--no-equation-normalize", dest="equation_normalize", action="store_false")
    p.add_argument("--system-mode", choices=["auto", "dense", "geometry-kernel", "geometry", "kernel", "projective-kernel", "lazy-feature", "lazy", "feature", "stream"], default="auto")
    p.add_argument("--dense-max-terms", type=int, default=250000)
    # KS-base geometry-kernel backend knobs (only used when the base itself is a large KS system).
    p.add_argument("--geometry-anchors", type=int, default=0)
    p.add_argument("--geometry-anchor-cap", type=int, default=4096)
    p.add_argument("--geometry-anchor-scales", default="0.25,0.5,1,2,4")
    p.add_argument("--geometry-dynamic-normalize", action="store_true", default=True)
    p.add_argument("--no-geometry-dynamic-normalize", dest="geometry_dynamic_normalize", action="store_false")
    p.add_argument("--geometry-self-normalize", action="store_true", default=True)
    p.add_argument("--no-geometry-self-normalize", dest="geometry_self_normalize", action="store_false")
    p.add_argument("--geometry-eval-block", type=int, default=128)
    p.add_argument("--lazy-features", type=int, default=0)
    p.add_argument("--lazy-feature-cap", type=int, default=8192)
    p.add_argument("--lazy-projective-normalize", action="store_true", default=True)
    p.add_argument("--no-lazy-projective-normalize", dest="lazy_projective_normalize", action="store_false")
    p.add_argument("--lazy-dynamic-normalize", action="store_true", default=True)
    p.add_argument("--no-lazy-dynamic-normalize", dest="lazy_dynamic_normalize", action="store_false")
    p.add_argument("--lazy-eval-block", type=int, default=128)
    # The local-jet geometry (the heart of 317).
    p.add_argument("--jet-quadratic", action="store_true", default=True, help="Also sample local curvature (second-order jet).")
    p.add_argument("--no-jet-quadratic", dest="jet_quadratic", action="store_false")
    p.add_argument("--jet-radius", type=float, default=1e-5, help="Relative hypercube radius for local-jet sampling.")
    p.add_argument("--jet-cache", action="store_true", default=True, help="Reuse local jets across nearby queries.")
    p.add_argument("--no-jet-cache", dest="jet_cache", action="store_false")
    p.add_argument("--jet-cache-decimals", type=int, default=9)
    p.add_argument("--jet-cache-cap", type=int, default=4096, help="[318-F] LRU bound on the jet cache (entries).")
    p.add_argument("--jet-sampling", choices=["paired", "cloud"], default="paired", help="[318-A] paired = central differences O(h^2) (2n+1 samples); cloud = legacy 317 random-sign lstsq.")
    # [318-B] vectorised swarm prefilter.
    p.add_argument("--swarm", action="store_true", default=True, help="[318-B] Run a batched damped-Newton prefilter over many atlas starts before sequential trials.")
    p.add_argument("--no-swarm", dest="swarm", action="store_false")
    p.add_argument("--swarm-size", type=int, default=0, help="Batch size of the swarm (0 = auto: min(pool, max(64, 16*count))).")
    p.add_argument("--swarm-iters", type=int, default=3, help="Batched Newton iterations in the swarm prefilter.")
    p.add_argument("--swarm-keep", type=int, default=0, help="Survivors promoted to priority starts (0 = auto: max(4*count, 8)).")
    p.add_argument("--swarm-line", default="1,0.5,0.25,0.1", help="Damping grid of the batched swarm line search.")
    p.add_argument("--swarm-sep", type=float, default=0.0, help="Pairwise diversity separation between kept survivors (0 = auto: 0.2*sqrt(n)).")
    # [318-C] ghost-trap escapes.
    p.add_argument("--ghost-escapes", type=int, default=2, help="[318-C] Phase/direction kicks attempted on the best iterate when a trial stalls unconverged (attracting non-root cycles exist).")
    p.add_argument("--ghost-kick", type=float, default=0.45, help="Relative amplitude of the ghost-escape kick.")
    p.add_argument("--ghost-gate", type=float, default=1e4, help="[318-C] Only spend ghost escapes when the best residual is below ghost_gate*accept (or fell 1000x from the start): hopeless trials fail fast instead.")
    # [318-E] early duplicate abort.
    p.add_argument("--early-dup-sep", type=float, default=1e-4, help="[318-E] Abort a trial whose iterate comes within this distance of an accepted root (0 = off).")
    # [320-A] Broyden-recycled jets.
    p.add_argument("--broyden", dest="broyden", action="store_true", default=True, help="[320-A] Recycle the local Jacobian across epochs with rank-1 Broyden secant updates; full paired jets only every --jet-refresh epochs or after a stall.")
    p.add_argument("--no-broyden", dest="broyden", action="store_false", help="Rebuild the full paired jet at every epoch (318 behaviour).")
    p.add_argument("--jet-refresh", type=int, default=4, help="[320-A] Epochs between full paired-jet rebuilds when --broyden is on.")
    # [320-B] deflated line search.
    p.add_argument("--deflation-alpha", type=float, default=0.15, help="[320-B] Root-repulsion radius of the deflated line-search merit r*prod(1+alpha/||y-y_k||) on the base chart (0 = off). Acceptance stays on the raw residual.")
    # [320-C] adaptive Levenberg-Marquardt.
    p.add_argument("--adaptive-lm", dest="adaptive_lm", action="store_true", default=True, help="[320-C] Live LM trust parameter: shrinks after strong descent, grows before a fresh-jet retry after a rejection.")
    p.add_argument("--no-adaptive-lm", dest="adaptive_lm", action="store_false", help="Static --lm-damping only (318 behaviour).")
    # [320-D] parabolic line refinement.
    p.add_argument("--parabolic-line", dest="parabolic_line", action="store_true", default=True, help="[320-D] Spend one extra oracle sample on the parabolic vertex of the line search.")
    p.add_argument("--no-parabolic-line", dest="parabolic_line", action="store_false")
    # [318-F] singular-direction seeds.
    p.add_argument("--eigen-starts", type=int, default=4, help="[318-F] Priority starts along the smallest singular directions of the origin jet (0 = off).")
    p.add_argument("--geometric-polish-steps", type=int, default=6)
    p.add_argument("--origin-seed", action="store_true", default=True)
    p.add_argument("--no-origin-seed", dest="origin_seed", action="store_false")
    p.add_argument("--origin-seed-h", type=float, default=1e-5)
    p.add_argument("--origin-seed-max-norm", type=float, default=0.0)
    p.add_argument("--origin-seed-period", type=int, default=0)
    p.add_argument("--linear-scale", type=float, default=1.0)
    p.add_argument("--count", type=int, default=8)
    p.add_argument("--pool", type=int, default=4096)
    p.add_argument("--epochs", type=int, default=24)
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--accept", "--residual-accept", type=float, default=1e-8)
    p.add_argument("--cluster-sep", type=float, default=1e-8)
    p.add_argument("--trial-timeout", type=float, default=0.0)
    p.add_argument("--line-search", type=int, default=12)
    p.add_argument("--hypercube-nodes", type=int, default=0)
    p.add_argument("--irp-layers", type=int, default=2)
    p.add_argument("--irp-inner-epochs", type=int, default=2)
    p.add_argument("--local-inner-epochs", type=int, default=3)
    p.add_argument("--lazy-direct-epochs", type=int, default=2, help="[320-A] Inner direct-corrector epochs per outer 312 epoch (320 default 2: the Broyden secant needs a compound drop before the IRP trigger judges it).")
    p.add_argument("--lazy-trigger-drop", type=float, default=0.82)
    p.add_argument("--lazy-trigger-after", type=int, default=1)
    p.add_argument("--lazy-bad-cond", type=float, default=1e10)
    p.add_argument("--lazy-log-energy", type=float, default=8.0)
    p.add_argument("--eager-irp", action="store_true", default=False)
    p.add_argument("--rescue-collapsed", action="store_true", default=False)
    # Regularised inverse-jet inversion (tames the giant near-singular step in high n).
    p.add_argument("--lm-damping", type=float, default=0.0, help="Levenberg-Marquardt damping on the inverse-jet inversion (0 = off, try 1e-2).")
    p.add_argument("--trust-radius", type=float, default=0.0, help="Cap the inverse-jet step to trust_radius*||y|| (0 = off, try 1.0).")
    p.add_argument("--trust-region", action="store_true", default=False, help="Convenience: enable LM damping + trust radius with sensible defaults.")
    p.add_argument("--irp-chart-top", type=int, default=2)
    p.add_argument("--irp-gains", default=None, help="Homothety gains for the IRP chart palette. Defaults to the Pandrosion ladder when --scale-ladder is enabled, otherwise legacy dyadic gains.")
    p.add_argument("--irp-phases", default="0,0.08,-0.08,0.19,-0.19")
    p.add_argument("--irp-top", type=int, default=14)
    p.add_argument("--irp-inversion", action="store_true", default=True)
    p.add_argument("--no-irp-inversion", dest="irp_inversion", action="store_false")
    p.add_argument("--collapse", action="store_true", default=True)
    p.add_argument("--no-collapse", dest="collapse", action="store_false")
    p.add_argument("--collapse-residual", type=float, default=1e-4)
    p.add_argument("--collapse-drop", type=float, default=0.42)
    p.add_argument("--collapse-rel-step", type=float, default=0.35)
    p.add_argument("--collapse-after", type=int, default=2)
    p.add_argument("--line-grid", default="1,0.75,0.5,0.35,0.25,0.18,0.125,0.09,0.0625,0.045,0.03125,0.02")
    p.add_argument("--powers", default=None)
    p.add_argument("--power-cap", type=float, default=1048576.0)
    # The 'complete logarithmic ladder' (Pandrosion proportional-means x^{k/p} scaling).
    p.add_argument("--scale-ladder", dest="scale_ladder", action="store_true", default=True, help="Use an equally-spaced logarithmic (proportional-means) ladder for the start homothety palette and IRP charts.")
    p.add_argument("--no-scale-ladder", dest="scale_ladder", action="store_false", help="Use the legacy dyadic/decadic start palette and dyadic IRP chart gains.")
    p.add_argument("--ladder-subdiv", type=int, default=3, help="Geometric means per octave (p): 3 inserts the cube-root means base^{k/3}; the x^{k/p} ladder.")
    p.add_argument("--ladder-octaves", type=int, default=12, help="Half-span of the start ladder in octaves of --ladder-base.")
    p.add_argument("--ladder-base", type=float, default=2.0, help="Octave base of the logarithmic ladder (>1).")
    p.add_argument("--angles", default=None)
    p.add_argument("--rays", default=None)
    p.add_argument("--startopt-steps", type=int, default=1)
    p.add_argument("--startopt-candidates", type=int, default=12)
    p.add_argument("--startopt-gains", default=None)
    p.add_argument("--startopt-micro-epochs", type=int, default=0)
    p.add_argument("--universal-cells", type=int, default=16)
    p.add_argument("--universal-shells", type=int, default=5)
    p.add_argument("--atlas-selection", choices=["diverse-shell", "compact-score"], default="diverse-shell", help="Automatic start selection. diverse-shell cycles deterministic logarithmic shells; compact-score keeps the older residual-min atlas cell selection.")
    p.add_argument("--atlas-bypass-startopt", action="store_true", default=True, help="Do not run StartOpt on diverse-shell atlas starts, preserving basin diversity for all-root extraction.")
    p.add_argument("--no-atlas-bypass-startopt", dest="atlas_bypass_startopt", action="store_false")
    p.add_argument("--out", default=None)
    p.add_argument("--outdir", default="/tmp/320_standalone_out")
    p.add_argument("--keep-trials", type=int, default=160)
    p.add_argument("--verbose-trials", action="store_true")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--system-source", choices=["ks", "kostlan", "polynomial", "poly", "expr", "expression"], default="ks")
    p.add_argument("--polys", "--poly", default=None)
    p.add_argument("--variables", default=None)
    p.add_argument("--starts", default=None)
    return p


def validate_args(args: argparse.Namespace) -> None:
    for case_raw in [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]:
        parse_case(case_raw)
    if int(args.count) < 0 or int(args.pool) < 0:
        raise ValueError("--count and --pool must be non-negative")
    if int(args.epochs) <= 0 or int(args.line_search) <= 0:
        raise ValueError("--epochs and --line-search must be positive")
    if not math.isfinite(float(args.accept)) or float(args.accept) <= 0.0:
        raise ValueError("--accept must be finite and positive")
    if not math.isfinite(float(args.tol)) or float(args.tol) < 0.0:
        raise ValueError("--tol must be finite and non-negative")
    if not math.isfinite(float(args.linear_scale)) or abs(float(args.linear_scale)) <= 1e-300:
        raise ValueError("--linear-scale must be finite and non-zero")
    if not math.isfinite(float(args.trial_timeout)) or float(args.trial_timeout) < 0.0:
        raise ValueError("--trial-timeout must be finite and non-negative")
    if not math.isfinite(float(args.cluster_sep)) or float(args.cluster_sep) <= 0.0:
        raise ValueError("--cluster-sep must be finite and positive")
    if float(args.early_dup_sep) < 0.0 or float(args.deflation_alpha) < 0.0:
        raise ValueError("duplicate separation and deflation alpha must be non-negative")
    if bool(args.scale_ladder) and (not math.isfinite(float(args.ladder_base)) or float(args.ladder_base) <= 1.0):
        raise ValueError("--ladder-base must be finite and > 1")
    if not parse_float_list(args.line_grid, [], positive=True):
        raise ValueError("--line-grid must contain at least one positive finite value")


def run_self_test_suite_320(args: argparse.Namespace) -> dict[str, Any]:
    """Deterministic regressions for extraction, deflation and singular roots."""
    specs = [
        ("quadratic-two-roots", "1,2", "x^2 - 3*x - 10", "-8,4", 2, 6),
        ("product-four-roots", "2,2", "x1^2 - 1; x2^2 - 1", "-2,-2;-2,2;2,-2;2,2", 4, 6),
        # Polishing is disabled on purpose: conditioning must still be estimated.
        ("multiple-root-conditioning", "1,2", "x^2", "0", 1, 0),
    ]
    records: list[dict[str, Any]] = []
    for name, case_raw, polys, starts, count, polish_steps in specs:
        a = argparse.Namespace(**vars(args))
        a.self_test = False
        a.system_source = "polynomial"
        a.system_mode = "auto"
        a.cases = case_raw
        a.polys = polys
        a.starts = starts
        a.variables = None
        a.count = count
        a.pool = max(count, 8)
        a.epochs = 24
        a.accept = min(float(a.accept), 1e-8)
        a.swarm = False
        a.eigen_starts = 0
        a.origin_seed = False
        a.startopt_steps = 0
        a.geometric_polish_steps = polish_steps
        result = run_case(a, case_raw)
        passed = bool(result["summary"]["success"] and len(result["roots"]) == count)
        if name == "multiple-root-conditioning":
            root = result["roots"][0] if result["roots"] else {}
            passed = bool(passed and root.get("root_near_multiple") and root.get("root_condition_status") == "singular" and root.get("cond") is None)
        records.append({"name": name, "passed": passed, "result": result})

    class _QuadraticTarget:
        def eval(self, y: Sequence[complex]) -> Any:
            x = np.asarray(y, dtype=np.complex128)[0]
            return np.asarray([x * x - 1.0], dtype=np.complex128)

        def eval_batch(self, Y: Any) -> Any:
            return np.vstack([self.eval(y) for y in np.asarray(Y, dtype=np.complex128)])

        def residual(self, y: Sequence[complex]) -> float:
            return finite_norm(self.eval(y))

    state: dict[str, Any] = {}
    qtarget = _QuadraticTarget()
    first = hypercube_inversejet_corrector(qtarget, [2.0 + 0.0j], 1, 1e-14, 1e-14, 0.0, line_grid=[1.0, 0.5, 0.25], parabolic=False, jet_state=state)
    second = hypercube_inversejet_corrector(qtarget, first["y"], 1, 1e-14, 1e-14, 0.0, line_grid=[1.0, 0.5, 0.25], parabolic=False, jet_state=state)
    broyden_passed = bool(second.get("broyden_updates", 0) >= 1 and second.get("hypercube_sampling") == "broyden-reuse")
    records.append({"name": "cross-call-anchored-broyden", "passed": broyden_passed, "first": first, "second": second})
    defl_idx, _, _ = deflated_line_choice(
        np.asarray([11.0, 9.0]), np.asarray([[1.0 + 0.0j], [0.5 + 0.0j]]),
        np.asarray([0.0 + 0.0j]), 10.0, 10.0, [np.asarray([0.5 + 0.0j])], 1.0,
    )
    records.append({"name": "acceptance-safe-deflation", "passed": bool(defl_idx == 1), "selected_index": defl_idx})
    return {"script": Path(__file__).name, "self_test": True, "passed": bool(all(r["passed"] for r in records)), "checks": records}


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if bool(args.trust_region):
        if float(args.lm_damping) <= 0.0:
            args.lm_damping = 1e-2
        if float(args.trust_radius) <= 0.0:
            args.trust_radius = 1.0
    validate_args(args)
    if bool(args.self_test):
        final = run_self_test_suite_320(args)
        out = Path(args.out) if args.out else Path(args.outdir) / "self_test_320.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(json_safe(final), indent=2, allow_nan=False), encoding="utf-8")
        print(f"320 self-test: {'PASS' if final['passed'] else 'FAIL'} ({sum(int(c['passed']) for c in final['checks'])}/{len(final['checks'])})", flush=True)
        print(f"out={out}", flush=True)
        return 0 if final["passed"] else 1
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    outputs = [run_case(args, c) for c in cases]
    final = outputs[0] if len(outputs) == 1 else {"script": Path(__file__).name, "autonomous": True, "cases": outputs}
    out = Path(args.out) if args.out else Path(args.outdir) / f"320_standalone_{cases[0].replace(',', 'x') if cases else 'case'}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(json_safe(final), indent=2, allow_nan=False), encoding="utf-8")
    print("=" * 120, flush=True)
    print("320 STANDALONE LOCAL-JET GEOMETRY -- FULLY GEOMETRIC SWARM + LAZY IRP / HYPERCUBE INVERSE-JET", flush=True)
    print("NumPy + stdlib only; no local flow imports. Broyden-recycled jets, deflated (root-repelling) line search,", flush=True)
    print("adaptive LM trust, parabolic line refinement, oracle-sample reuse, honest SVD conditioning estimate; residual is faithful (true oracle).", flush=True)
    print("=" * 120, flush=True)
    for r in outputs:
        s = r["summary"]
        g = r.get("geometry", {})
        sw = r.get("swarm", {})
        print(f"case={r.get('family')}({r['n']},{r['degree']}), base_backend={r.get('base_backend')}, geometry=local-jet-field (quad={g.get('use_quadratic')}, samples/jet={g.get('samples_per_jet')}), seed={r['seed']}", flush=True)
        print(f"roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} trials={s['trials_used']} duplicates={s['duplicates']} (early={s.get('duplicates_early', 0)}) failures={s['failures']} ghost_kicks={s.get('ghost_kicks_total', 0)}", flush=True)
        if sw.get("swarm_enabled"):
            print(f"swarm: size={sw.get('swarm_size')} iters={sw.get('swarm_iters')} survivors={sw.get('swarm_survivors')} kept={sw.get('swarm_kept')} full_jets={sw.get('swarm_full_jets')} broyden={sw.get('swarm_broyden_updates')} best_r={sw.get('swarm_best_residual'):.3e} ({sw.get('swarm_seconds', 0.0):.2f}s)", flush=True)
        print(f"seconds: extract={s['extract_seconds']:.2f}, total={s['total_seconds']:.2f}; oracle_samples={s['oracle_samples_total']} ({s['oracle_samples_per_trial']:.0f}/trial)", flush=True)
        if r.get("roots"):
            best = r["roots"][0]
            print(f"best_root: residual={float(best.get('geometric_residual', float('inf'))):.3e} (faithful), trial={best.get('trial')}", flush=True)
    print(f"out={out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
