#!/usr/bin/env python3
"""
305_pandrosion_highdim_universal_atlas_full_cubic_halley_numpy_engine.py

Autonomous ANCHORED FULL-CUBIC HALLEY Pandrosion root extractor.

301 adapts articles/001, Paper 0 bis: Cubic Pandrosion by Finite-Slope
Halley Acceleration, to a strict multivariate finite-slope cubic corrector.  It still uses the exact multivariate anchored
finite-slope matrix Q_G(a,b), but the accepted correction is never the raw
Steffensen/Pandrosion step.  Each local step is the full cubic Halley-style
finite-slope direction

    Q delta1 = -G(a)
    Q delta2 = -1/2 H_G(delta1,delta1)
    delta = delta1 + delta2

with H_G(delta1,delta1) estimated by symmetric finite probes.  Line search is
performed only along the cubic direction.  No analytic Jacobian, Hessian,
Newton fallback, Broyden update, homotopy path, SciPy, or imports from previous
flow scripts are used.

No imports from previous Pandrosion scripts.  This file contains its own:
  - Kostlan/dense polynomial generator for ks(n,d)
  - dense F evaluation engine + vectorized exact telescopic slope Q(a,b)
  - single geometric flow: heuristic Riemann/Mobius + unconditional Thales homothety
  - starting point optimization inspired by s0^opt = h(1)
  - local PURE Pandrosion corrector using vectorized exact monomial telescopic slope Q(a,b)
  - root clustering and JSON export

The goal is not to enumerate all Bezout paths.  It extracts useful complex roots
from high-degree systems using direct geometric starts.

Core flow per trial:

    u
      -> Riemann/Mobius chart with angle theta and pole p
      -> unconditional homothety Lambda
      -> startopt radial/geometric improvement
      -> y
      -> z = A y
      -> full cubic finite-slope Halley correction using Q_G(a,b)
      -> validate F(z) and deduplicate

Dependencies: Python stdlib + NumPy.  There is intentionally no dependency on
any other .py file.
"""
from __future__ import annotations

import argparse
import cmath
import dataclasses
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Iterable, List, Optional, Sequence, Tuple

def _bootstrap_numpy_path() -> None:
    """Make NumPy visible even when launched with python -S.

    This remains autonomous: it does not import any previous Pandrosion script; it
    only exposes normal site-packages directories that python -S hides.
    """
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
# Generic helpers
# ---------------------------------------------------------------------------

def now() -> float:
    return time.time()


def cjson(z: complex) -> list[float]:
    return [float(complex(z).real), float(complex(z).imag)]


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


def norm2(v: Sequence[complex]) -> float:
    return math.sqrt(sum(abs(complex(x)) ** 2 for x in v))


def ensure_numpy() -> None:
    if np is None:
        raise RuntimeError(
            "NumPy is required by this autonomous high-degree extractor. "
            f"Import error: {_NUMPY_IMPORT_ERROR!r}"
        )


# ---------------------------------------------------------------------------
# Polynomial system generation/evaluation
# ---------------------------------------------------------------------------

def compositions_leq(d: int, n: int) -> "np.ndarray":
    """All alpha in N^n with |alpha| <= d, shape (M,n), lexicographic-ish."""
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
    """sqrt(d! / ((d-|a|)! prod a_j!)) for total degree <= d."""
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


def stable_seed(n: int, d: int, seed_index: int = 0, salt: int = 0) -> int:
    return int(splitmix64(0x50414E44524F5349 + 1000003 * n + 9176 * d + 97 * seed_index + salt) & 0x7FFFFFFF)


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
        ensure_numpy()
        t0 = now()
        exps = compositions_leq(d, n)
        weights = multinomial_kostlan_weights(exps, d)
        seed = stable_seed(n, d, seed_index)
        rng = np.random.default_rng(seed)
        # Complex Gaussian Kostlan coefficients, normalized by sqrt(2) so E|a|^2=1 before weights.
        coeff = (rng.standard_normal((n, exps.shape[0])) + 1j * rng.standard_normal((n, exps.shape[0]))) / math.sqrt(2.0)
        coeff = coeff * weights[None, :]
        if equation_normalize:
            row_norm = np.linalg.norm(coeff, axis=1)
            row_norm = np.where(row_norm > 0, row_norm, 1.0)
            coeff = coeff / row_norm[:, None]
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
        # pow_tables[j][k] = z_j^k, k=0..d
        tables = []
        for zj in z:
            p = np.empty(self.d + 1, dtype=np.complex128)
            p[0] = 1.0 + 0.0j
            if self.d > 0:
                p[1] = zj
                for k in range(2, self.d + 1):
                    p[k] = p[k-1] * zj
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
        """Batched monomial table for row-vector points.

        Z has shape (B,n).  The return has shape (B,M), where M is the number
        of monomials.  This is used only to score geometric probes and
        line-search candidates; it does not change the finite-slope corrector.
        """
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
            if self.d > 0:
                p[:, 1] = ZZ[:, j]
                for k in range(2, self.d + 1):
                    p[:, k] = p[:, k - 1] * ZZ[:, j]
            mon *= p[:, self.exps[:, j]]
        return mon

    def eval_batch(self, Z: "np.ndarray") -> "np.ndarray":
        """Evaluate many points in one BLAS-backed block."""
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

    def _slope_power_table(self, aj: complex, bj: complex, pows_a_j: "np.ndarray", pows_b_j: "np.ndarray") -> "np.ndarray":
        """Vectorized pure Pandrosion power-slope table.

        slope[m] = sum_{r=0}^{m-1} bj^(m-1-r) aj^r, with slope[0]=0.

        This is the exact finite Thales/Pandrosion factor for z^m between
        aj and bj.  It uses the O(d) polynomial recurrence

            S_0 = 0,   S_m = bj^(m-1) + aj*S_(m-1)

        so it is not a derivative, not a divided-difference shortcut, and not
        a fallback algorithm.
        """
        slope = np.empty(self.d + 1, dtype=np.complex128)
        slope[0] = 0.0 + 0.0j
        acc = 0.0 + 0.0j
        ajc = complex(aj)
        for m in range(1, self.d + 1):
            acc = pows_b_j[m - 1] + ajc * acc
            slope[m] = acc
        return slope

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> "np.ndarray":
        """Vectorized exact derivative-free Pandrosion slope matrix Q(a,b).

        For every polynomial row F_i and two points a,b, this constructs Q so
        that

            F(b) - F(a) = Q(a,b) @ (b-a)

        using the telescopic monomial identity.  No Jacobian is formed and no
        derivative formula is called.

        Vectorization in 120:
          * index all monomial factors a_k^alpha_k and b_k^alpha_k once;
          * build prefix_b[j] = prod_{k<j} b_k^alpha_k as whole vectors;
          * build suffix_a[j] = prod_{k>=j} a_k^alpha_k as whole vectors;
          * build the coordinate slope table S_j[m] by the pure finite-sum
            recurrence, then index it as S_j[alpha_j];
          * compute every Q column as one BLAS-backed coeff @ term product.

        This preserves 115's pure Pandrosion formula while removing the slow
        Python loops over exponent masks and finite sums.  It is generic in n;
        there is no bivariate specialization.
        """
        ensure_numpy()
        t0 = now()
        aa = np.asarray(a, dtype=np.complex128)
        bb = np.asarray(b, dtype=np.complex128)
        pows_a = self._powers(aa)
        pows_b = self._powers(bb)
        M = self.terms_per_poly
        n = self.n

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

        # Materialize all Q-column monomial factors as one (M,n) block, then
        # use a single BLAS GEMM: Q = coeff @ terms.  This is faster than n
        # separate matrix-vector products and remains fully general in n.
        terms = np.empty((M, n), dtype=np.complex128)
        for j in range(n):
            slope_table = self._slope_power_table(aa[j], bb[j], pows_a[j], pows_b[j])
            terms[:, j] = prefix_b[j] * suffix_a[j + 1] * slope_table[self.exps[:, j]]
        Q = self.coeff @ terms

        self.slope_count = int(getattr(self, "slope_count", 0)) + 1
        self.seconds_slope = float(getattr(self, "seconds_slope", 0.0)) + (now() - t0)
        return Q

    def stats(self) -> dict[str, Any]:
        return {
            "eval_count": int(self.eval_count),
            "slope_count": int(getattr(self, "slope_count", 0)),
            "seconds_eval": float(self.seconds_eval),
            "seconds_slope": float(getattr(self, "seconds_slope", 0.0)),
            "terms_per_poly": self.terms_per_poly,
            "total_terms": self.total_terms,
        }


@dataclasses.dataclass
class LinearChart:
    A: Any
    Ainv: Any

    @classmethod
    def identity(cls, n: int, scale: float = 1.0) -> "LinearChart":
        ensure_numpy()
        A = np.eye(n, dtype=np.complex128) * complex(scale)
        Ainv = np.eye(n, dtype=np.complex128) / complex(scale)
        return cls(A=A, Ainv=Ainv)

    def z_from_y(self, y: Sequence[complex]) -> "np.ndarray":
        return self.A @ np.asarray(y, dtype=np.complex128)

    def y_from_z(self, z: Sequence[complex]) -> "np.ndarray":
        return self.Ainv @ np.asarray(z, dtype=np.complex128)


@dataclasses.dataclass
class TargetTrack:
    system: DenseKostlanSystem
    chart: LinearChart

    def eval(self, y: Sequence[complex]) -> "np.ndarray":
        z = self.chart.z_from_y(y)
        return self.system.eval(z)

    def eval_batch(self, Y: "np.ndarray") -> "np.ndarray":
        YY = np.asarray(Y, dtype=np.complex128)
        if YY.ndim == 1:
            return self.eval(YY)[None, :]
        # z = A @ y for column-vector y; for row batches this is Y @ A.T.
        Z = YY @ np.asarray(self.chart.A, dtype=np.complex128).T
        return self.system.eval_batch(Z)

    def slope_matrix(self, a_y: Sequence[complex], b_y: Sequence[complex]) -> "np.ndarray":
        a_z = self.chart.z_from_y(a_y)
        b_z = self.chart.z_from_y(b_y)
        return self.system.slope_matrix(a_z, b_z) @ self.chart.A

    def residual(self, y: Sequence[complex]) -> float:
        try:
            return float(np.linalg.norm(self.eval(y)))
        except Exception:
            return float("inf")

    def residuals_batch(self, Y: "np.ndarray") -> "np.ndarray":
        try:
            F = self.eval_batch(Y)
            rr = np.linalg.norm(F, axis=1)
            rr = np.asarray(rr, dtype=float)
            rr[~np.isfinite(rr)] = np.inf
            return rr
        except Exception:
            YY = np.asarray(Y, dtype=np.complex128)
            if YY.ndim == 1:
                return np.asarray([float("inf")], dtype=float)
            return np.full(int(YY.shape[0]), float("inf"), dtype=float)


# ---------------------------------------------------------------------------
# Geometric single flow: Mobius/Riemann + unconditional homothety + startopt
# ---------------------------------------------------------------------------

DEFAULT_POWERS = (
    [2.0 ** k for k in range(-20, 25)] +
    [3.0 * (2.0 ** k) for k in range(-16, 22)] +
    [5.0 * (2.0 ** k) for k in range(-14, 20)] +
    [10.0 ** k for k in range(-6, 7)]
)
DEFAULT_ANGLES_DEG = [0, 6, 12, 18, 24, 32, 40, 48, 56, 64, 72, 80, 86, 89, 90, 91, 94, 100, 108, 116, 128, 140, 152, 164, 172]
DEFAULT_RADII = [0.025, 0.04, 0.06, 0.08, 0.12, 0.18, 0.27, 0.40, 0.60, 0.85, 1.15, 1.55, 2.05, 2.75, 3.60, 4.80, 6.40]
DEFAULT_GAINS = [0.035, 0.055, 0.085, 0.12, 0.18, 0.27, 0.40, 0.58, 0.78, 1.0, 1.28, 1.65, 2.2, 3.0, 4.2, 6.0, 8.5, 12.0]


def raw_direction(n: int, trial: int, seed: int, normalize: bool = True) -> "np.ndarray":
    ensure_numpy()
    vals = []
    for j in range(n):
        h1 = splitmix64(seed + 0xD1A5E + 0x1000003 * trial + 0x9E37 * (j + 1))
        h2 = splitmix64(seed + 0xBADC0DE + 0x1000033 * trial + 0xC2B2 * (j + 1))
        ang = 2.0 * math.pi * u01(h1)
        amp = math.exp(0.45 * (2.0 * u01(h2) - 1.0))
        vals.append(amp * phase(ang))
    v = np.asarray(vals, dtype=np.complex128)
    if normalize:
        nm = float(np.linalg.norm(v))
        if nm > 0:
            v = v / nm * math.sqrt(n)
    return v


def mobius_homothety_start(
    n: int,
    trial: int,
    seed: int,
    powers: Sequence[float],
    angles: Sequence[float],
    radii: Sequence[float],
    cap: float,
    roots_found: int = 0,
    duplicates: int = 0,
    failures: int = 0,
    target_count: int = 1,
) -> tuple[Any, dict[str, Any]]:
    """Single universal Riemann/Mobius + Thales homothety start.

    This is the lightweight heuristic engine.  It is still one formula and one
    flow: the heuristic only changes the parameters (Lambda, theta, radius) of
    the same Mobius/homothety map.  No policy switch and no solver fallback.

    Formula per coordinate:
        y_j = Lambda * pole_j * (cos(theta_j)*(u_j/pole_j)+sin(theta_j))
                           / (-sin(theta_j)*(u_j/pole_j)+cos(theta_j))

    Heuristic idea:
      - Lambda is always drawn from a very wide Thales-power ladder.
      - When duplicates accumulate, Lambda is pushed upward to escape known
        basins; this is a formula-level amplification, not a separate branch.
      - theta is rotated through affine, Riemann, and near-infinity charts.
    """
    ensure_numpy()
    powers2 = sorted(set(min(max(float(x), 1e-300), float(cap)) for x in powers if float(x) > 0))
    if not powers2:
        powers2 = [1.0]
    Lp, La, Lr = len(powers2), max(1, len(angles)), max(1, len(radii))
    # Low-discrepancy permutation across the power ladder, with unconditional
    # high-power thrust.  This is intentionally not a fallback or multi-policy:
    # it is one parameter formula for Lambda.
    phi = 0.6180339887498948482
    q = (trial * phi + 0.071 * roots_found + 0.013 * duplicates) % 1.0
    power_index = (int(q * Lp) + 37 * trial + 11 * (trial // 7) + 5 * roots_found) % Lp
    base_power = powers2[power_index]
    dup_pressure = (duplicates + 1.0) / (roots_found + 1.0)
    fail_pressure = (failures + 1.0) / (trial + 1.0)
    progress = min(1.0, max(0.0, roots_found / max(1.0, float(target_count))))
    thrust_ladder = [1.0, 1.6, 2.5, 4.0, 6.5, 10.0, 16.0, 25.0, 40.0, 64.0, 100.0, 160.0, 256.0]
    thrust = thrust_ladder[(trial * 17 + roots_found * 3 + duplicates) % len(thrust_ladder)]
    # More roots found -> explore outer charts more; more duplicates -> stronger
    # homothety; more failures -> soften slightly but keep the same flow.
    amp = (thrust ** (0.18 + 0.82 * progress)) * ((1.0 + dup_pressure) ** 0.42) / ((1.0 + 0.25 * fail_pressure) ** 0.15)
    pwr = min(float(cap), max(1e-300, base_power * amp))

    theta0 = angles[(trial * 19 + roots_found * 7 + duplicates * 3) % La]
    theta_jitter = math.radians(4.0) * math.sin(1.324717957244746 * (trial + 1) + 0.31 * roots_found)
    radius0 = radii[(trial * 13 + failures * 5 + roots_found * 2) % Lr]
    radius = max(1e-300, float(radius0) * math.exp(0.22 * math.sin(0.754877666 * (trial + 1) + 0.17 * duplicates)))

    d = raw_direction(n, trial, seed, True)
    u = radius * d
    out = np.empty(n, dtype=np.complex128)
    poles = []
    theta_values = []
    for j in range(n):
        # Coordinate phase and tiny theta offsets are generic in n; no bivariate
        # special case is used anywhere.
        hj = splitmix64(seed + 0xA11CE + 982451653 * trial + 1009 * (j + 1))
        pole = phase(2.0 * math.pi * u01(hj))
        poles.append(pole)
        tj = theta0 + theta_jitter * math.cos(0.5 + j) + math.radians(2.0) * math.sin((j + 1) * (trial + 1) * 0.38196601125)
        theta_values.append(tj)
        c, s = math.cos(tj), math.sin(tj)
        w = u[j] / pole
        denom = (-s * w + c)
        if abs(denom) < 1e-12:
            # This is a deterministic chart regularization, not a fallback to a
            # different algorithm.  It keeps the same projective formula finite.
            denom += 1e-12 * phase(0.37 + j)
        out[j] = pwr * pole * ((c * w + s) / denom)
    meta = {
        "homothety": float(pwr),
        "base_homothety": float(base_power),
        "thales_thrust": float(thrust),
        "theta_deg": float(math.degrees(theta0)),
        "theta_jitter_deg": float(math.degrees(theta_jitter)),
        "theta_mean_deg": float(sum(math.degrees(t) for t in theta_values) / max(1, len(theta_values))),
        "base_radius": float(radius),
        "dup_pressure": float(dup_pressure),
        "fail_pressure": float(fail_pressure),
        "progress": float(progress),
        "chart": "single-flow/pure-mobius-riemann-thales",
    }
    return out, meta

def finite_residual(target: TargetTrack, y: Sequence[complex]) -> float:
    r = target.residual(y)
    return r if math.isfinite(r) else float("inf")


def _fold_highdim_to_complex(vals: list[complex], n: int, radius: float) -> "np.ndarray":
    ensure_numpy()
    out = np.zeros(n, dtype=np.complex128)
    for k, v in enumerate(vals):
        j = k % max(1, n)
        # deterministic complex phase while folding, to avoid aliasing high-D
        # opposite vertices into the same low-D direction.
        out[j] += complex(v) * phase(2.0 * math.pi * (((k + 1) * 0.6180339887498948 + (j + 1) * 0.41421356237) % 1.0))
    nm = max(1e-300, float(np.linalg.norm(out)))
    return out / nm * (float(radius) * math.sqrt(max(1, n)))


def _highdim_universal_complex_cell(n: int, idx: int, layer: int, radius: float) -> "np.ndarray":
    """High-dimensional universal geometry folded into C^n.

    The geometry is fixed and independent of F.  We generate vertices from
    hypercubes and regular simplexes in dimensions larger than n, then fold them
    into C^n by a deterministic phase projection.
    """
    ensure_numpy()
    dims = [max(4, n + 2), max(6, 2 * n + 2), max(8, 4 * n + 2), 12]
    D = dims[(idx + layer) % len(dims)]
    mode = (idx // len(dims)) % 4
    vals: list[complex] = []
    if mode == 0:
        # High-D hypercube vertex with complex quadrant signs.
        for k in range(D):
            h = splitmix64(0x305C0BE + 1000003 * idx + 9176 * layer + 193 * (k + 1))
            q = int(4 * u01(h)) % 4
            vals.append([1.0, 1j, -1.0, -1j][q])
    elif mode == 1:
        # High-D simplex-like: one dominant vertex against balanced tail.
        j0 = (idx + 3 * layer) % D
        for k in range(D):
            amp = math.sqrt(D) if k == j0 else -1.0 / math.sqrt(max(1, D))
            vals.append(amp * phase(2.0 * math.pi * (((k + 1) * 0.6180339887498948 + idx * 0.071) % 1.0)))
    elif mode == 2:
        # Cross-polytope pair: sparse high-D axis, folded densely into C^n.
        j0 = (idx * 5 + layer) % D
        for k in range(D):
            vals.append((1.0 if k == j0 else 0.18) * phase(2.0 * math.pi * (((k + 1) * (idx + 1) * 0.3027756377) % 1.0)))
    else:
        # Deterministic high-D projective spiral.
        for k in range(D):
            amp = math.exp(0.32 * math.sin((idx + 1) * (k + 1) * 1.324717957244746 + layer))
            vals.append(amp * phase(2.0 * math.pi * (((idx + 1) * 0.754877666 + (k + 1) * 0.569840291) % 1.0)))
    return _fold_highdim_to_complex(vals, n, radius)


def _universal_complex_cell(n: int, idx: int, layer: int, radius: float) -> "np.ndarray":
    """Fixed universal atlas cell in C^n, independent of F.

    Geometry mix:
      - simplex/coordinate vertices for axis visibility;
      - hypercube sign/phase corners for broad coverage;
      - golden projective spiral for generic directions;
      - shell radius supplied by a fixed radius ladder.
    """
    ensure_numpy()
    phi = 0.6180339887498948482
    vals = np.empty(n, dtype=np.complex128)
    mode = idx % 4
    if mode == 0:
        # simplex-like: one dominant coordinate plus a small balanced tail
        j0 = (idx // 4 + layer) % max(1, n)
        for j in range(n):
            amp = 1.0 if j == j0 else 1.0 / math.sqrt(max(1, n))
            ang = 2.0 * math.pi * ((idx + 1) * (j + 1) * phi + 0.137 * layer)
            vals[j] = amp * phase(ang)
    elif mode == 1:
        # hypercube phases: fixed low-discrepancy corners in complex signs
        for j in range(n):
            h = splitmix64(0x304C0BE + 65537 * idx + 4099 * layer + 193 * (j + 1))
            q = int(4 * u01(h)) % 4
            vals[j] = [1.0, 1j, -1.0, -1j][q]
    elif mode == 2:
        # projective golden spiral
        for j in range(n):
            amp = math.exp(0.35 * math.sin((idx + 1) * (j + 1) * 1.324717957244746 + layer))
            ang = 2.0 * math.pi * (((idx + 1) * phi + (j + 1) * phi * phi + 0.071 * layer) % 1.0)
            vals[j] = amp * phase(ang)
    else:
        # reciprocal/infinity-biased mixed shell: alternating large/small coordinates
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
    seed: int,
    powers: Sequence[float],
    angles: Sequence[float],
    radii: Sequence[float],
    cap: float,
    roots_found: int,
    duplicates: int,
    failures: int,
    target_count: int,
    universal_cells: int = 16,
    universal_shells: int = 5,
    cell_probe_radius: float = 0.14,
    cell_descent_min: float = 1.02,
    cell_equal_gap_min: float = 1e-10,
    cell_log_max: float = 80.0,
    universal_cycle: bool = True,
) -> tuple[Any, dict[str, Any]]:
    """305 fixed universal atlas.

    The atlas geometry is fixed for dimension n: shells times deterministic
    simplex/hypercube/projective cells.  It is not adapted to F and does not use
    a scalar score.  F only activates/deactivates cells by predicate tests.
    """
    ensure_numpy()
    # Fixed shell ladder: combine user radii/powers into bounded projective
    # scales.  The sequence is deterministic and independent of F.
    base_shells = [0.05, 0.12, 0.27, 0.60, 1.0, 1.7, 3.0, 5.5, 10.0, 18.0]
    shell_count = max(1, int(universal_shells))
    cell_count = max(1, int(universal_cells))
    admissible = []
    tested = 0
    for kk in range(cell_count):
        # Cycle deterministically through the universal atlas so each trial sees
        # a different window, but the atlas itself is fixed.
        atlas_idx = trial * cell_count + kk
        layer = atlas_idx // max(1, cell_count)
        radius = base_shells[layer % len(base_shells)]
        # Add a fixed homothety shell from DEFAULT_POWERS without using F.
        if powers:
            radius *= min(float(cap), max(1e-300, float(powers[(atlas_idx + 7 * layer) % len(powers)]))) ** 0.15

        cell_candidates: list[tuple[str, Any, dict[str, Any]]] = []
        cell_candidates.append((
            "highdim-hypercube-simplex-projective",
            _highdim_universal_complex_cell(n, atlas_idx, layer, radius),
            {"radius": float(radius), "theta_deg": None, "base_radius": float(radius)},
        ))
        cell_candidates.append((
            "native-lowdim-simplex-hypercube-projective",
            _universal_complex_cell(n, atlas_idx, layer, radius),
            {"radius": float(radius), "theta_deg": None, "base_radius": float(radius)},
        ))
        # Fixed Mobius/Riemann cell: same universal formula as 301, but with
        # pressures frozen to zero so the geometry is not adapted to F or to the
        # root history.  This adds projective/infinity coverage that plain
        # hypercube cells miss.
        try:
            my, mm = mobius_homothety_start(
                n, atlas_idx, seed + 0x304F1C + 65537 * kk,
                powers, angles, radii, cap,
                roots_found=0, duplicates=0, failures=0, target_count=max(1, target_count),
            )
            cell_candidates.append(("fixed-mobius-riemann", np.asarray(my, dtype=np.complex128), dict(mm)))
        except Exception:
            pass

        for geom_name, y, cmeta in cell_candidates:
            y = np.asarray(y, dtype=np.complex128)
            tested += 1
            try:
                f0 = target.eval(y)
                r0 = float(np.linalg.norm(f0))
            except Exception:
                continue
            if not math.isfinite(r0):
                continue
            log0 = _log_stability_energy_201(y)
            if not math.isfinite(log0) or log0 > float(cell_log_max):
                continue
            ynorm = max(1.0, float(np.linalg.norm(y)))
            # Fixed local cell probes: radial in/out, phase rotation, coordinate simplex.
            probes = [0.92 * y, 1.08 * y, y * phase(float(cell_probe_radius))]
            ej = np.zeros(n, dtype=np.complex128)
            ej[(atlas_idx + layer) % max(1, n)] = 1.0
            probes.append(y + float(cell_probe_radius) * ynorm * ej)
            probes.append(y - float(cell_probe_radius) * ynorm * ej)
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
            if descent >= float(cell_descent_min) and min_gap >= float(cell_equal_gap_min):
                # Predicate representative: first decreasing probe, not best probe.
                rep = y.copy()
                for ii, rr in enumerate(RP):
                    if math.isfinite(float(rr)) and float(rr) < r0 / float(cell_descent_min):
                        rep = P[int(ii)].copy()
                        break
                admissible.append({
                    "atlas_idx": int(atlas_idx), "layer": int(layer), "radius": float(cmeta.get("homothety", cmeta.get("radius", radius))),
                    "geometry_cell": str(geom_name), "cell_meta": cmeta,
                    "y": rep, "residual": float(r0), "log_energy": float(log0),
                    "descent_ratio": float(descent), "equal_gap_min": float(min_gap),
                    "representative_changed": bool(not np.allclose(rep, y)),
                })

    if admissible:
        pick = (trial + roots_found + 2 * duplicates + 3 * failures) % len(admissible) if bool(universal_cycle) else 0
        rec = admissible[int(pick)]
        meta = {
            "homothety": float(rec["radius"]),
            "base_homothety": float(rec["radius"]),
            "thales_thrust": 1.0,
            "theta_deg": None,
            "theta_jitter_deg": 0.0,
            "theta_mean_deg": None,
            "base_radius": float(rec["radius"]),
            "dup_pressure": float((duplicates + 1.0) / (roots_found + 1.0)),
            "fail_pressure": float((failures + 1.0) / (trial + 1.0)),
            "progress": float(min(1.0, max(0.0, roots_found / max(1.0, float(target_count))))),
            "chart": "305-highdim-universal-atlas/" + str(rec.get("geometry_cell", "hybrid")),
            "atlas_mode": "305-highdim-universal-atlas-no-score",
            "atlas_geometry": "highdim-hypercube+simplex+crosspolytope+projective-folded-to-Cn + fixed-mobius-shells",
            "atlas_selected_geometry_cell": str(rec.get("geometry_cell", "unknown")),
            "atlas_cells_tested": int(tested),
            "atlas_admissible_cells": int(len(admissible)),
            "atlas_selected_index": int(rec["atlas_idx"]),
            "atlas_selected_layer": int(rec["layer"]),
            "atlas_selection_rule": "cyclic-admissible-cell" if bool(universal_cycle) else "first-admissible-cell",
            "atlas_cell_residual": float(rec["residual"]),
            "atlas_cell_log_energy": float(rec["log_energy"]),
            "atlas_cell_descent_ratio": float(rec["descent_ratio"]),
            "atlas_cell_equal_gap_min": float(rec["equal_gap_min"]),
            "atlas_representative_changed": bool(rec["representative_changed"]),
        }
        return np.asarray(rec["y"], dtype=np.complex128).copy(), meta

    # Deterministic fixed-cell fallback if no cell passed predicates.  Not a
    # score winner; just the next universal cell.
    y = _universal_complex_cell(n, trial, trial // max(1, cell_count), 1.0)
    meta = {
        "homothety": 1.0, "base_homothety": 1.0, "thales_thrust": 1.0,
        "theta_deg": None, "theta_jitter_deg": 0.0, "theta_mean_deg": None,
        "base_radius": 1.0, "dup_pressure": float((duplicates + 1.0) / (roots_found + 1.0)),
        "fail_pressure": float((failures + 1.0) / (trial + 1.0)),
        "progress": float(min(1.0, max(0.0, roots_found / max(1.0, float(target_count))))),
        "chart": "305-highdim-universal-atlas/default-cell",
        "atlas_mode": "305-highdim-universal-atlas-no-score",
        "atlas_geometry": "highdim-hypercube+simplex+crosspolytope+projective-folded-to-Cn + fixed-mobius-shells",
            "atlas_selected_geometry_cell": str(rec.get("geometry_cell", "unknown")),
        "atlas_cells_tested": int(tested), "atlas_admissible_cells": 0,
        "atlas_selected_index": None, "atlas_selection_rule": "default-universal-cell-after-no-admissible-cell",
    }
    return y, meta


def startopt(target: TargetTrack, y0: Any, trial: int, seed: int, steps: int, candidates: int, gains: Sequence[float], micro_epochs: int) -> tuple[Any, dict[str, Any]]:
    """Batched start optimization in the same Thales chart.

    118 scored start candidates one-by-one.  120 keeps the same geometric pool
    but evaluates it as a NumPy batch, which is both faster and less biased by
    early candidate order.
    """
    ensure_numpy()
    best = np.asarray(y0, dtype=np.complex128).copy()
    best_r = finite_residual(target, best)
    initial = best_r
    evals = 1
    micro_total = 0
    chosen_gain = 1.0
    for step in range(max(0, steps)):
        base = best.copy()
        pool: list[tuple[float, Any]] = [(1.0, base)]
        for c in range(max(0, candidates - 1)):
            gain = float(gains[(trial + 3 * step + c) % len(gains)])
            cand = gain * base
            if c % 3 == 1:
                pert = []
                for j, val in enumerate(cand):
                    h1 = splitmix64(seed + 0x5157A47 + 65537 * trial + 4099 * c + 193 * (j + 1))
                    h2 = splitmix64(seed + 0x7150F7 + 104729 * trial + 8191 * c + 389 * (j + 1))
                    ph = 0.23 * (2.0 * u01(h1) - 1.0)
                    amp = math.exp(0.28 * (2.0 * u01(h2) - 1.0))
                    pert.append(val * amp * phase(ph))
                cand = np.asarray(pert, dtype=np.complex128)
            elif c % 3 == 2:
                fresh = raw_direction(len(base), trial + 31 * (step + 1) + c, seed, True)
                nm = max(1e-300, float(np.linalg.norm(base)))
                cand = 0.70 * cand + 0.30 * gain * nm * fresh / max(1e-300, float(np.linalg.norm(fresh)))
            pool.append((gain, np.asarray(cand, dtype=np.complex128)))

        Y = np.asarray([cand for _, cand in pool], dtype=np.complex128)
        R = target.residuals_batch(Y)
        evals += int(len(pool))
        scores = np.asarray([math.log1p(max(0.0, float(r))) + 1e-15 * math.log1p(float(np.linalg.norm(Y[idx]))) if math.isfinite(float(r)) else float("inf") for idx, r in enumerate(R)], dtype=float)
        idx_best = int(np.nanargmin(scores)) if np.any(np.isfinite(scores)) else 0
        cand_best = Y[idx_best].copy()
        r_best = float(R[idx_best]) if idx_best < len(R) else float("inf")
        if micro_epochs > 0:
            # Optional micro-corrector remains pure Pandrosion and is only run
            # on the best batched start candidate to keep 120 lightweight.
            loc = pandrosion_corrector(target, cand_best, max_epochs=micro_epochs, tol=1e-12, accept=0.0, trial_timeout=0.0, line_search=6)
            micro_total += int(loc.get("epochs", 0))
            if float(loc.get("residual", float("inf"))) < r_best:
                cand_best = np.asarray(loc["y"], dtype=np.complex128)
                r_best = float(loc["residual"])
        old_score = math.log1p(max(0.0, best_r)) + 1e-15 * math.log1p(float(np.linalg.norm(best))) if math.isfinite(best_r) else float("inf")
        new_score = math.log1p(max(0.0, r_best)) + 1e-15 * math.log1p(float(np.linalg.norm(cand_best))) if math.isfinite(r_best) else float("inf")
        if new_score < old_score:
            best = cand_best
            best_r = float(r_best)
            chosen_gain = float(pool[idx_best][0])
    return best, {
        "startopt_enabled": bool(steps > 0),
        "startopt_r0": float(initial),
        "startopt_r1": float(best_r),
        "startopt_ratio": (float(best_r / initial) if math.isfinite(best_r) and math.isfinite(initial) and initial > 0 else None),
        "startopt_steps": int(max(0, steps)),
        "startopt_evals": int(evals),
        "startopt_micro_epochs": int(micro_total),
        "startopt_gain": float(chosen_gain),
        "startopt_batch_numpy": True,
    }


# ---------------------------------------------------------------------------
# Corrector and root handling
# ---------------------------------------------------------------------------


def _log_stability_energy_201(y: Any, eps: float = 1e-12) -> float:
    """Projective/logarithmic scale energy inspired by Paper 0.

    This is not a solver and not a derivative. It is a scalar diagnostic:
        L(y) = || log(1 + |y_j|) ||_2.
    It penalizes probes that create unnecessary scale distortion.
    """
    ensure_numpy()
    yy = np.asarray(y, dtype=np.complex128).ravel()
    if yy.size == 0:
        return 0.0
    val = np.log1p(np.abs(yy) + float(eps))
    val[~np.isfinite(val)] = 0.0
    return float(np.linalg.norm(val) / math.sqrt(max(1, yy.size)))


def _log_stability_batch_201(Y: Any, eps: float = 1e-12) -> "np.ndarray":
    ensure_numpy()
    YY = np.asarray(Y, dtype=np.complex128)
    if YY.ndim == 1:
        YY = YY[None, :]
    vals = np.log1p(np.abs(YY) + float(eps))
    vals[~np.isfinite(vals)] = 0.0
    return np.linalg.norm(vals, axis=1) / math.sqrt(max(1, YY.shape[1]))


def _probe_endpoint_candidates(
    target: TargetTrack,
    y: "np.ndarray",
    f: "np.ndarray",
    residual: float,
    prev_delta: Optional["np.ndarray"],
    ep: int,
    direction_seed: int,
    probe_scale: float,
    probe_candidates: int,
    probe_radii: Sequence[float],
    include_self_probe: bool,
    condition_top: int = 2,
    condition_weight: float = 0.0025,
    axis_probes: bool = True,
    equal_value_weight: float = 0.015,
    equal_value_eps: float = 1e-12,
    curvature_top: int = 2,
    curvature_weight: float = 0.003,
    curvature_mid: float = 0.5,
    probe_log_weight: float = 0.0005,
    probe_log_delta_weight: float = 0.0015,
) -> tuple[Any, dict[str, Any]]:
    """Choose the finite-slope endpoint b by batched residual scoring.

    120 keeps the 118/120 theorem-guided principle b_* = argmin ||G(b)||, but uses
    a NumPy batch to score all probes at once.  Among the top residual probes it
    optionally adds penalties for the condition number of the exact finite-slope
    matrix Q_G(a,b), for near equal-value poles F(b)≈F(a), and for large
    finite-slope curvature measured by slope variation between a midpoint and b.
    This is not a Newton fallback: Q_G is the same telescopic Pandrosion slope
    used for the accepted step.
    """
    ensure_numpy()
    n = len(y)
    ynorm = max(1.0, float(np.linalg.norm(y)))
    radii = [float(r) for r in probe_radii if float(r) >= 0]
    if not radii:
        radii = [1.0]
    budget = max(1, int(probe_candidates))
    candidates: list[tuple[str, Any]] = []
    if include_self_probe:
        candidates.append(("self", y.copy()))

    if prev_delta is not None and np.all(np.isfinite(prev_delta)):
        pdn = max(1e-300, float(np.linalg.norm(prev_delta)))
        base = prev_delta / pdn * min(max(pdn, float(probe_scale) * ynorm), 2.5 * ynorm)
        candidates.append(("inertial", y + base))

    if axis_probes and budget > len(candidates):
        # Coordinate Thales probes help multivariate finite slopes see each
        # coordinate separately before generic spiral probes are used.
        rad0 = float(probe_scale) * ynorm * (radii[min(1, len(radii)-1)] if radii else 1.0)
        for j in range(n):
            if len(candidates) >= budget:
                break
            ej = np.zeros(n, dtype=np.complex128)
            ej[j] = 1.0 + 0.0j
            ph = phase(0.173 + 0.6180339887498948 * (ep + 1) + 1.41421356237 * (j + 1))
            candidates.append((f"axis-{j}", y + rad0 * ph * ej))

    k = 0
    while len(candidates) < budget:
        rad = float(probe_scale) * ynorm * radii[k % len(radii)]
        qdir = raw_direction(n, direction_seed + 104729 * (ep + 1) + 7919 * (k + 1), direction_seed ^ (0x120120 + 17 * k), True)
        qnorm = max(1e-300, float(np.linalg.norm(qdir)))
        qdir = qdir / qnorm * math.sqrt(max(1, n))
        ph = phase(0.6180339887498948 * (ep + 1) + 2.399963229728653 * (k + 1))
        step = rad * ph * qdir
        if float(np.linalg.norm(y)) > 0:
            step = step + (0.12 * rad) * y / ynorm * phase(0.38196601125 * (k + 1))
        tiny = 1e-12 * ynorm
        for j in range(n):
            if abs(step[j]) < tiny:
                step[j] += tiny * phase(0.17 + j + ep + k)
        candidates.append((f"geom-{k}", y + step))
        k += 1

    cand_slice = candidates[:budget]
    names = [name for name, _ in cand_slice]
    B = np.asarray([np.asarray(b, dtype=np.complex128) for _, b in cand_slice], dtype=np.complex128)
    FB = target.eval_batch(B)
    R = np.linalg.norm(FB, axis=1)
    D = np.linalg.norm(B - y[None, :], axis=1)
    scores = np.asarray([math.log1p(max(0.0, float(R[i]))) + 1e-14 * math.log1p(float(D[i])) if math.isfinite(float(R[i])) else float("inf") for i in range(len(cand_slice))], dtype=float)

    # 301 logarithmic/projective stability: probes that violently change
    # the logarithmic scale are penalized. This is the algorithmic version of
    # the Paper-0 rule "make |log y| small before correcting".
    log_y = _log_stability_energy_201(y)
    log_B = _log_stability_batch_201(B)
    log_delta = np.abs(log_B - log_y)
    if float(probe_log_weight) > 0:
        scores = scores + float(probe_log_weight) * np.log1p(log_B)
    if float(probe_log_delta_weight) > 0:
        scores = scores + float(probe_log_delta_weight) * np.log1p(log_delta)

    # Equal-value pole avoidance.  In the anchored map the poles satisfy
    # F(b)=F(a).  The self point b=a is a removable limit, so it is excluded
    # from the penalty; all other near-equal-value probes are discouraged.
    Fy_norm = float(np.linalg.norm(f))
    Fdiff = np.linalg.norm(FB - np.asarray(f, dtype=np.complex128)[None, :], axis=1)
    Fscale = R + Fy_norm + 1e-300
    equal_gap = Fdiff / Fscale
    nonself = D > (1e-12 * ynorm)
    pole_penalty = np.zeros_like(scores)
    if float(equal_value_weight) > 0:
        safe_gap = np.maximum(equal_gap, float(equal_value_eps))
        pole_penalty = np.where(nonself, float(equal_value_weight) * np.log1p(1.0 / safe_gap), 0.0)
        scores = scores + pole_penalty

    cond_best = None
    cond_scored = 0
    curvature_best = None
    curvature_scored = 0
    if int(condition_top) > 0 and len(cand_slice) > 0 and np.any(np.isfinite(scores)):
        order = np.argsort(scores)
        for idx in order[: max(0, min(int(condition_top), len(order)))]:
            try:
                Q = target.slope_matrix(y, B[int(idx)])
                cnd = float(np.linalg.cond(Q))
            except Exception:
                cnd = float("inf")
            cond_scored += 1
            if math.isfinite(cnd):
                scores[int(idx)] += float(condition_weight) * math.log1p(max(0.0, cnd))
            else:
                scores[int(idx)] = float("inf")

    # Tensor/curvature proxy using only finite slopes: compare Q(a,b) with
    # Q(a, a + mid*(b-a)).  Large variation means the averaged Jacobian is
    # changing rapidly along the probe segment, which is the multidimensional
    # analogue of a large curvature term.
    if int(curvature_top) > 0 and float(curvature_weight) > 0 and len(cand_slice) > 0 and np.any(np.isfinite(scores)):
        order = np.argsort(scores)
        curv_vals = []
        for idx in order[: max(0, min(int(curvature_top), len(order)))]:
            ii = int(idx)
            if not bool(nonself[ii]):
                continue
            try:
                bvec = B[ii]
                mid = y + float(curvature_mid) * (bvec - y)
                Q_full = target.slope_matrix(y, bvec)
                Q_mid = target.slope_matrix(y, mid)
                curv = float(np.linalg.norm(Q_full - Q_mid) / max(1e-300, np.linalg.norm(Q_full)))
            except Exception:
                curv = float("inf")
            curvature_scored += 1
            if math.isfinite(curv):
                scores[ii] += float(curvature_weight) * math.log1p(max(0.0, curv))
                curv_vals.append(curv)
            else:
                scores[ii] = float("inf")
        if curv_vals:
            curvature_best = float(min(curv_vals))

    if not np.any(np.isfinite(scores)):
        raise RuntimeError("no-finite-probe")
    best_idx = int(np.nanargmin(scores))
    best_name = names[best_idx]
    best_b = B[best_idx].copy()
    best_res = float(R[best_idx])
    best_distance = float(D[best_idx])
    if int(condition_top) > 0:
        try:
            cond_best = float(np.linalg.cond(target.slope_matrix(y, best_b)))
        except Exception:
            cond_best = None

    return best_b, {
        "probe_mode": "120-equal-value-pole-and-tensor-aware-batched-probe",
        "probe_name": best_name,
        "probe_candidates": int(len(cand_slice)),
        "probe_evals": int(len(cand_slice)),
        "probe_residual": float(best_res),
        "probe_distance": float(best_distance),
        "probe_improvement_proxy": (float(residual / best_res) if math.isfinite(best_res) and best_res > 0 and math.isfinite(residual) else None),
        "probe_self_enabled": bool(include_self_probe),
        "probe_axis_enabled": bool(axis_probes),
        "probe_condition_top": int(condition_top),
        "probe_condition_scored": int(cond_scored),
        "probe_condition_best": cond_best,
        "probe_equal_value_weight": float(equal_value_weight),
        "probe_equal_value_gap": float(equal_gap[best_idx]) if len(equal_gap) else None,
        "probe_equal_value_gap_min": float(np.nanmin(equal_gap)) if len(equal_gap) else None,
        "probe_equal_value_penalty": float(pole_penalty[best_idx]) if len(pole_penalty) else 0.0,
        "probe_curvature_top": int(curvature_top),
        "probe_curvature_scored": int(curvature_scored),
        "probe_curvature_best": curvature_best,
        "probe_log_energy_current": float(log_y),
        "probe_log_energy_best": float(log_B[best_idx]) if len(log_B) else None,
        "probe_log_delta_best": float(log_delta[best_idx]) if len(log_delta) else None,
        "probe_log_weight": float(probe_log_weight),
        "probe_log_delta_weight": float(probe_log_delta_weight),
        "probe_batch_numpy": True,
    }


def _line_lambdas(line_search: int, line_grid: Optional[Sequence[float]] = None) -> list[float]:
    vals: list[float] = []
    if line_grid is not None:
        for x in line_grid:
            try:
                xx = float(x)
                if math.isfinite(xx) and xx > 0:
                    vals.append(xx)
            except Exception:
                pass
    if not vals:
        vals = [1.0, 0.75, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03125, 0.02]
    k = 0
    while len(vals) < max(1, int(line_search)):
        vals.append(1.0 / (2.0 ** k))
        k += 1
    # keep order but remove near-duplicates
    out: list[float] = []
    for v in vals:
        vv = min(max(float(v), 1e-15), 4.0)
        if all(abs(vv - u) > 1e-15 for u in out):
            out.append(vv)
        if len(out) >= max(1, int(line_search)):
            break
    return out or [1.0]



def _halley_gate_201(residual: float, cond: Optional[float], gate_residual: float, cond_weight: float, min_gate: float, log_energy: float = 0.0, log_weight: float = 0.0, log_scale: float = 1.0) -> float:
    """Continuous gate for the cubic finite-slope Halley layer.

    Far from a root, this returns nearly zero and the engine behaves like 120/118.
    Near a root with an acceptable finite-slope condition, the gate approaches one.
    """
    try:
        r = max(0.0, float(residual))
    except Exception:
        r = 1e300
    gr = max(1e-300, float(gate_residual))
    rg = 1.0 / (1.0 + (r / gr) ** 2)
    if cond is None or not math.isfinite(float(cond)):
        cg = 0.0
    else:
        cg = 1.0 / (1.0 + float(cond_weight) * math.log1p(max(0.0, float(cond))))
    try:
        le = max(0.0, float(log_energy))
    except Exception:
        le = 0.0
    lg = 1.0 / (1.0 + float(log_weight) * le / max(1e-12, float(log_scale)))
    g = max(0.0, min(1.0, rg * cg * lg))
    return float(g if g >= float(min_gate) else 0.0)


def _gated_tensor_halley_delta_201(
    target: TargetTrack,
    y: "np.ndarray",
    f: "np.ndarray",
    Q: "np.ndarray",
    delta_raw: "np.ndarray",
    residual: float,
    cond: Optional[float],
    gate_residual: float = 0.25,
    probe_fraction: float = 0.50,
    cond_weight: float = 0.025,
    min_gate: float = 0.04,
    max_correction: float = 1.25,
    log_energy: float = 0.0,
    halley_log_weight: float = 0.03,
    halley_log_scale: float = 1.0,
    tensor_extra_directions: int = 0,
) -> tuple["np.ndarray", dict[str, Any]]:
    """Directional multivariate finite-slope Halley correction.

    Scalar 0-bis:
        H3 = a - 2 f D1 / (2 D1^2 - f D2)

    Vector form used here:
        Q1 delta1 = -F(a)
        Q1 delta2 = -1/2 H(delta1,delta1)
        delta_H = delta1 + delta2

    Q1 is the exact Pandrosion finite-slope matrix already used by 118/120.
    H(delta1,delta1) is estimated by a symmetric finite-slope probe along delta1.
    No analytic Jacobian or Hessian is formed.
    """
    ensure_numpy()
    gate = _halley_gate_201(residual, cond, gate_residual, cond_weight, min_gate, log_energy=log_energy, log_weight=halley_log_weight, log_scale=halley_log_scale)
    meta: dict[str, Any] = {
        "halley_gate": float(gate),
        "halley_used": False,
        "halley_probe_fraction": float(probe_fraction),
        "halley_delta2_norm": None,
        "halley_raw_norm": float(np.linalg.norm(delta_raw)),
        "halley_log_energy": float(log_energy),
        "halley_tensor_extra_directions": int(tensor_extra_directions),
    }
    if gate <= 0.0:
        return delta_raw, meta

    dnorm = float(np.linalg.norm(delta_raw))
    ynorm = max(1.0, float(np.linalg.norm(y)))
    if (not math.isfinite(dnorm)) or dnorm <= 1e-300:
        return delta_raw, meta

    frac = max(1e-6, min(1.0, float(probe_fraction)))
    hvec = frac * delta_raw

    try:
        # Direction 1: raw Pandrosion step.
        fp = target.eval(y + hvec)
        fm = target.eval(y - hvec)
        sec2_terms = [(fp - 2.0 * f + fm) / (frac * frac)]

        # 301 tensorial enrichment: add a small number of deterministic
        # transverse directions. They are cheap and only active when the gate is
        # already local, so they do not dominate the global search.
        extra = max(0, int(tensor_extra_directions))
        if extra > 0:
            dn = max(1e-300, float(np.linalg.norm(delta_raw)))
            for kk in range(extra):
                qdir = raw_direction(len(y), 0x201000 + 7919 * (kk + 1), 0x201777 + kk, True)
                qdir = np.asarray(qdir, dtype=np.complex128)
                qdir = qdir / max(1e-300, float(np.linalg.norm(qdir))) * dn
                if kk == 0:
                    qdir = 1j * delta_raw
                hvec2 = frac * qdir
                try:
                    fp2 = target.eval(y + hvec2)
                    fm2 = target.eval(y - hvec2)
                    sec2_terms.append((fp2 - 2.0 * f + fm2) / (frac * frac))
                except Exception:
                    pass
        sec2 = np.mean(np.asarray(sec2_terms, dtype=np.complex128), axis=0)
        delta2 = np.linalg.solve(Q, -0.5 * sec2)
        if not np.all(np.isfinite(delta2)):
            return delta_raw, meta

        d2norm = float(np.linalg.norm(delta2))
        max_d2 = float(max_correction) * max(dnorm, 1e-300)
        if d2norm > max_d2:
            delta2 = delta2 * (max_d2 / max(d2norm, 1e-300))
            d2norm = float(np.linalg.norm(delta2))

        if dnorm + d2norm > 20.0 * ynorm:
            scale = (20.0 * ynorm) / max(dnorm + d2norm, 1e-300)
            delta2 = delta2 * scale
            d2norm = float(np.linalg.norm(delta2))

        delta_halley = delta_raw + delta2
        delta_mix = (1.0 - gate) * delta_raw + gate * delta_halley
        meta.update({
            "halley_used": True,
            "halley_delta2_norm": float(d2norm),
            "halley_sec2_norm": float(np.linalg.norm(sec2)),
            "halley_tensor_terms": int(len(sec2_terms)),
            "halley_mix_norm": float(np.linalg.norm(delta_mix)),
        })
        return delta_mix, meta
    except Exception as exc:
        meta.update({"halley_error": type(exc).__name__})
        return delta_raw, meta


def _full_cubic_halley_delta_301(
    target: TargetTrack,
    y: "np.ndarray",
    f: "np.ndarray",
    Q: "np.ndarray",
    delta1: "np.ndarray",
    probe_fraction: float = 0.50,
    max_correction: float = 1.25,
    tensor_extra_directions: int = 0,
) -> tuple["np.ndarray", dict[str, Any]]:
    """Full cubic finite-slope Halley direction, no raw-step fallback.

    Articles/001 (Paper 0 bis) supplies the scalar finite-slope Halley formula.
    Multivariately, Q_G(a,b) is the exact finite-slope analogue of D1/J, and
    the symmetric second finite probe supplies the D2/Hessian defect.
    The cubic layer adds a symmetric finite-probe second defect and solves

        Q delta1 = -F(a),
        Q delta2 = -1/2 H(delta1, delta1),
        delta = delta1 + delta2.

    The returned delta is the cubic candidate itself.  If the tensor probe fails,
    the caller treats the epoch as a cubic failure rather than accepting delta1.
    """
    ensure_numpy()
    meta: dict[str, Any] = {
        "cubic_used": False,
        "cubic_probe_fraction": float(probe_fraction),
        "cubic_delta1_norm": float(np.linalg.norm(delta1)),
        "cubic_delta2_norm": None,
        "tensor_extra_directions": int(tensor_extra_directions),
    }
    dnorm = float(np.linalg.norm(delta1))
    ynorm = max(1.0, float(np.linalg.norm(y)))
    if (not math.isfinite(dnorm)) or dnorm <= 1e-300:
        raise RuntimeError("degenerate-cubic-delta1")
    frac = max(1e-6, min(1.0, float(probe_fraction)))
    sec2_terms = []
    hvec = frac * delta1
    fp = target.eval(y + hvec)
    fm = target.eval(y - hvec)
    sec2_terms.append((fp - 2.0 * f + fm) / (frac * frac))

    extra = max(0, int(tensor_extra_directions))
    if extra > 0:
        dn = max(1e-300, dnorm)
        for kk in range(extra):
            if kk == 0:
                qdir = 1j * delta1
            else:
                qdir = raw_direction(len(y), 0x301000 + 7919 * (kk + 1), 0x301777 + kk, True)
                qdir = np.asarray(qdir, dtype=np.complex128)
                qdir = qdir / max(1e-300, float(np.linalg.norm(qdir))) * dn
            hvec2 = frac * qdir
            fp2 = target.eval(y + hvec2)
            fm2 = target.eval(y - hvec2)
            sec2_terms.append((fp2 - 2.0 * f + fm2) / (frac * frac))

    sec2 = np.mean(np.asarray(sec2_terms, dtype=np.complex128), axis=0)
    delta2 = np.linalg.solve(Q, -0.5 * sec2)
    if not np.all(np.isfinite(delta2)):
        raise RuntimeError("nonfinite-cubic-delta2")

    d2norm = float(np.linalg.norm(delta2))
    max_d2 = float(max_correction) * max(dnorm, 1e-300)
    if d2norm > max_d2:
        delta2 = delta2 * (max_d2 / max(d2norm, 1e-300))
        d2norm = float(np.linalg.norm(delta2))
    if dnorm + d2norm > 20.0 * ynorm:
        scale = (20.0 * ynorm) / max(dnorm + d2norm, 1e-300)
        delta2 = delta2 * scale
        d2norm = float(np.linalg.norm(delta2))
    delta = delta1 + delta2
    if not np.all(np.isfinite(delta)):
        raise RuntimeError("nonfinite-cubic-delta")
    meta.update({
        "cubic_used": True,
        "cubic_delta2_norm": float(d2norm),
        "cubic_sec2_norm": float(np.linalg.norm(sec2)),
        "cubic_tensor_terms": int(len(sec2_terms)),
        "cubic_norm": float(np.linalg.norm(delta)),
    })
    return delta, meta


def pandrosion_corrector(
    target: TargetTrack,
    y0: Sequence[complex],
    max_epochs: int,
    tol: float,
    accept: float,
    trial_timeout: float,
    line_search: int = 12,
    probe_scale: float = 0.035,
    direction_seed: int = 0,
    probe_candidates: int = 8,
    probe_radii: Sequence[float] = (0.0, 0.5, 1.0, 2.0, 4.0),
    include_self_probe: bool = True,
    probe_condition_top: int = 2,
    probe_condition_weight: float = 0.0025,
    probe_axis: bool = True,
    probe_equal_value_weight: float = 0.015,
    probe_equal_value_eps: float = 1e-12,
    probe_curvature_top: int = 2,
    probe_curvature_weight: float = 0.003,
    probe_curvature_mid: float = 0.5,
    trust_cond_weight: float = 0.06,
    trust_cond_min: float = 4.0,
    line_grid: Optional[Sequence[float]] = None,
    halley_enabled: bool = True,
    halley_gate_residual: float = 0.25,
    halley_probe_fraction: float = 0.50,
    halley_cond_weight: float = 0.025,
    halley_min_gate: float = 0.04,
    halley_max_correction: float = 1.25,
    probe_log_weight: float = 0.0005,
    probe_log_delta_weight: float = 0.0015,
    halley_log_weight: float = 0.03,
    halley_log_scale: float = 1.0,
    tensor_extra_directions: int = 0,
) -> dict[str, Any]:
    """301 corrector: anchored divided-difference full cubic Halley only.

    The raw anchored finite-slope solve Q delta1=-G(a) is used only as the
    first ingredient of the cubic Halley formula.  The line search never tests
    or accepts the raw delta1 direction; it only tests delta1+delta2.
    """
    ensure_numpy()
    y = np.asarray(y0, dtype=np.complex128).copy()
    t0 = now()
    deadline = t0 + trial_timeout if trial_timeout and trial_timeout > 0 else None
    best_y = y.copy()
    best_r = finite_residual(target, y)
    ok = False
    status = "started"
    epochs = 0
    prev_delta = None
    last_cond = None
    last_probe_meta: dict[str, Any] = {}
    last_cubic_meta: dict[str, Any] = {}
    total_probe_evals = 0
    total_line_evals = 0
    total_cubic_evals = 0
    cubic_used_count = 0
    lambdas = _line_lambdas(line_search, line_grid)

    for ep in range(max(1, int(max_epochs))):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break
        try:
            f = target.eval(y)
            r = float(np.linalg.norm(f))
        except Exception as exc:
            status = f"eval-error:{type(exc).__name__}"
            break
        if math.isfinite(r) and r < best_r:
            best_r = r
            best_y = y.copy()
        if r <= max(float(tol), float(accept)) and (accept <= 0 or r < accept):
            ok = True
            status = "converged"
            break

        try:
            b, pmeta = _probe_endpoint_candidates(
                target=target, y=y, f=f, residual=r, prev_delta=prev_delta, ep=ep,
                direction_seed=direction_seed, probe_scale=float(probe_scale),
                probe_candidates=int(probe_candidates), probe_radii=probe_radii,
                include_self_probe=bool(include_self_probe), condition_top=int(probe_condition_top),
                condition_weight=float(probe_condition_weight), axis_probes=bool(probe_axis),
                equal_value_weight=float(probe_equal_value_weight), equal_value_eps=float(probe_equal_value_eps),
                curvature_top=int(probe_curvature_top), curvature_weight=float(probe_curvature_weight),
                curvature_mid=float(probe_curvature_mid), probe_log_weight=float(probe_log_weight),
                probe_log_delta_weight=float(probe_log_delta_weight),
            )
            last_probe_meta = pmeta
            total_probe_evals += int(pmeta.get("probe_evals", 0))
        except Exception as exc:
            status = f"probe-error:{type(exc).__name__}"
            break

        try:
            Q = target.slope_matrix(y, b)
            last_cond = float(np.linalg.cond(Q))
            delta1 = np.linalg.solve(Q, -f)
        except Exception as exc:
            status = f"anchored-slope-solve-error:{type(exc).__name__}"
            break
        if not np.all(np.isfinite(delta1)):
            status = "nonfinite-cubic-delta1"
            break

        ynorm = max(1.0, float(np.linalg.norm(y)))
        dnorm = float(np.linalg.norm(delta1))
        cap_factor = 18.0
        if float(trust_cond_weight) > 0 and last_cond is not None and math.isfinite(float(last_cond)):
            cap_factor = max(float(trust_cond_min), 18.0 / (1.0 + float(trust_cond_weight) * math.log1p(max(0.0, float(last_cond)))))
        if dnorm > cap_factor * ynorm:
            delta1 = delta1 * ((cap_factor * ynorm) / max(dnorm, 1e-300))

        try:
            delta, cmeta = _full_cubic_halley_delta_301(
                target=target, y=y, f=f, Q=Q, delta1=delta1,
                probe_fraction=float(halley_probe_fraction),
                max_correction=float(halley_max_correction),
                tensor_extra_directions=int(tensor_extra_directions),
            )
            last_cubic_meta = cmeta
            total_cubic_evals += 2 * int(cmeta.get("cubic_tensor_terms", 1))
            cubic_used_count += 1
        except Exception as exc:
            status = f"cubic-error:{type(exc).__name__}"
            break

        L = np.asarray(lambdas, dtype=float)
        Ytry = y[None, :] + L[:, None] * delta[None, :]
        Rtry = target.residuals_batch(Ytry)
        total_line_evals += int(len(Ytry))
        finite = np.isfinite(Rtry)
        better = finite & ((Rtry < r) | (Rtry < best_r))
        if not np.any(better):
            status = "no-cubic-decrease"
            break
        scores = np.where(better, Rtry * (1.0 + 1e-15 / np.maximum(L, 1e-15)), np.inf)
        idx = int(np.nanargmin(scores))
        lam = float(L[idx])
        rr = float(Rtry[idx])
        yy = Ytry[idx].copy()
        prev_delta = lam * delta
        y = yy
        if rr < best_r:
            best_y = yy.copy()
            best_r = rr
        epochs = ep + 1
        last_cubic_meta["cubic_chosen_lambda"] = lam
    else:
        status = "max-epochs"

    final_r = finite_residual(target, best_y)
    if final_r <= max(float(tol), float(accept)) and (accept <= 0 or final_r < accept):
        ok = True
        status = "converged"

    return {
        "accepted": bool(ok if accept <= 0 else (math.isfinite(final_r) and final_r < accept)),
        "ok": bool(ok),
        "status": status,
        "epochs": int(epochs),
        "residual": float(final_r),
        "y": best_y,
        "seconds": float(now() - t0),
        "slope_cond": last_cond,
        "corrector": "301-anchored-divided-difference-full-cubic-halley",
        "raw_pandrosion_step_accepted": False,
        "probe_total_evals": int(total_probe_evals),
        "line_search_batch_numpy": True,
        "line_search_evals": int(total_line_evals),
        "line_lambdas": [float(x) for x in lambdas],
        "trust_cond_weight": float(trust_cond_weight),
        "trust_cond_min": float(trust_cond_min),
        "halley_enabled": True,
        "halley_total_evals": int(total_cubic_evals),
        "halley_used_count": int(cubic_used_count),
        "halley_probe_fraction": float(halley_probe_fraction),
        "halley_max_correction": float(halley_max_correction),
        "tensor_extra_directions": int(tensor_extra_directions),
        **last_probe_meta,
        **last_cubic_meta,
    }


def cluster_index(roots: list[dict[str, Any]], z: Any, sep: float) -> Optional[int]:
    zz = np.asarray(z, dtype=np.complex128)
    zn = max(1.0, float(np.linalg.norm(zz)))
    for i, r in enumerate(roots):
        rz = np.asarray(r["z_complex"], dtype=np.complex128)
        dist = float(np.linalg.norm(zz - rz)) / max(zn, float(np.linalg.norm(rz)), 1.0)
        if dist <= sep:
            return i
    return None


def realness(z: Any) -> float:
    zz = np.asarray(z, dtype=np.complex128)
    return float(np.linalg.norm(zz.imag) / max(1e-300, np.linalg.norm(zz)))


def slope_condition_from_corrector(loc: dict[str, Any]) -> Optional[float]:
    c = loc.get("slope_cond")
    try:
        return float(c) if c is not None and math.isfinite(float(c)) else None
    except Exception:
        return None


def score_root(residual: float, realness_value: float, cond: Optional[float]) -> float:
    c = float(cond) if cond is not None and math.isfinite(float(cond)) else 1e300
    return float(math.log1p(max(0.0, residual)) + 0.01 * math.log1p(max(0.0, realness_value)) + 0.001 * math.log1p(max(0.0, c)))


# ---------------------------------------------------------------------------
# Main extractor
# ---------------------------------------------------------------------------

def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    ensure_numpy()
    t_case = now()
    n, d = parse_case(case_raw)
    system = DenseKostlanSystem.make(n, d, seed_index=int(args.seed_index), equation_normalize=bool(args.equation_normalize))
    # 200 keeps the 118/120 robust backbone; projective/homothetic flow remains the global geometry.
    chart = LinearChart.identity(n, scale=float(args.linear_scale))
    target = TargetTrack(system, chart)

    powers = sorted(set(round(float(x), 16) for x in parse_float_list(args.powers, DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in parse_float_list(args.angles, DEFAULT_ANGLES_DEG)]
    angles_deg = [math.degrees(x) for x in angles]
    radii = parse_float_list(args.rays, DEFAULT_RADII, positive=True)
    gains = parse_float_list(args.startopt_gains, DEFAULT_GAINS, positive=True)
    probe_radii = parse_float_list(args.probe_radii, [0.0, 0.35, 0.7, 1.0, 1.6, 2.6, 4.2], positive=False)
    probe_radii = [r for r in probe_radii if r >= 0] or [0.0, 1.0]

    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    duplicates = 0
    failures = 0
    t_extract = now()

    for trial in range(int(args.pool)):
        if len(roots) >= int(args.count):
            break
        y_raw, geom = universal_atlas_start(
            target, n, trial, system.seed + 0x113000, powers, angles, radii, float(args.power_cap),
            roots_found=len(roots), duplicates=duplicates, failures=failures, target_count=int(args.count),
            universal_cells=int(args.universal_cells), universal_shells=int(args.universal_shells),
            cell_probe_radius=float(args.cell_probe_radius), cell_descent_min=float(args.cell_descent_min),
            cell_equal_gap_min=float(args.cell_equal_gap_min), cell_log_max=float(args.cell_log_max),
            universal_cycle=bool(args.universal_cycle)
        )
        y0, smeta = startopt(target, y_raw, trial, system.seed + 0x112555, int(args.startopt_steps), int(args.startopt_candidates), gains, int(args.startopt_micro_epochs))
        loc = pandrosion_corrector(target, y0, max_epochs=int(args.epochs), tol=float(args.tol), accept=float(args.accept), trial_timeout=float(args.trial_timeout), line_search=int(args.line_search), probe_scale=float(getattr(args, "probe_scale", 0.035)), direction_seed=system.seed + 7919 * trial, probe_candidates=int(args.probe_candidates), probe_radii=probe_radii, include_self_probe=bool(args.probe_self), probe_condition_top=int(args.probe_condition_top), probe_condition_weight=float(args.probe_condition_weight), probe_axis=bool(args.probe_axis), probe_equal_value_weight=float(args.probe_equal_value_weight), probe_equal_value_eps=float(args.probe_equal_value_eps), probe_curvature_top=int(args.probe_curvature_top), probe_curvature_weight=float(args.probe_curvature_weight), probe_curvature_mid=float(args.probe_curvature_mid), trust_cond_weight=float(args.trust_cond_weight), trust_cond_min=float(args.trust_cond_min), line_grid=parse_float_list(args.line_grid, []), halley_enabled=bool(args.halley), halley_gate_residual=float(args.halley_gate_residual), halley_probe_fraction=float(args.halley_probe_fraction), halley_cond_weight=float(args.halley_cond_weight), halley_min_gate=float(args.halley_min_gate), halley_max_correction=float(args.halley_max_correction), probe_log_weight=float(args.probe_log_weight), probe_log_delta_weight=float(args.probe_log_delta_weight), halley_log_weight=float(args.halley_log_weight), halley_log_scale=float(args.halley_log_scale), tensor_extra_directions=int(args.tensor_extra_directions))
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
            "corrector": loc.get("corrector", "batch-probe-aware-pure-pandrosion-conditioned-exact-slope"),
            "slope_cond": loc.get("slope_cond"),
            "probe_mode": loc.get("probe_mode"),
            "probe_name": loc.get("probe_name"),
            "probe_candidates": loc.get("probe_candidates"),
            "probe_total_evals": loc.get("probe_total_evals"),
            "line_search_evals": loc.get("line_search_evals"),
            "probe_condition_best": loc.get("probe_condition_best"),
            "probe_condition_scored": loc.get("probe_condition_scored"),
            "probe_equal_value_gap": loc.get("probe_equal_value_gap"),
            "probe_equal_value_gap_min": loc.get("probe_equal_value_gap_min"),
            "probe_equal_value_penalty": loc.get("probe_equal_value_penalty"),
            "probe_curvature_best": loc.get("probe_curvature_best"),
            "probe_curvature_scored": loc.get("probe_curvature_scored"),
            "probe_axis_enabled": loc.get("probe_axis_enabled"),
            "probe_residual": loc.get("probe_residual"),
            "probe_distance": loc.get("probe_distance"),
                        "probe_improvement_proxy": loc.get("probe_improvement_proxy"),
            "halley_enabled": loc.get("halley_enabled"),
            "halley_gate": loc.get("halley_gate"),
            "halley_used": loc.get("halley_used"),
            "halley_used_count": loc.get("halley_used_count"),
            "halley_chosen_direction": loc.get("halley_chosen_direction"),
            "halley_delta2_norm": loc.get("halley_delta2_norm"),
            "halley_total_evals": loc.get("halley_total_evals"),
            **geom,
            **smeta,
        }
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
        rid = len(roots)
        cond = slope_condition_from_corrector(loc)
        realv = realness(z)
        roots.append({
            "id": rid,
            "source": "305-highdim-universal-atlas-full-cubic-halley",
            "trial": int(trial),
            "z_complex": np.asarray(z, dtype=np.complex128).copy(),
            "y_complex": np.asarray(loc["y"], dtype=np.complex128).copy(),
            "residual": float(r_orig),
            "realness": float(realv),
            "cond": cond,
            "score": score_root(float(r_orig), realv, cond),
            "epochs": int(loc.get("epochs", 0)),
            "seconds": float(loc.get("seconds", 0.0)),
            "corrector": loc.get("corrector", "batch-probe-aware-pure-pandrosion-conditioned-exact-slope"),
            "slope_cond": loc.get("slope_cond"),
            "probe_mode": loc.get("probe_mode"),
            "probe_name": loc.get("probe_name"),
            "probe_candidates": loc.get("probe_candidates"),
            "probe_total_evals": loc.get("probe_total_evals"),
            "line_search_evals": loc.get("line_search_evals"),
            "probe_condition_best": loc.get("probe_condition_best"),
            "probe_condition_scored": loc.get("probe_condition_scored"),
            "probe_equal_value_gap": loc.get("probe_equal_value_gap"),
            "probe_equal_value_gap_min": loc.get("probe_equal_value_gap_min"),
            "probe_equal_value_penalty": loc.get("probe_equal_value_penalty"),
            "probe_curvature_best": loc.get("probe_curvature_best"),
            "probe_curvature_scored": loc.get("probe_curvature_scored"),
            "probe_axis_enabled": loc.get("probe_axis_enabled"),
            "probe_residual": loc.get("probe_residual"),
            "probe_distance": loc.get("probe_distance"),
                        "probe_improvement_proxy": loc.get("probe_improvement_proxy"),
            "halley_enabled": loc.get("halley_enabled"),
            "halley_gate": loc.get("halley_gate"),
            "halley_used": loc.get("halley_used"),
            "halley_used_count": loc.get("halley_used_count"),
            "halley_chosen_direction": loc.get("halley_chosen_direction"),
            "halley_delta2_norm": loc.get("halley_delta2_norm"),
            "halley_total_evals": loc.get("halley_total_evals"),
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
        rr["z"] = root_to_json(zc)
        rr["y"] = root_to_json(yc)
        encoded_roots.append(rr)

    result = {
        "script": "305_pandrosion_highdim_universal_atlas_full_cubic_halley_numpy_engine.py",
        "autonomous": True,
        "dependencies": {"python_scripts": [], "numpy": bool(np is not None)},
        "mode": "305-highdim-universal-atlas-full-cubic-halley",
        "flow_formula": "u -> Mobius_Riemann(theta,pole)(u) -> y=Lambda*Mobius(u) -> StartOpt(y) -> anchored probe b -> Q_G(a,b) delta1=-G(a) -> symmetric finite-probe H(delta1,delta1) -> Q_G(a,b) delta2=-1/2 H -> full cubic delta=delta1+delta2 -> batched line search only on cubic direction; no raw Pandrosion-step acceptance, no analytic Newton/Jacobian/Hessian fallback",
        "case": f"{n},{d}",
        "family": "ks",
        "seed_index": int(args.seed_index),
        "seed": int(system.seed),
        "n": int(n),
        "degree": int(d),
        "terms_per_poly": system.terms_per_poly,
        "terms": system.total_terms,
        "bezout": system.bezout,
        "equation_normalize": bool(args.equation_normalize),
        "linear_A": [[cjson(chart.A[i, j]) for j in range(n)] for i in range(n)],
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
            "probe_candidates": int(args.probe_candidates),
            "probe_radii": probe_radii,
            "probe_self": bool(args.probe_self),
            "probe_axis": bool(args.probe_axis),
            "probe_condition_top": int(args.probe_condition_top),
            "probe_condition_weight": float(args.probe_condition_weight),
            "probe_equal_value_weight": float(args.probe_equal_value_weight),
            "probe_equal_value_eps": float(args.probe_equal_value_eps),
            "probe_curvature_top": int(args.probe_curvature_top),
            "probe_curvature_weight": float(args.probe_curvature_weight),
            "probe_curvature_mid": float(args.probe_curvature_mid),
            "trust_cond_weight": float(args.trust_cond_weight),
            "trust_cond_min": float(args.trust_cond_min),
            "line_grid": parse_float_list(args.line_grid, []),
            "halley": bool(args.halley),
            "halley_gate_residual": float(args.halley_gate_residual),
            "halley_probe_fraction": float(args.halley_probe_fraction),
            "halley_cond_weight": float(args.halley_cond_weight),
            "halley_min_gate": float(args.halley_min_gate),
            "halley_max_correction": float(args.halley_max_correction),
            "probe_log_weight": float(args.probe_log_weight),
            "probe_log_delta_weight": float(args.probe_log_delta_weight),
            "halley_log_weight": float(args.halley_log_weight),
            "halley_log_scale": float(args.halley_log_scale),
            "tensor_extra_directions": int(args.tensor_extra_directions),
            "powers": powers,
            "power_cap": float(args.power_cap),
            "angles_deg": angles_deg,
            "base_rays": radii,
            "startopt_steps": int(args.startopt_steps),
            "startopt_candidates": int(args.startopt_candidates),
            "startopt_gains": gains,
            "startopt_micro_epochs": int(args.startopt_micro_epochs),
            "universal_cells": int(args.universal_cells),
            "universal_shells": int(args.universal_shells),
            "cell_probe_radius": float(args.cell_probe_radius),
            "cell_descent_min": float(args.cell_descent_min),
            "cell_equal_gap_min": float(args.cell_equal_gap_min),
            "cell_log_max": float(args.cell_log_max),
            "universal_cycle": bool(args.universal_cycle),
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
            "extract_seconds": float(now() - t_extract),
            "total_seconds": float(now() - t_case),
            "eval_stats": system.stats(),
        },
    }
    return result


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="305 standalone NumPy Pandrosion engine: high-dimensional universal atlas plus anchored full cubic Halley.")
    p.add_argument("--cases", default="2,4", help="comma-separated case n,d; multiple cases can be separated by ';'")
    p.add_argument("--seed-index", type=int, default=0)
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--no-equation-normalize", dest="equation_normalize", action="store_false")
    p.add_argument("--linear-scale", type=float, default=1.0)
    p.add_argument("--count", "--thales-count", "--useful-count", type=int, default=8)
    p.add_argument("--pool", "--thales-pool", "--useful-pool", type=int, default=4096)
    p.add_argument("--epochs", "--thales-epochs", "--useful-epochs", type=int, default=24)
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--accept", "--residual-accept", type=float, default=1e-8)
    p.add_argument("--cluster-sep", type=float, default=1e-8)
    p.add_argument("--trial-timeout", "--thales-trial-timeout", "--useful-trial-timeout", type=float, default=0.0)
    p.add_argument("--line-search", type=int, default=12)
    p.add_argument("--probe-scale", type=float, default=0.035, help="base scale for finite-slope endpoint probes")
    p.add_argument("--probe-candidates", type=int, default=8, help="number of theorem-guided b-probes scored by ||F(b)|| per epoch")
    p.add_argument("--probe-radii", default="0,0.35,0.7,1,1.6,2.6,4.2", help="nonnegative multipliers for theorem-guided b-probes")
    p.add_argument("--probe-self", action="store_true", default=True, help="include b=a self-probe; Q is still computed by finite telescopic sums, not a Jacobian formula")
    p.add_argument("--no-probe-self", dest="probe_self", action="store_false")
    p.add_argument("--probe-axis", action="store_true", default=True, help="include coordinate-axis Thales probes in the finite-slope b-probe pool")
    p.add_argument("--no-probe-axis", dest="probe_axis", action="store_false")
    p.add_argument("--probe-condition-top", type=int, default=2, help="condition-number score for this many best residual b-probes; 0 disables")
    p.add_argument("--probe-condition-weight", type=float, default=0.0025, help="small log-cond penalty in b-probe scoring")
    p.add_argument("--probe-equal-value-weight", type=float, default=0.015, help="penalty for non-self probes with F(b) close to F(a), i.e. near equal-value poles")
    p.add_argument("--probe-equal-value-eps", type=float, default=1e-12, help="floor for equal-value pole gap")
    p.add_argument("--probe-curvature-top", type=int, default=2, help="score finite-slope curvature proxy on this many best probes; 0 disables")
    p.add_argument("--probe-curvature-weight", type=float, default=0.003, help="small penalty for finite-slope curvature variation along probe segment")
    p.add_argument("--probe-curvature-mid", type=float, default=0.5, help="midpoint fraction used in finite-slope curvature proxy")
    p.add_argument("--trust-cond-weight", type=float, default=0.06, help="condition-aware trust-region shrinkage; 0 disables")
    p.add_argument("--trust-cond-min", type=float, default=4.0, help="minimum trust-region factor in units of max(1,||y||)")
    p.add_argument("--line-grid", default="1,0.75,0.5,0.35,0.25,0.18,0.125,0.09,0.0625,0.045,0.03125,0.02", help="comma-separated lambdas for batched line search")
    p.add_argument("--halley", action="store_true", default=True, help="301: kept for CLI compatibility; full cubic Halley is always on")
    p.add_argument("--no-halley", dest="halley", action="store_false", help="301: ignored; raw Pandrosion fallback is intentionally not used")
    p.add_argument("--halley-gate-residual", type=float, default=0.25, help="301: kept for CLI compatibility; no gate is used")
    p.add_argument("--halley-probe-fraction", type=float, default=0.50, help="301: symmetric second-slope probe fraction along cubic delta1")
    p.add_argument("--halley-cond-weight", type=float, default=0.025, help="301: kept for CLI compatibility; no gate is used")
    p.add_argument("--halley-min-gate", type=float, default=0.04, help="301: kept for CLI compatibility; no gate is used")
    p.add_argument("--halley-max-correction", type=float, default=1.25, help="301: cap ||delta2|| relative to ||delta1||")
    p.add_argument("--probe-log-weight", type=float, default=0.0005, help="301: penalty for absolute logarithmic/projective scale energy of probes")
    p.add_argument("--probe-log-delta-weight", type=float, default=0.0015, help="301: penalty for probes that distort the current logarithmic scale")
    p.add_argument("--halley-log-weight", type=float, default=0.03, help="301: logarithmic stability penalty in the Halley gate")
    p.add_argument("--halley-log-scale", type=float, default=1.0, help="301: scale for logarithmic stability gate")
    p.add_argument("--tensor-extra-directions", type=int, default=0, help="301: extra symmetric finite-slope tensor directions near a root")
    p.add_argument("--self-test", action="store_true", help="run a small ks(2,2) smoke test and exit")

    p.add_argument("--powers", "--thales-powers", default=None)
    p.add_argument("--power-cap", "--thales-power-cap", type=float, default=1048576.0)
    p.add_argument("--angles", "--thales-angles", default=None, help="degrees, comma-separated")
    p.add_argument("--rays", "--thales-rays", default=None)
    p.add_argument("--startopt-steps", type=int, default=1)
    p.add_argument("--startopt-candidates", type=int, default=12)
    p.add_argument("--startopt-gains", default=None)
    p.add_argument("--startopt-micro-epochs", type=int, default=0)
    p.add_argument("--universal-cells", type=int, default=16, help="305: fixed universal atlas cells tested per trial")
    p.add_argument("--universal-shells", type=int, default=5, help="305: fixed shell layers for simplex/hypercube/projective atlas")
    p.add_argument("--cell-probe-radius", type=float, default=0.14, help="305: fixed local cell probe radius")
    p.add_argument("--cell-descent-min", type=float, default=1.02, help="305: required local descent ratio for admissible fixed cell")
    p.add_argument("--cell-equal-gap-min", type=float, default=1e-10, help="305: minimum equal-value gap for admissible fixed cell")
    p.add_argument("--cell-log-max", type=float, default=80.0, help="305: maximum log-scale energy for admissible fixed cell")
    p.add_argument("--universal-cycle", action="store_true", default=True, help="305: cycle through admissible universal cells")
    p.add_argument("--no-universal-cycle", dest="universal_cycle", action="store_false")
    p.add_argument("--out", "--thales-out", "--useful-out", default=None)
    p.add_argument("--outdir", default="/mnt/data/305_highdim_universal_atlas_out")
    p.add_argument("--keep-trials", type=int, default=160)
    p.add_argument("--verbose-trials", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(getattr(args, "self_test", False)):
        args.cases = "2,2"
        args.count = 2
        args.pool = min(int(args.pool), 128)
        args.epochs = min(int(args.epochs), 16)
        args.accept = min(float(args.accept), 1e-8)
        args.keep_trials = min(int(args.keep_trials), 20)
        args.out = args.out or "/mnt/data/305_highdim_universal_atlas_out/self_test_305.json"
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    outputs = [run_case(args, c) for c in cases]
    final: dict[str, Any]
    if len(outputs) == 1:
        final = outputs[0]
    else:
        final = {"script": "305_pandrosion_highdim_universal_atlas_full_cubic_halley_numpy_engine.py", "autonomous": True, "cases": outputs}
    if args.out:
        out = Path(args.out)
    else:
        first = cases[0].replace(",", "x") if cases else "case"
        out = Path(args.outdir) / f"305_highdim_universal_atlas_{first}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2), encoding="utf-8")

    print("=" * 120, flush=True)
    print("305 autonomous HIGH-DIM UNIVERSAL ATLAS + ANCHORED FULL-CUBIC HALLEY NumPy Pandrosion", flush=True)
    print("No dependency on previous Python scripts; anchored divided-difference Q plus full cubic finite-slope Halley only; no raw Pandrosion-step acceptance, no Newton/Jacobian fallback path.", flush=True)
    print("=" * 120, flush=True)
    for r in outputs:
        s = r["summary"]
        print(f"case=ks({r['n']},{r['degree']}), seed={r['seed']}, terms={r['terms']}, Bezout={r['bezout']}", flush=True)
        print(f"roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} trials={s['trials_used']} duplicates={s['duplicates']} failures={s['failures']}", flush=True)
        print(f"seconds: generation={s['generation_seconds']:.2f}, extract={s['extract_seconds']:.2f}, total={s['total_seconds']:.2f}", flush=True)
        if r.get("roots"):
            best = r["roots"][0]
            print(f"best_root: residual={float(best.get('residual', float('inf'))):.3e}, trial={best.get('trial')}, Lambda={best.get('homothety')}, theta={best.get('theta_deg')}, startopt_ratio={best.get('startopt_ratio')}", flush=True)
    print(f"out={out}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    finally:
        try:
            sys.stdout.flush(); sys.stderr.flush()
        except Exception:
            pass
