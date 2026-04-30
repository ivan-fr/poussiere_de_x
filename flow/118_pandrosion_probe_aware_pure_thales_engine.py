#!/usr/bin/env python3
"""
116_pandrosion_probe_aware_pure_thales_engine.py

Autonomous PROBE-AWARE vectorized PURE Pandrosion Thales root extractor.

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
      -> correct G(y)=F(Ay)=0 with Q_G(a,b), not Jacobian
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

        Vectorization in 116:
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

    def slope_matrix(self, a_y: Sequence[complex], b_y: Sequence[complex]) -> "np.ndarray":
        a_z = self.chart.z_from_y(a_y)
        b_z = self.chart.z_from_y(b_y)
        return self.system.slope_matrix(a_z, b_z) @ self.chart.A

    def residual(self, y: Sequence[complex]) -> float:
        try:
            return float(np.linalg.norm(self.eval(y)))
        except Exception:
            return float("inf")


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


def startopt(target: TargetTrack, y0: Any, trial: int, seed: int, steps: int, candidates: int, gains: Sequence[float], micro_epochs: int) -> tuple[Any, dict[str, Any]]:
    ensure_numpy()
    best = np.asarray(y0, dtype=np.complex128).copy()
    best_r = finite_residual(target, best)
    initial = best_r
    evals = 1
    micro_total = 0
    chosen_gain = 1.0
    for step in range(max(0, steps)):
        base = best.copy()
        pool = [(1.0, base)]
        for c in range(max(0, candidates - 1)):
            gain = float(gains[(trial + 3 * step + c) % len(gains)])
            cand = gain * base
            if c % 3 == 1:
                # phase/amplitude jitter inside the same chart
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
            pool.append((gain, cand))
        for gain, cand in pool:
            evals += 1
            yscore = cand
            r = finite_residual(target, yscore)
            if micro_epochs > 0:
                loc = pandrosion_corrector(target, yscore, max_epochs=micro_epochs, tol=1e-12, accept=0.0, trial_timeout=0.0, line_search=6)
                micro_total += int(loc.get("epochs", 0))
                if float(loc.get("residual", float("inf"))) < r:
                    yscore = loc["y"]
                    r = float(loc["residual"])
            # residual dominates; norm is only tie-breaking
            score = math.log1p(max(0.0, r)) + 1e-15 * math.log1p(float(np.linalg.norm(yscore))) if math.isfinite(r) else float("inf")
            old = math.log1p(max(0.0, best_r)) + 1e-15 * math.log1p(float(np.linalg.norm(best))) if math.isfinite(best_r) else float("inf")
            if score < old:
                best = np.asarray(yscore, dtype=np.complex128)
                best_r = float(r)
                chosen_gain = float(gain)
    return best, {
        "startopt_enabled": bool(steps > 0),
        "startopt_r0": float(initial),
        "startopt_r1": float(best_r),
        "startopt_ratio": (float(best_r / initial) if math.isfinite(best_r) and math.isfinite(initial) and initial > 0 else None),
        "startopt_steps": int(max(0, steps)),
        "startopt_evals": int(evals),
        "startopt_micro_epochs": int(micro_total),
        "startopt_gain": float(chosen_gain),
    }


# ---------------------------------------------------------------------------
# Corrector and root handling
# ---------------------------------------------------------------------------

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
) -> tuple[Any, dict[str, Any]]:
    """Choose the finite-slope probe b by the theorem-guided rule.

    The finite-slope theorem says that the contraction constant is proportional
    to ||b-zeta||.  Since zeta is unknown, 116 uses the residual ||G(b)|| as a
    computable proxy and selects, among one universal geometric probe family,

        b_* = argmin_b ||G(b)||.

    This does not change the solver architecture: every candidate is produced by
    the same finite Thales geometry around the current point.  There is no
    Newton fallback, no least-squares fallback, no homotopy path, and no
    bivariate special case.  The optional self-probe b=a is also evaluated by
    the same telescopic finite-slope formula; no Jacobian formula is formed.
    """
    ensure_numpy()
    n = len(y)
    ynorm = max(1.0, float(np.linalg.norm(y)))
    radii = [float(r) for r in probe_radii if float(r) >= 0]
    if not radii:
        radii = [1.0]
    candidates: list[tuple[str, Any]] = []
    if include_self_probe:
        candidates.append(("self", y.copy()))
    # Inertial probe from the previous accepted Pandrosion displacement.
    if prev_delta is not None and np.all(np.isfinite(prev_delta)):
        pdn = max(1e-300, float(np.linalg.norm(prev_delta)))
        base = prev_delta / pdn * min(max(pdn, float(probe_scale) * ynorm), 2.5 * ynorm)
        candidates.append(("inertial", y + base))
    # Deterministic residual-aware geometric probes.
    # The directions are complex and generic in dimension n.  The same rule is
    # used for every n, so this remains fully multivariate.
    budget = max(1, int(probe_candidates))
    k = 0
    while len(candidates) < budget:
        rad = float(probe_scale) * ynorm * radii[k % len(radii)]
        qdir = raw_direction(n, direction_seed + 104729 * (ep + 1) + 7919 * (k + 1), direction_seed ^ (0x116116 + 17 * k), True)
        qnorm = max(1e-300, float(np.linalg.norm(qdir)))
        qdir = qdir / qnorm * math.sqrt(max(1, n))
        # The phase spiral avoids repeated collinear probes while remaining a
        # single closed formula.
        ph = phase(0.6180339887498948 * (ep + 1) + 2.399963229728653 * (k + 1))
        step = rad * ph * qdir
        # Mix in a small radial component when y is nonzero; this reflects the
        # Thales start-optimization idea and helps make b a better root proxy.
        if float(np.linalg.norm(y)) > 0:
            step = step + (0.12 * rad) * y / ynorm * phase(0.38196601125 * (k + 1))
        # Avoid collapsed coordinates numerically; same formula, regularized.
        tiny = 1e-12 * ynorm
        for j in range(n):
            if abs(step[j]) < tiny:
                step[j] += tiny * phase(0.17 + j + ep + k)
        candidates.append((f"geom-{k}", y + step))
        k += 1

    best_name = ""
    best_b = None
    best_res = float("inf")
    evals = 0
    best_distance = 0.0
    for name, b in candidates[:budget]:
        try:
            rb = finite_residual(target, b)
        except Exception:
            rb = float("inf")
        evals += 1
        dist = float(np.linalg.norm(np.asarray(b, dtype=np.complex128) - y))
        # Residual is the theorem-guided proxy for ||b-zeta||.  Distance is only
        # a mild tie-breaker to avoid huge probes when residuals match.
        score = math.log1p(max(0.0, rb)) + 1e-14 * math.log1p(dist)
        old = math.log1p(max(0.0, best_res)) + 1e-14 * math.log1p(best_distance)
        if math.isfinite(score) and score < old:
            best_name = name
            best_b = np.asarray(b, dtype=np.complex128).copy()
            best_res = float(rb)
            best_distance = float(dist)
    if best_b is None:
        # Not a fallback algorithm: if no finite probe exists, the trial fails.
        raise RuntimeError("no-finite-probe")
    return best_b, {
        "probe_mode": "theorem-guided-residual-min",
        "probe_name": best_name,
        "probe_candidates": int(min(budget, len(candidates))),
        "probe_evals": int(evals),
        "probe_residual": float(best_res),
        "probe_distance": float(best_distance),
        "probe_improvement_proxy": (float(residual / best_res) if math.isfinite(best_res) and best_res > 0 and math.isfinite(residual) else None),
        "probe_self_enabled": bool(include_self_probe),
    }


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
) -> dict[str, Any]:
    """Probe-aware pure derivative-free Pandrosion corrector.

    Compared with 115, the finite-slope endpoint b is no longer a single
    predetermined perturbation.  The finite-slope theorem shows that the local
    contraction is controlled by the distance from b to the root.  Since the root
    is unknown, 116 chooses b by minimizing ||G(b)|| over a small universal family
    of geometric probes.

        b_* = argmin ||G(b)||,
        Q_G(a,b_*) delta = -G(a),
        a_next = a + lambda delta.

    No Jacobian is built, no derivative formula is called, and no alternative
    solver is used.  The only linear solve is the Pandrosion solve in the exact
    telescopic finite-slope matrix Q_G(a,b_*).
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
    total_probe_evals = 0
    n = len(y)
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
                target=target,
                y=y,
                f=f,
                residual=r,
                prev_delta=prev_delta,
                ep=ep,
                direction_seed=direction_seed,
                probe_scale=float(probe_scale),
                probe_candidates=int(probe_candidates),
                probe_radii=probe_radii,
                include_self_probe=bool(include_self_probe),
            )
            last_probe_meta = pmeta
            total_probe_evals += int(pmeta.get("probe_evals", 0))
        except Exception as exc:
            status = f"probe-error:{type(exc).__name__}"
            break
        try:
            Q = target.slope_matrix(y, b)
            last_cond = float(np.linalg.cond(Q))
            delta = np.linalg.solve(Q, -f)
        except Exception as exc:
            status = f"slope-solve-error:{type(exc).__name__}"
            break
        if not np.all(np.isfinite(delta)):
            status = "nonfinite-step"
            break
        ynorm = max(1.0, float(np.linalg.norm(y)))
        dnorm = float(np.linalg.norm(delta))
        if dnorm > 18.0 * ynorm:
            delta = delta * ((18.0 * ynorm) / max(dnorm, 1e-300))
        accepted = False
        base_r = r
        for k in range(max(1, int(line_search))):
            lam = 1.0 / (2.0 ** k)
            yy = y + lam * delta
            rr = finite_residual(target, yy)
            if math.isfinite(rr) and (rr < base_r or rr < best_r):
                prev_delta = lam * delta
                y = yy
                if rr < best_r:
                    best_y = yy.copy(); best_r = rr
                accepted = True
                break
        epochs = ep + 1
        if not accepted:
            status = "no-decrease"
            break
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
        "corrector": "probe-aware-pure-pandrosion-exact-telescopic-slope",
        "probe_total_evals": int(total_probe_evals),
        **last_probe_meta,
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
    # Fully autonomous 116 keeps A simple and robust; the projective/homothetic flow provides the main geometry.
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
        y_raw, geom = mobius_homothety_start(
            n, trial, system.seed + 0x113000, powers, angles, radii, float(args.power_cap),
            roots_found=len(roots), duplicates=duplicates, failures=failures, target_count=int(args.count)
        )
        y0, smeta = startopt(target, y_raw, trial, system.seed + 0x112555, int(args.startopt_steps), int(args.startopt_candidates), gains, int(args.startopt_micro_epochs))
        loc = pandrosion_corrector(target, y0, max_epochs=int(args.epochs), tol=float(args.tol), accept=float(args.accept), trial_timeout=float(args.trial_timeout), line_search=int(args.line_search), probe_scale=float(getattr(args, "probe_scale", 0.035)), direction_seed=system.seed + 7919 * trial, probe_candidates=int(args.probe_candidates), probe_radii=probe_radii, include_self_probe=bool(args.probe_self))
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
            "probe_mode": loc.get("probe_mode"),
            "probe_name": loc.get("probe_name"),
            "probe_candidates": loc.get("probe_candidates"),
            "probe_total_evals": loc.get("probe_total_evals"),
            "probe_residual": loc.get("probe_residual"),
            "probe_distance": loc.get("probe_distance"),
            "probe_improvement_proxy": loc.get("probe_improvement_proxy"),
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
            "source": "116-autonomous-probe-aware-vectorized-pure-pandrosion-single-flow",
            "trial": int(trial),
            "z_complex": np.asarray(z, dtype=np.complex128).copy(),
            "y_complex": np.asarray(loc["y"], dtype=np.complex128).copy(),
            "residual": float(r_orig),
            "realness": float(realv),
            "cond": cond,
            "score": score_root(float(r_orig), realv, cond),
            "epochs": int(loc.get("epochs", 0)),
            "seconds": float(loc.get("seconds", 0.0)),
            "corrector": loc.get("corrector", "pure-pandrosion-exact-telescopic-slope"),
            "slope_cond": loc.get("slope_cond"),
            "probe_mode": loc.get("probe_mode"),
            "probe_name": loc.get("probe_name"),
            "probe_candidates": loc.get("probe_candidates"),
            "probe_total_evals": loc.get("probe_total_evals"),
            "probe_residual": loc.get("probe_residual"),
            "probe_distance": loc.get("probe_distance"),
            "probe_improvement_proxy": loc.get("probe_improvement_proxy"),
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
        "script": "116_pandrosion_probe_aware_pure_thales_engine.py",
        "autonomous": True,
        "dependencies": {"python_scripts": [], "numpy": bool(np is not None)},
        "mode": "autonomous-probe-aware-vectorized-pure-pandrosion-single-flow/thales-power-riemann-startopt-exact-slope",
        "flow_formula": "u -> Mobius_Riemann(theta,pole)(u) -> y=Lambda*Mobius(u) -> StartOpt(y) -> z=A*y -> PURE Pandrosion Q_G(a,b) on F(Ay), no Jacobian",
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
            "powers": powers,
            "power_cap": float(args.power_cap),
            "angles_deg": angles_deg,
            "base_rays": radii,
            "startopt_steps": int(args.startopt_steps),
            "startopt_candidates": int(args.startopt_candidates),
            "startopt_gains": gains,
            "startopt_micro_epochs": int(args.startopt_micro_epochs),
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
    p = argparse.ArgumentParser(description="Autonomous vectorized pure Thales Pandrosion engine: one flow, general multivariate, no solver fallback.")
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
    p.add_argument("--powers", "--thales-powers", default=None)
    p.add_argument("--power-cap", "--thales-power-cap", type=float, default=1048576.0)
    p.add_argument("--angles", "--thales-angles", default=None, help="degrees, comma-separated")
    p.add_argument("--rays", "--thales-rays", default=None)
    p.add_argument("--startopt-steps", type=int, default=1)
    p.add_argument("--startopt-candidates", type=int, default=12)
    p.add_argument("--startopt-gains", default=None)
    p.add_argument("--startopt-micro-epochs", type=int, default=0)
    p.add_argument("--out", "--thales-out", "--useful-out", default=None)
    p.add_argument("--outdir", default="/mnt/data/116_probe_aware_pure_out")
    p.add_argument("--keep-trials", type=int, default=160)
    p.add_argument("--verbose-trials", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    outputs = [run_case(args, c) for c in cases]
    final: dict[str, Any]
    if len(outputs) == 1:
        final = outputs[0]
    else:
        final = {"script": "116_pandrosion_probe_aware_pure_thales_engine.py", "autonomous": True, "cases": outputs}
    if args.out:
        out = Path(args.out)
    else:
        first = cases[0].replace(",", "x") if cases else "case"
        out = Path(args.outdir) / f"116_probe_aware_pure_{first}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2), encoding="utf-8")

    print("=" * 120, flush=True)
    print("116 autonomous PROBE-AWARE VECTORIZED PURE Thales Pandrosion - single flow, no fallback", flush=True)
    print("No dependency on previous Python scripts; theorem-guided b-probes + vectorized pure exact-slope Q(a,b); no Newton/Jacobian fallback path.", flush=True)
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
