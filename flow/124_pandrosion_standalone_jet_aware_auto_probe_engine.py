"""
124_pandrosion_standalone_jet_aware_auto_probe_engine.py

Standalone FULL RIEMANN-CHART + JET-AWARE AUTO-PROBE Pandrosion engine.

This is the clean 124 release.  It is standalone like 123, but adds a
higher-order jet-aware local inverse-chart candidate.  The projective
Riemann/Mobius chart still handles global normalization; the new jet layer
uses analytic first/second directional jets of G(y)=F(Phi(y)) to propose a
locally flattened inverse-chart step, then keeps it only if residual line search
confirms improvement.  The exact finite-slope Pandrosion correction remains the
main correction and no previous engine is imported.

No imports from previous Pandrosion scripts, no hidden engine dependency. This file contains its own:
  - Kostlan/dense polynomial generator for ks(n,d)
  - dense F evaluation engine + vectorized exact telescopic slope Q(a,b)
  - explicit coordinate-wise Riemann/Mobius chart layer z = Phi(y)
  - inversion/homothety/return map handled by the chart
  - local PURE Pandrosion corrector Q_G(a,b_*) delta = -G(a)
  - structured nD probe geometries: star, simplex, axes, hypercube, simplex-cube, balanced, adaptive
  - Riemann-aware automatic probe-frame selection
  - optional jet-aware local inverse-chart predictor (order 0/1/2) scored by residual

Dependencies: Python stdlib + NumPy only.
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
            "NumPy is required by this standalone high-degree extractor. "
            f"Import error: {_NUMPY_IMPORT_ERROR!r}"
        )


def standalone_dependency_audit() -> dict[str, Any]:
    """Runtime certificate: this file has loaded no previous Pandrosion module."""
    forbidden = ("118_pandrosion", "119_pandrosion", "120_pandrosion", "121_pandrosion", "122_pandrosion")
    loaded = sorted(name for name in sys.modules if any(tok in name for tok in forbidden))
    return {
        "standalone": True,
        "python_script_dependencies": [],
        "forbidden_pandrosion_modules_loaded": loaded,
        "numpy": bool(np is not None),
    }


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

    def jacobian(self, z: Sequence[complex]) -> "np.ndarray":
        """Analytic dense Jacobian of the generated polynomial system.

        124 uses this only for the optional jet-aware inverse-chart candidate.
        The core Pandrosion finite-slope matrix remains the exact telescopic
        divided difference.  This routine is standalone and derivative formulas
        are computed directly from the stored monomial exponents.
        """
        ensure_numpy()
        t0 = now()
        zz = np.asarray(z, dtype=np.complex128)
        pows = self._powers(zz)
        M = self.terms_per_poly
        n = self.n
        J = np.empty((n, n), dtype=np.complex128)
        for j in range(n):
            term = np.ones(M, dtype=np.complex128)
            for k in range(n):
                ek = self.exps[:, k].astype(np.int64)
                if k == j:
                    fac = np.zeros(M, dtype=np.complex128)
                    mask = ek > 0
                    if np.any(mask):
                        fac[mask] = ek[mask].astype(np.float64) * pows[k][ek[mask] - 1]
                else:
                    fac = pows[k][ek]
                term *= fac
            J[:, j] = self.coeff @ term
        self.jac_count = int(getattr(self, "jac_count", 0)) + 1
        self.seconds_jac = float(getattr(self, "seconds_jac", 0.0)) + (now() - t0)
        return J

    def hessian_vec(self, z: Sequence[complex], v: Sequence[complex]) -> "np.ndarray":
        """Second directional derivative D^2F(z)[v,v].

        The implementation multiplies each monomial along the line z+t v and
        keeps only coefficients through t^2.  It avoids explicit n^2 Hessian
        materialization and remains practical for dense Kostlan systems.
        """
        ensure_numpy()
        t0 = now()
        zz = np.asarray(z, dtype=np.complex128)
        vv = np.asarray(v, dtype=np.complex128)
        pows = self._powers(zz)
        M = self.terms_per_poly
        prod0 = np.ones(M, dtype=np.complex128)
        prod1 = np.zeros(M, dtype=np.complex128)
        prod2 = np.zeros(M, dtype=np.complex128)
        for j in range(self.n):
            e = self.exps[:, j].astype(np.int64)
            a0 = pows[j][e]
            a1 = np.zeros(M, dtype=np.complex128)
            a2 = np.zeros(M, dtype=np.complex128)
            m1 = e > 0
            if np.any(m1):
                a1[m1] = e[m1].astype(np.float64) * pows[j][e[m1] - 1] * vv[j]
            m2 = e > 1
            if np.any(m2):
                a2[m2] = 0.5 * e[m2].astype(np.float64) * (e[m2] - 1).astype(np.float64) * pows[j][e[m2] - 2] * (vv[j] ** 2)
            new2 = prod2 * a0 + prod1 * a1 + prod0 * a2
            new1 = prod1 * a0 + prod0 * a1
            new0 = prod0 * a0
            prod0, prod1, prod2 = new0, new1, new2
        out = self.coeff @ (2.0 * prod2)
        self.hess_count = int(getattr(self, "hess_count", 0)) + 1
        self.seconds_hess = float(getattr(self, "seconds_hess", 0.0)) + (now() - t0)
        return out

    def stats(self) -> dict[str, Any]:
        return {
            "eval_count": int(self.eval_count),
            "slope_count": int(getattr(self, "slope_count", 0)),
            "jac_count": int(getattr(self, "jac_count", 0)),
            "hess_count": int(getattr(self, "hess_count", 0)),
            "seconds_eval": float(self.seconds_eval),
            "seconds_slope": float(getattr(self, "seconds_slope", 0.0)),
            "seconds_jac": float(getattr(self, "seconds_jac", 0.0)),
            "seconds_hess": float(getattr(self, "seconds_hess", 0.0)),
            "terms_per_poly": self.terms_per_poly,
            "total_terms": self.total_terms,
        }



@dataclasses.dataclass
class ProjectiveMobiusChart:
    """Coordinate-wise Riemann/Mobius chart used by 122.

    Each coordinate is mapped by

        z_j = Lambda * pole_j * (c_j*w_j + s_j) / (-s_j*w_j + c_j),
        w_j = y_j / pole_j,
        c_j = cos(theta_j), s_j = sin(theta_j).

    This is an explicit chart layer: the corrector runs in y-coordinates on
    G(y)=F(Phi(y)), then returns to the original variables by z=Phi(y).

    Special cases:
      theta = 0       -> affine homothety z = Lambda*y
      theta = pi/2    -> inversion chart z = -Lambda*pole^2/y
      general theta   -> Riemann/Mobius interpolation between affine and infinity
    """
    scale: complex
    poles: Any
    theta: Any
    denom_eps: float = 1e-12

    @classmethod
    def affine(cls, n: int, scale: float = 1.0) -> "ProjectiveMobiusChart":
        ensure_numpy()
        return cls(
            scale=complex(scale),
            poles=np.ones(n, dtype=np.complex128),
            theta=np.zeros(n, dtype=np.float64),
        )

    @property
    def n(self) -> int:
        return int(len(self.poles))

    def _safe_denom(self, y: Sequence[complex]) -> tuple[Any, Any, Any, Any]:
        yy = np.asarray(y, dtype=np.complex128)
        poles = np.asarray(self.poles, dtype=np.complex128)
        theta = np.asarray(self.theta, dtype=np.float64)
        c = np.cos(theta)
        s = np.sin(theta)
        w = yy / poles
        denom = -s * w + c
        small = np.abs(denom) < float(self.denom_eps)
        if np.any(small):
            # Stay in the same projective chart; only regularize the numerical
            # representation of a point too close to the chart pole.
            jitter = np.asarray([phase(0.173 + 0.619 * j) for j in range(len(yy))], dtype=np.complex128)
            denom = denom + small.astype(np.float64) * float(self.denom_eps) * jitter
        return w, c, s, denom

    def z_from_y(self, y: Sequence[complex]) -> "np.ndarray":
        ensure_numpy()
        yy = np.asarray(y, dtype=np.complex128)
        poles = np.asarray(self.poles, dtype=np.complex128)
        w, c, s, denom = self._safe_denom(yy)
        return complex(self.scale) * poles * ((c * w + s) / denom)

    def y_from_z(self, z: Sequence[complex]) -> "np.ndarray":
        ensure_numpy()
        zz = np.asarray(z, dtype=np.complex128)
        poles = np.asarray(self.poles, dtype=np.complex128)
        theta = np.asarray(self.theta, dtype=np.float64)
        c = np.cos(theta)
        s = np.sin(theta)
        t = zz / (complex(self.scale) * poles)
        denom = t * s + c
        small = np.abs(denom) < float(self.denom_eps)
        if np.any(small):
            jitter = np.asarray([phase(0.341 + 0.517 * j) for j in range(len(zz))], dtype=np.complex128)
            denom = denom + small.astype(np.float64) * float(self.denom_eps) * jitter
        w = (t * c - s) / denom
        return poles * w

    def slope_diag(self, a_y: Sequence[complex], b_y: Sequence[complex]) -> "np.ndarray":
        """Finite chart slope D so Phi(b)-Phi(a)=D@(b-a).

        Since Phi is coordinate-wise, D is diagonal.  For nonzero coordinate
        displacements we use the exact finite quotient.  For a collapsed
        coordinate we use the analytic local slope of the same Mobius chart,
        Lambda/(-s*w+c)^2, only as the removable-limit value of the finite
        quotient.
        """
        ensure_numpy()
        a = np.asarray(a_y, dtype=np.complex128)
        b = np.asarray(b_y, dtype=np.complex128)
        za = self.z_from_y(a)
        zb = self.z_from_y(b)
        dy = b - a
        out = np.empty(len(a), dtype=np.complex128)
        mask = np.abs(dy) > 1e-14 * np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
        out[mask] = (zb[mask] - za[mask]) / dy[mask]
        if np.any(~mask):
            w, c, s, denom = self._safe_denom(a)
            deriv = complex(self.scale) / (denom * denom)
            out[~mask] = deriv[~mask]
        return out

    def derivative_diag(self, y: Sequence[complex]) -> "np.ndarray":
        """Coordinate-wise analytic derivative Phi'(y)."""
        ensure_numpy()
        _w, _c, _s, denom = self._safe_denom(y)
        return complex(self.scale) / (denom * denom)

    def second_derivative_diag(self, y: Sequence[complex]) -> "np.ndarray":
        """Coordinate-wise analytic second derivative Phi''(y)."""
        ensure_numpy()
        poles = np.asarray(self.poles, dtype=np.complex128)
        _w, _c, s, denom = self._safe_denom(y)
        return 2.0 * complex(self.scale) * s / (poles * denom * denom * denom)

    def meta(self) -> dict[str, Any]:
        theta = np.asarray(self.theta, dtype=np.float64)
        poles = np.asarray(self.poles, dtype=np.complex128)
        pole_angles = [float(math.degrees(cmath.phase(complex(p)))) for p in poles[: min(4, len(poles))]]
        return {
            "chart_layer": "coordinate-projective-riemann-mobius",
            "chart_scale_abs": float(abs(complex(self.scale))),
            "chart_scale": cjson(complex(self.scale)),
            "chart_theta_mean_deg": float(np.mean(np.degrees(theta))) if len(theta) else 0.0,
            "chart_theta_min_deg": float(np.min(np.degrees(theta))) if len(theta) else 0.0,
            "chart_theta_max_deg": float(np.max(np.degrees(theta))) if len(theta) else 0.0,
            "chart_pole_angle_deg_first": pole_angles,
        }


@dataclasses.dataclass
class TargetTrack:
    system: DenseKostlanSystem
    chart: ProjectiveMobiusChart

    def eval(self, y: Sequence[complex]) -> "np.ndarray":
        z = self.chart.z_from_y(y)
        return self.system.eval(z)

    def slope_matrix(self, a_y: Sequence[complex], b_y: Sequence[complex]) -> "np.ndarray":
        a_z = self.chart.z_from_y(a_y)
        b_z = self.chart.z_from_y(b_y)
        qz = self.system.slope_matrix(a_z, b_z)
        dphi = self.chart.slope_diag(a_y, b_y)
        return qz @ np.diag(dphi)

    def jacobian(self, y: Sequence[complex]) -> "np.ndarray":
        """Analytic Jacobian of G(y)=F(Phi(y))."""
        ensure_numpy()
        yy = np.asarray(y, dtype=np.complex128)
        z = self.chart.z_from_y(yy)
        Jz = self.system.jacobian(z)
        d1 = self.chart.derivative_diag(yy)
        return Jz @ np.diag(d1)

    def hessian_vec(self, y: Sequence[complex], v: Sequence[complex]) -> "np.ndarray":
        """Second directional derivative D^2G(y)[v,v]."""
        ensure_numpy()
        yy = np.asarray(y, dtype=np.complex128)
        vv = np.asarray(v, dtype=np.complex128)
        z = self.chart.z_from_y(yy)
        d1 = self.chart.derivative_diag(yy)
        d2 = self.chart.second_derivative_diag(yy)
        vz = d1 * vv
        hz = self.system.hessian_vec(z, vz)
        Jz = self.system.jacobian(z)
        return hz + Jz @ (d2 * vv * vv)

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




# ---------------------------------------------------------------------------
# 119 structured n-dimensional probe geometries
# ---------------------------------------------------------------------------

def _norm_to_sqrt_n(v: Any) -> "np.ndarray":
    """Normalize a direction to Euclidean norm sqrt(n)."""
    ensure_numpy()
    arr = np.asarray(v, dtype=np.complex128)
    n = len(arr)
    nm = float(np.linalg.norm(arr))
    if nm <= 0 or not math.isfinite(nm):
        return arr
    return arr / nm * math.sqrt(max(1, n))


def _coord_phases(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """Deterministic coordinate phases to avoid real/axis degeneracy."""
    vals = []
    for j in range(n):
        h = splitmix64(seed + 0xC0FFEE + 4099 * (j + 1) + 65537 * (ep + 1) + 131071 * (k + 1))
        vals.append(phase(2.0 * math.pi * u01(h)))
    return np.asarray(vals, dtype=np.complex128)


def axis_direction(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """Cross-polytope direction: rotated +/- coordinate axes."""
    j = k % max(1, n)
    sign = 1.0 if ((k // max(1, n)) % 2 == 0) else -1.0
    v = np.zeros(n, dtype=np.complex128)
    v[j] = sign * _coord_phases(n, seed, ep, k)[j]
    return _norm_to_sqrt_n(v)


def simplex_direction(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """Regular-simplex-like direction, projected from R^{n+1} into C^n."""
    m = n + 1
    idx = k % m
    row = np.full(m, -1.0 / m, dtype=np.float64)
    row[idx] += 1.0
    v = row[:n].astype(np.complex128) * _coord_phases(n, seed, ep, k)
    return _norm_to_sqrt_n(v)


def hypercube_direction(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """Budgeted phase-rotated hypercube corner direction."""
    if n <= 12:
        idx = k % (1 << n)
    else:
        idx = splitmix64(seed + 0xABCDEF + 104729 * (ep + 1) + 9176 * (k + 1))
    signs = [1.0 if ((idx >> j) & 1) else -1.0 for j in range(n)]
    v = np.asarray(signs, dtype=np.complex128) * _coord_phases(n, seed, ep, k)
    return _norm_to_sqrt_n(v)


def spherical_direction(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """The 118-style generic complex star direction."""
    return raw_direction(n, seed + 104729 * (ep + 1) + 7919 * (k + 1), seed ^ (0x116116 + 17 * k), True)


def hybrid_direction(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """Axes + simplex + budgeted hypercube + spherical fill."""
    if k < 2 * n:
        return axis_direction(n, seed, ep, k)
    k2 = k - 2 * n
    if k2 < n + 1:
        return simplex_direction(n, seed, ep, k2)
    k3 = k2 - (n + 1)
    max_cube = min(2 ** min(n, 8), 2 * n + 4)
    if k3 < max_cube:
        return hypercube_direction(n, seed, ep, k3)
    return spherical_direction(n, seed, ep, k3 - max_cube)


def simplex_cube_direction(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """Alternate simplex stability and hypercube diagonal exploration."""
    if k % 2 == 0:
        return simplex_direction(n, seed, ep, k // 2)
    return hypercube_direction(n, seed, ep, k // 2)


def balanced_direction(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """Interleave simplex, axes, hypercube, and spherical directions."""
    block = k // 4
    kind = k % 4
    if kind == 0:
        return simplex_direction(n, seed, ep, block)
    if kind == 1:
        return axis_direction(n, seed, ep, block)
    if kind == 2:
        return hypercube_direction(n, seed, ep, block)
    return spherical_direction(n, seed, ep, block)


def adaptive_direction(n: int, seed: int, ep: int, k: int) -> "np.ndarray":
    """Epoch-adaptive ordering: stable early, cube-heavy later."""
    block = k // 4
    kind = k % 4
    if ep < 4:
        order = (simplex_direction, spherical_direction, axis_direction, hypercube_direction)
    else:
        order = (hypercube_direction, axis_direction, simplex_direction, spherical_direction)
    return order[kind](n, seed, ep, block)


DIRECTION_GEOMETRIES = {
    "star": spherical_direction,
    "spherical": spherical_direction,
    "axes": axis_direction,
    "simplex": simplex_direction,
    "hypercube": hypercube_direction,
    "hybrid": hybrid_direction,
    "simplex_cube": simplex_cube_direction,
    "simplex-cube": simplex_cube_direction,
    "balanced": balanced_direction,
    "adaptive": adaptive_direction,
}


def probe_direction(n: int, seed: int, ep: int, k: int, geometry: str) -> "np.ndarray":
    key = str(geometry or "simplex_cube").strip().lower()
    if key not in DIRECTION_GEOMETRIES:
        raise ValueError(f"unknown probe geometry {geometry!r}; choose from {sorted(DIRECTION_GEOMETRIES)}")
    return DIRECTION_GEOMETRIES[key](n, seed, ep, k)


# ---------------------------------------------------------------------------
# Automatic geometry tuning
# ---------------------------------------------------------------------------

AUTO_GEOMETRY_DEFAULT_CANDIDATES = ("star", "adaptive", "simplex_cube", "balanced")


def normalize_geometry_name(name: str) -> str:
    key = str(name or "adaptive").strip().lower().replace("-", "_")
    if key == "spherical":
        key = "star"
    if key == "simplex_cube":
        return "simplex_cube"
    return key


def parse_geometry_candidates(raw: Optional[str]) -> list[str]:
    if raw is None or str(raw).strip() == "":
        vals = list(AUTO_GEOMETRY_DEFAULT_CANDIDATES)
    else:
        vals = [normalize_geometry_name(x) for x in str(raw).replace(";", ",").split(",") if x.strip()]
    out: list[str] = []
    for g in vals:
        if g == "auto":
            continue
        if g in DIRECTION_GEOMETRIES and g not in out:
            out.append(g)
    return out or list(AUTO_GEOMETRY_DEFAULT_CANDIDATES)


def heuristic_auto_probe_geometry(n: int, d: int, candidates: Sequence[str] = AUTO_GEOMETRY_DEFAULT_CANDIDATES) -> str:
    """Riemann-aware cheap prior for the nD Pandrosion probe geometry.

    In 122 the chart layer already does the heavy global work:

        Riemann/Mobius chart -> homothety -> return map.

    The 118-vs-122 benchmark showed that after this global normalization, the
    fastest local endpoint frame is often a simple star, not necessarily the
    older simplex-cube/adaptive prior.  This rule keeps structured frames when
    they empirically help, but lets the projective chart layer carry the global
    geometric burden.
    """
    cset = {normalize_geometry_name(c) for c in candidates}
    ordered = [normalize_geometry_name(c) for c in candidates]
    def pick(g: str) -> str:
        g = normalize_geometry_name(g)
        if g in cset:
            return g
        return next((x for x in ordered if x in DIRECTION_GEOMETRIES), "star")

    n = int(n); d = int(d)
    # Low-dimensional cases: use the benchmarked Riemann-chart synergy.
    if n <= 2:
        if d <= 2:
            return pick("star")
        if d <= 5:
            return pick("simplex_cube")
        return pick("balanced")
    if n == 3:
        if d <= 2:
            return pick("adaptive")
        return pick("balanced")
    if n == 4:
        if d <= 2:
            return pick("star")
        return pick("balanced")
    # Higher n: avoid cube-heavy frames unless degree is large enough to need
    # more diagonal exploration.
    if n <= 7:
        return pick("balanced" if d <= 3 else "adaptive")
    return pick("balanced")

def _safe_log10_residual(r: float) -> float:
    try:
        rf = float(r)
        if not math.isfinite(rf):
            return 300.0
        return math.log10(max(1e-300, max(0.0, rf)))
    except Exception:
        return 300.0


def geometry_trial_score(residual: float, accepted: bool, epochs: int, probe_evals: int, seconds: float, cond: Optional[float]) -> float:
    """Lower is better.  Residual dominates; cost and conditioning break ties."""
    score = _safe_log10_residual(residual)
    if accepted:
        score -= 2.0
    else:
        score += 0.75
    score += 0.015 * max(0, int(epochs))
    score += 0.0008 * max(0, int(probe_evals))
    score += 0.05 * max(0.0, float(seconds))
    if cond is not None and math.isfinite(float(cond)):
        score += 0.004 * math.log1p(max(0.0, float(cond)))
    else:
        score += 0.25
    return float(score)


def choose_auto_probe_geometry(
    target: TargetTrack,
    chart: ProjectiveMobiusChart,
    system: DenseKostlanSystem,
    args: argparse.Namespace,
    n: int,
    d: int,
    powers: Sequence[float],
    angles: Sequence[float],
    radii: Sequence[float],
    gains: Sequence[float],
    probe_radii: Sequence[float],
) -> tuple[str, dict[str, Any]]:
    """Select adaptive/simplex_cube/balanced automatically from (n,d).

    The selected geometry is first guessed by a cheap (n,d) prior, then, unless
    disabled, refined by a tiny calibration run.  The calibration reuses the
    same Mobius--Thales starts for all candidate geometries and scores the
    finite-slope corrector output.  This is still Pandrosion: the only tuned
    object is the endpoint-probe figure used to build Q(a,b_*).
    """
    candidates = parse_geometry_candidates(getattr(args, "auto_geometry_candidates", None))
    heuristic = heuristic_auto_probe_geometry(n, d, candidates)
    cal_trials = max(0, int(getattr(args, "auto_calibration_trials", 0)))
    cal_epochs = max(1, int(getattr(args, "auto_calibration_epochs", max(1, min(6, int(getattr(args, "epochs", 6)))))))
    cal_epochs = min(cal_epochs, max(1, int(getattr(args, "epochs", cal_epochs))))
    cal_probe_candidates = max(1, int(getattr(args, "auto_calibration_probe_candidates", max(1, int(getattr(args, "probe_candidates", 8))))))
    t0 = now()

    if cal_trials <= 0 or len(candidates) <= 1:
        return heuristic, {
            "auto_geometry_enabled": True,
            "auto_geometry_mode": "heuristic-only",
            "selected": heuristic,
            "heuristic": heuristic,
            "candidates": candidates,
            "calibration_trials": 0,
            "calibration_epochs": 0,
            "results": [],
            "seconds": 0.0,
        }

    records: list[dict[str, Any]] = []
    # Use shared projective charts/starts so the comparison is about local probe geometry,
    # not about random starts or chart choices.
    for trial in range(cal_trials):
        y_raw, cal_chart, _geom = projective_chart_homothety_start(
            n, trial, system.seed + 0x120A00, powers, angles, radii, float(args.power_cap),
            roots_found=0, duplicates=0, failures=trial, target_count=max(1, int(getattr(args, "count", 1)))
        )
        cal_target = TargetTrack(system, cal_chart)
        y0, smeta = startopt(
            cal_target, y_raw, trial, system.seed + 0x120B00,
            int(getattr(args, "startopt_steps", 1)),
            int(getattr(args, "startopt_candidates", 8)),
            gains,
            int(getattr(args, "startopt_micro_epochs", 0)),
        )
        start_res = finite_residual(cal_target, y0)
        for geom_name in candidates:
            t1 = now()
            try:
                loc = pandrosion_corrector(
                    cal_target, y0,
                    max_epochs=cal_epochs,
                    tol=float(getattr(args, "tol", 1e-12)),
                    accept=float(getattr(args, "accept", 1e-8)),
                    trial_timeout=float(getattr(args, "trial_timeout", 0.0)),
                    line_search=min(int(getattr(args, "line_search", 10)), 8),
                    probe_scale=float(getattr(args, "probe_scale", 0.035)),
                    direction_seed=system.seed + 131071 * trial,
                    probe_candidates=cal_probe_candidates,
                    probe_radii=probe_radii,
                    include_self_probe=bool(getattr(args, "probe_self", True)),
                    probe_geometry=geom_name,
                )
                z = cal_chart.z_from_y(loc["y"])
                r_orig = float(np.linalg.norm(system.eval(z)))
                accepted = bool(math.isfinite(r_orig) and r_orig < float(getattr(args, "accept", 1e-8)))
                seconds = now() - t1
                score = geometry_trial_score(
                    r_orig, accepted,
                    int(loc.get("epochs", 0)),
                    int(loc.get("probe_total_evals", 0)),
                    seconds,
                    loc.get("slope_cond"),
                )
                records.append({
                    "trial": int(trial),
                    "geometry": geom_name,
                    "start_residual": float(start_res),
                    "residual": float(r_orig),
                    "accepted": bool(accepted),
                    "status": loc.get("status"),
                    "epochs": int(loc.get("epochs", 0)),
                    "probe_total_evals": int(loc.get("probe_total_evals", 0)),
                    "seconds": float(seconds),
                    "slope_cond": loc.get("slope_cond"),
                    "score": float(score),
                })
            except Exception as exc:
                records.append({
                    "trial": int(trial),
                    "geometry": geom_name,
                    "start_residual": float(start_res),
                    "residual": float("inf"),
                    "accepted": False,
                    "status": f"calibration-error:{type(exc).__name__}",
                    "epochs": 0,
                    "probe_total_evals": 0,
                    "seconds": float(now() - t1),
                    "slope_cond": None,
                    "score": 999.0,
                })

    by_geom: dict[str, list[dict[str, Any]]] = {g: [] for g in candidates}
    for r in records:
        by_geom.setdefault(str(r["geometry"]), []).append(r)
    summaries: list[dict[str, Any]] = []
    for idx, g in enumerate(candidates):
        vals = by_geom.get(g, [])
        if not vals:
            continue
        # Tiny tie-breaker: favor the heuristic prior when scores are effectively equal.
        tie = 0.0 if g == heuristic else 1e-4 * (idx + 1)
        mean_score = sum(float(v["score"]) for v in vals) / max(1, len(vals)) + tie
        finite_res = [float(v["residual"]) for v in vals if math.isfinite(float(v["residual"]))]
        summaries.append({
            "geometry": g,
            "mean_score": float(mean_score),
            "mean_log10_residual": float(sum(_safe_log10_residual(float(v["residual"])) for v in vals) / max(1, len(vals))),
            "accepted_rate": float(sum(1.0 if bool(v["accepted"]) else 0.0 for v in vals) / max(1, len(vals))),
            "mean_residual": float(sum(finite_res) / max(1, len(finite_res))) if finite_res else float("inf"),
            "mean_epochs": float(sum(float(v["epochs"]) for v in vals) / max(1, len(vals))),
            "mean_probe_evals": float(sum(float(v["probe_total_evals"]) for v in vals) / max(1, len(vals))),
            "mean_seconds": float(sum(float(v["seconds"]) for v in vals) / max(1, len(vals))),
            "runs": int(len(vals)),
        })
    summaries.sort(key=lambda r: (float(r["mean_score"]), 0 if str(r["geometry"]) == heuristic else 1))
    selected = heuristic
    best_score = None
    heuristic_score = None
    if summaries:
        best = summaries[0]
        best_score = float(best["mean_score"])
        for row in summaries:
            if str(row["geometry"]) == heuristic:
                heuristic_score = float(row["mean_score"])
                break
        if heuristic_score is None:
            heuristic_score = best_score
        # Cautious override: the micro-calibration is intentionally tiny and noisy.
        # It may overfit two or three starts, so it may override the (n,d) prior
        # only when it wins by a clear residual/cost margin.
        override_margin = float(getattr(args, "auto_override_margin", 0.75))
        if str(best["geometry"]) != heuristic and best_score < heuristic_score - override_margin:
            selected = str(best["geometry"])
    return selected, {
        "auto_geometry_enabled": True,
        "auto_geometry_mode": "heuristic+cautious-microcalibration",
        "selected": selected,
        "heuristic": heuristic,
        "best_calibrated": str(summaries[0]["geometry"]) if summaries else heuristic,
        "best_calibrated_score": best_score,
        "heuristic_score": heuristic_score,
        "override_margin": float(getattr(args, "auto_override_margin", 0.75)),
        "candidates": candidates,
        "calibration_trials": int(cal_trials),
        "calibration_epochs": int(cal_epochs),
        "calibration_probe_candidates": int(cal_probe_candidates),
        "results": summaries,
        "raw_trials": records[: min(len(records), int(getattr(args, "auto_keep_raw", 120)))],
        "seconds": float(now() - t0),
    }


def projective_chart_homothety_start(
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
) -> tuple[Any, ProjectiveMobiusChart, dict[str, Any]]:
    """Choose a chart Phi and a coordinate start y0 for the 122 flow.

    118/121 injected the Mobius--Thales geometry mostly into the starting point.
    122 promotes it to a true chart layer.  We choose a coordinate y0 and a
    projective map Phi, run Pandrosion on G(y)=F(Phi(y)), then return z=Phi(y).

    Formula per coordinate:
        Phi_j(y) = Lambda * pole_j * (cos(theta_j)*(y_j/pole_j)+sin(theta_j))
                              / (-sin(theta_j)*(y_j/pole_j)+cos(theta_j)).

    theta=0 gives affine homothety; theta≈pi/2 gives the inversion/Riemann
    chart.  Lambda is the Thales scale.  The finite-slope matrix for G remains
    exact because Phi is coordinate-wise and its finite quotient is diagonal.
    """
    ensure_numpy()
    powers2 = sorted(set(min(max(float(x), 1e-300), float(cap)) for x in powers if float(x) > 0))
    if not powers2:
        powers2 = [1.0]
    Lp, La, Lr = len(powers2), max(1, len(angles)), max(1, len(radii))
    phi = 0.6180339887498948482
    q = (trial * phi + 0.071 * roots_found + 0.013 * duplicates) % 1.0
    power_index = (int(q * Lp) + 37 * trial + 11 * (trial // 7) + 5 * roots_found) % Lp
    base_power = powers2[power_index]
    dup_pressure = (duplicates + 1.0) / (roots_found + 1.0)
    fail_pressure = (failures + 1.0) / (trial + 1.0)
    progress = min(1.0, max(0.0, roots_found / max(1.0, float(target_count))))
    thrust_ladder = [1.0, 1.6, 2.5, 4.0, 6.5, 10.0, 16.0, 25.0, 40.0, 64.0, 100.0, 160.0, 256.0]
    thrust = thrust_ladder[(trial * 17 + roots_found * 3 + duplicates) % len(thrust_ladder)]
    amp = (thrust ** (0.18 + 0.82 * progress)) * ((1.0 + dup_pressure) ** 0.42) / ((1.0 + 0.25 * fail_pressure) ** 0.15)
    Lambda = min(float(cap), max(1e-300, float(base_power) * amp))

    theta0 = angles[(trial * 19 + roots_found * 7 + duplicates * 3) % La]
    strict_single_chart = (len(angles) == 1 and (abs(theta0) < 1e-14 or abs(abs(theta0) - math.pi / 2.0) < 1e-14))
    theta_jitter = 0.0 if strict_single_chart else math.radians(4.0) * math.sin(1.324717957244746 * (trial + 1) + 0.31 * roots_found)
    radius0 = radii[(trial * 13 + failures * 5 + roots_found * 2) % Lr]
    radius = max(1e-300, float(radius0) * math.exp(0.22 * math.sin(0.754877666 * (trial + 1) + 0.17 * duplicates)))

    # y0 lives in chart coordinates.  The physical start is z0=Phi(y0).
    dvec = raw_direction(n, trial, seed, True)
    y0 = radius * dvec
    poles = np.empty(n, dtype=np.complex128)
    theta_values = np.empty(n, dtype=np.float64)
    for j in range(n):
        hj = splitmix64(seed + 0xA11CE + 982451653 * trial + 1009 * (j + 1))
        poles[j] = phase(2.0 * math.pi * u01(hj))
        if strict_single_chart:
            theta_values[j] = theta0
        else:
            theta_values[j] = theta0 + theta_jitter * math.cos(0.5 + j) + math.radians(2.0) * math.sin((j + 1) * (trial + 1) * 0.38196601125)
    chart = ProjectiveMobiusChart(scale=complex(Lambda), poles=poles, theta=theta_values)
    z0 = chart.z_from_y(y0)
    inv_pressure = float(np.mean(np.abs(np.sin(theta_values)))) if n else 0.0
    meta = {
        "homothety": float(Lambda),
        "base_homothety": float(base_power),
        "thales_thrust": float(thrust),
        "theta_deg": float(math.degrees(theta0)),
        "theta_jitter_deg": float(math.degrees(theta_jitter)),
        "theta_mean_deg": float(sum(math.degrees(t) for t in theta_values) / max(1, len(theta_values))),
        "base_radius": float(radius),
        "dup_pressure": float(dup_pressure),
        "fail_pressure": float(fail_pressure),
        "progress": float(progress),
        "inversion_pressure": float(inv_pressure),
        "chart": "124-coordinate-projective-mobius-thales-return-map",
        "physical_start_norm": float(np.linalg.norm(z0)),
        **chart.meta(),
    }
    return y0, chart, meta

# Backward-compatible alias name, but 122 returns (y0, chart, meta).
def mobius_homothety_start(*args: Any, **kwargs: Any) -> tuple[Any, ProjectiveMobiusChart, dict[str, Any]]:
    return projective_chart_homothety_start(*args, **kwargs)

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
    probe_geometry: str = "adaptive",
) -> tuple[Any, dict[str, Any]]:
    """Choose the finite-slope endpoint b with a structured nD probe figure.

    118 used an unstructured complex star of pseudo-random directions.  119 keeps
    the same exact Pandrosion step

        Q_G(a,b_*) delta = -G(a),

    but replaces the local probe figure by a selectable n-dimensional geometry.
    The default ``adaptive`` interleaves simplex/star/axis/cube probes; ``simplex_cube`` alternates regular-simplex probes, which give a
    stable affine frame, with budgeted hypercube diagonal probes, which explore
    several coordinates at once.  The selected endpoint remains

        b_* = argmin_b ||G(b)||,

    so this is not a Newton/Jacobian fallback and not a different solver: only
    the finite-slope endpoint geometry changes.
    """
    ensure_numpy()
    n = len(y)
    ynorm = max(1.0, float(np.linalg.norm(y)))
    radii = [float(r) for r in probe_radii if float(r) >= 0]
    if not radii:
        radii = [1.0]

    geometry_key = str(probe_geometry or "adaptive").strip().lower()
    candidates: list[tuple[str, Any]] = []
    if include_self_probe:
        candidates.append(("self", y.copy()))

    # Inertial probe from the previous accepted Pandrosion displacement.
    if prev_delta is not None and np.all(np.isfinite(prev_delta)):
        pdn = max(1e-300, float(np.linalg.norm(prev_delta)))
        base = prev_delta / pdn * min(max(pdn, float(probe_scale) * ynorm), 2.5 * ynorm)
        candidates.append(("inertial", y + base))

    budget = max(1, int(probe_candidates))
    k = 0
    while len(candidates) < budget:
        rad = float(probe_scale) * ynorm * radii[k % len(radii)]
        qdir = probe_direction(n, direction_seed, ep, k, geometry_key)
        qdir = _norm_to_sqrt_n(qdir)

        # Spiral the entire figure between epochs.  Structured figures would be
        # too rigid without this complex rotation; the rotation preserves the
        # figure while moving it through projective charts.
        ph = phase(0.6180339887498948 * (ep + 1) + 2.399963229728653 * (k + 1))
        step = rad * ph * qdir

        # Same radial Thales component as 118.  It couples the local probe figure
        # to the current chart scale without changing the Pandrosion formula.
        if float(np.linalg.norm(y)) > 0:
            step = step + (0.12 * rad) * y / ynorm * phase(0.38196601125 * (k + 1))

        # Avoid collapsed coordinates numerically; same finite geometry,
        # infinitesimally regularized.
        tiny = 1e-12 * ynorm
        for j in range(n):
            if abs(step[j]) < tiny:
                step[j] += tiny * phase(0.17 + j + ep + k)

        candidates.append((f"{geometry_key}-{k}", y + step))
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
        raise RuntimeError("no-finite-probe")
    return best_b, {
        "probe_mode": f"119-structured-{geometry_key}-residual-min",
        "probe_geometry": geometry_key,
        "probe_name": best_name,
        "probe_candidates": int(min(budget, len(candidates))),
        "probe_evals": int(evals),
        "probe_residual": float(best_res),
        "probe_distance": float(best_distance),
        "probe_improvement_proxy": (float(residual / best_res) if math.isfinite(best_res) and best_res > 0 and math.isfinite(residual) else None),
        "probe_self_enabled": bool(include_self_probe),
    }

def jet_inverse_chart_step(
    target: TargetTrack,
    y: Sequence[complex],
    f: Sequence[complex],
    order: int = 2,
) -> tuple[Optional["np.ndarray"], dict[str, Any]]:
    """Residual-gated inverse-chart jet predictor.

    The higher-order chart theory says that, in coordinates where G is locally
    flattened, the anchored multiplier loses its low-order terms.  Operationally
    this routine constructs the first nontrivial inverse-chart Taylor step at
    the current point.  For order=1 it is the affine inverse chart.  For
    order=2 it adds the quadratic inverse-jet correction

        u = u1 - 1/2 J^{-1} D^2G[u1,u1],   J u1 = -G(y).

    It is only a candidate: the corrector accepts it solely through residual
    line search, so bad local jets cannot force the orbit.
    """
    ensure_numpy()
    if int(order) <= 0:
        return None, {"jet_enabled": False, "jet_order": int(order), "jet_status": "disabled"}
    t0 = now()
    try:
        yy = np.asarray(y, dtype=np.complex128)
        ff = np.asarray(f, dtype=np.complex128)
        J = target.jacobian(yy)
        cond = float(np.linalg.cond(J))
        u1 = np.linalg.solve(J, -ff)
        if not np.all(np.isfinite(u1)):
            return None, {"jet_enabled": True, "jet_order": int(order), "jet_status": "nonfinite-affine", "jet_cond": cond, "jet_seconds": float(now()-t0)}
        step = u1
        corr_norm = 0.0
        if int(order) >= 2:
            h2 = target.hessian_vec(yy, u1)
            corr = np.linalg.solve(J, h2)
            if np.all(np.isfinite(corr)):
                step = u1 - 0.5 * corr
                corr_norm = float(np.linalg.norm(0.5 * corr))
        return np.asarray(step, dtype=np.complex128), {
            "jet_enabled": True,
            "jet_order": int(min(max(1, int(order)), 2)),
            "jet_status": "candidate",
            "jet_cond": cond,
            "jet_step_norm": float(np.linalg.norm(step)),
            "jet_affine_step_norm": float(np.linalg.norm(u1)),
            "jet_quadratic_correction_norm": float(corr_norm),
            "jet_seconds": float(now() - t0),
        }
    except Exception as exc:
        return None, {
            "jet_enabled": True,
            "jet_order": int(order),
            "jet_status": f"jet-error:{type(exc).__name__}",
            "jet_seconds": float(now() - t0),
        }


def _residual_line_search(
    target: TargetTrack,
    y: "np.ndarray",
    delta: "np.ndarray",
    base_r: float,
    best_r: float,
    line_search: int,
) -> tuple[bool, "np.ndarray", float, float, int]:
    """Shared residual line search for Pandrosion and jet candidates."""
    ensure_numpy()
    evals = 0
    best_local_y = y
    best_local_r = float("inf")
    best_lam = 0.0
    for k in range(max(1, int(line_search))):
        lam = 1.0 / (2.0 ** k)
        yy = y + lam * delta
        rr = finite_residual(target, yy)
        evals += 1
        if math.isfinite(rr) and rr < best_local_r:
            best_local_y = yy
            best_local_r = float(rr)
            best_lam = float(lam)
        if math.isfinite(rr) and (rr < base_r or rr < best_r):
            return True, yy, float(rr), float(lam), evals
    return False, best_local_y, float(best_local_r), float(best_lam), evals


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
    probe_geometry: str = "adaptive",
    jet_order: int = 2,
    jet_every: int = 1,
    jet_max_step: float = 18.0,
    jet_blend: bool = True,
) -> dict[str, Any]:
    """124 jet-aware probe Pandrosion corrector.

    The exact finite-slope Pandrosion correction is still built each epoch:

        b_* = argmin ||G(b)||,
        Q_G(a,b_*) delta_P = -G(a).

    124 adds a residual-gated inverse-chart jet candidate.  The theory from the
    higher-order chart papers says that a chart with flattened jets reduces the
    anchored multiplier.  Computationally we approximate that inverse chart at
    the current point through first/second analytic directional jets of
    G(y)=F(Phi(y)).  The resulting candidate is never trusted blindly: it is
    compared against the pure Pandrosion step by the same residual line search.
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
    last_jet_meta: dict[str, Any] = {"jet_enabled": bool(int(jet_order) > 0), "jet_order": int(jet_order), "jet_status": "not-run"}
    total_probe_evals = 0
    total_line_evals = 0
    jet_attempts = 0
    jet_accepts = 0
    pandrosion_accepts = 0
    blend_accepts = 0
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

        # 1) Exact finite-slope Pandrosion direction.
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
                probe_geometry=str(probe_geometry),
            )
            last_probe_meta = pmeta
            total_probe_evals += int(pmeta.get("probe_evals", 0))
        except Exception as exc:
            status = f"probe-error:{type(exc).__name__}"
            break
        try:
            Q = target.slope_matrix(y, b)
            last_cond = float(np.linalg.cond(Q))
            delta_p = np.linalg.solve(Q, -f)
        except Exception as exc:
            status = f"slope-solve-error:{type(exc).__name__}"
            break
        if not np.all(np.isfinite(delta_p)):
            status = "nonfinite-pandrosion-step"
            break

        ynorm = max(1.0, float(np.linalg.norm(y)))
        max_step = max(1e-12, float(jet_max_step)) * ynorm
        dnorm = float(np.linalg.norm(delta_p))
        if dnorm > max_step:
            delta_p = delta_p * (max_step / max(dnorm, 1e-300))

        candidates: list[tuple[str, np.ndarray, dict[str, Any]]] = [
            ("pandrosion", np.asarray(delta_p, dtype=np.complex128), {"candidate_kind": "pandrosion"})
        ]

        # 2) Optional inverse-chart jet candidate.  This is deliberately local
        # and residual-gated.  It often shortens the linear Pandrosion regime
        # once the Riemann chart has already normalized the global scale.
        if int(jet_order) > 0 and (int(jet_every) <= 1 or ep % max(1, int(jet_every)) == 0):
            jet_attempts += 1
            delta_j, jmeta = jet_inverse_chart_step(target, y, f, order=int(jet_order))
            last_jet_meta = jmeta
            if delta_j is not None and np.all(np.isfinite(delta_j)):
                jnorm = float(np.linalg.norm(delta_j))
                if jnorm > max_step:
                    delta_j = delta_j * (max_step / max(jnorm, 1e-300))
                    jmeta = dict(jmeta); jmeta["jet_step_capped"] = True
                candidates.append(("jet", np.asarray(delta_j, dtype=np.complex128), {"candidate_kind": "jet", **jmeta}))
                if bool(jet_blend):
                    # A safe bridge between the exact finite-slope step and the
                    # local inverse-chart step; useful when either one is a bit
                    # overconfident in a curved chart.
                    candidates.append(("blend", 0.5 * (np.asarray(delta_p) + np.asarray(delta_j)), {"candidate_kind": "blend", **jmeta}))
        else:
            last_jet_meta = {"jet_enabled": bool(int(jet_order) > 0), "jet_order": int(jet_order), "jet_status": "skipped"}

        # 3) Residual-gated selection among candidate steps.
        accepted = False
        base_r = r
        best_candidate = None
        best_candidate_data = None
        # Try pure Pandrosion first to preserve the original flow when it works.
        # Then test jet/blend candidates; the accepted step is the one with the
        # smallest residual among all successful line searches.
        for cname, delta, cmeta in candidates:
            if not np.all(np.isfinite(delta)):
                continue
            cnorm = float(np.linalg.norm(delta))
            if cnorm > max_step:
                delta = delta * (max_step / max(cnorm, 1e-300))
            ok_ls, yy, rr, lam, levals = _residual_line_search(
                target=target,
                y=y,
                delta=delta,
                base_r=base_r,
                best_r=best_r,
                line_search=int(line_search),
            )
            total_line_evals += int(levals)
            if ok_ls:
                record = (rr, cname, yy, lam, delta, cmeta)
                if best_candidate is None or rr < best_candidate[0]:
                    best_candidate = record
            elif best_candidate_data is None or rr < best_candidate_data[0]:
                best_candidate_data = (rr, cname, yy, lam, delta, cmeta)

        if best_candidate is not None:
            rr, cname, yy, lam, used_delta, cmeta = best_candidate
            prev_delta = lam * used_delta
            y = yy
            if rr < best_r:
                best_y = yy.copy(); best_r = rr
            accepted = True
            if cname == "jet":
                jet_accepts += 1
            elif cname == "blend":
                blend_accepts += 1
            else:
                pandrosion_accepts += 1
            last_probe_meta = dict(last_probe_meta)
            last_probe_meta.update({
                "accepted_candidate": cname,
                "accepted_lambda": float(lam),
                "accepted_residual": float(rr),
                "candidate_count": int(len(candidates)),
            })
            if cname in ("jet", "blend"):
                last_jet_meta = dict(last_jet_meta)
                last_jet_meta.update(cmeta)
                last_jet_meta["jet_last_accepted_as"] = cname
        epochs = ep + 1
        if not accepted:
            status = "no-decrease"
            if best_candidate_data is not None:
                rr, cname, yy, lam, _used_delta, _cmeta = best_candidate_data
                last_probe_meta = dict(last_probe_meta)
                last_probe_meta.update({
                    "best_rejected_candidate": cname,
                    "best_rejected_residual": float(rr),
                    "best_rejected_lambda": float(lam),
                    "candidate_count": int(len(candidates)),
                })
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
        "corrector": "124-jet-aware-probe-pandrosion-exact-telescopic-slope-plus-residual-gated-inverse-jet",
        "probe_total_evals": int(total_probe_evals),
        "line_search_evals": int(total_line_evals),
        "jet_attempts": int(jet_attempts),
        "jet_accepts": int(jet_accepts),
        "pandrosion_accepts": int(pandrosion_accepts),
        "blend_accepts": int(blend_accepts),
        **last_probe_meta,
        **last_jet_meta,
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


def select_chart_start(
    system: DenseKostlanSystem,
    args: argparse.Namespace,
    n: int,
    trial: int,
    powers: Sequence[float],
    angles: Sequence[float],
    radii: Sequence[float],
    gains: Sequence[float],
    roots_found: int,
    duplicates: int,
    failures: int,
    target_count: int,
) -> tuple[Any, TargetTrack, ProjectiveMobiusChart, dict[str, Any], dict[str, Any]]:
    """Select a projective chart/start pair by residual before correction.

    This is the missing 122 layer: instead of injecting one Mobius point and
    hoping its chart is good, Pandrosion tests a small chart portfolio
    (affine/Riemann/inversion angles already encoded by ``angles``), normalizes
    each one by Thales scaling, performs the same StartOpt in its own chart
    coordinates, and keeps the chart with the smallest residual of G(y)=F(Phi(y)).
    No Newton/Jacobian or external solver is used.
    """
    ensure_numpy()
    budget = max(1, int(getattr(args, "chart_candidates", 3)))
    best: Optional[tuple[float, Any, TargetTrack, ProjectiveMobiusChart, dict[str, Any], dict[str, Any]]] = None
    records: list[dict[str, Any]] = []
    for cc in range(budget):
        # Large deterministic offset: same outer trial, different chart sample.
        ttrial = int(trial) + cc * 1000003
        y_raw, chart, geom = projective_chart_homothety_start(
            n, ttrial, system.seed + 0x113000 + 7919 * cc, powers, angles, radii, float(args.power_cap),
            roots_found=roots_found, duplicates=duplicates, failures=failures + cc, target_count=target_count,
        )
        target = TargetTrack(system, chart)
        y0, smeta = startopt(
            target, y_raw, ttrial, system.seed + 0x112555 + 15485863 * cc,
            int(args.startopt_steps), int(args.startopt_candidates), gains, int(args.startopt_micro_epochs),
        )
        r = finite_residual(target, y0)
        z0 = chart.z_from_y(y0)
        # residual dominates; prefer moderate physical size when residual ties.
        score = math.log1p(max(0.0, r)) + 1e-12 * math.log1p(float(np.linalg.norm(z0))) if math.isfinite(r) else float("inf")
        rec = {
            "candidate": int(cc),
            "trial_index": int(ttrial),
            "score": float(score),
            "residual": float(r),
            "physical_norm": float(np.linalg.norm(z0)),
            "theta_mean_deg": geom.get("theta_mean_deg"),
            "homothety": geom.get("homothety"),
            "chart": geom.get("chart"),
        }
        records.append(rec)
        if best is None or score < best[0]:
            best = (float(score), y0, target, chart, geom, smeta)
    if best is None:
        raise RuntimeError("no-chart-candidate")
    score, y0, target, chart, geom, smeta = best
    geom = dict(geom)
    geom["chart_selection"] = "min-start-residual-over-projective-chart-portfolio"
    geom["chart_candidates"] = int(budget)
    geom["chart_selected_score"] = float(score)
    geom["chart_candidate_records"] = records[: min(len(records), int(getattr(args, "keep_chart_candidates", 12)))]
    return y0, target, chart, geom, smeta


# ---------------------------------------------------------------------------
# Main extractor
# ---------------------------------------------------------------------------

def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    ensure_numpy()
    t_case = now()
    n, d = parse_case(case_raw)
    system = DenseKostlanSystem.make(n, d, seed_index=int(args.seed_index), equation_normalize=bool(args.equation_normalize))
    # Base affine chart is used only for cheap geometry micro-calibration when requested.
    base_chart = ProjectiveMobiusChart.affine(n, scale=float(args.linear_scale))
    base_target = TargetTrack(system, base_chart)

    powers = sorted(set(round(float(x), 16) for x in parse_float_list(args.powers, DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    chart_layer = str(getattr(args, "chart_layer", "projective")).strip().lower()
    angles = [math.radians(x) for x in parse_float_list(args.angles, DEFAULT_ANGLES_DEG)]
    if chart_layer in ("affine", "linear"):
        angles = [0.0]
    elif chart_layer in ("inversion", "reciprocal"):
        angles = [math.pi / 2.0]
    angles_deg = [math.degrees(x) for x in angles]
    radii = parse_float_list(args.rays, DEFAULT_RADII, positive=True)
    gains = parse_float_list(args.startopt_gains, DEFAULT_GAINS, positive=True)
    probe_radii = parse_float_list(args.probe_radii, [0.0, 0.35, 0.7, 1.0, 1.6, 2.6, 4.2], positive=False)
    probe_radii = [r for r in probe_radii if r >= 0] or [0.0, 1.0]

    requested_probe_geometry = normalize_geometry_name(str(getattr(args, "probe_geometry", "auto")))
    if requested_probe_geometry in ("auto", "tuned", "auto_nd"):
        selected_probe_geometry, auto_geometry_meta = choose_auto_probe_geometry(
            target=base_target, chart=base_chart, system=system, args=args, n=n, d=d,
            powers=powers, angles=angles, radii=radii, gains=gains, probe_radii=probe_radii,
        )
    else:
        selected_probe_geometry = requested_probe_geometry
        auto_geometry_meta = {
            "auto_geometry_enabled": False,
            "auto_geometry_mode": "manual",
            "selected": selected_probe_geometry,
            "heuristic": heuristic_auto_probe_geometry(n, d, AUTO_GEOMETRY_DEFAULT_CANDIDATES),
            "candidates": [],
            "calibration_trials": 0,
            "calibration_epochs": 0,
            "results": [],
            "seconds": 0.0,
        }

    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    duplicates = 0
    failures = 0
    t_extract = now()

    for trial in range(int(args.pool)):
        if len(roots) >= int(args.count):
            break
        y0, target, chart, geom, smeta = select_chart_start(
            system=system,
            args=args,
            n=n,
            trial=trial,
            powers=powers,
            angles=angles,
            radii=radii,
            gains=gains,
            roots_found=len(roots),
            duplicates=duplicates,
            failures=failures,
            target_count=int(args.count),
        )
        loc = pandrosion_corrector(
            target, y0,
            max_epochs=int(args.epochs),
            tol=float(args.tol),
            accept=float(args.accept),
            trial_timeout=float(args.trial_timeout),
            line_search=int(args.line_search),
            probe_scale=float(getattr(args, "probe_scale", 0.035)),
            direction_seed=system.seed + 7919 * trial,
            probe_candidates=int(args.probe_candidates),
            probe_radii=probe_radii,
            include_self_probe=bool(args.probe_self),
            probe_geometry=str(selected_probe_geometry),
            jet_order=int(getattr(args, "jet_order", 2)),
            jet_every=int(getattr(args, "jet_every", 1)),
            jet_max_step=float(getattr(args, "jet_max_step", 18.0)),
            jet_blend=bool(getattr(args, "jet_blend", True)),
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
            "probe_mode": loc.get("probe_mode"),
            "probe_geometry": loc.get("probe_geometry"),
            "probe_name": loc.get("probe_name"),
            "probe_candidates": loc.get("probe_candidates"),
            "probe_total_evals": loc.get("probe_total_evals"),
            "probe_residual": loc.get("probe_residual"),
            "probe_distance": loc.get("probe_distance"),
            "probe_improvement_proxy": loc.get("probe_improvement_proxy"),
            "accepted_candidate": loc.get("accepted_candidate"),
            "accepted_lambda": loc.get("accepted_lambda"),
            "line_search_evals": loc.get("line_search_evals"),
            "jet_enabled": loc.get("jet_enabled"),
            "jet_order": loc.get("jet_order"),
            "jet_status": loc.get("jet_status"),
            "jet_cond": loc.get("jet_cond"),
            "jet_attempts": loc.get("jet_attempts"),
            "jet_accepts": loc.get("jet_accepts"),
            "pandrosion_accepts": loc.get("pandrosion_accepts"),
            "blend_accepts": loc.get("blend_accepts"),
            "jet_last_accepted_as": loc.get("jet_last_accepted_as"),
            "jet_step_norm": loc.get("jet_step_norm"),
            "jet_affine_step_norm": loc.get("jet_affine_step_norm"),
            "jet_quadratic_correction_norm": loc.get("jet_quadratic_correction_norm"),
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
            "source": "124-standalone-jet-aware-auto-probe-pandrosion",
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
            "probe_geometry": loc.get("probe_geometry"),
            "probe_name": loc.get("probe_name"),
            "probe_candidates": loc.get("probe_candidates"),
            "probe_total_evals": loc.get("probe_total_evals"),
            "probe_residual": loc.get("probe_residual"),
            "probe_distance": loc.get("probe_distance"),
            "probe_improvement_proxy": loc.get("probe_improvement_proxy"),
            "accepted_candidate": loc.get("accepted_candidate"),
            "accepted_lambda": loc.get("accepted_lambda"),
            "line_search_evals": loc.get("line_search_evals"),
            "jet_enabled": loc.get("jet_enabled"),
            "jet_order": loc.get("jet_order"),
            "jet_status": loc.get("jet_status"),
            "jet_cond": loc.get("jet_cond"),
            "jet_attempts": loc.get("jet_attempts"),
            "jet_accepts": loc.get("jet_accepts"),
            "pandrosion_accepts": loc.get("pandrosion_accepts"),
            "blend_accepts": loc.get("blend_accepts"),
            "jet_last_accepted_as": loc.get("jet_last_accepted_as"),
            "jet_step_norm": loc.get("jet_step_norm"),
            "jet_affine_step_norm": loc.get("jet_affine_step_norm"),
            "jet_quadratic_correction_norm": loc.get("jet_quadratic_correction_norm"),
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
        "script": "124_pandrosion_standalone_jet_aware_auto_probe_engine.py",
        "autonomous": True,
        "dependencies": standalone_dependency_audit(),
        "mode": "standalone-full-riemann-chart-jet-aware-auto-probe-pandrosion",
        "flow_formula": "choose Phi_{Lambda,theta,pole}; run StartOpt and PURE Pandrosion on G(y)=F(Phi(y)) with structured nD probes; return z=Phi(y); Q_G = Q_F(Phi(a),Phi(b)) diag((Phi(b)-Phi(a))/(b-a))",
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
        "chart_layer": "per-trial coordinate-wise Riemann/Mobius chart with homothety and return map",
        "parameters": {
            "count": int(args.count),
            "pool": int(args.pool),
            "accept": float(args.accept),
            "tol": float(args.tol),
            "cluster_sep": float(args.cluster_sep),
            "epochs": int(args.epochs),
            "trial_timeout": float(args.trial_timeout),
            "line_search": int(args.line_search),
            "jet_order": int(getattr(args, "jet_order", 2)),
            "jet_every": int(getattr(args, "jet_every", 1)),
            "jet_max_step": float(getattr(args, "jet_max_step", 18.0)),
            "jet_blend": bool(getattr(args, "jet_blend", True)),
            "probe_scale": float(args.probe_scale),
            "probe_candidates": int(args.probe_candidates),
            "probe_geometry_requested": str(getattr(args, "probe_geometry", "auto")),
            "probe_geometry_selected": str(selected_probe_geometry),
            "probe_geometry": str(selected_probe_geometry),
            "chart_layer": chart_layer,
            "chart_candidates": int(args.chart_candidates),
            "auto_geometry": auto_geometry_meta,
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
    p = argparse.ArgumentParser(description="Standalone 123 full Riemann-chart Pandrosion engine with Riemann-aware auto probe geometry, general multivariate, no solver fallback.")
    p.add_argument("--cases", default="2,4", help="comma-separated case n,d; multiple cases can be separated by ';'")
    p.add_argument("--seed-index", "--seed", type=int, default=0, help="deterministic case seed index")
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--no-equation-normalize", dest="equation_normalize", action="store_false")
    p.add_argument("--linear-scale", type=float, default=1.0, help="base affine scale used by the chart layer")
    p.add_argument("--chart-layer", default="projective", choices=["projective", "affine", "linear", "inversion", "reciprocal"], help="global chart layer: projective mixes affine/Riemann/inversion angles; affine forces theta=0; inversion forces theta=90 degrees")
    p.add_argument("--chart-candidates", type=int, default=3, help="number of projective chart/start candidates scored per trial before the Pandrosion corrector")
    p.add_argument("--keep-chart-candidates", type=int, default=8, help="number of chart-candidate records kept per trial/root in JSON")
    p.add_argument("--count", "--thales-count", "--useful-count", type=int, default=8)
    p.add_argument("--pool", "--thales-pool", "--useful-pool", type=int, default=4096)
    p.add_argument("--epochs", "--thales-epochs", "--useful-epochs", type=int, default=24)
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--accept", "--residual-accept", type=float, default=1e-8)
    p.add_argument("--cluster-sep", type=float, default=1e-8)
    p.add_argument("--trial-timeout", "--thales-trial-timeout", "--useful-trial-timeout", type=float, default=0.0)
    p.add_argument("--line-search", type=int, default=12)
    p.add_argument("--jet-order", type=int, default=2, choices=[0, 1, 2], help="local inverse-chart jet candidate order: 0=off, 1=affine inverse chart, 2=quadratic inverse-jet flattening")
    p.add_argument("--jet-every", type=int, default=1, help="attempt the jet candidate every N epochs")
    p.add_argument("--jet-max-step", type=float, default=18.0, help="cap candidate step norm by this multiple of max(1,||y||)")
    p.add_argument("--jet-blend", action="store_true", default=True, help="also test a 50/50 blend of Pandrosion and jet steps")
    p.add_argument("--no-jet-blend", dest="jet_blend", action="store_false")
    p.add_argument("--probe-scale", type=float, default=0.035, help="base scale for finite-slope endpoint probes")
    p.add_argument("--probe-candidates", type=int, default=8, help="number of theorem-guided b-probes scored by ||F(b)|| per epoch")
    p.add_argument("--probe-geometry", default="auto", choices=["auto", "star", "spherical", "axes", "simplex", "hypercube", "hybrid", "simplex_cube", "simplex-cube", "balanced", "adaptive"], help="structured nD geometry used to generate finite-slope endpoint probes; auto chooses among star/adaptive/simplex_cube/balanced with a Riemann-aware heuristic plus optional micro-calibration")
    p.add_argument("--auto-geometry-candidates", default="star,adaptive,simplex_cube,balanced", help="candidate geometries tested when --probe-geometry=auto")
    p.add_argument("--auto-calibration-trials", type=int, default=0, help="number of shared-start calibration trials used by auto geometry selection; 0 = heuristic-only n,d tuning (default)")
    p.add_argument("--auto-calibration-epochs", type=int, default=6, help="epochs per candidate during auto geometry calibration")
    p.add_argument("--auto-calibration-probe-candidates", type=int, default=8, help="probe candidates per epoch during auto geometry calibration")
    p.add_argument("--auto-override-margin", type=float, default=0.75, help="minimum score advantage required for calibration to override the (n,d) heuristic")
    p.add_argument("--auto-keep-raw", type=int, default=120, help="number of raw auto-calibration records kept in JSON")
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
    p.add_argument("--outdir", default="/mnt/data/124_jet_aware_auto_probe_out")
    p.add_argument("--keep-trials", type=int, default=160)
    p.add_argument("--verbose-trials", action="store_true")
    p.add_argument("--quiet", action="store_true", help="write JSON but suppress the human-readable console summary")
    p.add_argument("--self-test", action="store_true", help="run a small standalone smoke test by overriding cases/count/pool/epochs")
    p.add_argument("--dependency-audit", action="store_true", help="print standalone dependency audit and exit")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(getattr(args, "dependency_audit", False)):
        print(json.dumps(standalone_dependency_audit(), indent=2), flush=True)
        return 0
    if bool(getattr(args, "self_test", False)):
        # A fast, low-cost smoke test proving that the single file runs without
        # importing any previous Pandrosion script. It does not try to exhaust all roots.
        args.cases = "2,2"
        args.count = min(int(args.count), 1)
        args.pool = min(int(args.pool), 48)
        args.epochs = min(int(args.epochs), 8)
        args.probe_candidates = min(int(args.probe_candidates), 6)
        args.chart_candidates = min(int(args.chart_candidates), 3)
        args.auto_calibration_trials = 0
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    outputs = [run_case(args, c) for c in cases]
    final: dict[str, Any]
    if len(outputs) == 1:
        final = outputs[0]
    else:
        final = {"script": "124_pandrosion_standalone_jet_aware_auto_probe_engine.py", "standalone": True, "dependencies": standalone_dependency_audit(), "cases": outputs}
    if args.out:
        out = Path(args.out)
    else:
        first = cases[0].replace(",", "x") if cases else "case"
        out = Path(args.outdir) / f"124_jet_aware_auto_probe_{first}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2), encoding="utf-8")

    if not bool(getattr(args, "quiet", False)):
        print("=" * 120, flush=True)
        print("124 standalone FULL RIEMANN-CHART + JET-AWARE AUTO-PROBE Pandrosion", flush=True)
        print("No dependency on previous Python scripts; explicit chart Phi(y), homothety/inversion/return map + Riemann-aware automatic nD b-probes + residual-gated jet inverse-chart candidates.", flush=True)
        print("=" * 120, flush=True)
        for r in outputs:
            s = r["summary"]
            print(f"case=ks({r['n']},{r['degree']}), seed={r['seed']}, terms={r['terms']}, Bezout={r['bezout']}", flush=True)
            print(f"roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} trials={s['trials_used']} duplicates={s['duplicates']} failures={s['failures']}", flush=True)
            try:
                ag = r.get("parameters", {}).get("auto_geometry", {})
                print(f"auto_geometry: requested={r.get('parameters', {}).get('probe_geometry_requested')} selected={r.get('parameters', {}).get('probe_geometry_selected')} mode={ag.get('auto_geometry_mode')} heuristic={ag.get('heuristic')} calibration_seconds={float(ag.get('seconds', 0.0)):.3f}", flush=True)
            except Exception:
                pass
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
