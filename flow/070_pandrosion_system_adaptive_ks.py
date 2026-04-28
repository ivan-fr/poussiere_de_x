"""
FLOW 070 -- system-adaptive dD polynomial Pandrosion benchmark vs Lairez-style.

This file continues 068/069 with the idea that Pandrosion should not choose a
fixed ambient geometry (simplex/sparse/hypercube) blindly.  Instead it builds a
geometry from the system itself.

Core identity
-------------
For every local geometry used here we keep

    F(z) - F(a) = Q_G(a,z) (z-a),
    P_G(a,z)   = a - Q_G(a,z)^(-1) F(a).

The new geometry is `system`: it computes an edge slope frame and exactifies it
with a covector induced by the polynomial supports, coefficient magnitudes,
local Jacobian column energies, and the residual direction.  The rank-one
exactification guarantees the telescope identity exactly up to roundoff.

The file also contains:
  * exact Gauss-Legendre integral Pandrosion slope `integral_gl`;
  * exact convex blend `blend = theta*integral_gl + (1-theta)*system`;
  * rescue modes inherited from 068: simplex/sparse/path/support/cube;
  * gamma/phase retries and optional affine-chart recovery;
  * a Lairez-style baseline imported from 069 (gamma total-degree homotopy +
    analytic Newton corrector).  This is a local baseline, not the official
    Lairez package.

Typical commands
----------------
Smoke:
  python -S 070_pandrosion_system_adaptive_geometry.py --suite smoke --smoke

Quick benchmark:
  python -S 070_pandrosion_system_adaptive_geometry.py --suite quick --seeds 1 \
    --gamma-retries 1 --passes 1 --budget-070 15 --budget-lairez 15 \
    --csv /mnt/data/070_quick.csv --md /mnt/data/070_quick.md

High-degree KS stress:
  python -S 070_pandrosion_system_adaptive_geometry.py --family ks \
    --cases 2,8 2,10 3,4 --seeds 1 --gamma-retries 2 --passes 1 \
    --charts 1 --budget-070 30 --budget-lairez 30
"""
from __future__ import annotations

import argparse
import cmath
import csv
import importlib.util
import json
import math
import os
import random
import signal
import statistics
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

Complex = complex
Exponent = Tuple[int, ...]
Poly = Dict[Exponent, Complex]
System = List[Poly]
Vector = List[Complex]
Matrix = List[List[Complex]]

HERE = Path(__file__).resolve().parent
M069_PATH = HERE / "069_benchmark_068_vs_lairez_ks.py"


def load_069():
    spec = importlib.util.spec_from_file_location("flow069_for_070_adaptive", str(M069_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import 069 from {M069_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow069_for_070_adaptive"] = mod
    spec.loader.exec_module(mod)
    return mod


m069 = load_069()
m068 = m069.m068

MAX_Z = 1.0e8
MAX_TERM = 1.0e240

class TimeLimitExpired(Exception):
    pass


@contextmanager
def time_limit(seconds: float | None):
    if seconds is None or seconds <= 0:
        yield
        return
    old_handler = signal.getsignal(signal.SIGALRM)

    def handler(signum, frame):
        raise TimeLimitExpired()

    signal.signal(signal.SIGALRM, handler)
    signal.setitimer(signal.ITIMER_REAL, float(seconds))
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)


# -----------------------------------------------------------------------------
# Basic wrappers.
# -----------------------------------------------------------------------------

def finite(v: Complex) -> bool:
    z = complex(v)
    return math.isfinite(z.real) and math.isfinite(z.imag)


def finite_vec(z: Sequence[Complex]) -> bool:
    return all(finite(v) for v in z)


def safe_z(z: Sequence[Complex], bound: float = MAX_Z) -> bool:
    return finite_vec(z) and all(abs(v) <= bound for v in z)


def degree(poly: Poly) -> int:
    return m069.degree(poly)


def degrees(polys: System) -> List[int]:
    return m069.degrees(polys)


def bezout(polys: System) -> int:
    return m069.bezout(polys)


def system_degree(polys: System) -> int:
    return max((degree(p) for p in polys), default=1)


def F_eval(polys: System, z: Sequence[Complex]) -> Vector:
    return m069.F_eval(polys, z)


def residual_norm(polys: System, z: Sequence[Complex]) -> float:
    return m069.residual_norm(polys, z)


def residual_norm2(polys: System, z: Sequence[Complex]) -> float:
    return m069.residual_norm2(polys, z)


def jacobian_eval(polys: System, z: Sequence[Complex]) -> Matrix:
    return m069.jacobian_eval(polys, z)


def solve_linear(A: Matrix, b: Sequence[Complex], tol: float = 1e-13) -> Vector | None:
    return m069.solve_linear(A, b, tol=tol)


def scale_system(polys: System, gammas: Sequence[Complex]) -> System:
    return m069.scale_system(polys, gammas)


def homotopy_polys(target_gamma: System, start: System, t: float) -> System:
    return m069.homotopy_polys(target_gamma, start, t)


def timed_out(deadline: float | None) -> bool:
    return deadline is not None and time.time() > deadline


# -----------------------------------------------------------------------------
# Gauss-Legendre exact segment integral Q_int.
# -----------------------------------------------------------------------------

@lru_cache(maxsize=64)
def gauss_legendre01(m: int) -> Tuple[Tuple[float, ...], Tuple[float, ...]]:
    """m-point Gauss-Legendre nodes/weights on [0,1]."""
    m = max(1, int(m))
    eps = 3e-15
    x = [0.0] * m
    w = [0.0] * m
    half = (m + 1) // 2
    for i in range(half):
        z = math.cos(math.pi * (i + 0.75) / (m + 0.5))
        pp = 0.0
        while True:
            p1 = 1.0
            p2 = 0.0
            for j in range(1, m + 1):
                p3 = p2
                p2 = p1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
            pp = m * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / pp
            if abs(z - z1) <= eps:
                break
        wi = 2.0 / ((1.0 - z * z) * pp * pp)
        x[i] = 0.5 * (1.0 - z)
        x[m - 1 - i] = 0.5 * (1.0 + z)
        w[i] = 0.5 * wi
        w[m - 1 - i] = 0.5 * wi
    return tuple(x), tuple(w)


def quad_order(D: int, cap: int) -> int:
    # Jacobian entries have segment degree <= D-1; m nodes integrate degree 2m-1.
    # q=(D+1)//2 is exact.  cap may make it an intentional approximation for
    # very high degree, but benchmark defaults keep d<=24 exact.
    return max(1, min(int(cap), (max(1, int(D)) + 1) // 2))


def Q_integral_gl(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], quad_cap: int = 16) -> Matrix:
    n = len(z)
    q = quad_order(system_degree(polys), quad_cap)
    nodes, weights = gauss_legendre01(q)
    a = [complex(v) for v in anchor]
    delta = [complex(z[j]) - a[j] for j in range(n)]
    Q = [[0.0 + 0.0j for _ in range(n)] for _ in range(n)]
    for t, wt in zip(nodes, weights):
        y = [a[j] + t * delta[j] for j in range(n)]
        J = jacobian_eval(polys, y)
        if any(not finite_vec(row) for row in J):
            return [[complex("inf") for _ in range(n)] for _ in range(n)]
        for i in range(n):
            Qi = Q[i]
            Ji = J[i]
            for j in range(n):
                Qi[j] += wt * Ji[j]
    return Q


# -----------------------------------------------------------------------------
# New system-adaptive geometry.
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class SystemProfile:
    n: int
    degree: int
    coeff_rms: float
    coeff_max: float
    activity: Tuple[float, ...]
    support_count: int


def system_profile(polys: System) -> SystemProfile:
    n = len(polys)
    coeff2 = 0.0
    coeff_max = 0.0
    count = 0
    activity = [0.0] * n
    for p in polys:
        for alpha, coeff in p.items():
            c = abs(coeff)
            coeff2 += c * c
            coeff_max = max(coeff_max, c)
            count += 1
            for j, aj in enumerate(alpha):
                if aj:
                    activity[j] += c * aj
    rms = math.sqrt(coeff2 / max(1, count)) if count else 1.0
    m = max(activity) if activity else 0.0
    act = tuple((a / m if m > 0 else 1.0) for a in activity)
    return SystemProfile(n=n, degree=system_degree(polys), coeff_rms=rms or 1.0,
                         coeff_max=coeff_max or 1.0, activity=act, support_count=count)


def edge_frame(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], F_anchor: Sequence[Complex]) -> Matrix:
    n = len(z)
    B = [[0.0 + 0.0j for _ in range(n)] for _ in range(n)]
    for j in range(n):
        den = z[j] - anchor[j]
        if abs(den) < 1e-300:
            continue
        e = list(anchor)
        e[j] = z[j]
        Fe = F_eval(polys, e)
        if not finite_vec(Fe):
            continue
        for i in range(n):
            B[i][j] = (Fe[i] - F_anchor[i]) / den
    return B


def exactify(B: Matrix, Fa: Sequence[Complex], Fz: Sequence[Complex],
             anchor: Sequence[Complex], z: Sequence[Complex], w: Sequence[Complex]) -> Matrix:
    n = len(z)
    delta = [z[j] - anchor[j] for j in range(n)]
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [Fz[i] - Fa[i] - Bd[i] for i in range(n)]
    den = sum(w[j] * delta[j] for j in range(n))
    if abs(den) < 1e-28:
        w = [d.conjugate() for d in delta]
        den = sum(w[j] * delta[j] for j in range(n))
    if abs(den) < 1e-28:
        return [row[:] for row in B]
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / den
    return Q


def local_support_energy(poly: Poly, rho: Sequence[float], j: int, coeff_scale: float) -> float:
    # A coefficient/support-only proxy for |partial_j poly| on the current box.
    s = 0.0
    for alpha, coeff in poly.items():
        aj = alpha[j]
        if not aj:
            continue
        term = (abs(coeff) / coeff_scale) * aj
        for k, ak in enumerate(alpha):
            p = ak - 1 if k == j else ak
            if p:
                # Cap the local radius to prevent one huge variable from dominating everything.
                term *= min(50.0, max(1e-12, rho[k])) ** min(12, p)
        s += term
    return s


def system_covector(polys: System, anchor: Sequence[Complex], z: Sequence[Complex],
                    F_anchor: Sequence[Complex] | None = None) -> Vector:
    """Covector induced by the actual system.

    It combines four system quantities:
      1. support/degree activity per variable;
      2. local monomial support energy on the box joining a and z;
      3. Jacobian column energy at the midpoint;
      4. residual alignment J^*F(a), so the geometry sees the current equation.

    The final direction is conjugate(delta_j) times a positive scalar.  Therefore
    w·delta is positive real unless delta=0, which stabilizes exactification.
    """
    n = len(z)
    profile = system_profile(polys)
    delta = [z[j] - anchor[j] for j in range(n)]
    scale_delta = max(1e-15, max(abs(d) for d in delta))
    mid = [(anchor[j] + z[j]) * 0.5 for j in range(n)]
    rho = [1.0 + min(20.0, max(abs(anchor[j]), abs(z[j]), abs(mid[j]))) for j in range(n)]
    coeff_scale = max(profile.coeff_rms, 1e-300)

    # Jacobian metric and residual alignment at midpoint.
    J = jacobian_eval(polys, mid) if safe_z(mid) else [[0.0 + 0.0j for _ in range(n)] for _ in range(n)]
    col_energy = [0.0] * n
    align = [0.0] * n
    if any(not finite_vec(row) for row in J):
        J = [[0.0 + 0.0j for _ in range(n)] for _ in range(n)]
    if F_anchor is None or not finite_vec(F_anchor):
        F_anchor = [0.0 + 0.0j for _ in range(n)]
    Fn = max(1e-30, math.sqrt(sum(abs(v) ** 2 for v in F_anchor)))
    for j in range(n):
        col = [J[i][j] for i in range(n)]
        col_energy[j] = math.sqrt(sum(abs(c) ** 2 for c in col)) / coeff_scale
        align[j] = abs(sum(col[i].conjugate() * F_anchor[i] for i in range(n))) / (coeff_scale * Fn)

    support_energy = [0.0] * n
    for j in range(n):
        support_energy[j] = sum(local_support_energy(p, rho, j, coeff_scale) for p in polys)

    # Normalize components to avoid extreme anisotropy; keep relative structure.
    def normed(v: Sequence[float]) -> List[float]:
        m = max(v) if v else 0.0
        if m <= 0 or not math.isfinite(m):
            return [0.0 for _ in v]
        return [max(0.0, min(20.0, x / m)) for x in v]

    ce = normed(col_energy)
    ae = normed(align)
    se = normed(support_energy)
    w: Vector = []
    for j in range(n):
        # The 0.15 floor prevents dropping variables completely.  The activity
        # term is global, support/J terms are local, alignment sees the residual.
        metric = 0.15 + 0.35 * profile.activity[j] + 0.35 * se[j] + 0.70 * ce[j] + 0.35 * ae[j]
        # Prefer directions actually moved by z-a but do not overfit to one huge coordinate.
        dj = abs(delta[j]) / scale_delta
        metric *= 0.65 + 0.35 * min(5.0, dj)
        w.append(delta[j].conjugate() * metric)
    return w


def Q_system_adaptive(polys: System, anchor: Sequence[Complex], z: Sequence[Complex],
                      F_anchor: Sequence[Complex]) -> Matrix:
    B = edge_frame(polys, anchor, z, F_anchor)
    Fz = F_eval(polys, z)
    if not finite_vec(Fz):
        return B
    return exactify(B, F_anchor, Fz, anchor, z, system_covector(polys, anchor, z, F_anchor))


def blend_matrices(A: Matrix, B: Matrix, theta: float) -> Matrix:
    n = len(A)
    return [[theta * A[i][j] + (1.0 - theta) * B[i][j] for j in range(n)] for i in range(n)]


def Q_geometry(polys: System, anchor: Sequence[Complex], z: Sequence[Complex],
               F_anchor: Sequence[Complex], mode: str, quad_cap: int = 16) -> Matrix:
    mode = mode.lower()
    if mode in {"system", "adaptive", "system_adaptive"}:
        return Q_system_adaptive(polys, anchor, z, F_anchor)
    if mode in {"integral", "integral_gl", "gl"}:
        return Q_integral_gl(polys, anchor, z, quad_cap=quad_cap)
    if mode == "blend":
        Qs = Q_system_adaptive(polys, anchor, z, F_anchor)
        Qi = Q_integral_gl(polys, anchor, z, quad_cap=quad_cap)
        if any(not finite_vec(row) for row in Qi):
            return Qs
        if any(not finite_vec(row) for row in Qs):
            return Qi
        return blend_matrices(Qi, Qs, theta=0.55)
    if mode in {"simplex", "sparse", "support", "hypercube_cov", "path", "cube"}:
        return m068.Q_geometry(polys, anchor, z, F_anchor, mode)
    raise ValueError(f"unknown geometry mode {mode}")


# -----------------------------------------------------------------------------
# Pandrosion corrector.
# -----------------------------------------------------------------------------

def pandrosion_h(polys: System, anchor: Sequence[Complex], z: Sequence[Complex],
                 F_anchor: Sequence[Complex], mode: str, quad_cap: int = 16) -> Tuple[Vector, bool]:
    if not safe_z(anchor) or not safe_z(z) or not finite_vec(F_anchor):
        return list(z), False
    Q = Q_geometry(polys, anchor, z, F_anchor, mode, quad_cap=quad_cap)
    if any(not finite_vec(row) for row in Q):
        return list(z), False
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), False
    out = [anchor[i] + step[i] for i in range(len(z))]
    return (out, True) if safe_z(out) else (list(z), False)


def T2_step(polys: System, anchor: Sequence[Complex], z: Sequence[Complex],
            F_anchor: Sequence[Complex], mode: str, quad_cap: int = 16) -> Tuple[Vector, bool]:
    s1, ok = pandrosion_h(polys, anchor, z, F_anchor, mode, quad_cap=quad_cap)
    if not ok:
        return list(z), False
    s2, ok2 = pandrosion_h(polys, anchor, s1, F_anchor, mode, quad_cap=quad_cap)
    if not ok2:
        return s1, True
    out: Vector = []
    for k in range(len(z)):
        d0 = s1[k] - z[k]
        d2 = s2[k] - 2.0 * s1[k] + z[k]
        cand = s2[k] if abs(d2) < 1e-260 else z[k] - d0 * d0 / d2
        if abs(cand) > MAX_Z or abs(cand - z[k]) > 100.0 * max(1.0, abs(d0)):
            cand = s2[k]
        out.append(cand)
    return (out, True) if safe_z(out) else (s2, True)


def accept_damped(polys: System, z: Sequence[Complex], cand: Sequence[Complex],
                  best_r: float, best_z: Vector) -> Tuple[float, Vector]:
    if not safe_z(cand):
        return best_r, best_z
    direction = [cand[i] - z[i] for i in range(len(z))]
    for tau in (1.0, 0.72, 0.5, 0.35, 0.25, 0.125, 0.0625):
        trial = [z[i] + tau * direction[i] for i in range(len(z))]
        r = residual_norm(polys, trial)
        if math.isfinite(r) and r < best_r:
            best_r, best_z = r, trial
            if r < 1e-12:
                break
    return best_r, best_z


def deterministic_anchor(z: Sequence[Complex], epoch: int, strength: float = 0.070) -> Vector:
    scale = strength * max(1.0, max(abs(x) for x in z)) / (1.0 + 0.13 * epoch)
    out: Vector = []
    for j, zj in enumerate(z):
        th = 2.399963229728653 * (epoch + 1) * (j + 1)
        out.append(zj + scale * complex(math.cos(th), math.sin(1.618033988749895 * th))
                   + complex(0.004 * (j + 1), -0.003 * (j + 1)))
    return out


def parse_modes(s: str | Sequence[str]) -> Tuple[str, ...]:
    if isinstance(s, str):
        return tuple(x.strip() for x in s.split(",") if x.strip())
    return tuple(s)


def corrector(polys: System, z_init: Sequence[Complex], tol: float = 1e-9,
              max_epochs: int = 12, quad_cap: int = 16,
              modes: Sequence[str] = ("system", "integral_gl", "blend"),
              rescue_modes: Sequence[str] = (),
              deadline: float | None = None) -> Tuple[Vector, bool, int]:
    z = list(z_init)
    epochs_used = 0
    anchor = deterministic_anchor(z, 0)
    for epoch in range(max_epochs):
        if timed_out(deadline):
            return z, False, epochs_used
        epochs_used = epoch + 1
        rz = residual_norm(polys, z)
        if rz < tol:
            return z, True, epoch
        if not math.isfinite(rz):
            return z, False, epoch
        Fa = F_eval(polys, anchor)
        if finite_vec(Fa) and max(abs(v) for v in Fa) < tol:
            return list(anchor), True, epoch
        best_r, best_z = rz, list(z)

        # System geometry first: cheap and adapted to F.  Then GL integral if
        # system geometry is not enough.  Blend is exact and often stabilizes.
        tried = set()
        for mode in modes:
            if timed_out(deadline):
                break
            tried.add(mode)
            local_before = best_r
            cand, ok = pandrosion_h(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
            if ok:
                best_r, best_z = accept_damped(polys, z, cand, best_r, best_z)
            # T2 is useful only if the H step made at least weak progress.
            if best_r < 0.92 * local_before:
                cand, ok = T2_step(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                if ok:
                    best_r, best_z = accept_damped(polys, z, cand, best_r, best_z)
            if best_r < 0.22 * rz:
                break

        if best_r > 0.60 * rz:
            for mode in rescue_modes:
                if mode in tried or timed_out(deadline):
                    continue
                local_before = best_r
                cand, ok = pandrosion_h(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                if ok:
                    best_r, best_z = accept_damped(polys, z, cand, best_r, best_z)
                if best_r < 0.85 * local_before:
                    cand, ok = T2_step(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                    if ok:
                        best_r, best_z = accept_damped(polys, z, cand, best_r, best_z)
                if best_r < 0.38 * rz:
                    break

        if best_r > 0.985 * rz and not timed_out(deadline):
            cand, ok = m068.armijo_fallback(polys, z)
            if ok:
                best_r, best_z = accept_damped(polys, z, cand, best_r, best_z)

        if best_r >= rz:
            anchor = deterministic_anchor(z, epoch + 1, strength=0.050)
            continue
        z = list(best_z)
        anchor = list(z)  # next Pandrosion anchor follows the successful point
    return z, residual_norm(polys, z) < tol, epochs_used


# -----------------------------------------------------------------------------
# Homotopy tracker.
# -----------------------------------------------------------------------------

@dataclass
class PathResult:
    z: Vector
    ok: bool
    steps: int
    epochs: int
    t: float
    residual: float
    status: str


def tangent_predictor(target_gamma: System, start: System, z: Sequence[Complex], t: float,
                      dt: float, prev_z: Sequence[Complex] | None = None,
                      prev_t: float | None = None) -> Vector:
    n = len(z)
    try:
        Hcur = homotopy_polys(target_gamma, start, t)
        J = jacobian_eval(Hcur, z)
        dHdt = [fg - gs for fg, gs in zip(F_eval(target_gamma, z), F_eval(start, z))]
        dzdt = solve_linear(J, [-v for v in dHdt])
    except Exception:
        dzdt = None
    if dzdt is not None and finite_vec(dzdt):
        speed = max(abs(v) for v in dzdt)
        scale = max(1.0, max(abs(x) for x in z))
        if speed * abs(dt) < 18.0 * scale:
            pred = [z[j] + dt * dzdt[j] for j in range(n)]
            if safe_z(pred):
                return pred
    if prev_z is not None and prev_t is not None and t != prev_t:
        slope = [(z[j] - prev_z[j]) / (t - prev_t) for j in range(n)]
        pred = [z[j] + dt * slope[j] for j in range(n)]
        if safe_z(pred):
            return pred
    return list(z)


def track_one_070(target: System, start: System, target_gamma: System, z0: Sequence[Complex],
                  tol: float = 1e-9, max_steps: int = 420, max_epochs: int = 12,
                  quad_cap: int = 16, modes: Sequence[str] = ("system", "integral_gl", "blend"),
                  rescue_modes: Sequence[str] = (),
                  deadline: float | None = None) -> PathResult:
    z = list(z0)
    D = max(degrees(target))
    dt = min(0.010, max(0.0015, 0.060 / max(2, D * D)))
    dtmax = min(0.030, max(0.006, 0.18 / max(2, D)))
    t = 0.0
    steps = epochs = fails = 0
    prev_z: Vector | None = None
    prev_t: float | None = None
    status = "ok"
    while t < 1.0 - 1e-15 and steps < max_steps:
        if timed_out(deadline):
            status = "budget"
            break
        steps += 1
        tnext = min(1.0, t + dt)
        pred = tangent_predictor(target_gamma, start, z, t, tnext - t, prev_z=prev_z, prev_t=prev_t)
        Hnext = homotopy_polys(target_gamma, start, tnext)
        zc, ok, ep = corrector(Hnext, pred, tol=tol, max_epochs=max_epochs, quad_cap=quad_cap,
                                modes=modes, rescue_modes=rescue_modes, deadline=deadline)
        epochs += ep
        rh = residual_norm(Hnext, zc)
        if ok and math.isfinite(rh) and rh < 60.0 * tol:
            prev_z, prev_t = list(z), t
            z, t = list(zc), tnext
            dt = min(dtmax, dt * (1.22 if ep <= max(3, max_epochs // 3) else 1.08))
            fails = max(0, fails - 1)
            continue
        # softer rescue if very close but not inside the strict tolerance
        if math.isfinite(rh) and rh < 1e-6 and not timed_out(deadline):
            zr, okr, epr = corrector(Hnext, zc, tol=tol, max_epochs=max_epochs + 6, quad_cap=quad_cap,
                                      modes=("blend", "system", "integral_gl"),
                                      rescue_modes=rescue_modes,
                                      deadline=deadline)
            epochs += epr
            rr = residual_norm(Hnext, zr)
            if okr and rr < 60.0 * tol:
                prev_z, prev_t = list(z), t
                z, t = list(zr), tnext
                dt = min(dtmax, dt * 1.04)
                fails = max(0, fails - 1)
                continue
        fails += 1
        dt *= 0.5
        if dt < 5e-7 or fails > 80:
            status = "fail"
            break
    if t >= 1.0 - 1e-12 and not timed_out(deadline):
        zp, okp, ep = corrector(target, z, tol=tol, max_epochs=max_epochs + 10, quad_cap=quad_cap,
                                modes=("system", "integral_gl", "blend"),
                                rescue_modes=rescue_modes,
                                deadline=deadline)
        epochs += ep
        if okp and residual_norm(target, zp) < 1e-7:
            z = zp
    res = residual_norm(target, z)
    ok_final = t >= 1.0 - 1e-12 and res < 1e-7
    return PathResult(z=z, ok=ok_final, steps=steps, epochs=epochs, t=t,
                      residual=res, status=status if not ok_final else "ok")


# -----------------------------------------------------------------------------
# Starts, charts, clustering, algorithm runner.
# -----------------------------------------------------------------------------

def deterministic_gamma_vector(n: int, seed: int, retry_index: int = 0) -> List[Complex]:
    return m069.deterministic_gamma_vector(n, seed, retry_index)


def unique_append(found: List[Vector], z: Sequence[Complex], sep: float = 1e-6) -> bool:
    if not safe_z(z):
        return False
    n = len(z)
    for r in found:
        scale = max(1.0, max(abs(z[i]) for i in range(n)), max(abs(r[i]) for i in range(n)))
        if max(abs(z[i] - r[i]) for i in range(n)) <= sep * scale:
            return False
    found.append(list(z))
    return True


def polish_070(target: System, z: Sequence[Complex], tol: float, quad_cap: int,
               modes: Sequence[str], rescue_modes: Sequence[str], deadline: float | None) -> Vector | None:
    r = residual_norm(target, z)
    if math.isfinite(r) and r < 1e-7:
        return list(z)
    zp, ok, _ = corrector(target, z, tol=tol, max_epochs=20, quad_cap=quad_cap,
                          modes=modes, rescue_modes=rescue_modes,
                          deadline=deadline)
    if ok and residual_norm(target, zp) < 1e-7:
        return zp
    return None


@dataclass
class Result070:
    roots: List[Vector]
    bezout: int
    coverage: float
    paths_ok: int
    paths_fail: int
    paths: int
    steps: int
    epochs: int
    charts_used: int
    gamma_retries_used: int
    max_residual: float
    status: str


def run_chart_070(target: System, seed: int, passes: int, gamma_retries: int, tol: float,
                  max_steps: int, max_epochs: int, quad_cap: int, modes: Sequence[str],
                  rescue_modes: Sequence[str], found: List[Vector], stop_at: int,
                  deadline: float | None, path_timeout: float = 0.0) -> Tuple[int, int, int, int, int, int, str]:
    n = len(target)
    degs = degrees(target)
    B = math.prod(degs)
    okp = failp = paths = steps = epochs = gamma_used = 0
    status = "ok"
    start_cache: Dict[Tuple[int, int], Tuple[System, List[Vector]]] = {}
    for gi in range(max(1, gamma_retries)):
        if len(found) >= stop_at:
            break
        if timed_out(deadline):
            status = "budget"
            break
        gamma_used = gi + 1
        gammas = deterministic_gamma_vector(n, seed + 991 * B, gi)
        target_gamma = scale_system(target, gammas)
        for pi in range(max(1, passes)):
            if len(found) >= stop_at:
                break
            if timed_out(deadline):
                status = "budget"
                break
            key = (gi, pi)
            if key not in start_cache:
                phases = m068.deterministic_phases(n, pi, seed + 17 * B + 10007 * gi)
                start = m068.phase_start_system(degs, n, phases)
                roots0 = m068.phase_start_roots(degs, phases)
                start_cache[key] = (start, roots0)
            else:
                start, roots0 = start_cache[key]
            for z0 in roots0:
                if len(found) >= stop_at:
                    break
                if timed_out(deadline):
                    status = "budget"
                    break
                paths += 1
                # Hard per-path guard.  Some high-degree paths can enter a very
                # slow numerical corner before the soft deadline is checked.
                if path_timeout and path_timeout > 0:
                    ptimeout = path_timeout
                elif deadline is not None:
                    remaining = max(0.05, deadline - time.time())
                    ptimeout = min(remaining, max(0.15, remaining / max(1, stop_at - len(found))))
                else:
                    ptimeout = 0.0
                try:
                    with time_limit(ptimeout):
                        tr = track_one_070(target, start, target_gamma, z0, tol=tol,
                                           max_steps=max_steps, max_epochs=max_epochs,
                                           quad_cap=quad_cap, modes=modes,
                                           rescue_modes=rescue_modes, deadline=deadline)
                except TimeLimitExpired:
                    failp += 1
                    if timed_out(deadline):
                        status = "budget"
                        break
                    continue
                steps += tr.steps
                epochs += tr.epochs
                if tr.ok:
                    zp = polish_070(target, tr.z, tol, quad_cap, modes, rescue_modes, deadline)
                    if zp is not None:
                        unique_append(found, zp)
                    okp += 1
                else:
                    zp = polish_070(target, tr.z, tol, quad_cap, modes, rescue_modes, deadline) if tr.t > 0.88 else None
                    if zp is not None and unique_append(found, zp):
                        okp += 1
                    else:
                        failp += 1
                if tr.status == "budget":
                    status = "budget"
                    break
            if status == "budget":
                break
        if status == "budget":
            break
    return okp, failp, paths, steps, epochs, gamma_used, status


def run_070(target: System, seed: int, passes: int = 1, gamma_retries: int = 1,
            charts: int = 0, chart_passes: int = 1, chart_gamma_retries: int | None = None,
            tol: float = 1e-9, max_steps: int = 420, max_epochs: int = 12,
            quad_cap: int = 16, modes: Sequence[str] = ("system", "integral_gl", "blend"),
            rescue_modes: Sequence[str] = (),
            time_budget: float = 0.0, path_timeout: float = 0.0) -> Tuple[Result070, float]:
    t0 = time.time()
    deadline = t0 + time_budget if time_budget and time_budget > 0 else None
    B = bezout(target)
    found: List[Vector] = []
    okp = failp = paths = steps = epochs = 0
    gamma_used_total = 0
    charts_used = 1
    status = "ok"
    if chart_gamma_retries is None:
        chart_gamma_retries = max(1, gamma_retries)

    a, b, c, dsteps, e, gu, st = run_chart_070(target, seed, passes, gamma_retries, tol,
                                                max_steps, max_epochs, quad_cap, modes,
                                                rescue_modes, found, B, deadline, path_timeout=path_timeout)
    okp += a; failp += b; paths += c; steps += dsteps; epochs += e; gamma_used_total += gu
    if st == "budget":
        status = "budget"

    for ci in range(1, max(0, charts) + 1):
        if len(found) >= B or status == "budget":
            break
        if timed_out(deadline):
            status = "budget"
            break
        charts_used = ci + 1
        A, bvec = m068.deterministic_affine_chart(len(target), ci, 80000 + 101 * (len(target) + B))
        transformed = m068.compose_affine_system(target, A, bvec)
        chart_found: List[Vector] = []
        a, b, c, dsteps, e, gu, st = run_chart_070(transformed, seed + 1009 * ci,
                                                    chart_passes, chart_gamma_retries,
                                                    tol, max_steps, max_epochs,
                                                    quad_cap, modes, rescue_modes,
                                                    chart_found, B, deadline, path_timeout=path_timeout)
        okp += a; failp += b; paths += c; steps += dsteps; epochs += e; gamma_used_total += gu
        if st == "budget":
            status = "budget"
        for y in chart_found:
            z = m068.affine_map(A, bvec, y)
            zp = polish_070(target, z, tol, quad_cap, modes, rescue_modes, deadline)
            if zp is not None:
                unique_append(found, zp)
            if len(found) >= B or timed_out(deadline):
                break

    max_res = max((residual_norm(target, z) for z in found), default=float("inf"))
    result = Result070(roots=found, bezout=B, coverage=len(found) / max(1, B),
                       paths_ok=okp, paths_fail=failp, paths=paths, steps=steps,
                       epochs=epochs, charts_used=charts_used,
                       gamma_retries_used=gamma_used_total, max_residual=max_res,
                       status=status)
    return result, time.time() - t0


# -----------------------------------------------------------------------------
# Benchmark orchestration.
# -----------------------------------------------------------------------------

@dataclass
class BenchRow:
    family: str
    n: int
    d: int
    seed_index: int
    seed: int
    terms: int
    bezout: int
    alg: str
    roots: int
    coverage: float
    paths_ok: int
    paths_fail: int
    paths: int
    charts: int
    gamma_retries: int
    steps: int
    work: int
    max_residual: float
    seconds: float
    status: str


def parse_case(s: str) -> Tuple[int, int]:
    return m069.parse_case(s)


def suite_cases(name: str) -> List[Tuple[str, Tuple[int, int]]]:
    name = name.lower()
    if name == "smoke":
        return [("dense064", (2, 3)), ("ks", (2, 3)), ("ks", (3, 2))]
    if name == "quick":
        return [("ks", (2, 3)), ("ks", (2, 5)), ("ks", (3, 2)), ("dense064", (4, 2))]
    if name in {"high", "high_ks"}:
        return [("ks", (2, 8)), ("ks", (2, 10)), ("ks", (3, 4)), ("sparse_ks", (2, 10))]
    if name in {"full", "complete"}:
        return [
            ("dense064", (2, 3)), ("dense064", (3, 2)), ("dense064", (4, 2)),
            ("ks", (2, 3)), ("ks", (2, 5)), ("ks", (2, 8)),
            ("ks", (3, 2)), ("ks", (3, 3)), ("ks", (3, 4)),
            ("sparse_ks", (2, 8)), ("sparse_ks", (2, 10)), ("sparse_ks", (3, 4)),
        ]
    raise ValueError(f"unknown suite {name}")


def seed_for(family: str, n: int, d: int, seed_index: int) -> int:
    return m069.seed_for(family, n, d, seed_index)


def gen_system(family: str, n: int, d: int, seed: int) -> System:
    return m069.gen_system(family, n, d, seed)


def term_count(polys: System) -> int:
    return m069.term_count(polys)


def print_row(row: BenchRow) -> None:
    print(
        f"{row.family:>9} ({row.n:>2},{row.d:<2}) seed{row.seed_index:<2} "
        f"{row.alg:>12} | Bez={row.bezout:<4} terms={row.terms:<4} "
        f"roots={row.roots:<4} cov={100*row.coverage:5.1f}% "
        f"paths={row.paths_ok:>4}/{row.paths_fail:<4} tot={row.paths:<4} "
        f"charts={row.charts:<2} gamma={row.gamma_retries:<2} "
        f"steps={row.steps:<6} work={row.work:<6} "
        f"res={row.max_residual:.1e} time={row.seconds:7.2f}s {row.status}",
        flush=True,
    )


def row_from_lairez(family: str, n: int, d: int, si: int, seed: int, terms: int,
                    B: int, res: dict) -> BenchRow:
    return BenchRow(family=family, n=n, d=d, seed_index=si, seed=seed, terms=terms,
                    bezout=B, alg="lairez", roots=len(res.get("roots", [])),
                    coverage=float(res.get("coverage", 0.0)),
                    paths_ok=int(res.get("paths_ok", 0)), paths_fail=int(res.get("paths_fail", 0)),
                    paths=int(res.get("paths", 0)), charts=int(res.get("charts", 0)),
                    gamma_retries=0, steps=int(res.get("steps", 0)),
                    work=int(res.get("epochs_or_newton", 0)),
                    max_residual=float(res.get("max_residual", float("inf"))),
                    seconds=float(res.get("seconds", 0.0)), status=str(res.get("status", "ok")))


def summarize(rows: List[BenchRow]) -> str:
    lines: List[str] = []
    lines.append("\n" + "=" * 132)
    lines.append("SUMMARY BY FAMILY / ALGORITHM")
    lines.append("=" * 132)
    groups: Dict[Tuple[str, str], List[BenchRow]] = {}
    for r in rows:
        groups.setdefault((r.family, r.alg), []).append(r)
    lines.append(f"{'family':>10} {'alg':>12} {'cases':>5} {'avg cov':>8} {'full%':>7} {'avg time':>10} {'avg paths/Bez':>14}")
    lines.append("-" * 132)
    for (fam, alg), rs in sorted(groups.items()):
        avg_cov = statistics.mean(r.coverage for r in rs)
        full = sum(1 for r in rs if r.coverage >= 0.999999) / len(rs)
        avg_time = statistics.mean(r.seconds for r in rs)
        avg_pbez = statistics.mean(r.paths / max(1, r.bezout) for r in rs)
        lines.append(f"{fam:>10} {alg:>12} {len(rs):>5} {100*avg_cov:7.2f}% {100*full:6.1f}% {avg_time:9.2f}s {avg_pbez:13.2f}")
    lines.append("\n" + "=" * 132)
    lines.append("HEAD-TO-HEAD: 070-system minus Lairez-style")
    lines.append("=" * 132)
    by_key: Dict[Tuple[str, int, int, int], Dict[str, BenchRow]] = {}
    for r in rows:
        by_key.setdefault((r.family, r.n, r.d, r.seed_index), {})[r.alg] = r
    wins_cov = ties_cov = wins_time = compared = 0
    for key, dct in sorted(by_key.items()):
        if "070-system" not in dct or "lairez" not in dct:
            continue
        compared += 1
        a = dct["070-system"]
        l = dct["lairez"]
        dcov = a.coverage - l.coverage
        ratio = a.seconds / l.seconds if l.seconds > 0 else float("inf")
        if dcov > 1e-12:
            wins_cov += 1
        elif abs(dcov) <= 1e-12:
            ties_cov += 1
        if a.seconds < l.seconds:
            wins_time += 1
        fam, nn, dd, si = key
        lines.append(f"{fam:>9} ({nn},{dd}) seed{si:<2}: Δcov={100*dcov:+6.1f}%  time_ratio={ratio:7.2f}  paths_ratio={a.paths/max(1,l.paths):5.2f}")
    if compared:
        lines.append("-" * 132)
        lines.append(f"coverage: 070 wins {wins_cov}/{compared}, ties {ties_cov}/{compared}; time: 070 faster {wins_time}/{compared}")
    text = "\n".join(lines)
    print(text, flush=True)
    return text


def write_csv(rows: List[BenchRow], path: str) -> None:
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


def write_json(rows: List[BenchRow], path: str) -> None:
    with open(path, "w") as f:
        json.dump([asdict(r) for r in rows], f, indent=2)


def write_md(rows: List[BenchRow], summary: str, path: str) -> None:
    with open(path, "w") as f:
        f.write("# 070 benchmark: system-adaptive Pandrosion vs Lairez-style\n\n")
        f.write("`070-system` uses a geometry induced by the input system: supports, coefficients, local Jacobian metric and residual alignment. No Newton-ELS corrector is used.\n\n")
        f.write("| family | n | d | seed | alg | Bez | roots | cov | paths | steps | work | time | status |\n")
        f.write("|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for r in rows:
            f.write(f"| {r.family} | {r.n} | {r.d} | {r.seed_index} | {r.alg} | {r.bezout} | {r.roots} | {100*r.coverage:.1f}% | {r.paths} | {r.steps} | {r.work} | {r.seconds:.2f}s | {r.status} |\n")
        f.write("\n```text\n")
        f.write(summary)
        f.write("\n```\n")


def telescope_error(polys: System, mode: str, quad_cap: int = 16, seed: int = 1234) -> float:
    rng = random.Random(seed)
    n = len(polys)
    a = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) / 3 for _ in range(n)]
    z = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) / 3 for _ in range(n)]
    Fa = F_eval(polys, a)
    Fz = F_eval(polys, z)
    Q = Q_geometry(polys, a, z, Fa, mode, quad_cap=quad_cap)
    dlt = [z[j] - a[j] for j in range(n)]
    Qd = [sum(Q[i][j] * dlt[j] for j in range(n)) for i in range(n)]
    return max(abs(Qd[i] - (Fz[i] - Fa[i])) for i in range(n))


def main() -> None:
    parser = argparse.ArgumentParser(description="flow/070 system-adaptive Pandrosion vs Lairez-style")
    parser.add_argument("--suite", choices=["smoke", "quick", "full", "complete", "high", "high_ks"], default="quick")
    parser.add_argument("--family", default=None, help="Override suite with one family: dense064, ks, sparse_ks")
    parser.add_argument("--cases", nargs="*", type=parse_case, default=None)
    parser.add_argument("--seeds", type=int, default=1)
    parser.add_argument("--only", choices=["both", "070", "lairez"], default="both")
    parser.add_argument("--passes", type=int, default=1)
    parser.add_argument("--gamma-retries", type=int, default=1)
    parser.add_argument("--charts", type=int, default=0)
    parser.add_argument("--chart-passes", type=int, default=1)
    parser.add_argument("--chart-gamma-retries", type=int, default=None)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--max-steps", type=int, default=420)
    parser.add_argument("--max-epochs", type=int, default=12)
    parser.add_argument("--quad-cap", type=int, default=16)
    parser.add_argument("--modes", default="system,integral_gl,blend")
    parser.add_argument("--rescue-modes", default="", help="Optional fixed geometries, e.g. simplex,sparse,path. Default is none: system-built geometry only.")
    parser.add_argument("--lairez-max-steps", type=int, default=420)
    parser.add_argument("--lairez-newton-iters", type=int, default=12)
    parser.add_argument("--lairez-retries", type=int, default=2)
    parser.add_argument("--budget-070", type=float, default=0.0)
    parser.add_argument("--budget-lairez", type=float, default=0.0)
    parser.add_argument("--path-timeout", type=float, default=0.0, help="Hard timeout per 070 path; 0 derives it from --budget-070.")
    parser.add_argument("--max-bezout", type=int, default=160)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--csv", default=None)
    parser.add_argument("--json", default=None)
    parser.add_argument("--md", default=None)
    args = parser.parse_args()

    modes = parse_modes(args.modes)
    rescue_modes = parse_modes(args.rescue_modes)
    if args.family:
        if not args.cases:
            raise SystemExit("--family requires --cases")
        tasks = [(args.family, c) for c in args.cases]
    else:
        tasks = suite_cases(args.suite)

    print("=" * 132)
    print("070 -- system-adaptive dD polynomial Pandrosion vs Lairez-style gamma continuation")
    print("=" * 132)
    print(f"suite={args.suite}, seeds={args.seeds}, passes={args.passes}, gamma_retries={args.gamma_retries}, charts={args.charts}, tol={args.tol:g}")
    print(f"070 modes={','.join(modes)}; rescue={','.join(rescue_modes)}; no Newton-ELS corrector")
    print("Lairez baseline = gamma total-degree homotopy + analytic Newton corrector from 069.")
    print("=" * 132, flush=True)

    if args.smoke:
        for fam, (n, d) in [("dense064", (2, 3)), ("ks", (2, 5)), ("ks", (3, 3))]:
            polys = gen_system(fam, n, d, seed_for(fam, n, d, 0))
            es = telescope_error(polys, "system", quad_cap=args.quad_cap)
            ei = telescope_error(polys, "integral_gl", quad_cap=args.quad_cap)
            eb = telescope_error(polys, "blend", quad_cap=args.quad_cap)
            print(f"smoke {fam:>8} ({n},{d}) telescope: system={es:.2e} integral_gl={ei:.2e} blend={eb:.2e}", flush=True)

    rows: List[BenchRow] = []
    for family, (n, d) in tasks:
        B0 = d ** n
        if args.max_bezout and B0 > args.max_bezout:
            print(f"SKIP {family:>9} ({n},{d}) Bez={B0} > --max-bezout={args.max_bezout}", flush=True)
            continue
        for si in range(args.seeds):
            seed = seed_for(family, n, d, si)
            target = gen_system(family, n, d, seed)
            terms = term_count(target)
            B = bezout(target)
            if args.only in {"both", "070"}:
                t0 = time.time()
                try:
                    res070, sec070 = run_070(target, seed=70000 + seed, passes=args.passes,
                                             gamma_retries=args.gamma_retries, charts=args.charts,
                                             chart_passes=args.chart_passes,
                                             chart_gamma_retries=args.chart_gamma_retries,
                                             tol=args.tol, max_steps=args.max_steps,
                                             max_epochs=args.max_epochs, quad_cap=args.quad_cap,
                                             modes=modes, rescue_modes=rescue_modes,
                                             time_budget=args.budget_070, path_timeout=args.path_timeout)
                    row = BenchRow(family, n, d, si, seed, terms, B, "070-system",
                                   len(res070.roots), res070.coverage, res070.paths_ok,
                                   res070.paths_fail, res070.paths, res070.charts_used,
                                   res070.gamma_retries_used, res070.steps, res070.epochs,
                                   res070.max_residual, sec070, res070.status)
                except Exception as exc:
                    row = BenchRow(family, n, d, si, seed, terms, B, "070-system", 0, 0.0,
                                   0, B, B, args.charts + 1, args.gamma_retries,
                                   0, 0, float("inf"), time.time() - t0, f"error:{exc}")
                    print(f"ERROR 070 {family} ({n},{d}) seed{si}: {exc}", flush=True)
                rows.append(row)
                print_row(row)

            if args.only in {"both", "lairez"}:
                t0 = time.time()
                try:
                    resl = m069.run_lairez(target, seed=90000 + seed, tol=args.tol,
                                           max_steps=args.lairez_max_steps,
                                           max_newton_iter=args.lairez_newton_iters,
                                           retries=args.lairez_retries,
                                           time_budget=args.budget_lairez)
                    row = row_from_lairez(family, n, d, si, seed, terms, B, resl)
                    row.gamma_retries = args.lairez_retries
                except Exception as exc:
                    row = BenchRow(family, n, d, si, seed, terms, B, "lairez", 0, 0.0,
                                   0, B, B, 0, args.lairez_retries,
                                   0, 0, float("inf"), time.time() - t0, f"error:{exc}")
                    print(f"ERROR lairez {family} ({n},{d}) seed{si}: {exc}", flush=True)
                rows.append(row)
                print_row(row)

    if rows:
        summary = summarize(rows)
        if args.csv:
            write_csv(rows, args.csv)
            print(f"CSV written to {args.csv}", flush=True)
        if args.json:
            write_json(rows, args.json)
            print(f"JSON written to {args.json}", flush=True)
        if args.md:
            write_md(rows, summary, args.md)
            print(f"Markdown written to {args.md}", flush=True)
    print("=" * 132, flush=True)


if __name__ == "__main__":
    main()
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)
