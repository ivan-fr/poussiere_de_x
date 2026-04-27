"""
Multivariate Pandrosion vs Lairez 2017 — F: C^n -> C^n, find ALL D = d^n
Bezout solutions on Kostlan-Smale random systems.

Structure mirrored on the univariate ``pandrosion_vs_lairez.py``:

  1. PRIMITIVES                Schmidt slope-matrix M(a, z), evaluate F, DF
  2. PANDROSION-T_2 MV         Vector Aitken on two anchored Pandrosion steps
  3. PPH-MV (homotopy)         Cyclotomic-tensor start + T_2 corrector
  4. LAIREZ MV (homotopy)      Cyclotomic-tensor start + Newton corrector
  5. BENCHMARK                 Both methods on mvKS at small (n, d)

Multivariate Pandrosion difference quotient (Schmidt slope-matrix):
   M(a, z)_{i,j} = (F_i(z_1..z_j, a_{j+1}..a_n) - F_i(z_1..z_{j-1}, a_j..a_n))
                   / (z_j - a_j)
The matrix M satisfies F(z) - F(a) = M(a, z) (z - a).  At a=z, M -> DF(z)
(the Jacobian).  Pandrosion step: z_new = a - M(a, z)^{-1} F(a).
T_2 Aitken applied component-wise on two successive Pandrosion steps.

Multivariate homotopy (Lairez 2017):
   G_i(z) = z_i^d - 1            (start system, Bezout count d^n)
   H(z, t) = gamma (1-t) G(z) + t F(z),   gamma random complex
Anchors are tensor products of d-th roots of unity in C^n.
"""
from __future__ import annotations
import math
import time
import warnings
import sys
from itertools import product
from pathlib import Path

import numpy as np

_flow_dir = str(Path(__file__).resolve().parent.parent)
if _flow_dir not in sys.path:
    sys.path.insert(0, _flow_dir)
from importlib import import_module as _imp
_m010 = _imp("010_smale_multivariate")  # provides mvKS factory + evaluate_F + jacobian_F


# ---------------------------------------------------------------------------
# 1. PRIMITIVES
# ---------------------------------------------------------------------------

def evaluate_F(F, z):
    return _m010.evaluate_F(F, np.asarray(z, dtype=complex))


def jacobian_F(F, z, h=1e-7):
    return _m010.jacobian_F(F, np.asarray(z, dtype=complex), h=h)


def schmidt_slope_matrix(F, a, z, eps=1e-15):
    """Schmidt divided-difference matrix M(a, z) such that F(z) - F(a) = M (z - a).

    Component-wise definition: replace coordinates one by one.
        M[:,j] = (F(z_1..z_j, a_{j+1}..a_n) - F(z_1..z_{j-1}, a_j..a_n)) / (z_j - a_j)
    Falls back to Jacobian columns when a coordinate coincides.
    """
    a = np.asarray(a, dtype=complex)
    z = np.asarray(z, dtype=complex)
    n = len(a)
    M = np.zeros((n, n), dtype=complex)
    cur = a.copy()
    F_cur = evaluate_F(F, cur)
    for j in range(n):
        prev = cur.copy()
        cur[j] = z[j]
        if abs(z[j] - a[j]) < eps:
            # Use derivative column at this coordinate
            J = jacobian_F(F, cur)
            M[:, j] = J[:, j]
        else:
            F_next = evaluate_F(F, cur)
            M[:, j] = (F_next - F_cur) / (z[j] - a[j])
            F_cur = F_next
    return M


def pandrosion_step(F, a, z):
    """One multivariate Pandrosion step with anchor a:  z_new = a - M(a,z)^{-1} F(a)."""
    F_a = evaluate_F(F, a)
    M = schmidt_slope_matrix(F, a, z)
    try:
        delta = np.linalg.solve(M, F_a)
        return a - delta
    except np.linalg.LinAlgError:
        return z


def t2_pandrosion_mv(F, a, z):
    """Vector T_2 Aitken on two successive Pandrosion steps, anchored at a.

    s1 = z - M(a,z)^{-1} F(z)
    s2 = s1 - M(a,s1)^{-1} F(s1)
    T_2 = z - (s1 - z)^2 / (s2 - 2 s1 + z)   (component-wise division)
    Component-wise Aitken because the polynomial Aitken Delta^2 generalises
    coordinate-by-coordinate when M acts diagonally on the path tangent.
    """
    F_z = evaluate_F(F, z)
    M_az = schmidt_slope_matrix(F, a, z)
    try:
        s1 = z - np.linalg.solve(M_az, F_z)
    except np.linalg.LinAlgError:
        return z
    F_s1 = evaluate_F(F, s1)
    M_as1 = schmidt_slope_matrix(F, a, s1)
    try:
        s2 = s1 - np.linalg.solve(M_as1, F_s1)
    except np.linalg.LinAlgError:
        return s1
    denom = s2 - 2 * s1 + z
    safe = np.abs(denom) > 1e-15
    out = np.where(safe, z - (s1 - z) ** 2 / np.where(safe, denom, 1.0), s2)
    return out


# ---------------------------------------------------------------------------
# 2. PANDROSION-T_2 MV ORBIT (single-root local solver)
# ---------------------------------------------------------------------------

def pandrosion_orbit_mv(F, z0, max_epochs=80, tol=1e-10, eta=0.95):
    """T_2 Pandrosion with anchor = z (self-anchor -> reduces to Newton; the
    cached non-self anchor pattern is harder in the matrix case so we keep
    self-anchor T_2 with Armijo fallback on the Newton direction).
    """
    z = np.asarray(z0, dtype=complex).copy()
    for _ in range(max_epochs):
        F_z = evaluate_F(F, z)
        nz = np.linalg.norm(F_z)
        if nz < tol:
            return z, True
        z_cand = t2_pandrosion_mv(F, z, z)
        if np.linalg.norm(evaluate_F(F, z_cand)) <= eta * nz:
            z = z_cand
            continue
        # Armijo Newton fallback
        try:
            J = jacobian_F(F, z)
            direction = np.linalg.solve(J, F_z)
        except np.linalg.LinAlgError:
            return z, False
        accepted = False
        for j in range(20):
            tau = 2.0 ** (-j)
            z_new = z - tau * direction
            if np.linalg.norm(evaluate_F(F, z_new)) <= (1 - 0.1 * tau) * nz:
                z = z_new
                accepted = True
                break
        if not accepted:
            return z, False
    return z, np.linalg.norm(evaluate_F(F, z)) < tol


# ---------------------------------------------------------------------------
# 3. PPH-MV: cyclotomic-tensor start + Pandrosion-T_2 corrector + homotopy
# ---------------------------------------------------------------------------

def evaluate_G(z, d):
    return z ** d - 1.0


def jacobian_G(z, d):
    return np.diag(d * z ** (d - 1))


def evaluate_H(F, z, t, gamma, d):
    return gamma * (1 - t) * evaluate_G(z, d) + t * evaluate_F(F, z)


def jacobian_H(F, z, t, gamma, d):
    return gamma * (1 - t) * jacobian_G(z, d) + t * jacobian_F(F, z)


def H_t_partial(F, z, gamma, d):
    return -gamma * evaluate_G(z, d) + evaluate_F(F, z)


def schmidt_H(F, d, gamma, t, a, z):
    """Schmidt slope-matrix of H(., t):  M_H = gamma*(1-t)*M_G + t*M_F."""
    a = np.asarray(a, dtype=complex)
    z = np.asarray(z, dtype=complex)
    n = len(a)
    # G is separable diagonal: G_i depends only on z_i, so M_G is diagonal with
    # M_G[i,i] = (z_i^d - a_i^d)/(z_i - a_i) = sum_{k=0}^{d-1} z_i^k a_i^{d-1-k}.
    M_G = np.zeros((n, n), dtype=complex)
    for i in range(n):
        if abs(z[i] - a[i]) < 1e-15:
            M_G[i, i] = d * a[i] ** (d - 1)
        else:
            M_G[i, i] = (z[i] ** d - a[i] ** d) / (z[i] - a[i])
    M_F = schmidt_slope_matrix(F, a, z)
    return gamma * (1 - t) * M_G + t * M_F


def pph_corrector_mv(F, d, gamma, z, t, max_iter=10, tol=1e-10):
    """T_2-Pandrosion corrector with CACHED NON-SELF ANCHOR (Schmidt + Aitken).

    Univariate analog: T_2-fused-PURE uses Q(z_prev, z) (secant) instead of
    P'(z) (derivative).  Multivariate analog: M_H(z_prev, z) (Schmidt slope-
    matrix) instead of J_H(z) (Jacobian).  Order ~2.41 (golden) > Newton's 2.
    Aitken applied component-wise on two anchored Schmidt steps.
    """
    z = np.asarray(z, dtype=complex).copy()
    n = len(z)
    z_prev = z + 1e-6 * np.ones(n, dtype=complex)
    for _ in range(max_iter):
        Hv = evaluate_H(F, z, t, gamma, d)
        if np.linalg.norm(Hv) < tol:
            return z, True
        # Two anchored Pandrosion steps using Schmidt slope-matrix
        M1 = schmidt_H(F, d, gamma, t, z_prev, z)
        try:
            s1 = z - np.linalg.solve(M1, Hv)
        except np.linalg.LinAlgError:
            return z, False
        H_s1 = evaluate_H(F, s1, t, gamma, d)
        M2 = schmidt_H(F, d, gamma, t, z_prev, s1)
        try:
            s2 = s1 - np.linalg.solve(M2, H_s1)
        except np.linalg.LinAlgError:
            z_prev, z = z, s1
            continue
        # Component-wise Aitken Delta^2
        denom = s2 - 2 * s1 + z
        safe = np.abs(denom) > 1e-15
        z_new = np.where(
            safe,
            z - (s1 - z) ** 2 / np.where(safe, denom, 1.0),
            s2,
        )
        z_prev = z
        z = z_new
    return z, np.linalg.norm(evaluate_H(F, z, t, gamma, d)) < tol


def pph_track_mv(F, d, z0, gamma, dt_init=0.05, dt_min=1e-6,
                  max_steps=400, tol=1e-10):
    """One homotopy path z(t) from t=0 to t=1, Euler predictor + Pandrosion corrector."""
    z = np.asarray(z0, dtype=complex).copy()
    t = 0.0
    dt = dt_init
    fails = 0
    for _ in range(max_steps):
        if t >= 1.0 - 1e-9:
            return z, True
        h = min(dt, 1.0 - t)
        try:
            J = jacobian_H(F, z, t, gamma, d)
            H_t = H_t_partial(F, z, gamma, d)
            dzdt = -np.linalg.solve(J, H_t)
        except np.linalg.LinAlgError:
            return z, False
        z_pred = z + h * dzdt
        z_corr, ok = pph_corrector_mv(F, d, gamma, z_pred, t + h, tol=tol)
        if ok:
            z = z_corr
            t = t + h
            dt = min(dt * 1.5, 0.2)
            fails = 0
        else:
            fails += 1
            if fails >= 5 or dt < dt_min:
                return z, False
            dt *= 0.5
    return z, t >= 1.0 - 1e-9


def find_all_pph_mv(F, n, d, gamma=None, tol=1e-10):
    """Find all D=d^n solutions via PPH-MV homotopy.

    Anchors: tensor product of d-th roots of unity. Each path tracked
    independently with Euler predictor + Pandrosion corrector.
    """
    if gamma is None:
        gamma = complex(0.6, 0.8)
    omega = np.exp(2j * np.pi * np.arange(d) / d)
    found = []
    for k_tuple in product(range(d), repeat=n):
        z0 = np.array([omega[k] for k in k_tuple], dtype=complex)
        z_end, ok = pph_track_mv(F, d, z0, gamma, tol=tol)
        if not ok:
            for g_alt in [complex(0.7, 0.5), complex(-0.5, 0.6)]:
                z_end, ok = pph_track_mv(F, d, z0, g_alt, tol=tol)
                if ok:
                    break
        found.append(z_end)
    return found


# ---------------------------------------------------------------------------
# 4. LAIREZ MV (cyclotomic-tensor start + Newton corrector, baseline)
# ---------------------------------------------------------------------------

def lairez_corrector_mv(F, d, gamma, z, t, max_iter=10, tol=1e-10):
    for _ in range(max_iter):
        Hv = evaluate_H(F, z, t, gamma, d)
        if np.linalg.norm(Hv) < tol:
            return z, True
        try:
            J = jacobian_H(F, z, t, gamma, d)
            delta = np.linalg.solve(J, Hv)
        except np.linalg.LinAlgError:
            return z, False
        z = z - delta
    return z, np.linalg.norm(evaluate_H(F, z, t, gamma, d)) < tol


def lairez_track_mv(F, d, z0, gamma, dt_init=0.05, dt_min=1e-6,
                     max_steps=400, tol=1e-10):
    z = np.asarray(z0, dtype=complex).copy()
    t = 0.0
    dt = dt_init
    fails = 0
    for _ in range(max_steps):
        if t >= 1.0 - 1e-9:
            return z, True
        h = min(dt, 1.0 - t)
        try:
            J = jacobian_H(F, z, t, gamma, d)
            H_t = H_t_partial(F, z, gamma, d)
            dzdt = -np.linalg.solve(J, H_t)
        except np.linalg.LinAlgError:
            return z, False
        z_pred = z + h * dzdt
        z_corr, ok = lairez_corrector_mv(F, d, gamma, z_pred, t + h, tol=tol)
        if ok:
            z = z_corr
            t = t + h
            dt = min(dt * 1.5, 0.2)
            fails = 0
        else:
            fails += 1
            if fails >= 5 or dt < dt_min:
                return z, False
            dt *= 0.5
    return z, t >= 1.0 - 1e-9


def find_all_lairez_mv(F, n, d, gamma=None, tol=1e-10):
    if gamma is None:
        gamma = complex(0.6, 0.8)
    omega = np.exp(2j * np.pi * np.arange(d) / d)
    found = []
    for k_tuple in product(range(d), repeat=n):
        z0 = np.array([omega[k] for k in k_tuple], dtype=complex)
        z_end, ok = lairez_track_mv(F, d, z0, gamma, tol=tol)
        if not ok:
            for g_alt in [complex(0.7, 0.5), complex(-0.5, 0.6)]:
                z_end, ok = lairez_track_mv(F, d, z0, g_alt, tol=tol)
                if ok:
                    break
        found.append(z_end)
    return found


# ---------------------------------------------------------------------------
# 5. BENCHMARK
# ---------------------------------------------------------------------------

def max_residual_mv(F, roots):
    if not roots:
        return float("inf")
    return max(float(np.linalg.norm(evaluate_F(F, np.asarray(r)))) for r in roots)


def count_distinct(roots, tol=1e-3):
    seen = []
    for r in roots:
        r = np.asarray(r)
        if all(np.linalg.norm(r - s) > tol for s in seen):
            seen.append(r)
    return len(seen)


def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    print("=" * 96)
    print("Multivariate Pandrosion-T_2 homotopy vs Lairez Newton homotopy")
    print("on Kostlan-Smale random systems F: C^n -> C^n")
    print("=" * 96)

    rng = np.random.default_rng(2026)

    configs = [
        (2, 2),  # D = 4
        (2, 3),  # D = 9
        (3, 2),  # D = 8
        (2, 4),  # D = 16
        (3, 3),  # D = 27
    ]
    n_systems = 3

    print(f"\n{'(n,d)':>7} {'D':>4} {'method':<32} "
          f"{'time (s)':>10} {'roots':>8} {'distinct':>10} {'max ||F(r)||':>15}")
    print("-" * 96)

    for n, d in configs:
        D = d ** n
        polys = [_m010.random_mvks_system(n, d, rng) for _ in range(n_systems)]

        for name, fn in [
            ("PPH-MV (Pandrosion homotopy)", lambda F: find_all_pph_mv(F, n, d)),
            ("Lairez MV (Newton homotopy)",   lambda F: find_all_lairez_mv(F, n, d)),
        ]:
            t_total = 0.0
            n_roots_avg = 0
            n_distinct_avg = 0
            max_res = 0.0
            for F in polys:
                t0 = time.perf_counter()
                try:
                    roots = fn(F)
                except Exception:
                    roots = []
                t_total += time.perf_counter() - t0
                n_roots_avg += len(roots)
                n_distinct_avg += count_distinct(roots)
                if roots:
                    res = max_residual_mv(F, roots)
                    if res > max_res:
                        max_res = res
            n_roots_avg /= n_systems
            n_distinct_avg /= n_systems
            t_avg = t_total / n_systems
            print(f"{str((n, d)):>7} {D:>4} {name:<32} {t_avg:>10.3f} "
                  f"{n_roots_avg:>4.1f}/{D:<3} {n_distinct_avg:>4.1f}/{D:<3}    "
                  f"{max_res:>15.2e}")
        print()

    print("=" * 96)
    print("Done.")


if __name__ == "__main__":
    main()
