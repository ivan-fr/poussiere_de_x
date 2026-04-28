"""
PAPER: 062
TITLE: T_2 with K=1 (paper 8 empirical optimum) + universal hypercube/sparse/simplex
       + Armijo fallback + soft repulsion for full Bezout coverage attempt
STATUS: empirical, distinct from flow/061 by accelerator and reanchor period

WHAT IS T_2 AT K=1 (paper 8 §3, Proposition 3.1)
================================================
Hierarchy T_n (paper 8 lines 90-96):
    s_1 = h(z),  s_2 = h(s_1)
    T_2(z) = z - (s_1 - z)^2 / (s_2 - 2 s_1 + z)         # Aitken delta-squared
    T_n(z) = T_{n-1}(z) - (s_n - T_{n-1}(z)) / (lambda_hat - 1)   # Richardson, n >= 3

K is the anchor reanchor period (paper 1 Def 4.1):
    K = infinity : fixed-anchor Pandrosion
    K = 3        : reanchor every 3 steps (used by flow/057, 060, 061)
    K = 1        : reanchor every step                                <-- here

Paper 8 line 120-122 (key clarification):
    "T_2 at K=1 is built entirely from h and Aitken's Delta^2.
     The denominator Q(z_0, .) is evaluated at distinct points
     (z_0, z) and (z_0, s_1), never as a derivative limit, so the
     pole-avoidance mechanism is preserved. In particular, T_2 at
     K=1 does NOT reduce to Newton's method."

Concretely: with K=1, the anchor is the PREVIOUS iterate (lag-1 secant).
At iteration n: anchor a_n = z_{n-1}, iterate z_n. Q(a_n, z_n) is the
divided difference between two consecutive iterates - secant-like, NOT
derivative.

Paper 8 Proposition 3.1 empirical scaling (univariate, single-root):
    T_2 mean cost / (d log d) = 0.336    (best)
    T_3 mean cost / (d log d) = 0.413    (used in 057, 060, 061)
    T_4, T_8, T_12 progressively worse.

DIFFERENCE FROM flow/061
========================
    flow/061 uses T_3 + K=3, full reanchor every 3 steps.
    flow/062 uses T_2 + K=1, lag-1 secant anchor.

The accelerator and reanchor structure are different. The four universal
Q_F geometries (hypercube/path/cube/simplex/sparse) are kept identical,
and so is the Armijo fallback (paper 7) and the soft repulsion deflation.

EXPECTATION
===========
T_2 K=1 is faster per orbit (-20% cost) but the Bezout coverage problem
is about basin geometry, not orbit speed. So I expect coverage similar
to 061 unless the lag-1 secant dynamic happens to reach different basins.

This script is the empirical test.
"""

from __future__ import annotations

import cmath
import itertools
import math
import random
import time
from itertools import product as iprod


# ============================================================================
# Polynomial primitives
# ============================================================================
def eval_poly(poly, z):
    total = 0.0 + 0.0j
    try:
        for exp, coeff in poly.items():
            term = complex(coeff)
            for zi, ai in zip(z, exp):
                term *= zi ** ai
            total += term
    except (OverflowError, ZeroDivisionError):
        return complex("inf")
    if not (math.isfinite(total.real) and math.isfinite(total.imag)):
        return complex("inf")
    return total


def F_eval(polys, z):
    return [eval_poly(p, z) for p in polys]


def is_finite_vec(z):
    return all(math.isfinite(complex(v).real) and math.isfinite(complex(v).imag) for v in z)


def residual_norm(polys, z):
    if not is_finite_vec(z):
        return float("inf")
    vals = F_eval(polys, z)
    if any(not math.isfinite(v.real) or not math.isfinite(v.imag) for v in vals):
        return float("inf")
    return max(abs(v) for v in vals)


def degree(poly):
    return max((sum(e) for e in poly), default=0)


def bezout_estimate(polys):
    out = 1
    for p in polys:
        out *= max(degree(p), 1)
    return out


# ============================================================================
# Linear algebra
# ============================================================================
def solve_linear(A, b):
    n = len(A)
    M = [list(row) + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-14:
            return None
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0.0 + 0.0j] * n
    for i in range(n - 1, -1, -1):
        rhs = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = rhs / M[i][i]
    return x


# ============================================================================
# Universal Q_F geometries (same four as flow/061)
# ============================================================================
def schmidt_path(polys, anchor, z, order):
    n = len(z)
    Q = [[0.0 + 0.0j] * n for _ in range(n)]
    cur = list(z)
    F_cur = F_eval(polys, cur)
    for j in order:
        nxt = list(cur)
        nxt[j] = anchor[j]
        F_next = F_eval(polys, nxt)
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300:
            continue
        for i in range(n):
            Q[i][j] = (F_cur[i] - F_next[i]) / denom
        cur, F_cur = nxt, F_next
    return Q


def schmidt_cube(polys, anchor, z, max_orders=8):
    n = len(z)
    perms = list(itertools.permutations(range(n)))[:max_orders]
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for order in perms:
        Q = schmidt_path(polys, anchor, z, order)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Q[i][j] / len(perms)
    return Q_avg


def edge_frame(polys, anchor, z, F_anchor):
    n = len(z)
    B = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        edge = list(anchor)
        edge[j] = z[j]
        F_edge = F_eval(polys, edge)
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300:
            continue
        for i in range(n):
            B[i][j] = (F_edge[i] - F_anchor[i]) / denom
    return B


def schmidt_simplex(polys, anchor, z, F_anchor):
    n = len(z)
    delta = [z[i] - anchor[i] for i in range(n)]
    B = edge_frame(polys, anchor, z, F_anchor)
    norm2 = sum(abs(v)**2 for v in delta)
    if norm2 < 1e-28:
        return B
    F_z = F_eval(polys, z)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [v.conjugate() for v in delta]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def schmidt_sparse(polys, anchor, z, F_anchor):
    n = len(z)
    activity = [0.0] * n
    for poly in polys:
        for exp, coeff in poly.items():
            mag = abs(coeff)
            for j, e in enumerate(exp):
                activity[j] += mag * e
    m = max(activity) if activity else 0.0
    activity = [a if a > 0 else max(1.0, 0.1 * m) for a in activity]
    delta = [z[i] - anchor[i] for i in range(n)]
    B = edge_frame(polys, anchor, z, F_anchor)
    F_z = F_eval(polys, z)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [activity[j] * delta[j].conjugate() for j in range(n)]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def Q_geometry(polys, anchor, z, F_anchor, mode):
    if mode == "path":
        return schmidt_path(polys, anchor, z, tuple(range(len(z))))
    if mode == "cube":
        return schmidt_cube(polys, anchor, z)
    if mode == "simplex":
        return schmidt_simplex(polys, anchor, z, F_anchor)
    if mode == "sparse":
        return schmidt_sparse(polys, anchor, z, F_anchor)
    raise ValueError(mode)


# ============================================================================
# Pandrosion T_2 step (paper 8 §3): Aitken Delta-squared, no Richardson
# ============================================================================
def pandrosion_h(polys, anchor, z, F_anchor, mode):
    """h(z) = anchor - Q^{-1} F(anchor) where Q = Q_F(anchor, z) for the chosen mode."""
    Q = Q_geometry(polys, anchor, z, F_anchor, mode)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if not is_finite_vec(out):
        return list(z), False
    return out, True


def T2_step(polys, anchor, z, F_anchor, mode):
    """T_2: Aitken Delta-squared on (z, s_1, s_2). Coordinate-wise.

    s_1 = h(z),  s_2 = h(s_1),
    T_2(z)_i = z_i - (s_1_i - z_i)^2 / (s_2_i - 2 s_1_i + z_i).

    This is paper 8 line 93 in multivariate (per coordinate)."""
    s1, ok = pandrosion_h(polys, anchor, z, F_anchor, mode)
    if not ok:
        return list(z), False
    s2, ok = pandrosion_h(polys, anchor, s1, F_anchor, mode)
    if not ok:
        return s1, True
    n = len(z)
    out = []
    for k in range(n):
        d0 = s1[k] - z[k]
        d2 = s2[k] - 2 * s1[k] + z[k]
        out.append(s2[k] if abs(d2) < 1e-300 else z[k] - d0 * d0 / d2)
    return out, True


# ============================================================================
# Newton-ELS (paper 9-mv) for fallback when T_2 K=1 stalls
# ============================================================================
def newton_ELS_step(polys, z):
    n = len(z)
    F_z = F_eval(polys, z)
    if not is_finite_vec(F_z):
        return list(z), False
    h = 1e-7
    J = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        zp = list(z)
        zp[j] += h
        F_p = F_eval(polys, zp)
        for i in range(n):
            J[i][j] = (F_p[i] - F_z[i]) / h
    direction = solve_linear(J, [-v for v in F_z])
    if direction is None:
        return list(z), False
    r0 = sum(abs(v)**2 for v in F_z)
    best_z, best_r = list(z), r0
    for k in range(-6, 7):
        tau = 2.0 ** k
        cand = [z[i] + tau * direction[i] for i in range(n)]
        if not is_finite_vec(cand):
            continue
        F_c = F_eval(polys, cand)
        r = sum(abs(v)**2 for v in F_c) if is_finite_vec(cand) else float("inf")
        if r < best_r:
            best_z, best_r = cand, r
    return best_z, best_r < r0


# ============================================================================
# Armijo non-holomorphic fallback (paper 7 Definition 2.1)
# ============================================================================
def armijo_fallback(polys, z0, sigma=0.1, beta=0.5, j_max=8):
    n = len(z0)
    F0 = F_eval(polys, z0)
    if not is_finite_vec(F0):
        return list(z0), False
    F0_norm2 = sum(abs(v)**2 for v in F0)
    d_eff = max(sum(degree(p) for p in polys), 1)
    eps = max(1e-7, F0_norm2 ** (1.0 / d_eff) / (2 * d_eff))
    best_dir = None
    best_score = -float("inf")
    for j in range(n):
        for u in [1.0, -1.0, 1j, -1j]:
            zp = list(z0)
            zp[j] = z0[j] + eps * u
            F_p = F_eval(polys, zp)
            if not is_finite_vec(F_p):
                continue
            Q_steff = [(F_p[i] - F0[i]) / (eps * u) for i in range(n)]
            score = sum((F0[i].conjugate() * Q_steff[i]).real for i in range(n))
            if score > best_score and any(abs(q) > 1e-12 for q in Q_steff):
                best_score = score
                best_dir = (j, Q_steff)
    if best_dir is None:
        return list(z0), False
    j, Q_steff = best_dir
    if abs(Q_steff[j]) < 1e-12:
        return list(z0), False
    step_full = [0.0 + 0.0j] * n
    step_full[j] = -F0[j] / Q_steff[j]
    for jb in range(j_max + 1):
        b = beta ** jb
        cand = [z0[i] + b * step_full[i] for i in range(n)]
        F_c = F_eval(polys, cand)
        if not is_finite_vec(F_c):
            continue
        r2 = sum(abs(v)**2 for v in F_c)
        if r2 <= F0_norm2 * (1 - 2 * sigma * b * 0.05):
            return cand, True
    return list(z0), False


# ============================================================================
# Combined orbit with T_2 + K=1 (lag-1 secant anchor)
# ============================================================================
def combined_orbit_T2K1(polys, z_init, mode="path", tol=1e-9, max_epochs=80):
    """Orbit with T_2 (Aitken Delta-squared) and K=1 reanchoring.

    K=1 means: at each step, the anchor is the previous iterate.
    Lag-1 secant: a_n = z_{n-1}, iterate z_n. Q is the divided difference
    between consecutive iterates -- distinct points, derivative-free.
    """
    n = len(z_init)
    z = list(z_init)
    # initial anchor: small perturbation so anchor != z at step 0
    anchor = [zi + complex(0.13 * (i + 1), 0.07) for i, zi in enumerate(z_init)]
    for epoch in range(max_epochs):
        F_anchor = F_eval(polys, anchor)
        r_anchor = max(abs(v) for v in F_anchor) if is_finite_vec(F_anchor) else float("inf")
        if r_anchor < tol:
            return anchor, True
        r_z = residual_norm(polys, z)
        if r_z < tol:
            return z, True
        # T_2 step with current geometry mode
        z_t2, ok_t2 = T2_step(polys, anchor, z, F_anchor, mode)
        r_t2 = residual_norm(polys, z_t2) if ok_t2 else float("inf")
        # Newton-ELS fallback
        z_nels, ok_nels = newton_ELS_step(polys, z)
        r_nels = residual_norm(polys, z_nels) if ok_nels else float("inf")
        candidates = [(r_z, z)]
        if ok_t2 and math.isfinite(r_t2):
            candidates.append((r_t2, z_t2))
        if ok_nels and math.isfinite(r_nels):
            candidates.append((r_nels, z_nels))
        best_r, best_z = min(candidates, key=lambda x: x[0])
        if best_r > 0.99 * r_z:
            z_arm, ok_arm = armijo_fallback(polys, z)
            if ok_arm:
                r_arm = residual_norm(polys, z_arm)
                if r_arm < best_r:
                    best_r, best_z = r_arm, z_arm
        if best_r >= r_z:                       # stalled: shake anchor
            anchor = [a + complex(0.05 * random.gauss(0, 1), 0.05 * random.gauss(0, 1)) for a in anchor]
            continue
        # K=1: the new anchor is the OLD iterate (lag-1 secant structure).
        anchor = list(z)
        z = best_z
    return z, residual_norm(polys, z) < tol


# ============================================================================
# Strategy B' multivariate (paper 9-mv, same as flow/061)
# ============================================================================
def gen_starts_Bprime_mv(n, count, R=2.0, q=0.7, seed=20260427):
    rng = random.Random(seed)
    starts = []
    for k in range(count):
        u = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
        norm = sum(abs(v)**2 for v in u) ** 0.5
        if norm < 1e-12:
            u = [complex(1.0, 0.0) for _ in range(n)]
            norm = math.sqrt(n)
        u = [v / norm for v in u]
        rho = R * q ** k
        if rho < 0.05:
            rho = 0.05 + 0.5 * rng.random()
        starts.append([rho * v for v in u])
    return starts


def is_new_root(z, found, tol=1e-4):
    return all(max(abs(z[i] - r[i]) for i in range(len(z))) > tol for r in found)


# ============================================================================
# Soft repulsion deflation (same as flow/061)
# ============================================================================
def F_eval_deflated(polys, z, found_roots, sigma=0.5):
    F = F_eval(polys, z)
    if not found_roots:
        return F
    penalty = 1.0
    for r in found_roots:
        d2 = sum(abs(z[i] - r[i])**2 for i in range(len(z)))
        penalty *= (1.0 - math.exp(-d2 / (sigma**2)))
    if penalty < 1e-10:
        return [complex("inf")] * len(F)
    return [v / penalty for v in F]


def newton_ELS_deflated(polys, z, found_roots, sigma=0.5):
    n = len(z)
    F_z = F_eval_deflated(polys, z, found_roots, sigma)
    if not is_finite_vec(F_z):
        return list(z), False
    h = 1e-6
    J = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        zp = list(z)
        zp[j] += h
        F_p = F_eval_deflated(polys, zp, found_roots, sigma)
        for i in range(n):
            J[i][j] = (F_p[i] - F_z[i]) / h
    direction = solve_linear(J, [-v for v in F_z])
    if direction is None:
        return list(z), False
    r0 = sum(abs(v)**2 for v in F_z)
    best_z, best_r = list(z), r0
    for k in range(-6, 7):
        tau = 2.0 ** k
        cand = [z[i] + tau * direction[i] for i in range(n)]
        if not is_finite_vec(cand):
            continue
        F_c = F_eval_deflated(polys, cand, found_roots, sigma)
        r = sum(abs(v)**2 for v in F_c) if is_finite_vec(cand) else float("inf")
        if r < best_r:
            best_z, best_r = cand, r
    return best_z, best_r < r0


def deflated_orbit(polys, z_init, found_roots, sigma=0.4, tol=1e-9, max_epochs=60):
    z = list(z_init)
    for _ in range(max_epochs):
        if residual_norm(polys, z) < tol:
            return z, True
        z_new, ok = newton_ELS_deflated(polys, z, found_roots, sigma)
        if not ok:
            return z, False
        if max(abs(z_new[i] - z[i]) for i in range(len(z))) < 1e-12:
            return z, residual_norm(polys, z) < tol
        z = z_new
    return z, residual_norm(polys, z) < tol


# ============================================================================
# Full coverage attempt with T_2 K=1
# ============================================================================
def full_coverage_T2K1(polys, n_factor=3, n_min=12, tol=1e-9):
    n = len(next(iter(polys[0])))
    bez = bezout_estimate(polys)
    n_starts = max(n_min, n_factor * bez)
    starts = gen_starts_Bprime_mv(n, n_starts)
    found = []
    modes_used = {"path": 0, "simplex": 0, "sparse": 0, "cube": 0, "deflation": 0}
    diverged = 0
    geom_modes = ["path", "simplex", "sparse"]
    if n <= 4:
        geom_modes.append("cube")
    for z0 in starts:
        progress = False
        for mode in geom_modes:
            z, ok = combined_orbit_T2K1(polys, z0, mode, tol=tol)
            if ok and is_new_root(z, found):
                found.append(z)
                modes_used[mode] += 1
                progress = True
                break
        if not progress:
            diverged += 1
    if len(found) < bez:
        rng = random.Random(20260428)
        n_extra = max(20, bez * 2)
        for trial in range(n_extra):
            if len(found) >= bez:
                break
            z0 = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
            z, ok = deflated_orbit(polys, z0, found, sigma=0.4, tol=tol)
            if ok and is_new_root(z, found):
                found.append(z)
                modes_used["deflation"] += 1
    return {
        "roots": found,
        "n_starts": n_starts,
        "bezout": bez,
        "diverged": diverged,
        "coverage": len(found) / max(bez, 1),
        "modes": modes_used,
    }


# ============================================================================
# Test systems (same generator as 057, 059, 060, 061)
# ============================================================================
def gen_random_poly_system(n, d, seed):
    rng = random.Random(seed)
    polys = []
    for i in range(n):
        poly = {}
        for alpha in iprod(*[range(d + 1) for _ in range(n)]):
            if sum(alpha) > d:
                continue
            if rng.random() < 0.7:
                poly[tuple(alpha)] = complex(rng.gauss(0.0, 1.0), 0.15 * rng.gauss(0.0, 1.0))
        diag = tuple(d if k == i else 0 for k in range(n))
        poly[diag] = poly.get(diag, 0.0 + 0.0j) + 1.0
        m = max(abs(v) for v in poly.values())
        polys.append({k: v / m for k, v in poly.items()})
    return polys


def main():
    print("=" * 116)
    print("flow/062 -- T_2 K=1 (paper 8 empirical optimum) + universal geometries + Armijo + deflation")
    print("=" * 116)
    print("Distinct from 061: T_2 (not T_3), K=1 reanchor (not K=3), lag-1 secant anchor.")
    print("Goal: assess if T_2 K=1 changes Bezout coverage vs 061's T_3 K=3.")
    print()
    # Same seed pattern as 061 to allow direct comparison.
    cases = [(2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 2), (3, 3), (4, 2)]
    print(f"  {'(n,d)':>8} {'Bez':>5} {'starts':>7} | "
          f"{'roots':>5} {'cov%':>6} {'modes used':>40} {'time':>7}")
    print("-" * 116)
    for n, d in cases:
        polys = gen_random_poly_system(n, d, seed=61000 + 100 * n + d)  # SAME seed as 061
        t0 = time.time()
        result = full_coverage_T2K1(polys, n_factor=3, n_min=12)
        elapsed = time.time() - t0
        cov = 100.0 * result["coverage"]
        modes = " ".join(f"{k}:{v}" for k, v in result["modes"].items() if v > 0)
        print(f"  ({n:>2},{d:>2}) {result['bezout']:>5} {result['n_starts']:>7} | "
              f"{len(result['roots']):>5} {cov:>5.1f}% {modes:>40} {elapsed:>6.1f}s")
    print()
    print("=" * 116)
    print("VERDICT (compare with 061 same seed)")
    print("=" * 116)
    print("  If 062 coverage > 061 by >5pp on most cells: T_2 K=1 helps coverage.")
    print("  If 062 coverage similar to 061: T_2 K=1 affects orbit speed only,")
    print("                                 not Bezout coverage (basin geometry).")
    print("  If 062 coverage < 061: lag-1 secant has more aggressive dynamics,")
    print("                         may overshoot/escape some basins.")
    print("=" * 116)


if __name__ == "__main__":
    main()
