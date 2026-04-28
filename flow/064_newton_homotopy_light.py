"""
PAPER: 064
TITLE: Newton homotopy light (no gamma rotation) with Pandrosion T_2 K=1 corrector
STATUS: empirical; this IS homotopy, but minimal (single path, no Lairez safety)

DIFFERENCE FROM flow/059
========================
flow/059 uses gamma * F (random rotation) and TOTAL-DEGREE start system
G_i = z_i^{d_i} - 1.  The gamma rotation provides the LAIREZ SAFETY:
under generic gamma, all paths H_t = (1-t)G + t gamma F are smooth and
disjoint for t in [0, 1).

flow/064 (this script) drops the gamma rotation.  It uses
    H_t = (1 - t) G + t F,
i.e., the straight-line homotopy from G to F.  No gamma safety.

WHY DROP GAMMA?
===============
The user's question is whether a Pandrosion-pure approach (with paper 7
fallback + universal geometries) can get full coverage.  Including
gamma is "doing Lairez homotopy"; including only the straight-line
H_t = (1-t)G + t F is a more minimal homotopy (Newton-homotopy /
classical continuation) which may FAIL on a measure-zero set of F where
the path crosses a singularity (path-jumping).

In practice, on random KS systems, the singular set is generically
avoided by the straight-line homotopy.  We expect 100% coverage on
most cases, but possibly some path failures where the gamma trick
would have rescued.

CORRECTOR
=========
Per path, we use the Pandrosion T_2 K=1 portfolio (4 geometries +
Armijo fallback) of flow/062.  No symbolic Jacobian unless the
fallback fires.

EXPECTATION
===========
- Coverage close to 100% on (2,2)..(4,2).
- Marginal failures (1-2 paths out of D) when straight-line homotopy
  hits a singularity.  These would be rescued by gamma rotation
  (= flow/059), at the cost of admitting "we did Lairez".
- This is HONESTLY homotopy, just with a less safe starting point
  than 059.  The flow/064 result is a measure of "what does even
  minimal homotopy buy us?".
"""

from __future__ import annotations

import cmath
import itertools
import math
import random
import time
from itertools import product as iprod


# ============================================================================
# Polynomial primitives (same as 061-063)
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
# Linear algebra and 4 geometries (same as 062)
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


def schmidt_path(polys, anchor, z, order):
    n = len(z); Q = [[0.0 + 0.0j] * n for _ in range(n)]
    cur = list(z); F_cur = F_eval(polys, cur)
    for j in order:
        nxt = list(cur); nxt[j] = anchor[j]
        F_next = F_eval(polys, nxt)
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300: continue
        for i in range(n):
            Q[i][j] = (F_cur[i] - F_next[i]) / denom
        cur, F_cur = nxt, F_next
    return Q


def schmidt_cube(polys, anchor, z, max_orders=8):
    n = len(z); perms = list(itertools.permutations(range(n)))[:max_orders]
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for order in perms:
        Q = schmidt_path(polys, anchor, z, order)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Q[i][j] / len(perms)
    return Q_avg


def edge_frame(polys, anchor, z, F_anchor):
    n = len(z); B = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        edge = list(anchor); edge[j] = z[j]
        F_edge = F_eval(polys, edge); denom = z[j] - anchor[j]
        if abs(denom) < 1e-300: continue
        for i in range(n):
            B[i][j] = (F_edge[i] - F_anchor[i]) / denom
    return B


def schmidt_simplex(polys, anchor, z, F_anchor):
    n = len(z); delta = [z[i] - anchor[i] for i in range(n)]
    B = edge_frame(polys, anchor, z, F_anchor)
    if sum(abs(v)**2 for v in delta) < 1e-28: return B
    F_z = F_eval(polys, z)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [v.conjugate() for v in delta]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28: return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def schmidt_sparse(polys, anchor, z, F_anchor):
    n = len(z); activity = [0.0] * n
    for poly in polys:
        for exp, coeff in poly.items():
            for j, e in enumerate(exp):
                activity[j] += abs(coeff) * e
    m = max(activity) if activity else 0.0
    activity = [a if a > 0 else max(1.0, 0.1 * m) for a in activity]
    delta = [z[i] - anchor[i] for i in range(n)]
    B = edge_frame(polys, anchor, z, F_anchor)
    F_z = F_eval(polys, z)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [activity[j] * delta[j].conjugate() for j in range(n)]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28: return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def Q_geometry(polys, anchor, z, F_anchor, mode):
    if mode == "path": return schmidt_path(polys, anchor, z, tuple(range(len(z))))
    if mode == "cube": return schmidt_cube(polys, anchor, z)
    if mode == "simplex": return schmidt_simplex(polys, anchor, z, F_anchor)
    if mode == "sparse": return schmidt_sparse(polys, anchor, z, F_anchor)
    raise ValueError(mode)


def pandrosion_h(polys, anchor, z, F_anchor, mode):
    Q = Q_geometry(polys, anchor, z, F_anchor, mode)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None: return list(z), False
    out = [anchor[i] + step[i] for i in range(len(z))]
    return (out, True) if is_finite_vec(out) else (list(z), False)


def T2_step(polys, anchor, z, F_anchor, mode):
    s1, ok = pandrosion_h(polys, anchor, z, F_anchor, mode)
    if not ok: return list(z), False
    s2, ok = pandrosion_h(polys, anchor, s1, F_anchor, mode)
    if not ok: return s1, True
    n = len(z); out = []
    for k in range(n):
        d0 = s1[k] - z[k]; d2 = s2[k] - 2 * s1[k] + z[k]
        out.append(s2[k] if abs(d2) < 1e-300 else z[k] - d0 * d0 / d2)
    return out, True


def newton_ELS_step(polys, z):
    n = len(z); F_z = F_eval(polys, z)
    if not is_finite_vec(F_z): return list(z), False
    h = 1e-7
    J = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        zp = list(z); zp[j] += h; F_p = F_eval(polys, zp)
        for i in range(n):
            J[i][j] = (F_p[i] - F_z[i]) / h
    direction = solve_linear(J, [-v for v in F_z])
    if direction is None: return list(z), False
    r0 = sum(abs(v)**2 for v in F_z); best_z, best_r = list(z), r0
    for k in range(-6, 7):
        tau = 2.0 ** k
        cand = [z[i] + tau * direction[i] for i in range(n)]
        if not is_finite_vec(cand): continue
        F_c = F_eval(polys, cand)
        r = sum(abs(v)**2 for v in F_c) if is_finite_vec(cand) else float("inf")
        if r < best_r: best_z, best_r = cand, r
    return best_z, best_r < r0


def armijo_fallback(polys, z0, sigma=0.1, beta=0.5, j_max=8):
    n = len(z0); F0 = F_eval(polys, z0)
    if not is_finite_vec(F0): return list(z0), False
    F0_norm2 = sum(abs(v)**2 for v in F0)
    d_eff = max(sum(degree(p) for p in polys), 1)
    eps = max(1e-7, F0_norm2 ** (1.0 / d_eff) / (2 * d_eff))
    best_dir, best_score = None, -float("inf")
    for j in range(n):
        for u in [1.0, -1.0, 1j, -1j]:
            zp = list(z0); zp[j] = z0[j] + eps * u
            F_p = F_eval(polys, zp)
            if not is_finite_vec(F_p): continue
            Q_steff = [(F_p[i] - F0[i]) / (eps * u) for i in range(n)]
            score = sum((F0[i].conjugate() * Q_steff[i]).real for i in range(n))
            if score > best_score and any(abs(q) > 1e-12 for q in Q_steff):
                best_score, best_dir = score, (j, Q_steff)
    if best_dir is None: return list(z0), False
    j, Q_steff = best_dir
    if abs(Q_steff[j]) < 1e-12: return list(z0), False
    step_full = [0.0 + 0.0j] * n
    step_full[j] = -F0[j] / Q_steff[j]
    for jb in range(j_max + 1):
        b = beta ** jb
        cand = [z0[i] + b * step_full[i] for i in range(n)]
        F_c = F_eval(polys, cand)
        if not is_finite_vec(F_c): continue
        r2 = sum(abs(v)**2 for v in F_c)
        if r2 <= F0_norm2 * (1 - 2 * sigma * b * 0.05):
            return cand, True
    return list(z0), False


# ============================================================================
# Newton-homotopy primitives
# ============================================================================
def add_scaled(out, poly, scale):
    for exp, coeff in poly.items():
        out[exp] = out.get(exp, 0.0 + 0.0j) + scale * coeff


def clean_dict(d, eps=1e-14):
    return {k: v for k, v in d.items() if abs(v) > eps}


def homotopy_polys(target, start, t):
    """H_t = (1 - t) * G + t * F.  No gamma rotation -- straight-line homotopy."""
    polys = []
    for f, g in zip(target, start):
        h = {}
        add_scaled(h, g, 1.0 - t)
        add_scaled(h, f, t)
        polys.append(clean_dict(h))
    return polys


def total_degree_start_system(degrees, n):
    """G_i(z) = z_i^{d_i} - 1.  Roots are products of d_i-th roots of unity."""
    polys = []
    zero = tuple([0] * n)
    for i, d in enumerate(degrees):
        exp = tuple(d if j == i else 0 for j in range(n))
        polys.append({exp: 1.0 + 0.0j, zero: -1.0 + 0.0j})
    return polys


def roots_of_unity(d):
    return [cmath.exp(2j * math.pi * k / d) for k in range(d)]


def start_roots(degrees):
    return [list(root) for root in iprod(*[roots_of_unity(d) for d in degrees])]


# ============================================================================
# Path tracker with Pandrosion T_2 K=1 + 4-geometry portfolio + Armijo
# ============================================================================
def correct_T2K1_portfolio(polys, z_init, tol=1e-9, max_epochs=30):
    """Combined T_2 K=1 corrector with portfolio of 4 geometries + Armijo."""
    n = len(z_init); z = list(z_init)
    anchor = [zi + complex(0.13 * (i + 1), 0.07) for i, zi in enumerate(z_init)]
    geom_modes = ["path", "simplex", "sparse"]
    if n <= 4: geom_modes.append("cube")
    for epoch in range(max_epochs):
        F_anchor = F_eval(polys, anchor)
        if is_finite_vec(F_anchor) and max(abs(v) for v in F_anchor) < tol:
            return anchor, True
        r_z = residual_norm(polys, z)
        if r_z < tol:
            return z, True
        # Try each geometry mode in T_2 step
        best_r = r_z; best_z = z
        for mode in geom_modes:
            z_t2, ok_t2 = T2_step(polys, anchor, z, F_anchor, mode)
            if ok_t2:
                r_t2 = residual_norm(polys, z_t2)
                if math.isfinite(r_t2) and r_t2 < best_r:
                    best_r, best_z = r_t2, z_t2
        # Newton-ELS as backup
        z_n, ok_n = newton_ELS_step(polys, z)
        if ok_n:
            r_n = residual_norm(polys, z_n)
            if math.isfinite(r_n) and r_n < best_r:
                best_r, best_z = r_n, z_n
        # Armijo if no descent
        if best_r > 0.99 * r_z:
            z_a, ok_a = armijo_fallback(polys, z)
            if ok_a:
                r_a = residual_norm(polys, z_a)
                if r_a < best_r: best_r, best_z = r_a, z_a
        if best_r >= r_z:
            anchor = [a + complex(0.05*random.gauss(0,1), 0.05*random.gauss(0,1)) for a in anchor]
            continue
        anchor = list(z); z = best_z
    return z, residual_norm(polys, z) < tol


def track_one_path(target, start, z0, tol=1e-9):
    """Track a single homotopy path from t=0 (z=z0 root of start system)
    to t=1 (z = root of target).  Predictor + corrector."""
    t = 0.0; dt = 0.005; z = list(z0); n = len(z)
    prev_z, prev_t = None, None
    failures = 0
    while t < 1.0 - 1e-15:
        trial_dt = min(dt, 1.0 - t)
        t_next = t + trial_dt
        # Predictor: linear extrapolation in (z, t)
        if prev_z is None:
            pred = list(z)
        else:
            slope = [(z[i] - prev_z[i]) / (t - prev_t) for i in range(n)]
            pred = [z[i] + trial_dt * slope[i] for i in range(n)]
        # Corrector
        H = homotopy_polys(target, start, t_next)
        z_corr, ok = correct_T2K1_portfolio(H, pred, tol=tol, max_epochs=20)
        if ok and residual_norm(H, z_corr) < 10.0 * tol:
            prev_z, prev_t = list(z), t
            z, t = z_corr, t_next
            dt = min(0.02, dt * 1.1)
        else:
            failures += 1
            dt *= 0.5
            if dt < 1e-5 or failures > 60:
                return z, False
    final_r = residual_norm(target, z)
    return z, final_r < 1e-7


def track_all_roots(target, tol=1e-9):
    n = len(target)
    degrees = [max(1, degree(p)) for p in target]
    start_sys = total_degree_start_system(degrees, n)
    z0_list = start_roots(degrees)
    found = []
    paths_ok, paths_fail = 0, 0
    for z0 in z0_list:
        z, ok = track_one_path(target, start_sys, z0, tol=tol)
        if ok:
            paths_ok += 1
            if all(max(abs(z[i] - r[i]) for i in range(n)) > 1e-4 for r in found):
                found.append(z)
        else:
            paths_fail += 1
    bez = math.prod(degrees)
    return {"roots": found, "bezout": bez, "paths_ok": paths_ok,
            "paths_fail": paths_fail, "coverage": len(found) / max(bez, 1)}


def gen_random_poly_system(n, d, seed):
    rng = random.Random(seed); polys = []
    for i in range(n):
        poly = {}
        for alpha in iprod(*[range(d + 1) for _ in range(n)]):
            if sum(alpha) > d: continue
            if rng.random() < 0.7:
                poly[tuple(alpha)] = complex(rng.gauss(0.0, 1.0), 0.15 * rng.gauss(0.0, 1.0))
        diag = tuple(d if k == i else 0 for k in range(n))
        poly[diag] = poly.get(diag, 0.0 + 0.0j) + 1.0
        m = max(abs(v) for v in poly.values())
        polys.append({k: v / m for k, v in poly.items()})
    return polys


def main():
    print("=" * 116)
    print("flow/064 -- Newton homotopy light (no gamma) with Pandrosion T_2 K=1 corrector")
    print("=" * 116)
    print("Straight-line H_t = (1-t)*G + t*F, G_i = z_i^{d_i} - 1, no Lairez gamma.")
    print("Honest about being homotopy: this IS path tracking, just the minimal version.")
    print()
    cases = [(2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 2), (3, 3), (4, 2)]
    print(f"  {'(n,d)':>8} {'Bez':>5} | {'roots':>5} {'cov%':>6} "
          f"{'paths ok/fail':>14} {'time':>7}")
    print("-" * 116)
    for n, d in cases:
        polys = gen_random_poly_system(n, d, seed=61000 + 100 * n + d)  # SAME seed
        t0 = time.time()
        result = track_all_roots(polys)
        elapsed = time.time() - t0
        cov = 100.0 * result["coverage"]
        print(f"  ({n:>2},{d:>2}) {result['bezout']:>5} | "
              f"{len(result['roots']):>5} {cov:>5.1f}% "
              f"{result['paths_ok']:>4}/{result['paths_fail']:<4} {elapsed:>6.1f}s")
    print()
    print("=" * 116)
    print("VERDICT")
    print("=" * 116)
    print("  Compare with:")
    print("    flow/057 multistart only:                 30-70%")
    print("    flow/061 T3 K=3 + everything:             56-100%")
    print("    flow/062 T2 K=1 + everything:             56-100%")
    print("    flow/063 T2 K=1 + symmetric targeting:    63-100%")
    print("    flow/064 (this) Newton-homotopy + T2 K=1: ?-100% (expect near 100%)")
    print("    flow/059 full gamma-homotopy + 051:       100% always")
    print()
    print("  064 IS homotopy.  If it reaches 100%, that confirms: homotopy structure")
    print("  (not the corrector or the gamma rotation) is what gives Bezout coverage.")
    print("=" * 116)


if __name__ == "__main__":
    main()
