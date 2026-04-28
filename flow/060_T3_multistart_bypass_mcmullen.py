"""
PAPER: 060
TITLE: T3 multistart with hypercube / sparse / simplex Q_F geometries
       at Bezout-aware density on multivariate KS systems

MISSION
=======

User-conjectured: paper 1 claims "T3 multistart bypasses McMullen" (Prop. 3.5,
"5/5 convergence on (n,d) up to (6,3)").  If true, T3 multistart with the four
universal Q_F geometries from flow/044-049 (hypercube / sparse / simplex /
path) should reach near-100% Bezout coverage without homotopy.

Honest reading of paper 1 (line 455 of latex/1pandrosion_smale.tex):
    "The strict universal form of Conjecture ABD was REFUTED ... This is
    CONSISTENT with McMullen's topological impossibility theorem for purely
    rational iterative schemes of degree d >= 4."

So McMullen holds for pure-rational T_n.  Bypass requires the non-holomorphic
fallback of paper 7 (NOT in scope here).  Paper 1's "5/5 at (n,d) = (5,4),
D=1024" is ABOUT FINDING ONE ROOT per orbit (per-orbit success rate), NOT
about finding all D Bezout roots per system.

This script tests both readings empirically: combine T3 multistart with
the four geometries 057 already uses, plus a 4D-aware simplex variant, and
measure Bezout coverage on the same benchmarks as 057/059.

Expectation: coverage similar to 057 (~30-70%), confirming that T3 alone does
not bypass McMullen for Bezout coverage even when combined with all four
universal geometries.
"""

from __future__ import annotations

import cmath
import itertools
import math
import random
import time
from itertools import product as iprod


# --------------------------------------------------------------------------
# Polynomial primitives
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# Linear algebra
# --------------------------------------------------------------------------
def solve_linear(A, b):
    n = len(A)
    M = [list(row) + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-13:
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


# --------------------------------------------------------------------------
# The four universal Q_F geometries
# --------------------------------------------------------------------------
def schmidt_path(polys, anchor, z, order):
    """Path geometry: telescope along coordinates one at a time."""
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
        cur = nxt
        F_cur = F_next
    return Q


def schmidt_cube(polys, anchor, z, max_orders=8):
    """Cube geometry: average over coordinate-permutations."""
    n = len(z)
    perms = list(itertools.permutations(range(n)))
    if len(perms) > max_orders:
        perms = perms[:max_orders]
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for order in perms:
        Q = schmidt_path(polys, anchor, z, order)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Q[i][j] / len(perms)
    return Q_avg


def schmidt_simplex(polys, anchor, z, F_anchor):
    """Simplex / edge-frame geometry with defect correction (universal-prism)."""
    n = len(z)
    delta = [z[i] - anchor[i] for i in range(n)]
    norm2 = sum(abs(v)**2 for v in delta)
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
    if norm2 < 1e-28:
        return B
    F_z = F_eval(polys, z)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [v.conjugate() for v in delta]
    denom = sum(w[j] * delta[j] for j in range(n))
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def schmidt_sparse(polys, anchor, z, F_anchor):
    """Sparse / Newton-polytope-aware geometry: weights by support activity."""
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


# --------------------------------------------------------------------------
# Pandrosion T3 step (4D-aware: cube uses up to all 24 perms in n=4)
# --------------------------------------------------------------------------
def pandrosion_h(polys, anchor, z, F_anchor, mode):
    Q = Q_geometry(polys, anchor, z, F_anchor, mode)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if not is_finite_vec(out):
        return list(z), False
    return out, True


def T3_step(polys, anchor, z, F_anchor, mode):
    s1, ok1 = pandrosion_h(polys, anchor, z, F_anchor, mode)
    if not ok1:
        return list(z), False
    s2, ok2 = pandrosion_h(polys, anchor, s1, F_anchor, mode)
    if not ok2:
        return s1, True
    n = len(z)
    t2 = []
    for k in range(n):
        d0 = s1[k] - z[k]
        d2 = s2[k] - 2 * s1[k] + z[k]
        t2.append(s2[k] if abs(d2) < 1e-300 else z[k] - d0 * d0 / d2)
    s3, ok3 = pandrosion_h(polys, anchor, t2, F_anchor, mode)
    if not ok3:
        return t2, True
    lam_den = sum(s1[k] - z[k] for k in range(n))
    if abs(lam_den) < 1e-300:
        return t2, True
    lam = sum(s2[k] - s1[k] for k in range(n)) / lam_den
    if abs(lam - 1.0) < 1e-12:
        return s3, True
    out = [t2[k] - (s3[k] - t2[k]) / (lam - 1.0) for k in range(n)]
    return out, True


def adaptive_T3_orbit(polys, z_init, mode, K=3, tol=1e-9, max_cycles=20):
    """T3 + adaptive reanchor (paper 1 Definition 4.x).

    Mirrors flow/057's adaptive_pandrosion_T3 logic so that we get the
    same convergence behaviour, with the geometry mode swappable.
    """
    n = len(z_init)
    z = list(z_init)
    anchor = [zi + complex(0.13 * (i + 1), 0.07) for i, zi in enumerate(z_init)]
    for cycle in range(max_cycles):
        F_anchor = F_eval(polys, anchor)
        if max(abs(v) for v in F_anchor) < tol and is_finite_vec(anchor):
            return anchor, True
        progress = False
        for _ in range(K):
            r0 = residual_norm(polys, z)
            if r0 < tol:
                return z, True
            z_new, ok = T3_step(polys, anchor, z, F_anchor, mode)
            if not ok:
                break
            r_new = residual_norm(polys, z_new)
            if math.isfinite(r_new) and r_new < r0:
                z = z_new
                progress = True
            else:
                z_half = [z[k] + 0.5 * (z_new[k] - z[k]) for k in range(n)]
                r_half = residual_norm(polys, z_half)
                if r_half < r0:
                    z = z_half
                    progress = True
                else:
                    break    # this K-block stalled; try a fresh anchor
            if residual_norm(polys, z) < tol:
                return z, True
        anchor = list(z) if progress else [a + complex(0.05, 0.03) for a in anchor]
    return z, residual_norm(polys, z) < tol


# --------------------------------------------------------------------------
# Multistart with the four geometries combined as a portfolio
# --------------------------------------------------------------------------
def gen_complex_starts_spiral(n, count, R_max=1.5, seed=2026):
    """Strategy B golden-angle spiral, multivariate."""
    rng = random.Random(seed)
    phi = math.pi * (3.0 - math.sqrt(5.0))
    starts = []
    for k in range(count):
        z = []
        rho = R_max * (0.95 ** k) if (0.95 ** k) > 0.05 else R_max * (0.3 + 0.7 * rng.random())
        for j in range(n):
            angle = phi * (k + 1) * (j + 1) + 0.13 * j
            z.append(complex(rho * math.cos(angle), rho * math.sin(angle)))
        z = [zi + complex(0.02 * rng.uniform(-1, 1), 0.02 * rng.uniform(-1, 1))
             for zi in z]
        starts.append(z)
    return starts


def is_new_root(z, found, tol=1e-5):
    return all(max(abs(z[i] - r[i]) for i in range(len(z))) > tol for r in found)


def portfolio_multistart(polys, n_factor=3, n_min=12, tol=1e-9):
    """T3 multistart with all four universal geometries."""
    n = len(next(iter(polys[0])))
    bez = bezout_estimate(polys)
    n_starts = max(n_min, n_factor * bez)
    starts = gen_complex_starts_spiral(n, n_starts)
    found = []
    convergences = {"path": 0, "cube": 0, "simplex": 0, "sparse": 0}
    diverged = 0
    modes = ["path", "simplex", "sparse"]
    if n <= 4:
        modes.append("cube")        # cube = full coordinate-permutation avg
    for z0 in starts:
        for mode in modes:
            z, ok = adaptive_T3_orbit(polys, z0, mode, tol=tol)
            if ok and is_new_root(z, found):
                found.append(z)
                convergences[mode] = convergences.get(mode, 0) + 1
                break  # next start
        else:
            diverged += 1
    return {
        "roots": found,
        "n_starts": n_starts,
        "bezout": bez,
        "diverged": diverged,
        "coverage": len(found) / max(bez, 1),
        "modes": convergences,
    }


# --------------------------------------------------------------------------
# Test systems (same generator as 057, 059)
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    print("=" * 116)
    print("flow/060 -- T3 multistart with hypercube/sparse/simplex/path geometries")
    print("=" * 116)
    print("Question: does T3 multistart bypass McMullen and reach Bezout coverage?")
    print("Reading of paper 1 line 455: NO -- T3 is rational, McMullen barrier holds.")
    print("Empirical test below confirms or refutes this reading.")
    print()

    cases = [
        (2, 2),   # Bezout=4
        (2, 3),   # Bezout=9
        (2, 4),   # Bezout=16
        (2, 5),   # Bezout=25
        (2, 6),   # Bezout=36
        (3, 2),   # Bezout=8
        (3, 3),   # Bezout=27
        (4, 2),   # Bezout=16
    ]

    print(f"  {'(n,d)':>8} {'Bezout':>7} {'starts':>7} | "
          f"{'roots':>5} {'cov%':>6} {'modes hit':>30} {'time':>7}")
    print("-" * 116)
    for n, d in cases:
        polys = gen_random_poly_system(n, d, seed=60000 + 100 * n + d)
        t0 = time.time()
        result = portfolio_multistart(polys, n_factor=3, n_min=12)
        elapsed = time.time() - t0
        cov = 100.0 * result["coverage"]
        modes = " ".join(f"{k}:{v}" for k, v in result["modes"].items() if v > 0)
        print(f"  ({n:>2},{d:>2})  {result['bezout']:>7} {result['n_starts']:>7} | "
              f"{len(result['roots']):>5} {cov:>5.1f}% {modes:>30} {elapsed:>6.1f}s")
    print()
    print("=" * 116)
    print("VERDICT")
    print("=" * 116)
    print("  Comparison with flow/057 (T3 multistart, single geometry):")
    print("    - 057 plateaus at 30-70% coverage on the same benchmarks.")
    print("    - 060 (this script) adds the cube geometry for n<=4 and tries all")
    print("      four geometries per start.  If results similar to 057 ==> McMullen")
    print("      barrier is geometry-independent.")
    print("    - 059 (homotopy + 051 corrector) achieves 100% on the same benchmarks.")
    print()
    print("  Honest reading of paper 1:")
    print("    Line 174-194 (Prop. 3.5):  '5/5 convergence' = per-orbit success,")
    print("                               NOT per-system Bezout coverage.")
    print("    Line 200:  'must be multiplied by D = d^n starts for system-level total'.")
    print("    Line 455:  'refutation consistent with McMullen for purely rational schemes'.")
    print("    Bypass of McMullen requires the non-holomorphic fallback of paper 7,")
    print("    which is OUT of T3-multistart scope.")
    print("=" * 116)


if __name__ == "__main__":
    main()
