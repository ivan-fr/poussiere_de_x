"""
PAPER: 006 (acceleration piste 1)
TITLE: Asymmetric Pandrosion dD --- Multi-Target Roots via Brezinski-Vector
STATUS: empirical test; compares vs n separate Pandrosion 2D + Steffensen

PISTE 1 (acceleration via shared evaluations)
=============================================
Compute n simultaneous roots x_1^{1/p}, ..., x_n^{1/p} via the asymmetric
Pandrosion dD iteration (no rank-1 collapse: x_i distinct breaks the
diagonal symmetry).  Each step shares ONE evaluation of S_p^{dD}, then
each component s_i is updated separately:

    s_{i, n+1} = 1 - (x_i - 1) / (x_i * S_p^{dD}(s_n))

After n iterations, apply Brezinski-vector-Aitken (epsilon_2 = vector
Aitken) to accelerate the coupled vector sequence (s_1, ..., s_{d-1}).

Baseline: solve each x_i^{1/p} independently with Pandrosion 2D + Steffensen,
which costs ~5 evaluations per root, so 5n total.  Pandrosion dD asymmetric
+ Brezinski costs ~5 evaluations of S_p^{dD} (shared).

Expected gain: x(n-1)/n on TOTAL cost if the asymmetric coupling converges
at comparable rate.

Verifies:
  1. Asymmetric dD iteration converges to the right roots.
  2. Brezinski-vector-Aitken accelerates beyond linear.
  3. Wall-clock comparison: dD-asym + Brezinski vs 2D x n + Steffensen.
"""
from __future__ import annotations
import math
import time
from math import comb

import numpy as np


def S_dD(s_vec, p: int) -> complex:
    s_vec = [complex(s) for s in s_vec]
    n = len(s_vec)
    total = 0.0 + 0.0j
    def gen(rem, depth):
        if depth == 0:
            yield (); return
        for k in range(rem + 1):
            for r in gen(rem - k, depth - 1):
                yield (k,) + r
    for alpha in gen(p - 1, n):
        term = 1.0 + 0.0j
        for i, a in enumerate(alpha):
            if a:
                term *= s_vec[i] ** a
        total += term
    return total


def pandrosion_2d_steffensen(x: float, p: int = 2, max_iter=20, tol=1e-13):
    """Baseline: Pandrosion 2D + Steffensen for single x^{1/p}."""
    def g(s):
        S = sum(s**k for k in range(p))
        return 1 - (x - 1) / (x * S)
    s = 0.5
    n_evals = 0
    for _ in range(max_iter):
        g1 = g(s); n_evals += 1
        g2 = g(g1); n_evals += 1
        denom = g2 - 2*g1 + s
        if abs(denom) < 1e-300:
            return g2, n_evals
        s_new = s - (g1 - s)**2 / denom
        if abs(s_new - s) < tol:
            return s_new, n_evals
        s = s_new
    return s, n_evals


def samelson_inverse(v):
    norm_sq = float(np.vdot(v, v).real)
    if norm_sq < 1e-30:
        return np.zeros_like(v)
    return np.conj(v) / norm_sq


def brezinski_vector_step(s_minus, s_zero, s_plus):
    """Brezinski-Aitken vector Delta^2 = Wynn epsilon_2."""
    ds = s_zero - s_minus
    ds_p = s_plus - s_zero
    diff = ds_p - ds
    norm_sq = float(np.vdot(diff, diff).real)
    if norm_sq < 1e-30:
        return s_plus
    lambda_v = np.vdot(diff, ds_p).real / norm_sq
    return s_plus - lambda_v * ds_p


def pandrosion_dD_asym_brezinski(x_vec, p: int, max_iter=30, tol=1e-13):
    """Asymmetric Pandrosion dD with Brezinski-vector acceleration.

    Note: with x_i distinct, NO rank-1 collapse; the iteration is genuinely
    multivariate.  We track 3 successive iterates and apply Brezinski every
    3 steps.
    """
    n = len(x_vec)
    s = np.array([0.5 + 0.0j] * n, dtype=complex)
    n_evals = 0
    for outer in range(max_iter):
        # Three Pandrosion dD iterations (one S_p^{dD} eval each)
        s_minus = s.copy()
        Spsn = S_dD(s, p); n_evals += 1
        s_zero = np.array([1 - (x_vec[i] - 1) / (x_vec[i] * Spsn) for i in range(n)])
        Spsn = S_dD(s_zero, p); n_evals += 1
        s_plus = np.array([1 - (x_vec[i] - 1) / (x_vec[i] * Spsn) for i in range(n)])
        # Brezinski vector acceleration
        s_acc = brezinski_vector_step(s_minus, s_zero, s_plus)
        if np.max(np.abs(s_acc - s_plus)) < tol:
            return s_acc, n_evals
        s = s_acc
    return s, n_evals


def pandrosion_dD_asym_raw(x_vec, p: int, max_iter=300, tol=1e-13):
    """Asymmetric Pandrosion dD without acceleration (linear baseline)."""
    n = len(x_vec)
    s = np.array([0.5 + 0.0j] * n, dtype=complex)
    n_evals = 0
    for _ in range(max_iter):
        Spsn = S_dD(s, p); n_evals += 1
        s_new = np.array([1 - (x_vec[i] - 1) / (x_vec[i] * Spsn) for i in range(n)])
        if np.max(np.abs(s_new - s)) < tol:
            return s_new, n_evals
        s = s_new
    return s, n_evals


def main():
    print("=" * 76)
    print("PISTE 1 -- Asymmetric Pandrosion dD + Brezinski-vector multi-target")
    print("=" * 76)
    test_targets = [
        (2.0, 3.0, 5.0),                             # n = 3
        (1.5, 2.5, 3.5, 4.5),                        # n = 4
        (2.0, 3.0, 5.0, 7.0, 11.0),                  # n = 5
        (1.2, 1.5, 2.0, 3.0, 5.0, 8.0),              # n = 6
    ]
    for x_tuple in test_targets:
        n = len(x_tuple)
        # Baseline: n separate 2D + Steffensen
        t0 = time.perf_counter()
        n_evals_baseline = 0
        roots_2d = []
        for x in x_tuple:
            s_star, n_e = pandrosion_2d_steffensen(x, p=2)
            roots_2d.append(s_star)
            n_evals_baseline += n_e
        t_baseline = time.perf_counter() - t0
        # Method dD asymmetric raw (linear)
        t0 = time.perf_counter()
        s_dD_raw, ne_raw = pandrosion_dD_asym_raw(x_tuple, p=2)
        t_dD_raw = time.perf_counter() - t0
        # Method dD asymmetric + Brezinski
        t0 = time.perf_counter()
        s_dD_brez, ne_brez = pandrosion_dD_asym_brezinski(x_tuple, p=2)
        t_dD_brez = time.perf_counter() - t0
        print(f"\nTargets x = {x_tuple}, n = {n}")
        print(f"  baseline 2D + Steffensen, n separate roots:")
        print(f"    {n_evals_baseline} S_p evals total, time {t_baseline*1e6:.1f} us")
        print(f"  asymmetric dD (raw, linear):")
        print(f"    {ne_raw} S_p^{{dD}} evals, time {t_dD_raw*1e6:.1f} us")
        print(f"  asymmetric dD + Brezinski:")
        print(f"    {ne_brez} S_p^{{dD}} evals, time {t_dD_brez*1e6:.1f} us")
        # Check correctness: each s_i should give x_i^{-1/2} (since p=2, fixed point of paper-0 is 1/sqrt(x))
        # But asymmetric dD has different fixed point structure -- let's just check residual on the dD eq
        Sp = S_dD(s_dD_brez, 2)
        max_residual = max(abs(x_tuple[i] * Sp * (1 - s_dD_brez[i]) - (x_tuple[i] - 1))
                            for i in range(n))
        print(f"    max dD residual: {max_residual:.2e}")
        # Practical question: do the dD roots correspond to the 2D roots?
        # 2D Pandrosion solves x*(1+s)*(1-s) = x-1 -> s = 1/sqrt(x)
        # dD asymmetric solves x_i * S_p^{dD}(s) * (1-s_i) = x_i - 1
        # These are DIFFERENT problems!  The asymmetric dD does NOT compute x_i^{1/p}.
        print(f"    NOTE: dD asym roots are NOT x_i^{{1/p}}; they solve a different system.")


if __name__ == "__main__":
    main()
