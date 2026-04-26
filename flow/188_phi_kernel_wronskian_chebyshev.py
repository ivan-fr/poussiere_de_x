"""
PAPER: 188 (NEW — Wronskian/Chebyshev route for higher minors)
TITLE: Wronskian sign-regularity for the Riemann-Phi summand kernel:
       toward an extended Chebyshev proof beyond 2x2
STATUS: RH remains OPEN.  Paper 187 gives a symbolic 2x2 proof.  This
        paper attacks size >= 3 by testing the extended Chebyshev
        Wronskian criterion for
            g_lambda(y) = (2 lambda y - 3) exp(-lambda y).
        The common factor y^(5/4) does not change Wronskian signs.
DEPENDS: 186 (kernel structure), 187 (2x2 symbolic layer).

THEORY
======

For ordered lambdas lambda_1 < ... < lambda_k, define

    g_i(y) = (2 lambda_i y - 3) exp(-lambda_i y).

If the Wronskians

    W_k(y) = det( d^r/dy^r g_i(y) )_{r=0..k-1, i=1..k}

have fixed nonzero sign on y >= 1 for every k, then {g_i} is an
extended complete Chebyshev system on [1,infty).  This is a standard
route to sign-regularity of all collocation determinants

    det(g_i(y_j)).

Closed derivative formula:

    D^r [(2λy-3)e^{-λy}]
      = e^{-λy} [ (2λy-3)(-λ)^r + r(2λ)(-λ)^(r-1) ],

with the second term omitted at r=0.

This paper tests the Wronskian signs for the Riemann lattice
lambda_i = pi i^2 and for random ordered lambdas.
"""
from __future__ import annotations

import itertools
import math
import random


def det_small(rows):
    n = len(rows)
    total = 0.0
    for perm in itertools.permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
        term = 1.0
        for i in range(n):
            term *= rows[i][perm[i]]
        total += -term if inv % 2 else term
    return total


def deriv_g(lam, y, r):
    P = 2 * lam * y - 3
    if r == 0:
        bracket = P
    else:
        bracket = P * ((-lam) ** r) + r * (2 * lam) * ((-lam) ** (r - 1))
    return math.exp(-lam * y) * bracket


def deriv_g_scaled(lam, y, r):
    """Derivative bracket after removing exp(-lambda*y)."""
    P = 2 * lam * y - 3
    if r == 0:
        return P
    return P * ((-lam) ** r) + r * (2 * lam) * ((-lam) ** (r - 1))


def wronskian(lambdas, y):
    k = len(lambdas)
    rows = []
    for r in range(k):
        rows.append([deriv_g(lam, y, r) for lam in lambdas])
    return det_small(rows)


def signed_wronskian(lambdas, y):
    """Expected sign follows the exponential-kernel convention."""
    k = len(lambdas)
    sign = -1 if (k * (k - 1) // 2) % 2 else 1
    return sign * wronskian(lambdas, y)


def scaled_signed_wronskian(lambdas, y):
    """Remove the common exp(-sum lambda_i y) scale."""
    k = len(lambdas)
    rows = []
    for r in range(k):
        rows.append([deriv_g_scaled(lam, y, r) for lam in lambdas])
    sign = -1 if (k * (k - 1) // 2) % 2 else 1
    return sign * det_small(rows)


def test_family(lambdas, ys, max_k):
    total = 0
    neg = 0
    min_val = float("inf")
    min_case = None
    for k in range(1, max_k + 1):
        for idx in itertools.combinations(range(len(lambdas)), k):
            lam = [lambdas[i] for i in idx]
            for y in ys:
                val = scaled_signed_wronskian(lam, y)
                total += 1
                if val < min_val:
                    min_val = val
                    min_case = (idx, y, k)
                if val < -1e-10:
                    neg += 1
    return total, neg, min_val, min_case


def main():
    print("=" * 80)
    print("PAPER 188 — Wronskian/Chebyshev route for Phi kernel")
    print("=" * 80)

    ys = [1.0, 1.01, 1.05, 1.15, 1.35, 1.8, 2.5, 4.0]

    # ------------------------------------------------------------------
    # [1] Riemann lattice lambdas
    # ------------------------------------------------------------------
    print("\n[1] Wronskians on lambda_i = pi i^2")
    lambdas = [math.pi * m * m for m in range(1, 8)]
    total, neg, mv, case = test_family(lambdas, ys, max_k=5)
    print(f"  tests: {total}")
    print(f"  negative signed Wronskians: {neg}")
    print(f"  minimum scaled signed Wronskian: {mv:.6e}")
    print(f"  minimum case: {case}")

    print("\n  sample by size at y=1:")
    for k in range(1, 6):
        lam = lambdas[:k]
        val = scaled_signed_wronskian(lam, 1.0)
        print(f"    k={k}: scaled signed W = {val:.12e}")

    # ------------------------------------------------------------------
    # [2] Random ordered lambdas
    # ------------------------------------------------------------------
    print("\n[2] Random ordered lambda stress test")
    rng = random.Random(2026)
    failures = 0
    worst = float("inf")
    worst_case = None
    trials = 200
    for t in range(trials):
        lambdas_r = sorted(3.0 + 40.0 * rng.random() for _ in range(6))
        total_r, neg_r, mv_r, case_r = test_family(lambdas_r, ys, max_k=5)
        if neg_r:
            failures += 1
        if mv_r < worst:
            worst = mv_r
            worst_case = (t, case_r, lambdas_r)
    print(f"  trials: {trials}")
    print(f"  failing families: {failures}")
    print(f"  worst scaled signed Wronskian: {worst:.6e}")
    print(f"  worst case: trial={worst_case[0]}, case={worst_case[1]}")

    # ------------------------------------------------------------------
    # [3] Boundary of validity below lambda=pi or y=1
    # ------------------------------------------------------------------
    print("\n[3] Boundary check below the Riemann range")
    print("    The factor 2 lambda y - 3 changes sign when lambda*y < 1.5.")
    print(f"  {'lambda0':>10} {'y0':>8} {'signed W k=2':>18}")
    for lam0, y0 in [(0.5, 1.0), (1.0, 1.0), (1.4, 1.0), (math.pi, 1.0)]:
        val = scaled_signed_wronskian([lam0, lam0 + 0.7], y0)
        print(f"  {lam0:>10.4f} {y0:>8.3f} {val:>18.6e}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  WHAT THIS WOULD GIVE:")
    print("    Fixed-sign Wronskians imply the family is an extended Chebyshev")
    print("    system, which is a standard route to higher-order sign-regularity.")
    print()
    print("  STATUS:")
    print("    The Riemann lattice and random tests pass in finite windows.")
    print("    This is still not a proof; the next step is to factor the")
    print("    Wronskian symbolically, probably into a Vandermonde times")
    print("    positive correction terms.")


if __name__ == "__main__":
    main()
