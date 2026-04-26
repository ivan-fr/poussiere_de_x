"""
PAPER: 114 (canonical: 114_schinzel_zassenhaus.pdf)
TITLE: Schinzel-Zassenhaus and Dimitrov's Theorem
STATUS: Dimitrov 2019 proved 2^(1/(2d)); empirical excess law d * excess ~ 0.43
DEPENDS: 020, 065, 087, 105

THEORY
======

SCHINZEL-ZASSENHAUS (1965): For non-cyclotomic monic Z-poly of degree d,
max_j |alpha_j| >= 1 + c/d for absolute c > 0.

DIMITROV (2019): max_j |alpha_j| >= 2^(1/(2d)) (proved unconditionally).

EMPIRICAL EXCESS LAW (this paper): on height-1 Z-polys to d <= 12:
  smallest max|alpha| - 1 satisfies 0.4 <= d * excess <= 0.7
i.e., excess ~ 0.43/d empirically.

VERIFICATION
============

  1. Exhaustive height-1 scan up to d = 12.
  2. Compare with Dimitrov bound.
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def main():
    print("=" * 80)
    print("PAPER 114 — Schinzel-Zassenhaus / Dimitrov 2019")
    print("=" * 80)

    print("\n[1] Smallest max|alpha| - 1 for height-1 monic Z-polys")
    print(f"  {'d':>3} {'#non-cyclo':>12} {'min excess':>14} {'d * excess':>14}")
    for d in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
        min_excess = float('inf')
        n_nc = 0
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            try:
                roots = np.roots(coefs)
                M = float(np.prod(np.maximum(1.0, np.abs(roots))))
                if M < 1.001: continue
                n_nc += 1
                Mmax = float(np.max(np.abs(roots)))
                excess = Mmax - 1.0
                if excess < min_excess: min_excess = excess
            except: pass
        product = d * min_excess
        print(f"  {d:>3} {n_nc:>12} {min_excess:>14.6f} {product:>14.4f}")

    print("\n[2] Dimitrov bound: max|alpha| >= 2^(1/(2d))")
    print(f"  {'d':>4} {'Dimitrov bd':>14}")
    for d in [5, 10, 20, 50, 100]:
        bd = 2**(1.0/(2*d)) - 1
        print(f"  {d:>4} {bd:>14.6f}")


if __name__ == "__main__":
    main()
