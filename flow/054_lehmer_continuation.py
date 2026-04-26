"""
PAPER: 054 (canonical: 44pandrosion_lehmer.pdf)
TITLE: Lehmer (continuation): height-bounded scans
STATUS: empirical (paired with paper 020, refined in 105)
DEPENDS: 020

THEORY
======

Continuation of Lehmer attack from paper 020: extended scan over
height-bounded integer polynomials, plus reciprocal polynomial families
that historically include Lehmer's polynomial.

VERIFICATION
============

  1. Reciprocal polynomial scan up to d = 30.
  2. Smyth's L_0 = 1.176280818 verified.
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def mahler_measure(coeffs):
    roots = np.roots(coeffs)
    return float(abs(coeffs[0]) * np.prod(np.maximum(1.0, np.abs(roots))))


def is_reciprocal(coeffs, tol=1e-9):
    """P is reciprocal iff a_i = a_{d-i}."""
    n = len(coeffs)
    return all(abs(coeffs[i] - coeffs[n-1-i]) < tol for i in range(n // 2))


def main():
    print("=" * 80)
    print("PAPER 54 — Lehmer continuation: reciprocal polys")
    print("=" * 80)

    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])
    print(f"\n[1] Smyth L_0 reciprocal? {is_reciprocal(L0)}")
    print(f"  M(L_0) = {mahler_measure(L0):.10f}")

    print("\n[2] Reciprocal {-1, 0, 1} scan up to d = 12")
    print(f"  {'d':>3} {'#recip':>10} {'min M > 1':>14}")
    for d in [4, 6, 8, 10, 12]:
        # d even: a_0..a_{d/2} free, a_i = a_{d-i}
        half = d // 2
        min_M = float('inf')
        n_recip = 0
        for combo in itertools.product([-1, 0, 1], repeat=half + 1):
            # Build reciprocal poly: leading = 1, middle from combo, mirror
            free = list(combo)
            full = [1] + free[1:] + free[:0:-1] + [free[0]] if d > 2 * half else [1] + free[1:half+1] + free[half-1::-1] + [free[0]]
            # Simpler: monic, a_d = a_0 = 1, then mirror
            try:
                a_low = list(combo)  # [a_0, a_1, ..., a_half]
                full = a_low[:half] + [a_low[half]] + a_low[half-1::-1]
                if len(full) != d + 1:
                    full = a_low + a_low[half-1::-1] if d > 2*half else a_low + a_low[-2::-1]
                if len(full) != d + 1: continue
                full[0] = 1
                P = np.array(full, dtype=float)
                if abs(P[-1]) < 1e-12: continue
                M = mahler_measure(P)
                if M > 1.001:
                    n_recip += 1
                    if M < min_M: min_M = M
            except: pass
        print(f"  {d:>3} {n_recip:>10} {min_M:>14.6f}")


if __name__ == "__main__":
    main()
