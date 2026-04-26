"""
PAPER: 190 (NEW — lattice Q root barrier)
TITLE: Riemann-lattice positivity of Q_k as a univariate root-barrier
       problem in y
STATUS: RH remains OPEN.  Paper 189 reduced the Wronskian sign problem to
        positivity of Q_k(2*pi*m_i^2*y).  This paper turns each fixed
        integer set m_1<...<m_k into a univariate polynomial in y and
        checks whether all positive roots lie below y=1.
DEPENDS: 189 (Wronskian Vandermonde factor).

THEORY
======

For fixed integers M={m_1<...<m_k}, define

    q_M(y) := Q_k(2*pi*m_1^2*y, ..., 2*pi*m_k^2*y).

This is a degree-k polynomial:

    q_M(y) = sum_{r=0}^k (-1)^(k-r) (2(k-r)+1)!!
             e_r(2*pi*m_1^2,...,2*pi*m_k^2) y^r.

The Wronskian program needs

    q_M(y) > 0  for all y >= 1.

A concrete finite diagnostic for a given M:
  - compute all roots of q_M;
  - verify q_M(1)>0;
  - verify every positive real root is <1.

This does not prove the all-k theorem, but it converts the remaining
analytic target into a root-location problem for explicit polynomials.
"""
from __future__ import annotations

import itertools
import math


def odd_double_factorial(n):
    out = 1.0
    for m in range(n, 0, -2):
        out *= m
    return out


def elementary_symmetric(xs, r):
    if r == 0:
        return 1.0
    return sum(math.prod(c) for c in itertools.combinations(xs, r))


def q_coeffs_for_ms(ms):
    """Return coefficients low-to-high for q_M(y)."""
    k = len(ms)
    base = [2 * math.pi * m * m for m in ms]
    coeffs = []
    for r in range(k + 1):
        coeffs.append(((-1) ** (k - r))
                      * odd_double_factorial(2 * (k - r) + 1)
                      * elementary_symmetric(base, r))
    return coeffs


def poly_eval(coeffs, y):
    total = 0.0
    for c in reversed(coeffs):
        total = total * y + c
    return total


def main():
    print("=" * 80)
    print("PAPER 190 — lattice Q root barrier")
    print("=" * 80)

    try:
        import numpy as np
    except ImportError:
        print("\n  [numpy required]")
        return

    # ------------------------------------------------------------------
    # [1] Consecutive initial blocks
    # ------------------------------------------------------------------
    print("\n[1] Consecutive initial blocks M={1,...,k}")
    print(f"  {'k':>3} {'q_M(1)':>18} {'largest pos root':>18}"
          f" {'# pos roots >=1':>16}")
    for k in range(2, 21):
        ms = tuple(range(1, k + 1))
        coeffs = q_coeffs_for_ms(ms)
        roots = np.roots(list(reversed(coeffs)))
        pos_roots = sorted(r.real for r in roots if abs(r.imag) < 1e-7 and r.real > 0)
        largest = max(pos_roots) if pos_roots else float("nan")
        n_ge1 = sum(1 for r in pos_roots if r >= 1 - 1e-8)
        print(f"  {k:>3} {poly_eval(coeffs, 1.0):>18.6e}"
              f" {largest:>18.8f} {n_ge1:>16}")

    # ------------------------------------------------------------------
    # [2] Exhaustive subsets up to bounded M
    # ------------------------------------------------------------------
    print("\n[2] Exhaustive subset root barrier")
    print(f"  {'M_max':>5} {'k':>3} {'subsets':>9} {'bad':>6}"
          f" {'worst root':>14} {'worst set':>18}")
    for M_max in [8, 10, 12]:
        for k in range(2, min(7, M_max) + 1):
            bad = 0
            total = 0
            worst_root = -float("inf")
            worst_set = None
            for ms in itertools.combinations(range(1, M_max + 1), k):
                coeffs = q_coeffs_for_ms(ms)
                roots = np.roots(list(reversed(coeffs)))
                pos_roots = [r.real for r in roots if abs(r.imag) < 1e-7 and r.real > 0]
                local = max(pos_roots) if pos_roots else -float("inf")
                if local > worst_root:
                    worst_root = local
                    worst_set = ms
                if poly_eval(coeffs, 1.0) <= 0 or any(r >= 1 - 1e-8 for r in pos_roots):
                    bad += 1
                total += 1
            print(f"  {M_max:>5} {k:>3} {total:>9} {bad:>6}"
                  f" {worst_root:>14.8f} {str(worst_set):>18}")

    # ------------------------------------------------------------------
    # [3] Root pattern for the hardest small set
    # ------------------------------------------------------------------
    print("\n[3] Root pattern for initial block k=10")
    ms = tuple(range(1, 11))
    coeffs = q_coeffs_for_ms(ms)
    roots = np.roots(list(reversed(coeffs)))
    roots_sorted = sorted(roots, key=lambda z: (z.real, z.imag))
    for r in roots_sorted:
        print(f"  root: {r.real: .10f} {r.imag:+.3e}i")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  PROGRESS:")
    print("    For each fixed lattice subset M, Wronskian positivity on y>=1")
    print("    reduces to a concrete root-location certificate for q_M(y).")
    print()
    print("  WHAT THE TESTS INDICATE:")
    print("    The consecutive initial blocks and bounded subsets tested have")
    print("    no positive roots >= 1.")
    print()
    print("  REMAINING THEOREM:")
    print("    Prove a uniform root barrier: every positive real root of q_M")
    print("    lies below 1 for all finite M subset positive integers.")


if __name__ == "__main__":
    main()
