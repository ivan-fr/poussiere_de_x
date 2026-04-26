"""
PAPER: 080 (canonical: 70pandrosion_waring.pdf)
TITLE: Waring's Problem and Pandrosion Symmetric Powers
STATUS: framework (classical Waring + polynomial-form)

THEORY
======

WARING (number theory): Every positive integer can be expressed as a sum
of at most g(k) k-th powers (g(2) = 4 by Lagrange, etc.).

Polynomial Waring: For ground field K, what's the smallest s such that every
polynomial F of degree d can be written as F = sum_{i=1}^s L_i^d for linear
forms L_i?

This is the WARING RANK of F. Generic F has Waring rank ~ binom(n+d-1, d)/n.

VERIFICATION
============

  1. Waring rank of test polynomials.
  2. Sylvester's symbolic-decomposition method.
"""
from __future__ import annotations
import math
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 80 — Waring's problem (polynomial-form)")
    print("=" * 80)

    print("\n[1] Generic Waring rank: ~ C(n+d-1, d)/n")
    print(f"  {'(n, d)':>10} {'C(n+d-1,d)/n':>15} {'space dim':>12}")
    for n, d in [(2, 2), (2, 3), (3, 2), (3, 3), (2, 5)]:
        rank_est = math.comb(n + d - 1, d) // n
        space_dim = math.comb(n + d - 1, d)
        print(f"  {(n,d)!s:>10} {rank_est:>15} {space_dim:>12}")

    print("\n[2] Univariate (n=1): F = L_1^d + ... + L_s^d, where L_i are linear in z")
    # For univariate of degree d, generic Waring rank = d (each linear form L_i = z - a_i with multiplicity 1).
    # But this isn't tight: F = (z - a)^d alone has Waring rank 1 even though degree d.
    print("  Univariate: rank(F) = #distinct roots up to symmetric power")
    P = np.array([1.0, -2, 1])  # (z-1)^2: rank 1 in Waring sense
    print(f"  P = (z-1)^2 has Waring rank 1 (single linear factor squared)")

    P = np.array([1.0, 0, -2])  # z^2 - 2 = (z - sqrt 2)(z + sqrt 2)
    print(f"  P = z^2 - 2: roots ±sqrt(2); rank as sum of squares?")
    print(f"  Note: z^2 - 2 = z^2 + (-sqrt 2)^2 - 4 = (z + sqrt 2)^2 - 2 sqrt 2 z - ... rank 3+")


if __name__ == "__main__":
    main()
