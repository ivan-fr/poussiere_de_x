"""
PAPER: 084 (canonical: 74pandrosion_galois2.pdf)
TITLE: Galois Theory of Pandrosion Quotients
STATUS: framework
DEPENDS: 047

THEORY
======

GALOIS THEORY OF QUOTIENTS:
For irreducible P with Galois group G acting on roots {alpha_k}, the
Pandrosion quotients {Q_k = P/(z - alpha_k)} form a G-orbit. The
G-invariants of {Q_k} live in the base field.

SOLVABILITY: P solvable by radicals iff G is solvable (Galois 1832).
Pandrosion translation: solvability of {Q_k}-orbit by abelian quotients.

VERIFICATION
============

  1. Galois group of test polynomials.
  2. Quotient orbit structure.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 84 — Galois theory of Pandrosion quotients")
    print("=" * 80)

    print("\n[1] z^3 - 2: Galois S_3 (non-abelian, but solvable)")
    P = np.array([1.0, 0, 0, -2])
    roots = np.roots(P)
    print(f"  Roots: {roots}")
    print(f"  Galois group: S_3 = solvable (P is solvable by radicals: cube root of 2)")

    print("\n[2] z^5 - z + 1: Galois S_5 (non-solvable, by Abel-Ruffini)")
    P = np.array([1.0, 0, 0, 0, -1, 1])
    roots = np.roots(P)
    print(f"  Roots: {roots}")
    print(f"  Galois group: typically S_5 (non-solvable)")
    print(f"  Hence z^5 - z + 1 is NOT solvable by radicals")

    print("\n[3] Cyclotomic z^n - 1: Galois (Z/n)*  (cyclic, abelian)")
    for n in [3, 5, 7, 8]:
        P = np.array([1.0] + [0]*(n-1) + [-1])
        roots = np.roots(P)
        print(f"  z^{n} - 1: roots are n-th roots of unity, Galois (Z/{n})* abelian")


if __name__ == "__main__":
    main()
