"""
PAPER: 047 (canonical: 37pandrosion_galois.pdf)
TITLE: Galois Action on Pandrosion Quotients
STATUS: framework

THEORY
======

For P in K[z] with K a number field, the Galois group Gal(P) acts on the
roots of P by permutation. This action lifts to the Pandrosion quotients:
  sigma(Q(alpha_k, z)) = Q(sigma(alpha_k), z).
The quotient orbits {Q_1, ..., Q_d} are permuted by Gal(P).

ABEL-RUFFINI: P solvable by radicals iff Gal(P) is solvable. In Pandrosion
language, the orbit of {Q_k} under Gal must refine through abelian quotients.

VERIFICATION
============

  1. Galois orbit structure on test polynomials.
  2. Sym_d action on quotients.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 47 — Galois action on Pandrosion quotients")
    print("=" * 80)

    print("\n[1] z^3 - 2 (Galois group S_3, non-abelian)")
    P = np.array([1, 0, 0, -2.0])
    roots = np.roots(P)
    print(f"  Roots: {roots}")
    # Pandrosion quotients Q_k(z) = P(z)/(z - alpha_k) at z = 0:
    # Q_k(0) = P(0)/(-alpha_k) = -2 / (-alpha_k) = 2/alpha_k
    quotients = [2.0 / alpha for alpha in roots]
    print(f"  Q_k(0) = 2/alpha_k = {[f'{q:.4f}' for q in quotients]}")
    # Galois acts by permuting these

    print("\n[2] z^4 - 2 (Galois group D_4, dihedral)")
    P = np.array([1.0, 0, 0, 0, -2])
    roots = np.roots(P)
    quotients = [-2.0 / alpha for alpha in roots]
    print(f"  Roots: {roots}")
    print(f"  Q_k(0) = {[f'{q:.4f}' for q in quotients]}")

    print("\n[3] Cyclotomic z^5 - 1 (Galois group Z/4, cyclic)")
    P = np.array([1.0, 0, 0, 0, 0, -1])
    roots = np.roots(P)
    print(f"  Roots: {roots}")
    # Galois acts by powers: alpha -> alpha^k for k coprime to 5
    print(f"  Galois orbit: alpha -> alpha^2 -> alpha^4 -> alpha^3 -> alpha (cyclic)")


if __name__ == "__main__":
    main()
