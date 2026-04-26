"""
PAPER: 090 (canonical: 80pandrosion_kronecker.pdf)
TITLE: Kronecker's Theorem on Roots of Unity
STATUS: proved (Kronecker 1857)

THEORY
======

KRONECKER (1857): If P is a monic polynomial in Z[z] with all roots in
{|z| <= 1}, then either:
  - P has a root at 0, or
  - All roots of P are roots of unity.

EQUIVALENT (cyclotomic characterization): M(P) = 1 iff P is a product of
cyclotomic polynomials (and ± z^k factors).

PANDROSION CONNECTION: Kronecker boundary case for Lehmer's conjecture
(papers 020, 105). For Lehmer's L_0, M = 1.176 just above the cyclotomic
boundary.

VERIFICATION
============

  1. Mahler measure 1 = cyclotomic.
  2. Smyth L_0 just above 1 (1.176280).
"""
from __future__ import annotations
import numpy as np


def mahler_measure(P):
    return float(abs(P[0]) * np.prod(np.maximum(1.0, np.abs(np.roots(P)))))


def main():
    print("=" * 80)
    print("PAPER 90 — Kronecker's theorem")
    print("=" * 80)

    print("\n[1] Cyclotomic polys: M = 1")
    cyclos = [
        ("z - 1", np.array([1.0, -1])),
        ("z + 1", np.array([1.0, 1])),
        ("z^2 + 1", np.array([1.0, 0, 1])),
        ("z^2 + z + 1", np.array([1.0, 1, 1])),
        ("z^4 + z^3 + z^2 + z + 1", np.array([1.0, 1, 1, 1, 1])),
    ]
    for name, P in cyclos:
        M = mahler_measure(P)
        print(f"  {name}: M = {M:.6f}")

    print("\n[2] Non-cyclotomic with all roots in |z| <= 1: must have z = 0 root")
    P = np.array([1.0, 0, 0])  # z^2 (root at 0)
    M = mahler_measure(P)
    print(f"  P = z^2: M = {M}, root at 0 trivially")

    print("\n[3] Smyth L_0: just above cyclotomic boundary")
    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])
    M = mahler_measure(L0)
    print(f"  L_0: M = {M:.6f} (just above 1, Lehmer's conjecture: M >= 1.17628 universally)")


if __name__ == "__main__":
    main()
