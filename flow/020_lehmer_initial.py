"""
PAPER: 020 (canonical: 10pandrosion_lehmer.pdf)
TITLE: Lehmer's Conjecture: Initial Pandrosion Energy Approach
STATUS: empirical (refined in papers 044, 105, 112)
DEPENDS: 011

THEORY
======

LEHMER (1933): For monic non-cyclotomic integer polynomial P, the Mahler
measure M(P) = |a_d| prod max(1, |alpha_j|) satisfies M(P) >= L for some
universal L > 1.

PANDROSION-ENERGY ATTACK:
For P with roots alpha_j, the residue energy is
  E(P) = sum 1/|P'(alpha_j)|^2.
By the discriminant identity (paper 22): prod |P'(alpha_j)| = |Disc(P)|.
For non-cyclotomic Z-polys, |Disc(P)| >= 1 (integer).

PARTIAL: This bounds the harmonic energy but does NOT give M(P) > 1.

VERIFICATION
============

  1. Mahler measure of Smyth's L_0 = 1.17628.
  2. Smallest M > 1 in {-1,0,1} coeffs scan.
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def mahler_measure(coeffs):
    roots = np.roots(coeffs)
    return float(abs(coeffs[0]) * np.prod(np.maximum(1.0, np.abs(roots))))


def main():
    print("=" * 80)
    print("PAPER 20 — Lehmer initial Pandrosion energy")
    print("=" * 80)

    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], dtype=float)
    M_L0 = mahler_measure(L0)
    print(f"\n[1] Smyth's L_0 (Lehmer's polynomial): M(L_0) = {M_L0:.10f}")
    print(f"    Literature: 1.176280818...")

    print("\n[2] Plastic number z^3 - z - 1")
    plastic = np.array([1, 0, -1, -1.0])
    print(f"  M = {mahler_measure(plastic):.10f}  (plastic = 1.32472)")

    print("\n[3] Exhaustive scan height-1 coeffs: smallest M > 1")
    print(f"  {'d':>3} {'#total':>10} {'min M > 1':>14}")
    for d in [3, 4, 5, 6, 7, 8, 9, 10]:
        min_M = float('inf')
        n_total = 0
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            n_total += 1
            try:
                roots = np.roots(coefs)
                M = float(np.prod(np.maximum(1.0, np.abs(roots))))
                if M > 1.001 and M < min_M:
                    min_M = M
            except: pass
        print(f"  {d:>3} {n_total:>10} {min_M:>14.6f}")

    print("\n[4] Pandrosion energy: sum 1/P'(alpha_j)^2")
    for name, P in [("L_0", L0), ("plastic", plastic), ("z^3 - 2", np.array([1.0, 0, 0, -2]))]:
        roots = np.roots(P)
        Pp = np.polyder(P)
        E = sum(1.0 / abs(np.polyval(Pp, ak))**2 for ak in roots)
        print(f"  {name}: E = sum 1/|P'(a_k)|^2 = {E:.6f}")


if __name__ == "__main__":
    main()
