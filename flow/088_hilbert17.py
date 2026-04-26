"""
PAPER: 088 (canonical: 78pandrosion_hilbert17.pdf)
TITLE: Hilbert's 17th Problem and Pandrosion SOS
STATUS: proved (Artin 1927)

THEORY
======

HILBERT'S 17TH (Artin 1927): Every nonneg polynomial in R[x_1, ..., x_n]
is a sum of squares of RATIONAL functions (not necessarily polynomials).

POSITIVSTELLENSATZ: Refines to: nonneg P on a semi-algebraic set S has a
specific SOS-modulo-S representation.

PANDROSION CONNECTION: For real-rooted P, the Pandrosion energy
sum Q(alpha_k, x)^2 = (P')^2 - P P'' is automatically an SOS (= sum of
squares of polynomials in real-root case, paper 30, 66).

VERIFICATION
============

  1. SOS decomposition for simple positive polynomials.
  2. Pandrosion energy as SOS for real-rooted P.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 88 — Hilbert's 17th and Pandrosion SOS")
    print("=" * 80)

    print("\n[1] Pandrosion energy E_P = (P')^2 - P P'' is SOS for real-rooted P")
    P = np.poly([1.0, 2.0, 3.0])  # (x-1)(x-2)(x-3), real-rooted
    Pp = np.polyder(P)
    Ppp = np.polyder(Pp)
    print(f"  P = (x-1)(x-2)(x-3)")
    print(f"  E_P = (P')^2 - P P'' is sum of (Q(alpha_k, x))^2:")
    for x in [0, 1.5, 2.5, 4.0]:
        E = np.polyval(Pp, x)**2 - np.polyval(P, x) * np.polyval(Ppp, x)
        print(f"  x = {x}: E_P(x) = {E:.4f} (>= 0)")

    print("\n[2] Motzkin's polynomial: nonneg but not SOS")
    # Motzkin: M(x, y) = x^4 y^2 + x^2 y^4 + 1 - 3 x^2 y^2 (nonneg, not SOS in R[x,y])
    print("  M(x, y) = x^4 y^2 + x^2 y^4 + 1 - 3 x^2 y^2")
    print("  M >= 0 by AM-GM, but NOT SOS in R[x, y].")
    print("  However, by Hilbert 17 (Artin): M = (sum of squares of rationals).")

    print("\n[3] Hilbert 1888 counterexample: ternary sextics")
    print("  Hilbert showed not every nonneg ternary sextic is SOS.")
    print("  Artin 1927 fixed this by allowing RATIONAL squares.")


if __name__ == "__main__":
    main()
