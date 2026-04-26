"""
PAPER: 089 (canonical: 79pandrosion_kadison.pdf)
TITLE: Kadison-Singer, Interlacing Families, and Pandrosion Barrier
STATUS: proved (MSS 2015)
DEPENDS: 011

THEORY
======

KADISON-SINGER (1959, proved by Marcus-Spielman-Srivastava 2015):
Every state on a maximal abelian subalgebra of B(H) extends uniquely.

MSS BARRIER FUNCTION: Phi(x) = sum 1/(x - alpha_j)^2 = -F_P'(x).
This is the Pandrosion-CURVATURE (paper 86).

INTERLACING FAMILIES: A method to bound max root of expected characteristic
polynomial.

PANDROSION VIEW: F_P' is the gradient of the energy / barrier, controlling
where the max root sits.

VERIFICATION
============

  1. MSS barrier Phi(x) = -F_P'(x) = sum 1/(x - alpha_j)^2.
  2. Phi positive on R for real-rooted P.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 89 — MSS / Kadison-Singer: Pandrosion barrier")
    print("=" * 80)

    print("\n[1] Phi(x) = sum 1/(x - alpha_j)^2 = -F_P'(x)")
    P = np.poly([1.0, 2.0, 3.0, 4.0])
    roots = sorted(np.roots(P).real)
    Pp = np.polyder(P)
    Ppp = np.polyder(Pp)
    for x in [0.0, 1.5, 2.5, 5.0]:
        # Phi via direct sum
        Phi_direct = sum(1.0 / (x - r)**2 for r in roots)
        # Phi via -F_P': F_P = P'/P, F_P' = (P''P - (P')^2)/P^2 = -E/P^2
        E = np.polyval(Pp, x)**2 - np.polyval(P, x) * np.polyval(Ppp, x)
        Phi_via_E = E / np.polyval(P, x)**2
        print(f"  x = {x}: Phi direct = {Phi_direct:.4f}, Phi = E/P^2 = {Phi_via_E:.4f}")

    print("\n[2] MSS interlacing families: max-root bound")
    print("  Idea: family of polys p_S(x) interlace; max root of expected")
    print("  poly bounded by max root of any single one.")
    print("\n  In Pandrosion terms: barrier Phi(x) controls the rightmost zero.")


if __name__ == "__main__":
    main()
