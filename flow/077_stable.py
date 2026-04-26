"""
PAPER: 077 (canonical: 67pandrosion_stable.pdf)
TITLE: Stable and Hyperbolic Polynomials: Unified Pandrosion Classification
STATUS: proved (classical, unified Pandrosion-form)
DEPENDS: 029, 067

THEORY
======

UNIFIED CLASSIFICATION via Pandrosion field F_P = P'/P:

  Stability class    Root region        Pandrosion criterion
  ----------------- ----------------    ------------------------------------
  Hyperbolic         alpha_j in R       (P')^2 - P P'' >= 0 on R (Turán)
  Hurwitz            Re(alpha_j) < 0    Re F_P(iy) > 0 for all y in R
  Upper HP           Im(alpha_j) > 0    Im F_P(x) > 0 for all x in R
  Schur              |alpha_j| < 1      Re(z F_P(z)) > 0 on |z| = 1

VERIFICATION
============

  1. Test each stability class.
  2. Unified Pandrosion field criterion.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 77 — Unified stability via Pandrosion field")
    print("=" * 80)

    print("\n[1] Hurwitz: roots in left half-plane, Re F_P(iy) > 0")
    P_hurwitz = np.poly([-1, -2, -3.0])
    Pp = np.polyder(P_hurwitz)
    for y in [-2, 0, 2]:
        F = np.polyval(Pp, 1j * y) / np.polyval(P_hurwitz, 1j * y)
        print(f"  P=(z+1)(z+2)(z+3), y={y}: Re F = {F.real:+.4f}")

    print("\n[2] Upper HP: roots Im > 0, Im F_P(x) > 0 on R")
    P_uhp = np.poly([1+1j, 2+0.5j, 0.5+2j])
    Pp = np.polyder(P_uhp)
    for x in [-2.0, 0.0, 1.5]:
        F = np.polyval(Pp, x) / np.polyval(P_uhp, x)
        print(f"  P with Im(roots)>0, x={x}: Im F = {F.imag:+.4f}")

    print("\n[3] Schur: |roots| < 1, Re(z F_P) > 0 on |z| = 1")
    P_schur = np.poly([0.5, -0.3, 0.2 + 0.4j])
    Pp = np.polyder(P_schur)
    for theta in [0, np.pi/3, np.pi/2, np.pi]:
        z = np.exp(1j * theta)
        zF = z * np.polyval(Pp, z) / np.polyval(P_schur, z)
        print(f"  theta = {theta:.3f}: Re(z F) = {zF.real:+.4f}")

    print("\n[4] Hyperbolic: real roots, T = (P')^2 - P P'' >= 0 on R")
    P_hyp = np.poly([1.0, 2.0, 3.0])
    Pp = np.polyder(P_hyp)
    Ppp = np.polyder(Pp)
    for x in [0.0, 1.5, 2.5, 4.0]:
        T = np.polyval(Pp, x)**2 - np.polyval(P_hyp, x) * np.polyval(Ppp, x)
        print(f"  P=(z-1)(z-2)(z-3), x={x}: T = {T:+.4f}")


if __name__ == "__main__":
    main()
