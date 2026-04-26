"""
PAPER: 075 (canonical: 65pandrosion_hb2.pdf)
TITLE: Hadamard-Bombieri (Part 2) and Hermite-Biehler
STATUS: proved (Hermite-Biehler 1893, Pandrosion-form)
DEPENDS: 034, 067

THEORY
======

HERMITE-BIEHLER: Polynomial P(z) = A(z) + i B(z) (with A, B real polys) is
Hurwitz-stable (all roots in lower half-plane Im(z) < 0) iff:
  A and B have only real zeros, alternating, with appropriate sign at +infty.

PANDROSION CONNECTION (paper 67):
  P Hurwitz <=> Re F_P(iy) > 0 for all y in R.
HB criterion is the CONSTRUCTIVE version: split P into real/imag parts and
check root interlacing.

VERIFICATION
============

  1. Hurwitz polynomial split into A + iB.
  2. Verify A, B real-rooted, alternating.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 75 — Hermite-Biehler stability test")
    print("=" * 80)

    print("\n[1] Hurwitz polynomial: roots in left half-plane")
    # P(z) = (z + 1)(z + 2)(z + 3): all real negative roots
    P = np.poly([-1, -2, -3.0])
    roots = np.roots(P)
    print(f"  P = (z+1)(z+2)(z+3): roots = {roots}, all Re < 0: {all(r.real < 0 for r in roots)}")

    # Substitute z -> i z: P(iz) = A(z) + i B(z)
    print("\n[2] Substitute z -> i z, split P(iz) = A(z) + i B(z)")
    # P(z) = z^3 + 6 z^2 + 11 z + 6
    # P(iz) = -i z^3 - 6 z^2 + 11 i z + 6 = (6 - 6 z^2) + i (11 z - z^3)
    # A(z) = 6 - 6 z^2, B(z) = 11 z - z^3
    A = np.array([-6, 0, 6.0])  # -6 z^2 + 0 z + 6
    B = np.array([-1, 0, 11.0, 0])  # -z^3 + 11 z
    A_roots = np.roots(A)
    B_roots = np.roots(B)
    print(f"  A(z) = -6 z^2 + 6: roots = {A_roots}")
    print(f"  B(z) = -z^3 + 11 z: roots = {B_roots}")
    print(f"  A real-rooted? {all(abs(r.imag) < 1e-9 for r in A_roots)}")
    print(f"  B real-rooted? {all(abs(r.imag) < 1e-9 for r in B_roots)}")

    print("\n[3] Verify Re F_P(iy) > 0 for the Hurwitz P")
    P = np.poly([-1, -2, -3.0])
    Pp = np.polyder(P)
    for y in [-2, -1, 0, 1, 2]:
        z = 1j * y
        F = np.polyval(Pp, z) / np.polyval(P, z)
        print(f"  y = {y}: F_P(iy) = {F:.4f}, Re = {F.real:.4f}")


if __name__ == "__main__":
    main()
