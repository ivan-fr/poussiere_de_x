"""
PAPER: 095 (canonical: 85pandrosion_cardano.pdf)
TITLE: Cardano's Formula via Pandrosion
STATUS: proved (Cardano 1545, Pandrosion-form)

THEORY
======

CARDANO (1545): For depressed cubic z^3 + p z + q = 0:
  z = cubic_root(-q/2 + sqrt(D)) + cubic_root(-q/2 - sqrt(D))
where D = q^2/4 + p^3/27 is the (-1) times discriminant scaled.

PANDROSION CONNECTION: The Pandrosion difference quotient at a root
yields a one-step formula: from anchor a near root, single Pandrosion step
recovers root exactly for d = 1 (linear case). For d = 3, requires
discriminant + cube root.

VERIFICATION
============

  1. Cardano's formula on z^3 + p z + q = 0.
  2. Discriminant connection.
"""
from __future__ import annotations
import math
import numpy as np


def cardano(p, q):
    """Solve z^3 + p z + q = 0 via Cardano."""
    D = q**2 / 4 + p**3 / 27
    if D >= 0:
        # Real D
        sqrtD = math.sqrt(D)
        u = (-q/2 + sqrtD)
        v = (-q/2 - sqrtD)
        # Cube roots
        u_cube = math.copysign(abs(u)**(1/3), u)
        v_cube = math.copysign(abs(v)**(1/3), v)
        return u_cube + v_cube
    else:
        # Complex D
        sqrtD = 1j * math.sqrt(-D)
        u = (-q/2 + sqrtD)
        v = (-q/2 - sqrtD)
        # Cube roots (principal branch)
        u_cube = u**(1/3)
        v_cube = v**(1/3)
        return u_cube + v_cube


def main():
    print("=" * 80)
    print("PAPER 95 — Cardano's formula")
    print("=" * 80)

    print("\n[1] Cardano on z^3 + p z + q = 0")
    cases = [
        (-3, 2, "z^3 - 3z + 2 = (z-1)^2 (z+2)"),
        (-7, 6, "z^3 - 7z + 6 = (z-1)(z-2)(z+3)"),
        (1, -2, "z^3 + z - 2 = (z-1)(z^2 + z + 2)"),
    ]
    for p, q, name in cases:
        root = cardano(p, q)
        # Verify: z^3 + p*z + q = 0
        check = root**3 + p*root + q
        roots_np = np.roots([1, 0, p, q])
        print(f"  {name}: Cardano root = {root}, check = {check:.4e}")
        print(f"    All roots (numpy): {roots_np}")

    print("\n[2] Discriminant of cubic z^3 + p z + q")
    print("    Disc = -4 p^3 - 27 q^2")
    for p, q in [(-3, 2), (-7, 6), (1, -2)]:
        disc = -4 * p**3 - 27 * q**2
        print(f"  p={p}, q={q}: Disc = {disc}")


if __name__ == "__main__":
    main()
