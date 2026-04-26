"""
PAPER: 040 (canonical: 30pandrosion_zeta.pdf)
TITLE: Pandrosion-Zeta Function
STATUS: framework (zeta_P(s) = sum 1/P'(alpha_k)^s)
DEPENDS: 003, 011

THEORY
======

PANDROSION-ZETA FUNCTION:
  zeta_P(s) := sum_{k=1}^d 1 / P'(alpha_k)^s.

Special values:
  zeta_P(1) = sum 1/P'(alpha_k) = 0  (Lagrange-Sylvester, paper 3).
  zeta_P(2) = sum 1/P'(alpha_k)^2 = -2/disc(P)  (for monic Z-polys, related).

VERIFICATION
============

  1. zeta_P(1) = 0 (Lagrange-Sylvester identity).
  2. zeta_P(2) and connection to discriminant.
  3. Smyth's L_0: compute zeta(s) for s = 1, 2, 3.
"""
from __future__ import annotations
import numpy as np


def zeta_P(P, s):
    roots = np.roots(P)
    Pp = np.polyder(P)
    return sum(1.0 / np.polyval(Pp, ak)**s for ak in roots)


def main():
    print("=" * 80)
    print("PAPER 40 — Pandrosion-zeta function")
    print("=" * 80)

    test_polys = [
        ("z^3 - z - 1 (plastic)", np.array([1.0, 0, -1, -1])),
        ("z^4 + z + 1", np.array([1.0, 0, 0, 1, 1])),
        ("Smyth L_0", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])),
    ]

    print("\n[1] zeta_P(1) = sum 1/P'(alpha) = 0 (Lagrange-Sylvester)")
    for name, P in test_polys:
        z1 = zeta_P(P, 1)
        print(f"  {name}: zeta_P(1) = {z1:.4e}")

    print("\n[2] zeta_P(2) and discriminant")
    for name, P in test_polys:
        z2 = zeta_P(P, 2)
        # disc via product
        roots = np.roots(P)
        d = len(roots)
        disc = float(np.prod([(roots[i] - roots[j])**2
                              for i in range(d) for j in range(i+1, d)]).real)
        print(f"  {name}: zeta_P(2) = {z2:.4f}, disc = {disc:.4e}")

    print("\n[3] zeta_P(s) for s = 1, 2, 3, 4 on Smyth L_0")
    P = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])
    for s in [1, 2, 3, 4]:
        z = zeta_P(P, s)
        print(f"  zeta_P({s}) = {z}")


if __name__ == "__main__":
    main()
