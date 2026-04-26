"""
PAPER: 072 (canonical: 62pandrosion_ek.pdf)
TITLE: Edelman-Kostlan (Part 2): Eneström-Kakeya / unit-disk
STATUS: proved (Eneström 1893, Kakeya 1912)
DEPENDS: 026

THEORY
======

ENESTRÖM-KAKEYA: For polynomial P(z) = a_0 + a_1 z + ... + a_d z^d with
0 < a_0 <= a_1 <= ... <= a_d (real positive non-decreasing coeffs),
all roots satisfy |alpha| <= 1.
Conversely, decreasing coefficients give |alpha| >= 1.

PANDROSION CONNECTION: Roots in unit disk iff Re(z F_P(z)) > 0 on |z| = 1
(paper 67). Eneström-Kakeya gives a sufficient coefficient condition.

VERIFICATION
============

  1. Increasing coefficients: roots in unit disk.
  2. Eneström-Kakeya bounds tight at z = -1 and z = 1.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 72 — Eneström-Kakeya")
    print("=" * 80)

    print("\n[1] Increasing coeffs: max|alpha| <= 1")
    test_polys = [
        ([1, 2, 3], "1 + 2z + 3z^2"),
        ([1, 2, 4, 5, 7], "1 + 2z + 4z^2 + 5z^3 + 7z^4"),
        ([1, 1.5, 2.0, 2.5], "1 + 1.5 z + 2 z^2 + 2.5 z^3"),
    ]
    for coefs, name in test_polys:
        # numpy convention: high-to-low. Reverse.
        P = list(reversed(coefs))
        roots = np.roots(P)
        max_r = max(abs(r) for r in roots)
        print(f"  {name}: max|root| = {max_r:.4f}  ({'OK' if max_r <= 1 + 1e-9 else 'EXCEEDS'})")

    print("\n[2] Decreasing coefficients: min|alpha| >= 1")
    for coefs, name in [
        ([3, 2, 1], "3 + 2z + z^2"),
        ([7, 5, 4, 2, 1], "7 + 5z + 4z^2 + 2z^3 + z^4"),
    ]:
        P = list(reversed(coefs))
        roots = np.roots(P)
        min_r = min(abs(r) for r in roots)
        print(f"  {name}: min|root| = {min_r:.4f}  ({'OK' if min_r >= 1 - 1e-9 else 'BELOW'})")

    print("\n[3] Tight case: cyclotomic-like with constant coefficients")
    P = list(reversed([1, 1, 1, 1]))  # 1 + z + z^2 + z^3
    roots = np.roots(P)
    print(f"  1 + z + z^2 + z^3: roots = {roots}, all |.| = {[f'{abs(r):.4f}' for r in roots]}")


if __name__ == "__main__":
    main()
