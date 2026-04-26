"""
PAPER: 099 (canonical: 89pandrosion_perron.pdf)
TITLE: Perron-Frobenius via Pandrosion Companion
STATUS: proved (Perron 1907, Frobenius 1912)

THEORY
======

PERRON-FROBENIUS: A non-negative matrix has a non-neg real dominant
eigenvalue (Perron root), with non-neg eigenvector.

PANDROSION CONNECTION: For companion matrix C_P (whose char poly is P with
non-neg coefficients), the Perron root is the rightmost pole of F_P =
P'/P on R^+, detected by sign change.

VERIFICATION
============

  1. Perron root of non-neg matrix.
  2. Sign change of F_P at Perron root.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 99 — Perron-Frobenius via Pandrosion companion")
    print("=" * 80)

    print("\n[1] Non-negative matrix: dominant eigenvalue is real, positive")
    A = np.array([[2, 1, 1], [1, 3, 1], [1, 1, 4.0]])
    eigs = np.linalg.eigvals(A)
    print(f"  A = non-neg 3x3:")
    print(f"  Eigenvalues: {sorted(eigs, key=lambda e: e.real, reverse=True)}")
    perron = max(eigs, key=lambda e: e.real)
    print(f"  Perron root: {perron.real:.4f}")

    print("\n[2] Companion poly of A: F_P sign change at Perron root")
    P = np.poly(eigs)
    roots = sorted(np.roots(P).real, reverse=True)
    print(f"  Char poly P, real roots: {roots}")
    Pp = np.polyder(P)
    perron_real = roots[0]
    eps = 0.01
    F_left = np.polyval(Pp, perron_real - eps) / np.polyval(P, perron_real - eps)
    F_right = np.polyval(Pp, perron_real + eps) / np.polyval(P, perron_real + eps)
    print(f"  F_P({perron_real - eps}) = {F_left:.4f}")
    print(f"  F_P({perron_real + eps}) = {F_right:.4f}")
    print(f"  Sign change at Perron root: {F_left * F_right < 0}")


if __name__ == "__main__":
    main()
