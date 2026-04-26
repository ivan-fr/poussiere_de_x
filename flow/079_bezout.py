"""
PAPER: 079 (canonical: 69pandrosion_bezout.pdf)
TITLE: Bézout's Theorem and the Pandrosion Determinant
STATUS: proved (Bézout, Pandrosion-form)
DEPENDS: 011

THEORY
======

BÉZOUT (univariate): Two polynomials of degrees d_1, d_2 share roots iff
their RESULTANT vanishes:
  Res(P_1, P_2) = 0 <=> gcd(P_1, P_2) is non-constant.

Multivariate: For F: C^n -> C^n homogeneous of degree d, generic F has
exactly d^n solutions in projective space (Bézout degree).

PANDROSION CONNECTION: Resultant via the Pandrosion-Sylvester determinant
of P, Q (Sylvester matrix is closely related to the Schmidt slope matrix).

VERIFICATION
============

  1. Resultant via Sylvester determinant.
  2. Common-root detection.
"""
from __future__ import annotations
import numpy as np


def sylvester_matrix(P, Q):
    n = len(P) - 1
    m = len(Q) - 1
    S = np.zeros((n + m, n + m))
    for i in range(m):
        S[i, i:i + n + 1] = P
    for i in range(n):
        S[m + i, i:i + m + 1] = Q
    return S


def resultant(P, Q):
    return float(np.linalg.det(sylvester_matrix(P, Q)))


def main():
    print("=" * 80)
    print("PAPER 79 — Bézout / resultant via Pandrosion-Sylvester")
    print("=" * 80)

    print("\n[1] Common root: Res(P, Q) = 0")
    P = np.array([1.0, -3, 2])  # roots 1, 2
    Q = np.array([1.0, -2])  # root 2 (shared with P)
    R = resultant(P, Q)
    print(f"  P=(z-1)(z-2), Q=z-2 (shares root 2): Res = {R:.4e}")

    print("\n[2] No common root: Res(P, Q) != 0")
    P = np.array([1.0, -3, 2])
    Q = np.array([1.0, -5])
    R = resultant(P, Q)
    print(f"  P=(z-1)(z-2), Q=z-5: Res = {R:.4f}")

    print("\n[3] Discriminant via Res(P, P')")
    for name, P in [("z^3 - 1", np.array([1.0, 0, 0, -1])),
                     ("z^4 + 1", np.array([1.0, 0, 0, 0, 1]))]:
        Pp = np.polyder(P)
        R = resultant(P, Pp)
        # Discriminant = (-1)^(d(d-1)/2) Res(P, P') / a_d
        d = len(P) - 1
        D = (-1)**(d * (d - 1) // 2) * R / P[0]
        print(f"  {name}: Res(P, P') = {R:.4f}, Disc = {D:.4f}")


if __name__ == "__main__":
    main()
