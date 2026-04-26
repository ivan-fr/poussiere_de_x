"""
PAPER: 068 (canonical: 58pandrosion_lc.pdf)
TITLE: Linear Combination Identities for Pandrosion
STATUS: framework

THEORY
======

LINEAR COMBINATION: For polynomials P, Q and constants a, b:
  (a P + b Q)' / (a P + b Q) is NOT linear in F_P, F_Q in general.
But the residue STRUCTURE is preserved: each root of aP + bQ contributes
a pole.

For special cases (P = Q): trivially F = F_P = F_Q.
For (P, Q = P'): a F_P + b F_{P'} relates to log(P + small P').

VERIFICATION
============

  1. Linear combinations of Pandrosion fields.
  2. Logarithmic derivative of products.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 68 — Linear combination identities")
    print("=" * 80)

    P = np.array([1.0, 0, 0, -1])  # z^3 - 1
    Q = np.array([1.0, 0, 0, 0, 0])  # z^4
    z = 0.5 + 0.3j

    F_P = np.polyval(np.polyder(P), z) / np.polyval(P, z)
    F_Q = np.polyval(np.polyder(Q), z) / np.polyval(Q, z)
    print(f"\n[1] F_P({z}) = {F_P:.4f}")
    print(f"  F_Q({z}) = {F_Q:.4f}")

    # F_{PQ} = F_P + F_Q (logarithmic derivative of product)
    PQ = np.convolve(P, Q)
    F_PQ = np.polyval(np.polyder(PQ), z) / np.polyval(PQ, z)
    print(f"  F_{{PQ}}({z}) = {F_PQ:.4f}, F_P + F_Q = {F_P + F_Q:.4f}")
    print(f"  diff = {abs(F_PQ - (F_P + F_Q)):.2e}  (should be ~0: log of product)")

    print("\n[2] Sum a P + b Q: roots are different from P, Q individually")
    a, b = 0.5, 1.0
    # Pad to same length
    max_len = max(len(P), len(Q))
    P_pad = np.concatenate([np.zeros(max_len - len(P)), P])
    Q_pad = np.concatenate([np.zeros(max_len - len(Q)), Q])
    R = a * P_pad + b * Q_pad
    F_R = np.polyval(np.polyder(R), z) / np.polyval(R, z)
    print(f"  a={a}, b={b}: F_R({z}) = {F_R:.4f}")
    print(f"  This is NOT a linear combination of F_P and F_Q.")


if __name__ == "__main__":
    main()
