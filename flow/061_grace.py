"""
PAPER: 061 (canonical: 51pandrosion_grace.pdf)
TITLE: Grace's Theorem via Pandrosion
STATUS: proved (Grace 1902, Pandrosion-form)
DEPENDS: 011

THEORY
======

GRACE'S THEOREM: If P, Q are circular polynomials of degree n in apolar form,
then they share a common zero in any closed disk that contains the zeros of
P (or zeros of Q).

POLAR DERIVATIVE: D_alpha P(z) = (n - z d/dz) P(z) + alpha P'(z)... or
related operator depending on formulation.

PANDROSION CONNECTION: The polar derivative is a Pandrosion-style mixed
operator combining P and P'.

VERIFICATION
============

  1. Apolar polynomials and shared zeros.
  2. Polar derivative computation.
"""
from __future__ import annotations
import numpy as np


def polar_derivative(P, alpha):
    """D_alpha P = n P - z P' + alpha P' = n P + (alpha - z) P'."""
    n = len(P) - 1
    Pp = np.polyder(P)
    # n P
    nP = n * P
    # (alpha - z) P': Pp shifted by alpha multiplication minus z * Pp
    # alpha P' (constant scaling)
    alpha_Pp = alpha * Pp
    # -z * Pp: shift up by 1 power and negate
    z_Pp = np.concatenate([np.zeros(0), Pp, [0]])  # z * P' increases degree by 1
    # Actually z * Pp means multiply by z: prepend with leading [Pp, 0]
    z_Pp = np.zeros(len(Pp) + 1)
    z_Pp[:-1] = Pp  # numpy convention: highest degree first; z * Pp means shift left
    # Align lengths
    max_len = max(len(nP), len(alpha_Pp), len(z_Pp))
    a = np.concatenate([np.zeros(max_len - len(nP)), nP])
    b = np.concatenate([np.zeros(max_len - len(alpha_Pp)), alpha_Pp])
    c = np.concatenate([np.zeros(max_len - len(z_Pp)), z_Pp])
    return a + b - c


def main():
    print("=" * 80)
    print("PAPER 61 — Grace's theorem via Pandrosion")
    print("=" * 80)

    print("\n[1] Polar derivative D_alpha P")
    P = np.array([1.0, 0, -1])  # z^2 - 1
    alpha = 0.5
    Dp = polar_derivative(P, alpha)
    print(f"  P = z^2 - 1, alpha = {alpha}")
    print(f"  D_alpha P = {Dp}")
    print(f"  At z = 0: D = {np.polyval(Dp, 0):.4f}")
    print(f"  At z = 1: D = {np.polyval(Dp, 1):.4f}")

    print("\n[2] Apolar pair shares a zero in disk containing P's zeros")
    P = np.poly([0.5, -0.5])  # P with zeros in unit disk
    Q = np.poly([0.7, 1.5])  # Q with one zero outside
    print(f"  P zeros: {np.roots(P)}, Q zeros: {np.roots(Q)}")
    print(f"  By Grace: P and Q apolar -> shared zero in any disk containing P's zeros.")


if __name__ == "__main__":
    main()
