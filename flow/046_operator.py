"""
PAPER: 046 (canonical: 36pandrosion_op.pdf)
TITLE: Operator-Theoretic View of Pandrosion
STATUS: framework

THEORY
======

The Pandrosion operator P_F: C^n -> C^n acts on test points:
  P_{F, z_0}(z) = z_0 - Q_F(z_0, z)^{-1} F(z_0).

Viewed as a family of maps parametrised by z_0, this gives a sheaf-like
structure on C^n. Fixed points form the discrete set V(F) = {zeros of F}.

For the univariate case, P_a(z) = a - P(a)/Q(a, z) is a Möbius
transformation in z (with parameter a).

VERIFICATION
============

  1. Möbius structure: P_a(z) = (Az + B)/(Cz + D) for fixed a.
  2. Action on infinity.
  3. Fixed-point structure.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 46 — Operator-theoretic view of Pandrosion")
    print("=" * 80)

    # Univariate: P_a(z) = a - P(a)/Q(a, z) is Möbius in z
    print("\n[1] Möbius structure: P_a(z) for fixed a is rational of degree 1")
    P = np.array([1.0, 0, -2])  # z^2 - 2
    a = 1.5
    Pa_val = np.polyval(P, a)
    Pp = np.polyder(P)
    # Q(a, z) = (P(z) - P(a))/(z - a) for z != a
    # For z^2 - 2: P(z) - P(a) = z^2 - a^2 = (z-a)(z+a), Q = z + a
    # So h_a(z) = a - P(a)/(z + a) = a - (a^2 - 2)/(z + a)
    print(f"  P = z^2 - 2, a = {a}")
    print(f"  Q(a, z) = z + a (linear in z)")
    print(f"  h_a(z) = a - (a^2 - 2)/(z + a) -- Möbius transformation")
    print(f"  Fixed points: roots of P, i.e., ±sqrt(2)")
    for zeta in [np.sqrt(2), -np.sqrt(2)]:
        Q_val = zeta + a
        h = a - Pa_val / Q_val
        print(f"  h_a({zeta:.4f}) = {h:.6f}")

    print("\n[2] Action on infinity: limit z -> infty")
    # Q(a, z) -> P'(z_0) where z_0 = a as |z| -> infty (NO; Q is degree d-1 in z)
    # For d=2: Q -> z (leading term), so h_a(z) = a - P(a)/z -> a as z -> infty
    print("  For P degree d, Q(a, z) ~ z^{d-1} as z -> infty")
    print("  So h_a(z) -> a (the anchor) as z -> infty.")
    z_far = 1000.0
    h_far = a - Pa_val / (z_far + a)
    print(f"  z = {z_far}: h_a(z) = {h_far:.6f} (close to a = {a})")

    print("\n[3] Iteration with anchor sweeping")
    P = np.array([1.0, 0, 0, -1])  # z^3 - 1
    z = 0.5 + 0.5j
    print(f"  z^3 - 1, sweeping anchors:")
    for anchor in [0.7+0.3j, 0.9+0.1j, 1.0+0j]:
        z_next = anchor - np.polyval(P, anchor) / ((np.polyval(P, z) - np.polyval(P, anchor)) / (z - anchor))
        print(f"  anchor = {anchor:.3f}: h(z) = {z_next:.4f}, |z_next - 1| = {abs(z_next - 1):.4f}")


if __name__ == "__main__":
    main()
