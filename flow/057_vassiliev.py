"""
PAPER: 057 (canonical: 47pandrosion_vassiliev.pdf)
TITLE: Vassiliev Invariants from the Pandrosion Quotient
STATUS: framework
DEPENDS: 048

THEORY
======

VASSILIEV (finite-type) invariants:
  v_n(K) = type-n invariants of knot K.
For Alexander polynomial Delta_K, the n-th derivative at 1:
  v_n = Delta_K^(n)(1) / n!  (Pandrosion Q-chain at t = 1).

Energy polynomial: E_K(t) = (Delta_K')^2 - Delta_K * Delta_K'' in Z[t].
Knot zeta: zeta_K(s) = sum 1/Delta_K'(alpha_k)^s.
By paper 003: zeta_K(1) = 0 universally.

VERIFICATION
============

  1. Compute v_2 for trefoil, figure-8.
  2. Knot zeta function zeta_K(s).
  3. Energy polynomial E_K.
"""
from __future__ import annotations
import numpy as np
import math


def main():
    print("=" * 80)
    print("PAPER 57 — Vassiliev knot invariants via Pandrosion")
    print("=" * 80)

    knots = [
        ("3_1 (trefoil)", np.array([1.0, -1, 1])),
        ("4_1 (figure-8)", np.array([1.0, -3, 1])),
        ("5_1", np.array([1.0, -1, 1, -1, 1])),
        ("5_2", np.array([2.0, -3, 2])),
        ("6_1", np.array([2.0, -5, 2])),
    ]

    print("\n[1] v_2 = Delta_K''(1)/2 (second Vassiliev)")
    print(f"  {'knot':>20} {'v_2':>8}")
    for name, P in knots:
        Pp = np.polyder(P)
        Ppp = np.polyder(Pp)
        v2 = np.polyval(Ppp, 1.0) / 2
        print(f"  {name:>20} {v2:>8.2f}")

    print("\n[2] Knot zeta zeta_K(s) = sum 1/Delta'(alpha_k)^s")
    for name, P in knots[:3]:
        roots = np.roots(P)
        Pp = np.polyder(P)
        z1 = sum(1.0 / np.polyval(Pp, ak) for ak in roots)
        z2 = sum(1.0 / np.polyval(Pp, ak)**2 for ak in roots)
        print(f"  {name}: zeta(1) = {abs(z1):.2e} (~0), zeta(2) = {z2:.4f}")

    print("\n[3] Energy polynomial E_K(t) = (Delta')^2 - Delta * Delta''")
    for name, P in knots[:3]:
        Pp = np.polyder(P)
        Ppp = np.polyder(Pp)
        Pp_sq = np.convolve(Pp, Pp)
        PPpp = np.convolve(P, Ppp)
        max_len = max(len(Pp_sq), len(PPpp))
        a = np.concatenate([np.zeros(max_len - len(Pp_sq)), Pp_sq])
        b = np.concatenate([np.zeros(max_len - len(PPpp)), PPpp])
        E = a - b
        E_at_1 = np.polyval(E, 1.0)
        print(f"  {name}: E_K(1) = {E_at_1:.4f}, deg(E) = {len(E) - 1}")


if __name__ == "__main__":
    main()
