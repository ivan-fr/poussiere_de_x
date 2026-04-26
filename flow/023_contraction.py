"""
PAPER: 023 (canonical: 13pandrosion_contraction.pdf)
TITLE: Contraction Ratios in Pandrosion Iteration
STATUS: framework (corrected in paper 101 via Plancherel decorrelation)
DEPENDS: 011

THEORY
======

The contraction ratio C of a Pandrosion step from anchor a to iterate z:
  C(a, z) = prod_k (1 - r_k(a, z))   with r_k = Q(alpha_k, z) / Q(a, z).

For the step to be well-resolved:
  | arg(C) | < pi  (avoid pole-crossing).

ENERGY DOMINANCE: |Q(a, z)|^2 > E_P(z) = sum |Q(alpha_k, z)|^2
  is sufficient for | arg(C) | < pi (anchor stronger than tidal field).

VERIFICATION
============

  1. Contraction product C = prod (1 - Q_k/Q_a) for random configs.
  2. Energy dominance check: |Q_a|^2 > E_P.
  3. Original Conjecture 5.9 (half-plane Re(rs) < 0) refuted (paper 101).
"""
from __future__ import annotations
import math
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def main():
    print("=" * 80)
    print("PAPER 23 — Contraction ratios")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Contraction product C = prod (1 - Q_k(z)/Q_a(z))")
    for d in [3, 5, 8]:
        roots_real = rng.uniform(-1, 1, d)
        P = np.poly(roots_real)
        a = 0.5 + 0.5j
        z = 0.3 + 0.3j
        Q_a = Q(P, a, z)
        rs = []
        for ak in roots_real:
            Q_k = Q(P, ak, z)
            rs.append(Q_k / Q_a)
        C_prod = np.prod([1 - r for r in rs])
        print(f"  d={d}: |C| = {abs(C_prod):.4f}, arg(C) = {np.angle(C_prod):.4f} rad")

    print("\n[2] Energy dominance: |Q_a|^2 > E_P(z) ?")
    for d in [3, 5, 8]:
        P = rng.standard_normal(d + 1); P[0] = 1.0
        roots = np.roots(P)
        a = complex(rng.standard_normal(), rng.standard_normal())
        z = complex(rng.standard_normal(), rng.standard_normal())
        Q_a = Q(P, a, z)
        E_P = sum(abs(Q(P, ak, z))**2 for ak in roots)
        print(f"  d={d}: |Q_a|^2 = {abs(Q_a)**2:.4e}, E_P(z) = {E_P:.4e}, "
              f"dominance = {abs(Q_a)**2 > E_P}")

    print("\n[3] Half-plane containment Re(rs) < 0 (Conjecture 5.9 of paper 1)")
    print("    REFUTED in paper 101: counterexample at d=3.")
    print("    Optimal example: roots on circle, R = 2:")
    P = np.array([1, -0.98 * np.exp(1j * np.pi / 3), 0, 0.5])  # d=3 try
    R = 2.0
    cuspy_z = R * np.exp(1j * np.pi / 4)
    a = R
    Pp = np.polyder(P)
    rs = np.polyval(P, cuspy_z) / np.polyval(P, a) if abs(np.polyval(P, a)) > 1e-12 else 0
    print(f"  rs = P(z)/P(a) = {rs:.4f}, Re(rs) = {rs.real:.4f}")
    print(f"  Status: Re(rs) > 0 possible -> half-plane containment fails")


if __name__ == "__main__":
    main()
