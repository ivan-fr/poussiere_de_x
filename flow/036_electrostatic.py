"""
PAPER: 036 (canonical: 26pandrosion_electro.pdf)
TITLE: Stieltjes Electrostatics and the Pandrosion Quotient
STATUS: framework (complete dictionary)
DEPENDS: 011, 023

THEORY
======

ELECTROSTATIC DICTIONARY:
  Force from charge alpha_k:    F_k = 1/(z - alpha_k) = Q(alpha_k, z)/P(z)
  Total field:                  F = P'/P = sum F_k
  Equilibrium (F = 0):          P'(beta) = 0  (critical point)
  Influence weight:             lambda_k = |F_k|^2 / sum |F_j|^2
  Tidal field:                  T = sum F_k^2 = E_P/P^2
  k-body interaction:           e_k(F) = P^(k)/(k! P)
  Self-energy at root:          E_P(alpha_k) = P'(alpha_k)^2

GAUSS-LUCAS = equilibrium containment.
LAGUERRE = T >= 0 on R (Pandrosion energy nonneg).

VERIFICATION
============

  1. Gauss-Lucas = electrostatic equilibrium (paper 24 link).
  2. Laguerre = T >= 0 on R (paper 30 link).
  3. Tidal field at critical points.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 36 — Stieltjes electrostatic dictionary")
    print("=" * 80)
    rng = np.random.default_rng(0)

    P = np.array([1, 0, 0, -1.0])  # z^3 - 1
    roots = np.roots(P)
    crits = np.roots(np.polyder(P))

    print(f"\n[1] z^3 - 1: roots {roots}, critical points {crits}")

    print("\n[2] Equilibrium check at critical points: sum F_k(beta) = 0")
    for beta in crits:
        F_total = sum(1.0 / (beta - ak) for ak in roots)
        print(f"  beta = {beta:.6f}: |sum F_k| = {abs(F_total):.2e}")

    print("\n[3] Influence weights at critical points")
    for beta in crits:
        Fs = [1.0 / (beta - ak) for ak in roots]
        weights = [abs(F)**2 for F in Fs]
        s = sum(weights)
        weights = [w / s for w in weights]
        print(f"  beta = {beta:.4f}: weights = {[f'{w:.4f}' for w in weights]}")

    print("\n[4] Tidal field T(z) = sum F_k(z)^2 on real line for real-rooted P")
    P_real = np.poly([1.0, 2.0, 3.0])  # real-rooted
    roots_real = np.roots(P_real)
    print(f"  P with roots {roots_real}")
    for x in [0.0, 1.5, 2.5, 4.0]:
        T = sum(1.0 / (x - ak)**2 for ak in roots_real)
        print(f"  x = {x}: T(x) = {T:.4f}  (>= 0 on R: Laguerre)")


if __name__ == "__main__":
    main()
