"""
PAPER: 086 (canonical: 76pandrosion_atiyah.pdf)
TITLE: Polynomial Index Theorem via Pandrosion Field
       (also: Atiyah configurations on S^2 - initial work)
STATUS: proved (argument principle = polynomial Atiyah-Singer index)
DEPENDS: 016, 011

THEORY
======

POLYNOMIAL ATIYAH-SINGER (argument principle):
  ind(P, gamma) = (1 / 2pi i) integral_gamma F_P(z) dz = #(roots inside gamma).

CURVATURE: F_P'(z) = -T(z)/P(z)^2 where T = (P')^2 - P P'' is the Turán form.

ROUCHÉ'S THEOREM as Pandrosion perturbation: small perturbation of P
preserves the index.

VERIFICATION
============

  1. Argument principle: integral of F_P = root count.
  2. Curvature identity F_P' = -T/P^2.
  3. Rouché perturbation stability.
"""
from __future__ import annotations
import numpy as np


def winding(P, center, R, n_pts=512):
    thetas = 2 * np.pi * np.arange(n_pts) / n_pts
    z = center + R * np.exp(1j * thetas)
    Pp = np.polyder(P)
    F = np.polyval(Pp, z) / np.polyval(P, z)
    dz = 1j * R * np.exp(1j * thetas)
    integral = np.mean(F * dz) * 2 * np.pi
    return integral / (2 * np.pi * 1j)


def main():
    print("=" * 80)
    print("PAPER 86 — Polynomial Atiyah-Singer index via Pandrosion field")
    print("=" * 80)

    print("\n[1] Argument principle: ind(P, gamma) = root count")
    P = np.array([1.0, 0, 0, -1])  # z^3 - 1, roots 1, ω, ω^2
    for R in [0.5, 1.5, 2.0]:
        w = winding(P, 0, R)
        n_inside = sum(1 for r in np.roots(P) if abs(r) < R)
        print(f"  R = {R}: winding = {w.real:.4f}, roots inside = {n_inside}")

    print("\n[2] Curvature: F_P'(z) = -T(z) / P(z)^2")
    P = np.array([1.0, 0, 0, -1])
    Pp = np.polyder(P)
    Ppp = np.polyder(Pp)
    for z_val in [0.5, 0.7 + 0.3j, 2.0]:
        # F_P = P'/P, F_P' = (P'' P - (P')^2) / P^2 = -((P')^2 - P P'') / P^2 = -T/P^2
        F = np.polyval(Pp, z_val) / np.polyval(P, z_val)
        # Direct F'
        T = np.polyval(Pp, z_val)**2 - np.polyval(P, z_val) * np.polyval(Ppp, z_val)
        F_prime_via_T = -T / np.polyval(P, z_val)**2
        # Numerical F'
        h = 1e-6
        F_p = np.polyval(Pp, z_val + h) / np.polyval(P, z_val + h)
        F_m = np.polyval(Pp, z_val - h) / np.polyval(P, z_val - h)
        F_prime_num = (F_p - F_m) / (2 * h)
        print(f"  z = {z_val}: -T/P^2 = {F_prime_via_T:.4f}, F' numerical = {F_prime_num:.4f}, "
              f"diff = {abs(F_prime_via_T - F_prime_num):.2e}")

    print("\n[3] Rouché stability: small perturbation preserves index")
    Q = np.array([1.0, 0, 0, -1])  # z^3 - 1
    print("  Q = z^3 - 1: index in |z|=2 = 3")
    for eps in [0.01, 0.1, 0.5]:
        # P = Q + eps * (z^2)
        P_pert = np.array([1.0, eps, 0, -1])
        w = winding(P_pert, 0, 2.0)
        print(f"  P = z^3 + {eps} z^2 - 1: winding = {w.real:.4f}")


if __name__ == "__main__":
    main()
