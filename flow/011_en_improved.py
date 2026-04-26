"""
PAPER: 011 (canonical: 01pandrosion_en_improved.pdf)
TITLE: The Pandrosion Energy: Improved Estimates
STATUS: proved (refinements of Pandrosion energy bounds)
DEPENDS: 000

THEORY
======

For monic P(z) = prod (z - alpha_j), the Pandrosion ENERGY is:
  E_P(z) := sum_{k=1}^d Q(alpha_k, z)^2 = (P'(z))^2 - P(z) P''(z).

KEY IDENTITY (re-affirmed): E_P = (P')^2 - P P'' (Turán-form).

IMPROVED LOWER BOUNDS:
  E_P(x) >= 0 on R if all roots are real (Turán inequality, paper 56).
  E_P(alpha_k) = P'(alpha_k)^2 = (residue inverse squared).
  prod E_P(alpha_k) = D(P)^2 (paper 22 discriminant identity).

VERIFICATION
============

  1. Energy identity E_P = (P')^2 - P P''.
  2. Energy at roots = squared derivative.
  3. Product over roots = discriminant squared.
"""
from __future__ import annotations
import math
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def energy(P, z):
    """E_P(z) = (P'(z))^2 - P(z) P''(z)."""
    Pp = np.polyder(P)
    Ppp = np.polyder(Pp)
    return np.polyval(Pp, z)**2 - np.polyval(P, z) * np.polyval(Ppp, z)


def energy_via_quotients(P, z):
    """E_P(z) = sum Q(alpha_k, z)^2."""
    roots = np.roots(P)
    return sum(Q(P, ak, z)**2 for ak in roots)


def discriminant(P):
    """|disc(P)| = prod_{i<j} |alpha_i - alpha_j|^2."""
    roots = np.roots(P)
    d = len(roots)
    log_disc = sum(2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
                   for i in range(d) for j in range(i+1, d))
    return math.exp(log_disc)


def main():
    print("=" * 80)
    print("PAPER 11 — Pandrosion Energy: improved estimates")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Energy identity: E_P(z) = (P')^2 - P P'' = sum Q(a_k, z)^2")
    for d in [3, 5, 8]:
        P = rng.standard_normal(d + 1)
        P[0] = 1.0
        z = complex(rng.standard_normal(), rng.standard_normal())
        e1 = energy(P, z)
        e2 = energy_via_quotients(P, z)
        print(f"  d={d}: E via formula = {e1:.6f}, E via Q-sum = {e2:.6f}, "
              f"diff = {abs(e1-e2):.2e}")

    print("\n[2] E_P(alpha_k) = (P'(alpha_k))^2")
    for d in [3, 5, 8]:
        P = rng.standard_normal(d + 1); P[0] = 1.0
        roots = np.roots(P)
        Pp = np.polyder(P)
        max_err = max(abs(energy(P, ak) - np.polyval(Pp, ak)**2) for ak in roots)
        print(f"  d={d}: max |E(a_k) - P'(a_k)^2| = {max_err:.2e}")

    print("\n[3] prod_k E_P(alpha_k) = |disc(P)|^2")
    for d in [3, 4, 5, 6]:
        P = rng.standard_normal(d + 1); P[0] = 1.0
        roots = np.roots(P)
        prod_e = float(np.prod([energy(P, ak).real for ak in roots]))
        disc_sq = discriminant(P)**2
        print(f"  d={d}: prod E = {prod_e:.4e}, |disc|^2 = {disc_sq:.4e}, "
              f"ratio = {prod_e/disc_sq:.6f}")


if __name__ == "__main__":
    main()
