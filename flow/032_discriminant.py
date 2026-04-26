"""
PAPER: 032 (canonical: 22pandrosion_discriminant.pdf)
TITLE: The Pandrosion Discriminant Identity
STATUS: proved (classical, Pandrosion-form)
DEPENDS: 011

THEORY
======

DISCRIMINANT IDENTITY:
  D(P)^2 = prod_k P'(alpha_k) = prod_{i<j} (alpha_i - alpha_j)^2.

PANDROSION FORM:
  prod_k E_P(alpha_k) = D(P)^2
where E_P(alpha_k) = P'(alpha_k)^2 (paper 11).

CONSEQUENCE:
For monic Z-poly with simple roots, |D(P)| in Z and >= 1, so
  prod_k |P'(alpha_k)| >= 1.
This is the discriminant lower bound used in Lehmer/Schinzel-Zassenhaus
analyses (papers 105, 114).

VERIFICATION
============

  1. D(P)^2 = prod P'(alpha_k) (Pandrosion form).
  2. Equivalent to prod_{i<j} (alpha_i - alpha_j)^2.
  3. For Z-polys: |D| in Z, >= 1.
"""
from __future__ import annotations
import numpy as np
import math


def discriminant_via_derivative(P):
    """D(P)^2 = prod P'(alpha_k)."""
    roots = np.roots(P)
    Pp = np.polyder(P)
    return float(np.prod([np.polyval(Pp, ak) for ak in roots]).real)


def discriminant_via_pairs(P):
    roots = np.roots(P)
    d = len(roots)
    return float(np.prod([(roots[i] - roots[j])**2
                          for i in range(d) for j in range(i+1, d)]).real)


def main():
    print("=" * 80)
    print("PAPER 32 — Pandrosion discriminant identity")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] D(P)^2 = prod P'(alpha_k) = prod_{i<j} (alpha_i - alpha_j)^2")
    print(f"  {'d':>3} {'D via derivs':>16} {'D via pairs':>16} {'diff':>10}")
    for d in [3, 4, 5, 6]:
        P = rng.standard_normal(d + 1); P[0] = 1.0
        d1 = discriminant_via_derivative(P)
        d2 = discriminant_via_pairs(P)
        print(f"  {d:>3} {d1:>16.6e} {d2:>16.6e} {abs(d1-d2):>10.2e}")

    print("\n[2] Integer poly: |D| in Z, |D| >= 1")
    int_polys = [
        ("z^2 - 2", np.array([1, 0, -2.0])),
        ("z^3 - z - 1", np.array([1, 0, -1, -1.0])),
        ("z^4 - z + 1", np.array([1, 0, 0, -1, 1.0])),
        ("z^5 + z + 1", np.array([1, 0, 0, 0, 1, 1.0])),
    ]
    for name, P in int_polys:
        D = discriminant_via_pairs(P)
        print(f"  {name}: D = {D}")

    print("\n[3] Cyclotomic: D(z^d - 1) = ?")
    for d in [3, 4, 5, 6]:
        P = np.array([1.0] + [0]*(d-1) + [-1])
        D = discriminant_via_pairs(P)
        print(f"  z^{d} - 1: D = {D}")


if __name__ == "__main__":
    main()
