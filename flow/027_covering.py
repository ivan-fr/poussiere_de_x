"""
PAPER: 027 (canonical: 17pandrosion_covering.pdf)
TITLE: Covering Numbers for Polynomial Roots via Pandrosion
STATUS: framework

THEORY
======

The COVERING NUMBER N(P, eps) is the minimum number of disks of radius eps
needed to cover all roots of P. By Smale gamma-theory:
  N(P, eps) <= d * (gamma_max / eps)^k
for some k depending on configuration.

Pandrosion connection: gamma(P, z) = sup_{k>=2} |P^(k)(z)/(k! P'(z))|^{1/(k-1)}
appears in basin radius eta_F = c/(gamma * |P|).

VERIFICATION
============

  1. Cover roots of various polynomials with eps-balls.
  2. Compare with gamma-theoretic prediction.
"""
from __future__ import annotations
import math
import numpy as np


def covering_number(roots, eps):
    """Greedy: count balls of radius eps that cover all roots."""
    n = len(roots)
    covered = [False] * n
    n_balls = 0
    while not all(covered):
        # Pick first uncovered
        for i in range(n):
            if not covered[i]:
                center = roots[i]
                break
        # Mark all roots within eps of center
        for j in range(n):
            if not covered[j] and abs(roots[j] - center) <= eps:
                covered[j] = True
        n_balls += 1
    return n_balls


def smale_gamma(P, z):
    Pp_z = np.polyval(np.polyder(P), z)
    if abs(Pp_z) < 1e-15: return float('inf')
    d = len(P) - 1
    gammas = []
    for k in range(2, d + 1):
        Pk = np.polyval(np.polyder(P, k), z)
        v = abs(Pk / (math.factorial(k) * Pp_z))
        if v > 0: gammas.append(v**(1.0 / (k - 1)))
    return max(gammas) if gammas else 0


def main():
    print("=" * 80)
    print("PAPER 27 — Covering numbers via Pandrosion gamma")
    print("=" * 80)

    # Various polynomials
    test_polys = [
        ("z^d - 1, d=8", np.array([1.0] + [0]*7 + [-1])),
        ("z^d - 1, d=12", np.array([1.0] + [0]*11 + [-1])),
        ("Wilkinson d=10", np.poly([1, 2, 3, 4, 5, 6, 7, 8, 9, 10.0])),
    ]
    for name, P in test_polys:
        roots = np.roots(P)
        print(f"\n  {name}")
        for eps in [0.1, 0.3, 0.5, 1.0]:
            N = covering_number(roots, eps)
            print(f"    eps = {eps}: N(P, eps) = {N}  /  d = {len(roots)}")

    print("\n[2] Smale gamma at random points")
    rng = np.random.default_rng(0)
    for d in [5, 10, 20]:
        P = rng.standard_normal(d + 1); P[0] = 1.0
        for _ in range(3):
            z = complex(rng.standard_normal(), rng.standard_normal())
            g = smale_gamma(P, z)
            print(f"  d={d}, z={z:.3f}: gamma = {g:.4f}")


if __name__ == "__main__":
    main()
