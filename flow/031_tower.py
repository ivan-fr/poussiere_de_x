"""
PAPER: 031 (canonical: 21pandrosion_tower.pdf)
TITLE: Pandrosion Tower / k-body interactions
STATUS: framework

THEORY
======

The Pandrosion TOWER constructs higher-derivative interactions:
  e_k(F_P) = P^(k) / (k! * P)
where F_P = P'/P. Each level e_k corresponds to a k-body interaction in
the electrostatic analogy.

Identity: e_2(F) - F^2 = -F'/2 (paper 26 Stieltjes electrostatics).

VERIFICATION
============

  1. Compute e_k(F_P) = P^(k)/(k! P) for various k.
  2. Verify the relation e_2 - F^2 = -F'/2.
"""
from __future__ import annotations
import math
import numpy as np


def e_k(P, k, z):
    """e_k(F_P)(z) = P^(k)(z) / (k! * P(z))."""
    Pk = np.polyder(P, k)
    return np.polyval(Pk, z) / (math.factorial(k) * np.polyval(P, z))


def main():
    print("=" * 80)
    print("PAPER 31 — Pandrosion tower e_k = P^(k)/(k! P)")
    print("=" * 80)
    rng = np.random.default_rng(0)

    P = rng.standard_normal(8)
    P[0] = 1.0
    z = 0.5 + 0.5j

    print("\n[1] e_k for k = 1, 2, 3, 4")
    for k in range(1, 5):
        ek = e_k(P, k, z)
        print(f"  e_{k} = P^({k})(z)/({k}! P(z)) = {ek}")

    print("\n[2] Identity: e_2 - F^2 = -F'/2 (where F = e_1 = P'/P)")
    F = e_k(P, 1, z)
    e2 = e_k(P, 2, z)
    # F' = (P'/P)' = (P''P - P'^2)/P^2 = e_2/(1/2) ... let me compute directly
    # F'(z) = P''(z)/P(z) - (P'(z)/P(z))^2 = 2 e_2 - F^2
    # So e_2 - F^2 = -F'/2 ⟺ 2 e_2 - 2 F^2 = -F' ⟺ F' = 2 F^2 - 2 e_2
    Pp = np.polyval(np.polyder(P), z)
    Ppp = np.polyval(np.polyder(P, 2), z)
    Pz = np.polyval(P, z)
    F_prime = Ppp/Pz - (Pp/Pz)**2
    lhs = e2 - F**2
    rhs = -F_prime / 2
    print(f"  e_2 - F^2 = {lhs}")
    print(f"  -F'/2     = {rhs}")
    print(f"  diff      = {abs(lhs - rhs):.2e}")

    print("\n[3] Tower for d=5 polynomial: e_1, ..., e_5")
    P = rng.standard_normal(6); P[0] = 1.0
    z = complex(rng.standard_normal(), rng.standard_normal())
    for k in range(1, 6):
        print(f"  k={k}: e_{k}(z) = {e_k(P, k, z)}")


if __name__ == "__main__":
    main()
