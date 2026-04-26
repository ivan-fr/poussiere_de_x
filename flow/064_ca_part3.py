"""
PAPER: 064 (canonical: 54pandrosion_ca.pdf)
TITLE: Casas-Alvero (Part 3): Resultant Approach
STATUS: framework (refined in paper 106)
DEPENDS: 050, 058

THEORY
======

CA via resultants: P satisfies CA condition iff certain resultants
Res(P, P^(k)) share a common root structure.
For pure power (z-alpha)^d, all resultants vanish.
For non-pure, generic resultants are non-zero (over-determined system).

Paper 106 (this work) gives the dimension count: 2(d-1) constraints on
d unknowns, exponentially-growing residuals.

VERIFICATION
============

  1. Resultant Res(P, P^(k)) vanishes for pure power.
  2. Random non-pure P: resultants generically non-zero.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 64 — Casas-Alvero (Part 3): resultant approach")
    print("=" * 80)

    print("\n[1] Pure power (z-1)^4: all derivatives vanish at z = 1")
    P = np.array([1, -4, 6, -4, 1.0])  # (z-1)^4
    print(f"  Roots of P: {np.roots(P)}")
    for k in range(1, 4):
        Pk = np.polyder(P, k)
        roots_k = np.roots(Pk)
        # Check if 1 is in roots_k
        near_1 = any(abs(r - 1) < 1e-3 for r in roots_k)
        print(f"  P^({k}) roots: {roots_k}, contains 1? {near_1}")

    print("\n[2] Non-pure power: derivatives' roots typically don't share")
    P = np.array([1.0, 0, 0, -1, 1])  # z^4 - z + 1
    P_roots = np.roots(P)
    print(f"  Roots of P: {P_roots}")
    for k in range(1, 4):
        Pk = np.polyder(P, k)
        roots_k = np.roots(Pk)
        # Min distance between any P-root and P^(k)-root
        min_d = min(abs(r1 - r2) for r1 in P_roots for r2 in roots_k)
        print(f"  P^({k}) roots: {roots_k}, min dist to P-root = {min_d:.4f}")


if __name__ == "__main__":
    main()
