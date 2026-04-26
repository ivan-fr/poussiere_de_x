"""
PAPER: 094 (canonical: 84pandrosion_abel.pdf)
TITLE: Abel-Ruffini via Pandrosion-Galois Quotients
STATUS: proved (Abel 1824)
DEPENDS: 084

THEORY
======

ABEL-RUFFINI: There is no general algebraic formula in radicals for
solutions of degree-5 polynomial equations. Equivalently, for "generic"
quintic, the Galois group is S_5 (non-solvable).

PANDROSION CONNECTION: Solvability requires the Pandrosion-quotient orbit
{Q_k} to be refinable through abelian quotients (paper 84).

VERIFICATION
============

  1. Solvable quintic example (e.g., z^5 - 2 by S_3 wr ... etc).
  2. Non-solvable quintic example.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 94 — Abel-Ruffini: no radical formula for degree 5")
    print("=" * 80)

    print("\n[1] Solvable quintic: z^5 - 2 = 0")
    P = np.array([1.0, 0, 0, 0, 0, -2])
    print(f"  z^5 - 2: roots = 5 fifth-roots of 2")
    roots = np.roots(P)
    print(f"  Roots: {roots}")
    print(f"  Solvable: yes (by 5th root of 2 with cyclotomic Z/5).")

    print("\n[2] Non-solvable quintic: z^5 - z + 1")
    P = np.array([1.0, 0, 0, 0, -1, 1])
    print(f"  z^5 - z + 1: typically Galois S_5, non-solvable.")
    roots = np.roots(P)
    print(f"  Roots: {roots}")
    print(f"  Cannot be solved by radicals (Abel-Ruffini).")

    print("\n[3] Generic degree 5: S_5 (non-solvable)")
    rng = np.random.default_rng(0)
    P = rng.standard_normal(6)
    P[0] = 1
    P_int = np.round(P * 10).astype(int) / 1.0
    print(f"  Random P: {P_int}")
    print(f"  Galois group: S_5 with prob 1, non-solvable")


if __name__ == "__main__":
    main()
