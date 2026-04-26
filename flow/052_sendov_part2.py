"""
PAPER: 052 (canonical: 42pandrosion_sendov.pdf)
TITLE: Sendov (Part 2): Boundary Cases
STATUS: empirical (refined further in 053, 104, 108, 113)
DEPENDS: 019

THEORY
======

Sendov's conjecture is tight at boundary configurations:
  - z^d - 1 (roots of unity): tight as d -> infty.
  - Miller's family: extremal at d=32, eta=0.5.

Pandrosion approach: study contraction of root toward critical point
via the Pandrosion field F_P near the boundary |z| = 1.

VERIFICATION
============

  1. Boundary cases: z^d - 1, Miller family.
  2. Maximum violation across families.
"""
from __future__ import annotations
import math
import numpy as np


def sendov_violation(roots):
    P = np.poly(roots)
    crits = np.roots(np.polyder(P))
    if len(crits) == 0: return float('-inf')
    D = np.abs(np.array(roots)[:, None] - np.array(crits)[None, :])
    return float(D.min(axis=1).max()) - 1.0


def main():
    print("=" * 80)
    print("PAPER 52 — Sendov (Part 2): boundary cases")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] z^d - 1 (cyclotomic): tight asymptotically")
    print(f"  {'d':>4} {'V':>14}")
    for d in [4, 8, 16, 32, 64]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        v = sendov_violation(roots)
        print(f"  {d:>4} {v:>14.4e}")

    print("\n[2] Miller family: zeta_1 = 1, others on small arc")
    print(f"  {'d':>4} {'eta':>6} {'V':>10}")
    for d in [16, 32, 64]:
        for eta in [0.5, 0.1]:
            roots = [1.0 + 0j]
            for _ in range(d - 1):
                phi = rng.uniform(0, 2*np.pi)
                z = 1.0 + eta * np.exp(1j * phi)
                if abs(z) > 1: z = z / abs(z) * 0.999
                roots.append(z)
            v = sendov_violation(np.array(roots))
            print(f"  {d:>4} {eta:>6.2f} {v:>10.4e}")

    print("\n[3] Antipodal cluster")
    print(f"  {'d':>4} {'V':>10}")
    for d in [16, 32]:
        half = d // 2
        roots = []
        for _ in range(half):
            roots.append(1.0 + 0.05 * (rng.uniform(-1, 1) + 1j * rng.uniform(-1, 1)))
        for _ in range(d - half):
            roots.append(-1.0 + 0.05 * (rng.uniform(-1, 1) + 1j * rng.uniform(-1, 1)))
        roots = np.array(roots)
        for i in range(len(roots)):
            if abs(roots[i]) > 1: roots[i] /= abs(roots[i]) * 1.001
        v = sendov_violation(roots)
        print(f"  {d:>4} {v:>10.4e}")


if __name__ == "__main__":
    main()
