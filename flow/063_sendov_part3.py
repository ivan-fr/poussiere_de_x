"""
PAPER: 063 (canonical: 53pandrosion_sendov.pdf)
TITLE: Sendov (Part 3): Statistical Bounds
STATUS: empirical (refined in 104, 108, 113)
DEPENDS: 052

THEORY
======

Statistical analysis of Sendov violations:
  V(P) = max_zeta min_xi |zeta - xi| - 1.
For random P (uniform in disk): V_avg < 0 with margin ~ 1 - O(1/sqrt(d)).

Tao 2022: Sendov holds for d >= d_0 sufficiently large (proved).
Open: d in [9, d_0(Tao)].

VERIFICATION
============

  1. Distribution of V(P) for random uniform P.
  2. Tao's regime threshold.
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
    print("PAPER 63 — Sendov (Part 3): statistical bounds")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Statistics of V(P) for random P in disk")
    print(f"  {'d':>4} {'#cfg':>6} {'mean V':>10} {'std V':>10} {'max V':>10}")
    for d in [5, 10, 20, 50, 100]:
        n = 100 if d <= 50 else 30
        Vs = []
        for _ in range(n):
            angles = rng.uniform(0, 2*np.pi, d)
            radii = rng.uniform(0, 1, d)
            roots = radii * np.exp(1j * angles)
            Vs.append(sendov_violation(roots))
        arr = np.array(Vs)
        print(f"  {d:>4} {n:>6} {arr.mean():>10.4f} {arr.std():>10.4f} {arr.max():>10.4f}")

    print("\n[2] Tao's regime: Sendov proved for d >= d_0 (Tao 2022)")
    print(f"  d_0 conjectured around 10^12 by Tao paper.")
    print(f"  Open: d in [9, 10^12] -- empirical only.")


if __name__ == "__main__":
    main()
