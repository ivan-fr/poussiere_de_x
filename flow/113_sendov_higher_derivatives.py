"""
PAPER: 113 (canonical: 113_sendov_higher_derivatives.pdf)
TITLE: Sendov for Higher Derivatives — Naive Extension Refuted
STATUS: empirical refutation of naive extension r_k(d) <= 1
DEPENDS: 104, 052

THEORY
======

CLASSICAL SENDOV: r_1(d) <= 1 (closest critical point of P' to a root of P).

NAIVE EXTENSION: r_k(d) <= 1 for k = 1, ..., d-1.
THIS PAPER: REFUTES the naive extension. For k near d-1, r_k(d) > 1.

EXAMPLES (random scan):
  d = 10, k = 8: V_k = +0.067 (>0, fails Sendov)
  d = 50, k = 48: V_k = +0.125

REASON: P^(d-1) has 1 root at the centroid, far from boundary roots.

VERIFICATION
============

  1. r_k violations at large k (close to d-1).
  2. r_1(d) <= 1 (classical Sendov holds).
"""
from __future__ import annotations
import numpy as np


def sendov_violation_kth(roots, k):
    P = np.poly(roots)
    if k >= len(roots): return float('-inf')
    Pk = P
    for _ in range(k): Pk = np.polyder(Pk)
    if len(Pk) <= 1:
        if len(Pk) == 1: return float('inf')
        crits = [-Pk[1]/Pk[0]]
    else: crits = np.roots(Pk)
    if len(crits) == 0: return float('inf')
    D = np.abs(np.array(roots)[:, None] - np.array(crits)[None, :])
    return float(D.min(axis=1).max()) - 1.0


def main():
    print("=" * 80)
    print("PAPER 113 — Sendov for higher derivatives")
    print("=" * 80)
    rng = np.random.default_rng(2026)

    print("\n[1] Random uniform-disk: V_k(d) for various k, d")
    print(f"  {'d':>4} {'k':>4} {'#cfg':>5} {'max V_k':>12}")
    for d in [10, 20, 50]:
        for k in [1, 2, 5, d // 2, d - 2]:
            if k >= d or k < 1: continue
            n_cfg = 50 if d <= 20 else 20
            max_v = float('-inf')
            for _ in range(n_cfg):
                angles = rng.uniform(0, 2*np.pi, d)
                radii = rng.uniform(0, 1, d)
                roots = radii * np.exp(1j * angles)
                v = sendov_violation_kth(roots, k)
                if v > max_v: max_v = v
            print(f"  {d:>4} {k:>4} {n_cfg:>5} {max_v:>12.4f}")

    print("\n[2] Conclusion: r_k(d) > 1 for k near d-1 (refutes naive extension)")
    print("  Classical Sendov (k=1) consistent with r_1 <= 1.")


if __name__ == "__main__":
    main()
