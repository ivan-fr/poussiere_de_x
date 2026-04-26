"""
PAPER: 065 (canonical: 55pandrosion_sz.pdf)
TITLE: Schinzel-Zassenhaus Initial (Pandrosion-form)
STATUS: framework (refined in paper 114)
DEPENDS: 020

THEORY
======

SCHINZEL-ZASSENHAUS (1965): For non-cyclotomic monic Z-poly of degree d,
max_j |alpha_j| >= 1 + c/d for absolute c > 0.

DIMITROV (2019): max_j |alpha_j|^d >= sqrt(2), so log max >= log 2 / (2d).

Pandrosion-Hadamard connection (paper 77): the discriminant identity
det G = |Disc|^2 with Disc(P) integer >= 1 for distinct roots.

VERIFICATION
============

  1. Smallest |alpha|_max for non-cyclotomic small-height polys.
  2. Compare with Dimitrov bound.
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def main():
    print("=" * 80)
    print("PAPER 65 — Schinzel-Zassenhaus initial")
    print("=" * 80)

    print("\n[1] Smallest max|alpha| for height-1 monic Z-polys")
    print(f"  {'d':>3} {'#non-cyclo':>12} {'min max|alpha| - 1':>22}")
    for d in [4, 6, 8, 10]:
        min_excess = float('inf')
        n_nc = 0
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            roots = np.roots(coefs)
            M = float(np.prod(np.maximum(1.0, np.abs(roots))))
            if M < 1.001: continue  # cyclotomic
            n_nc += 1
            Mmax = float(np.max(np.abs(roots)))
            excess = Mmax - 1.0
            if excess < min_excess: min_excess = excess
        print(f"  {d:>3} {n_nc:>12} {min_excess:>22.6f}")

    print("\n[2] Dimitrov bound: max|alpha| >= 2^(1/(2d))")
    for d in [5, 10, 20, 50]:
        bd = 2**(1.0 / (2 * d)) - 1
        print(f"  d = {d}: Dimitrov bound = {bd:.6f}")


if __name__ == "__main__":
    main()
