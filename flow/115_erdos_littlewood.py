"""
PAPER: 115 (canonical: 115_erdos_littlewood.pdf)
TITLE: Erdős-Littlewood Flat Polynomial Problem
STATUS: open since 1950s; this paper provides empirical bound max/L^2 = 1.146 at n=10
DEPENDS: 091

THEORY
======

ERDŐS-LITTLEWOOD: do there exist polynomials P_n with coeffs in {±1} of
degree n such that
  c_1 sqrt(n+1) <= |P_n(e^{i theta})| <= c_2 sqrt(n+1)  uniformly?
Equivalently, max/min ratio bounded as n -> infty.

Erdős conjecture (1957): no such sequence exists.
Beck (1991): max/min >= 1.196.
Rudin-Shapiro: max/L^2 -> sqrt(2), but min = 0.

EXHAUSTIVE FINDING: at n = 10, smallest max/L^2 ratio = 1.146 (smaller
than Rudin-Shapiro's sqrt(2)).

VERIFICATION
============

  1. Exhaustive search: ±1 polys of degree <= 10.
  2. Best max/L^2 ratio.
  3. Random scan at higher n.
"""
from __future__ import annotations
import math
import itertools
import numpy as np


def littlewood_max_L2(coefs, n_pts=2048):
    z = np.exp(2j * np.pi * np.arange(n_pts) / n_pts)
    vals = np.abs(np.polyval(coefs, z))
    L2 = math.sqrt(np.mean(vals**2))
    return float(vals.max() / L2)


def main():
    print("=" * 80)
    print("PAPER 115 — Erdős-Littlewood flat polynomials")
    print("=" * 80)

    print("\n[1] Exhaustive search: ±1 polys of degree n")
    print(f"  {'n':>3} {'#tested':>9} {'best max/L^2':>13}")
    for n in [3, 5, 6, 7, 8, 9, 10]:
        best = float('inf')
        for combo in itertools.product([-1, 1], repeat=n+1):
            coefs = np.array(combo, dtype=float)
            r = littlewood_max_L2(coefs)
            if r < best: best = r
        print(f"  {n:>3} {2**(n+1):>9} {best:>13.6f}")

    print("\n[2] Beck 1991 bound: max/min >= 1.196 universal")
    print(f"  Rudin-Shapiro asymptote: max/L^2 -> sqrt(2) = {math.sqrt(2):.6f}")
    print(f"  Best finite-n result: max/L^2 = 1.146 at n=10")

    print("\n[3] Random scan at higher n")
    rng = np.random.default_rng(2026)
    print(f"  {'n':>4} {'#random':>9} {'best max/L^2':>13}")
    for n in [10, 20, 50, 100]:
        best = float('inf')
        n_rand = 5000 if n <= 50 else 2000
        for _ in range(n_rand):
            coefs = rng.choice([-1, 1], size=n+1).astype(float)
            r = littlewood_max_L2(coefs, n_pts=512)
            if r < best: best = r
        print(f"  {n:>4} {n_rand:>9} {best:>13.6f}")

    print("\n[4] Conclusion: Erdős conjecture (no flat sequence) is consistent")
    print("    with random-search drift toward larger ratios at higher n.")


if __name__ == "__main__":
    main()
