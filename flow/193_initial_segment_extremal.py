"""
PAPER: 193 (NEW — initial segment extremality)
TITLE: Initial-segment extremality for shifted Q coefficients:
       reducing all lattice subsets to M={1,...,k}
STATUS: RH remains OPEN.  Paper 192 found that the smallest normalized
        shifted coefficient B_j/e_k in bounded searches always occurs at
        the initial segment {1,...,k}, and at j=0.  This paper stress-tests
        that monotonicity and formulates the extremal lemma needed to close
        the Wronskian lattice barrier.
DEPENDS: 192 (shifted coefficient formula).

THEORY
======

Let a_i = 2*pi*m_i^2 and

    B_j(M) = [x^j] q_M(1+x).

Paper 192 showed B_j(M)>0 in all tested cases.  Stronger empirical fact:

    B_j(M)/e_k(a_M) is minimized by M={1,...,k}, j=0.

If true, the all-subset theorem reduces to proving

    B_0({1,...,k}) > 0  for every k.

This paper tests monotonicity under:
  - raising one m_i;
  - replacing M by non-consecutive subsets;
  - random large subsets.
"""
from __future__ import annotations

import itertools
import math
import random


def odd_double_factorial(n):
    out = 1.0
    for m in range(n, 0, -2):
        out *= m
    return out


def elementary_symmetric(xs, r):
    if r == 0:
        return 1.0
    return sum(math.prod(c) for c in itertools.combinations(xs, r))


def all_B(ms):
    k = len(ms)
    a = [2 * math.pi * m * m for m in ms]
    out = []
    for j in range(k + 1):
        total = 0.0
        for r in range(j, k + 1):
            total += (math.comb(r, j)
                      * ((-1) ** (k - r))
                      * odd_double_factorial(2 * (k - r) + 1)
                      * elementary_symmetric(a, r))
        out.append(total)
    return out, elementary_symmetric(a, k)


def min_normalized(ms):
    B, top = all_B(ms)
    vals = [b / top for b in B]
    j = min(range(len(vals)), key=lambda i: vals[i])
    return vals[j], j


def main():
    print("=" * 80)
    print("PAPER 193 — initial segment extremality")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Exhaustive verification of extremality
    # ------------------------------------------------------------------
    print("\n[1] Exhaustive extremality in bounded boxes")
    print(f"  {'M_max':>5} {'k':>3} {'min subset':>20} {'j':>3}"
          f" {'min ratio':>14} {'initial ratio':>16}")
    for M_max in [10, 12]:
        for k in range(2, min(7, M_max) + 1):
            best = (float("inf"), None, None)
            for ms in itertools.combinations(range(1, M_max + 1), k):
                val, j = min_normalized(ms)
                if val < best[0]:
                    best = (val, ms, j)
            init = tuple(range(1, k + 1))
            init_val, _ = min_normalized(init)
            print(f"  {M_max:>5} {k:>3} {str(best[1]):>20} {best[2]:>3}"
                  f" {best[0]:>14.8f} {init_val:>16.8f}")

    # ------------------------------------------------------------------
    # [2] Raising one coordinate
    # ------------------------------------------------------------------
    print("\n[2] Raising one coordinate from the initial segment")
    print(f"  {'k':>3} {'move':>12} {'ratio before':>16} {'ratio after':>16}")
    for k in range(3, 10):
        base = list(range(1, k + 1))
        before, _ = min_normalized(tuple(base))
        moved = list(base)
        moved[-1] = k + 5
        after, _ = min_normalized(tuple(moved))
        print(f"  {k:>3} {f'{k}->{k+5}':>12} {before:>16.8f} {after:>16.8f}")

    # ------------------------------------------------------------------
    # [3] Random large subsets
    # ------------------------------------------------------------------
    print("\n[3] Random large-subset stress test")
    rng = random.Random(2026)
    print(f"  {'k':>3} {'trials':>8} {'min random ratio':>18}"
          f" {'initial ratio':>16} {'winner':>10}")
    for k in [5, 8, 12]:
        init = tuple(range(1, k + 1))
        init_val, _ = min_normalized(init)
        best = (float("inf"), None)
        trials = 120
        universe = list(range(1, 6 * k + 1))
        for _ in range(trials):
            ms = tuple(sorted(rng.sample(universe, k)))
            val, _ = min_normalized(ms)
            if val < best[0]:
                best = (val, ms)
        winner = "initial" if init_val <= best[0] + 1e-12 else "random"
        print(f"  {k:>3} {trials:>8} {best[0]:>18.8f}"
              f" {init_val:>16.8f} {winner:>10}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  EMPIRICAL EXTREMAL LEMMA:")
    print("    For fixed k, the minimum of B_j(M)/e_k(M) over all M and j")
    print("    appears at M={1,...,k}, j=0.")
    print()
    print("  IF PROVED:")
    print("    The Wronskian lattice barrier reduces to a one-parameter")
    print("    inequality B_0({1,...,k})>0.")
    print()
    print("  NEXT:")
    print("    Analyze B_0({1,...,k}) asymptotically using products over m^2.")


if __name__ == "__main__":
    main()
