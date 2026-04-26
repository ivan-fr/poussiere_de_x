"""
PAPER: 083 (canonical: 73pandrosion_egz.pdf)
TITLE: Erdős-Ginzburg-Ziv Theorem and Pandrosion Generators
STATUS: proved (Erdős-Ginzburg-Ziv 1961)

THEORY
======

EGZ (1961): Among any 2n-1 integers, there exist n whose sum is divisible
by n.

POLYNOMIAL VERSION (via generating functions):
prod (1 + x^{a_i}) has properties tied to the additive structure.

VERIFICATION
============

  1. EGZ for small n.
  2. Counterexamples to "2n-2" (tightness).
"""
from __future__ import annotations
import itertools


def egz_check(numbers, n):
    """Check if there exists a subset of size n summing to 0 mod n."""
    for combo in itertools.combinations(numbers, n):
        if sum(combo) % n == 0:
            return True, list(combo)
    return False, None


def main():
    print("=" * 80)
    print("PAPER 83 — Erdős-Ginzburg-Ziv theorem")
    print("=" * 80)

    print("\n[1] EGZ for various n: 2n-1 numbers, n with sum 0 mod n")
    print(f"  {'n':>3} {'numbers':>30} {'subset':>20}")
    cases = [
        (3, [1, 1, 2, 2, 3]),
        (4, [1, 2, 3, 4, 5, 6, 7]),
        (5, [1, 1, 1, 1, 1, 2, 3, 4, 5]),
    ]
    for n, nums in cases:
        ok, subset = egz_check(nums, n)
        print(f"  {n:>3} {str(nums):>30} {str(subset):>20}")

    print("\n[2] EGZ tightness: only 2n-2 numbers can fail")
    # Example: n=3, take {0, 0, 1, 1}: 4 = 2n-2 numbers, no 3 summing to 0 mod 3
    print("  n=3: 4 = 2*3 - 2 numbers {0,0,1,1}: 3-subsets are (0,0,1) sum 1, (0,1,1) sum 2: no 0 mod 3")
    nums = [0, 0, 1, 1]
    ok, subset = egz_check(nums, 3)
    print(f"  {nums}: 3-subset summing to 0 mod 3? {ok}")


if __name__ == "__main__":
    main()
