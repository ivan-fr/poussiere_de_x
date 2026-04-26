"""
PAPER: 192 (NEW — shifted coefficient formula)
TITLE: Explicit coefficient formula for q_M(1+x): reducing the Wronskian
       root barrier to elementary symmetric inequalities
STATUS: RH remains OPEN.  Paper 191 found that q_M(1+x) has positive
        coefficients in tested Riemann-lattice cases.  This paper derives
        the exact coefficient formula and identifies the inequalities that
        would prove the barrier uniformly.
DEPENDS: 191 (coefficient criteria).

THEORY
======

Let a_i = 2*pi*m_i^2 and

    q_M(y) = sum_{r=0}^k (-1)^(k-r) c_{k-r} e_r(a) y^r,
    c_j = (2j+1)!!.

The coefficient of x^j in q_M(1+x) is

    B_j(M) = sum_{r=j}^k binom(r,j) (-1)^(k-r)
             c_{k-r} e_r(a).

So paper 191's certificate is exactly:

    B_j(M) > 0   for every j=0,...,k.

Equivalently, define the differential operator theta_j picking the j-th
Taylor coefficient at y=1:

    B_j = q_M^(j)(1)/j!.

The remaining theorem is no longer about roots; it is a finite family of
elementary symmetric inequalities on the separated lattice a_i=2*pi*m_i^2.
"""
from __future__ import annotations

import itertools
import math


def odd_double_factorial(n):
    out = 1.0
    for m in range(n, 0, -2):
        out *= m
    return out


def elementary_symmetric(xs, r):
    if r == 0:
        return 1.0
    return sum(math.prod(c) for c in itertools.combinations(xs, r))


def shifted_B(ms, j):
    k = len(ms)
    a = [2 * math.pi * m * m for m in ms]
    total = 0.0
    for r in range(j, k + 1):
        total += (math.comb(r, j)
                  * ((-1) ** (k - r))
                  * odd_double_factorial(2 * (k - r) + 1)
                  * elementary_symmetric(a, r))
    return total


def all_B(ms):
    return [shifted_B(ms, j) for j in range(len(ms) + 1)]


def main():
    print("=" * 80)
    print("PAPER 192 — shifted coefficient formula")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Coefficient profiles
    # ------------------------------------------------------------------
    print("\n[1] Coefficient profiles B_j for initial blocks")
    for k in [3, 5, 8, 12]:
        ms = tuple(range(1, k + 1))
        B = all_B(ms)
        ratios = [B[j] / B[j + 1] for j in range(len(B) - 1)]
        print(f"  k={k}, min B={min(B):.6e}, max B={max(B):.6e}")
        print("    B_j/B_{j+1}:", " ".join(f"{r:.3e}" for r in ratios[:8]))

    # ------------------------------------------------------------------
    # [2] Worst coefficient index in bounded search
    # ------------------------------------------------------------------
    print("\n[2] Worst normalized shifted coefficient")
    print("    Normalize by top coefficient B_k=e_k(a).")
    print(f"  {'M_max':>5} {'k':>3} {'min B_j/B_k':>16} {'j':>4} {'set':>18}")
    for M_max in [8, 10, 12, 14]:
        for k in range(2, min(8, M_max) + 1):
            best = (float("inf"), None, None)
            for ms in itertools.combinations(range(1, M_max + 1), k):
                B = all_B(ms)
                top = B[-1]
                for j, bj in enumerate(B):
                    val = bj / top
                    if val < best[0]:
                        best = (val, j, ms)
            print(f"  {M_max:>5} {k:>3} {best[0]:>16.6e}"
                  f" {best[1]:>4} {str(best[2]):>18}")

    # ------------------------------------------------------------------
    # [3] Sufficient dominance inequalities
    # ------------------------------------------------------------------
    print("\n[3] Alternating-tail dominance check")
    print("    For each j, compare positive last term to absolute sum of earlier negatives.")
    print(f"  {'k':>3} {'worst margin':>16} {'j':>4}")
    for k in range(2, 16):
        ms = tuple(range(1, k + 1))
        a = [2 * math.pi * m * m for m in ms]
        worst = float("inf")
        worst_j = None
        for j in range(k + 1):
            terms = []
            for r in range(j, k + 1):
                term = (math.comb(r, j)
                        * ((-1) ** (k - r))
                        * odd_double_factorial(2 * (k - r) + 1)
                        * elementary_symmetric(a, r))
                terms.append(term)
            pos = sum(t for t in terms if t > 0)
            neg = -sum(t for t in terms if t < 0)
            margin = (pos - neg) / pos if pos else -float("inf")
            if margin < worst:
                worst = margin
                worst_j = j
        print(f"  {k:>3} {worst:>16.6e} {worst_j:>4}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  REDUCTION:")
    print("    The Wronskian lattice barrier is implied by B_j(M)>0 for all j.")
    print("    Each B_j is an explicit alternating sum of elementary symmetric")
    print("    functions of 2*pi*m_i^2.")
    print()
    print("  NEXT THEOREM:")
    print("    Prove these alternating sums are positive using separation")
    print("    m_1^2 < ... < m_k^2.  The coefficient test is stronger and more")
    print("    elementary than root-location.")


if __name__ == "__main__":
    main()
