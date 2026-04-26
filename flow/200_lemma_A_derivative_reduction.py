"""
PAPER: 200 (NEW -- derivative reduction for Lemma A)
TITLE: Reducing initial-segment extremality to a generalized lobe
       positivity derivative condition
STATUS: RH remains OPEN.

This paper attacks the open Lemma A from paper 193.

Let a_i = 2*pi*m_i^2 and x_i = 1/a_i.  After normalizing by e_k(a),
the shifted coefficient from paper 192 is

    F_{k,j}(x) := B_j(M)/e_k(a)
      = sum_{s=0}^{k-j} (-1)^s binom(k-s,j) (2s+1)!! e_s(x).

The initial-segment extremality conjecture says that, for every lattice
subset M of size k and every j,

    F_{k,j}(x_M) >= F_{k,0}(alpha_1,...,alpha_k),
    alpha_i = 1/(2*pi*i^2).

KEY OBSERVATION
===============

The binomial factor has a clean subset expansion:

    binom(k-s,j) e_s(x_1,...,x_k)
      = sum_{|J|=j} e_s(x_{[k]\\J}).

Therefore

    F_{k,j}(x) = sum_{|J|=j} G_{k-j}(x_{[k]\\J}),

where

    G_n(u) = sum_{s=0}^n (-1)^s (2s+1)!! e_s(u)
           = E[ Z^2 prod_i (1 - u_i Z^2) ].

Thus all j>0 coefficients are sums of B_0-type objects of smaller size.

DERIVATIVE REDUCTION
====================

For G_n,

    d/dx_i G_n(x) = - H_{n-1}(x without i),

where

    H_r(u) = sum_{s=0}^r (-1)^s (2s+3)!! e_s(u)
           = E[ Z^4 prod_i (1 - u_i Z^2) ].

So Lemma A follows from the generalized derivative positivity statement:

    H_r(x_M) >= 0 for every admissible Riemann-lattice subset M.

If H_r >= 0, then G_n decreases as each x_i increases.  Since every
ordered lattice subset satisfies m_i >= i, hence x_i <= alpha_i, the
minimum occurs at the initial segment x_i = alpha_i.  The j>0 case then
follows by the subset-sum identity and the already decreasing sequence
R_n = G_n(alpha_1,...,alpha_n).

This paper verifies the algebraic identities and stress-tests H_r >= 0.
It does NOT prove H_r >= 0.
"""
from __future__ import annotations

import itertools
import math
import random


def odd_double_factorial(n: int) -> float:
    out = 1.0
    for m in range(n, 0, -2):
        out *= m
    return out


def elementary_coeffs(xs: list[float]) -> list[float]:
    coeffs = [1.0] + [0.0] * len(xs)
    for x in xs:
        for r in range(len(xs), 0, -1):
            coeffs[r] += coeffs[r - 1] * x
    return coeffs


def reciprocal_xs(ms: tuple[int, ...]) -> list[float]:
    return [1.0 / (2.0 * math.pi * m * m) for m in ms]


def G_from_xs(xs: list[float]) -> float:
    es = elementary_coeffs(xs)
    return sum(
        (-1.0) ** s * odd_double_factorial(2 * s + 1) * es[s]
        for s in range(len(xs) + 1)
    )


def H_from_xs(xs: list[float]) -> float:
    es = elementary_coeffs(xs)
    return sum(
        (-1.0) ** s * odd_double_factorial(2 * s + 3) * es[s]
        for s in range(len(xs) + 1)
    )


def F_kj_from_xs(xs: list[float], j: int) -> float:
    k = len(xs)
    es = elementary_coeffs(xs)
    return sum(
        (-1.0) ** s
        * math.comb(k - s, j)
        * odd_double_factorial(2 * s + 1)
        * es[s]
        for s in range(k - j + 1)
    )


def F_kj_by_complements(xs: list[float], j: int) -> float:
    total = 0.0
    indices = range(len(xs))
    for keep_out in itertools.combinations(indices, j):
        keep_out_set = set(keep_out)
        complement = [x for i, x in enumerate(xs) if i not in keep_out_set]
        total += G_from_xs(complement)
    return total


def R_initial(k: int) -> float:
    return G_from_xs(reciprocal_xs(tuple(range(1, k + 1))))


def min_H_exhaustive(k: int, m_max: int) -> tuple[float, tuple[int, ...]]:
    best = (float("inf"), ())
    for ms in itertools.combinations(range(1, m_max + 1), k):
        val = H_from_xs(reciprocal_xs(ms))
        if val < best[0]:
            best = (val, ms)
    return best


def min_F_exhaustive(k: int, m_max: int) -> tuple[float, tuple[int, ...], int]:
    best = (float("inf"), (), -1)
    for ms in itertools.combinations(range(1, m_max + 1), k):
        xs = reciprocal_xs(ms)
        for j in range(k + 1):
            val = F_kj_from_xs(xs, j)
            if val < best[0]:
                best = (val, ms, j)
    return best


def main() -> None:
    print("=" * 80)
    print("PAPER 200 -- Lemma A derivative reduction")
    print("=" * 80)

    print("\n[1] Subset-sum identity check")
    print("    F_{k,j}(x) = sum_{|J|=j} G_{k-j}(x_{[k]\\J})")
    print(f"  {'k':>3} {'max abs residual':>18}")
    rng = random.Random(200)
    for k in range(2, 9):
        worst = 0.0
        for _ in range(40):
            ms = tuple(sorted(rng.sample(range(1, 5 * k + 1), k)))
            xs = reciprocal_xs(ms)
            for j in range(k + 1):
                direct = F_kj_from_xs(xs, j)
                expanded = F_kj_by_complements(xs, j)
                worst = max(worst, abs(direct - expanded))
        print(f"  {k:>3} {worst:>18.3e}")

    print("\n[2] Derivative identity check")
    print("    dG/dx_i = -H(x without i), verified by finite differences.")
    print(f"  {'k':>3} {'max abs residual':>18}")
    eps = 1.0e-7
    for k in range(2, 9):
        worst = 0.0
        for _ in range(30):
            ms = tuple(sorted(rng.sample(range(1, 5 * k + 1), k)))
            xs = reciprocal_xs(ms)
            for i in range(k):
                up = list(xs)
                down = list(xs)
                up[i] += eps
                down[i] -= eps
                finite = (G_from_xs(up) - G_from_xs(down)) / (2 * eps)
                exact = -H_from_xs(xs[:i] + xs[i + 1 :])
                worst = max(worst, abs(finite - exact))
        print(f"  {k:>3} {worst:>18.3e}")

    print("\n[3] Exhaustive H_r positivity stress test")
    print("    H_r(M) = E[Z^4 prod(1 - x_i Z^2)].")
    print(f"  {'M_max':>5} {'r':>3} {'min H_r':>16} {'witness':>22}")
    for m_max in [10, 12, 14]:
        for r in range(1, min(8, m_max) + 1):
            val, ms = min_H_exhaustive(r, m_max)
            print(f"  {m_max:>5} {r:>3} {val:>16.10f} {str(ms):>22}")

    print("\n[4] Lemma A full coefficient spot-check")
    print("    Minimum of F_{k,j}=B_j/e_k over subsets and j.")
    print(f"  {'M_max':>5} {'k':>3} {'min F':>16} {'witness':>22} {'j':>3} {'R_k':>16}")
    for m_max in [10, 12, 14]:
        for k in range(2, min(8, m_max) + 1):
            val, ms, j = min_F_exhaustive(k, m_max)
            print(
                f"  {m_max:>5} {k:>3} {val:>16.10f}"
                f" {str(ms):>22} {j:>3} {R_initial(k):>16.10f}"
            )

    print("\n[5] CONSEQUENCE IF H_r >= 0 IS PROVED")
    print("  1. G_n decreases in every reciprocal variable x_i.")
    print("  2. Since ordered lattice subsets satisfy m_i >= i, we have")
    print("     x_i <= 1/(2*pi*i^2).")
    print("  3. Therefore G_n(M) >= G_n({1,...,n}) = R_n.")
    print("  4. The subset-sum identity gives F_{k,j}(M) as a sum of such G terms.")
    print("  5. Since R_n decreases to exp(-pi/4), every F_{k,j}(M) >= R_k.")

    print("\n[6] HONEST ASSESSMENT")
    print("  NEW:")
    print("    Lemma A is reduced to a clean generalized lobe-positivity")
    print("    statement H_r(M) >= 0 for every lattice subset M.")
    print()
    print("  VERIFIED HERE:")
    print("    Algebraic subset identity, derivative identity, and exhaustive")
    print("    positivity checks up to the displayed boxes.")
    print()
    print("  STILL OPEN:")
    print("    A proof of H_r(M) >= 0 for all lattice subsets.")
    print("    This is close to paper 195's T_k positivity, but now for")
    print("    arbitrary subsets rather than only the initial segment.")


if __name__ == "__main__":
    main()
