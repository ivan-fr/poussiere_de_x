"""
PAPER: 201 (NEW -- anchor bound candidate for generalized H positivity)
TITLE: Evidence that the generalized H_r barrier is minimized by a single
       active lattice point m=1 and all other points sent to infinity
STATUS: RH remains OPEN.

Paper 200 reduced Lemma A to the generalized derivative positivity

    H_M := E[ Z^4 prod_{m in M} (1 - Z^2/(2*pi*m^2)) ] >= 0.

This paper probes the sharper candidate bound

    H_M >= H_{ {1} } = 3 - 15/(2*pi) = 0.6126758536...

for every finite Riemann-lattice subset M.

If true, the derivative condition in paper 200 is closed with a clean
margin, and Lemma A follows.

The key empirical surprise is that H_M is NOT minimized by the initial
segment.  It is minimized, in finite boxes, by keeping m=1 and pushing all
other selected lattice points as far right as possible:

    M = {1, M_max-r+2, ..., M_max}.

That points to an "anchor + harmless tail" theorem rather than the same
initial-segment extremality used for R_k.
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


def moment_product(xs: list[float], base_power: int) -> float:
    """E[Z^base_power prod(1 - x_i Z^2)] for even base_power."""
    es = elementary_coeffs(xs)
    return sum(
        (-1.0) ** s
        * odd_double_factorial(base_power + 2 * s - 1)
        * es[s]
        for s in range(len(xs) + 1)
    )


def H(ms: tuple[int, ...]) -> float:
    return moment_product(reciprocal_xs(ms), 4)


def derivative_H_wrt_member(ms: tuple[int, ...], index: int) -> float:
    """Derivative of H with respect to x_index, not m_index."""
    xs = reciprocal_xs(ms)
    rest = xs[:index] + xs[index + 1 :]
    return -moment_product(rest, 6)


def exhaustive_min(r: int, m_max: int) -> tuple[float, tuple[int, ...]]:
    best = (float("inf"), ())
    for ms in itertools.combinations(range(1, m_max + 1), r):
        val = H(ms)
        if val < best[0]:
            best = (val, ms)
    return best


def random_min(r: int, m_max: int, trials: int, rng: random.Random) -> tuple[float, tuple[int, ...]]:
    best = (float("inf"), ())
    for _ in range(trials):
        ms = tuple(sorted(rng.sample(range(1, m_max + 1), r)))
        val = H(ms)
        if val < best[0]:
            best = (val, ms)
    return best


def right_tail_candidate(r: int, m_max: int) -> tuple[int, ...]:
    if r == 1:
        return (1,)
    return (1,) + tuple(range(m_max - r + 2, m_max + 1))


def main() -> None:
    print("=" * 80)
    print("PAPER 201 -- generalized H anchor bound")
    print("=" * 80)

    h_anchor = H((1,))
    print("\n[1] Proposed universal lower bound")
    print(f"  H_{{1}} = 3 - 15/(2*pi) = {h_anchor:.12f}")
    print(f"  positive margin = {h_anchor:.12f}")

    print("\n[2] Single-point H values")
    print(f"  {'m':>4} {'H_{m}':>16}")
    for m in [1, 2, 3, 4, 5, 10, 20, 50]:
        print(f"  {m:>4} {H((m,)):>16.12f}")

    print("\n[3] Exhaustive finite-box minima")
    print(f"  {'M_max':>5} {'r':>3} {'min H':>16} {'minus H1':>14} {'witness':>28}")
    for m_max in [10, 12, 14, 16]:
        for r in range(1, min(8, m_max) + 1):
            val, ms = exhaustive_min(r, m_max)
            print(
                f"  {m_max:>5} {r:>3} {val:>16.12f}"
                f" {val - h_anchor:>14.3e} {str(ms):>28}"
            )

    print("\n[4] Right-tail candidate convergence")
    print("    Candidate: {1, M_max-r+2, ..., M_max}.")
    print(f"  {'M_max':>6} {'r':>3} {'H(candidate)':>16} {'minus H1':>14}")
    for m_max in [20, 50, 100, 250, 1000]:
        for r in [2, 4, 8, 12]:
            ms = right_tail_candidate(r, m_max)
            val = H(ms)
            print(f"  {m_max:>6} {r:>3} {val:>16.12f} {val - h_anchor:>14.3e}")

    print("\n[5] Random large-box search")
    rng = random.Random(201)
    print(f"  {'M_max':>6} {'r':>3} {'trials':>8} {'min H':>16} {'minus H1':>14} {'witness':>24}")
    for m_max, r, trials in [(100, 8, 2000), (100, 14, 2000), (500, 8, 2000), (500, 20, 2000)]:
        val, ms = random_min(r, m_max, trials, rng)
        print(
            f"  {m_max:>6} {r:>3} {trials:>8} {val:>16.12f}"
            f" {val - h_anchor:>14.3e} {str(ms[:6]) + ('...' if len(ms) > 6 else ''):>24}"
        )

    print("\n[6] Derivative signs at right-tail candidates")
    print("    The anchor derivative is negative: the minimum wants x_1 as large")
    print("    as allowed, i.e. m=1.  Tail derivatives are positive: the minimum")
    print("    wants tail reciprocal variables smaller, i.e. farther right.")
    print(f"  {'M':>32} {'anchor dH':>14} {'min tail dH':>16} {'max tail dH':>16}")
    for r, m_max in [(3, 50), (5, 50), (8, 100), (12, 250)]:
        ms = right_tail_candidate(r, m_max)
        derivs = [derivative_H_wrt_member(ms, i) for i in range(len(ms))]
        tail_derivs = derivs[1:]
        print(
            f"  {str(ms):>32} {derivs[0]:>14.10f}"
            f" {min(tail_derivs):>16.10f} {max(tail_derivs):>16.10f}"
        )

    print("\n[7] HONEST ASSESSMENT")
    print("  NEW CANDIDATE THEOREM:")
    print("    H_M >= H_{1} = 3 - 15/(2*pi) > 0.")
    print()
    print("  EVIDENCE:")
    print("    Exhaustive finite boxes and random large boxes find no violations.")
    print("    The minimizing pattern is not initial segment; it is one anchor")
    print("    m=1 plus a tail sent toward infinity.")
    print()
    print("  STILL OPEN:")
    print("    Prove the anchor-tail inequality.  A likely route is to split")
    print("    cases m_1=1 and m_1>=2, then show the tail contribution is")
    print("    nonnegative in the first case and harmlessly small in the second.")


if __name__ == "__main__":
    main()
