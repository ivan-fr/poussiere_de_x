"""
PAPER: 194 (NEW — one-parameter initial-block inequality)
TITLE: The final one-parameter barrier candidate:
       positivity of B_0({1,...,k})/e_k for all k
STATUS: RH remains OPEN.  Paper 193 suggests that all shifted coefficient
        inequalities reduce to the initial segment M={1,...,k}, j=0.
        This paper studies the resulting one-parameter sequence.
DEPENDS: 193 (initial segment extremality).

THEORY
======

Let a_m = 2*pi*m^2 and M_k={1,...,k}.  The candidate final inequality is

    R_k := B_0(M_k)/e_k(M_k) > 0.

Since e_r/e_k over M_k equals the elementary symmetric function of
reciprocals 1/a_m of degree k-r, we get

    R_k = sum_{s=0}^k (-1)^s (2s+1)!! E_s(k),

where

    E_s(k) = e_s( 1/(2*pi*1^2), ..., 1/(2*pi*k^2) ).

As k -> infinity, E_s(k) increases to e_s over the infinite sequence.
The limiting generating function is

    prod_{m>=1} (1 + z/(2*pi*m^2))
      = sinh(sqrt(pi*z/2)) / sqrt(pi*z/2).

So the limiting candidate is

    R_infty = sum_s (-1)^s (2s+1)!! [z^s] sinh(sqrt(pi*z/2))/sqrt(pi*z/2)
            = exp(-pi/4).

This paper computes R_k, estimates R_infty, and searches for a closed
form / positivity certificate.
"""
from __future__ import annotations

import math


def odd_double_factorial(n):
    out = 1.0
    for m in range(n, 0, -2):
        out *= m
    return out


def reciprocal_es(k):
    """Elementary symmetric coefficients of x_m=1/(2*pi*m^2)."""
    coeffs = [1.0] + [0.0] * k
    for m in range(1, k + 1):
        x = 1.0 / (2 * math.pi * m * m)
        for s in range(k, 0, -1):
            coeffs[s] += coeffs[s - 1] * x
    return coeffs


def R_k(k):
    E = reciprocal_es(k)
    return sum(((-1) ** s) * odd_double_factorial(2 * s + 1) * E[s]
               for s in range(k + 1))


def R_infty_trunc(S=80):
    """Use E_s(infty) from product sinh(sqrt(2*pi*z))/sqrt(2*pi*z).

    coefficient E_s = (pi/2)^s / (2s+1)!.
    """
    total = 0.0
    for s in range(S + 1):
        E_s = (math.pi / 2) ** s / math.factorial(2 * s + 1)
        total += ((-1) ** s) * odd_double_factorial(2 * s + 1) * E_s
    return total


def main():
    print("=" * 80)
    print("PAPER 194 — initial-block asymptotic")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] R_k values
    # ------------------------------------------------------------------
    print("\n[1] R_k = B_0({1,...,k})/e_k")
    print(f"  {'k':>4} {'R_k':>18} {'delta to limit est':>20}")
    Rinf = R_infty_trunc(80)
    for k in list(range(1, 21)) + [30, 40, 60, 80, 120]:
        rk = R_k(k)
        print(f"  {k:>4} {rk:>18.12f} {rk - Rinf:>20.12e}")
    print(f"  R_infty estimate = {Rinf:.12f}")

    # ------------------------------------------------------------------
    # [2] Closed-form guess
    # ------------------------------------------------------------------
    print("\n[2] Closed-form diagnostics")
    print("  R_infty series:")
    print("    sum_s (-1)^s (2s+1)!! (pi/2)^s/(2s+1)!")
    print("  Since (2s+1)!!/(2s+1)! = 1/(2^s s!),")
    print("    R_infty = sum_s (-1)^s (pi/4)^s/s! = exp(-pi/4).")
    print(f"  exp(-pi/4) = {math.exp(-math.pi / 4):.12f}")
    print(f"  difference estimate - exp(-pi/4) = {Rinf - math.exp(-math.pi / 4):.3e}")

    # ------------------------------------------------------------------
    # [3] Monotonicity check
    # ------------------------------------------------------------------
    print("\n[3] Monotonicity toward exp(-pi/4)")
    prev = None
    violations = 0
    for k in range(1, 150):
        rk = R_k(k)
        if prev is not None and rk > prev + 1e-14:
            violations += 1
        prev = rk
    print(f"  monotone decreasing violations over k<=149: {violations}")
    print(f"  lower target exp(-pi/4) positive: {math.exp(-math.pi / 4):.12f}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  MAJOR REDUCTION:")
    print("    If paper 193's extremal lemma is proved, the whole shifted")
    print("    coefficient barrier follows from R_k >= exp(-pi/4)>0.")
    print()
    print("  CLOSED FORM:")
    print("    The infinite limit is exactly exp(-pi/4), via Euler's sinh product.")
    print()
    print("  REMAINING THEOREM:")
    print("    Prove R_k decreases to exp(-pi/4) from above, and prove the")
    print("    initial-segment extremal lemma from paper 193.")


if __name__ == "__main__":
    main()
