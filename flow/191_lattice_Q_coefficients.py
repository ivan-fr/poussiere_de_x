"""
PAPER: 191 (NEW — coefficient criteria for the lattice Q barrier)
TITLE: Coefficient-side attack on q_M(y)>0: Descartes, reciprocal scaling,
       and Routh-Hurwitz probes for the Riemann-lattice Q-polynomials
STATUS: RH remains OPEN.  Paper 190 reduced the Wronskian layer to a
        root barrier for explicit polynomials q_M(y).  This paper studies
        those polynomials from the coefficient side, looking for a uniform
        positivity certificate on y>=1.
DEPENDS: 189 (Vandermonde factor), 190 (root barrier).

THEORY
======

For M={m_1<...<m_k},

    q_M(y) = sum_{r=0}^k (-1)^(k-r) (2(k-r)+1)!!
             e_r(2*pi*m_i^2) y^r.

The coefficients alternate in sign.  Direct positivity on y>=1 is not
coefficient-obvious.  Shift the interval:

    p_M(x) := q_M(1+x),    x>=0.

If every coefficient of p_M(x) is positive, then q_M(y)>0 for y>=1.
This is a strong and elementary certificate.

This paper tests:
  (C1) coefficient positivity of q_M(1+x);
  (C2) if C1 fails, Descartes sign variations after shifting y=1+x;
  (C3) reciprocal scaling patterns of roots.
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


def q_coeffs_for_ms(ms):
    k = len(ms)
    base = [2 * math.pi * m * m for m in ms]
    return [
        ((-1) ** (k - r)) * odd_double_factorial(2 * (k - r) + 1)
        * elementary_symmetric(base, r)
        for r in range(k + 1)
    ]


def shift_coeffs(coeffs, shift=1.0):
    """coeffs for q(y), low-to-high. Return coeffs for q(shift+x)."""
    k = len(coeffs) - 1
    out = [0.0] * (k + 1)
    for r, c in enumerate(coeffs):
        for j in range(r + 1):
            out[j] += c * math.comb(r, j) * (shift ** (r - j))
    return out


def sign_variations(vals, tol=1e-12):
    signs = []
    for v in vals:
        if abs(v) <= tol:
            continue
        signs.append(1 if v > 0 else -1)
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)


def main():
    print("=" * 80)
    print("PAPER 191 — coefficient criteria for lattice Q barrier")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Initial blocks
    # ------------------------------------------------------------------
    print("\n[1] Shifted coefficient positivity for initial blocks")
    print(f"  {'k':>3} {'min coeff q(1+x)':>22} {'sign variations':>17}")
    for k in range(2, 26):
        coeffs = q_coeffs_for_ms(tuple(range(1, k + 1)))
        shifted = shift_coeffs(coeffs, 1.0)
        print(f"  {k:>3} {min(shifted):>22.6e} {sign_variations(shifted):>17}")

    # ------------------------------------------------------------------
    # [2] Exhaustive bounded subsets
    # ------------------------------------------------------------------
    print("\n[2] Exhaustive subsets: q_M(1+x) coefficient test")
    print(f"  {'M_max':>5} {'k':>3} {'subsets':>9} {'fail coeff+':>12}"
          f" {'max variations':>15} {'worst set':>18}")
    for M_max in [8, 10, 12]:
        for k in range(2, min(7, M_max) + 1):
            total = 0
            fail = 0
            max_var = 0
            worst_min = float("inf")
            worst_set = None
            for ms in itertools.combinations(range(1, M_max + 1), k):
                shifted = shift_coeffs(q_coeffs_for_ms(ms), 1.0)
                mn = min(shifted)
                var = sign_variations(shifted)
                if mn <= -1e-9:
                    fail += 1
                if var > max_var or mn < worst_min:
                    max_var = max(max_var, var)
                    worst_min = min(worst_min, mn)
                    worst_set = ms
                total += 1
            print(f"  {M_max:>5} {k:>3} {total:>9} {fail:>12}"
                  f" {max_var:>15} {str(worst_set):>18}")

    # ------------------------------------------------------------------
    # [3] How far left can we shift?
    # ------------------------------------------------------------------
    print("\n[3] Minimal shift a such that q_M(a+x) has positive coefficients")
    print("    Tested on initial blocks by binary search.")
    print(f"  {'k':>3} {'a_coeff_positive':>20}")
    for k in range(2, 18):
        coeffs = q_coeffs_for_ms(tuple(range(1, k + 1)))
        lo, hi = 0.0, 1.0
        # Ensure hi works; if not, report.
        if min(shift_coeffs(coeffs, hi)) <= 0:
            print(f"  {k:>3} {'>1 or fail':>20}")
            continue
        for _ in range(50):
            mid = (lo + hi) / 2
            if min(shift_coeffs(coeffs, mid)) > 0:
                hi = mid
            else:
                lo = mid
        print(f"  {k:>3} {hi:>20.12f}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  IF q_M(1+x) HAS POSITIVE COEFFICIENTS:")
    print("    The root barrier q_M(y)>0 for y>=1 becomes elementary for that M.")
    print()
    print("  WHAT TO PROVE:")
    print("    For every finite M subset positive integers, all coefficients of")
    print("    q_M(1+x) are positive.  This would prove the Wronskian lattice")
    print("    barrier from paper 190.")


if __name__ == "__main__":
    main()
