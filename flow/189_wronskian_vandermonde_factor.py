"""
PAPER: 189 (NEW — Wronskian Vandermonde factor)
TITLE: Vandermonde factorisation of the Phi-kernel Wronskians and the
       remaining positivity polynomial
STATUS: RH remains OPEN.  Paper 188 reduced higher sign-regularity to
        fixed-sign Wronskians.  This paper factors those Wronskians:
        the signed Wronskian equals a positive Vandermonde factor times
        an explicit symmetric polynomial Q_k.  The tempting theorem
        Q_k(t_1,...,t_k)>0 for all t_i >= 2*pi is FALSE on diagonals
        for some larger odd k.  The next theorem must exploit the
        spacing of the Riemann lattice lambda_i=pi i^2.
DEPENDS: 188 (Wronskian/Chebyshev route).

THEORY
======

Let
    g_lambda(y) = (2 lambda y - 3) exp(-lambda y).

After removing exp(-sum lambda_i y), the signed Wronskian for
lambda_1 < ... < lambda_k is

    SW_k = Vandermonde(lambda_1,...,lambda_k) * Q_k(t_1,...,t_k),
    t_i := 2 lambda_i y.

Symbolic computation gives

    Q_k(t) = sum_{r=0}^k (-1)^(k-r) (2(k-r)+1)!! e_r(t),

where e_r is the elementary symmetric polynomial.

Examples:
    Q_2 = e_2 - 3 e_1 + 15
    Q_3 = e_3 - 3 e_2 + 15 e_1 - 105
    Q_4 = e_4 - 3 e_3 + 15 e_2 - 105 e_1 + 945

in the t variables.  Since Riemann's range has

    t_i = 2 lambda_i y >= 2*pi > 6,

one might hope the lower bound alone gives Q_k>0.  The diagonal test
below disproves that.  The real target is the separated-lattice form

    Q_k(2*pi*m_1^2*y, ..., 2*pi*m_k^2*y) > 0

for 1 <= m_1 < ... < m_k and y >= 1.
"""
from __future__ import annotations

import itertools
import math
import random


def odd_double_factorial(n):
    out = 1
    for m in range(n, 0, -2):
        out *= m
    return out


def elementary_symmetric(xs, r):
    if r == 0:
        return 1.0
    return sum(math.prod(c) for c in itertools.combinations(xs, r))


def Q(xs):
    k = len(xs)
    total = 0.0
    for r in range(k + 1):
        total += ((-1) ** (k - r)) * odd_double_factorial(2 * (k - r) + 1) * elementary_symmetric(xs, r)
    return total


def Q_lower_bound_crude(xs):
    """A crude sufficient bound using only min(t_i).

    If all t_i >= T, then e_r(t) >= C(k,r) T^r.  This is not a proof of
    Q>0 because of alternating signs, but it gives a quick diagnostic.
    """
    k = len(xs)
    T = min(xs)
    return sum(((-1) ** (k - r)) * odd_double_factorial(2 * (k - r) + 1)
               * math.comb(k, r) * (T ** r)
               for r in range(k + 1))


def main():
    print("=" * 80)
    print("PAPER 189 — Wronskian Vandermonde factor")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Symbolic factorisation for k <= 5
    # ------------------------------------------------------------------
    print("\n[1] Symbolic Q_k pattern")
    try:
        import sympy as sp

        y = sp.symbols("y")
        for k in range(2, 6):
            ls = sp.symbols(f"l1:{k+1}")

            def br(lam, r):
                P = 2 * lam * y - 3
                if r == 0:
                    return P
                return sp.expand(P * ((-lam) ** r) + r * (2 * lam) * ((-lam) ** (r - 1)))

            M = sp.Matrix([[br(lam, r) for lam in ls] for r in range(k)])
            D = sp.factor(M.det())
            V = sp.prod(ls[j] - ls[i] for i in range(k) for j in range(i + 1, k))
            signed = (-1) ** (k * (k - 1) // 2)
            Q_lam = sp.factor(signed * D / V)
            # Convert to t variables mentally by t_i=2 lambda_i y.
            print(f"  k={k}: signed_det/V =")
            print(f"    {Q_lam}")
    except ImportError:
        print("  [sympy unavailable]")

    # ------------------------------------------------------------------
    # [2] Positivity of Q on sampled lower-bound range
    # ------------------------------------------------------------------
    print("\n[2] Numerical positivity of Q_k(t), sampled t_i >= 2*pi")
    rng = random.Random(2026)
    print(f"  {'k':>3} {'trials':>8} {'min Q':>16} {'min crude bound':>18}")
    for k in range(2, 11):
        min_q = float("inf")
        min_b = float("inf")
        trials = 3000
        for _ in range(trials):
            xs = sorted(2 * math.pi + 80 * rng.random() for _ in range(k))
            q = Q(xs)
            b = Q_lower_bound_crude(xs)
            min_q = min(min_q, q)
            min_b = min(min_b, b)
        print(f"  {k:>3} {trials:>8} {min_q:>16.6e} {min_b:>18.6e}")

    # ------------------------------------------------------------------
    # [2b] Riemann lattice positivity
    # ------------------------------------------------------------------
    print("\n[2b] Riemann lattice Q_k(2*pi*m_i^2*y)")
    print(f"  {'k':>3} {'min over subsets/y':>22} {'case':>24}")
    y_grid = [1.0, 1.01, 1.05, 1.2, 1.7, 2.5, 4.0]
    for k in range(2, 11):
        best = (float("inf"), None)
        for ms in itertools.combinations(range(1, 13), k):
            for y in y_grid:
                xs = [2 * math.pi * m * m * y for m in ms]
                q = Q(xs)
                if q < best[0]:
                    best = (q, (ms, y))
        print(f"  {k:>3} {best[0]:>22.6e} {str(best[1]):>24}")

    # ------------------------------------------------------------------
    # [3] Diagonal counter-test
    # ------------------------------------------------------------------
    print("\n[3] Diagonal counter-test Q_k(T,...,T)")
    print(f"  {'k':>3} {'Q_k(2pi)':>18} {'first integer T with Q>0':>28}")
    for k in range(2, 15):
        q2 = Q([2 * math.pi] * k)
        first = None
        for T in range(1, 30):
            if Q([float(T)] * k) > 0:
                first = T
                break
        print(f"  {k:>3} {q2:>18.6e} {str(first):>28}")

    # ------------------------------------------------------------------
    # [4] Probabilistic representation
    # ------------------------------------------------------------------
    print("\n[4] Probabilistic representation")
    print("  Since E[Z^(2m)] = (2m-1)!! for Z~N(0,1),")
    print("  Q_k(t) = E[ Z^2 * prod_i (t_i - Z^2) ].")
    print("  This representation is exact but not immediately positive,")
    print("  because the integrand changes sign when Z^2 exceeds some t_i.")

    # ------------------------------------------------------------------
    # [5] Honest assessment
    # ------------------------------------------------------------------
    print("\n[5] HONEST ASSESSMENT")
    print("  PROGRESS:")
    print("    Higher Wronskians reduce to Vandermonde positivity times Q_k.")
    print("    The Riemann lattice t_i=2*pi*m_i^2*y is positive in tested ranges.")
    print()
    print("  IMPORTANT CORRECTION:")
    print("    The simpler theorem 'Q_k>0 for all t_i>=2pi' is false:")
    print("    diagonal tests are negative for some larger odd k.  Separation")
    print("    of the lambda_i matters.")
    print()
    print("  REMAINING THEOREM:")
    print("    Prove Q_k(2*pi*m_1^2*y,...,2*pi*m_k^2*y)>0 for all k,")
    print("    y>=1, and 1<=m_1<...<m_k.")
    print("    That would give the extended Chebyshev part of the Phi-kernel")
    print("    sign-regularity program.")


if __name__ == "__main__":
    main()
