"""
PAPER: 005 (canonical: 5pandrosion_smale.pdf)
TITLE: The Quartic Case of Smale's Conjecture in Pandrosion Form --
       A Reduction to a Resultant Inequality
STATUS: d=4 case is classical (Tischler 1989); this paper records the
        Pandrosion-form reduction.
DEPENDS: 004

THEORY
======

Setup: P(z) = z^4 + a*z^3 + b*z^2 + z, monic quartic with P(0) = 0, P'(0) = 1.
Companion: R(z) = z^3 + a*z^2 + b*z + 1.
Smale ratio: S(c) = |R(c)| = |P(c)/c|.

Critical points satisfy P'(c) = 4c^3 + 3a*c^2 + 2b*c + 1 = 0.

PANDROSION REDUCTION (Theorem 2.1):
At each critical point c:
  R(c) = tilde_q(c)/4,    tilde_q(c) = a*c^2 + 2b*c + 3.

CRITICAL SIMPLIFICATION (Theorem 2.2):
  tilde_q(c) = 2(1 - c^2(a + 2c)).

COROLLARY 2.3 (Smale d=4 reformulation):
For all (a, b) in C^2:
  exists c with 4c^3 + 3a*c^2 + 2b*c + 1 = 0 AND |1 - c^2(a + 2c)| <= 3/2.

EXTREMAL POLYNOMIAL (Theorem 3.1):
For P_0(z) = z^4 - z (a = b = 0): all 3 critical points satisfy S(c) = 3/4
EXACTLY. Critical equation: 4c^3 + 1 = 0, c^3 = -1/4. Then tilde_q(c) =
2(1 - 2c^3) = 2(1 + 1/2) = 3. S(c) = 3/4.

PRODUCT FORMULA (Theorem 4.1):
  prod_{k=1}^{3} tilde_q(c_k)  =  4a^3 - a^2 b^2 - 18 ab + 4 b^3 + 27.

This is the discriminant of the depressed cubic 4z^4 + 3az^2 + 2bz + 1
related to substitution z -> z - a/4 of the critical cubic.
At extremal (0, 0): prod tilde_q = 27 = 3^3.

VERIFICATION
============

This script verifies:
  1. The critical simplification tilde_q(c) = 2(1 - c^2(a + 2c)) at critical c.
  2. The extremal P_0(z) = z^4 - z achieves min |tilde_q| = 4/3 = 3/4 * 4.
  3. The product formula prod tilde_q(c_k) = N(a, b) = 4a^3 - a^2 b^2 - 18 ab + 4 b^3 + 27.
  4. Random samples: min over critical points of |tilde_q/4| <= 3/4.
"""
from __future__ import annotations
import math
import cmath
import numpy as np


def quartic_critical_points(a, b):
    """Roots of P'(c) = 4c^3 + 3 a c^2 + 2 b c + 1 = 0."""
    return np.roots([4, 3*a, 2*b, 1])


def tilde_q_quartic(a, b, c):
    """tilde_q(c) = a c^2 + 2 b c + 3 (raw form)."""
    return a*c**2 + 2*b*c + 3


def tilde_q_critical_simplified(a, c):
    """At critical c: tilde_q(c) = 2(1 - c^2(a + 2c))."""
    return 2 * (1 - c**2 * (a + 2*c))


def smale_ratio_quartic(a, b, c):
    """S(c) = |R(c)| = |c^3 + a c^2 + b c + 1|."""
    return abs(c**3 + a*c**2 + b*c + 1)


def product_formula_predicted(a, b):
    """N(a, b) = 4a^3 - a^2 b^2 - 18 a b + 4 b^3 + 27."""
    return 4*a**3 - a**2 * b**2 - 18*a*b + 4*b**3 + 27


def main():
    print("=" * 80)
    print("PAPER 5 — Quartic case (d=4) of Smale's MVC, Pandrosion form")
    print("=" * 80)

    # 1. Critical simplification
    print("\n[1] Critical simplification: tilde_q(c) = 2(1 - c^2(a + 2c)) at critical c")
    rng = np.random.default_rng(0)
    max_err = 0.0
    for _ in range(20):
        a = complex(rng.standard_normal(), rng.standard_normal())
        b = complex(rng.standard_normal(), rng.standard_normal())
        for c in quartic_critical_points(a, b):
            tq_raw = tilde_q_quartic(a, b, c)
            tq_simp = tilde_q_critical_simplified(a, c)
            err = abs(tq_raw - tq_simp)
            if err > max_err: max_err = err
    print(f"  max |raw - simplified| = {max_err:.2e} (expected ~0)")

    # 2. Extremal polynomial P_0 = z^4 - z (a=b=0): all critical S = 3/4
    print("\n[2] Extremal P_0(z) = z^4 - z (a=b=0): S(c) = 3/4 at all critical points")
    a, b = 0.0, 0.0
    crits = quartic_critical_points(a, b)
    print(f"  Critical points: {crits}")
    for c in crits:
        # For a=b=0: critical eq 4c^3 + 1 = 0, so c^3 = -1/4
        # tilde_q(c) = 2(1 - 2c^3) = 2(1 + 1/2) = 3
        # S(c) = |R(c)| = |c^3 + 1| = |-1/4 + 1| = 3/4
        S = smale_ratio_quartic(a, b, c)
        tq = tilde_q_critical_simplified(a, c)
        ratio = abs(tq) / 4
        print(f"    c = {c:.6f}, S(c) = {S:.6f}, |tilde_q|/4 = {ratio:.6f} "
              f"(target 3/4 = {0.75})")

    # 3. Product formula
    print("\n[3] Product formula: prod tilde_q(c_k) = 4a^3 - a^2 b^2 - 18 ab + 4 b^3 + 27")
    print(f"  {'a':>10} {'b':>10} {'product':>14} {'predicted':>14} {'diff':>10}")
    for _ in range(8):
        a = complex(rng.standard_normal(), rng.standard_normal())
        b = complex(rng.standard_normal(), rng.standard_normal())
        prod = 1.0 + 0j
        for c in quartic_critical_points(a, b):
            prod *= tilde_q_critical_simplified(a, c)
        pred = product_formula_predicted(a, b)
        print(f"  {a.real:>10.3f} {b.real:>10.3f} {prod.real:>14.4f} "
              f"{pred.real:>14.4f} {abs(prod - pred):>10.2e}")

    # 4. Random sampling: min |tilde_q/4| <= 3/4 (Smale d=4)
    print("\n[4] Random samples: min over critical points of |tilde_q/4| <= 3/4")
    n_total = 1000
    violations = 0
    max_min = 0.0
    for _ in range(n_total):
        a = complex(rng.standard_normal()*2, rng.standard_normal()*2)
        b = complex(rng.standard_normal()*2, rng.standard_normal()*2)
        ratios = [abs(tilde_q_critical_simplified(a, c))/4 for c in quartic_critical_points(a, b)]
        min_r = min(ratios)
        if min_r > 0.75 + 1e-9: violations += 1
        if min_r > max_min: max_min = min_r
    print(f"  {n_total} random (a, b): {violations} violations, max min ratio = {max_min:.6f}")


if __name__ == "__main__":
    main()
