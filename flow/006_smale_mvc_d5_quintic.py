"""
PAPER: 006 (canonical: 6pandrosion_smale.pdf)
TITLE: Toward the Quintic Case of Smale's Mean Value Conjecture --
       A Centered Pandrosion Slice and a Compactness Program
STATUS: proof candidate (not full proof for d=5)
DEPENDS: 005

THEORY
======

Setup: centered two-parameter quintic slice
  P(z) = z^5 + a*z^3 + b*z^2 + z,    (a, b) in C^2.
This slice contains the symmetric extremal P_0(z) = z^5 + z (a = b = 0).

Critical points: P'(c) = 5c^4 + 3a*c^2 + 2b*c + 1 = 0.

PANDROSION REDUCTION:
  tilde_q(c) = (1 - 3c^4 - a*c^2)/2     at each critical point c.

Smale ratio: S(c) = |tilde_q(c)|.

Smale's MVC for d=5:
  for all (a, b), min_{c critical} |tilde_q(c)| <= 4/5.

EXTREMAL POLYNOMIAL P_0(z) = z^5 + z (Theorem 2.1):
At (a, b) = (0, 0), P_0'(c) = 5c^4 + 1 = 0, so c^4 = -1/5, c^2 = ±i/sqrt(5).
  tilde_q(c) = (1 - 3 * (-1/5) - 0)/2 = (1 + 3/5)/2 = 4/5  (identically).

So all 4 critical points have |tilde_q(c)| = 4/5 (the conjectural maximum).

PRODUCT FORMULA (Theorem 1.1):
  prod_{j=1}^{4} tilde_q(c_j(a, b)) = N(a, b)/625
where N(a, b) = 16a^4 - 4a^3 b^2 - 128a^2 + 144ab^2 - 27b^4 + 256.
At extremal (0, 0): N(0, 0) = 256 = 4^4, prod = 256/625 = (4/5)^4.

THE COMPACTNESS PROGRAM (Section 5):
Three ingredients:
  (i) Extremal P_0 achieves |tilde_q| = 4/5 at all 4 critical points (rigorous).
  (ii) Local-descent certificate at (0, 0) (mixed symbolic/numerical).
  (iii) Asymptotic descent from compact parameter sets (partial).
Combining gives a proof candidate; gaps:
  - Mixed second-order cone estimate around delta a_re axis.
  - Discriminant-locus case requires Lojasiewicz analysis.
  - Reducing general d=5 quintics to the centered slice is separate.

VERIFICATION
============

This script verifies:
  1. The extremal P_0(z) = z^5 + z achieves |tilde_q| = 4/5 at all critical points.
  2. The product formula prod tilde_q = N(a, b)/625.
  3. Numerical verification: 1000 random (a, b), no violation of |tilde_q| <= 4/5.
"""
from __future__ import annotations
import math
import cmath
import numpy as np


def quintic_critical_points(a, b):
    """Roots of P'(c) = 5c^4 + 3a c^2 + 2b c + 1 = 0."""
    return np.roots([5, 0, 3*a, 2*b, 1])


def tilde_q_quintic(a, c):
    """tilde_q(c) = (1 - 3c^4 - a c^2) / 2 at critical c."""
    return (1 - 3*c**4 - a*c**2) / 2


def smale_ratio_quintic(a, b, c):
    """S(c) = |R(c)| where R(z) = z^4 + a z^2 + b z + 1."""
    return abs(c**4 + a*c**2 + b*c + 1)


def product_predicted(a, b):
    """N(a, b)/625 where N = 16a^4 - 4a^3 b^2 - 128a^2 + 144ab^2 - 27b^4 + 256."""
    return (16*a**4 - 4*a**3 * b**2 - 128*a**2 + 144*a*b**2 - 27*b**4 + 256) / 625


def main():
    print("=" * 80)
    print("PAPER 6 — Quintic case (d=5) of Smale's MVC")
    print("=" * 80)

    # 1. Extremal P_0(z) = z^5 + z (a=b=0): all critical S = 4/5
    print("\n[1] Extremal P_0(z) = z^5 + z: |tilde_q| = 4/5 at all critical points")
    a, b = 0.0, 0.0
    crits = quintic_critical_points(a, b)
    print(f"  Critical points: {[f'{c:.4f}' for c in crits]}")
    for c in crits:
        tq = tilde_q_quintic(a, c)
        S = smale_ratio_quintic(a, b, c)
        print(f"    c = {c:.4f}: |tilde_q(c)| = {abs(tq):.6f}, S = {S:.6f}")

    # 2. Product formula
    print("\n[2] Product formula: prod tilde_q(c_j) = N(a, b)/625")
    print(f"  At (0,0): N(0,0)/625 = 256/625 = (4/5)^4 = {(4/5)**4:.6f}")
    rng = np.random.default_rng(0)
    print(f"\n  {'a':>10} {'b':>10} {'product':>14} {'predicted':>14}")
    for _ in range(8):
        a = complex(rng.standard_normal(), rng.standard_normal())
        b = complex(rng.standard_normal(), rng.standard_normal())
        prod = 1.0 + 0j
        for c in quintic_critical_points(a, b):
            prod *= tilde_q_quintic(a, c)
        pred = product_predicted(a, b)
        print(f"  {a.real:>10.3f} {b.real:>10.3f} {prod.real:>14.4f} {pred.real:>14.4f}")

    # 3. Random sampling: min |tilde_q| <= 4/5 (Smale d=5)
    print("\n[3] 10000 random (a, b): is min |tilde_q(c)| always <= 4/5?")
    n_total = 10000
    violations = 0
    max_min = 0.0
    for _ in range(n_total):
        a = complex(rng.standard_normal()*2, rng.standard_normal()*2)
        b = complex(rng.standard_normal()*2, rng.standard_normal()*2)
        crits = quintic_critical_points(a, b)
        ratios = [abs(tilde_q_quintic(a, c)) for c in crits]
        min_r = min(ratios)
        if min_r > 0.8 + 1e-6: violations += 1
        if min_r > max_min: max_min = min_r
    print(f"  Tested {n_total} random (a, b): {violations} violations, max min = {max_min:.6f}")
    print(f"  Conjecture (open): min <= 4/5 = 0.800000")


if __name__ == "__main__":
    main()
