"""
PAPER: 187 (NEW — symbolic 2x2 sign-regularity)
TITLE: Symbolic 2x2 positivity for the Riemann-Phi summand kernel
       K(lambda,y)=y^(5/4)(2 lambda y - 3) exp(-lambda y)
STATUS: RH remains OPEN.  This paper proves/checks the first nontrivial
        sign-regular layer for the kernel isolated in paper 186: all 2x2
        signed minors are positive for lambda_1 < lambda_2 and
        1 <= y_1 < y_2.  Higher minors remain open.
DEPENDS: 186 (exponential-kernel structure).

THEORY
======

The positive column factor y^(5/4) does not affect determinant signs, so
consider

    L(lambda,y) = (2 lambda y - 3) exp(-lambda y).

For a<b and x<y, the signed 2x2 determinant is

    D = - det [[L(a,x), L(a,y)],
               [L(b,x), L(b,y)]].

Factoring exp(-a y - b x) gives

    D = (2 a y - 3)(2 b x - 3)
        - exp(-(b-a)(y-x)) (2 a x - 3)(2 b y - 3).

For the Riemann range a >= pi, x >= 1, every factor (2 lambda y - 3)
is positive.  The target inequality is therefore

    R := ((2 a y - 3)(2 b x - 3)) /
         ((2 a x - 3)(2 b y - 3))
       > exp(-(b-a)(y-x)).

This paper verifies the stronger monotonicity mechanism numerically and
symbolically differentiates log R against the rectangle variables.  The
2x2 layer is a real proof target; size >= 3 needs a separate theorem.
"""
from __future__ import annotations

import itertools
import math


def signed_2x2_expr(a, b, x, y):
    """Signed 2x2 determinant after removing common positive exponent."""
    return ((2 * a * y - 3) * (2 * b * x - 3)
            - math.exp(-(b - a) * (y - x))
            * (2 * a * x - 3) * (2 * b * y - 3))


def log_ratio_minus_exp_bound(a, b, x, y):
    R = ((2 * a * y - 3) * (2 * b * x - 3)) / (
        (2 * a * x - 3) * (2 * b * y - 3)
    )
    return math.log(R) + (b - a) * (y - x)


def full_signed_det(a, b, x, y):
    def K(lam, yy):
        return (yy ** 1.25) * (2 * lam * yy - 3) * math.exp(-lam * yy)

    det = K(a, x) * K(b, y) - K(a, y) * K(b, x)
    return -det


def main():
    print("=" * 80)
    print("PAPER 187 — symbolic 2x2 positivity for Phi kernel")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Sympy symbolic derivatives of log R
    # ------------------------------------------------------------------
    print("\n[1] Symbolic derivative structure")
    try:
        import sympy as sp

        a, b, x, y = sp.symbols("a b x y", positive=True)
        logR = (sp.log(2 * a * y - 3) + sp.log(2 * b * x - 3)
                - sp.log(2 * a * x - 3) - sp.log(2 * b * y - 3))
        d_y = sp.factor(sp.diff(logR, y))
        d_x = sp.factor(sp.diff(logR, x))
        print("  d/dy log R =")
        print(f"    {d_y}")
        print("  d/dx log R =")
        print(f"    {d_x}")
        print("  Therefore for H=log(R)+(b-a)(y-x):")
        print("    H(x,x)=0 and dH/dy > 0 for b>a, y>=x>=1.")
        print("    Hence H>0 when y>x, proving the 2x2 factored inequality.")
    except ImportError:
        print("  [sympy unavailable]")

    # ------------------------------------------------------------------
    # [2] Grid verification of the factored inequality
    # ------------------------------------------------------------------
    print("\n[2] Factored inequality check on Riemann range")
    lambdas = [math.pi * m * m for m in range(1, 10)]
    ys = [1.0, 1.01, 1.05, 1.15, 1.4, 2.0, 3.5, 5.0]
    min_factored = float("inf")
    min_full = float("inf")
    min_case = None
    violations = 0
    total = 0
    for i, j in itertools.combinations(range(len(lambdas)), 2):
        a = lambdas[i]
        b = lambdas[j]
        for p, q in itertools.combinations(range(len(ys)), 2):
            x = ys[p]
            y = ys[q]
            factored = signed_2x2_expr(a, b, x, y)
            margin = log_ratio_minus_exp_bound(a, b, x, y)
            full = full_signed_det(a, b, x, y)
            total += 1
            # Use the factored/logarithmic form for correctness.  The full
            # determinant underflows for large lambda*y in double precision.
            if factored <= 0 or margin <= 0:
                violations += 1
            if margin < min_factored:
                min_factored = margin
                min_full = full
                min_case = (a, b, x, y, factored)
    print(f"  rectangles tested: {total}")
    print(f"  violations: {violations}")
    print(f"  minimum log-margin: {min_factored:.6e}")
    print(f"  corresponding full signed determinant: {min_full:.6e}")
    print(f"  min case (a,b,x,y,factored): {min_case}")

    # ------------------------------------------------------------------
    # [3] Near-boundary asymptotic
    # ------------------------------------------------------------------
    print("\n[3] Near-boundary behaviour y=x+h")
    a = math.pi
    b = 4 * math.pi
    x = 1.0
    print(f"  fixed a=pi, b=4pi, x=1")
    print(f"  {'h':>10} {'signed det':>16} {'det / h':>16}")
    for h in [1e-1, 1e-2, 1e-3, 1e-4]:
        d = full_signed_det(a, b, x, x + h)
        print(f"  {h:>10.1e} {d:>16.8e} {d / h:>16.8e}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  PROGRESS:")
    print("    The 2x2 layer has a clean factored inequality.  This is the")
    print("    first analytic foothold for the kernel sign-regularity program.")
    print()
    print("  STILL OPEN:")
    print("    Size >= 3 minors.  The 2x2 proof does not imply total positivity.")
    print()
    print("  NEXT MOVE:")
    print("    Try to represent K(lambda,y) as a composition of known")
    print("    sign-regular kernels or derive a Wronskian criterion for the")
    print("    family y -> K(lambda_i,y).")


if __name__ == "__main__":
    main()
