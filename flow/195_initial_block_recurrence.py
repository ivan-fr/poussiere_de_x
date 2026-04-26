"""
PAPER: 195 (NEW — recurrence for the initial-block barrier)
TITLE: Recurrence and Gaussian-product form for
       R_k = B_0({1,...,k})/e_k({1,...,k})
STATUS: RH remains OPEN.  Paper 194 reduced one branch of the Wronskian
        barrier to R_k > 0 and observed R_k decreases to exp(-pi/4).
        This paper derives a recurrence/expectation form and tests the
        positivity of the increments R_k - R_{k+1}.
DEPENDS: 194 (initial-block asymptotic).

THEORY
======

Let
    x_m = 1/(2*pi*m^2),
    P_k(z) = prod_{m=1}^k (1 - x_m z).

Then

    R_k = sum_s (-1)^s (2s+1)!! e_s(x_1,...,x_k)
        = E[ Z^2 P_k(Z^2) ],        Z ~ N(0,1).

The infinite product is

    P_infty(z) = prod_m (1 - z/(2*pi*m^2))
               = sin(sqrt(pi*z/2)) / sqrt(pi*z/2),

so

    R_infty = E[ Z^2 sin(sqrt(pi/2)Z)/(sqrt(pi/2)Z) ] = exp(-pi/4).

The finite-step difference is

    R_k - R_{k+1}
      = x_{k+1} E[ Z^4 P_k(Z^2) ].

Therefore monotone decrease R_k >= R_{k+1} follows from

    T_k := E[ Z^4 P_k(Z^2) ] >= 0.

This paper tests T_k and nearby analogues.  Higher analogues are not all
positive; the theorem needed for monotonicity is specifically T_k >= 0.
"""
from __future__ import annotations

import math


def odd_double_factorial(n):
    out = 1.0
    for m in range(n, 0, -2):
        out *= m
    return out


def product_coeffs_x(k):
    """P_k(z)=prod(1-x_m z), low-to-high coefficients."""
    coeffs = [1.0]
    for m in range(1, k + 1):
        x = 1.0 / (2 * math.pi * m * m)
        coeffs.append(0.0)
        for j in range(len(coeffs) - 1, 0, -1):
            coeffs[j] -= x * coeffs[j - 1]
    return coeffs


def gaussian_even_moment(power):
    """E[Z^power] for even power."""
    if power % 2:
        return 0.0
    if power == 0:
        return 1.0
    return odd_double_factorial(power - 1)


def gaussian_product_expectation(k, z_power):
    """E[Z^z_power P_k(Z^2)]."""
    coeffs = product_coeffs_x(k)
    return sum(c * gaussian_even_moment(z_power + 2 * j)
               for j, c in enumerate(coeffs))


def R_k(k):
    return gaussian_product_expectation(k, 2)


def main():
    print("=" * 80)
    print("PAPER 195 — recurrence for initial-block barrier")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Recurrence check
    # ------------------------------------------------------------------
    print("\n[1] Recurrence R_k - R_{k+1} = x_{k+1} T_k")
    print(f"  {'k':>4} {'R_k':>16} {'R_k-R_{k+1}':>16}"
          f" {'x*T_k':>16} {'T_k':>16}")
    for k in list(range(1, 16)) + [25, 40, 80]:
        rk = R_k(k)
        diff = rk - R_k(k + 1)
        xnext = 1.0 / (2 * math.pi * (k + 1) ** 2)
        Tk = gaussian_product_expectation(k, 4)
        print(f"  {k:>4} {rk:>16.10f} {diff:>16.8e}"
              f" {xnext * Tk:>16.8e} {Tk:>16.8e}")

    # ------------------------------------------------------------------
    # [2] Positivity of T_k and higher moments
    # ------------------------------------------------------------------
    print("\n[2] Gaussian-product expectations E[Z^(2r) P_k(Z^2)]")
    print("    r=2 gives monotonicity of R_k; higher r need not be positive.")
    print(f"  {'k max':>6} {'r':>3} {'min value':>16} {'arg k':>8}")
    for r in range(1, 7):
        min_val = float("inf")
        arg = None
        for k in range(1, 120):
            val = gaussian_product_expectation(k, 2 * r)
            if val < min_val:
                min_val = val
                arg = k
        print(f"  {119:>6} {r:>3} {min_val:>16.8e} {arg:>8}")

    # ------------------------------------------------------------------
    # [3] Scaling of monotone tail
    # ------------------------------------------------------------------
    print("\n[3] Tail scaling R_k - exp(-pi/4)")
    target = math.exp(-math.pi / 4)
    print(f"  {'k':>5} {'R_k-target':>16} {'k*(gap)':>16}")
    for k in [10, 20, 40, 80, 120]:
        gap = R_k(k) - target
        print(f"  {k:>5} {gap:>16.8e} {k * gap:>16.8e}")

    # ------------------------------------------------------------------
    # [3b] Limit of T_k
    # ------------------------------------------------------------------
    print("\n[3b] Limit of T_k")
    c2 = math.pi / 2
    T_limit = (3 - c2) * math.exp(-c2 / 2)
    print("  P_infty(Z^2)=sin(sqrt(pi/2)Z)/(sqrt(pi/2)Z)")
    print("  lim T_k = E[Z^4 P_infty(Z^2)]")
    print("          = (3 - pi/2) exp(-pi/4)")
    print(f"          = {T_limit:.12f}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  NEW REDUCTION:")
    print("    R_k monotonicity is equivalent to T_k=E[Z^4 P_k(Z^2)] >= 0.")
    print()
    print("  OBSERVED:")
    print("    T_k is positive in the tested range and tends to")
    print("    (3 - pi/2) exp(-pi/4)>0.  Higher moments are not all positive,")
    print("    so the monotonicity theorem must target T_k specifically.")
    print()
    print("  NEXT THEOREM:")
    print("    Prove E[Z^4 prod_{m<=k}(1-Z^2/(2*pi*m^2))] >= 0 for all k.")
    print("    This would prove R_k decreases to exp(-pi/4)>0.")


if __name__ == "__main__":
    main()
