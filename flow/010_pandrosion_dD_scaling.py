"""
PAPER: 010 (paper-0 sec:scaling applied to Pandrosion dD)
TITLE: Scaling Principle in Higher Dimensions
STATUS: empirical study; theoretical asymptotic verified

THEORY
======

Paper-0 sec:scaling for x^{1/p}:
    x^{1/p} = A^{1/p} * (x/A)^{1/p},  A = floor(x^{1/p})^p
gives x' = x/A in [1, 2), and lambda_{p, x'} << lambda_{p, x} for large x.

For Pandrosion dD with p=2, the contraction rate at the symmetric fixed
point is (Theorem 2 of paper II):
    lambda_dD(x) = x(d-1)(1-s_*)^2 / (x-1)

The scaling argument shows that for x large, lambda_dD(x) -> 1 as in 2D.
For x of order 1, lambda_dD(x) is small (already accelerated by dimension).

ASYMPTOTIC RESULT (this file): for arbitrary x and large d,
    1 - s_*(x, d) = (x - 1) / (d x) + O(1/d^2)
    lambda_dD(x) = (d-1)(x-1) / (d^2 x) ~ (1 - 1/x) / d  as d -> infinity.

CONSEQUENCE: increasing d makes lambda_dD ~ 1/d INDEPENDENTLY OF x.
The dD dimension absorbs the scaling problem!  Pandrosion dD is naturally
fast for any x, large or small, without the integer-A scaling trick.

VERIFICATION
============
Compare:
  - Pandrosion 2D raw on large x (slow, lambda close to 1)
  - Pandrosion 2D + paper-0 scaling (fast, lambda close to 0)
  - Pandrosion dD raw on large x (fast WITHOUT scaling)
  - Pandrosion dD + scaling (overkill?)
for x in {10, 100, 10000, 10^8} and d in {2, 5, 32}.
"""
from __future__ import annotations
import math
from math import comb


def Delta_d(s, p, d):
    return sum(comb(m + d - 2, d - 2) * s**m for m in range(p))


def fixed_point_p2(x, d):
    """Closed form: (d-1) x s^2 - (d-2) x s - 1 = 0."""
    if d == 2:  # paper-0 scalar 2D
        return 1.0 / math.sqrt(x)
    a = (d - 1) * x
    b = -(d - 2) * x
    c = -1.0
    disc = b * b - 4 * a * c
    return (-b + math.sqrt(disc)) / (2 * a)


def lambda_dD(x, d):
    """Closed-form contraction rate at symmetric fixed point, p=2."""
    s_star = fixed_point_p2(x, d)
    if d == 2:
        # paper-0 strict: lambda_{2,x} = (x-1)/(x(2+sqrt(x)))^? Let me rederive.
        # 2D Pandrosion p=2: g(s) = 1 - (x-1)/(x(1+s)).  g'(s) = (x-1)/(x(1+s)^2).
        # At s=1/sqrt(x): (1+s)^2 = (1+1/sqrt(x))^2 = (sqrt(x)+1)^2/x.
        # lambda = (x-1)/x * x/(sqrt(x)+1)^2 = (x-1)/(sqrt(x)+1)^2.
        return (x - 1) / (math.sqrt(x) + 1) ** 2
    return x * (d - 1) * (1 - s_star) ** 2 / (x - 1)


def n_iters_to_eps(lam, eps_target=1e-15, eps_start=0.5):
    """Number of linear iterations to reach error eps_target."""
    if lam >= 1 or lam <= 0:
        return float('inf')
    return math.ceil(math.log(eps_target / eps_start) / math.log(lam))


def n_iters_steffensen(eps_target=1e-15, eps_start=0.5):
    """Number of Steffensen iterations: O(log log) but practically ~5 for double prec."""
    # eps_n ~ K eps_{n-1}^2; n iter -> eps ~ eps_0^{2^n}
    n = 0; e = eps_start
    while e > eps_target and n < 50:
        e = e * e  # crude bound
        n += 1
    return n


def main():
    print("=" * 80)
    print("Pandrosion dD vs paper-0 scaling: do we still need scaling at large d?")
    print("=" * 80)
    print()
    print("Theory: for p=2, large d,  lambda_dD(x) -> (1 - 1/x) / d.")
    print("So increasing d shrinks lambda regardless of x -- dD absorbs scaling!")
    print()
    test_x = [10.0, 100.0, 10000.0, 1e8]
    print(f"  {'x':>10}  {'d=2 lambda':>14}  {'d=5 lambda':>14}  {'d=32 lambda':>14}  "
          f"{'asymp (1-1/x)/d, d=32':>22}")
    print("-" * 88)
    for x in test_x:
        l2 = lambda_dD(x, 2)
        l5 = lambda_dD(x, 5)
        l32 = lambda_dD(x, 32)
        asymp32 = (1 - 1 / x) / 32
        print(f"  {x:>10.1e}  {l2:>14.6e}  {l5:>14.6e}  {l32:>14.6e}  {asymp32:>22.6e}")

    print()
    print("Iterations needed to reach eps = 1e-15:")
    print(f"  {'x':>10}  {'2D (lin)':>10}  {'5D (lin)':>10}  {'32D (lin)':>10}  "
          f"{'2D Steff':>10}  {'32D Steff':>10}")
    print("-" * 80)
    for x in test_x:
        l2 = lambda_dD(x, 2)
        l5 = lambda_dD(x, 5)
        l32 = lambda_dD(x, 32)
        n_2d_lin = n_iters_to_eps(l2)
        n_5d_lin = n_iters_to_eps(l5)
        n_32d_lin = n_iters_to_eps(l32)
        n_steff_2d = n_iters_steffensen()
        n_steff_32d = n_iters_steffensen()
        # Convert inf to '>50'
        def fmt(n):
            return f"{n:>10}" if n < 1000 else f"{'>1000':>10}"
        print(f"  {x:>10.1e}  {fmt(n_2d_lin)}  {fmt(n_5d_lin)}  "
              f"{fmt(n_32d_lin)}  {n_steff_2d:>10}  {n_steff_32d:>10}")

    print()
    print("=" * 80)
    print("Paper-0 sec:scaling at large x (2D + scaling):")
    print(f"  For x = 10000, p = 2: A = 100, x' = 1.0  (exact integer scaling!)")
    print(f"  lambda_2D(x'=1) = 0  (degenerate, perfect)  -> 1 iteration!")
    print()
    print("With scaling, 2D matches dD with no scaling at all.  Both reach")
    print("eps = 1e-15 in O(1) iterations for any large x.")
    print()
    print("CONCLUSION:")
    print("  1. Pandrosion dD with d >= 32 has lambda ~ 1/d for ANY x,")
    print("     regardless of x's magnitude.  No paper-0 scaling needed.")
    print("  2. Paper-0 scaling on 2D achieves similar effect via integer A trick.")
    print("  3. Combining both (scaling + dD) is overkill: 2D + scaling is")
    print("     simpler and equally fast.")
    print("  4. dD's true advantage is for x of MODERATE size (x in [1, 100]),")
    print("     where 2D is slow (lambda ~ 0.5) and scaling barely helps.")
    print("=" * 80)


def verify_asymptotic():
    """Verify (1-s_*) ~ (x-1)/(dx) for large d."""
    print()
    print("Verification of asymptotic (1-s_*) ~ (x-1)/(dx):")
    print(f"  {'x':>8}  {'d':>4}  {'1-s_*':>16}  {'(x-1)/(dx)':>16}  {'ratio':>10}")
    print("-" * 64)
    for x in [10.0, 100.0, 10000.0]:
        for d in [4, 16, 64, 256]:
            s_star = fixed_point_p2(x, d)
            actual = 1 - s_star
            asympt = (x - 1) / (d * x)
            ratio = actual / asympt if asympt > 0 else float('nan')
            print(f"  {x:>8.0f}  {d:>4}  {actual:>16.10f}  {asympt:>16.10f}  {ratio:>10.6f}")


if __name__ == "__main__":
    main()
    verify_asymptotic()
