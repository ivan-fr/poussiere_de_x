"""
Numerical verification of the reviewer's technical complaints.
Checks:
  (1) limit of lambda_{p,x} as p -> infinity with x fixed
  (2) monotonicity claim "lambda -> 1 as p -> infinity" (Lean paper)
  (3) |lambda_{p,x}| < 1 for x < 1 (Prop 11.1 of expository paper)
"""

from __future__ import annotations
import math


def lam(p: int, x: float) -> float:
    alpha = x ** (1.0 / p)
    num = (alpha - 1.0) * sum(k * alpha ** (p - 1 - k) for k in range(1, p))
    den = sum(alpha ** k for k in range(p))
    return num / den


# ---------------------------------------------------------------------------
# (1) Limit as p -> infinity
# ---------------------------------------------------------------------------
print("=" * 72)
print("(1) lim_{p -> oo} lambda_{p,x} for fixed x > 1")
print()
print("Two candidate limits:")
print("  (a) ln(x) / (1 + ln(x))   [article 1, Prop 9.1(3)]")
print("  (b) 1 - ln(x) / (x - 1)   [reviewer-implied analytic limit]")
print()
print("  x     |  lam(p=50) lam(p=200) lam(p=1000)  | (a) cand    (b) cand")
print("-" * 72)
for x in [1.5, 2.0, 3.0, 10.0, 100.0]:
    cand_a = math.log(x) / (1.0 + math.log(x))
    cand_b = 1.0 - math.log(x) / (x - 1.0)
    l50 = lam(50, x)
    l200 = lam(200, x)
    l1000 = lam(1000, x)
    print(f"  {x:<5} |  {l50:.5f}    {l200:.5f}    {l1000:.5f}     | "
          f"{cand_a:.5f}    {cand_b:.5f}")

# ---------------------------------------------------------------------------
# (2) Does lambda -> 1 as p -> infinity with x fixed?
# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(2) Does lambda -> 1 as p -> oo with x fixed? (Lean paper claim)")
print()
for x in [2.0, 10.0]:
    vals = [lam(p, x) for p in [10, 100, 1000, 10000]]
    print(f"  x = {x:<4}: lam(p=10) = {vals[0]:.5f}, "
          f"lam(p=100) = {vals[1]:.5f}, "
          f"lam(p=1000) = {vals[2]:.5f}, "
          f"lam(p=10000) = {vals[3]:.5f}")
    print(f"             1 - ln(x)/(x-1) = {1 - math.log(x)/(x-1):.5f}")

# ---------------------------------------------------------------------------
# (3) |lambda_{p,x}| < 1 for x < 1?  (Prop 11.1 of article 1)
# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(3) Is |lambda_{p,x}| < 1 for all 0 < x < 1?  (Prop 11.1 claim)")
print()
print("  (p, x)             | alpha=x^(1/p)  | lambda          | |lambda| < 1 ?")
print("-" * 72)
# Cases the paper lists (Table 6) are safe; try cases with x very small.
cases = [
    (2, 0.5), (3, 0.5), (3, 0.125),                # paper's table
    (3, 0.05), (3, 0.01), (3, 0.001),              # small x
    (5, 0.01), (10, 0.001),
    (3, 1e-6),                                     # extreme
]
for p, x in cases:
    alpha = x ** (1.0 / p)
    l = lam(p, x)
    ok = "YES" if abs(l) < 1.0 else "***  NO  ***"
    print(f"  (p={p}, x={x:<8g}) | alpha={alpha:.5f}  | "
          f"lambda={l:+.5f}  | {ok}")

# ---------------------------------------------------------------------------
# (4) Threshold alpha* for p = 3 such that |lambda| < 1 iff alpha > alpha*
# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(4) For p=3, threshold alpha* where |lambda_{3,alpha^3}| = 1")
print()
# |lambda| = (1-a)(2+a)/(1+a+a^2) = 1  =>  2a^2 + 2a - 1 = 0
# a = (-2 + sqrt(12))/4 = (-1 + sqrt(3))/2.
a_star = (-1 + math.sqrt(3.0)) / 2.0
print(f"  Analytic threshold: alpha* = (-1 + sqrt(3))/2 = {a_star:.6f}")
print(f"  Corresponding x*   = alpha*^3              = {a_star**3:.6f}")
print(f"  lambda_{{3,x}} at alpha just above/below:")
for a in [a_star - 0.001, a_star, a_star + 0.001]:
    x_ = a ** 3
    l = lam(3, x_)
    print(f"    alpha = {a:.6f}, x = {x_:.6f}, lambda = {l:+.6f}")
