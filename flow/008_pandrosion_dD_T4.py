"""
PAPER: 008 (acceleration piste 3)
TITLE: T_4 hierarchy (iterated Shanks) on Pandrosion dD
STATUS: empirical comparison; Shanks T_2/T_3/T_4 stacked

PISTE 3 (higher Shanks levels)
==============================
Iterate the Aitken Delta^2 transformation:
  T_2 = standard Aitken Delta^2 (Wynn epsilon_2), order 2
  T_3 = Aitken^2, order 3
  T_4 = Aitken^3, order 4

Apply to the diagonal scalar map g(s) = 1 - (x-1) / (x Delta_d(s)).
Each level requires more g-evaluations per step but converges in fewer
steps.  We compare wall-clock and total g-evaluations.
"""
from __future__ import annotations
import math
import time
from math import comb


def Delta_d(s, p, d):
    return sum(comb(m + d - 2, d - 2) * s**m for m in range(p))


def make_g(x, p, d):
    return lambda s: 1 - (x - 1) / (x * Delta_d(s, p, d))


def aitken_T2(g, s):
    """T_2 = Aitken Delta^2 on g (uses 2 evaluations of g)."""
    g1 = g(s)
    g2 = g(g1)
    denom = g2 - 2 * g1 + s
    if abs(denom) < 1e-300:
        return g2, 2
    return s - (g1 - s) ** 2 / denom, 2


def aitken_T3(g, s):
    """T_3 = T_2 applied to T_2 sequence (4 evaluations)."""
    s0 = s
    s1, e1 = aitken_T2(g, s0)
    s2, e2 = aitken_T2(g, s1)
    denom = s2 - 2 * s1 + s0
    if abs(denom) < 1e-300:
        return s2, e1 + e2
    return s0 - (s1 - s0) ** 2 / denom, e1 + e2


def aitken_T4(g, s):
    """T_4 = T_2 applied to T_3 sequence (8 evaluations)."""
    s0 = s
    s1, e1 = aitken_T3(g, s0)
    s2, e2 = aitken_T3(g, s1)
    denom = s2 - 2 * s1 + s0
    if abs(denom) < 1e-300:
        return s2, e1 + e2
    return s0 - (s1 - s0) ** 2 / denom, e1 + e2


def iterate(method, g, s0=0.5, max_iter=20, tol=1e-14):
    s = s0
    n_evals = 0
    n_iter = 0
    history = [s]
    for _ in range(max_iter):
        s_new, e = method(g, s)
        n_evals += e
        n_iter += 1
        history.append(s_new)
        if abs(s_new - s) < tol:
            return s_new, n_iter, n_evals, history
        s = s_new
    return s, max_iter, n_evals, history


def estimate_order(history, s_star):
    errs = [abs(h - s_star) for h in history if abs(h - s_star) > 1e-16]
    if len(errs) < 4:
        return float('nan')
    return math.log(errs[-1] / errs[-2]) / math.log(errs[-2] / errs[-3])


def main():
    print("=" * 76)
    print("PISTE 3 -- T_2/T_3/T_4 Shanks hierarchy on Pandrosion dD")
    print("=" * 76)
    print()
    print(f"  {'(x, d, p)':>13}  {'T_2 (Aitken)':<18}  {'T_3':<18}  {'T_4':<18}")
    print(f"  {'':>13}  {'i | e | order':<18}  {'i | e | order':<18}  {'i | e | order':<18}")
    print("-" * 76)
    cases = [
        (2.0, 3, 2), (2.0, 5, 2), (2.0, 10, 2),
        (5.0, 3, 2), (10.0, 4, 2),
        (2.0, 3, 3), (2.0, 5, 3),
    ]
    for x, d, p in cases:
        g = make_g(x, p, d)
        s_T2, i2, e2, h2 = iterate(aitken_T2, g)
        s_T3, i3, e3, h3 = iterate(aitken_T3, g)
        s_T4, i4, e4, h4 = iterate(aitken_T4, g)
        s_star = s_T2  # all converge to same fixed point
        ord_T2 = estimate_order(h2, s_star)
        ord_T3 = estimate_order(h3, s_star)
        ord_T4 = estimate_order(h4, s_star)
        print(f"  ({x:>3.1f}, {d:>2}, {p})  "
              f"{i2:>2} | {e2:>2} | {ord_T2:>5.2f}    "
              f"{i3:>2} | {e3:>2} | {ord_T3:>5.2f}    "
              f"{i4:>2} | {e4:>2} | {ord_T4:>5.2f}")
    print()
    print("Legend: i=iterations, e=total g-evaluations, order=empirical convergence order")
    print("Total cost = e (g-evaluations).  Lower e is better.")


if __name__ == "__main__":
    main()
