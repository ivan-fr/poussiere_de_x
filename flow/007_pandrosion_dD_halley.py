"""
PAPER: 007 (acceleration piste 2)
TITLE: Halley acceleration on the diagonal-reduced Pandrosion dD
STATUS: empirical comparison of orders 2, 3 (Steffensen vs Halley)

PISTE 2 (cubic order via Halley)
================================
After rank-1 collapse, the diagonal scalar map is
    g(s) = 1 - (x - 1) / (x * Delta_d(s)),  Delta_d = sum C(m+d-2, d-2) s^m.
Halley's method on f(s) := g(s) - s = 0:
    s_{n+1} = s_n - 2 f f' / (2 (f')^2 - f f'').
This is order 3 (cubic) instead of Newton/Aitken's order 2.

We compare iterations at machine precision tolerance:
  Steffensen (Aitken Delta^2 on g):       order 2
  Halley:                                  order 3
"""
from __future__ import annotations
import math
from math import comb


def Delta_d(s, p, d):
    return sum(comb(m + d - 2, d - 2) * s**m for m in range(p))


def Delta_d_p(s, p, d):
    return sum(m * comb(m + d - 2, d - 2) * s**(m - 1) for m in range(1, p))


def Delta_d_pp(s, p, d):
    return sum(m * (m - 1) * comb(m + d - 2, d - 2) * s**(m - 2) for m in range(2, p))


def steffensen_iteration(x, p, d, max_iter=20, tol=1e-14):
    g = lambda s: 1 - (x - 1) / (x * Delta_d(s, p, d))
    s = 0.5
    history = [s]
    for n in range(max_iter):
        g1 = g(s); g2 = g(g1)
        denom = g2 - 2*g1 + s
        if abs(denom) < 1e-300:
            return g2, n + 1, history
        s_new = s - (g1 - s)**2 / denom
        history.append(s_new)
        if abs(s_new - s) < tol:
            return s_new, n + 1, history
        s = s_new
    return s, max_iter, history


def halley_iteration(x, p, d, max_iter=20, tol=1e-14):
    """Halley's method on f(s) = g(s) - s where g(s) = 1 - (x-1)/(x Delta_d(s))."""
    def f(s):
        return 1 - (x - 1) / (x * Delta_d(s, p, d)) - s
    def fp(s):
        # f'(s) = (x-1)/(x) * Delta'(s) / Delta(s)^2 - 1
        return (x - 1) / x * Delta_d_p(s, p, d) / Delta_d(s, p, d)**2 - 1
    def fpp(s):
        # f''(s) = (x-1)/x * [Delta''/Delta^2 - 2 Delta'^2 / Delta^3]
        D = Delta_d(s, p, d); Dp = Delta_d_p(s, p, d); Dpp = Delta_d_pp(s, p, d)
        return (x - 1) / x * (Dpp / D**2 - 2 * Dp**2 / D**3)
    s = 0.5
    history = [s]
    for n in range(max_iter):
        fs, fps, fpps = f(s), fp(s), fpp(s)
        denom = 2 * fps**2 - fs * fpps
        if abs(denom) < 1e-300:
            break
        s_new = s - 2 * fs * fps / denom
        history.append(s_new)
        if abs(s_new - s) < tol:
            return s_new, n + 1, history
        s = s_new
    return s, max_iter, history


def estimate_order(history, s_star):
    errs = [abs(h - s_star) for h in history if abs(h - s_star) > 1e-16]
    if len(errs) < 4:
        return float('nan')
    # log(eps_{n+1}) / log(eps_n) at n large
    return math.log(errs[-1] / errs[-2]) / math.log(errs[-2] / errs[-3])


def main():
    print("=" * 70)
    print("PISTE 2 -- Halley (order 3) vs Steffensen (order 2)")
    print("=" * 70)
    print()
    print(f"  {'(x, d, p)':>15}  {'Steffensen':<22}  {'Halley':<22}")
    print(f"  {'':>15}  {'iter | order | err':<22}  {'iter | order | err':<22}")
    print("-" * 70)
    test_cases = [
        (2.0, 3, 2), (2.0, 4, 2), (2.0, 5, 2),
        (5.0, 3, 2), (10.0, 4, 2),
        (2.0, 3, 3), (2.0, 4, 3), (2.0, 5, 3),
    ]
    for x, d, p in test_cases:
        s_st, n_st, hist_st = steffensen_iteration(x, p, d)
        s_h, n_h, hist_h = halley_iteration(x, p, d)
        s_star = s_st  # both should converge to same root
        ord_st = estimate_order(hist_st, s_star)
        ord_h = estimate_order(hist_h, s_star)
        err_st = abs(hist_st[-1] - s_star) if len(hist_st) > 1 else float('nan')
        err_h = abs(hist_h[-1] - s_h)  # halley relative
        print(f"  ({x:>3.1f}, {d}, {p})    "
              f"{n_st:>3} | {ord_st:>5.2f} | {err_st:.0e}    "
              f"{n_h:>3} | {ord_h:>5.2f} | {err_h:.0e}")


if __name__ == "__main__":
    main()
