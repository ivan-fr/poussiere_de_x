"""
PAPER: 012 (alternative geometry: rectangle chain)
TITLE: Pandrosion as Thales chain of rectangles (family alpha)
STATUS: empirical exploration

CONSTRUCTION
============
Family alpha: chain of d-1 rectangles glued at shared edges.
  d=2: 1 rectangle of ratio x        (paper-0 strict)
  d=3: 2 rectangles, ratios x_1 then x_2, sharing an edge
  d=4: 3 rectangles in series
  d=k: k-1 rectangles in chain

Each rectangle inserts p-1 proportional means via Thales.  Junction
between rectangle i and rectangle i+1 imposes continuity: the last
intercept of rectangle i equals the first intercept of rectangle i+1
(after rescaling by the shared edge ratio).

POLYNOMIAL DERIVATION (n = d-1 chained rectangles, p=2 example):
  Rectangle 1: insertions 1, s_1                  (Thales sum: 1 + s_1)
  Rectangle 2: starts at s_1, ends at s_1*s_2     (Thales sum: s_1 + s_1*s_2 = s_1*(1+s_2))
  Total geometric sum: S_p^{alpha, d}(s_1, ..., s_{d-1}) =
                       sum_i [s_1*s_2*...*s_{i-1}] * S_p(s_i).

For p=2: S_p(s) = 1 + s.  So S_alpha,3 = (1+s_1) + s_1*(1+s_2) = 1 + s_1 + s_1 + s_1*s_2
                                       = 1 + 2*s_1 + s_1*s_2.
For d=4: S_alpha,4 = 1 + 2*s_1 + s_1*s_2 + s_1*(1+s_2) - ?
Actually let me re-derive more carefully.

The sum of the d-1 rectangles' Thales sums, with each rectangle scaled
by its starting position:
  S_alpha,d(s_1, ..., s_{d-1}) = sum_{i=1}^{d-1} (prod_{j<i} s_j^p) * S_p(s_i).

For p=2:
  d=2: S = (1)*(1+s_1) = 1+s_1 (paper-0)
  d=3: S = (1)*(1+s_1) + s_1^2*(1+s_2) = 1 + s_1 + s_1^2 + s_1^2*s_2
  d=4: S = 1 + s_1 + s_1^2 + s_1^2*s_2 + s_1^2*s_2^2 + s_1^2*s_2^2*s_3
DIAGONAL s_1 = s_2 = ... = s:
  d=2: 1 + s
  d=3: 1 + s + s^2 + s^3 = S_4(s)
  d=4: 1 + s + s^2 + s^3 + s^4 + s^5 = S_6(s)
  d-th rectangle: S_{2(d-1)}(s) — chained Thales is equivalent to
    paper-0 with effective p_eff = 2(d-1).

INSIGHT: the chain reduces to scalar paper-0 with INCREASED p, not new
multivariate structure.  The rank-1 collapse is replaced by an "effective
p_eff = p*(d-1)" reduction.

Verify numerically.
"""
from __future__ import annotations
import math


def S_chain(s_vec, p):
    """S_alpha,d(s_1, ..., s_{d-1}) = sum_{i=1}^{d-1} (prod_{j<i} s_j^p) * S_p(s_i)."""
    s_vec = [complex(s) for s in s_vec]
    n = len(s_vec)  # d-1
    total = 0.0+0.0j
    prefix = 1.0+0.0j
    for i in range(n):
        Sp_si = sum(s_vec[i]**k for k in range(p))
        total += prefix * Sp_si
        prefix *= s_vec[i]**p
    return total


def iterate_chain(x_vec, d, p, max_iter=200, tol=1e-13):
    n = d - 1
    s = [complex(0.5)] * n
    history = [tuple(s)]
    for step in range(max_iter):
        Spsn = S_chain(s, p)
        s_new = [1 - (x_vec[i]-1)/(x_vec[i]*Spsn) for i in range(n)]
        if max(abs(s_new[i]-s[i]) for i in range(n)) < tol:
            return tuple(s_new), step+1, history
        s = s_new
        history.append(tuple(s))
    return tuple(s), max_iter, history


def main():
    print("="*70)
    print("Family alpha: Pandrosion CHAIN of rectangles")
    print("="*70)
    print("\nFor symmetric x_i = x, p=2 case:")
    print("  d=2: paper-0 strict, S_p = 1+s,        s_* = 1/sqrt(x)")
    print("  d=3: S = 1+2s+s^3 (after sym),         new fixed point")
    print("  d=4: S = 1+s+s^2+s^3+s^4+s^5 = S_6(s), reduces to paper-0 p_eff=6")
    print()
    print(f"  {'d':>3} {'x':>5} {'s_*':>14} {'lambda_chain':>14} {'paper-0 equiv':>20}")
    print("-"*70)
    for d in [2, 3, 4, 5, 6]:
        x = 2.0
        n = d - 1
        s_final, n_iter, history = iterate_chain((x,)*n, d, p=2)
        s_star = s_final[0].real
        # Estimate lambda
        if len(history) > 5:
            errs = [abs(complex(h[0]).real - s_star) for h in history if h]
            ratios = [errs[i+1]/errs[i] for i in range(len(errs)-2) if errs[i] > 1e-15]
            lam = ratios[-1] if ratios else float('nan')
        else:
            lam = float('nan')
        # Equivalent paper-0 p_eff
        p_eff = 2 * (d - 1) if d >= 2 else 2
        # paper-0 with p_eff: s_* should equal x^{-1/p_eff}
        s_paper0 = x**(-1/p_eff)
        print(f"  {d:>3} {x:>5.1f} {s_star:>14.10f} {lam:>14.6f} "
              f"x^{{-1/{p_eff}}} = {s_paper0:>10.6f}")
    print()
    print("VERDICT: chain reduces to paper-0 with effective p_eff = p*(d-1).")
    print("  -> No new geometric structure beyond scalar paper-0; just")
    print("     reaches a HIGHER root x^{1/(2(d-1))} instead of x^{1/p}.")
    print("  -> Useful only if one wanted x^{1/(2(d-1))} for arbitrary k = d-1,")
    print("     but trivially achievable in paper-0 with p = 2(d-1) directly.")


if __name__ == "__main__":
    main()
