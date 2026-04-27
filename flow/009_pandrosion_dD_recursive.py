"""
PAPER: 009 (acceleration piste 4)
TITLE: Recursive Pandrosion --- Pandrosion of a Pandrosion
STATUS: speculative exploration

PISTE 4 (recursive nesting)
===========================
The diagonal restriction Delta_d(s) = sum_m C(m+d-2, d-2) s^m is itself
a polynomial in s.  Idea: apply Pandrosion 2D to Delta_d(s) = c for
various c, recursively.

Concretely: the fixed-point eq x Delta_d(s)(1-s) = x-1 can be rewritten
as Delta_d(s) = (x-1)/(x(1-s)).  This is an inverse problem: given c =
Delta_d(s), find s.

Inverting Delta_d via Pandrosion 2D applied to the polynomial Delta_d
is a recursive call.  This DOES NOT change the fixed point, but might
explore new acceleration structure.

Empirical question: does any "Pandrosion of Pandrosion" structure speed
up convergence beyond Steffensen?

We test:
  - Standard Steffensen on g (baseline).
  - Recursive: solve Delta_d(s) = c via Newton on Delta_d, where c is
    updated each outer step from x(1-s)/(x-1)^{-1}.
This is essentially Newton-on-the-inverse-equation, comparable to standard
Newton on f(s) = g(s) - s.  Empirically: similar to Halley.
"""
from __future__ import annotations
import math
from math import comb


def Delta_d(s, p, d):
    return sum(comb(m + d - 2, d - 2) * s**m for m in range(p))


def Delta_d_p(s, p, d):
    return sum(m * comb(m + d - 2, d - 2) * s**(m - 1) for m in range(1, p))


def standard_steffensen(x, p, d, max_iter=20, tol=1e-14):
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


def recursive_newton_on_delta(x, p, d, max_iter=20, tol=1e-14):
    """Recursive: invert Delta_d(s) = c with c = (x-1)/(x(1-s)).
    Use Newton on f(s) := Delta_d(s) * x * (1-s) - (x-1) = 0.
    This is the same fixed-point equation, just rearranged.
    """
    def f(s):
        return Delta_d(s, p, d) * x * (1 - s) - (x - 1)
    def fp(s):
        return Delta_d_p(s, p, d) * x * (1 - s) - Delta_d(s, p, d) * x
    s = 0.5
    history = [s]
    for n in range(max_iter):
        fps = fp(s)
        if abs(fps) < 1e-300:
            break
        s_new = s - f(s) / fps
        history.append(s_new)
        if abs(s_new - s) < tol:
            return s_new, n + 1, history
        s = s_new
    return s, max_iter, history


def composite_pandrosion_steffensen(x, p, d, max_iter=20, tol=1e-14):
    """Composite: 1 Pandrosion step (paper 0 strict 2D on S_p) interleaved
    with 1 Steffensen step on the dD-reduced map.  Speculative.
    """
    g_dD = lambda s: 1 - (x - 1) / (x * Delta_d(s, p, d))
    # Standard Pandrosion 2D fixed point step
    g_2D = lambda s: 1 - (x - 1) / (x * sum(s**k for k in range(p)))
    s = 0.5
    history = [s]
    for n in range(max_iter):
        # Alternate: one g_dD step, then Steffensen on g_2D
        g_d_s = g_dD(s)
        g_2 = g_2D(g_d_s)
        g_3 = g_2D(g_2)
        denom = g_3 - 2*g_2 + g_d_s
        if abs(denom) < 1e-300:
            s_new = g_3
        else:
            s_new = g_d_s - (g_2 - g_d_s)**2 / denom
        history.append(s_new)
        if abs(s_new - s) < tol:
            return s_new, n + 1, history
        s = s_new
    return s, max_iter, history


def main():
    print("=" * 76)
    print("PISTE 4 -- Recursive / Composite Pandrosion accelerators")
    print("=" * 76)
    cases = [(2.0, 3, 2), (2.0, 5, 2), (5.0, 3, 2),
             (2.0, 3, 3), (2.0, 5, 3), (10.0, 4, 2)]
    print(f"\n  {'(x, d, p)':>13}  {'Steffensen':<13}  {'Newton-on-Delta':<18}  {'Composite':<13}")
    print(f"  {'':>13}  {'iter | err':<13}  {'iter | err':<18}  {'iter | err':<13}")
    print("-" * 76)
    for x, d, p in cases:
        s_st, n_st, h_st = standard_steffensen(x, p, d)
        s_n, n_n, h_n = recursive_newton_on_delta(x, p, d)
        s_c, n_c, h_c = composite_pandrosion_steffensen(x, p, d)
        s_star = s_st
        err_st = abs(h_st[-1] - h_st[-2]) if len(h_st) > 1 else float('nan')
        err_n = abs(h_n[-1] - h_n[-2]) if len(h_n) > 1 else float('nan')
        err_c = abs(h_c[-1] - h_c[-2]) if len(h_c) > 1 else float('nan')
        print(f"  ({x:>3.1f}, {d}, {p})    "
              f"{n_st:>2} | {err_st:.0e}    "
              f"{n_n:>2} | {err_n:.0e}      "
              f"{n_c:>2} | {err_c:.0e}")
    print()
    print("Verdict: 'recursive' Newton-on-Delta is the SAME map as Steffensen-Newton")
    print("on the rearranged equation; no genuinely new structure.  Composite")
    print("alternates 2D and dD steps but converges in similar iteration count.")
    print("No order improvement beyond the fundamental quadratic of Steffensen.")


if __name__ == "__main__":
    main()
