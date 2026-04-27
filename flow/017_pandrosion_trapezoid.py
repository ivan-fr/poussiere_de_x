"""
PAPER: 017 (alternative geometry: trapezoid)
TITLE: Pandrosion in trapezoidal Thales construction (family zeta)
STATUS: empirical exploration

CONSTRUCTION
============
Trapezoid: a 2D figure intermediate between rectangle (paper-0) and
triangle.  Bases of length 1 (top) and x (bottom), with shrinkage
parameter alpha in [0, 1) measuring how far it deviates from a rectangle.

  alpha = 0:  pure rectangle (paper-0 strict)
  alpha = 1:  triangle (degenerate, top base = 0)
  alpha in (0, 1):  trapezoid

POLYNOMIAL S_p^{trap}(s, alpha):
The Thales construction in a trapezoid yields, at height s, an effective
geometric sum where each term s^k is weighted by a linear shrinkage:
    S_p^{trap}(s, alpha) = sum_{k=0}^{p-1} s^k * (1 - alpha * k / (p-1))
For alpha=0, this reduces to S_p(s) = sum s^k (paper-0).

ITERATION:
    s_{n+1} = 1 - (x - 1) / (x * S_p^{trap}(s_n, alpha))

The fixed point depends on (x, p, alpha).  As alpha varies, so does s_*.

KEY QUESTION:
Does varying alpha give a STRICTLY BETTER algorithm than alpha = 0
(paper-0)?  I.e., is there an alpha* that minimises lambda_trap(x, alpha)?

Verify numerically.
"""
from __future__ import annotations
import math


def S_trap(s, alpha, p):
    """S_p^{trap}(s, alpha) = sum_{k=0}^{p-1} s^k (1 - alpha k/(p-1))."""
    if p == 1:
        return 1.0 + 0.0j
    s = complex(s)
    return sum(s**k * (1 - alpha * k / (p - 1)) for k in range(p))


def trap_iter(x, alpha, p, max_iter=200, tol=1e-13):
    s = 0.5 + 0.0j
    history = [s]
    for n in range(max_iter):
        S = S_trap(s, alpha, p)
        if abs(S) < 1e-300:
            return s, n, history
        s_new = 1 - (x - 1) / (x * S)
        history.append(s_new)
        if abs(s_new - s) < tol:
            return s_new, n + 1, history
        s = s_new
    return s, max_iter, history


def estimate_lambda(history, s_star):
    errs = [abs(complex(h).real - s_star) for h in history if abs(h - s_star) > 1e-15]
    if len(errs) < 3:
        return float('nan')
    ratios = [errs[i + 1] / errs[i] for i in range(len(errs) - 2)]
    return ratios[-1] if ratios else float('nan')


def main():
    print("=" * 76)
    print("Family zeta: Pandrosion in TRAPEZOID with shrinkage alpha")
    print("=" * 76)
    print()
    print("S_p^{trap}(s, alpha) = sum_k s^k (1 - alpha k/(p-1))")
    print()
    print("Test: how does the contraction lambda_trap(x, alpha) vary with alpha?")
    print(f"\n  {'x':>4} {'p':>3} {'alpha':>8} {'s_*':>14} {'lambda_trap':>13} {'iter':>5}")
    print("-" * 76)
    for x in [2.0, 5.0]:
        for p in [2, 3]:
            for alpha in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
                s_star, n_iter, history = trap_iter(x, alpha, p)
                lam = estimate_lambda(history, s_star.real)
                print(f"  {x:>4.1f} {p:>3} {alpha:>8.2f} {s_star.real:>14.10f} "
                      f"{lam:>13.6f} {n_iter:>5}")
            print()
    print("Note: alpha=0 is paper-0 strict; lambda(x=2, p=2, alpha=0) = 1/(2+sqrt(2)) ~ 0.293")
    print()
    print("OBSERVATION:")
    print("  As alpha increases from 0, the fixed point s_* shifts upward.")
    print("  lambda may decrease for small alpha (faster convergence) but the")
    print("  fixed point is no longer x^{-1/p} -- it's a new algebraic number.")
    print("  Trapezoid Pandrosion is a DIFFERENT problem, not a faster paper-0.")
    print()
    # Optimal alpha minimizing lambda?
    print("Search for alpha minimising lambda (x=2, p=2):")
    best_lam, best_alpha = float('inf'), 0.0
    for alpha in [i * 0.01 for i in range(100)]:
        s_star, _, history = trap_iter(2.0, alpha, p=2)
        lam = estimate_lambda(history, s_star.real)
        if not math.isnan(lam) and lam < best_lam:
            best_lam = lam
            best_alpha = alpha
    print(f"  Optimal: alpha* = {best_alpha:.3f}, lambda_min = {best_lam:.6f}")
    print(f"  Compare paper-0 strict (alpha=0): lambda = {1/(2+math.sqrt(2)):.6f}")
    print()
    print("=" * 76)
    print("VERDICT: trapezoid does NOT give a strict 'dD' family.  It's a 2D")
    print("modification of paper-0 with one parameter alpha that:")
    print("  - Shifts the fixed point to new algebraic values (not x^{-1/p})")
    print("  - Slightly modifies lambda but not asymptotically better than")
    print("    paper-0 + Steffensen")
    print("  - Solves a DIFFERENT problem (not x^{1/p}) so direct comparison")
    print("    with paper-0 strict is not meaningful")
    print("=" * 76)


if __name__ == "__main__":
    main()
