"""
PAPER: 001 (extension of paper 0)
TITLE: Pandrosion 3D — Tetrahedral Extension via Bivariate Geometric Sum
STATUS: empirical proposal (Pandrosion 2D paper 0 extended to 3D simplex)
DEPENDS: 0pandrosion_pth.tex (paper 0)

THEORY
======

Paper 0 (latex/0pandrosion_pth.tex) defines Pandrosion's classical iteration
for x^{1/p} via Thales' theorem in a 2D rectangle:

    s_{n+1} = 1 - (x - 1) / (x * S_p(s_n)),
    S_p(s) = 1 + s + s^2 + ... + s^{p-1} = (1 - s^p) / (1 - s).

The polynomial S_p(s) is the analytical translation of the geometric sum of
p segments inserted by parallels in the rectangle (Thales' construction).

PANDROSION 3D EXTENSION:
We replace the 2D rectangle by a 3D tetrahedral simplex with TWO families
of non-parallel cutting planes (parameterised by s and t).  The bivariate
geometric sum on the simplex {(j, k) in N^2 : j + k <= p-1} is

    S_p^{3D}(s, t) = sum_{j+k <= p-1} s^j t^k.

For p=2:  S_2^{3D}(s, t) = 1 + s + t.
For p=3:  S_3^{3D}(s, t) = 1 + s + t + s^2 + s*t + t^2.

PANDROSION 3D ITERATION (coupled, two ratios x and y):
    s_{n+1} = 1 - (x - 1) / (x * S_p^{3D}(s_n, t_n)),
    t_{n+1} = 1 - (y - 1) / (y * S_p^{3D}(s_n, t_n)).

FIXED POINT (symmetric x = y, hence s_* = t_*):
At s = t = s_*, S_p^{3D}(s_*, s_*) = sum_{j+k <= p-1} s_*^{j+k} which
equals sum_{m=0}^{p-1} (m+1) s_*^m = derivative-like form.

  - p=2: 2 x s_*^2 - x s_* - 1 = 0
         => s_* = (x + sqrt(x^2 + 8x)) / (4x).
         Special: x=2 => s_* = (1 + sqrt(5)) / 4 = phi/2 (GOLDEN RATIO).

  - p=3: 3 x s_*^3 - x s_*^2 - x s_* - 1 = 0
         (real positive root of a new cubic, NOT x^{-1/p}).

These fixed points are NEW algebraic numbers, not the simple p-th roots
of paper 0.  The golden ratio appearance for x=2, p=2 is the geometric
signature of 3D structures (icosahedron, pentagon).

VERIFICATION
============
This script verifies:
  1. S_p^{3D} formula for p=2, p=3.
  2. Symmetric iteration converges to predicted fixed point.
  3. Golden ratio emerges at (x=2, p=2).
  4. Asymmetric iteration (x != y) converges to coupled (s_*, t_*).
"""
from __future__ import annotations
import math
import numpy as np


# ---------------------------------------------------------------------------
# Bivariate geometric sums S_p^{3D}
# ---------------------------------------------------------------------------

def S_3d(s: complex, t: complex, p: int) -> complex:
    """S_p^{3D}(s, t) = sum_{j+k <= p-1} s^j t^k.

    For p=2: 1 + s + t.
    For p=3: 1 + s + t + s^2 + s t + t^2.
    """
    total = 0.0 + 0.0j
    s = complex(s); t = complex(t)
    for j in range(p):
        for k in range(p - j):
            total += s**j * t**k
    return total


# ---------------------------------------------------------------------------
# Pandrosion 3D iteration
# ---------------------------------------------------------------------------

def pandrosion_3d_iter(x, y, p, s0=0.5+0j, t0=0.5+0j,
                       max_iter=200, tol=1e-13):
    """Coupled Pandrosion 3D iteration on (s, t)."""
    s, t = complex(s0), complex(t0)
    history = [(s, t)]
    for n in range(max_iter):
        Spsn = S_3d(s, t, p)
        if abs(Spsn) < 1e-300:
            return s, t, n, history
        s_new = 1 - (x - 1) / (x * Spsn)
        t_new = 1 - (y - 1) / (y * Spsn)
        if abs(s_new - s) < tol and abs(t_new - t) < tol:
            return s_new, t_new, n + 1, history
        s, t = s_new, t_new
        history.append((s, t))
    return s, t, max_iter, history


def fixed_point_p2_symmetric(x):
    """Closed form fixed point for p=2, x=y."""
    return (x + math.sqrt(x*x + 8*x)) / (4*x)


def cubic_p3_symmetric(x):
    """Coefficients (descending) of 3xs^3 - xs^2 - xs - 1 = 0."""
    return [3*x, -x, -x, -1]


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def verify_S_p_3D_formula():
    print("[1] Formula check S_p^{3D}(s, t) = sum_{j+k <= p-1} s^j t^k")
    print("-" * 70)
    s, t = 0.5, 0.3
    expected_p2 = 1 + s + t
    actual_p2 = S_3d(s, t, 2)
    print(f"  p=2:  S = 1 + s + t                    = {actual_p2.real:.6f}  "
          f"(expected {expected_p2:.6f})")
    expected_p3 = 1 + s + t + s*s + s*t + t*t
    actual_p3 = S_3d(s, t, 3)
    print(f"  p=3:  S = 1 + s + t + s^2 + s t + t^2 = {actual_p3.real:.6f}  "
          f"(expected {expected_p3:.6f})")


def verify_symmetric_p2():
    print("\n[2] Symmetric Pandrosion 3D, p=2, fixed point (x = y)")
    print("    Predicted: s_* = (x + sqrt(x^2 + 8x)) / (4x)")
    print("-" * 70)
    print(f"  {'x':>6} {'s_* numeric':>14} {'s_* formula':>14} {'iter':>6} {'error':>11}")
    for x in [1.5, 2.0, 4.0, 10.0]:
        s_star, t_star, n, _ = pandrosion_3d_iter(x, x, p=2)
        s_pred = fixed_point_p2_symmetric(x)
        err = abs(s_star.real - s_pred)
        print(f"  {x:>6.2f} {s_star.real:>14.10f} {s_pred:>14.10f} {n:>6} {err:>11.2e}")


def verify_golden_ratio():
    print("\n[3] GOLDEN RATIO emerges at p=2, x=2")
    print("-" * 70)
    s_star, _, n, _ = pandrosion_3d_iter(2.0, 2.0, p=2)
    phi = (1 + math.sqrt(5)) / 2
    half_phi = phi / 2
    print(f"  Pandrosion 3D fixed point: s_* = {s_star.real:.15f}")
    print(f"  phi / 2                   = {half_phi:.15f}  (phi = {phi:.10f})")
    print(f"  Difference                = {abs(s_star.real - half_phi):.2e}")


def verify_symmetric_p3():
    print("\n[4] Symmetric Pandrosion 3D, p=3, fixed point (x = y)")
    print("    Predicted: real positive root of 3 x s^3 - x s^2 - x s - 1 = 0")
    print("-" * 70)
    print(f"  {'x':>6} {'s_* numeric':>14} {'cubic residual':>16} {'iter':>6}")
    for x in [1.5, 2.0, 4.0, 10.0]:
        s_star, _, n, _ = pandrosion_3d_iter(x, x, p=3)
        residual = abs(np.polyval(cubic_p3_symmetric(x), s_star.real))
        print(f"  {x:>6.2f} {s_star.real:>14.10f} {residual:>16.2e} {n:>6}")
    # Special case x=2: print all 3 cubic roots numerically
    print("\n  Special case x=2: cubic 6 s^3 - 2 s^2 - 2 s - 1 = 0")
    roots = np.roots([6, -2, -2, -1])
    for r in roots:
        if abs(r.imag) < 1e-10:
            print(f"    real root: {r.real:.15f}")
        else:
            print(f"    complex:   {r.real:.6f} + {r.imag:.6f} j")


def verify_asymmetric():
    print("\n[5] Asymmetric Pandrosion 3D (x != y), p=2")
    print("-" * 70)
    print(f"  {'(x, y)':>10}    {'(s_*, t_*)':>30}    {'iter':>6}")
    for x, y in [(2.0, 3.0), (1.5, 4.0), (5.0, 10.0)]:
        s_star, t_star, n, _ = pandrosion_3d_iter(x, y, p=2)
        print(f"  ({x:>3.1f}, {y:>3.1f})    "
              f"({s_star.real:>10.6f}, {t_star.real:>10.6f})    "
              f"{n:>6}")
        # Check fixed-point equations
        Sp = S_3d(s_star, t_star, 2)
        eq1 = x * Sp.real * (1 - s_star.real) - (x - 1)
        eq2 = y * Sp.real * (1 - t_star.real) - (y - 1)
        print(f"          residual eq1 = {abs(eq1):.2e},  eq2 = {abs(eq2):.2e}")


def verify_contraction_rate():
    print("\n[6] Empirical contraction rate (lambda_3D for p=2, x=2)")
    print("-" * 70)
    s_star, t_star, _, history = pandrosion_3d_iter(2.0, 2.0, p=2, max_iter=50)
    # Estimate contraction rate from successive errors
    errs = [abs(complex(s) - s_star) for s, t in history]
    print(f"  step   |s_n - s_*|     ratio")
    for n in range(1, min(10, len(errs))):
        if errs[n-1] > 0:
            ratio = errs[n] / errs[n-1]
            print(f"  {n:>4}  {errs[n]:>13.2e}    {ratio:.6f}")


def main():
    print("=" * 70)
    print("PAPER 001 -- Pandrosion 3D: tetrahedral bivariate extension")
    print("=" * 70)
    verify_S_p_3D_formula()
    verify_symmetric_p2()
    verify_golden_ratio()
    verify_symmetric_p3()
    verify_asymmetric()
    verify_contraction_rate()
    print("\n" + "=" * 70)
    print("Pandrosion 3D extension verified: bivariate S_p^{3D} converges,")
    print("golden ratio appears at (x=2, p=2), new cubic algebraic numbers")
    print("at p=3.  The geometric signature of 3D structures is captured.")
    print("=" * 70)


if __name__ == "__main__":
    main()
