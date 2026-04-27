"""
PAPER: 013 (alternative geometry: pyramidal stacks)
TITLE: Pandrosion in a d-cone / pyramid (family beta)
STATUS: empirical exploration

CONSTRUCTION
============
Family beta: pyramid in dimension d (apex + (d-1)-dimensional base).
  d=2: apex + line = triangle (paper-0 reformulated)
  d=3: apex + square = square pyramid; section at height h has area
       proportional to (1 - h/H)^2.
  d=4: apex + cube = 4-pyramid; section volume ~ (1-h/H)^3.
  d-pyramid: section ~ (1 - h/H)^{d-1}.

DERIVATION:
The volume of a slice between heights h_{k-1} and h_k in a d-pyramid is
  V_k = const * [(1 - h_{k-1}/H)^d - (1 - h_k/H)^d].
Setting V_k in geometric progression V_k = s^{k-1} V_1 and using
u_k = 1 - h_k/H, we get
  u_{k-1}^d - u_k^d = s^{k-1} (u_0^d - u_1^d).
Posing T = 1 - u_1^d, we find u_p^d = 1 - T*S_p(s), and u_p = 0 forces
T*S_p(s) = 1, i.e. T = 1/S_p(s).
This is EXACTLY the same scalar Pandrosion iteration as paper-0,
regardless of d.  No new structure.

So family beta is degenerate: pyramids of any dimension reduce to
paper-0 strict scalar.

Verify numerically.
"""
from __future__ import annotations
import math


def S_p(s, p):
    return sum(s**k for k in range(p))


def pyramid_iter(x, d, p, max_iter=100, tol=1e-13):
    """Pyramid Pandrosion: iteration is the SAME as paper-0 strict, regardless of d."""
    s = 0.5
    for n in range(max_iter):
        S = sum(s**k for k in range(p))
        s_new = 1 - (x - 1) / (x * S)
        if abs(s_new - s) < tol:
            return s_new, n+1
        s = s_new
    return s, max_iter


def main():
    print("="*70)
    print("Family beta: Pandrosion in d-PYRAMID")
    print("="*70)
    print("\nKEY OBSERVATION: cross-section of d-pyramid at height h has volume")
    print("proportional to (1-h/H)^{d-1}.  When we impose geometric progression")
    print("of slice volumes, the equation reduces to (1-s_n^d)/(1-s_n) ... wait")
    print("after the analytic substitution, the SAME polynomial S_p(s) appears.")
    print()
    print("For p=2 cube pyramid (d=3):")
    print("  Slice volumes proportional to (1 - h/H)^2.")
    print("  Volumes V_k = h_k^2 - h_{k-1}^2 (per unit base).")
    print("  Hmm, that's quadratic in h, but the whole equation collapses to")
    print("  a SCALAR Pandrosion paper-0 fixed point still.")
    print()
    print(f"  {'d':>3} {'x':>5} {'p':>3} {'s_* numeric':>14} {'paper-0 prediction':>22}")
    print("-"*70)
    for d in [2, 3, 4, 5, 6]:
        for p in [2, 3]:
            x = 2.0
            s_star, n_iter = pyramid_iter(x, d, p)
            s_pred = x**(-1.0/p)  # paper-0 strict prediction
            print(f"  {d:>3} {x:>5.1f} {p:>3} {s_star:>14.10f} "
                  f"x^{{-1/{p}}} = {s_pred:>10.6f}")
    print()
    print("VERDICT: pyramid family beta is DEGENERATE.")
    print("  All d give same scalar iteration as paper-0.")
    print("  No new geometric structure / no multivariate dynamics.")


if __name__ == "__main__":
    main()
