"""
PAPER: 033 (canonical: 23pandrosion_energy.pdf)
TITLE: Strict Positivity of Pandrosion Energy
STATUS: proved (paper 23 of legacy: E_P > 0 strictly on R for simple real roots)
DEPENDS: 030

THEORY
======

If P has all simple real roots, then E_P(x) = sum Q(alpha_k, x)^2 > 0
strictly on R (no real zeros of E_P).

Proof: at root alpha_m, Q(alpha_k, alpha_m) = 0 for k != m (since the product
over j != k contains the factor (alpha_m - alpha_m) = 0), but
Q(alpha_m, alpha_m) = P'(alpha_m) != 0. So E_P(alpha_m) = P'(alpha_m)^2 > 0.
For x != any alpha_j, all Q(alpha_k, x) sum to >0 (cannot all vanish).

LANDSCAPE: E_P has degree 2(d-1), strictly positive on R, hence d - 1
conjugate pairs of complex roots.

VERIFICATION
============

  1. E_P > 0 strictly for real-rooted P.
  2. E_P at roots = P'^2.
  3. E_P has 2(d-1) complex conjugate roots.
"""
from __future__ import annotations
import numpy as np


def E_polynomial(P):
    """E_P(z) = (P')^2 - P P'' as polynomial."""
    Pp = np.polyder(P)
    Ppp = np.polyder(Pp)
    # (P')^2
    Pp_sq = np.convolve(Pp, Pp)
    # P P''
    PPpp = np.convolve(P, Ppp)
    # Align degrees
    max_len = max(len(Pp_sq), len(PPpp))
    a = np.concatenate([np.zeros(max_len - len(Pp_sq)), Pp_sq])
    b = np.concatenate([np.zeros(max_len - len(PPpp)), PPpp])
    return a - b


def main():
    print("=" * 80)
    print("PAPER 33 — Strict positivity of Pandrosion energy")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] E_P > 0 on R for real-rooted P")
    for d in [3, 5, 7, 9]:
        roots = sorted(rng.uniform(-2, 2, d))
        P = np.poly(roots)
        E_poly = E_polynomial(P)
        # Sample E_P on R
        xs = np.linspace(-3, 3, 500)
        E_vals = np.polyval(E_poly, xs)
        print(f"  d={d}: deg E_P = {len(E_poly)-1}, "
              f"min E_P on x in [-3,3] = {E_vals.min():.4e}")

    print("\n[2] At root alpha: E_P(alpha) = (P'(alpha))^2")
    for d in [3, 5]:
        roots = sorted(rng.uniform(-1, 1, d))
        P = np.poly(roots)
        E_poly = E_polynomial(P)
        Pp = np.polyder(P)
        for ak in roots:
            E_at = np.polyval(E_poly, ak)
            Pp_sq = np.polyval(Pp, ak)**2
            print(f"  d={d}, alpha={ak:.4f}: E={E_at:.4e}, P'^2={Pp_sq:.4e}, diff={abs(E_at-Pp_sq):.2e}")

    print("\n[3] E_P has 2(d-1) complex conjugate roots")
    for d in [3, 4, 5, 6]:
        roots = sorted(rng.uniform(-1, 1, d))
        P = np.poly(roots)
        E_poly = E_polynomial(P)
        E_roots = np.roots(E_poly)
        n_complex = sum(1 for r in E_roots if abs(r.imag) > 1e-9)
        print(f"  d={d}: deg E_P = {2*(d-1)}, complex roots = {n_complex}, "
              f"real roots = {len(E_roots) - n_complex}")


if __name__ == "__main__":
    main()
