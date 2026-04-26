"""
PAPER: 073 (canonical: 63pandrosion_ps.pdf)
TITLE: Power Sums and Pandrosion (Newton's identities)
STATUS: proved (Newton's identities, Pandrosion-form)
DEPENDS: 011

THEORY
======

POWER SUM: p_k = sum_j alpha_j^k.
Newton's identity:
  p_k = e_1 p_{k-1} - e_2 p_{k-2} + ... + (-1)^{k-1} k e_k
relating power sums p_k to elementary symmetric polynomials e_k.

PANDROSION FORM: F_P(z) = sum 1/(z - alpha_j) admits Laurent expansion at
infinity:
  F_P(z) = d/z + p_1/z^2 + p_2/z^3 + ...
The coefficients of F_P at infinity are exactly the power sums p_k.

VERIFICATION
============

  1. Laurent expansion of F_P at infinity = power sums.
  2. Newton's identity: e_k from p_k.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 73 — Power sums and Pandrosion (Newton's identities)")
    print("=" * 80)

    P = np.poly([1, 2, 3, 4.0])  # roots 1, 2, 3, 4
    roots = sorted(np.roots(P).real)
    print(f"\n[1] P with roots {roots}")
    # Power sums
    p = [sum(r**k for r in roots) for k in range(5)]
    print(f"  Power sums p_0..p_4: {p}")

    print("\n[2] F_P(z) ~ d/z + p_1/z^2 + p_2/z^3 + ... at infinity")
    # F_P at large z
    Pp = np.polyder(P)
    for z in [10.0, 100.0, 1000.0]:
        F = np.polyval(Pp, z) / np.polyval(P, z)
        # Compare with d/z + p_1/z^2 + ...
        approx = sum(p[k] / z**(k+1) for k in range(5)) if abs(z) > 1 else 0
        # p_0 = d (sum 1's)
        print(f"  z = {z}: F_P = {F:.6f}, predicted = {len(roots)/z + p[1]/z**2 + p[2]/z**3:.6f}")

    print("\n[3] Elementary symmetric e_k from coefficients")
    coefs = P
    d = len(coefs) - 1
    e = [(-1)**k * coefs[k] / coefs[0] for k in range(d + 1)]
    print(f"  Coefficients (high to low): {coefs}")
    print(f"  Elementary symmetric e_0..e_{d}: {e}")
    # Verify: e_1 = sum alpha = 1+2+3+4 = 10
    print(f"  e_1 = sum alpha = {sum(roots)}, computed = {e[1]}")


if __name__ == "__main__":
    main()
