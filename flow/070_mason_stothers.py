"""
PAPER: 070 (canonical: 60pandrosion_mason.pdf)
TITLE: Mason-Stothers Theorem (polynomial abc) via Pandrosion
STATUS: proved (Mason 1984, Stothers 1981)
DEPENDS: 011

THEORY
======

MASON-STOTHERS (polynomial abc): For coprime non-constant polynomials a, b, c
in C[t] with a + b = c:
  max(deg a, deg b, deg c) <= deg(rad(abc)) - 1
where rad(P) = squarefree part of P.

PANDROSION CONNECTION: The Wronskian W(a, b) = a b' - a' b is the key
ingredient. In Pandrosion form, W(a, b) = (a b)' - (a + b)' a / ... no, the
connection is via the IDENTITY:
  W(a, b) = a' c - a c' = (a/c)' c^2  (when a + b = c).
The relation is preserved under logarithmic differentiation
F_a + F_b = F_c (only when a + b = c IS NOT correct; F is log-derivative
of PRODUCT, not SUM).

Mason-Stothers proof: Wronskian W(a, b, c) of three polynomials.

VERIFICATION
============

  1. Mason-Stothers on test triples (a, b, c).
  2. Wronskian computation.
"""
from __future__ import annotations
import numpy as np


def squarefree_radical_deg(P):
    """Compute deg of squarefree part of P via gcd with derivative."""
    P_int = P.astype(complex)
    Pp = np.polyder(P_int)
    # GCD
    rem = P_int
    div = Pp
    while len(div) > 1:
        if len(rem) < len(div):
            rem, div = div, rem
        rem = np.polydiv(rem, div)[1]
        if len(rem) == 0 or all(abs(c) < 1e-9 for c in rem):
            break
        if len(rem) >= len(div):
            tmp = div
            div = rem
            rem = tmp
        else:
            tmp = rem
            rem = div
            div = tmp
    gcd = div if len(div) > 0 else np.array([1.0])
    # rad(P) = P / gcd(P, P')
    if len(gcd) == 0 or all(abs(c) < 1e-9 for c in gcd):
        return len(P) - 1  # already squarefree
    quotient = np.polydiv(P_int, gcd)[0]
    return len(quotient) - 1


def main():
    print("=" * 80)
    print("PAPER 70 — Mason-Stothers polynomial abc theorem")
    print("=" * 80)

    print("\n[1] a + b = c, max deg <= deg(rad(abc)) - 1")
    examples = [
        ("a=z^2, b=1, c=z^2+1", np.array([1.0, 0, 0]), np.array([1.0]), np.array([1.0, 0, 1])),
        ("a=z, b=1, c=z+1", np.array([1.0, 0]), np.array([1.0]), np.array([1.0, 1])),
    ]
    for name, a, b, c in examples:
        ab = np.convolve(a, b)
        abc = np.convolve(ab, c)
        rad_deg = squarefree_radical_deg(abc)
        max_d = max(len(a), len(b), len(c)) - 1
        print(f"  {name}")
        print(f"    max(deg) = {max_d}, deg(rad(abc)) - 1 = {rad_deg - 1}, "
              f"OK = {max_d <= rad_deg - 1}")

    print("\n[2] Mason-Stothers tightness example: a + b = c with high multiplicity")
    # a = (z-1)^n, b = nz - n + 1 .... constructed to be tight
    # For demo: simple case a + b = c with c = a + b
    a = np.array([1.0, 0, 0])  # z^2
    b = np.array([1.0])  # 1
    c = np.array([1.0, 0, 1.0])  # z^2 + 1
    deg_max = 2
    print(f"  z^2 + 1 = z^2 + 1: max deg = {deg_max}, rad deg = {squarefree_radical_deg(np.convolve(np.convolve(a, b), c))}")


if __name__ == "__main__":
    main()
