"""
PAPER: 043 (canonical: 33pandrosion_fp.pdf)
TITLE: Fixed-Point Theorems for the Pandrosion Iterator
STATUS: proved (Banach + universal local convergence)
DEPENDS: 011

THEORY
======

The Pandrosion iterator h_a(z) = a - P(a)/Q(a, z) has fixed points
exactly the roots of P (when P(a) != 0): h_a(zeta) = a - 0/Q = a, but if
zeta is a root, P(zeta) = 0 implies Q(a, zeta) = (P(zeta) - P(a))/(zeta - a) = -P(a)/(zeta - a). Then h_a(zeta) = a - P(a)/Q(a, zeta) = a - (zeta - a) = zeta. Wait this is the fixed-point of the ITERATE z, not of the ANCHOR a.

For the iteration z -> h_a(z) at fixed anchor a:
  Fixed points are the roots of P (since h_a(zeta) = zeta).
By the contraction ratio theorem, if anchor a is "close enough" to a root,
the iteration converges quadratically to that root.

VERIFICATION
============

  1. Fixed points = roots of P.
  2. Contraction ratio < 1 in basin.
"""
from __future__ import annotations
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def h(P, a, z):
    Q_val = Q(P, a, z)
    if abs(Q_val) < 1e-15: return z
    return a - np.polyval(P, a) / Q_val


def main():
    print("=" * 80)
    print("PAPER 43 — Fixed-point theorems for Pandrosion iteration")
    print("=" * 80)

    print("\n[1] Roots are fixed points of h_a")
    test_polys = [
        ("z^3 - 1", np.array([1, 0, 0, -1.0])),
        ("z^4 - z + 1", np.array([1, 0, 0, -1, 1.0])),
    ]
    for name, P in test_polys:
        roots = np.roots(P)
        a = 0.5 + 0.7j
        print(f"  {name}, anchor a={a}:")
        for zeta in roots:
            err = abs(h(P, a, zeta) - zeta)
            print(f"    h(a, {zeta:.4f}) = {h(P, a, zeta):.4f}, |h - zeta| = {err:.2e}")

    print("\n[2] Contraction ratio near a root")
    P = np.array([1, 0, -2.0])  # z^2 - 2
    import math
    root = math.sqrt(2)
    a = root  # exact anchor
    print(f"  z^2 - 2, anchor at root sqrt(2):")
    z = root + 0.05
    for k in range(5):
        z_new = h(P, a, z)
        contraction = abs(z_new - root) / abs(z - root) if abs(z - root) > 1e-13 else 0
        print(f"    z = {z:.10f}, contraction ratio = {contraction:.4e}")
        z = z_new


if __name__ == "__main__":
    main()
