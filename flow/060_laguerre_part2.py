"""
PAPER: 060 (canonical: 50pandrosion_laguerre.pdf)
TITLE: Laguerre's Theorem (Part 2): Disk Stability
STATUS: proved (Laguerre 1880, Pandrosion-form)
DEPENDS: 030

THEORY
======

LAGUERRE'S THEOREM: If P has all roots in disk D = {|z| <= R}, then so do
all derivatives P', P'', ..., P^(d-1).

Combined with Marden (paper 25) for d=3: critical points are foci of
inscribed Steiner ellipse.

PANDROSION FORM: For roots in D, the convex hull is in D. Gauss-Lucas
gives critical points in convex hull, hence in D. Iterating gives all
derivatives' roots in D.

VERIFICATION
============

  1. Roots of P in disk -> roots of P^(k) in disk.
  2. Tightness for cyclotomic z^d - 1 (boundary).
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 60 — Laguerre theorem: disk stability of derivatives")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Random polys with roots in disk: derivatives stay in disk")
    print(f"  {'d':>4} {'R':>4} {'#trials':>9} {'#OK (all derivs in D)':>22}")
    for d in [4, 6, 8, 10]:
        R = 1.0
        n_test = 30
        n_ok = 0
        for _ in range(n_test):
            angles = rng.uniform(0, 2*np.pi, d)
            radii = rng.uniform(0, R, d)
            roots = radii * np.exp(1j * angles)
            P = np.poly(roots)
            ok = True
            for k in range(1, d):
                Pk = np.polyder(P, k)
                roots_k = np.roots(Pk)
                if any(abs(r) > R + 1e-9 for r in roots_k):
                    ok = False
                    break
            if ok: n_ok += 1
        print(f"  {d:>4} {R:>4.1f} {n_test:>9} {n_ok:>22}")

    print("\n[2] z^d - 1: boundary case (roots on |z| = 1)")
    for d in [4, 6, 8]:
        P = np.array([1.0] + [0]*(d-1) + [-1])
        Pp = np.polyder(P)
        roots_p = np.roots(Pp)
        max_r = max(abs(r) for r in roots_p)
        print(f"  z^{d} - 1: max|root P'| = {max_r:.4f} (roots P' should be at 0)")


if __name__ == "__main__":
    main()
