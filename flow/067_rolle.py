"""
PAPER: 067 (canonical: 57pandrosion_rolle.pdf)
TITLE: Rolle's Theorem via the Pandrosion Field
STATUS: proved (classical, Pandrosion-form interpretation)
DEPENDS: 011

THEORY
======

ROLLE: For real-rooted P with simple roots alpha_1 < alpha_2 < ... < alpha_d,
between any two consecutive roots there's at least one root of P'.

PANDROSION INTERPRETATION: The Pandrosion field F_P = P'/P has poles at
{alpha_j} (negative residue jump). Between consecutive roots, F_P decreases
monotonically from +infty (just right of alpha_j) to -infty (just left of
alpha_{j+1}), crossing 0 at a root of P'.

VERIFICATION
============

  1. F_P sign change between consecutive roots.
  2. Critical point in each interval (alpha_j, alpha_{j+1}).
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 67 — Rolle's theorem via Pandrosion field")
    print("=" * 80)

    P = np.poly([1.0, 2.0, 4.0, 7.0])  # roots 1, 2, 4, 7
    roots = sorted(np.roots(P).real)
    crits = sorted(np.roots(np.polyder(P)).real)
    print(f"\n  Roots of P: {roots}")
    print(f"  Critical points P': {crits}")

    print("\n[1] One critical point in each (alpha_j, alpha_{j+1})")
    for j in range(len(roots) - 1):
        a, b = roots[j], roots[j+1]
        cs_in = [c for c in crits if a < c < b]
        print(f"  ({a}, {b}): {len(cs_in)} critical point(s) at {cs_in}")

    print("\n[2] Pandrosion field sign change between roots")
    Pp = np.polyder(P)
    for j in range(len(roots) - 1):
        a, b = roots[j], roots[j+1]
        # Sample F_P just inside the interval
        eps = 0.01 * (b - a)
        F_left = np.polyval(Pp, a + eps) / np.polyval(P, a + eps)
        F_right = np.polyval(Pp, b - eps) / np.polyval(P, b - eps)
        print(f"  F_P({a + eps:.3f}) = {F_left:+.4f}, F_P({b - eps:.3f}) = {F_right:+.4f}, "
              f"sign change: {F_left * F_right < 0}")


if __name__ == "__main__":
    main()
