"""
PAPER: 078 (canonical: 68pandrosion_borwein.pdf)
TITLE: Borwein Conjectures and Pandrosion Positivity
STATUS: framework (Borwein conjectures partially proven)
DEPENDS: 011

THEORY
======

BORWEIN CONJECTURES (1990s): For Pi_n(q) = prod_{k=0}^{n-1}(1 - q^{3k+1})(1 - q^{3k+2}),
the coefficients in the q-expansion alternate in sign in a specific pattern
(Borwein A, B, C conjectures). Some cases proven (Wang, Krattenthaler), but
general open.

PANDROSION CONNECTION: Coefficient positivity ~ Pólya-Schur preserver
applied to specific generating functions.

VERIFICATION
============

  1. Compute first coefficients of Pi_n(q).
  2. Check sign pattern.
"""
from __future__ import annotations
import numpy as np


def Pi_n(n, max_deg):
    """Pi_n(q) = prod_{k=0}^{n-1}(1 - q^{3k+1})(1 - q^{3k+2}), truncated to deg max_deg."""
    P = np.zeros(max_deg + 1)
    P[0] = 1.0  # constant 1
    P = P[::-1]  # high-to-low convention for numpy
    # Multiply factors
    for k in range(n):
        for r in [3*k + 1, 3*k + 2]:
            if r > max_deg: continue
            # Factor (1 - q^r) in low-to-high: [1, 0, ..., 0, -1] (length r+1)
            factor_low = [0]*(r+1)
            factor_low[0] = 1
            factor_low[r] = -1
            factor = factor_low[::-1]  # high-to-low
            P = np.convolve(P, factor)
            # Trim to max_deg + 1
            if len(P) > max_deg + 1:
                P = P[-(max_deg + 1):]
    return P


def main():
    print("=" * 80)
    print("PAPER 78 — Borwein conjectures")
    print("=" * 80)

    for n in [1, 2, 3]:
        P = Pi_n(n, 25)
        # Reverse to low-to-high for display
        c_low = P[::-1]
        # First few non-zero coefficients
        print(f"\n  Pi_{n}(q) low-deg coeffs: {[int(c) for c in c_low[:12]]}")

    print("\n[2] Borwein A conjecture: signs of c_k for Pi_n alternate by k mod 3")
    P = Pi_n(2, 30)
    c_low = [int(c) for c in P[::-1]]
    # Group by k mod 3
    print(f"  k=0,3,6,...: {[c_low[k] for k in range(0, len(c_low), 3)][:6]}")
    print(f"  k=1,4,7,...: {[c_low[k] for k in range(1, len(c_low), 3)][:6]}")
    print(f"  k=2,5,8,...: {[c_low[k] for k in range(2, len(c_low), 3)][:6]}")


if __name__ == "__main__":
    main()
