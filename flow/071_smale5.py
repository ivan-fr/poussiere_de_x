"""
PAPER: 071 (canonical: 61pandrosion_smale5.pdf)
TITLE: Smale-Pandrosion (Part 5): Final adjustments
STATUS: framework
DEPENDS: 053, 062

THEORY
======

Final tweaks to the Pandrosion-Smale machinery: refined gamma-bounds,
multi-scale anchor selection, hybrid univariate/multivariate convergence.

VERIFICATION
============

  1. Refined gamma estimate.
  2. Multi-scale anchors on test polynomials.
"""
from __future__ import annotations
import math
import numpy as np


def gamma_max(P, z, k_max=8):
    Pp_z = np.polyval(np.polyder(P), z)
    if abs(Pp_z) < 1e-15: return float('inf')
    d = len(P) - 1
    gammas = []
    for k in range(2, min(d, k_max) + 1):
        Pk = np.polyval(np.polyder(P, k), z)
        v = abs(Pk / (math.factorial(k) * Pp_z))
        if v > 0: gammas.append(v**(1.0 / (k - 1)))
    return max(gammas) if gammas else 0


def main():
    print("=" * 80)
    print("PAPER 71 — Smale-Pandrosion (Part 5): final adjustments")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Refined gamma at near-root vs far points")
    P = np.array([1.0, 0, 0, 0, -1])  # z^4 - 1
    for z in [1.0, 1.1, 1.5, 2.0, 5.0]:
        g = gamma_max(P, z)
        print(f"  z = {z}: gamma = {g:.4f}")

    print("\n[2] Multi-scale anchors: choose anchor with smallest gamma")
    P = rng.standard_normal(8); P[0] = 1.0
    anchors = [complex(rng.standard_normal(), rng.standard_normal()) for _ in range(5)]
    z = 0.5
    print(f"  Test point z = {z}, gamma values at anchors:")
    for a in anchors:
        g = gamma_max(P, a)
        print(f"    anchor {a:.3f}: gamma = {g:.4f}")


if __name__ == "__main__":
    main()
