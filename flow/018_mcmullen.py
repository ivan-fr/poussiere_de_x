"""
PAPER: 018 (canonical: 08pandrosion_mcmullen.pdf)
TITLE: McMullen's Barrier and the Pandrosion Bypass
STATUS: framework (McMullen 1987 obstruction + non-holomorphic bypass)
DEPENDS: 007, 017

THEORY
======

McMULLEN'S THEOREM (1987):
There exists no purely-rational iterator (a single rational map z -> r(z)
applied uniformly) that converges to a root of every degree-d polynomial
for d >= 4.

PANDROSION BYPASS:
A fixed-anchor Pandrosion iteration is a rational map, so it inherits
McMullen's obstruction. The bypass is to:
  - Use ADAPTIVE anchor (re-anchoring), making the iteration not a single
    rational map but a sequence of them.
  - Combine with NON-HOLOMORPHIC fallback (|P(z)| comparisons) — paper 7.

The non-holomorphic check is a real-valued residual test; in BSS model,
branching on inequalities costs O(1) operations.

VERIFICATION
============

  1. Demonstrate McMullen's barrier: a fixed-anchor T_2 stagnates on z^4 + 1.
  2. Adaptive anchor + Armijo fallback escapes the trap.
"""
from __future__ import annotations
import math
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def fixed_anchor_T2(P, a, z, max_iter=30):
    """Fixed-anchor T_2 (Aitken) — pure rational iteration, McMullen-blocked."""
    orbit = [z]
    for _ in range(max_iter):
        s1 = z - np.polyval(P, z) / Q(P, a, z)
        s2 = s1 - np.polyval(P, s1) / Q(P, a, s1)
        denom = s2 - 2*s1 + z
        if abs(denom) < 1e-15: break
        z_new = z - (s1 - z)**2 / denom
        if abs(z_new - z) < 1e-13: break
        z = z_new
        orbit.append(z)
    return orbit


def adaptive_T2_with_fallback(P, z0, max_iter=50, eta=0.95):
    """Re-anchor each step + Armijo fallback if no descent."""
    z = z0
    orbit = [z]
    Pz = abs(np.polyval(P, z))
    for _ in range(max_iter):
        if Pz < 1e-13: break
        # Adaptive: anchor = z
        a = z
        s1 = z - np.polyval(P, z) / Q(P, a, z)
        s2 = s1 - np.polyval(P, s1) / Q(P, a, s1)
        denom = s2 - 2*s1 + z
        z_cand = s2 if abs(denom) < 1e-15 else z - (s1 - z)**2 / denom
        Pz_cand = abs(np.polyval(P, z_cand))
        if Pz_cand <= eta * Pz:
            z = z_cand
            Pz = Pz_cand
            orbit.append(z)
            continue
        # Fallback: Armijo on Newton direction
        Pp_z = np.polyval(np.polyder(P), z)
        if abs(Pp_z) < 1e-15: break
        direction = np.polyval(P, z) / Pp_z
        success = False
        for j in range(20):
            tau = 2**(-j)
            z_new = z - tau * direction
            if abs(np.polyval(P, z_new)) <= (1 - 0.1 * tau) * Pz:
                z = z_new
                Pz = abs(np.polyval(P, z))
                orbit.append(z)
                success = True
                break
        if not success: break
    return orbit


def main():
    print("=" * 80)
    print("PAPER 18 — McMullen's barrier and Pandrosion bypass")
    print("=" * 80)

    # Polynomial z^4 + 1 — known McMullen-style trap candidate
    P = np.array([1.0, 0, 0, 0, 1.0])
    print(f"\n  P = z^4 + 1, roots: {np.roots(P)}")

    print("\n[1] Fixed-anchor T_2 (rational map) — McMullen-blocked")
    z0 = 0.7  # generic start
    a = 0.5  # poor anchor
    orbit = fixed_anchor_T2(P, a, z0, max_iter=20)
    final = abs(np.polyval(P, orbit[-1]))
    print(f"  Anchor a = {a}, z_0 = {z0}: {len(orbit)} iters, |P(z_final)| = {final:.4e}")

    print("\n[2] Adaptive T_2 + Armijo fallback (McMullen bypass)")
    orbit2 = adaptive_T2_with_fallback(P, z0, max_iter=50)
    final2 = abs(np.polyval(P, orbit2[-1]))
    print(f"  z_0 = {z0}: {len(orbit2)} iters, |P(z_final)| = {final2:.4e}")

    print("\n[3] Many starts: convergence rate")
    rng = np.random.default_rng(0)
    n_total = 50
    n_conv_fixed = 0
    n_conv_adapt = 0
    for _ in range(n_total):
        z0 = complex(rng.uniform(-2, 2), rng.uniform(-2, 2))
        a = complex(rng.uniform(-2, 2), rng.uniform(-2, 2))
        # Fixed anchor
        orb = fixed_anchor_T2(P, a, z0, max_iter=20)
        if abs(np.polyval(P, orb[-1])) < 1e-6: n_conv_fixed += 1
        # Adaptive
        orb2 = adaptive_T2_with_fallback(P, z0, max_iter=50)
        if abs(np.polyval(P, orb2[-1])) < 1e-6: n_conv_adapt += 1
    print(f"  Fixed-anchor convergence: {n_conv_fixed}/{n_total}")
    print(f"  Adaptive + fallback:      {n_conv_adapt}/{n_total}")


if __name__ == "__main__":
    main()
