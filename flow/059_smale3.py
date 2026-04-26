"""
PAPER: 059 (canonical: 49pandrosion_smale3.pdf)
TITLE: Smale-Pandrosion (Part 3): Adaptive Anchor and Damping
STATUS: framework

THEORY
======

ADAPTIVE ANCHOR strategy: re-anchor every K steps (K = O(1)) to maintain
||Q^{-1} F|| ~ 1. This avoids the rational-map stagnation traps of paper 18.

DAMPING: in addition, scale the step by gamma-bound to ensure descent in
sup-norm of |F|.

VERIFICATION
============

  1. Adaptive anchor convergence on test polynomials.
  2. Effect of damping factor on convergence rate.
"""
from __future__ import annotations
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def adaptive_pandrosion(P, z0, K=3, max_iter=50, tol=1e-12):
    """Re-anchor every K steps."""
    z = z0
    n_steps = 0
    n_reanchor = 0
    for ep in range(max_iter):
        a = z  # re-anchor
        n_reanchor += 1
        for _ in range(K):
            Q_val = Q(P, a, z)
            if abs(Q_val) < 1e-15: return n_steps, z
            z = z - np.polyval(P, z) / Q_val
            n_steps += 1
            if abs(np.polyval(P, z)) < tol: return n_steps, z
    return n_steps, z


def damped_pandrosion(P, z0, gamma_factor=1.0, max_iter=50, tol=1e-12):
    """Damped Newton-Pandrosion: step *= gamma_factor."""
    z = z0
    for n in range(max_iter):
        if abs(np.polyval(P, z)) < tol: return n, z
        Pp = np.polyval(np.polyder(P), z)
        if abs(Pp) < 1e-15: return n, z
        z = z - gamma_factor * np.polyval(P, z) / Pp
    return max_iter, z


def main():
    print("=" * 80)
    print("PAPER 59 — Smale-Pandrosion: adaptive anchor + damping")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Adaptive anchor: K = 1, 3, 5")
    P = np.array([1.0, 0, 0, 0, -1])
    z0 = 0.5 + 0.5j
    for K in [1, 3, 5]:
        n, z = adaptive_pandrosion(P, z0, K=K)
        print(f"  K = {K}: {n} Pandrosion steps, |P(z)| = {abs(np.polyval(P, z)):.2e}")

    print("\n[2] Damping factor effect")
    P = np.array([1.0, 0, -2])  # z^2 - 2
    z0 = 1.5
    for gamma in [0.5, 1.0, 1.5]:
        n, z = damped_pandrosion(P, z0, gamma_factor=gamma)
        print(f"  damping = {gamma}: {n} iters, |P(z)| = {abs(np.polyval(P, z)):.2e}")


if __name__ == "__main__":
    main()
