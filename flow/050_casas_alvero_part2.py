"""
PAPER: 050 (canonical: 40pandrosion_ca.pdf)
TITLE: Casas-Alvero (Part 2): Pandrosion-Spectrum Refinement
STATUS: empirical (refined further in paper 106)
DEPENDS: 022

THEORY
======

Continuation of paper 22 (Casas-Alvero): the Pandrosion-spectrum
{Q(alpha_j, z)}_{j=1}^d is studied as a discrete object.

Casas-Alvero is equivalent to: this spectrum is non-injective in z only at
the joint vanishing locus of {P^(k)} and only when P is a pure power.

Strengthening: the system of constraints
  P(alpha_k) = P^(k)(alpha_k) = 0,  k = 1, ..., d-1
is over-determined by d-2 equations (paper 106 structural attack).

VERIFICATION
============

  1. Spectrum Q(alpha_k, z) at various z.
  2. Random scan: no counterexamples to Casas-Alvero.
"""
from __future__ import annotations
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def main():
    print("=" * 80)
    print("PAPER 50 — Casas-Alvero (Part 2): spectrum refinement")
    print("=" * 80)

    P = np.array([1, -4, 6, -4, 1.0])  # (z-1)^4
    roots = np.roots(P)
    print(f"\n[1] P = (z-1)^4: roots {roots}")
    for z in [0.0, 0.5, 1.5, 2.0]:
        spec = [Q(P, ak, z) for ak in roots]
        print(f"  z = {z}: spectrum = {[f'{s:.4f}' for s in spec]}")

    print("\n[2] Generic d=5 polynomial: spectrum at random z")
    rng = np.random.default_rng(0)
    P = rng.standard_normal(6); P[0] = 1.0
    roots = np.roots(P)
    z = 0.7 + 0.3j
    spec = [Q(P, ak, z) for ak in roots]
    print(f"  P degree 5, z = {z}: |spectrum| = {[f'{abs(s):.4f}' for s in spec]}")

    print("\n[3] Pure-power detection: spectrum constant iff P = (z-alpha)^d")
    for name, P in [("(z-1)^4", np.array([1, -4, 6, -4, 1.0])),
                     ("z^4 - 1", np.array([1.0, 0, 0, 0, -1])),
                     ("z^4 - z + 1", np.array([1.0, 0, 0, -1, 1]))]:
        roots = np.roots(P)
        z = 0.5
        spec = [Q(P, ak, z) for ak in roots]
        std = float(np.std(np.abs(spec)))
        print(f"  {name}: std(|spec|) = {std:.4f}  ({'pure power' if std < 1e-6 else 'mixed'})")


if __name__ == "__main__":
    main()
