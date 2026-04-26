"""
PAPER: 022 (canonical: 12pandrosion_casas_alvero.pdf)
TITLE: Casas-Alvero Conjecture: Initial Pandrosion-Spectrum Approach
STATUS: empirical scan (refined in papers 040, 048, 054, 106)
DEPENDS: 011

THEORY
======

CASAS-ALVERO (2001): If P in C[z] is monic of degree d and shares a root
with each of P', P'', ..., P^{(d-1)}, then P(z) = (z - alpha)^d.

PANDROSION-SPECTRUM REFORMULATION:
The Pandrosion-quotient spectrum at z is
  Spec_P(z) = {Q(alpha_j, z) : j = 1, ..., d}.
Casas-Alvero is equivalent to: this spectrum is injective in z at no point of
the joint vanishing locus of {P^(k)}_{k=1}^{d-1} unless P is a pure power.

VERIFICATION
============

  1. (z-1)^d satisfies Casas-Alvero condition vacuously (single root).
  2. Random perturbations: no counterexample.
  3. Dimension count: 2(d-1) constraints on d variables (over-determined).
"""
from __future__ import annotations
import math
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 22 — Casas-Alvero: initial spectrum approach")
    print("=" * 80)

    print("\n[1] (z-1)^d: P^(k)(1) = 0 for all k < d")
    for d in [3, 5, 8]:
        base = np.array([1.0])
        for _ in range(d):
            base = np.convolve(base, np.array([1, -1]))
        derivs = [base.copy()]
        P = base.copy()
        for _ in range(d - 1):
            P = np.polyder(P)
            derivs.append(P)
        # All evaluate to ~0 at z = 1
        evals = [abs(np.polyval(D, 1.0)) for D in derivs]
        print(f"  d={d}: |P(1)|, |P'(1)|, ..., |P^(d-1)(1)| = {[f'{e:.2e}' for e in evals[:5]]}")

    print("\n[2] Dimension count: 2(d-1) equations vs d unknowns")
    print(f"  {'d':>3} {'#equations':>12} {'#unknowns':>11} {'over-det by':>14}")
    for d in [4, 6, 8, 10, 12, 16]:
        eqs = 2 * (d - 1)
        vars_ = d
        print(f"  {d:>3} {eqs:>12} {vars_:>11} {eqs - vars_:>14}")

    print("\n[3] Random perturbation of (z-1)^d: still satisfies CA?")
    rng = np.random.default_rng(0)
    print(f"  {'d':>3} {'#trials':>9} {'#exact CA non-power':>22}")
    for d in [4, 6, 8]:
        n_test = 200
        n_ca_non_power = 0
        for _ in range(n_test):
            base = np.array([1.0])
            for _ in range(d):
                base = np.convolve(base, np.array([1, -1]))
            base = base.astype(complex)
            perturb = 0.01 * (rng.standard_normal(len(base)) + 1j * rng.standard_normal(len(base)))
            perturb[0] = 0
            P = base + perturb
            # Check: does P share a root with each P^(k)?
            satisfies = True
            roots = np.roots(P)
            for k in range(1, d):
                Pk = np.polyder(P, k)
                roots_k = np.roots(Pk)
                shared = any(min(abs(r - rk) for rk in roots_k) < 1e-3 for r in roots)
                if not shared:
                    satisfies = False
                    break
            if satisfies:
                # Pure power test: all roots equal?
                if max(abs(r - roots[0]) for r in roots) > 1e-3:
                    n_ca_non_power += 1
        print(f"  {d:>3} {n_test:>9} {n_ca_non_power:>22}")


if __name__ == "__main__":
    main()
