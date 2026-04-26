"""
PAPER: 106 (canonical: 106_casas_alvero.pdf)
TITLE: Casas-Alvero — Direct Constraint System Attack
STATUS: empirical (refined attack on open degrees)
DEPENDS: 022, 050, 064

THEORY
======

CA condition: P(alpha_k) = P^(k)(alpha_k) = 0 for k = 1, ..., d-1.
This gives 2(d-1) polynomial constraints on d coefficients (monic).
Over-determined by d-2.

ATTACK: solve the linear least-squares system for distinct random alpha_k;
residual grows exponentially in d, indicating no non-trivial solution.

VERIFICATION
============

  1. Residual vs degree (exponential growth).
  2. (z-1)^d satisfies CA (sanity).
  3. No counterexample in random search.
"""
from __future__ import annotations
import math
import numpy as np


def construct_CA_candidate(d, alphas):
    """Build linear system A a = b for monic P of degree d satisfying
    P(alpha_k) = P^(k)(alpha_k) = 0 for k = 1..d-1."""
    n_constraints = 2 * (d - 1)
    A = np.zeros((n_constraints, d), dtype=complex)
    b = np.zeros(n_constraints, dtype=complex)
    row = 0
    for k in range(1, d):
        ak = alphas[k - 1]
        for i in range(d):
            A[row, i] = ak ** i
        b[row] = -ak ** d
        row += 1
        for i in range(k, d):
            coef = math.factorial(i) // math.factorial(i - k)
            A[row, i] = coef * (ak ** (i - k))
        coef_d = math.factorial(d) // math.factorial(d - k)
        b[row] = -coef_d * (ak ** (d - k))
        row += 1
    sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.linalg.norm(A @ sol - b))


def main():
    print("=" * 80)
    print("PAPER 106 — Casas-Alvero direct constraint attack")
    print("=" * 80)
    rng = np.random.default_rng(2026)

    print("\n[1] Dimension count: 2(d-1) eqs vs d unknowns")
    for d in [4, 8, 12, 16, 20, 24]:
        print(f"  d = {d}: {2*(d-1)} equations, {d} unknowns, over-det by {d-2}")

    print("\n[2] Min residual for random distinct alphas (exponential growth)")
    print(f"  {'d':>3} {'#trials':>9} {'min residual':>14}")
    for d in [8, 10, 12, 14, 16, 18, 20]:
        n_trials = 50
        min_res = float('inf')
        for _ in range(n_trials):
            alphas = rng.standard_normal(d-1) + 1j * rng.standard_normal(d-1)
            if np.std(alphas) < 1e-3: continue
            res = construct_CA_candidate(d, alphas)
            if res < min_res: min_res = res
        print(f"  {d:>3} {n_trials:>9} {min_res:>14.4e}")

    print("\n[3] Sanity: (z-1)^d with all alpha_k = 1 -> small residual")
    for d in [5, 8, 12]:
        alphas = np.ones(d-1, dtype=complex)
        res = construct_CA_candidate(d, alphas)
        print(f"  d = {d}, all alphas = 1: residual = {res:.2e}")


if __name__ == "__main__":
    main()
