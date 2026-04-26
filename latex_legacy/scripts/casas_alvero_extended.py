"""
Extended Casas-Alvero attack: scan large degrees with proper tolerance handling.

Conjecture: P shares root with each P^(k), 1 <= k <= d-1, ⟹ P = (z-alpha)^d.

Key idea: rather than perturb (z-1)^d (which always satisfies CA at perturbed
levels), we directly seek polynomials matching the CA condition that AREN'T
pure powers.

Approach:
  1. For each k = 1, ..., d-1, pick a random alpha_k.
  2. Construct P satisfying P(alpha_k) = 0 AND P^(k)(alpha_k) = 0 for each k.
     This is d-1 constraints on the d coefficients of P (modulo monic).
  3. If a solution exists with NOT all alpha_k equal, we have a counterexample.
"""
from __future__ import annotations
import math
import numpy as np


def derivatives_all(coeffs):
    """Return [P, P', P'', ..., P^(d-1)]."""
    derivs = [coeffs.copy()]
    P = coeffs.copy()
    d = len(coeffs) - 1
    for _ in range(d - 1):
        P = np.polyder(P)
        derivs.append(P)
    return derivs


def evaluates_to_zero_jointly(coeffs, alphas, tol=1e-6):
    """Test: does P satisfy P(alpha_0)=0 and P^(k)(alpha_k)=0 for each k?

    alphas: list of d-1 complex numbers, alphas[k-1] is the shared root with P^(k).
    """
    d = len(coeffs) - 1
    derivs = derivatives_all(coeffs)
    for k in range(1, d):
        val = abs(np.polyval(derivs[k], alphas[k-1]))
        # Also need P(alpha_k) = 0: P shares root with P^(k)
        Pval = abs(np.polyval(coeffs, alphas[k-1]))
        if val > tol or Pval > tol:
            return False
    return True


def construct_CA_candidate(d, alphas, rng):
    """Construct a monic poly of degree d such that for each k,
    P(alpha_k) = 0 and P^(k)(alpha_k) = 0.

    This is 2(d-1) constraints, more than the d coefficients available
    (modulo monic = d-1 free coeffs). So generically NO solution exists,
    confirming Casas-Alvero.

    But: if alpha_1 = alpha_2 = ... = alpha_{d-1} = alpha, then constraints
    collapse to "P has root of multiplicity d at alpha", giving P = (z-alpha)^d.

    Test: for distinct alphas, no solution should exist (Casas-Alvero).
    """
    # Linear system: vector of coefs [a_0, a_1, ..., a_{d-1}] (a_d = 1, monic)
    # Constraints:
    #   For each k = 1, ..., d-1:
    #     P(alpha_k) = 0  AND  P^(k)(alpha_k) = 0
    #
    # Total: 2(d-1) equations in d unknowns. Over-determined for d >= 3.
    # If solution exists with distinct alphas, that's a counterexample.

    n_constraints = 2 * (d - 1)
    n_unknowns = d  # a_0, ..., a_{d-1}, with a_d = 1 fixed

    A = np.zeros((n_constraints, n_unknowns), dtype=complex)
    b = np.zeros(n_constraints, dtype=complex)

    row = 0
    for k in range(1, d):
        alpha_k = alphas[k - 1]
        # Constraint 1: P(alpha_k) = 0
        # P(z) = sum_{i=0}^{d} a_i z^i, with a_d = 1
        # P(alpha_k) = sum_{i=0}^{d} a_i alpha_k^i = 0
        for i in range(d):
            A[row, i] = alpha_k ** i
        b[row] = -alpha_k ** d  # since a_d = 1
        row += 1

        # Constraint 2: P^(k)(alpha_k) = 0
        # P^(k)(z) = sum_{i=k}^{d} a_i * i*(i-1)*...*(i-k+1) * z^{i-k}
        for i in range(k, d):
            coef = math.factorial(i) // math.factorial(i - k)
            A[row, i] = coef * (alpha_k ** (i - k))
        # a_d term:
        coef_d = math.factorial(d) // math.factorial(d - k)
        b[row] = -coef_d * (alpha_k ** (d - k))
        row += 1

    # Solve via least squares (over-determined)
    sol, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)
    res_norm = float(np.linalg.norm(A @ sol - b))
    coeffs = np.array([1.0+0j] + list(sol[::-1]))  # high to low order
    return coeffs, res_norm


def search_CA_counterexample(d, n_trials, rng, alpha_range=2.0):
    """Search for distinct-alpha solutions of the CA constraints."""
    successes = 0
    smallest_residual = float('inf')
    smallest_residual_polly = None
    smallest_residual_alphas = None
    for trial in range(n_trials):
        # Generate d-1 distinct random alpha values
        alphas = rng.standard_normal(d-1) + 1j * rng.standard_normal(d-1)
        alphas = alphas * alpha_range
        # Make sure they're not all equal
        if np.std(alphas) < 1e-3:
            continue
        coeffs, res = construct_CA_candidate(d, alphas, rng)
        if res < smallest_residual:
            smallest_residual = res
            smallest_residual_polly = coeffs
            smallest_residual_alphas = alphas
        if res < 1e-6:
            successes += 1
    return successes, smallest_residual, smallest_residual_polly, smallest_residual_alphas


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 95, flush=True)
    print("CASAS-ALVERO ATTACK: solve constraints directly", flush=True)
    print("=" * 95, flush=True)
    print("\nStrategy: construct P satisfying P(alpha_k) = P^(k)(alpha_k) = 0 for each k.")
    print("If solution exists with distinct alphas, that's a counterexample to CA.")

    print(f"\n{'d':>3} {'#trials':>9} {'#exact':>10} {'min residual':>15} {'status':>20}",
          flush=True)
    rng = np.random.default_rng(20260603)
    open_degrees = [8, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24]
    for d in open_degrees:
        n_trials = 500 if d < 16 else 200
        successes, min_res, _, _ = search_CA_counterexample(d, n_trials, rng)
        status = "PASSES (no counterexample)" if successes == 0 else f"FOUND ({successes})"
        print(f"{d:>3} {n_trials:>9} {successes:>10} {min_res:>15.4e} {status:>20}",
              flush=True)

    # Sanity: pure power should give residual = 0 with all alphas equal
    print("\n[Sanity] (z-1)^d with all alphas = 1:")
    for d in [5, 8, 12]:
        alphas = np.ones(d-1, dtype=complex)
        coeffs, res = construct_CA_candidate(d, alphas, rng)
        print(f"  d = {d}: residual = {res:.4e}, polynomial coeffs[0] = {coeffs[0]}, "
              f"coeffs[-1] = {coeffs[-1]:.4f}",
              flush=True)

    # Test: the system is over-determined (2(d-1) eqs > d unknowns).
    # The "no counterexample" finding is then expected from generic dimension counting,
    # which is essentially Casas-Alvero (the generic case).
    # The deep open question: are there special non-generic solutions for non-prime-power d?
    print("\n[Dimension count]")
    for d in [4, 8, 10, 12, 14, 15, 16]:
        n_eqs = 2 * (d - 1)
        n_vars = d
        print(f"  d = {d}: equations = {n_eqs}, variables = {n_vars}, "
              f"over-determined by {n_eqs - n_vars}",
              flush=True)


if __name__ == "__main__":
    main()
