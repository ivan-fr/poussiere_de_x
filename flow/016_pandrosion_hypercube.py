"""
PAPER: 016 (alternative geometry: d-hypercube)
TITLE: Pandrosion in d-HYPERCUBE (family epsilon)
STATUS: empirical exploration; expected to be SEPARABLE

CONSTRUCTION
============
d-hypercube of sides 1 x 1 x ... x 1 x x  (one expanded direction).
Insert p-1 hyperplanes parallel to each axis, parameterised by s_i.
The d insertions are independent: each direction is orthogonal.

POLYNOMIAL:
    S_p^{cube, d}(s_1, ..., s_d) = prod_{i=1}^{d} S_p(s_i)
                                 = prod_i (1 + s_i + s_i^2 + ... + s_i^{p-1}).

ITERATION:
    s_{i, n+1} = 1 - (x_i - 1) / (x_i * S_p^{cube,d}(s_n))

CORRECTED VERDICT (after numerical test):
The product structure S = prod_i S_p(s_i) does NOT separate.
At the symmetric fixed point x_i = x, s_i = s, we get the equation
    x * S_p(s)^d * (1 - s) = x - 1.
For p=2, this becomes a polynomial of degree d+1 in s, with a new
algebraic fixed point that is NOT x^{-1/p}.

Example: d=2, p=2, x=2: 2 (1+s)^2 (1-s) = 1 -> 2s^3 + 2s^2 - 2s - 1 = 0
         s_* approx 0.8546 (NEW algebraic number, not 1/sqrt(2)).

For asymmetric x_i, the system is genuinely coupled multivariate.
The hypercube/product structure is THEREFORE a NEW geometry, not a
trivial product of paper-0 instances.

Verify numerically.
"""
from __future__ import annotations
import math


def S_p(s, p):
    return sum(s**k for k in range(p))


def S_hypercube(s_vec, p):
    """S_cube,d(s_1, ..., s_d) = prod_i S_p(s_i)."""
    total = 1.0+0.0j
    for s in s_vec:
        total *= sum(complex(s)**k for k in range(p))
    return total


def iterate_hypercube(x_vec, d, p, max_iter=200, tol=1e-13):
    n = d
    s = [complex(0.5)] * n
    history = [tuple(s)]
    for step in range(max_iter):
        Spsn = S_hypercube(s, p)
        s_new = [1 - (x_vec[i]-1)/(x_vec[i]*Spsn) for i in range(n)]
        if max(abs(s_new[i]-s[i]) for i in range(n)) < tol:
            return tuple(s_new), step+1, history
        s = s_new
        history.append(tuple(s))
    return tuple(s), max_iter, history


def main():
    print("="*70)
    print("Family epsilon: Pandrosion in d-HYPERCUBE")
    print("="*70)
    print()
    print("Polynomial: S_cube,d(s_1, ..., s_d) = prod_i S_p(s_i)  (SEPARABLE)")
    print()
    # Symmetric case: NEW algebraic fixed points
    print("Symmetric x_i = x = 2: NEW fixed points (NOT paper-0).")
    print(f"  Equation: x * S_p(s)^d * (1-s) = x - 1")
    print(f"  {'d':>3} {'p':>3} {'s_*':>14} {'cubic/poly residual':>22} {'iter':>5}")
    print("-"*70)
    import numpy as np
    for d in [2, 3, 4, 5, 6]:
        for p in [2, 3]:
            x = 2.0
            s_final, n_iter, _ = iterate_hypercube((x,)*d, d, p)
            s_star = s_final[0].real
            # Verify numerically: check x * S_p(s_*)^d * (1 - s_*) = x - 1
            S_val = sum(s_star**k for k in range(p))
            residual = abs(x * S_val**d * (1 - s_star) - (x - 1))
            print(f"  {d:>3} {p:>3} {s_star:>14.10f} {residual:>22.2e} {n_iter:>5}")

    # Specific cubic verification at d=2, p=2, x=2
    print()
    print("Special case d=2, p=2, x=2: cubic 2s^3 + 2s^2 - 2s - 1 = 0")
    coeffs = [2, 2, -2, -1]
    roots = np.roots(coeffs)
    print(f"  numeric roots: {roots}")
    real_root = [r.real for r in roots if abs(r.imag) < 1e-9 and r.real > 0]
    print(f"  positive real root: {real_root[0]:.15f}")

    print()
    print("Asymmetric x_i different per direction:")
    print(f"  d=3, x=(2, 3, 5), p=2")
    s_final, n_iter, _ = iterate_hypercube((2.0, 3.0, 5.0), 3, 2)
    print(f"  iterations: {n_iter}")
    Sp_prod = 1.0
    for s in s_final:
        Sp_prod *= sum(s.real**k for k in range(2))
    for i, s in enumerate(s_final):
        x = (2.0, 3.0, 5.0)[i]
        residual = abs(x * Sp_prod * (1 - s.real) - (x - 1))
        print(f"    s_{i} = {s.real:.10f}   eq residual = {residual:.2e}")

    print()
    print("="*70)
    print("CORRECTED VERDICT: hypercube is GENUINELY MULTIVARIATE.")
    print("  - The product polynomial S_p^{cube,d} = prod_i S_p(s_i) couples")
    print("    variables non-trivially: at the symmetric diagonal, the fixed")
    print("    point is a new algebraic number, NOT x^{-1/p}.")
    print("  - For d=2, p=2, x=2: s_* ~ 0.8546 (cubic root, paper-0 gives 0.707).")
    print("  - Asymmetric x_i: genuinely coupled system, NOT d independent")
    print("    paper-0 instances.")
    print("  - This is a TRUE alternative geometry that escapes the simplex's")
    print("    rank-1 collapse and produces new algebraic fixed points.")
    print("="*70)


if __name__ == "__main__":
    main()
