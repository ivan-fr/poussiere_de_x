"""
PAPER: 038 (canonical: 28pandrosion_descartes.pdf)
TITLE: Descartes' Rule of Signs via Pandrosion
STATUS: proved (Descartes 1637, sign-pattern reformulation)
DEPENDS: 011

THEORY
======

DESCARTES' RULE: The number of positive real roots of P (counted with
multiplicity) is at most the number of sign changes in the coefficient
sequence (a_0, a_1, ..., a_d), and differs from it by an even number.

PANDROSION FORM: For P(z) = sum c_k z^k, F_P(x) = P'(x)/P(x) on x > 0:
  F_P(x) = (sum k c_k x^{k-1}) / (sum c_k x^k).
The number of zeros of P on (0, infty) is bounded by the sign changes in
{c_0, c_1, ..., c_d}.

VERIFICATION
============

  1. Descartes rule on test polynomials.
  2. Tightness: equality for polynomials with all positive real roots.
"""
from __future__ import annotations
import numpy as np


def sign_changes(coefs):
    """Count sign changes in non-zero coefficients."""
    nonzero = [c for c in coefs if abs(c) > 1e-12]
    sgns = [1 if c > 0 else -1 for c in nonzero]
    return sum(1 for i in range(len(sgns)-1) if sgns[i] * sgns[i+1] < 0)


def main():
    print("=" * 80)
    print("PAPER 38 — Descartes' rule of signs")
    print("=" * 80)

    print("\n[1] Number of positive real roots vs sign changes")
    print(f"  {'P':>40} {'#pos roots':>12} {'sign changes':>14}")
    polys = [
        ("z^2 - 3z + 2 = (z-1)(z-2)", [1, -3, 2]),
        ("z^3 - 6z^2 + 11z - 6 = (z-1)(z-2)(z-3)", [1, -6, 11, -6]),
        ("z^3 + 1", [1, 0, 0, 1]),
        ("z^4 - 1", [1, 0, 0, 0, -1]),
        ("z^5 - z^4 + z^3 - z^2 + z - 1", [1, -1, 1, -1, 1, -1]),
    ]
    for name, P in polys:
        # Reverse coefs to low-to-high for Descartes (P = sum c_k x^k)
        c_low = list(reversed(P))
        sc = sign_changes(c_low)
        roots = np.roots(P)
        n_pos = sum(1 for r in roots if abs(r.imag) < 1e-9 and r.real > 1e-9)
        print(f"  {name:>40} {n_pos:>12} {sc:>14}")

    print("\n[2] Equality case: all positive real roots")
    P = np.poly([1, 2, 3, 4.0])  # all positive roots
    coefs = list(reversed(P))
    sc = sign_changes(coefs)
    n_pos = sum(1 for r in np.roots(P) if r.real > 0 and abs(r.imag) < 1e-9)
    print(f"  P = (z-1)(z-2)(z-3)(z-4): #pos roots = {n_pos}, sign changes = {sc}")


if __name__ == "__main__":
    main()
