"""
PAPER: 100 (canonical: 90pandrosion_capstone.pdf)
TITLE: Vieta ↔ Newton Unified — One Identity, Ninety Papers
STATUS: synthesis (capstone of Pandrosion legacy series)
DEPENDS: 011, 037, 073, 086

THEORY
======

THE FUNDAMENTAL PANDROSION IDENTITY:
  F_P(z) * P(z) = P'(z),    F_P(z) = sum 1/(z - alpha_j).

This single identity is the load-bearing structure of the entire 90-paper
legacy series. Each operation on F_P unfolds a different theorem:

  Operation               Result
  -----                   ------
  Poles                   Roots
  Residues                Vieta formulas
  Laurent at infinity     Newton power sums
  F_P * P                 P' (definition)
  F_P'                    -T/P^2 (Turán curvature)
  Contour integral        Argument principle
  Real part on S^1        Mahler measure (Jensen)
  Stability of integral   Rouché
  Rightmost pole          Perron-Frobenius
  |F_P|^2 on S^1          Szegő/Toeplitz

VERIFICATION
============

The fundamental identity F_P * P = P' for many random polynomials.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 100 — Capstone: F_P * P = P' (one identity, ninety papers)")
    print("=" * 80)

    print("\n[1] Fundamental identity: F_P(z) * P(z) = P'(z)")
    rng = np.random.default_rng(0)
    max_err = 0.0
    n_tests = 500
    for _ in range(n_tests):
        d = rng.integers(2, 8)
        P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
        z = complex(rng.standard_normal(), rng.standard_normal())
        Pp = np.polyder(P)
        roots = np.roots(P)
        # F_P(z) = sum 1/(z - alpha_j)
        F = sum(1.0 / (z - r) for r in roots)
        lhs = F * np.polyval(P, z)
        rhs = np.polyval(Pp, z)
        err = abs(lhs - rhs)
        if err > max_err: max_err = err
    print(f"  {n_tests} random tests: max |F_P * P - P'| = {max_err:.2e}")

    print("\n[2] Operations on F_P unfold the theory:")
    operations = [
        ("Poles of F_P", "Roots of P (paper 11)"),
        ("Residues of F_P at alpha_k", "Multiplicities (paper 11)"),
        ("Laurent at infinity", "Power sums p_k (paper 73)"),
        ("F_P * P", "P' (the identity itself)"),
        ("F_P'", "-T/P^2 (Turán curvature, papers 30, 86)"),
        ("Contour integral oint F_P", "Root count (argument principle, paper 86)"),
        ("Re F_P on imaginary axis", "Hurwitz stability (paper 75)"),
        ("Re F_P > 0 on |z| = 1", "Schur stability (paper 29)"),
        ("|F_P|^2 on circle", "Szego / Toeplitz (paper 92)"),
        ("Rightmost pole on R+", "Perron-Frobenius (paper 99)"),
    ]
    for op, result in operations:
        print(f"  {op:>30}  ->  {result}")

    print("\n[3] One identity, ninety papers. Pandrosion.")


if __name__ == "__main__":
    main()
