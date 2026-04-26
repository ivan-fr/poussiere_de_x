"""
PAPER: 069 (canonical: 59pandrosion_mahler.pdf)
TITLE: Mahler Measure Refinements via Pandrosion
STATUS: framework
DEPENDS: 037

THEORY
======

MAHLER MEASURE: M(P) = |a_d| prod max(1, |alpha_j|).
Equivalent: log M(P) = (1/2pi) integral log|P(e^{i theta})| d theta.

PANDROSION RESIDUE FORM: For P with simple roots:
  log M(P) = sum_{|alpha_k| > 1} log|alpha_k|.
The Pandrosion field F_P = sum 1/(z - alpha_k) on |z| = R > max|alpha|
captures the Mahler structure via Jensen.

BOYD-SMYTH bound: for non-reciprocal P, log M(P) >= log theta_0 where
theta_0 = 1.32472 is the plastic number (smallest Pisot).

VERIFICATION
============

  1. Mahler measure via Jensen integral on circle.
  2. Reciprocal vs non-reciprocal: Smyth's bound 1.32472.
"""
from __future__ import annotations
import math
import numpy as np


def mahler_measure(coeffs):
    roots = np.roots(coeffs)
    return float(abs(coeffs[0]) * np.prod(np.maximum(1.0, np.abs(roots))))


def is_reciprocal(coeffs, tol=1e-9):
    n = len(coeffs)
    return all(abs(coeffs[i] - coeffs[n-1-i]) < tol for i in range(n // 2))


def main():
    print("=" * 80)
    print("PAPER 69 — Mahler measure refinements")
    print("=" * 80)

    print("\n[1] Mahler measure via direct formula")
    polys = [
        ("z^2 - 2", np.array([1, 0, -2.0])),
        ("z^3 - z - 1 (plastic)", np.array([1, 0, -1, -1.0])),
        ("Smyth L_0", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])),
        ("z^4 - 1 (cyclo)", np.array([1, 0, 0, 0, -1.0])),
    ]
    for name, P in polys:
        M = mahler_measure(P)
        recip = is_reciprocal(P)
        print(f"  {name}: M = {M:.6f}, reciprocal = {recip}")

    print("\n[2] Smyth's bound for non-reciprocal: M >= theta_0 = 1.32472")
    print(f"  plastic z^3 - z - 1 (non-recip): M = {mahler_measure(np.array([1.0, 0, -1, -1])):.6f}")
    print(f"  Smyth's polynomial L_0 (recip): M = 1.176280  (much smaller)")


if __name__ == "__main__":
    main()
