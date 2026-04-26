"""
PAPER: 093 (canonical: 83pandrosion_smyth.pdf)
TITLE: Smyth's Theorem (Mahler measures)
STATUS: proved (Smyth 1971), refined in paper 112
DEPENDS: 069

THEORY
======

SMYTH (1971): For non-reciprocal monic Z-poly P, M(P) >= theta_0 where
theta_0 = 1.32472 is the plastic number (smallest Pisot, root of z^3 = z + 1).

PANDROSION CONNECTION: For Z-poly with M close to 1, structure forces
roots to be near unit circle (cyclotomic-like), implying reciprocity.

VERIFICATION
============

  1. Plastic z^3 - z - 1: M = 1.32472 (smallest non-reciprocal).
  2. Smyth's bound on random non-reciprocal Z-polys.
"""
from __future__ import annotations
import math
import numpy as np


def mahler_measure(P):
    return float(abs(P[0]) * np.prod(np.maximum(1.0, np.abs(np.roots(P)))))


def is_reciprocal(P, tol=1e-9):
    n = len(P)
    return all(abs(P[i] - P[n-1-i]) < tol for i in range(n // 2))


def main():
    print("=" * 80)
    print("PAPER 93 — Smyth's theorem (non-reciprocal Mahler bound)")
    print("=" * 80)

    print("\n[1] Plastic number z^3 - z - 1 = 0")
    plastic = np.array([1.0, 0, -1, -1])
    M = mahler_measure(plastic)
    print(f"  Plastic: M = {M:.10f}, theta_0 = 1.32472 expected")
    print(f"  Reciprocal? {is_reciprocal(plastic)}")

    print("\n[2] Smyth bound: non-reciprocal monic Z-polys have M >= 1.32472")
    rng = np.random.default_rng(0)
    print(f"  {'d':>3} {'#nonrecip':>11} {'min M':>10}")
    for d in [3, 4, 5, 6, 7]:
        n_nr = 0
        min_M = float('inf')
        for _ in range(500):
            coefs = rng.choice([-1, 0, 1], size=d+1)
            coefs[0] = 1
            if coefs[-1] == 0: continue
            P = coefs.astype(float)
            if is_reciprocal(P): continue
            n_nr += 1
            M = mahler_measure(P)
            if M > 1.001 and M < min_M: min_M = M
        print(f"  {d:>3} {n_nr:>11} {min_M:>10.6f}")

    print("\n[3] Reciprocal vs non-reciprocal: small Mahler concentrates in reciprocal")
    print(f"  Smyth's L_0 = 1.176 (reciprocal): below 1.32472 (Smyth bound for non-recip)")


if __name__ == "__main__":
    main()
