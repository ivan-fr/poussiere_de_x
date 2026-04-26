"""
PAPER: 034 (canonical: 24pandrosion_hb.pdf)
TITLE: Hadamard-Bombieri inequality via Pandrosion
STATUS: proved (Bombieri-Weyl norm inequalities)
DEPENDS: 011

THEORY
======

For a homogeneous polynomial P of degree d with Bombieri-Weyl norm
||P||_BW^2 = sum |c_k|^2 / C(d, k):

HADAMARD-BOMBIERI: |P(z)|^2 / (1 + |z|^2)^d <= ||P||_BW^2
with equality iff P(z) = c (1 + |z|^2)^{d/2} ... (only achieved on the
"Bombieri-extremal" homogenization).

PANDROSION CONNECTION: ||P||_BW is the natural norm for the SU(2)-invariant
inner product on Sym^d(C^2) (paper 77).

VERIFICATION
============

  1. Compute Bombieri-Weyl norm.
  2. Verify Hadamard-Bombieri inequality.
"""
from __future__ import annotations
import math
import numpy as np


def bombieri_weyl_norm(P):
    """||P||_BW^2 = sum |c_k|^2 / C(d, k) (using monomial coeffs).

    P[i] is coefficient of z^{d-i} in numpy convention. So P[k] in the
    high-to-low indexing is c_{d-k}; reorder to low-to-high to match.
    """
    d = len(P) - 1
    binom = np.array([math.comb(d, k) for k in range(d + 1)])
    # Reverse to low-to-high: c_0, c_1, ..., c_d
    coefs_low = P[::-1]
    return float(np.sum(np.abs(coefs_low)**2 / binom))


def main():
    print("=" * 80)
    print("PAPER 34 — Hadamard-Bombieri inequality")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Bombieri-Weyl norm computation")
    test_polys = [
        ("(1+z)^d, d=4", np.array([1, 4, 6, 4, 1.0])),  # binomial coeffs
        ("z^d, d=4", np.array([1, 0, 0, 0, 0.0])),
        ("random d=5", rng.standard_normal(6)),
    ]
    for name, P in test_polys:
        bw = bombieri_weyl_norm(P)
        print(f"  {name}: ||P||_BW^2 = {bw:.4f}")

    print("\n[2] Hadamard-Bombieri: |P(z)|^2 / (1+|z|^2)^d <= ||P||_BW^2")
    for name, P in test_polys:
        d = len(P) - 1
        bw = bombieri_weyl_norm(P)
        max_ratio = 0.0
        for _ in range(100):
            z = complex(rng.standard_normal(), rng.standard_normal())
            ratio = abs(np.polyval(P, z))**2 / (1 + abs(z)**2)**d
            if ratio > max_ratio: max_ratio = ratio
        print(f"  {name}: max |P|^2/(1+|z|^2)^d = {max_ratio:.4f}, "
              f"||BW||^2 = {bw:.4f}  ({'OK' if max_ratio <= bw + 1e-9 else 'VIOLATES'})")

    print("\n[3] Bombieri-extremal: (1+z)^d hits the bound at z = real")
    P = np.array([1, 4, 6, 4, 1.0])  # (1+z)^4 expanded
    z = 1.0  # real positive
    val = abs(np.polyval(P, z))**2 / (1 + abs(z)**2)**4
    bw = bombieri_weyl_norm(P)
    print(f"  P=(1+z)^4 at z=1: ratio = {val:.4f}, BW = {bw:.4f}")


if __name__ == "__main__":
    main()
