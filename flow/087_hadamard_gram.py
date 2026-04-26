"""
PAPER: 087 (canonical: 77pandrosion_hadamard.pdf)
TITLE: Hadamard's Inequality and the Pandrosion Gram Matrix
STATUS: proved (Pandrosion-Hadamard identity: det G = |Disc|^2)
DEPENDS: 011, 032

THEORY
======

PANDROSION GRAM MATRIX on the unit circle:
  G_{ij} = (1/2pi) integral Q_i(e^{i theta}) conj(Q_j(e^{i theta})) d theta

KEY IDENTITY (Pandrosion-Hadamard):
  det G = |Disc(P)|^2 = prod_{i<j} |alpha_i - alpha_j|^2.

ROOTS OF UNITY: Q_k of P(z) = z^n - 1 are mutually orthogonal.

VERIFICATION
============

  1. det G = |Disc|^2 to machine precision.
  2. Roots of unity: Q_k orthogonal.
"""
from __future__ import annotations
import math
import numpy as np


def Q_poly(P, k):
    """Q_k(z) = P(z)/(z - alpha_k)."""
    roots = np.roots(P)
    other = [r for i, r in enumerate(roots) if i != k]
    return np.poly(other)


def gram_circle(P, n_pts=512):
    d = len(P) - 1
    thetas = 2 * np.pi * np.arange(n_pts) / n_pts
    z = np.exp(1j * thetas)
    Q_vals = np.zeros((d, n_pts), dtype=complex)
    for k in range(d):
        Qk = Q_poly(P, k)
        Q_vals[k] = np.polyval(Qk, z)
    G = (Q_vals @ Q_vals.conj().T) / n_pts
    return G


def main():
    print("=" * 80)
    print("PAPER 87 — Pandrosion-Hadamard: det G = |Disc|^2")
    print("=" * 80)

    print("\n[1] det G = |Disc|^2 (machine precision)")
    test_polys = [
        ("z^2 - 1", np.array([1.0, 0, -1])),
        ("z^3 - z - 1", np.array([1.0, 0, -1, -1])),
        ("z^4 + 1", np.array([1.0, 0, 0, 0, 1])),
        ("Smyth L_0", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])),
    ]
    for name, P in test_polys:
        G = gram_circle(P)
        sign, log_det = np.linalg.slogdet(G)
        # Discriminant
        roots = np.roots(P)
        d = len(roots)
        log_disc = sum(2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
                      for i in range(d) for j in range(i+1, d))
        print(f"  {name}: log det G = {log_det.real:.4f}, log |Disc|^2 = {log_disc:.4f}, "
              f"diff = {abs(log_det.real - log_disc):.2e}")

    print("\n[2] Roots of unity z^n - 1: Q_k orthogonal")
    for n in [3, 4, 5]:
        P = np.array([1.0] + [0]*(n-1) + [-1])
        G = gram_circle(P)
        # Off-diagonal magnitudes
        max_off = max(abs(G[i, j]) for i in range(n) for j in range(n) if i != j)
        print(f"  z^{n} - 1: max off-diagonal |G_ij| = {max_off:.4e} (should ~ 0)")


if __name__ == "__main__":
    main()
