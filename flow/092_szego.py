"""
PAPER: 092 (canonical: 82pandrosion_szego.pdf)
TITLE: Szegő's Theorem and Pandrosion-Mahler
STATUS: proved (Szegő 1915)
DEPENDS: 037, 069

THEORY
======

SZEGŐ'S THEOREM: For Toeplitz determinants D_n associated to a positive
function w on |z| = 1:
  D_n^{1/n} -> G(w) := exp((1/2pi) integral log w d theta).

For w = |P|^2 with P a polynomial, G(w) = M(P)^2 (Mahler measure squared).
This is the Pandrosion-Mahler manifestation of Szegő.

VERIFICATION
============

  1. Toeplitz determinant convergence.
  2. Limit = M(P)^2 for w = |P|^2.
"""
from __future__ import annotations
import math
import numpy as np


def toeplitz_from_function(w_func, n, n_fourier=512):
    """Compute Toeplitz determinant D_n from Fourier coefficients of w."""
    # Fourier coefficients
    thetas = 2 * np.pi * np.arange(n_fourier) / n_fourier
    w_vals = w_func(thetas)
    c = np.fft.fft(w_vals) / n_fourier
    # T_{ij} = c_{i-j}
    T = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            k = (i - j) % n_fourier
            T[i, j] = c[k]
    return T


def main():
    print("=" * 80)
    print("PAPER 92 — Szegő's theorem and Pandrosion-Mahler")
    print("=" * 80)

    P = np.array([1.0, 0, -2])  # z^2 - 2
    roots = np.roots(P)
    M = float(np.prod(np.maximum(1.0, np.abs(roots))))
    print(f"\n[1] P = z^2 - 2, M(P) = {M:.4f}")

    # Toeplitz D_n for w = |P|^2
    def w_func(theta):
        z = np.exp(1j * np.array(theta))
        return np.abs(np.polyval(P, z))**2

    print(f"\n[2] Toeplitz D_n^{{1/n}} -> M(P)^2 = {M**2:.4f}")
    for n in [3, 5, 8, 12]:
        T = toeplitz_from_function(w_func, n)
        sign, log_det = np.linalg.slogdet(T)
        D_root = math.exp(log_det.real / n)
        print(f"  n = {n}: D_n^(1/n) = {D_root:.4f}")


if __name__ == "__main__":
    main()
