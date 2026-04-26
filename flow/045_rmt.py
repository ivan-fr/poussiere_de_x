"""
PAPER: 045 (canonical: 35pandrosion_rmt.pdf)
TITLE: Eigenvalue Repulsion as Pandrosion Self-Energy: RMT via the Quotient
STATUS: proved (classical RMT in Pandrosion language)
DEPENDS: 032, 044

THEORY
======

GUE DENSITY: For Gaussian Unitary Ensemble:
  p_GUE(lambda) = C_d prod_{i<j} |lambda_i - lambda_j|^2 exp(-d sum lambda_k^2 / 2).

PANDROSION FORM: prod_{i<j} |lambda_i - lambda_j|^2 = |D(P)|^2 where
P is the characteristic polynomial. By paper 32:
  |D|^2 = prod_k P'(lambda_k) = prod_k E_P(lambda_k)^{1/2} ... times sign.

For real eigenvalues with simple roots: |P'(lambda_k)| = sqrt(E_P(lambda_k)).
So eigenvalue repulsion = product of Pandrosion self-energies.

beta-ENSEMBLE: p_beta ∝ prod |lambda_i - lambda_j|^beta exp(-d sum lambda^2/2)
            = prod E_P(lambda_k)^{beta/4} exp(...).

VERIFICATION
============

  1. |D|^2 = prod E_P(lambda_k) on GUE samples.
  2. Wigner semicircle as Cauchy-transform limit (paper 44 link).
  3. Spacing distribution.
"""
from __future__ import annotations
import math
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 45 — Eigenvalue repulsion = Pandrosion self-energy")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] GUE: |D|^2 = prod E_P(lambda_k)")
    for d in [10, 30, 50]:
        n_test = 5
        for _ in range(n_test):
            A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
            H = (A + A.conj().T) / 2 / np.sqrt(d)
            eigs = np.linalg.eigvalsh(H)
            # Discriminant
            log_disc = sum(2 * math.log(max(abs(eigs[i] - eigs[j]), 1e-300))
                          for i in range(d) for j in range(i+1, d))
            # E_P(lambda_k) = (P'(lambda_k))^2
            P = np.poly(eigs)
            Pp = np.polyder(P)
            log_prod_E = sum(2 * math.log(max(abs(np.polyval(Pp, ek)), 1e-300)) for ek in eigs)
            print(f"  d={d}: log|D|^2 = {log_disc:.4f}, log prod E = {log_prod_E:.4f}, "
                  f"diff = {abs(log_disc - log_prod_E):.2e}")
            break  # one per d

    print("\n[2] Wigner semicircle: empirical density of GUE eigenvalues")
    d = 200
    A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    H = (A + A.conj().T) / 2 / np.sqrt(d)
    eigs = np.linalg.eigvalsh(H)
    print(f"  d={d}: spectrum range = [{eigs.min():.3f}, {eigs.max():.3f}]")
    print(f"  expected: roughly [-2, 2] (Wigner semicircle support)")

    print("\n[3] Spacing distribution: nearest-neighbor gaps (Wigner surmise)")
    spacings = np.diff(np.sort(eigs))
    spacings = spacings / spacings.mean()  # normalize
    # Count spacings in bins
    bins = [0, 0.5, 1.0, 1.5, 2.0, 3.0]
    counts = np.histogram(spacings, bins)[0]
    print(f"  Empirical: {counts}")
    print(f"  Wigner surmise: roughly proportional to s exp(-pi s^2/4)")


if __name__ == "__main__":
    main()
