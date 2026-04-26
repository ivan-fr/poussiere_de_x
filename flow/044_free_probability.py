"""
PAPER: 044 (canonical: 34pandrosion_free.pdf)
TITLE: Free Probability and the Pandrosion Cauchy Transform
STATUS: framework

THEORY
======

For random matrix H, the Cauchy transform of its spectrum is
  G_mu(z) = E[Tr (zI - H)^{-1}] / d.

PANDROSION CONNECTION: For characteristic poly P(z) = det(zI - H):
  Tr (zI - H)^{-1} = P'(z) / P(z) = F_P(z).
So G_mu(z) = F_P(z) / d. The Cauchy transform IS the Pandrosion field
divided by d.

In free probability, R-transform R(z) = G_mu^{-1}(z) - 1/z is the additive
counterpart, and Voiculescu's S-transform handles multiplicative free
convolution.

VERIFICATION
============

  1. Cauchy transform = F_P / d for characteristic polys.
  2. Wigner semicircle as Cauchy transform limit.
"""
from __future__ import annotations
import math
import numpy as np


def cauchy_transform(P, z):
    """G(z) = Tr(zI - H)^{-1} / d = (P'/P)(z) / d."""
    Pp = np.polyder(P)
    d = len(P) - 1
    return np.polyval(Pp, z) / np.polyval(P, z) / d


def main():
    print("=" * 80)
    print("PAPER 44 — Free probability via Pandrosion Cauchy transform")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Cauchy transform G(z) = F_P(z)/d for char poly")
    # GUE: random Hermitian matrix
    d = 50
    A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    H = (A + A.conj().T) / 2 / np.sqrt(d)
    eigs = np.linalg.eigvalsh(H)
    P = np.poly(eigs)

    z = 3.0  # away from spectrum
    G_via_inverse = float(np.real(np.trace(np.linalg.inv(z * np.eye(d) - H)) / d))
    G_via_pandrosion = float(np.real(cauchy_transform(P, z)))
    print(f"  d={d}, z={z}:")
    print(f"  G via Tr(zI-H)^-1/d = {G_via_inverse:.6f}")
    print(f"  G via F_P/d         = {G_via_pandrosion:.6f}")
    print(f"  diff = {abs(G_via_inverse - G_via_pandrosion):.2e}")

    print("\n[2] Wigner semicircle limit")
    # Cauchy transform of semicircle on [-2, 2]: G(z) = (z - sqrt(z^2 - 4))/2
    print(f"  z       G_empirical (d=50)     G_semicircle")
    for z_test in [3.0, 2.5, 5.0]:
        G_emp = float(np.real(cauchy_transform(P, z_test)))
        G_sc = (z_test - math.sqrt(z_test**2 - 4)) / 2
        print(f"  {z_test}: {G_emp:.6f}            {G_sc:.6f}")


if __name__ == "__main__":
    main()
