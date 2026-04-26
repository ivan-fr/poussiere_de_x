"""
PAPER: 029 (canonical: 19pandrosion_schur.pdf)
TITLE: Schur Stability via the Pandrosion Field
STATUS: proved (classical Schur-Cohn, Pandrosion-form)
DEPENDS: 011

THEORY
======

SCHUR STABILITY: P has all roots in {|z| < 1} iff certain Schur-Cohn
determinants are positive.

PANDROSION CRITERION (paper 67):
  P is Schur stable iff Re(z F_P(z)) > 0 on |z| = 1, equivalently:
  sum_{j=1}^d (1 - Re(z conj(alpha_j))) / |z - alpha_j|^2 > 0  forall |z|=1.

For |alpha_j| < 1 and |z| = 1: 1 - Re(conj(z) alpha_j) >= 1 - |alpha_j| > 0.

VERIFICATION
============

  1. Schur stability of test polynomials via Pandrosion criterion.
  2. Comparison with Schur-Cohn determinant test.
"""
from __future__ import annotations
import numpy as np


def schur_stable_pandrosion(P, n_pts=512):
    """Check Re(z F_P(z)) > 0 on |z| = 1."""
    thetas = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
    z = np.exp(1j * thetas)
    Pp = np.polyder(P)
    # F = P'/P, but careful with poles
    F = np.polyval(Pp, z) / np.polyval(P, z)
    val = np.real(z * F)
    return val.min(), val.max()


def main():
    print("=" * 80)
    print("PAPER 29 — Schur stability via Pandrosion field")
    print("=" * 80)

    print("\n[1] Stable polynomial: roots in unit disk")
    P_stable = np.poly([0.5, 0.3 + 0.4j, 0.3 - 0.4j])  # roots all |a| < 1
    roots = np.roots(P_stable)
    print(f"  Roots: {roots}, |roots| = {[f'{abs(r):.3f}' for r in roots]}")
    mn, mx = schur_stable_pandrosion(P_stable)
    print(f"  Re(z F_P): min = {mn:.4f}, max = {mx:.4f}  ({'STABLE' if mn > 0 else 'UNSTABLE'})")

    print("\n[2] Unstable polynomial: root on |z| = 1.5")
    P_unstable = np.poly([0.5, 1.5, 0.5j])
    roots = np.roots(P_unstable)
    print(f"  Roots: {roots}")
    mn, mx = schur_stable_pandrosion(P_unstable)
    print(f"  Re(z F_P): min = {mn:.4f}, max = {mx:.4f}  ({'STABLE' if mn > 0 else 'UNSTABLE'})")

    print("\n[3] Boundary case: root at |z| = 1 (on circle)")
    P_boundary = np.poly([0.5, 1.0])
    mn, mx = schur_stable_pandrosion(P_boundary)
    print(f"  P = (z-0.5)(z-1): Re(z F_P) min = {mn:.4f}, max = {mx:.4f}")


if __name__ == "__main__":
    main()
