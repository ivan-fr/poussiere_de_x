"""
PAPER: 021 (canonical: 11pandrosion_smale_mv.pdf)
TITLE: Pandrosion-Smale Multivariate: Initial Treatment
STATUS: framework (refined in paper 010)
DEPENDS: 010

THEORY
======

Multivariate Pandrosion operator on F: C^n -> C^n:
  P_{F, z_0}(z) := z_0 - Q_F(z_0, z)^{-1} F(z_0)
with Schmidt slope matrix Q_F satisfying telescoping
  F(z) - F(z_0) = Q_F(z_0, z) (z - z_0).

When z_0 = z, Q_F = DF (Jacobian) and the iteration is multivariate Newton.

VERIFICATION
============

  1. Schmidt slope matrix Q_F via coordinate-by-coordinate finite differences.
  2. Telescoping identity F(z) - F(z_0) = Q_F (z - z_0).
  3. Newton recovery: Q_F(z, z) = DF(z).
"""
from __future__ import annotations
import numpy as np


def Q_F(F, z0, z):
    n = len(z)
    z0 = np.array(z0, dtype=complex)
    z = np.array(z, dtype=complex)
    Q = np.zeros((n, n), dtype=complex)
    z_curr = z0.copy()
    F_curr = np.array(F(z_curr), dtype=complex)
    for j in range(n):
        z_next = z_curr.copy()
        z_next[j] = z[j]
        F_next = np.array(F(z_next), dtype=complex)
        if abs(z[j] - z0[j]) < 1e-15:
            h = 1e-7
            z_h = z_curr.copy()
            z_h[j] += h
            F_h = np.array(F(z_h), dtype=complex)
            Q[:, j] = (F_h - F_curr) / h
        else:
            Q[:, j] = (F_next - F_curr) / (z[j] - z0[j])
        z_curr = z_next
        F_curr = F_next
    return Q


def main():
    print("=" * 80)
    print("PAPER 21 — Multivariate Pandrosion: Schmidt slope matrix")
    print("=" * 80)

    def F(z):
        return [z[0]**2 + z[1]**2 - 1, z[0] * z[1] - 0.5]

    rng = np.random.default_rng(0)

    print("\n[1] Telescoping identity verification")
    max_err = 0.0
    for _ in range(20):
        z0 = rng.standard_normal(2) + 1j * rng.standard_normal(2)
        z = rng.standard_normal(2) + 1j * rng.standard_normal(2)
        Q = Q_F(F, z0, z)
        F_z = np.array(F(z))
        F_z0 = np.array(F(z0))
        err = np.linalg.norm(F_z - F_z0 - Q @ (z - z0))
        if err > max_err: max_err = err
    print(f"  max ||F(z) - F(z_0) - Q (z - z_0)|| = {max_err:.2e}")

    print("\n[2] Newton recovery: Q(z, z) = DF(z)")
    z = np.array([0.7 + 0.5j, 0.3 - 0.2j])
    Q_diag = Q_F(F, z, z)
    # Analytic Jacobian: [[2x, 2y], [y, x]]
    DF = np.array([[2*z[0], 2*z[1]], [z[1], z[0]]])
    err = np.linalg.norm(Q_diag - DF)
    print(f"  ||Q(z,z) - DF(z)|| = {err:.2e}")

    print("\n[3] Pandrosion-Newton iteration")
    z = np.array([0.95, 0.55])
    for k in range(8):
        Q = Q_F(F, z, z)
        F_z = np.array(F(z))
        z = z - np.linalg.solve(Q, F_z)
        print(f"  iter {k+1}: z = {z}, |F| = {np.linalg.norm(F(z)):.2e}")


if __name__ == "__main__":
    main()
