"""
PAPER: 049 (canonical: 39pandrosion_quantum.pdf)
TITLE: Pandrosion in Quantum Mechanics
STATUS: framework

THEORY
======

For a Hamiltonian H with characteristic polynomial P:
  Resolvent R(z) = (zI - H)^{-1}, Tr R(z) = F_P(z) (paper 44).
  Spectral measure: roots of P = energy eigenvalues.

GREEN'S FUNCTION: G(z) = (z - E_n)^{-1} integrated over n-th eigenstate.
The Pandrosion identity sum 1/(z - E_n) = F_P(z) relates spectrum to
Pandrosion field.

VERIFICATION
============

  1. Resolvent = Pandrosion field for finite Hamiltonian.
  2. Spectral resolution.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 49 — Pandrosion in quantum mechanics: resolvent")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Resolvent Tr(zI - H)^{-1} = F_P(z) for char poly P")
    d = 5
    A = rng.standard_normal((d, d))
    H = (A + A.T) / 2  # Hermitian (real symmetric here)
    eigs = np.linalg.eigvalsh(H)
    P = np.poly(eigs)
    Pp = np.polyder(P)
    print(f"  H eigenvalues: {eigs}")
    for z_test in [10.0, 5.0, 1.0]:
        # Avoid eigenvalues
        if any(abs(z_test - e) < 0.1 for e in eigs): z_test += 0.5
        Tr_R = float(np.real(np.trace(np.linalg.inv(z_test * np.eye(d) - H))))
        F_P = np.polyval(Pp, z_test) / np.polyval(P, z_test)
        print(f"  z = {z_test}: Tr R = {Tr_R:.6f}, F_P = {F_P:.6f}")

    print("\n[2] Spectral measure: poles at eigenvalues")
    P = np.array([1.0, -3, 2])  # eigenvalues 1, 2
    print(f"  P = z^2 - 3z + 2 (eigenvalues 1, 2)")
    Pp = np.polyder(P)
    for z in [0.5, 1.5, 2.5]:
        F = np.polyval(Pp, z) / np.polyval(P, z)
        print(f"  F_P({z}) = {F:.4f} = 1/(z - 1) + 1/(z - 2) = {1/(z-1) + 1/(z-2):.4f}")


if __name__ == "__main__":
    main()
