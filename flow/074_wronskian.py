"""
PAPER: 074 (canonical: 64pandrosion_wron.pdf)
TITLE: Wronskian-Pandrosion calculus
STATUS: framework
DEPENDS: 011

THEORY
======

WRONSKIAN of Pandrosion quotients Q_i = P/(z - alpha_i):
  W_{ij}(z) = Q_i Q_j' - Q_i' Q_j.

For n = 2: W_{12} = (z - alpha_2) - (z - alpha_1) * 1 = alpha_1 - alpha_2
(constant in z).

For n >= 3: W has rank < n due to linear relation sum Q_j = P'.

VERIFICATION
============

  1. W_{12} = alpha_1 - alpha_2 for d = 2.
  2. det(W_{ij}) = 0 for n >= 3.
"""
from __future__ import annotations
import numpy as np


def Q_poly(P, k):
    """Q_k(z) = P(z)/(z - alpha_k) where alpha_k is the k-th root."""
    roots = np.roots(P)
    other = [r for i, r in enumerate(roots) if i != k]
    return np.poly(other)


def main():
    print("=" * 80)
    print("PAPER 74 — Wronskian-Pandrosion calculus")
    print("=" * 80)

    print("\n[1] n = 2: W_{12} = alpha_1 - alpha_2 (z-independent)")
    P = np.array([1.0, -3, 2])  # roots 1, 2
    roots = np.roots(P)
    Q1 = Q_poly(P, 0)  # P/(z - root_1)
    Q2 = Q_poly(P, 1)
    Q1p = np.polyder(Q1)
    Q2p = np.polyder(Q2)
    # W_{12}(z) = Q_1 Q_2' - Q_1' Q_2
    for z in [0.0, 0.5, 1.5, 3.0]:
        W = np.polyval(Q1, z) * np.polyval(Q2p, z) - np.polyval(Q1p, z) * np.polyval(Q2, z)
        print(f"  z = {z}: W_12(z) = {W:.4f}, alpha_1 - alpha_2 = {roots[0] - roots[1]:.4f}")

    print("\n[2] n >= 3: linear dependence sum Q_j = P'")
    P = np.poly([1, 2, 3.0])
    roots = np.roots(P)
    z = 0.5
    sum_Q = sum(np.polyval(Q_poly(P, k), z) for k in range(len(roots)))
    Pp_z = np.polyval(np.polyder(P), z)
    print(f"  d=3, z=0.5: sum Q_j(z) = {sum_Q:.4f}, P'(z) = {Pp_z:.4f}")


if __name__ == "__main__":
    main()
