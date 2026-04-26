"""
PAPER: 056 (canonical: 46pandrosion_hp.pdf)
TITLE: Hardy-Pandrosion (HP) inequalities
STATUS: framework

THEORY
======

Hardy-Pandrosion: bounds on power-sum / coefficient norms of polynomials
via the Pandrosion field. For P(z) = sum c_k z^k:
  ||P||_2^2 = sum |c_k|^2  (Parseval on |z|=1)
  ||P||_inf <= sqrt(d+1) ||P||_2  (Hardy)

VERIFICATION
============

  1. Parseval identity for ||P||_2.
  2. Hardy bound ||P||_inf <= sqrt(d+1) ||P||_2.
"""
from __future__ import annotations
import math
import numpy as np


def L2_norm_circle(P, n_pts=512):
    z = np.exp(2j * np.pi * np.arange(n_pts) / n_pts)
    return float(math.sqrt(np.mean(np.abs(np.polyval(P, z))**2)))


def Linf_norm_circle(P, n_pts=512):
    z = np.exp(2j * np.pi * np.arange(n_pts) / n_pts)
    return float(np.max(np.abs(np.polyval(P, z))))


def main():
    print("=" * 80)
    print("PAPER 56 — Hardy-Pandrosion inequalities")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Parseval: ||P||_2 = sqrt(sum |c_k|^2)")
    for d in [3, 5, 8]:
        P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
        L2 = L2_norm_circle(P)
        coef_norm = math.sqrt(sum(abs(c)**2 for c in P))
        print(f"  d={d}: L2 = {L2:.4f}, sqrt(sum|c|^2) = {coef_norm:.4f}")

    print("\n[2] Hardy: ||P||_inf <= sqrt(d+1) ||P||_2")
    for d in [3, 5, 10]:
        P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
        L2 = L2_norm_circle(P)
        Linf = Linf_norm_circle(P)
        bound = math.sqrt(d + 1) * L2
        print(f"  d={d}: L_inf = {Linf:.4f}, sqrt(d+1) L_2 = {bound:.4f}, "
              f"ratio = {Linf/bound:.4f}")


if __name__ == "__main__":
    main()
