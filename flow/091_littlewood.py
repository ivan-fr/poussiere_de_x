"""
PAPER: 091 (canonical: 81pandrosion_littlewood.pdf)
TITLE: Littlewood's Polynomials
STATUS: framework (open Erdős-Littlewood)
DEPENDS: 011

THEORY
======

LITTLEWOOD POLYNOMIAL: P with coefficients in {+1, -1}.
Erdős-Littlewood problem: do flat Littlewood polynomials exist?
  c_1 sqrt(d+1) <= |P_n(e^{i theta})| <= c_2 sqrt(d+1)  uniformly in theta.

Beck (1991): max/min ratio >= 1.196.
Rudin-Shapiro: max/L^2 -> sqrt 2, but min = 0.

VERIFICATION
============

  1. Rudin-Shapiro polynomial recursion.
  2. max/L^2 ratio approaches sqrt 2.
"""
from __future__ import annotations
import math
import numpy as np


def rudin_shapiro(k):
    """Build Rudin-Shapiro polynomials via P_{n+1} = P_n + z^{2^n} Q_n,
    Q_{n+1} = P_n - z^{2^n} Q_n. Returns P_k."""
    P = np.array([1, 1])
    Q = np.array([1, -1])
    for _ in range(k):
        Pn = np.concatenate([P, Q])
        Qn = np.concatenate([P, -Q])
        P, Q = Pn, Qn
    return P


def main():
    print("=" * 80)
    print("PAPER 91 — Littlewood polynomials and Rudin-Shapiro")
    print("=" * 80)

    print("\n[1] Rudin-Shapiro: max/L^2 -> sqrt(2)")
    print(f"  {'k':>3} {'deg':>5} {'max':>10} {'L^2':>10} {'max/L^2':>10}")
    for k in range(1, 8):
        P = rudin_shapiro(k)
        n_pts = 2048
        z = np.exp(2j * np.pi * np.arange(n_pts) / n_pts)
        vals = np.abs(np.polyval(P[::-1], z))
        max_v = float(vals.max())
        L2 = math.sqrt(np.mean(vals**2))
        ratio = max_v / L2
        print(f"  {k:>3} {len(P)-1:>5} {max_v:>10.4f} {L2:>10.4f} {ratio:>10.4f}")
    print(f"  sqrt(2) = {math.sqrt(2):.4f}")

    print("\n[2] Beck 1991: max/min >= 1.196 for any Littlewood P")
    P = rudin_shapiro(3)  # k=3, degree 15
    z = np.exp(2j * np.pi * np.arange(2048) / 2048)
    vals = np.abs(np.polyval(P[::-1], z))
    print(f"  Rudin-Shapiro k=3 (deg 15): max/min = {vals.max() / max(vals.min(), 1e-9):.4f}")


if __name__ == "__main__":
    main()
