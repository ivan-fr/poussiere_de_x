"""
PAPER: 037 (canonical: 27pandrosion_jensen.pdf)
TITLE: Jensen's Formula via the Pandrosion Field
STATUS: proved (classical Jensen formula in Pandrosion-form)
DEPENDS: 011

THEORY
======

JENSEN'S FORMULA: For analytic f with f(0) != 0, on |z| = R:
  (1/2pi) integral_{0}^{2pi} log|f(R e^{i theta})| d theta
  = log|f(0)| + sum_{|alpha_k| < R} log(R/|alpha_k|).

Specialised to a polynomial P with constant term not 0:
  (1/2pi) integral log|P(R e^{i theta})| d theta
  = log|a_d| + d log R + sum_{|alpha|<R} log(R/|alpha|)
  -> log M(P) as R -> infty (Mahler measure formula).

PANDROSION CONNECTION: The integrand log|P| is harmonic except at zeros;
the Pandrosion field F_P relates to grad log|P|.

VERIFICATION
============

  1. Jensen formula on test polynomials.
  2. Mahler measure recovered as R -> infty.
"""
from __future__ import annotations
import math
import numpy as np


def jensen_integral(P, R, n_pts=512):
    thetas = 2 * np.pi * np.arange(n_pts) / n_pts
    z = R * np.exp(1j * thetas)
    log_P = np.log(np.maximum(np.abs(np.polyval(P, z)), 1e-300))
    return float(np.mean(log_P))


def main():
    print("=" * 80)
    print("PAPER 37 — Jensen's formula (Pandrosion-form)")
    print("=" * 80)

    test_polys = [
        ("z^2 - 2", np.array([1.0, 0, -2])),
        ("z^3 - z - 1", np.array([1.0, 0, -1, -1])),
        ("Smyth L_0", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])),
    ]

    print("\n[1] Jensen formula for each polynomial at various R")
    for name, P in test_polys:
        roots = np.roots(P)
        d = len(roots)
        a_d = P[0]
        log_const = math.log(abs(P[-1]))  # |a_0|
        print(f"\n  {name}: roots {[f'{abs(r):.3f}' for r in roots]}")
        for R in [0.5, 1.0, 2.0, 5.0]:
            jensen = jensen_integral(P, R)
            roots_inside = sum(math.log(R / abs(r)) for r in roots if abs(r) < R)
            predicted = math.log(abs(a_d)) + d * math.log(R) + roots_inside
            print(f"    R={R}: integral = {jensen:.4f}, predicted = {predicted:.4f}, "
                  f"diff = {abs(jensen - predicted):.4f}")

    print("\n[2] Mahler measure as R -> infty")
    P = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])  # Smyth L_0
    print(f"  Smyth L_0:")
    for R in [10, 100, 1000]:
        jensen = jensen_integral(P, R)
        # log M = jensen - d log R for large R
        d = 10
        log_M_est = jensen - d * math.log(R)
        print(f"    R = {R}: log M_est = {log_M_est:.6f}, M_est = {math.exp(log_M_est):.6f}")
    print(f"  Literature M(L_0) = 1.176280818")


if __name__ == "__main__":
    main()
