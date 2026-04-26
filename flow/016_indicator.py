"""
PAPER: 016 (canonical: 06pandrosion_indicator.pdf)
TITLE: The Pandrosion Indicator: Topological Tracking of Roots
STATUS: framework (winding-number / argument-principle based)
DEPENDS: 011

THEORY
======

The Pandrosion indicator I_P(gamma) for a closed curve gamma not passing
through any root:
  I_P(gamma) = (1 / 2 pi) * variation of arg(P) around gamma
             = number of roots inside gamma (argument principle).

In Pandrosion language: I_P(gamma) = (1/2pi i) * integral of F_P over gamma
where F_P = P'/P is the Pandrosion field.

VERIFICATION
============

  1. Argument principle: I_P matches root count for various contours.
  2. Pandrosion-field integration recovers winding number.
"""
from __future__ import annotations
import math
import cmath
import numpy as np


def winding_number(P, center, radius, n_pts=512):
    """Estimate (1/2pi i) integral_{|z-center|=R} P'/P dz."""
    thetas = 2 * math.pi * np.arange(n_pts) / n_pts
    zs = center + radius * np.exp(1j * thetas)
    Pp = np.polyder(P)
    F_vals = np.polyval(Pp, zs) / np.polyval(P, zs)
    # Integrate: trapezoidal sum * d_theta with i R e^{i theta}
    dz_dtheta = 1j * radius * np.exp(1j * thetas)
    integrand = F_vals * dz_dtheta
    integral = np.mean(integrand) * 2 * math.pi
    return integral / (2 * math.pi * 1j)


def main():
    print("=" * 80)
    print("PAPER 16 — Pandrosion indicator: argument principle")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Argument principle: count roots inside |z| = R")
    P = np.array([1.0, 0, 0, 0, -1.0])  # z^4 - 1, roots at 1, -1, i, -i
    print("  P = z^4 - 1, roots at ±1, ±i")
    for R in [0.5, 1.5, 2.0]:
        wind = winding_number(P, 0, R)
        # Count actual roots inside
        roots = np.roots(P)
        true_count = sum(1 for r in roots if abs(r) < R)
        print(f"  R = {R}: winding = {wind.real:.4f}, true count = {true_count}")

    print("\n[2] Random polynomial: winding around a small disk")
    for d in [3, 5, 8]:
        P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
        P[0] = 1.0
        roots = np.roots(P)
        # Pick a small disk around the first root
        center = roots[0]
        R = min(abs(roots[0] - r) for r in roots[1:]) / 2
        wind = winding_number(P, center, R)
        print(f"  d={d}: small disk around root, winding = {wind.real:.4f} (expect 1)")

    print("\n[3] Multiplicity detection: (z-1)^3 has triple root at 1")
    P = np.poly([1, 1, 1, 2.0])
    wind = winding_number(P, 1.0, 0.3)
    print(f"  P = (z-1)^3 (z-2): winding around z=1 = {wind.real:.4f} (expect 3)")


if __name__ == "__main__":
    main()
