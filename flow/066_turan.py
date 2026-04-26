"""
PAPER: 066 (canonical: 56pandrosion_turan.pdf)
TITLE: Turán's Inequality as Pandrosion Squares
STATUS: proved (Turán 1939, Pandrosion-form proof)
DEPENDS: 030

THEORY
======

TURÁN (1939): For real-rooted polynomial P with simple roots:
  (P'(x))^2 - P(x) P''(x) > 0   for all real x.

PANDROSION FORM: This is the Pandrosion ENERGY inequality
  E_P(x) = sum Q(alpha_k, x)^2 > 0   on R (paper 30, 33).

Equivalently: T(z) = E_P/P^2 = (sum_k 1/(z-alpha_k)^2) > 0 (real line).

TURÁN-CRITERION (converse): if (P')^2 - P P'' >= 0 on R, then P has only
real roots (real-rootedness criterion).

VERIFICATION
============

  1. Turán inequality on real-rooted polynomials.
  2. Counterexample for non-real-rooted P.
  3. Tightness near double root.
"""
from __future__ import annotations
import numpy as np


def turan_form(P, x):
    Pp = np.polyder(P)
    Ppp = np.polyder(Pp)
    return np.polyval(Pp, x)**2 - np.polyval(P, x) * np.polyval(Ppp, x)


def main():
    print("=" * 80)
    print("PAPER 66 — Turán's inequality as Pandrosion squares")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Real-rooted P: T(x) > 0 on R")
    for d in [3, 5, 8]:
        roots = sorted(rng.uniform(-2, 2, d))
        P = np.poly(roots)
        xs = np.linspace(-3, 3, 500)
        Ts = [turan_form(P, x) for x in xs]
        print(f"  d={d}: min T(x) on x in [-3,3] = {min(Ts):.4e}")

    print("\n[2] Polynomial with complex roots: T can be negative")
    P = np.poly([1+1j, 1-1j, 0.5])
    xs = np.linspace(-3, 3, 500)
    Ts = [turan_form(P, x).real for x in xs]
    print(f"  P with complex roots: min T = {min(Ts):.4f}, max T = {max(Ts):.4f}")

    print("\n[3] Double root: T = 0 at double root")
    P = np.poly([1, 1, 2.0])  # double at 1
    print(f"  P = (z-1)^2 (z-2): T(1) = {turan_form(P, 1.0):.4e}, T(2) = {turan_form(P, 2.0):.4e}")


if __name__ == "__main__":
    main()
