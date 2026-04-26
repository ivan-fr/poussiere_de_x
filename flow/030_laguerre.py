"""
PAPER: 030 (canonical: 20pandrosion_laguerre.pdf)
TITLE: Laguerre's Inequality from Pandrosion Energy
STATUS: proved (Laguerre 1885, Pandrosion-energy form)
DEPENDS: 011, 056

THEORY
======

LAGUERRE'S INEQUALITY: For real-rooted polynomial P (with simple roots),
  (P'(x))^2 - P(x) P''(x) >= 0  for all real x.

PANDROSION FORM: This is exactly E_P(x) = sum Q(alpha_k, x)^2 >= 0,
the strict positivity of Pandrosion energy on R.

Reverse direction: if E_P(x) >= 0 for all real x, then P has only real roots
(Turán-Pandrosion criterion, paper 56).

VERIFICATION
============

  1. E_P(x) >= 0 on R for real-rooted P.
  2. E_P(x) can be negative for non-real-rooted P.
"""
from __future__ import annotations
import math
import numpy as np


def E(P, x):
    Pp = np.polyder(P)
    Ppp = np.polyder(Pp)
    return np.polyval(Pp, x)**2 - np.polyval(P, x) * np.polyval(Ppp, x)


def main():
    print("=" * 80)
    print("PAPER 30 — Laguerre's inequality (Pandrosion-energy form)")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Real-rooted polynomials: E_P >= 0 on R")
    for d in [3, 5, 8]:
        roots = sorted(rng.uniform(-2, 2, d))
        P = np.poly(roots)
        xs = np.linspace(-3, 3, 100)
        Es = [E(P, x).real for x in xs]
        min_E = min(Es)
        print(f"  d={d}, roots in [-2,2]: min E_P(x) on x in [-3,3] = {min_E:.4e}  "
              f"({'OK' if min_E >= -1e-9 else 'NEGATIVE!'})")

    print("\n[2] Polynomial with complex roots: E_P can be negative")
    P = np.poly([1+1j, 1-1j, 0.5])  # 1 real, 2 complex conjugate
    xs = np.linspace(-3, 3, 200)
    Es = [E(P, x).real for x in xs]
    min_E = min(Es)
    print(f"  P with 1 real + complex pair: min E_P = {min_E:.4f}  "
          f"({'NEG (expected)' if min_E < 0 else 'pos'})")

    print("\n[3] At a real root alpha: E_P(alpha) = (P'(alpha))^2 >= 0")
    for d in [3, 5, 7]:
        roots = rng.uniform(-1, 1, d)
        P = np.poly(roots)
        for alpha in roots[:2]:
            E_at = E(P, alpha)
            Pp_sq = np.polyval(np.polyder(P), alpha)**2
            print(f"  d={d}, alpha={alpha:.3f}: E(alpha)={E_at:.4e}, P'(alpha)^2={Pp_sq:.4e}")


if __name__ == "__main__":
    main()
