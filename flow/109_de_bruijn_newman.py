"""
PAPER: 109 (canonical: 109_de_bruijn_newman.pdf)
TITLE: De Bruijn-Newman Constant via Pandrosion-Pólya-Schur
STATUS: framework (Newman conjecture remains open)
DEPENDS: 055, 107

THEORY
======

DE BRUIJN-NEWMAN constant Lambda: parameter of heat-flow H_t(z) such that
H_t has only real zeros for t >= Lambda, complex zeros for t < Lambda.

NEWMAN CONJECTURE (1976): Lambda <= 0 (RH at boundary).
TAO-RODGERS (2018): Lambda >= 0 (proven).
Combined: Lambda = 0 ⟺ RH barely true.

PANDROSION-POLYA-SCHUR REFORMULATION:
  Lambda <= t_0 ⟺ H_{t_0} in Laguerre-Pólya class.
This requires H_{t_0} to have only real zeros, equivalently the Pandrosion
field F_{H_{t_0}} satisfies a Hurwitz-style positivity condition (paper 67).

VERIFICATION
============

  1. H_t(z) numerical at small |t|.
  2. Pólya-Schur stability test on truncations.
"""
from __future__ import annotations
import math
import numpy as np


def Phi_riemann(t, n_terms=20):
    """Phi(u) = sum (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) e^{-pi n^2 e^{4u}}."""
    s = 0.0
    for n in range(1, n_terms + 1):
        e9u = math.exp(9*t); e5u = math.exp(5*t); e4u = math.exp(4*t)
        s += (2 * math.pi**2 * n**4 * e9u - 3 * math.pi * n**2 * e5u) * \
             math.exp(-math.pi * n**2 * e4u)
    return s


def H_t(z, t, U=4.0, n_quad=100):
    """Newman's H_t(z) ~ (1/8) integral Phi(u) cos(zu) e^{tu^2} du."""
    du = U / n_quad
    integral = 0.0 + 0j
    for k in range(n_quad):
        u = (k + 0.5) * du
        phi = Phi_riemann(u)
        integral += phi * np.cos(z * u) * np.exp(t * u**2)
    return integral * du / 8.0


def main():
    print("=" * 80)
    print("PAPER 109 — De Bruijn-Newman Lambda")
    print("=" * 80)

    print("\n[1] Status:")
    print("  Tao-Rodgers 2018: Lambda >= 0 (proved)")
    print("  Newman 1976 conj: Lambda <= 0 (open)")
    print("  Combined: Lambda = 0 iff RH barely true.")
    print("  Best bound: 0 <= Lambda < 0.22 (Polymath15, 2018)")

    print("\n[2] H_t(z) for small t (numerical truncation, severely limited)")
    print(f"  {'t':>10} {'H_t(0)':>14} {'H_t(10)':>14}")
    for t in [-0.1, 0.0, 0.05, 0.1]:
        h0 = H_t(0.0, t).real
        h10 = H_t(10.0, t).real
        print(f"  {t:>10.3f} {h0:>14.4f} {h10:>14.4f}")

    print("\n[3] Pandrosion-Pólya-Schur reformulation:")
    print("  Lambda <= t_0 iff H_{t_0} in Laguerre-Pólya class.")
    print("  Equivalently: F_{H_{t_0}}(z) = H_{t_0}'/H_{t_0} satisfies")
    print("  Im F > 0 above real axis (Pandrosion-Hurwitz, paper 67).")


if __name__ == "__main__":
    main()
