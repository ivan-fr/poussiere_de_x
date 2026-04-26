"""
PAPER: 055 (canonical: 45pandrosion_riemann.pdf)
TITLE: Riemann zeta and Pandrosion (initial)
STATUS: framework (refined in papers 107, 110)
DEPENDS: 011

THEORY
======

Initial framework: Riemann xi(s) is entire of order 1, even after the change
of variable s -> 1/2 + iz. Its zeros (= Riemann zeros) are conjecturally
all on the critical line.

PANDROSION FORMULATION:
  F_xi(z) = xi'(z)/xi(z) = sum_rho 1/(z - rho)
  T_xi(z) = (xi'(z))^2 - xi(z) xi''(z)  (Turán form)

By Pandrosion-Turán (paper 30, 56): if xi has only real zeros, then
T_xi(z) >= 0 on R.

VERIFICATION
============

  1. xi(1/2 + i t) at the first few zeros.
  2. Turán form T_xi at sample points.
"""
from __future__ import annotations
import math


def main():
    print("=" * 80)
    print("PAPER 55 — Riemann zeta-Pandrosion (initial)")
    print("=" * 80)

    try:
        from mpmath import mp, mpc, pi, gamma, zeta, zetazero
        mp.dps = 30
    except ImportError:
        print("\n  [mpmath not available]")
        return

    print("\n[1] First 5 Riemann zeros (Im part on critical line)")
    for k in range(1, 6):
        rho = zetazero(k)
        print(f"  rho_{k} = {complex(rho)}")

    print("\n[2] xi(1/2 + i*Im(rho)) ~ 0 at zeros")
    def xi(s):
        s_mp = mpc(complex(s))
        return complex(0.5 * s_mp * (s_mp - 1) * pi**(-s_mp/2) * gamma(s_mp/2) * zeta(s_mp))

    for k in [1, 2, 3]:
        rho = zetazero(k)
        s = 0.5 + 1j * float(rho.imag)
        v = xi(s)
        print(f"  xi(1/2 + i*{float(rho.imag):.4f}) = {v:.4e}")

    print("\n[3] Turán form: T_xi(t) = (xi')^2 - xi * xi'' on critical line")

    def xi_critical(t):
        return xi(0.5 + 1j * t).real  # should be real on critical line

    def turan_xi(t, h=1e-3):
        f0 = xi_critical(t)
        fp = xi_critical(t + h)
        fm = xi_critical(t - h)
        f1 = (fp - fm) / (2 * h)
        f2 = (fp - 2*f0 + fm) / (h*h)
        return f1**2 - f0 * f2

    for t in [0.0, 5.0, 10.0, 14.0, 20.0]:
        T = turan_xi(t)
        print(f"  T_xi({t}) = {T:.4e}  (should be >= 0 if RH)")


if __name__ == "__main__":
    main()
