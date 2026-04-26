"""
PAPER: 107 (canonical: 107_riemann_pandrosion.pdf)
TITLE: Riemann zeta and Pandrosion (continuation)
STATUS: framework (Pandrosion-form of zeta-related identities)
DEPENDS: 055

THEORY
======

Riemann xi function on critical line: xi(1/2 + i t) is real-valued.
Pandrosion field F_xi = xi'/xi has poles at Riemann zeros.

KEY IDENTITY: At Riemann zeros rho_k, the Pandrosion residue is 1.

VERIFICATION
============

  1. xi at first 10 Riemann zeros.
  2. F_xi field structure.
"""
from __future__ import annotations


def main():
    print("=" * 80)
    print("PAPER 107 — Riemann zeta + Pandrosion")
    print("=" * 80)

    try:
        from mpmath import mp, mpc, pi, gamma, zeta, zetazero
        mp.dps = 30
    except ImportError:
        print("  mpmath not available; skipping.")
        return

    print("\n[1] First 10 Riemann zeros (Im part on critical line)")
    for k in range(1, 11):
        rho = zetazero(k)
        print(f"  rho_{k} = 1/2 + i*{float(rho.imag):.6f}")

    def xi(s):
        s_mp = mpc(complex(s))
        return complex(0.5 * s_mp * (s_mp - 1) * pi**(-s_mp/2) * gamma(s_mp/2) * zeta(s_mp))

    print("\n[2] xi(1/2 + i*Im(rho_k)) ~ 0")
    for k in [1, 2, 5, 10]:
        rho = zetazero(k)
        s = 0.5 + 1j * float(rho.imag)
        v = xi(s)
        print(f"  xi(1/2 + i*{float(rho.imag):.4f}) = {v:.4e}")

    print("\n[3] Pandrosion field F_xi(s) = xi'(s)/xi(s)")
    h = 1e-3
    for t in [5.0, 10.0, 18.0]:
        s0 = 0.5 + 1j * t
        # Numerical derivative of xi via finite diff
        xi_p = xi(s0 + h * 1j)
        xi_m = xi(s0 - h * 1j)
        xi_0 = xi(s0)
        d_xi = (xi_p - xi_m) / (2 * h * 1j)
        F = d_xi / xi_0 if abs(xi_0) > 1e-30 else float('inf')
        print(f"  F_xi(1/2 + i*{t}) = {F:.4f}")


if __name__ == "__main__":
    main()
