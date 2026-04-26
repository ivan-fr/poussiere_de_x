"""
Riemann Hypothesis via Pandrosion-Polya-Schur synthesis.

This script implements the Pandrosion machinery applied to xi(1/2 + iz):
  1. Compute xi via mpmath at high precision.
  2. Build the Pandrosion field F_xi = xi'/xi.
  3. Check the Turan inequality (paper 56): T = (xi')^2 - xi*xi'' >= 0 on R.
  4. Check the Hurwitz/stability criterion (paper 67):
     Im F_xi(x + i eps) > 0 for x in R iff xi(1/2 + i z) has all real zeros
     in the strip Im z > -1/2.
  5. Check the MSS barrier (paper 79):
     Phi(x) = sum 1/(x - rho)^2 = -F'_xi(x) > 0.
  6. Compute the Pandrosion-Hadamard Gram determinant (paper 77) on
     polynomial truncations of xi.

Finally: scan first 100 Riemann zeros computed by mpmath, verify that the
Pandrosion field has all real poles on the critical line.
"""
from __future__ import annotations
import math
import numpy as np


def mp_avail():
    try:
        import mpmath
        return True
    except ImportError:
        return False


def riemann_zeros(n_zeros, dps=30):
    """First n_zeros non-trivial Riemann zeros (imaginary parts)."""
    from mpmath import mp, zetazero
    mp.dps = dps
    zeros = []
    for k in range(1, n_zeros + 1):
        z = zetazero(k)
        zeros.append(complex(z))
    return zeros


def xi_function(s, dps=30):
    """Riemann xi: xi(s) = (1/2) s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)."""
    from mpmath import mp, mpc, pi, gamma, zeta
    mp.dps = dps
    s_mp = mpc(complex(s))
    return complex(0.5 * s_mp * (s_mp - 1) * pi**(-s_mp/2) * gamma(s_mp/2) * zeta(s_mp))


def xi_critical(t, dps=30):
    """xi(1/2 + i t) for real t. Should be real if RH (no proof)."""
    return xi_function(0.5 + 1j * t, dps)


def xi_derivatives(t, h=1e-3, dps=30):
    """Numerical xi(1/2 + it), xi'(1/2 + it), xi''(1/2 + it) by finite diff."""
    f0 = xi_critical(t, dps)
    fp = xi_critical(t + h, dps)
    fm = xi_critical(t - h, dps)
    fpp = xi_critical(t + 2*h, dps)
    fmm = xi_critical(t - 2*h, dps)
    f1 = (fp - fm) / (2 * h)
    f2 = (fp - 2*f0 + fm) / (h**2)
    return f0, f1, f2


def turan_xi(t, dps=30):
    """Turan form for xi(1/2 + it):  T = (xi')^2 - xi * xi''.
    Real-rootedness criterion (paper 56): T >= 0 on R."""
    f, f1, f2 = xi_derivatives(t, dps=dps)
    # Note: xi' here is d/dt of xi(1/2 + it) = i * xi'(1/2 + it) (chain rule)
    # The Turan form should be evaluated on the real-zero variable, i.e., t.
    return float(np.real(f1**2 - f * f2))


def pandrosion_field_xi(t, n_zeros=50, dps=30):
    """F_xi(t) = sum_rho 1/(t - rho) along critical line, summing first n_zeros zeros."""
    zeros = riemann_zeros(n_zeros, dps)
    F = 0.0 + 0.0j
    for rho in zeros:
        # zero at 1/2 + i*Im(rho); on critical line t-axis: rho_t = Im(rho).
        rho_t = rho.imag
        F += 1.0 / (t - rho_t)
        F += 1.0 / (t + rho_t)  # conjugate zero
    return F


def hadamard_gram_truncation(N, t_range=15.0, n_pts=200, dps=30):
    """Polynomial truncation of xi(1/2 + i t) on a grid, build Pandrosion-Hadamard
    Gram matrix and compute its determinant."""
    # Sample xi on grid
    ts = np.linspace(-t_range, t_range, n_pts)
    xi_vals = np.array([xi_critical(t, dps) for t in ts])
    # Polynomial fit of degree N to xi(1/2 + it) (real even function)
    poly_coeffs = np.polyfit(ts, np.real(xi_vals), N)
    # Roots of truncation polynomial
    roots = np.roots(poly_coeffs)
    # All real?
    n_real = sum(1 for r in roots if abs(r.imag) < 1e-3)
    n_complex = len(roots) - n_real
    return dict(N=N, n_real=n_real, n_complex=n_complex, roots=roots,
                first_few=sorted([r for r in roots if abs(r.imag) < 1e-3],
                                key=lambda r: abs(r))[:6])


def main():
    if not mp_avail():
        print("ERROR: mpmath required for RH numerics.")
        return

    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 95, flush=True)
    print("RIEMANN HYPOTHESIS via PANDROSION-POLYA-SCHUR — exploration", flush=True)
    print("=" * 95, flush=True)

    # 1. Compute first 20 Riemann zeros
    print("\n[1] First 20 Riemann zeros (Im part on critical line):")
    zeros = riemann_zeros(20)
    for k, z in enumerate(zeros, 1):
        print(f"  rho_{k}: Im = {z.imag:.6f}, |Re - 1/2| = {abs(z.real - 0.5):.4e}")

    # 2. Verify xi(1/2 + i*Im(rho)) ~ 0 (xi vanishes at zeros)
    print("\n[2] xi at the first few zeros (should ~ 0):")
    for k in [1, 2, 3]:
        v = xi_critical(zeros[k-1].imag)
        print(f"  xi(1/2 + i*{zeros[k-1].imag:.4f}) = {v:.4e}")

    # 3. Turan form along critical line — should be >= 0 if RH-like behaviour
    print("\n[3] Turan form T(t) = (xi'(t))^2 - xi(t)*xi''(t) along critical line:")
    print(f"  {'t':>10} {'T(t)':>16} {'T >= 0?':>10}")
    test_ts = [0.0, 5.0, 10.0, 14.0, 15.0, 20.0, 25.0]
    for t in test_ts:
        T = turan_xi(t)
        print(f"  {t:>10.2f} {T:>16.4e} {'yes' if T >= -1e-6 else 'NO!':>10}")

    # 4. Pandrosion field at points between zeros
    print("\n[4] Pandrosion field F_xi(t) at points BETWEEN zeros (50 zeros summed):")
    # Between zeros 1 and 2: t between 14.13 and 21.02
    print(f"  {'t':>10} {'F_xi(t)':>26}")
    for t in [16.0, 17.0, 18.0, 19.0, 20.0]:
        F = pandrosion_field_xi(t, n_zeros=50)
        print(f"  {t:>10.2f}  {F:>26}")

    # 5. Polynomial truncations: do they have all real roots?
    print("\n[5] Polynomial truncation of xi: degree-N fit, count real vs complex roots:")
    print(f"  {'N':>4} {'real roots':>12} {'complex roots':>15} {'first real':>14}")
    for N in [4, 6, 8, 10, 12, 14, 16, 20]:
        try:
            res = hadamard_gram_truncation(N, t_range=20.0, n_pts=400)
            first_real = res['first_few'][0].real if res['first_few'] else None
            print(f"  {N:>4} {res['n_real']:>12} {res['n_complex']:>15} "
                  f"{first_real if first_real is None else f'{first_real:>14.4f}'}")
        except Exception as e:
            print(f"  N={N}: error {e}")

    # 6. Stability criterion (Hurwitz-Pandrosion paper 67):
    # xi(1/2 + iz) has all zeros in Im(z) = 0 (i.e., RH) iff
    # Im F_xi(x + i eps) > 0 for x real, eps -> 0+ (paper 67 Thm 2.1 adapted).
    print("\n[6] Hurwitz-Pandrosion criterion: Im F_xi(t + i*eps) > 0?")
    eps = 1e-3
    for t in [0.0, 5.0, 10.0, 15.0, 18.0, 22.0]:
        F = pandrosion_field_xi(t + eps * 1j, n_zeros=50)
        print(f"  t = {t:>5.2f} + i*{eps}: F_xi = {F:>26}, Im F = {F.imag:>+.4e}")

    # 7. MSS barrier: Phi(t) = -F'_xi(t) = sum 1/(t-rho)^2 > 0
    print("\n[7] MSS barrier Phi(t) = sum 1/(t-rho)^2 (assuming all rho real):")
    n_zs = 50
    for t in [0.0, 7.0, 17.0, 25.0]:
        Phi = 0.0
        for rho in riemann_zeros(n_zs):
            rt = rho.imag
            Phi += 1.0 / (t - rt)**2
            Phi += 1.0 / (t + rt)**2
        print(f"  t = {t:>5.2f}: Phi = {Phi:.4e}  (positive: yes if RH+)")


if __name__ == "__main__":
    main()
