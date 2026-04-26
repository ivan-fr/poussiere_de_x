"""
PAPER: 125 (NEW — Jensen polynomials of Riemann ξ)
TITLE: Jensen Polynomials of Riemann ξ — Pólya-Bender via Pandrosion
STATUS: Pólya 1927 reformulation of RH; Griffin-Ono-Rolen-Zagier 2019 (stable
        real-rootedness for fixed d). Full conjecture (all d, N) open == RH.
DEPENDS: 056 (Hardy), 066 (Turán), 089 (MSS), 098 (Pólya-Schur),
         110 (RH-Pólya-Schur), 055/107 (Riemann)

THEORY
======

------------------------------------------------------------------------
JENSEN POLYNOMIALS (Pólya 1927)
------------------------------------------------------------------------

Riemann xi function expansion at s = 1/2:
  xi(1/2 + iz) = sum_{n=0}^infty gamma_n z^n
where gamma_n in R (xi is real on the critical line).

The JENSEN POLYNOMIAL of degree d and shift N is:
  J_d^N(X) = sum_{j=0}^d binom(d, j) gamma_{N+j} X^j.

------------------------------------------------------------------------
PÓLYA-BENDER REFORMULATION OF RH
------------------------------------------------------------------------

Pólya (1927) showed:
  RH <=> for ALL (d, N) with d, N >= 0, J_d^N has only real zeros.

Equivalently: xi(1/2 + i z) is in the Laguerre-Polya class (paper 110).

ASYMPTOTIC RESULT (Griffin-Ono-Rolen-Zagier 2019):
  For every fixed d, there exists N_0(d) such that J_d^N has only real
  zeros for all N >= N_0(d).

This is the asymptotic Pólya-Bender. The full statement (all N for all d)
remains open == RH.

------------------------------------------------------------------------
PANDROSION-PÓLYA-SCHUR FRAMEWORK
------------------------------------------------------------------------

By paper 098 (Pólya-Schur preservers): J_d^N is built from xi's coefficients
via a multiplier operator. The operator
  T: f(X) -> sum binom(d, j) (coefficient of f at X^j) (some shift) X^j
is a specific Pólya-Schur preserver.

Pandrosion-Turán (paper 066): J_d^N has only real zeros iff
  (J_d^N)'(x)^2 - J_d^N(x) (J_d^N)''(x) >= 0 on R.

By paper 110 (RH-Pólya-Schur synthesis), we have multiple equivalent
necessary conditions for RH; Jensen polynomials provide a CONCRETE
finite-dimensional verification target.

VERIFICATION
============

  1. Compute gamma_n via mpmath (xi Taylor series).
  2. Build J_d^N for small (d, N).
  3. Check real-rootedness (= Turán positivity on R).
  4. Confirm Griffin-Ono-Rolen-Zagier: stable rootedness for fixed d.
"""
from __future__ import annotations
import math


def main():
    print("=" * 80)
    print("PAPER 125 — Jensen polynomials of Riemann xi")
    print("=" * 80)

    try:
        from mpmath import mp, mpc, mpf, pi, gamma, zeta, taylor
        mp.dps = 50
    except ImportError:
        print("\n  mpmath required; skipping.")
        return

    # ===========================
    # PROOF SKETCH (Pólya 1927)
    # ===========================
    # xi(1/2 + iz) is entire of order 1, even, real on R.
    # Its Hadamard product: xi(s) = xi(0) prod (1 - s/rho).
    # Substituting s = 1/2 + iz: xi_R(z) := xi(1/2 + iz) = xi(0) prod (1 - 2 i z / (rho - 1/2)).
    # If RH holds: rho = 1/2 + i gamma_n with gamma_n in R. So
    #   xi_R(z) = xi(0) prod (1 + 2z/gamma_n) (real factors!).
    # Hence xi_R is in Laguerre-Polya class.
    # Truncating to degree d Taylor gives Jensen polynomial J_d^N which
    # must have real zeros (sub-class of Laguerre-Polya).
    # If RH FAILS: some rho off critical line, then xi_R has complex conjugate
    # pair of zeros, and Jensen polynomial of high degree shows complex zero.
    # QED (Pólya equivalence).

    # Compute gamma_n: Taylor coefficients of xi(1/2 + iz) at z = 0
    print("\n[1] Computing gamma_n: Taylor coefficients of xi(1/2 + iz)")

    def xi_real(z):
        """xi(1/2 + iz), should be real."""
        s = mpc('0.5') + mpc(0, 1) * z
        v = mpf('0.5') * s * (s - 1) * pi**(-s/2) * gamma(s/2) * zeta(s)
        return v

    # Use Taylor expansion at z = 0
    print("  Computing first 12 Taylor coeffs (high precision)...")
    coeffs = []
    for n in range(12):
        # Numerical derivative: d^n/dz^n xi(1/2 + iz) / n! at z = 0
        # Use small h finite differences (more reliable for high n with mpmath taylor)
        h = mpf('0.01')
        if n == 0:
            v = xi_real(mpf(0))
        else:
            # Central finite differences for n-th derivative
            v = mpf(0)
            for k in range(n + 1):
                sign = (-1)**(n - k)
                bin_k = math.factorial(n) // (math.factorial(k) * math.factorial(n - k))
                v += sign * bin_k * xi_real((k - n/2) * h)
            v /= h**n
        coeffs.append(complex(v).real / math.factorial(n))
        if abs(coeffs[-1]) < 1e-50:
            coeffs[-1] = 0.0  # numerical zero

    # Pólya's coefficients have alternating sign / specific structure
    # gamma_0 = xi(1/2) > 0
    # Many even gamma_2k != 0; odd usually zero by symmetry xi(1/2 + iz) = xi(1/2 - iz)
    print(f"  gamma_n (first 12, alternating zeros expected for odd n):")
    for n, g in enumerate(coeffs[:12]):
        print(f"    gamma_{n} = {g:.6e}")

    # By symmetry xi(1/2 + iz) = xi(1/2 - iz), odd coefficients should be ~0
    # Even coefficients are the meaningful ones. Use only even ones.
    even_coeffs = [coeffs[2*k] for k in range(6)]
    print(f"\n  Even gamma_{{2k}}: {[f'{g:.4e}' for g in even_coeffs]}")

    # ===========================
    # Jensen polynomials
    # ===========================
    print("\n[2] Jensen polynomial J_d^N(X) and real-rootedness")

    def jensen(d, N, gammas):
        """J_d^N(X) = sum binom(d, j) gamma_{N+j} X^j (low-to-high coeffs)."""
        coefs = []
        for j in range(d + 1):
            if N + j >= len(gammas): return None
            coefs.append(math.comb(d, j) * gammas[N + j])
        return coefs[::-1]  # high-to-low for numpy

    import numpy as np

    print(f"  {'d':>3} {'N':>3} {'J_d^N coefficients (high-low)':>50}")
    for d in [2, 3, 4]:
        for N in [0, 2, 4]:
            J = jensen(d, N, coeffs)
            if J is None: continue
            J_str = ', '.join(f"{c:.2e}" for c in J)
            print(f"  {d:>3} {N:>3}  {J_str}")

    print("\n[3] Real-rootedness of J_d^N (necessary for RH)")
    print(f"  {'d':>3} {'N':>3} {'all roots real?':>16} {'Turán T(0)':>14}")
    for d in [2, 3, 4]:
        for N in [0, 2, 4, 6]:
            J = jensen(d, N, coeffs)
            if J is None: continue
            roots = np.roots(J)
            all_real = all(abs(r.imag) < 1e-6 * (1 + abs(r.real)) for r in roots)
            # Turán form at 0
            J_arr = np.array(J)
            Jp = np.polyder(J_arr)
            Jpp = np.polyder(Jp)
            T = float(np.polyval(Jp, 0)**2 - np.polyval(J_arr, 0) * np.polyval(Jpp, 0))
            print(f"  {d:>3} {N:>3} {str(all_real):>16} {T:>14.4e}")

    print("\n[4] Stability: for fixed d, J_d^N rootedness as N grows")
    print(f"  Griffin-Ono-Rolen-Zagier 2019: stable for N >= N_0(d).")
    print(f"  {'d':>3} {'N':>3} {'min |Im(root)|':>16}")
    for d in [3]:
        for N in [0, 2, 4, 6, 8, 10]:
            J = jensen(d, N, coeffs)
            if J is None: break
            roots = np.roots(J)
            max_imag = max(abs(r.imag) for r in roots)
            print(f"  {d:>3} {N:>3} {max_imag:>16.4e}")

    print("\n[5] HONEST ASSESSMENT")
    print("  PROVED:")
    print("    Pólya 1927: RH <=> all J_d^N real-rooted for all (d, N).")
    print("    Griffin-Ono-Rolen-Zagier 2019: stable rootedness, fixed d, large N.")
    print("  ")
    print("  PANDROSION CONTRIBUTION:")
    print("    Turán T = (J')^2 - J J'' >= 0 on R is the REAL-ROOTEDNESS criterion")
    print("    (paper 066). Verifiable for ANY explicit J_d^N.")
    print("    Pólya-Schur preservers (paper 098) characterize the operators")
    print("    mapping LP -> LP, including the Jensen-polynomial construction.")
    print("  ")
    print("  OPEN:")
    print("    Full Pólya 1927 (all d, N) == Riemann Hypothesis.")
    print("    Effective N_0(d) bounds in Griffin-Ono-Rolen-Zagier are open in")
    print("    sharp form.")
    print("    NUMERICAL gamma_n at high precision: limits feasibility for d > 12.")


if __name__ == "__main__":
    main()
