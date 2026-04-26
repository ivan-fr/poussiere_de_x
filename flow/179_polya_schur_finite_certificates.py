"""
PAPER: 179 (NEW — finite Pólya-Schur/Jensen certificates)
TITLE: Riemann via finite Pólya-Schur certificates: Jensen-polynomial
       hyperbolicity, Hermite matrices, and Pandrosion-Turán margins
STATUS: RH remains OPEN.  This paper builds a finite, reproducible
        certificate pipeline for the Pólya-Schur/Laguerre-Pólya door:
        every certified Jensen polynomial is evidence for RH, not a proof.
DEPENDS: 098 (Pólya-Schur), 110 (RH-Pólya-Schur synthesis),
         125 (Jensen polynomials), 135 (Pandrosion-Turán quantitative).

THEORY
======

------------------------------------------------------------------------
THE DOOR
------------------------------------------------------------------------

Let
    Xi(t) := xi(1/2 + i t),
where xi is Riemann's entire xi-function.  Xi is real and even on R.

Pólya's reformulation:
    RH  <=>  Xi lies in the Laguerre-Pólya class
        <=>  every Jensen polynomial built from Xi's Taylor coefficients
             has only real zeros.

For Taylor coefficients
    Xi(t) = sum_{n >= 0} gamma_n t^n,
the Jensen polynomial is
    J_d^N(X) = sum_{j=0}^d binom(d, j) gamma_{N+j} X^j.

Because Xi is even, odd gamma_n vanish.  Therefore the raw all-index
Jensen polynomials are frequently degenerate.  The useful finite tests use
the even subsequence
    A_m := (-1)^m gamma_{2m}
so that
    Xi(t) = sum_m (-1)^m A_m t^{2m}.

Then
    J_d^N(X) = sum_{j=0}^d binom(d, j) A_{N+j} X^j.

------------------------------------------------------------------------
FINITE CERTIFICATES
------------------------------------------------------------------------

For a real polynomial P:

  (C1) root certificate:
       all roots have imaginary part <= tolerance.

  (C2) Hermite matrix certificate:
       H(P) = (s_{i+j})_{0 <= i,j < d}, where s_k is the k-th power
       sum of the roots.  P has all real roots iff H(P) is positive
       semidefinite and full rank (numerically: min eigenvalue >= 0).

  (C3) Pandrosion-Turán grid certificate:
       T_P(x) := P'(x)^2 - P(x) P''(x) >= 0 on a real grid.
       This is necessary for real-rootedness and exact for polynomials
       with all real simple roots.

PAPER 179 contribution:
  - compute A_m from Xi derivatives;
  - build finite Jensen polynomials;
  - emit a compact certificate table (roots, Hermite min eigenvalue,
    Turán grid margin);
  - identify the next obstruction: turning finite certificates into
    uniform-in-(d,N) inequalities.
"""
from __future__ import annotations

import math


def _poly_eval_low(coeffs, x):
    """Evaluate a polynomial stored low-to-high."""
    total = 0.0
    for c in reversed(coeffs):
        total = total * x + c
    return total


def _poly_derivative_low(coeffs):
    return [j * coeffs[j] for j in range(1, len(coeffs))]


def _jensen_low(A, d, N):
    """J_d^N(X) = sum binom(d,j) A[N+j] X^j, low-to-high."""
    if N + d >= len(A):
        return None
    return [math.comb(d, j) * A[N + j] for j in range(d + 1)]


def _trim_low(coeffs, tol=1e-80):
    out = list(coeffs)
    while len(out) > 1 and abs(out[-1]) < tol:
        out.pop()
    return out


def _turan_margin(coeffs, radius=20.0, n_grid=2001):
    """Minimum of P'^2 - P P'' on a symmetric grid."""
    p = _trim_low(coeffs)
    p1 = _poly_derivative_low(p)
    p2 = _poly_derivative_low(p1)
    min_T = float("inf")
    argmin = 0.0
    for i in range(n_grid):
        x = -radius + 2 * radius * i / (n_grid - 1)
        v = _poly_eval_low(p1, x) ** 2 - _poly_eval_low(p, x) * _poly_eval_low(p2, x)
        if v < min_T:
            min_T = v
            argmin = x
    return min_T, argmin


def main():
    print("=" * 80)
    print("PAPER 179 — finite Pólya-Schur/Jensen certificates")
    print("=" * 80)

    try:
        import numpy as np
        from mpmath import gamma, mp, mpc, pi, zeta
    except ImportError:
        print("\n  [numpy + mpmath required]")
        return

    mp.dps = 50

    def xi(s):
        s = mpc(s)
        return mp.mpf("0.5") * s * (s - 1) * pi ** (-s / 2) * gamma(s / 2) * zeta(s)

    def Xi(t):
        return mp.re(xi(mpc("0.5", t)))

    # ------------------------------------------------------------------
    # [1] Coefficients A_m = (-1)^m gamma_{2m}
    # ------------------------------------------------------------------
    print("\n[1] Taylor coefficients A_m = (-1)^m gamma_{2m}")
    print("    Fast Cauchy/FFT extraction from Xi(r exp(i theta)).")

    # Cauchy extraction is much faster than high-order automatic
    # differentiation here.  For Xi(t) = sum c_n t^n, sample on a circle
    # and use the discrete Cauchy formula c_n ~= mean Xi(r w_j) w_j^{-n}/r^n.
    max_m = 14
    n_coeff = 2 * max_m
    n_fft = 96
    radius = mp.mpf("0.55")
    samples = []
    for j in range(n_fft):
        theta = 2 * pi * j / n_fft
        w = mp.e ** (mpc(0, theta))
        samples.append(xi(mpc("0.5") + mpc(0, 1) * radius * w))

    coeffs = []
    for n in range(n_coeff + 1):
        total = mpc(0)
        for j, val in enumerate(samples):
            theta = 2 * pi * j / n_fft
            total += val * (mp.e ** (-mpc(0, n * theta)))
        coeffs.append(total / n_fft / (radius ** n))

    A = []
    for m in range(max_m + 1):
        gamma_2m = mp.re(coeffs[2 * m])
        A_m = ((-1) ** m) * gamma_2m
        A.append(float(A_m))
        print(f"  A_{m:02d} = {float(A_m): .12e}", flush=True)

    print("\n  Sanity: A_m should be positive for the tested prefix.")
    print(f"  min A_m = {min(A):.6e}")

    # ------------------------------------------------------------------
    # [2] Finite Jensen hyperbolicity table
    # ------------------------------------------------------------------
    print("\n[2] Jensen hyperbolicity certificates")
    print(f"  {'d':>3} {'N':>3} {'max |Im root|':>15} {'min eig H':>14}"
          f" {'min Turan':>14} {'status':>10}")

    passed = 0
    total = 0
    worst_imag = 0.0
    worst_eig = float("inf")
    worst_turan = float("inf")

    for d in range(2, 8):
        for N in range(0, 8):
            coeffs_low = _jensen_low(A, d, N)
            if coeffs_low is None:
                continue
            coeffs_low = _trim_low(coeffs_low)
            if len(coeffs_low) != d + 1:
                continue

            coeffs_high = np.array(list(reversed(coeffs_low)), dtype=float)
            roots = np.roots(coeffs_high)
            max_imag = float(max(abs(r.imag) for r in roots))

            # Hermite matrix from numerical roots: H_ij = sum r^(i+j).
            # The raw matrix is badly conditioned when the roots are large
            # (Jensen roots routinely sit at scale 10^2-10^4).  Scaling roots
            # by their max modulus preserves real-rootedness and makes the
            # PSD diagnostic meaningful in double precision.
            scale = max(1.0, float(max(abs(r) for r in roots)))
            roots_scaled = roots / scale
            power_sums = []
            for k in range(2 * d - 1):
                power_sums.append(float(np.real(np.sum(roots_scaled ** k))))
            H = np.array([[power_sums[i + j] for j in range(d)] for i in range(d)])
            eigs = np.linalg.eigvalsh((H + H.T) / 2.0)
            min_eig = float(eigs[0])

            # Turán margin on a radius adapted to root scale.
            root_radius = max(1.0, float(max(abs(r.real) for r in roots)))
            min_T, _ = _turan_margin(coeffs_low, radius=1.5 * root_radius)

            ok = max_imag < 1e-7 and min_eig > -1e-7 and min_T > -1e-9
            status = "CERT" if ok else "CHECK"
            total += 1
            passed += int(ok)
            worst_imag = max(worst_imag, max_imag)
            worst_eig = min(worst_eig, min_eig)
            worst_turan = min(worst_turan, min_T)
            print(f"  {d:>3} {N:>3} {max_imag:>15.3e} {min_eig:>14.3e}"
                  f" {min_T:>14.3e} {status:>10}")

    # ------------------------------------------------------------------
    # [3] A concrete certificate polynomial
    # ------------------------------------------------------------------
    print("\n[3] Example certificate object: J_6^2")
    d, N = 6, 2
    coeffs_low = _jensen_low(A, d, N)
    if coeffs_low is not None:
        roots = np.roots(np.array(list(reversed(coeffs_low)), dtype=float))
        print("  coefficients low-to-high:")
        print("   ", " ".join(f"{c:.8e}" for c in coeffs_low))
        print("  roots:")
        for r in sorted(roots, key=lambda z: z.real):
            print(f"    {r.real: .10e} {r.imag:+.2e}i")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print(f"  finite certificates passed: {passed}/{total}")
    print(f"  worst max imaginary part:   {worst_imag:.3e}")
    print(f"  worst Hermite min eig:      {worst_eig:.3e}")
    print(f"  worst Turan grid margin:    {worst_turan:.3e}")
    print()
    print("  WHAT THIS DOES:")
    print("    - Produces reproducible finite Pólya-Schur/Jensen certificates.")
    print("    - Gives a concrete target for Lean/interval-arithmetic upgrades.")
    print("    - Fits the RH door from paper 110 exactly.")
    print()
    print("  WHAT REMAINS OPEN:")
    print("    - Uniform proof for all d and N, equivalent to RH.")
    print("    - Rigorous interval certificates replacing floating-point roots.")
    print("    - Effective bounds extending Griffin-Ono-Rolen-Zagier from")
    print("      fixed-d asymptotics to a usable all-(d,N) theorem.")


if __name__ == "__main__":
    main()
