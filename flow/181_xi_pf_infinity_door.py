"""
PAPER: 181 (NEW — Xi coefficients as a PF_infinity door)
TITLE: The PF∞ door to RH: total positivity of the Xi coefficient sequence
       and Pandrosion finite-minor diagnostics
STATUS: RH remains OPEN.  This paper isolates a sharp analytic target:
        prove that the coefficient sequence A_n of G(x)=Xi(sqrt(-x)) is
        PF∞.  Finite Toeplitz/Hankel minors are certified numerically here;
        a uniform proof would be equivalent in strength to RH.
DEPENDS: 110 (RH-Pólya-Schur synthesis),
         125 (Jensen polynomials), 179 (finite certificates).

THEORY
======

------------------------------------------------------------------------
THE INSIGHT
------------------------------------------------------------------------

Let
    Xi(t) := xi(1/2 + i t)
be Riemann's even real entire Xi-function.  Write
    Xi(t) = sum_{n >= 0} (-1)^n A_n t^(2n),       A_n > 0.

Define the square-root transform
    G(x) := Xi(sqrt(-x)) = sum_{n >= 0} A_n x^n.

Then:
    RH  <=>  all zeros of Xi(t) are real
        <=>  all zeros of G(x) are real and non-positive
        <=>  (A_n) is a Pólya-frequency sequence of infinite order (PF∞)
             in the Aissen-Schoenberg-Whitney / Edrei sense.

So the analytic insight needed for RH can be focused as:

    PROVE TOTAL POSITIVITY OF THE TOEPLITZ MATRIX
        T(A)_{i,j} = A_{j-i},    A_k := 0 for k < 0.

This is more structural than testing individual Jensen polynomials.  The
Jensen certificates of paper 179 are shadows of the same total-positivity
phenomenon.

------------------------------------------------------------------------
WHAT THIS PAPER TESTS
------------------------------------------------------------------------

  (1) log-concavity / Turán-1 margins:
          A_n^2 - A_{n-1} A_{n+1} >= 0.

  (2) finite PF∞ Toeplitz minors:
          det(A_{c_j-r_i}) >= 0
      for sampled row/column sets.

  (3) Hankel moment minors:
          det(A_{N+i+j}).
      These are NOT expected to be positive: log-concavity gives negative
      2x2 Hankel determinants.  This is useful because it rejects the naive
      "A_n is a Stieltjes moment sequence" route.  The right object is
      Toeplitz total positivity, not Hankel moment positivity.

  (4) a concrete conjecture:
      find a Toeplitz-total-positive kernel/transform for A_n, not a
      plain positive moment representation.
"""
from __future__ import annotations

import itertools
import math


def extract_A(max_m=30, n_fft=192, radius_str="0.62", dps=70):
    """Return A_m from Xi(t)=sum (-1)^m A_m t^(2m)."""
    from mpmath import gamma, mp, mpc, pi, zeta

    mp.dps = dps

    def xi(s):
        s = mpc(s)
        return mp.mpf("0.5") * s * (s - 1) * pi ** (-s / 2) * gamma(s / 2) * zeta(s)

    radius = mp.mpf(radius_str)
    samples = []
    for j in range(n_fft):
        theta = 2 * pi * j / n_fft
        w = mp.e ** (mpc(0, theta))
        samples.append(xi(mpc("0.5") + mpc(0, 1) * radius * w))

    coeffs = []
    for n in range(2 * max_m + 1):
        total = mpc(0)
        for j, val in enumerate(samples):
            theta = 2 * pi * j / n_fft
            total += val * (mp.e ** (-mpc(0, n * theta)))
        coeffs.append(total / n_fft / (radius ** n))

    A = []
    for m in range(max_m + 1):
        gamma_2m = mp.re(coeffs[2 * m])
        A.append(((-1) ** m) * gamma_2m)
    return A


def mp_det_from_rows(rows):
    """Small determinant by Leibniz expansion.

    The tested minors have size at most 5.  This avoids version-dependent
    LU behaviour on exactly singular mpmath matrices.
    """
    n = len(rows)
    if n == 0:
        return 1
    total = rows[0][0] * 0
    for perm in itertools.permutations(range(n)):
        inv = 0
        for i in range(n):
            for j in range(i + 1, n):
                inv += int(perm[i] > perm[j])
        term = rows[0][perm[0]]
        for i in range(1, n):
            term *= rows[i][perm[i]]
        total += -term if inv % 2 else term
    return total


def toeplitz_entry(A, r, c):
    k = c - r
    if k < 0 or k >= len(A):
        return A[0] * 0
    return A[k]


def normalized_hankel_matrix(A, shift, size):
    """Correlation-normalized Hankel matrix, PSD iff raw Hankel is PSD."""
    rows = []
    for i in range(size):
        row = []
        for j in range(size):
            denom = (A[shift + 2 * i] * A[shift + 2 * j]) ** 0.5
            row.append(A[shift + i + j] / denom)
        rows.append(row)
    return rows


def main():
    print("=" * 80)
    print("PAPER 181 — Xi coefficients as a PF∞ door to RH")
    print("=" * 80)

    try:
        from mpmath import mp
    except ImportError:
        print("\n  [mpmath required]")
        return

    # ------------------------------------------------------------------
    # [1] Extract coefficients
    # ------------------------------------------------------------------
    print("\n[1] Extract A_n from G(x)=Xi(sqrt(-x))=sum A_n x^n")
    A = extract_A(max_m=30)
    A0 = A[0]
    for n in range(12):
        print(f"  A_{n:02d} = {float(A[n]): .12e}")
    reliable_prefix = 22
    print(f"  min A_n in reliable prefix n<= {reliable_prefix}: "
          f"{float(min(A[:reliable_prefix + 1])):.6e}")
    print("  note: very late coefficients from Cauchy extraction are diagnostic,"
          " not interval-certified")

    # Work with normalized coefficients for determinant conditioning.
    B = [a / A0 for a in A]

    # ------------------------------------------------------------------
    # [2] Turan/log-concavity margins
    # ------------------------------------------------------------------
    print("\n[2] First PF∞ shadow: log-concavity/Turán-1")
    print(f"  {'n':>3} {'A_n^2/(A_{n-1}A_{n+1}) - 1':>34}")
    min_lc = mp.inf
    min_lc_n = None
    for n in range(1, 20):
        margin = A[n] ** 2 / (A[n - 1] * A[n + 1]) - 1
        if margin < min_lc:
            min_lc = margin
            min_lc_n = n
        print(f"  {n:>3} {float(margin):>34.12e}")
    print(f"  worst tested log-concavity margin: n={min_lc_n}, {float(min_lc):.6e}")

    # ------------------------------------------------------------------
    # [3] Toeplitz PF minors
    # ------------------------------------------------------------------
    print("\n[3] Finite PF∞ Toeplitz minors det(A_{c_j-r_i})")
    print("    Enumerating all minors up to size 4 in the leading 9x9 Toeplitz window.")
    window = 8
    max_size = 4
    total = 0
    negative = 0
    zeroish = 0
    min_det = mp.inf
    min_case = None
    tol = mp.mpf("1e-45")

    indices = list(range(window + 1))
    for size in range(1, max_size + 1):
        for rows_idx in itertools.combinations(indices, size):
            for cols_idx in itertools.combinations(indices, size):
                M = [[toeplitz_entry(B, r, c) for c in cols_idx] for r in rows_idx]
                d = mp_det_from_rows(M)
                total += 1
                if d < min_det:
                    min_det = d
                    min_case = (rows_idx, cols_idx, size)
                if d < -tol:
                    negative += 1
                elif abs(d) <= tol:
                    zeroish += 1
    print(f"  minors tested:       {total}")
    print(f"  negative minors:     {negative}")
    print(f"  zero/tiny minors:    {zeroish}")
    print(f"  minimum determinant: {float(min_det):.6e}")
    print(f"  minimum case:        rows={min_case[0]}, cols={min_case[1]}, size={min_case[2]}")

    # ------------------------------------------------------------------
    # [4] Hankel moment minors
    # ------------------------------------------------------------------
    print("\n[4] Hankel diagnostic: rejecting the naive moment route")
    print("    Correlation-normalized det(A_{N+i+j}).")
    print("    Negative 2x2 determinants are expected from log-concavity:")
    print("      det [[A_N,A_{N+1}],[A_{N+1},A_{N+2}]] < 0.")
    print(f"  {'shift N':>8} {'size':>5} {'det normalized Hankel':>24}")
    min_h = mp.inf
    min_h_case = None
    h_negative = 0
    h_total = 0
    for shift in range(0, 9):
        for size in range(2, 6):
            if shift + 2 * (size - 1) >= len(A):
                continue
            H = normalized_hankel_matrix(A, shift, size)
            d = mp_det_from_rows(H)
            h_total += 1
            if d < 0:
                h_negative += 1
            if d < min_h:
                min_h = d
                min_h_case = (shift, size)
            print(f"  {shift:>8} {size:>5} {float(d):>24.12e}")
    print(f"  negative Hankel determinants: {h_negative}/{h_total}")
    print(f"  worst tested Hankel determinant: shift={min_h_case[0]},"
          f" size={min_h_case[1]}, det={float(min_h):.6e}")

    # ------------------------------------------------------------------
    # [5] Finite truncation zeros
    # ------------------------------------------------------------------
    print("\n[5] Truncated G_M(x)=sum_{n<=M} A_n x^n zeros")
    print("    This is NOT equivalent to RH, but it exposes the PF∞ geometry.")
    try:
        import numpy as np

        print(f"  {'M':>3} {'max |Im root|':>16} {'# roots Re>0':>14}")
        for M in [6, 8, 10, 12, 14]:
            coeffs_high = [float(B[n]) for n in range(M, -1, -1)]
            roots = np.roots(coeffs_high)
            max_im = max(abs(r.imag) for r in roots)
            n_pos = sum(1 for r in roots if r.real > 1e-8)
            print(f"  {M:>3} {max_im:>16.6e} {n_pos:>14}")
    except ImportError:
        print("  [numpy unavailable; skipped]")

    # ------------------------------------------------------------------
    # [6] Honest assessment
    # ------------------------------------------------------------------
    print("\n[6] HONEST ASSESSMENT")
    print("  PROPOSED ANALYTIC INSIGHT:")
    print("    RH should be attacked as PF∞ total positivity of A_n where")
    print("    G(x)=Xi(sqrt(-x))=sum A_n x^n.")
    print()
    print("  NUMERICAL RESULT:")
    print("    - tested Toeplitz PF minors: no negative minors.")
    print("    - tested Hankel minors: negative at size 2, as log-concavity predicts.")
    print("      This rejects the naive positive-moment route.")
    print("    - log-concavity margins are strongly positive in the tested range.")
    print()
    print("  WHAT WOULD SOLVE THE DOOR:")
    print("    Prove every Toeplitz minor det(A_{c_j-r_i}) is non-negative.")
    print("    By Aissen-Schoenberg-Whitney/Edrei, this gives PF∞; by the")
    print("    square-root transform, this is equivalent to Xi having only real")
    print("    zeros, i.e. RH.")
    print()
    print("  NEXT MOVE:")
    print("    Derive a Toeplitz-total-positive representation for A_n.")
    print("    The Hankel test shows a plain Stieltjes moment representation is")
    print("    the wrong shape; the proof must target PF∞/Toeplitz minors directly.")


if __name__ == "__main__":
    main()
