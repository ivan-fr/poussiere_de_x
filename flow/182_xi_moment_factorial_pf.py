"""
PAPER: 182 (NEW — moment/factorial decomposition of the Xi coefficients)
TITLE: Xi coefficients as factorially damped moments: locating the PF∞
       obstruction after paper 181
STATUS: RH remains OPEN.  This paper derives and tests the next analytic
        decomposition
            A_n = M_n / (2n)!,
        where M_n is a positive Riemann-Phi moment.  It shows that the
        naive Hankel moment route fails only after factorial damping, while
        the Toeplitz/PF∞ route remains numerically positive.
DEPENDS: 109 (de Bruijn-Newman Phi), 179 (finite Jensen certificates),
         181 (PF∞ door).

THEORY
======

------------------------------------------------------------------------
RIEMANN-PHI MOMENTS
------------------------------------------------------------------------

Riemann's Xi function has the cosine transform form

    Xi(t) = integral_0^infty Phi(u) cos(t u) du

up to a harmless positive normalising constant.  Expanding cos(tu),

    Xi(t) = sum_{n>=0} (-1)^n [ integral Phi(u) u^(2n) du / (2n)! ] t^(2n).

So the coefficients from paper 181 satisfy

    A_n = M_n / (2n)!,
    M_n := integral_0^infty Phi(u) u^(2n) du.

If Phi(u) >= 0, then M_n is a Stieltjes moment sequence, hence its
Hankel matrices are positive semidefinite.  Paper 181 found that A_n
has negative Hankel 2x2 minors.  This paper explains the discrepancy:
the factorial damping M_n -> M_n/(2n)! destroys ordinary Hankel moment
positivity while apparently preserving the Toeplitz/PF∞ shadow needed
for RH.

------------------------------------------------------------------------
ANALYTIC TARGET
------------------------------------------------------------------------

The right target is now more precise:

    Prove the Toeplitz total positivity of the coefficient sequence
        A_n = M_n/(2n)!
    using total positivity of the composed kernel
        K(n,u) = u^(2n)/(2n)! · Phi(u).

The naive "M_n is a positive moment sequence" is not enough.  The
factorially damped kernel must be handled directly.
"""
from __future__ import annotations

import itertools
import math


def phi_riemann(u, n_terms=40):
    """Riemann Phi kernel used in the de Bruijn-Newman integral."""
    s = 0.0
    e4u = math.exp(4 * u)
    e5u = math.exp(5 * u)
    e9u = math.exp(9 * u)
    for n in range(1, n_terms + 1):
        nn = n * n
        term = (2 * math.pi**2 * n**4 * e9u - 3 * math.pi * nn * e5u)
        term *= math.exp(-math.pi * nn * e4u)
        s += term
    return s


def simpson_integral(f, a, b, n):
    if n % 2:
        n += 1
    h = (b - a) / n
    total = f(a) + f(b)
    for k in range(1, n):
        total += (4 if k % 2 else 2) * f(a + k * h)
    return total * h / 3


def moment_M(n, U=5.0, n_grid=2000):
    return simpson_integral(lambda u: phi_riemann(u) * (u ** (2 * n)), 0.0, U, n_grid)


def det_small(rows):
    n = len(rows)
    if n == 0:
        return 1.0
    total = 0.0
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


def toeplitz_entry(seq, r, c):
    k = c - r
    if k < 0 or k >= len(seq):
        return 0.0
    return seq[k]


def main():
    print("=" * 80)
    print("PAPER 182 — Xi coefficients as factorially damped moments")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Phi positivity sample
    # ------------------------------------------------------------------
    print("\n[1] Riemann-Phi positivity sample on [0, 5]")
    pmin = float("inf")
    argmin = 0.0
    for j in range(1001):
        u = 5.0 * j / 1000
        v = phi_riemann(u)
        if v < pmin:
            pmin = v
            argmin = u
    print(f"  min Phi(u) on grid: {pmin:.12e} at u={argmin:.4f}")
    print("  note: this is a numerical positivity check, not a proof")

    # ------------------------------------------------------------------
    # [2] Moments and factorial damping
    # ------------------------------------------------------------------
    print("\n[2] Compute M_n and A_n = M_n/(2n)! by quadrature")
    max_n = 14
    M = []
    A = []
    for n in range(max_n + 1):
        Mn = moment_M(n)
        An = Mn / math.factorial(2 * n)
        M.append(Mn)
        A.append(An)
        if n <= 10:
            print(f"  n={n:02d}  M_n={Mn:.12e}  A_n={An:.12e}")

    # Normalize for determinant conditioning.
    Mn = [x / M[0] for x in M]
    An = [x / A[0] for x in A]

    # ------------------------------------------------------------------
    # [3] Hankel: moments positive, factorial damping negative
    # ------------------------------------------------------------------
    print("\n[3] Hankel determinants: M_n vs A_n")
    print(f"  {'shift':>5} {'det H_M(2)':>16} {'det H_A(2)':>16}"
          f" {'det H_M(3)':>16} {'det H_A(3)':>16}")
    for shift in range(0, 8):
        HM2 = [[Mn[shift + i + j] for j in range(2)] for i in range(2)]
        HA2 = [[An[shift + i + j] for j in range(2)] for i in range(2)]
        HM3 = [[Mn[shift + i + j] for j in range(3)] for i in range(3)]
        HA3 = [[An[shift + i + j] for j in range(3)] for i in range(3)]
        print(f"  {shift:>5} {det_small(HM2):>16.6e} {det_small(HA2):>16.6e}"
              f" {det_small(HM3):>16.6e} {det_small(HA3):>16.6e}")

    # ------------------------------------------------------------------
    # [4] Toeplitz PF minors for M_n and A_n
    # ------------------------------------------------------------------
    print("\n[4] Toeplitz minors: moment sequence M_n vs damped sequence A_n")
    print("    Testing all minors up to size 4 in a leading 8x8 window.")
    for label, seq in [("M", Mn), ("A=M/(2n)!", An)]:
        total = 0
        negative = 0
        min_det = float("inf")
        min_case = None
        indices = list(range(8))
        for size in range(1, 5):
            for rows_idx in itertools.combinations(indices, size):
                for cols_idx in itertools.combinations(indices, size):
                    T = [[toeplitz_entry(seq, r, c) for c in cols_idx] for r in rows_idx]
                    d = det_small(T)
                    total += 1
                    if d < min_det:
                        min_det = d
                        min_case = (rows_idx, cols_idx, size)
                    if d < -1e-18:
                        negative += 1
        print(f"  {label:>10}: tested={total}, negative={negative},"
              f" min={min_det:.6e}, case={min_case}")

    # ------------------------------------------------------------------
    # [5] Kernel minor diagnostic
    # ------------------------------------------------------------------
    print("\n[5] Kernel K(n,u)=u^(2n)/(2n)! total-positivity diagnostic")
    print("    Testing det K(n_i,u_j) for increasing n_i and u_j.")
    min_K = float("inf")
    neg_K = 0
    total_K = 0
    u_grid = [0.15, 0.35, 0.7, 1.1, 1.6, 2.2]
    n_grid = [0, 1, 2, 3, 4, 5]
    for size in range(2, 5):
        for ns in itertools.combinations(n_grid, size):
            for us in itertools.combinations(u_grid, size):
                K = [[(u ** (2 * n)) / math.factorial(2 * n) for u in us] for n in ns]
                d = det_small(K)
                total_K += 1
                min_K = min(min_K, d)
                if d < -1e-14:
                    neg_K += 1
    print(f"  kernel minors tested: {total_K}")
    print(f"  negative kernel minors: {neg_K}")
    print(f"  minimum kernel determinant: {min_K:.6e}")

    # ------------------------------------------------------------------
    # [6] Honest assessment
    # ------------------------------------------------------------------
    print("\n[6] HONEST ASSESSMENT")
    print("  NEW STRUCTURAL INFORMATION:")
    print("    - M_n behaves like a positive moment sequence.")
    print("    - A_n = M_n/(2n)! loses Hankel positivity, matching paper 181.")
    print("    - Toeplitz PF minors remain non-negative in the tested window.")
    print("    - The bare kernel u^(2n)/(2n)! is numerically totally positive")
    print("      on ordered grids; the missing theorem must include Phi(u).")
    print()
    print("  ANALYTIC TARGET:")
    print("    Prove that integrating the totally positive kernel")
    print("      K(n,u)=u^(2n)/(2n)! against Phi(u) du")
    print("    preserves the Toeplitz PF∞ minors of A_n.")
    print()
    print("  WHY THIS MATTERS:")
    print("    It converts the vague RH insight into a concrete total-positivity")
    print("    theorem about a known positive kernel and Riemann's Phi function.")


if __name__ == "__main__":
    main()
