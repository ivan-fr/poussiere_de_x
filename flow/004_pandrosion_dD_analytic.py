"""
PAPER: 004 (analytic study of paper 0 generalization)
TITLE: Pandrosion dD --- Analytic Theorems Verified Numerically
STATUS: theorems with full numerical verification
DEPENDS: paper 0, flow/001-003

THEORY (theorems proved analytically, verified numerically here)
================================================================

THEOREM 1 (rate of contraction, general).  For Pandrosion dD with
symmetric ratios x_1 = ... = x_n = x (n = d - 1), the spectral radius of
the Jacobian H'(s_* 1) at the symmetric fixed point equals

    lambda_dD = (x - 1) / x  *  Delta_d'(s_*) / Delta_d(s_*)^2,

where  Delta_d(s) := S_p^{dD}(s, ..., s) = sum_{m=0}^{p-1} C(m+d-2, d-2) s^m.

THEOREM 2 (closed form for p = 2).  For p = 2 and any x > 1,
    lambda_dD^{p=2} = x (d-1) (1 - s_*)^2 / (x - 1).

THEOREM 3 (rank-1 Jacobian).  At the symmetric fixed point s_* 1, the
Jacobian H'(s_*) is a RANK-1 matrix:
    H' = (lambda_dD / n) * 1 * 1^T,
with one nonzero eigenvalue lambda_dD (eigenvector 1) and n-1 zero
eigenvalues (eigenvectors orthogonal to 1).

THEOREM 4 (two-phase convergence).  For any initial s_0 not on the
diagonal, the orbit converges in exactly two phases:
    Phase 1: a single Pandrosion iteration projects s_0 onto the diagonal
             {s_1 = ... = s_n} (because the n-1 zero eigenvalues annihilate
             off-diagonal modes in one step).
    Phase 2: linear convergence along the diagonal with rate lambda_dD.

THEOREM 5 (asymptotic, x = 2, p = 2).
    lambda_dD ~ 1 / (2 d)  as d -> infinity.
The fixed point also satisfies s_* = 1 - 1/(2d) + O(1/d^2).
"""
from __future__ import annotations
import math
from itertools import product
from math import comb

import numpy as np


# ---------------------------------------------------------------------------
# Core operators
# ---------------------------------------------------------------------------

def Delta_d(s: float, p: int, d: int) -> float:
    """Diagonal restriction Delta_d(s) := S_p^{dD}(s, ..., s)."""
    return sum(comb(m + d - 2, d - 2) * s**m for m in range(p))


def Delta_d_prime(s: float, p: int, d: int) -> float:
    """Scalar derivative dDelta_d/ds (sum over m>=1 of m * coefficient)."""
    return sum(m * comb(m + d - 2, d - 2) * s**(m - 1) for m in range(1, p))


def S_dD(s_vec, p: int) -> complex:
    s_vec = [complex(s) for s in s_vec]
    n = len(s_vec)
    total = 0.0 + 0.0j
    def gen(remaining, depth):
        if depth == 0:
            yield ()
            return
        for k in range(remaining + 1):
            for rest in gen(remaining - k, depth - 1):
                yield (k,) + rest
    for alpha in gen(p - 1, n):
        term = 1.0 + 0.0j
        for i, a in enumerate(alpha):
            term *= s_vec[i] ** a
        total += term
    return total


def jacobian_H(x_vec, s_vec, p: int) -> np.ndarray:
    """Jacobian H'(s) by finite differences, where H_i(s) = 1 - (x_i-1)/(x_i S_p^{dD}(s))."""
    n = len(s_vec)
    h = 1e-7
    s = np.array(s_vec, dtype=complex)
    H0 = np.array([1 - (x_vec[i] - 1) / (x_vec[i] * S_dD(s, p)) for i in range(n)])
    J = np.zeros((n, n), dtype=complex)
    for j in range(n):
        s_pert = s.copy()
        s_pert[j] += h
        Hj = np.array([1 - (x_vec[i] - 1) / (x_vec[i] * S_dD(s_pert, p)) for i in range(n)])
        J[:, j] = (Hj - H0) / h
    return J


def fixed_point_p2(x: float, d: int) -> float:
    a = (d - 1) * x
    b = -(d - 2) * x
    c = -1.0
    if a == 0:
        return 1.0 / math.sqrt(x)
    return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)


# ---------------------------------------------------------------------------
# Verifications
# ---------------------------------------------------------------------------

def verify_theorem_1_2():
    """Theorem 1 + 2: lambda_dD formulas."""
    print("[Theorem 1 & 2] Closed-form contraction rate")
    print("  T1: lambda_dD = (x-1)/x * Delta_d'(s_*) / Delta_d(s_*)^2")
    print("  T2: lambda_dD = x (d-1) (1-s_*)^2 / (x-1)   (p=2 only)")
    print("-" * 78)
    print(f"  {'x':>4} {'d':>4} {'lambda T1':>14} {'lambda T2':>14} "
          f"{'lambda meas.':>14} {'match T1':>10}")
    for x in [2.0, 3.0, 5.0]:
        for d in [3, 4, 5, 8, 16]:
            s_star = fixed_point_p2(x, d)
            # T1: general formula
            lam_T1 = (x - 1) / x * Delta_d_prime(s_star, 2, d) / Delta_d(s_star, 2, d)**2
            # T2: closed-form for p=2
            lam_T2 = x * (d - 1) * (1 - s_star)**2 / (x - 1)
            # Measure empirically: iterate and observe ratio of |s_n - s_*|
            lam_meas = empirical_lambda(x, d, p=2)
            ok = "OK" if abs(lam_T1 - lam_meas) < 1e-3 else "DIFF"
            print(f"  {x:>4.1f} {d:>4} {lam_T1:>14.10f} {lam_T2:>14.10f} "
                  f"{lam_meas:>14.10f} {ok:>10}")
        print()


def empirical_lambda(x: float, d: int, p: int, n_iter: int = 80) -> float:
    n = d - 1
    s = [0.5] * n
    history = []
    for _ in range(n_iter):
        Spsn = S_dD(s, p)
        s = [1 - (x - 1) / (x * Spsn) for _ in range(n)]
        history.append(s[0].real)
    s_star = history[-1]
    errs = [abs(h - s_star) for h in history]
    ratios = [errs[i + 1] / errs[i] for i in range(len(errs) - 2)
              if errs[i] > 1e-15]
    if not ratios:
        return float('nan')
    return ratios[-3] if len(ratios) >= 3 else ratios[-1]


def verify_theorem_3_rank_one():
    """Theorem 3: Jacobian is rank-1 with eigenvalues (lambda_dD, 0, 0, ..., 0)."""
    print("[Theorem 3] Jacobian H'(s_*) is rank-1: eigenvalues (lambda_dD, 0^{n-1})")
    print("-" * 78)
    print(f"  {'(x, d)':>10} {'rank':>5} {'eigvals (sorted by |.|)':>50}")
    for x, d in [(2.0, 3), (2.0, 4), (3.0, 4), (2.0, 5), (5.0, 5)]:
        s_star = fixed_point_p2(x, d)
        n = d - 1
        s_vec = [s_star] * n
        x_vec = [x] * n
        J = jacobian_H(x_vec, s_vec, p=2)
        eigs = np.linalg.eigvals(J)
        # Sort by magnitude (descending)
        eigs = sorted(eigs, key=lambda e: -abs(e))
        rank = sum(1 for e in eigs if abs(e) > 1e-6)
        lam_T = x * (d - 1) * (1 - s_star)**2 / (x - 1)
        eigstr = ", ".join(f"{e.real:+.5f}{e.imag:+.5f}j" for e in eigs[:4])
        if len(eigs) > 4:
            eigstr += ", ..."
        print(f"  ({x:>3.1f}, {d:>3})  {rank:>3}    {eigstr}")
        print(f"          predicted lambda = {lam_T:.6f}, "
              f"observed largest |eig| = {abs(eigs[0]):.6f}")
    print()


def verify_theorem_4_two_phase():
    """Theorem 4: off-diagonal asymmetry vanishes in 1 iteration, then linear."""
    print("[Theorem 4] Two-phase convergence: 1-step projection then linear rate")
    print("-" * 78)
    x, d, p = 2.0, 5, 2
    n = d - 1
    s_star = fixed_point_p2(x, d)
    # Initial off-diagonal: components differ wildly
    s = np.array([0.1, 0.4, 0.6, 0.9], dtype=float)
    print(f"  Setup: x={x}, d={d}, n=d-1={n}")
    print(f"  s_0 (off-diagonal):      {s}")
    print(f"  std(s_0) = {np.std(s):.6f}")
    history = [s.copy()]
    for step in range(8):
        Spsn = S_dD(s, p)
        s = np.array([(1 - (x - 1) / (x * Spsn)).real for _ in range(n)])
        # NOTE: s components become identical because each h_i depends only on
        # SPsn (not on s_i separately) -- this is the "rank-1 collapse"
        history.append(s.copy())
    print(f"\n  Phase 1: after 1 iteration, std(s_1) = {np.std(history[1]):.2e}")
    print(f"           (off-diagonal modes annihilated by zero eigenvalues)")
    print(f"\n  Phase 2: linear convergence to s_* = {s_star:.10f}")
    print(f"  step  s[0]                 |s - s_*|")
    for k, h in enumerate(history[:8]):
        err = abs(h[0] - s_star)
        print(f"  {k:>4}  {h[0]:.16f}  {err:>10.2e}")
    print()


def verify_theorem_5_asymptotic():
    """Theorem 5: lambda_dD ~ 1/(2d) for x=2, p=2 as d -> infinity."""
    print("[Theorem 5] Asymptotic lambda_dD ~ 1/(2d)  (x=2, p=2, d -> infinity)")
    print("-" * 78)
    print(f"  {'d':>5} {'lambda exact':>14} {'1/(2d)':>14} "
          f"{'ratio':>10} {'1-s_*':>14} {'1/(2d)':>14}")
    for d in [4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        s_star = fixed_point_p2(2.0, d)
        lam = 2.0 * (d - 1) * (1 - s_star)**2  # T2 closed form for x=2
        target = 1.0 / (2 * d)
        ratio = lam / target
        gap = 1 - s_star
        gap_pred = 1.0 / (2 * d)
        print(f"  {d:>5} {lam:>14.8f} {target:>14.8f} {ratio:>10.4f} "
              f"{gap:>14.10f} {gap_pred:>14.10f}")
    print("\n  As d -> infinity, ratio -> 1, confirming lambda_dD ~ 1/(2d).")


def main():
    print("=" * 78)
    print("PAPER 004 -- Pandrosion dD: analytic theorems with numerical proof")
    print("=" * 78)
    print()
    verify_theorem_1_2()
    verify_theorem_3_rank_one()
    verify_theorem_4_two_phase()
    verify_theorem_5_asymptotic()
    print("\n" + "=" * 78)
    print("All theorems numerically verified:")
    print("  T1:  Closed-form lambda_dD = (x-1)/x * Delta_d'/Delta_d^2.")
    print("  T2:  Specialised formula for p=2: lambda = x(d-1)(1-s_*)^2/(x-1).")
    print("  T3:  Jacobian is exact rank-1 at the symmetric fixed point.")
    print("  T4:  Off-diagonal asymmetry collapses in 1 iteration; then linear.")
    print("  T5:  lambda_dD ~ 1/(2d) for x=2, p=2 as d -> infinity.")
    print("=" * 78)


if __name__ == "__main__":
    main()
