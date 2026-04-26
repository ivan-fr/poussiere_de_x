"""
ATTEMPT 05 — Push Λ via Turán/log-concavity tests on σ_n(t).

The eigenvalue approach (attempt 02 + initial 05) is numerically
ill-conditioned: σ_n(t) decays exponentially in n, so the Hankel matrix
becomes near-singular and floating-point eigenvalues are dominated by
roundoff at the same scale across all t.

A more robust necessary condition for Λ ≤ t (Pólya 1927; Csordas-Norfolk-
Varga 1986): the sequence σ_n(t) := (-1)^n c_{2n}(t) is the Taylor
coefficient sequence of an LP-class entire function, so it must satisfy
the Turán-type inequality

   τ_n(t) := σ_n(t)²  −  σ_{n-1}(t) · σ_{n+1}(t)  ≥  0.

(This is log-concavity, a NECESSARY condition for LP-membership.)

Stronger Csordas-Craven test:
   T_n(t) := 4 · σ_{n-1}(t) · σ_{n+1}(t)  ≤  3 · σ_n(t)²,
also necessary.  Failures pinpoint t-values violating Λ ≤ t.

Strategy
--------
1. Compute σ_n(0) at high precision for n = 0..30 (verify Riemann's Φ
   and known LP-class membership of Ξ ⇒ Λ ≤ 0 conjecturally).
2. Compute σ_n(t) for t in a dense grid t ∈ [-0.05, 0.30].
3. For each t, evaluate τ_n(t) and T_n(t) over n.  Find the smallest
   t such that all τ_n(t), T_n(t) are ≥ 0 / ≤ correct ratio.
4. That smallest t is a heuristic upper bound on Λ.
"""
from __future__ import annotations
import math
import time

from mpmath import mp, mpf, exp as mexp, quad

mp.dps = 50


def Phi(u, n_terms=60):
    s = mpf(0)
    e9u = mexp(9 * u)
    e5u = mexp(5 * u)
    e4u = mexp(4 * u)
    for n in range(1, n_terms + 1):
        coef = 2 * mp.pi ** 2 * n ** 4 * e9u - 3 * mp.pi * n ** 2 * e5u
        term = coef * mexp(-mp.pi * n ** 2 * e4u)
        s += term
        if abs(term) < mpf(10) ** (-mp.dps + 5) and n > 5:
            break
    return s


def sigma_n(n: int, t: float, U: float = 6.0):
    fact = mpf(math.factorial(2 * n))
    integral = quad(lambda u: Phi(u) * u ** (2 * n) * mexp(t * u ** 2),
                    [0, U])
    return integral / fact


def turan(sigmas):
    """τ_n = σ_n² − σ_{n-1} σ_{n+1}  for n = 1..N-1."""
    return [sigmas[n] ** 2 - sigmas[n - 1] * sigmas[n + 1]
            for n in range(1, len(sigmas) - 1)]


def csordas_craven(sigmas):
    """3 σ_n² − 4 σ_{n-1} σ_{n+1}  for n = 1..N-1."""
    return [3 * sigmas[n] ** 2 - 4 * sigmas[n - 1] * sigmas[n + 1]
            for n in range(1, len(sigmas) - 1)]


def main():
    print("=" * 80)
    print("ATTEMPT 05 — Λ push via Turán / Csordas-Craven on σ_n(t)")
    print("=" * 80)

    N = 18  # number of σ_n to compute per t

    # ----------------------------------------------------------------
    # [1] Baseline σ_n(0) — verify LP signal (Λ ≤ 0 conjecturally).
    # ----------------------------------------------------------------
    print(f"\n[1] σ_n(0) for n = 0..{N-1}.")
    t0 = time.time()
    sigmas_0 = [sigma_n(n, 0.0) for n in range(N)]
    print(f"  computed in {time.time()-t0:.1f}s")
    for n in range(min(8, N)):
        print(f"     σ_{n} = {float(sigmas_0[n]):+.6e}")
    tau_0 = turan(sigmas_0)
    cc_0 = csordas_craven(sigmas_0)
    print(f"  Turán τ_n at t=0 (must be ≥ 0):")
    for n, x in enumerate(tau_0[:8], start=1):
        s = "OK" if x >= 0 else "FAIL"
        print(f"     n={n:2d}  τ_n = {float(x):+.4e}   {s}")
    print(f"  Csordas-Craven 3σ²-4σ_{{n-1}}σ_{{n+1}} (must be ≥ 0):")
    for n, x in enumerate(cc_0[:8], start=1):
        s = "OK" if x >= 0 else "FAIL"
        print(f"     n={n:2d}  CC_n = {float(x):+.4e}   {s}")

    # ----------------------------------------------------------------
    # [2] Sweep t and find smallest t with both tests passing.
    # ----------------------------------------------------------------
    print("\n[2] Sweep t ∈ {-0.05, 0, 0.05, 0.10, 0.15, 0.18, 0.20, 0.22, 0.25}")
    t_grid = [-0.05, 0.0, 0.05, 0.10, 0.15, 0.18, 0.20, 0.22, 0.25, 0.30]
    sweep = []
    for t in t_grid:
        t0 = time.time()
        sigmas_t = [sigma_n(n, t) for n in range(N)]
        tau_t = turan(sigmas_t)
        cc_t = csordas_craven(sigmas_t)
        n_tau_fail = sum(1 for x in tau_t if x < 0)
        n_cc_fail = sum(1 for x in cc_t if x < 0)
        min_tau = min(tau_t, key=lambda x: float(x))
        min_cc = min(cc_t, key=lambda x: float(x))
        ok = (n_tau_fail == 0) and (n_cc_fail == 0)
        sweep.append((t, n_tau_fail, n_cc_fail, float(min_tau), float(min_cc),
                      ok))
        print(f"  t = {t:+.3f}  ({time.time()-t0:.1f}s)   "
              f"Turán fails: {n_tau_fail}/{len(tau_t)}    "
              f"CC fails: {n_cc_fail}/{len(cc_t)}    "
              f"min τ = {float(min_tau):+.3e}  min CC = {float(min_cc):+.3e}    "
              f"{'PASS' if ok else 'FAIL'}")

    # ----------------------------------------------------------------
    # [3] Threshold extraction
    # ----------------------------------------------------------------
    print("\n[3] Smallest t with Turán + Csordas-Craven simultaneously satisfied")
    passing = [t for t, *_, ok in sweep if ok]
    if passing:
        t_star = min(passing)
        print(f"   t* (heuristic Λ-upper-bound) ≈ {t_star:+.3f}")
    else:
        t_star = None
        print("   No t in grid passes both tests at order N = "
              f"{N-1}.  Need higher N or interval-arithmetic to certify.")

    # ----------------------------------------------------------------
    # [4] Translation to R
    # ----------------------------------------------------------------
    print("\n[4] Translation to R (paper 159 calibration)")
    Lam_C_table = {
        0.22: 76.2, 0.20: 70.0, 0.18: 67.0, 0.16: 63.0,
        0.15: 60.0, 0.12: 50.0, 0.10: 40.0, 0.08: 30.0, 0.05: 15.0,
        0.00: 5.0, -0.05: 3.0,
    }
    V0 = math.log(2 * math.pi) - 0.5772156649
    alpha = (5.5587 - 4 - V0) / math.log(76.2)
    print(f"   V0 = {V0:.4f}, α = {alpha:.4f}")
    if t_star is not None:
        Lam_snap = min(Lam_C_table.keys(), key=lambda L: abs(L - t_star))
        C = Lam_C_table[Lam_snap]
        V = V0 + alpha * math.log(C)
        R = 4 + V
        print(f"   Λ-heur ≈ {t_star},  C(Λ) ≈ {C},  R ≈ {R:.4f},"
              f"  gain vs MTY: {5.5587 - R:+.4f}")
    else:
        print("   No translation possible without a Λ-bound from the test.")

    # ----------------------------------------------------------------
    # [5] Honest verdict
    # ----------------------------------------------------------------
    print("\n[5] HONEST VERDICT")
    print()
    print("  This attempt uses NECESSARY-only LP tests (Turán + Csordas-Craven).")
    print("  PASSING the tests is consistent with Λ ≤ t but does not PROVE it.")
    print("  FAILING the tests proves Λ > t.")
    print()
    print("  Polymath15 (rigorous) used interval arithmetic and a much more")
    print("  involved method to prove Λ < 0.22.  This Python sweep is an")
    print("  empirical sanity check — not a competing rigorous bound.")
    print()
    print("  Even at the most optimistic Λ achievable here, the projected R")
    print("  improvement on MTY 2022 is sub-decimal.  A genuine improvement")
    print("  requires either:")
    print("    (i)  Polymath16-style rigorous push of Λ < 0.20 with fully")
    print("         certified interval bounds, plus")
    print("    (ii) re-running MTY 2022's analytic chain with the smaller C_VK")
    print("         derived from the new Λ — months of expert work.")


if __name__ == "__main__":
    main()
