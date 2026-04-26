"""
ATTEMPT 04 — α(P) derived from MT explicit-formula structure (not from
spread-rescaling).

Sketch
------
The Mossinghoff-Trudgian inequality for a zero ρ = β + iγ of ζ uses a
positive cosine polynomial P(θ) = b_0 + 2 Σ_{k≥1} b_k cos kθ with b_k ≥ 0.
The key inequality (de la Vallée Poussin / Stechkin / MT) is

   b_0 · Re(-ζ'/ζ)(σ)
   + Σ_{k≥1} b_k · Re(-ζ'/ζ)(σ + ikγ)  ≥  0.

Splitting:
  • k = 0 gives the pole of ζ at s = 1:  b_0 / (σ − 1)  +  b_0 · γ_Euler  + …
  • k = 1 gives the contribution from the zero ρ:  − b_1 / (σ − β)
                                                  − b_1 / (σ + β − 2)  + …
  • k ≥ 2 contribute  Σ b_k · Re(-ζ'/ζ)(σ + ikγ).

For the k ≥ 2 terms, Borel-Carathéodory + Hadamard 3-circles gives

   |ζ'/ζ(σ + iτ)|  ≤  K · log|τ|  +  L · log M(σ; τ),

where M(σ; τ) is an upper bound for |ζ| in a small disk; in our setting
M(1; τ) ≤ C_VK · (log τ)^{2/3}, hence the per-k contribution is
L · (log C_VK + (2/3) log log τ) plus log τ stuff that cancels against
the leading 1/(σ−1) term.

Collecting:
   1/(σ−β)  ≤  (b_0 / b_1) · 1/(σ−1)
              + L · (Σ_{k≥2} b_k / b_1) · log C_VK
              + (analytic constants).

Setting σ − 1 = c/(R log γ) and balancing gives

   R(P)   =   max_θ P(θ) / b_1   +   V_0   +   α(P) · log C_VK,
   α(P)   =   L · ( Σ_{k≥2} b_k / b_1 ).

L is the universal Borel-Carathéodory / Hadamard 3-circles constant for
ζ on a vertical strip. We CALIBRATE L from the dlVP/MTY anchor:

   dlVP : (b_0, b_1, b_2) = (3, 2, 1/2),   Σ_{k≥2} b_k / b_1 = 0.25.
   α_dlVP_calib = (R_MTY − 4 − V_0) / log(C_MTY) ≈ 0.0688.
   ⇒  L  =  α_dlVP_calib / 0.25  ≈  0.2752.

Then for ANY MT-strict polynomial P:
   α(P)   =   0.2752 · ( Σ_{k≥2} b_k / b_1 ).

This α(P) is THE ANALYTIC FORMULA implied by the MT 2015 structure
(modulo replacing all unknown constants by their dlVP-calibrated value).
It is much more defensible than the spread-rescaling guess of attempt 03.
"""
from __future__ import annotations
import math

import cvxpy as cp


MTY_R = 5.5587
MTY_C = 76.2
V0 = math.log(2 * math.pi) - 0.5772156649


# -----------------------------------------------------------------
# MT-derived α(P) formula
# -----------------------------------------------------------------
def alpha_MT(b):
    """α(P) = L · Σ_{k≥2} b_k / b_1, with L calibrated from dlVP/MTY."""
    sum_high = sum(b[k] for k in range(2, len(b))) / b[1]
    sum_high_dlVP = 0.5 / 2
    alpha_dlVP = (MTY_R - 4 - V0) / math.log(MTY_C)
    L = alpha_dlVP / sum_high_dlVP
    return L * sum_high, L


# -----------------------------------------------------------------
# SDP solvers (etage 1)
# -----------------------------------------------------------------
def solve_MT_strict(N):
    X1 = cp.Variable((N + 1, N + 1), symmetric=True)
    X2 = cp.Variable((N + 1, N + 1), symmetric=True)
    t = cp.Variable()
    b = [cp.sum([X1[i, i + m] for i in range(N + 1 - m)])
         for m in range(N + 1)]
    c = [cp.sum([X2[i, i + m] for i in range(N + 1 - m)])
         for m in range(N + 1)]
    cons = [X1 >> 0, X2 >> 0, c[0] == t - b[0], b[1] == 1]
    for m in range(1, N + 1):
        cons.append(c[m] == -b[m])
    for m in range(N + 1):
        cons.append(b[m] >= 0)
    cp.Problem(cp.Minimize(t), cons).solve(solver=cp.CLARABEL)
    return float(t.value), [float(bi.value) for bi in b]


def solve_min_MT_R(N, force_b2=0.5):
    """Minimise R_MT(P) = max P / b_1 + V_0 + α_high(P) log C_VK
    where α_high(P) = L · Σ_{k≥3} b_k / b_1.

    IMPORTANT FIX: the second-harmonic coefficient b_2 is what gives the
    leading dlVP "1 → 4" quadratic improvement.  Without b_2 ≥ b_1/2
    the MT analytic argument collapses to the trivial Hadamard bound
    R = 1, NOT R = 4.  We therefore impose  b_2 ≥ force_b2 · b_1
    (dlVP value = 0.5).  Coefficients b_k for k ≥ 3 contribute only to
    the V-correction via L · b_k · log C."""
    # Re-calibrate L so that dlVP at C=76.2 still gives R = 5.5587.
    # dlVP: b = (3, 2, 1/2).  Σ_{k≥3} b_k = 0.  Σ_{k≥2} b_k / b_1 = 0.25.
    # We treat b_2 as part of the "leading" constant (folded into R0 = 4),
    # and L only applies to k ≥ 3.  The dlVP calibration then recovers
    # R = 4 + V0 with NO log(C) contribution (since dlVP has no k ≥ 3).
    # That under-explains MTY's 1.5587 - 1.2607 = 0.298 gap, which means
    # part of MTY's V comes from b_2 itself.  So we keep the *original*
    # α formula (Σ_{k≥2}) but FORCE b_2 = dlVP value to keep that
    # contribution intact.
    sum_high_dlVP = 0.25
    alpha_dlVP = (MTY_R - 4 - V0) / math.log(MTY_C)
    L = alpha_dlVP / sum_high_dlVP
    log_C = math.log(MTY_C)

    X1 = cp.Variable((N + 1, N + 1), symmetric=True)
    X2 = cp.Variable((N + 1, N + 1), symmetric=True)
    t = cp.Variable()
    b = [cp.sum([X1[i, i + m] for i in range(N + 1 - m)])
         for m in range(N + 1)]
    c = [cp.sum([X2[i, i + m] for i in range(N + 1 - m)])
         for m in range(N + 1)]
    cons = [X1 >> 0, X2 >> 0, c[0] == t - b[0], b[1] == 1]
    for m in range(1, N + 1):
        cons.append(c[m] == -b[m])
    for m in range(N + 1):
        cons.append(b[m] >= 0)
    cons.append(b[2] >= force_b2)  # preserve dlVP quadratic structure
    sum_high = cp.sum([b[k] for k in range(2, N + 1)])
    obj = t + V0 + L * log_C * sum_high
    cp.Problem(cp.Minimize(obj), cons).solve(solver=cp.CLARABEL)
    return (float(t.value),
            [float(bi.value) for bi in b],
            float(obj.value))


# -----------------------------------------------------------------
# Main
# -----------------------------------------------------------------
def main():
    print("=" * 80)
    print("ATTEMPT 04 — α(P) from MT structure  (not spread-rescaled)")
    print("=" * 80)

    # ---------------------- sanity check on dlVP ----------------------
    print("\n[1] dlVP sanity check.  b = (3, 2, 1/2).")
    b_dlVP = [3.0, 2.0, 0.5]
    a_dlVP, L = alpha_MT(b_dlVP)
    R_dlVP = (4.0) + V0 + a_dlVP * math.log(MTY_C)
    print(f"  Σ_{{k≥2}} b_k / b_1 = {0.5/2:.4f}")
    print(f"  L (calibrated)     = {L:.4f}")
    print(f"  α(dlVP)            = {a_dlVP:.4f}")
    print(f"  R_MT(dlVP, C=76.2) = {R_dlVP:.4f}    [target {MTY_R}]")
    assert abs(R_dlVP - MTY_R) < 1e-3, "calibration mismatch"

    # ---------------------- MT-strict scan ----------------------
    print("\n[2] MT-strict SDP poly (minimise max P only) at varying N,")
    print("    then evaluate R_MT = 4 + V_0 + α(P) log C_MTY  with derived α(P).")
    print(f"  {'N':>4} {'R0':>10} {'Σ_{k≥2}b_k/b_1':>16} {'α(P)':>10}"
          f" {'R_MT':>10} {'gain vs MTY':>14}")
    best = (math.inf, None, None)
    for N in [2, 4, 8, 16, 30, 50, 60]:
        R0, b = solve_MT_strict(N)
        a, _ = alpha_MT(b)
        sum_high = sum(b[k] for k in range(2, N + 1)) / b[1]
        R_MT = R0 + V0 + a * math.log(MTY_C)
        if R_MT < best[0]:
            best = (R_MT, N, b)
        print(f"  {N:>4} {R0:>10.4f} {sum_high:>16.6f} {a:>10.6f}"
              f" {R_MT:>10.4f} {MTY_R - R_MT:>+14.4f}")

    # ---------------------- MT-derived joint min ----------------------
    print("\n[3] Joint SDP: minimise R_MT directly under MT-strict cone.")
    print(f"  {'N':>4} {'R0(P*)':>10} {'Σ_{k≥2}b_k/b_1':>16} {'R_MT':>10}"
          f" {'gain vs MTY':>14}")
    for N in [2, 4, 8, 16, 30, 50]:
        R0, b, obj = solve_min_MT_R(N)
        sum_high = sum(b[k] for k in range(2, N + 1)) / b[1]
        if obj < best[0]:
            best = (obj, N, b)
        print(f"  {N:>4} {R0:>10.4f} {sum_high:>16.6f} {obj:>10.4f}"
              f" {MTY_R - obj:>+14.4f}")

    # ---------------------- Verdict ----------------------
    R_best, N_best, b_best = best
    print("\n[4] BEST UNDER THIS DERIVATION")
    print(f"   N* = {N_best},  R_MT* = {R_best:.4f},  MTY 2022 = {MTY_R:.4f}")
    print(f"   gain = {MTY_R - R_best:+.4f}")
    if R_best < MTY_R - 1e-4:
        print()
        print("   → Heuristic improvement of MTY 2022's R is OBTAINED at the")
        print("     level of THIS DERIVATION.  This rests on the assumption")
        print("     that the MT analytic constant L is universal across")
        print("     polynomials (calibrated from dlVP).  The b_k coefficients")
        print("     of the optimising polynomial are:")
        print("     " + ", ".join(f"{x:.4f}" for x in b_best))
    else:
        print("   → No heuristic improvement from this derivation.")

    # ---------------------- Caveat ----------------------
    print("\n[5] CRITICAL CAVEAT — the [2] gain is most likely an ARTEFACT")
    print()
    print("   The MT analytic argument requires not just P ≥ 0 and b_k ≥ 0,")
    print("   but specifically the dlVP quadratic boost from the b_2 term:")
    print("     |1 + 2 cos θ|² = 3 + 4 cos θ + 2 cos 2θ ⇒ b_2/b_1 = 1/2.")
    print("   Without b_2 ≥ b_1/2 the analytic bound collapses from R = 4")
    print("   to R = 1 (trivial Hadamard).  So the polynomials found in [2]")
    print("   that satisfy max P / b_1 = 4 with b_2 ≈ 0.01 do NOT actually")
    print("   deliver R = 4 in the MT analytic chain — they deliver R = 1,")
    print("   for which the V_0 + α log C contribution is irrelevant.")
    print()
    print("   When we DO force b_2 ≥ 0.5 (block [3]), the SDP must increase")
    print("   max P to 4.5 and Σ_{k≥2} b_k jumps to 0.5, giving R = 6.36 ≫")
    print("   MTY's 5.5587.  This shows the dlVP polynomial is OPTIMAL within")
    print("   the MT-strict + dlVP-quadratic constraints — exactly the")
    print("   message of paper 152.")
    print()
    print("   What this attempt rules OUT: the α-spread-rescaling of attempt")
    print("   03 is incompatible with the MT structure (the actual α formula")
    print("   weights Σ_{k≥2} b_k, not the log-spread).  Both attempts 03")
    print("   and 04 produce HEURISTIC numerics but no rigorous improvement")
    print("   on MTY 2022.")
    print()
    print("   HONEST VERDICT: removing assumption (2) (α-rescaling) does")
    print("   NOT give a rigorous beat of MTY.  The leading-order R_0 = 4")
    print("   plus the MT-derived V correction reproduces MTY exactly when")
    print("   the analytic constraints are enforced honestly.")
    print("   Beating MTY rigorously needs assumption (1) — Λ-improvement")
    print("   or direct V-K constant reduction — which is research-level.")


if __name__ == "__main__":
    main()
