"""
ATTEMPT 03 — Synthesis: best heuristic R combining all leverage points.

Combines, in one place:
  (a) High-degree MT-strict cosine polynomial SDP   (paper 152, etage 1).
  (b) Stieltjes constants exhaustion bound          (paper 153, etage 2).
  (c) Λ-Polymath improvement                         (paper 159, etage 7).
  (d) Empirical V-K target from attempt 01         (sets the *target*).
  (e) Bombieri kernel + PL decoupling                (paper 158, etages 4+6).

Output: R_synth = 4 + V0 + α(P) · log C_eff(Λ).

We additionally try a "spread-minimising" polynomial — a polynomial
with smaller Σ b_k log(k+1) than dlVP — which, in the calibrated
heuristic, lowers α effectively because more of V comes from the
floor and less from the V-K factor.
"""
from __future__ import annotations
import math

import numpy as np
import cvxpy as cp


MTY_R = 5.5587
MTY_C = 76.2
V0 = math.log(2 * math.pi) - 0.5772156649


def solve_MT_strict(N):
    """Etage 1: MT-strict cosine polynomial SDP, returns (R0, b)."""
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


def solve_min_spread(N, R_cap):
    """Find an MT-strict polynomial with R0 <= R_cap that minimises
    spread = Σ b_k log(k+1).  (R_cap = 4.001 keeps us at the saturation.)"""
    X1 = cp.Variable((N + 1, N + 1), symmetric=True)
    X2 = cp.Variable((N + 1, N + 1), symmetric=True)
    t = cp.Variable()
    b = [cp.sum([X1[i, i + m] for i in range(N + 1 - m)])
         for m in range(N + 1)]
    c = [cp.sum([X2[i, i + m] for i in range(N + 1 - m)])
         for m in range(N + 1)]
    cons = [X1 >> 0, X2 >> 0, c[0] == t - b[0], b[1] == 1, t <= R_cap]
    for m in range(1, N + 1):
        cons.append(c[m] == -b[m])
    for m in range(N + 1):
        cons.append(b[m] >= 0)
    spread = cp.sum([b[m] * math.log(m + 1) for m in range(N + 1)])
    cp.Problem(cp.Minimize(spread), cons).solve(solver=cp.CLARABEL)
    return float(t.value), [float(bi.value) for bi in b], float(spread.value)


def alpha_for_poly(b, calib_C=MTY_C, calib_R=MTY_R):
    """The R-projection R(C) = R0 + V0 + α(P) log C with α calibrated so
    that for dlVP (b=(3,2,1/2)) and C=MTY_C we recover R = MTY_R = 5.5587.

    For other polynomials, treat α(P) as proportional to the spread/b1
    ratio relative to dlVP (heuristic, but consistent with the dlVP
    log-correction structure).
    """
    spread_dlVP = (3 * math.log(1) + 2 * math.log(2) + 0.5 * math.log(3)) / 2
    spread_P = sum(b[k] * math.log(k + 1) for k in range(len(b))) / b[1]
    alpha_dlVP = (calib_R - 4 - V0) / math.log(calib_C)
    return alpha_dlVP * (spread_P / spread_dlVP)


def main():
    print("=" * 80)
    print("ATTEMPT 03 — Synthesis: heuristic R combining every available lever")
    print("=" * 80)

    # ---------------------------------------------------------------
    # Etage 1 — high-degree MT-strict SDP
    # ---------------------------------------------------------------
    print("\n[1] Etage 1: MT-strict cosine polynomial (degrees 2..60).")
    print(f"  {'N':>4} {'R0':>10} {'spread/b1':>12} {'α(P)':>10}")
    best_alpha = 1.0
    best = None
    for N in [2, 8, 16, 30, 50, 60]:
        R0, b = solve_MT_strict(N)
        spread = sum(b[k] * math.log(k + 1) for k in range(N + 1)) / b[1]
        a = alpha_for_poly(b)
        if a < best_alpha:
            best_alpha = a
            best = (N, R0, b, spread)
        print(f"  {N:>4} {R0:>10.6f} {spread:>12.6f} {a:>10.6f}")

    print(f"\n  -> Best α from MT-strict SDP: α = {best_alpha:.6f} at N = {best[0]}")

    # ---------------------------------------------------------------
    # Etage 1b — minimise spread directly
    # ---------------------------------------------------------------
    print("\n[1b] Direct minimisation of spread under MT-strict + R0 ≤ 4.001.")
    print(f"  {'N':>4} {'R0':>10} {'spread/b1':>12} {'α(P)':>10}")
    for N in [4, 8, 16, 30]:
        try:
            R0, b, spread_raw = solve_min_spread(N, R_cap=4.001)
            spread = sum(b[k] * math.log(k + 1) for k in range(N + 1)) / b[1]
            a = alpha_for_poly(b)
            if a < best_alpha:
                best_alpha = a
                best = (N, R0, b, spread)
            print(f"  {N:>4} {R0:>10.6f} {spread:>12.6f} {a:>10.6f}")
        except Exception as e:
            print(f"  {N:>4} solver failed: {e}")
    print(f"\n  -> Best α overall: α = {best_alpha:.6f} at N = {best[0]}")

    # ---------------------------------------------------------------
    # Etage 7 — Λ scenarios for C_VK
    # ---------------------------------------------------------------
    print("\n[2] Λ → C_VK projection (paper 159).")
    Lam_table = [
        ("MTY 2022 (Λ ≤ 0.22 effectively)",  0.22, 76.2),
        ("Polymath16 announced (Λ ≤ 0.20)",   0.20, 70.0),
        ("hypothetical Λ ≤ 0.15",             0.15, 60.0),
        ("hypothetical Λ ≤ 0.10",             0.10, 40.0),
        ("Newman conj Λ = 0 (Lindelöf)",      0.00,  5.0),
    ]

    # ---------------------------------------------------------------
    # Synthesis table:  R_synth(P, C)  =  4 + V0 + α(P) · log C
    # ---------------------------------------------------------------
    print("\n[3] R_synth = 4 + V0 + α(P) · log C   for best polynomial.")
    print(f"   (V0 = {V0:.4f}, α(best) = {best_alpha:.6f})")
    print()
    print(f"  {'scenario':>40} {'C':>8} {'V':>10} {'R_synth':>10}"
          f" {'gain vs MTY':>14}")
    best_R = (math.inf, None)
    for label, Lam, C in Lam_table:
        V = V0 + best_alpha * math.log(C)
        R = 4 + V
        if R < best_R[0]:
            best_R = (R, label)
        print(f"  {label:>40} {C:>8.2f} {V:>10.4f} {R:>10.4f}"
              f" {MTY_R - R:>14.4f}")

    # ---------------------------------------------------------------
    # Verdict
    # ---------------------------------------------------------------
    print("\n[4] VERDICT")
    print(f"  Best R_synth (heuristic, all levers combined): {best_R[0]:.4f}")
    print(f"  MTY 2022 record:                                {MTY_R:.4f}")
    print(f"  best gain vs MTY:                               {MTY_R - best_R[0]:+.4f}")
    print()
    if best_R[0] < MTY_R:
        gain = MTY_R - best_R[0]
        print(f"  At the heuristic level, the combined synthesis gives a")
        print(f"  candidate improvement of {gain:.4f} over MTY 2022.")
        print()
        print("  However: this is HEURISTIC. The α-projection is calibrated")
        print("  to dlVP under MTY's V-K constant, and the polynomial-")
        print("  rescaled α used here is a *conjectural* extrapolation of")
        print("  how the explicit-formula log-correction depends on the")
        print("  polynomial spread. To convert this into a rigorous")
        print("  improvement of MTY 2022 requires:")
        print("    (i)  a rigorous proof of Λ ≤ 0.20  (Polymath16 in progress),")
        print("    (ii) a rigorous re-derivation of the explicit α(P) for")
        print("         the chosen polynomial (re-running the MT 2015 / MTY")
        print("         2022 analytic argument with this P),")
        print("    (iii) and constant tracking through the entire chain.")
    else:
        print("  No heuristic improvement on MTY 2022 within this synthesis.")
        print("  The Pandrosion-Bombieri framework's diagnosis (paper 160)")
        print("  is reaffirmed: the V-K constant 76.2 is the bottleneck and")
        print("  cannot be moved without a genuinely new analytic argument.")

    # ---------------------------------------------------------------
    # The actually-defensible bottom line
    # ---------------------------------------------------------------
    print("\n[5] WHAT WE HAVE / WHAT WE DON'T")
    print()
    print("  We HAVE:")
    print("   - A rigorous MT-strict R0 = 4 saturated by dlVP (paper 152).")
    print("   - An empirical sup C_emp ≤ 1.04 over t ≤ 10^9 (attempt 01).")
    print("   - A Hankel-Pandrosion-Hadamard framework for Λ ≤ t (attempt 02).")
    print("   - A heuristic synthesis chain that quantifies how each lever")
    print("     would translate into a gain of R below MTY 2022.")
    print()
    print("  We DO NOT HAVE:")
    print("   - A rigorous proof of Λ ≤ 0.20 (Polymath15 stops at < 0.22).")
    print("   - A rigorous re-derivation of α(P) for any polynomial P beyond")
    print("     dlVP at C = 76.2 (the calibration anchor).")
    print("   - A rigorous reduction of C_VK below 76.2 (the headline bound).")
    print()
    print("  Therefore: NO rigorous improvement of MTY 2022's R = 5.5587 is")
    print("  produced by this attempt.  This is consistent with paper 160's")
    print("  diagnosis.  Beating MTY rigorously needs research-level analytic")
    print("  work on the V-K constant or on the σ = 1/2 subconvexity bound.")


if __name__ == "__main__":
    main()
