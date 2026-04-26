"""
PAPER: 199 (NEW — Wronskian barrier synthesis: Lemma A + Lemma B)
TITLE: Combining Lemma A (initial-segment extremality) and Lemma B
       (R_k > exp(−π/4) > 0) to close the Wronskian root-barrier
       reduction toward RH
STATUS: RH remains OPEN.
        This paper assembles the chain papers 189 → 197 into a single
        statement: IF Lemma A holds (paper 193, conjectural extremality
        of the initial segment), then Lemma B (papers 194-198, almost-
        rigorous via the lobe decay bound) closes the Wronskian
        positivity barrier for ALL admissible lattice subsets M.
DEPENDS: 189 (Wronskian Vandermonde), 190 (lattice Q barrier),
         193 (Lemma A), 194 (asymptotic R_k → exp(−π/4)),
         195 (recurrence), 196 (lobe decomposition),
         197 (lobe decay), 198 (Mehler-Hermite closed form).

THEORY
======

------------------------------------------------------------------------
THE WRONSKIAN BARRIER PROGRAM (papers 189-199 outline)
------------------------------------------------------------------------

For RH (or its weak Wronskian / Schmidt reformulation):

  Step 1 (paper 189, 190):
    Reduce sign analysis on the Riemann ξ Wronskian to positivity of
       q_M(y)  :=  Q_k(2π m_1² y, ..., 2π m_k² y)     for all M, y ≥ 1.

  Step 2 (paper 192):
    Express q_M(y) =  Σ_{r=0}^k (−1)^{k−r} (2(k−r)+1)!! e_r(2π m_i²) y^r
    with shifted coefficient
       B_j(M)  :=  [x^j] q_M(1 + x).
    Then q_M(y) > 0 for y ≥ 1 follows from B_j(M) ≥ 0 ∀ j and B_0(M) > 0.

  Step 3 (Lemma A, paper 193):
    Initial-segment extremality:
       B_j(M) / e_k(M)  ≥  B_0(M_k) / e_k(M_k),     M_k := {1,...,k}.
    Reduces to a one-parameter inequality.

  Step 4 (Lemma B, papers 194-197):
    R_k = B_0(M_k) / e_k(M_k) > 0.
    R_k decreases monotonically to exp(−π/4) ≈ 0.4559.
    Reduced to T_k = E[Z⁴ P_k(Z²)] ≥ 0  (paper 195).
    T_k ≥ 0 reduced to Gaussian-lobe dominance (paper 196).
    Lobe decay L_{j+1}/L_j ≤ ρ < 1 with ρ ≤ 0.378 numerical, < 0.6 rigorous
    (paper 197).
    Hence T_k ≥ (1 − ρ) L_0 > 0.

  Step 5 (this paper):
    Combine A + B.  Then q_M(y) > 0 for all M, y ≥ 1.
    The Wronskian barrier holds.

------------------------------------------------------------------------
WHAT THIS DOES (AND DOESN'T) IMPLY
------------------------------------------------------------------------

The Wronskian barrier is a NECESSARY (but probably not sufficient) step
for RH via the Pandrosion-Wronskian-Schmidt reduction.  Closing it would
remove one structural obstruction.

It does NOT directly imply RH.  Even with full Wronskian positivity,
RH requires additional spectral information (zero distribution).

This paper documents the chain status; it does not claim RH.

------------------------------------------------------------------------
PANDROSION CONTRIBUTION (this paper)
------------------------------------------------------------------------

(W1) Assemble the full chain 189 → 198.
(W2) Identify what's RIGOROUS vs CONJECTURAL:
     - Lemma A (paper 193): conjectural, supported by exhaustive search.
     - Lemma B (papers 194-197): almost rigorous, with explicit bounds.
     - Step 5 (this paper): conditional on A.
(W3) State the remaining open problems.

VERIFICATION
============

  1. Verify R_k > exp(−π/4) numerically for k up to 50.
  2. Verify Lemma A on extended search (already in paper 193).
  3. Combined bound:  q_M(1) > 0 for all tested M.
"""
from __future__ import annotations
import math
import itertools


def odd_double_factorial(n):
    out = 1.0
    for m in range(n, 0, -2):
        out *= m
    return out


def reciprocal_es(k):
    coeffs = [1.0] + [0.0] * k
    for m in range(1, k + 1):
        x = 1.0 / (2 * math.pi * m * m)
        for s in range(k, 0, -1):
            coeffs[s] += coeffs[s - 1] * x
    return coeffs


def R_k(k):
    """R_k = B_0({1,...,k}) / e_k({1,...,k}) = sum (-1)^s (2s+1)!! E_s(k)."""
    coeffs = reciprocal_es(k)
    return sum((-1) ** s * odd_double_factorial(2 * s + 1) * coeffs[s]
               for s in range(k + 1))


def main():
    print("=" * 80)
    print("PAPER 199 — Wronskian barrier synthesis")
    print("=" * 80)

    # ---------------------------------------------------------------------
    # [1] R_k → exp(−π/4) verification
    # ---------------------------------------------------------------------
    R_inf = math.exp(-math.pi / 4)
    print(f"\n[1] R_k decreases to exp(−π/4) ≈ {R_inf:.10f}")
    print(f"  {'k':>4} {'R_k':>16} {'R_k − exp(−π/4)':>20}")
    for k in [1, 2, 3, 5, 10, 20, 30, 50]:
        Rk = R_k(k)
        print(f"  {k:>4} {Rk:>16.10f} {Rk - R_inf:>20.4e}")

    # ---------------------------------------------------------------------
    # [2] Lemma A spot-check (paper 193 territory)
    # ---------------------------------------------------------------------
    print("\n[2] Lemma A spot-check: B_0({1..k})/e_k({1..k}) is the minimum")
    print("    over all M ⊂ {1, ..., M_max} with |M| = k.")
    M_max = 8
    print(f"  M_max = {M_max}")
    print(f"  {'k':>3} {'min over M':>14} {'M-witness':>22} {'init-segment ratio':>20}")
    for k in range(2, 6):
        best = (float('inf'), None)
        init_val = None
        for ms in itertools.combinations(range(1, M_max + 1), k):
            a = [2 * math.pi * m * m for m in ms]
            # B_0 = sum_{r=0}^k (-1)^{k-r} (2(k-r)+1)!! e_r(a)
            B0 = 0.0
            for r in range(k + 1):
                e_r = sum(math.prod(c)
                          for c in itertools.combinations(a, r))
                B0 += (-1) ** (k - r) * odd_double_factorial(
                    2 * (k - r) + 1) * e_r
            e_k_val = math.prod(a)
            ratio = B0 / e_k_val
            if ms == tuple(range(1, k + 1)):
                init_val = ratio
            if ratio < best[0]:
                best = (ratio, ms)
        print(f"  {k:>3} {best[0]:>14.6f} {str(best[1]):>22}"
              f" {init_val:>20.6f}")

    # ---------------------------------------------------------------------
    # [3] Combined synthesis
    # ---------------------------------------------------------------------
    print("\n[3] Combined synthesis: Lemma A + Lemma B")
    print("    For every admissible lattice M of size k:")
    print("      Lemma A:  B_0(M)/e_k(M) ≥ B_0({1..k})/e_k({1..k}) = R_k.")
    print(f"      Lemma B:  R_k ≥ {R_inf:.4f}  (decreases to exp(−π/4)).")
    print("    Hence  B_0(M)/e_k(M) ≥ exp(−π/4) > 0  for every M.")
    print("    This implies  q_M(1) > 0  for every admissible lattice M.")

    # ---------------------------------------------------------------------
    # [4] Status table
    # ---------------------------------------------------------------------
    print("\n[4] Status of each step in the chain")
    rows = [
        ("189: Wronskian → q_M(y) > 0",        "rigorous reduction"),
        ("190: q_M(y) > 0 for y ≥ 1",          "rigorous reduction"),
        ("192: B_j(M) ≥ 0",                    "rigorous reduction"),
        ("193: Lemma A (initial extremality)", "CONJECTURAL — verified"),
        ("                                  ", "exhaustively for M_max ≤ 12"),
        ("194: R_k → exp(−π/4) limit",         "rigorous (algebraic)"),
        ("195: R_k − R_{k+1} = x_{k+1} T_k",    "rigorous"),
        ("196: T_k = lobe sum",                 "rigorous"),
        ("197: lobe decay L_{j+1}/L_j ≤ 0.6",  "rigorous (loose) /"),
        ("                                  ", "tighter numerical < 0.4"),
        ("198: Mehler-Hermite closed form",     "computational"),
        ("199: A + B ⟹ q_M(1) > 0",          "conditional on A"),
    ]
    for step, status in rows:
        print(f"  {step:<40} {status}")

    # ---------------------------------------------------------------------
    # [5] HONEST ASSESSMENT
    # ---------------------------------------------------------------------
    print("\n[5] HONEST ASSESSMENT")
    print("  RH itself: remains OPEN.")
    print()
    print("  RIGOROUS (papers 189-198):")
    print("    R_k decreases to exp(−π/4) > 0 (algebraic reduction).")
    print("    T_k ≥ (1 − ρ) L_0 > 0 with ρ ≤ 0.6 (lobe decay).")
    print("    Hence Lemma B is RIGOROUS to within the loose ρ bound.")
    print()
    print("  CONJECTURAL:")
    print("    Lemma A (initial-segment extremality, paper 193).")
    print("    Verified exhaustively for M_max ≤ 12 (paper 193).")
    print("    The full Wronskian barrier is conditional on A.")
    print()
    print("  STILL OPEN:")
    print("    - Lemma A as a rigorous theorem.")
    print("    - Connection of the Wronskian barrier to RH itself")
    print("      (additional spectral conditions needed).")
    print()
    print("  PATH FORWARD:")
    print("    - Paper 200: rigorous Lemma A via Schur-Horn / majorization")
    print("      argument (the e_r symmetric functions form a totally")
    print("      positive sequence under the M-ordering).")
    print("    - Paper 201: connection of Wronskian barrier to ξ-zero")
    print("      distribution.")


if __name__ == "__main__":
    main()
