"""
PAPER: 204 (NEW — Term-by-term bound attempt for anchor-tail series)
TITLE: Why generic Newton-Maclaurin bounds fail to close H_M ≥ H_{1}:
       term s = 1 and term s = 2 are quantitatively comparable, forcing
       the proof to use the specific x_m = 1/(2π m²) structure.
STATUS: RH remains OPEN.  This paper does NOT close Lemma A and is NOT
        a proof of RH.  It identifies the analytic obstruction:
        generic positive-sequence bounds (Newton, Maclaurin, AM-GM) are
        insufficient because terms s = 1 and s = 2 of the alternating
        anchor-tail series have COMPARABLE magnitude.
DEPENDS: 200, 201, 202, 203 (anchor-tail framework).

THEORY
======

------------------------------------------------------------------------
THE OBSTRUCTION
------------------------------------------------------------------------

For T ⊂ {2, 3, ...} and the anchor-tail series

   H_M − H_{1}  =  Σ_{s ≥ 1}  (−1)^s · e_s(x_T) · A_s,

the term-by-term magnitudes for x_m = 1/(2π m²), T = {2, ..., k}:

   |term s = 1|   =   e_1 · |A_1|     ≈   0.0925 · 1.713  ≈  0.158
   |term s = 2|   =   e_2 · |A_2|     ≈   0.0034 · 45.43  ≈  0.155
   |term s = 3|   =   e_3 · |A_3|     ≈   ... ≈ 0.10
   |term s = 4|   =   ...

i.e. consecutive terms are of the SAME ORDER OF MAGNITUDE.  A "generic"
bound  Σ_{s ≥ 2} |term s| ≤ ρ · |term s = 1|  with ρ < 1 fails.

Concretely, by Maclaurin's inequality for positive x_T:
   e_s(x_T)  ≤  e_1(x_T)^s  /  s!.
And |A_s| = (2s+3)!! · ((2s+5)/(2π) − 1).  So
   |term s|  ≤  e_1^s  /  s!  ·  (2s+3)!!  ·  ((2s+5)/(2π) − 1).

This bound, summed over s ≥ 2, EXCEEDS |term s = 1| for e_1 ≥ 0.07.
Since e_1 can reach π/12 − 1/(2π) ≈ 0.103 for T = {2, 3, ...}, the
generic bound DOES NOT close.

------------------------------------------------------------------------
WHY THE NET SUM IS STILL POSITIVE
------------------------------------------------------------------------

Empirically (paper 202, 203):  H_M − H_{1} > 0  always.
The cancellation between alternating terms is large; the residual
is small but POSITIVE.

So the proof must use the SPECIFIC SEQUENCE  x_m = 1/(2π m²),  not
generic positive reals.  The key is the DIOPHANTINE / number-theoretic
structure:  e_s(x_T) for x_m = 1/(2π m²) is related to evaluations of
ζ(2), ζ(4), ... (Faulhaber-type identities).

------------------------------------------------------------------------
PANDROSION CONTRIBUTION (this paper)
------------------------------------------------------------------------

(P1) Quantify why generic bounds fail: term s = 1 ≈ term s = 2.
(P2) Identify the structural avenue: x_m = 1/(2π m²) gives e_s in terms
     of ζ-values and Bernoulli numbers.
(P3) Numerical confirmation: net positivity comes from CANCELLATION.

VERIFICATION
============

  1. Term magnitudes |term s| computed for T = {2, ..., k}, k up to 12.
  2. Generic Maclaurin bound vs actual.
  3. Confirm: rigorous closure needs structural input.
"""
from __future__ import annotations
import math
import itertools


def odd_double_factorial(n):
    out = 1.0
    for q in range(n, 0, -2):
        out *= q
    return out


def A_s(s):
    return odd_double_factorial(2 * s + 3) - odd_double_factorial(2 * s + 5) / (2 * math.pi)


def elementary_symmetric_xs(xs, s):
    if s == 0:
        return 1.0
    if s > len(xs):
        return 0.0
    total = 0.0
    for c in itertools.combinations(xs, s):
        total += math.prod(c)
    return total


def main():
    print("=" * 80)
    print("PAPER 204 — Term-by-term bound attempt for anchor-tail series")
    print("=" * 80)

    pi = math.pi
    H1 = 3 - 15 / (2 * pi)
    print(f"\n  H_{{1}} = {H1:.10f}")

    # ---------------------------------------------------------------------
    # [1] Term magnitudes for T = {2, ..., k}
    # ---------------------------------------------------------------------
    print("\n[1] Term magnitudes |term s| = e_s · |A_s| for T = {2, ..., k}")
    print(f"  {'k':>3} {'s=1':>10} {'s=2':>10} {'s=3':>10} {'s=4':>10}"
          f" {'s=5':>10} {'partial sum':>14}")
    for k in [2, 4, 6, 8, 10, 12]:
        T = list(range(2, k + 1))
        xs_T = [1 / (2 * pi * m * m) for m in T]
        terms = []
        for s in range(1, min(6, len(T) + 1)):
            e_s = elementary_symmetric_xs(xs_T, s)
            terms.append(e_s * abs(A_s(s)))
        # Partial sum (all s)
        partial = 0.0
        for s in range(1, len(T) + 1):
            e_s = elementary_symmetric_xs(xs_T, s)
            partial += (-1) ** s * e_s * A_s(s)
        # pad terms list to 5
        while len(terms) < 5:
            terms.append(0.0)
        print(f"  {k:>3} {terms[0]:>10.5f} {terms[1]:>10.5f} {terms[2]:>10.5f}"
              f" {terms[3]:>10.5f} {terms[4]:>10.5f} {partial:>14.6f}")

    # ---------------------------------------------------------------------
    # [2] Generic Maclaurin upper bound vs actual (T = {2, ..., 8})
    # ---------------------------------------------------------------------
    print("\n[2] Generic Maclaurin upper bound  e_s ≤ e_1^s / s!")
    print("    Comparison to actual e_s for T = {2, ..., 8}.")
    T = list(range(2, 9))
    xs_T = [1 / (2 * pi * m * m) for m in T]
    e1 = sum(xs_T)
    print(f"  e_1 = {e1:.6f}")
    print(f"  {'s':>3} {'e_s actual':>14} {'e_1^s/s!  bound':>18} {'ratio':>10}")
    for s in range(1, 8):
        actual = elementary_symmetric_xs(xs_T, s)
        bound = e1 ** s / math.factorial(s)
        ratio = actual / bound if bound > 0 else 0.0
        print(f"  {s:>3} {actual:>14.6e} {bound:>18.6e} {ratio:>10.4f}")

    # ---------------------------------------------------------------------
    # [3] Bound on Σ_{s ≥ 2} |term s| via Maclaurin
    # ---------------------------------------------------------------------
    print("\n[3] Maclaurin-bounded tail:  Σ_{s ≥ 2} (e_1^s / s!) · |A_s|")
    print(f"  {'k':>3} {'e_1':>10} {'|term 1|':>12} {'Maclaurin tail':>16}"
          f" {'tail/|term 1|':>16}")
    for k in [4, 6, 8, 10, 12]:
        T = list(range(2, k + 1))
        xs_T = [1 / (2 * pi * m * m) for m in T]
        e1 = sum(xs_T)
        term1 = e1 * abs(A_s(1))
        tail = 0.0
        for s in range(2, 12):
            tail += e1 ** s / math.factorial(s) * abs(A_s(s))
        ratio = tail / term1 if term1 > 0 else float('inf')
        print(f"  {k:>3} {e1:>10.6f} {term1:>12.6f} {tail:>16.6f}"
              f" {ratio:>16.4f}")

    # ---------------------------------------------------------------------
    # [4] Why net is still positive: cancellation
    # ---------------------------------------------------------------------
    print("\n[4] Cancellation: cumulative partial sums")
    print(f"    For T = {{2, ..., 12}}, partial sums truncated at S")
    T = list(range(2, 13))
    xs_T = [1 / (2 * pi * m * m) for m in T]
    cum = 0.0
    print(f"  {'S':>3} {'partial sum':>16} {'sign':>5}")
    for s in range(1, len(T) + 1):
        e_s = elementary_symmetric_xs(xs_T, s)
        cum += (-1) ** s * e_s * A_s(s)
        print(f"  {s:>3} {cum:>16.10f} {('+' if cum > 0 else '−'):>5}")

    # ---------------------------------------------------------------------
    # [5] Structural input: e_s(x_T) and ζ-values
    # ---------------------------------------------------------------------
    print("\n[5] Structural input: e_s for x_m = 1/(2π m²) and ζ(2k)")
    print("    For T = {2, 3, ..., ∞}: e_1(x_T) = (1/(2π))(ζ(2) − 1) = π/12 − 1/(2π)")
    e1_inf = math.pi / 12 - 1 / (2 * pi)
    print(f"    e_1(∞) = {e1_inf:.6f}")
    print("    Higher e_s relate to multiple zeta values (Apéry-type sums).")
    print("    Specifically e_s involves the s-th elementary symmetric")
    print("    function of {1/m² : m ≥ 2}.")
    print("    These are EXPLICITLY computable but the link to the cancel-")
    print("    lation Σ (−1)^s e_s A_s ≥ 0 is non-trivial.")

    # ---------------------------------------------------------------------
    # [6] HONEST ASSESSMENT
    # ---------------------------------------------------------------------
    print("\n[6] HONEST ASSESSMENT")
    print()
    print("  THIS PAPER IS NOT A PROOF OF RH.  RH remains OPEN.")
    print()
    print("  WHAT THIS PAPER DOES:")
    print("    Quantifies the obstruction to a generic-Newton-Maclaurin")
    print("    proof of the anchor-tail inequality H_M ≥ H_{1}.")
    print()
    print("    Specifically: for T = {2, ..., k} and the natural sequence")
    print("    x_m = 1/(2π m²), the terms in the alternating series")
    print("       H_M − H_{1} = Σ (−1)^s · e_s · A_s")
    print("    satisfy  |term s = 1| ≈ |term s = 2|  in magnitude.")
    print()
    print("    Hence Σ_{s ≥ 2} |term s| > |term s = 1|,  defeating any")
    print("    proof strategy that relies solely on bounding the tail by")
    print("    the leading term.  The net positivity comes from CANCEL-")
    print("    LATION, not domination.")
    print()
    print("  WHAT WOULD CLOSE THE INEQUALITY:")
    print("    A structural identity tied to the Diophantine sequence")
    print("    x_m = 1/(2π m²).  For instance:")
    print("    - A Riemann ζ-functional identity.")
    print("    - A Pandrosion-Polya-Schur preserver argument.")
    print("    - A Mehler-Hermite expansion with explicit signs.")
    print()
    print("  STILL OPEN:")
    print("    - Anchor-tail inequality (rigorous).")
    print("    - Lemma A → Wronskian barrier → RH.")
    print("    - RH itself: unsolved Millennium Problem.")
    print()
    print("  PATH FORWARD:")
    print("    - Paper 205+: investigate Riemann ζ-structural input on")
    print("      e_s(1/(2π m²)) sums.")
    print("    - This is research-level; no quick path to closure.")


if __name__ == "__main__":
    main()
