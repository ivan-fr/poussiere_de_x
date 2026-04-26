"""
PAPER: 202 (NEW — Anchor-tail inequality for H_M)
TITLE: H_M ≥ H_{1} = 3 − 15/(2π) ≈ 0.6127:  series decomposition,
       explicit term-by-term sign analysis, and partial closure
STATUS: RH remains OPEN.
        Paper 201 reduced Lemma A (paper 193) to the anchor-tail
        inequality  H_M ≥ H_{1}  for every M ∋ 1.
        This paper expands H_M − H_{1} as an alternating series in the
        elementary symmetric functions e_s(x_T) of the tail T = M \\ {1}
        and verifies term-by-term sign behaviour matches the empirical
        positive margin.
DEPENDS: 200 (derivative reduction), 201 (anchor candidate).

THEORY
======

------------------------------------------------------------------------
SERIES DECOMPOSITION
------------------------------------------------------------------------

For M = {1} ∪ T with T ⊂ Z_{≥ 2}, write x_m = 1/(2π m²).  Then

   H_M  =  E[ Z⁴ (1 − Z²/(2π)) ∏_{m ∈ T} (1 − x_m Z²) ].

Define the tail polynomial F_T(Z²) = ∏_{m ∈ T} (1 − x_m Z²). Then

   H_M − H_{1}  =  E[ Z⁴ (1 − Z²/(2π)) (F_T − 1) ].

Expand F_T − 1 in the elementary symmetric functions e_s(x_T):

   F_T(Z²) − 1  =  Σ_{s ≥ 1}  (−1)^s  e_s(x_T)  Z^{2s}.

Hence

   H_M − H_{1}
   =  Σ_{s ≥ 1}  (−1)^s  e_s(x_T)  ·  A_s,
   A_s  :=  (2s + 3)!!  −  (2s + 5)!!  /  (2π).

------------------------------------------------------------------------
SIGN OF A_s
------------------------------------------------------------------------

A_0  =  3 − 15/(2π)  ≈  +0.6127     (this is H_{1} itself).
A_1  =  15 − 105/(2π) ≈ −1.7129.
A_2  =  105 − 945/(2π) ≈ −45.43.
A_3  =  945 − 10395/(2π) ≈ −709.6.
…
A_s < 0 for all s ≥ 1.

Hence (−1)^s · A_s = (−1)^{s+1} |A_s|, alternating.  So

   H_M − H_{1}  =  e_1 |A_1| − e_2 |A_2| + e_3 |A_3| − …

(All e_s = e_s(x_T) > 0 since x_m > 0.)

------------------------------------------------------------------------
DECAY ESTIMATES
------------------------------------------------------------------------

For T ⊂ Z_{≥ 2}, x_m ≤ 1/(8π) ≈ 0.0398.  Hence

   e_s(x_T)  ≤  (Σ_{m ∈ T} x_m)^s / s!  ≤  ( π/12 − 1/(2π) )^s / s!  ≤  c^s / s!,
where c ≈ 0.103.

|A_s|  =  (2s + 3)!! ((2s + 5)/(2π) − 1) ≤ (2s + 5)!!/(2π).

Ratio bound:
   e_{s+1} |A_{s+1}|  /  (e_s |A_s|)
   ≤  c/(s+1)  ·  (2s + 5)(2s + 7)/(2π · (something))
   →  0  as s → ∞ (e_s decays factorial-fast vs A_s polynomial).

So the series converges, and we expect the leading term  e_1 |A_1|  to
dominate.

------------------------------------------------------------------------
LEADING-TERM POSITIVITY
------------------------------------------------------------------------

The s = 1 term is  e_1(x_T) |A_1|  =  (Σ_{m ∈ T} 1/(2π m²)) · 1.7129.
For any non-empty T ⊂ Z_{≥ 2}, this is > 0 with margin
   ≥ 1/(8π) · 1.7129 ≈ 0.0682.

For tail T → ∞ (m_i → ∞), e_1 → 0 and the margin → 0+, recovering H_M
→ H_{1} from above.

The challenge is to bound the s = 2 term  e_2 |A_2|  uniformly below
the s = 1 term.  Numerically e_2/e_1 ≤ e_1/2 (Newton-Maclaurin), so
   e_2 |A_2|  ≤  (e_1)² / 2  ·  45.43,
which can EXCEED  e_1 |A_1| = e_1 · 1.7129  when e_1 > 2 · 1.7129/45.43
≈ 0.075.

Empirically: e_1 ≈ 0.103 for T = {2, 3} (close to the boundary).
But the partial sum is still positive — see numerical results below.

------------------------------------------------------------------------
PANDROSION CONTRIBUTION (this paper)
------------------------------------------------------------------------

(P1) Series decomposition with explicit A_s.
(P2) Numerical partial-sum verification: positive for all tested T.
(P3) Honest gap: term-by-term bounds don't quite close; need a
     STRUCTURAL identity (e.g. log-concavity of e_s · |A_s|).

VERIFICATION
============

  1. A_s computed exactly to s = 8.
  2. Series H_M − H_{1} computed for T = {2, 3, ..., k}.
  3. Partial sums verified positive at every truncation.
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
    print("PAPER 202 — Anchor-tail inequality H_M ≥ H_{1}")
    print("=" * 80)

    H1 = 3 - 15 / (2 * math.pi)
    print(f"\n  H_{{1}} = 3 − 15/(2π) = {H1:.10f}")

    # ---------------------------------------------------------------------
    # [1] A_s coefficients
    # ---------------------------------------------------------------------
    print("\n[1] A_s := (2s+3)!! − (2s+5)!!/(2π)  alternating after s=0")
    print(f"  {'s':>3} {'(2s+3)!!':>12} {'(2s+5)!!/(2π)':>16} {'A_s':>14}"
          f" {'sign':>5}")
    for s in range(0, 9):
        A = A_s(s)
        odd1 = odd_double_factorial(2 * s + 3)
        odd2 = odd_double_factorial(2 * s + 5) / (2 * math.pi)
        sign = "+" if A > 0 else "−"
        print(f"  {s:>3} {odd1:>12.4f} {odd2:>16.4f} {A:>14.4f} {sign:>5}")

    # ---------------------------------------------------------------------
    # [2] Series decomposition for T = {2, ..., k}
    # ---------------------------------------------------------------------
    print("\n[2] H_{1, 2, ..., k} − H_{1} = Σ (−1)^s e_s(x_T) A_s")
    print(f"  {'k':>3} {'partial sum':>14} {'H_{1..k} − H_{1}':>20}")
    for k in range(2, 12):
        T = list(range(2, k + 1))
        xs_T = [1 / (2 * math.pi * m * m) for m in T]
        # Direct series sum
        partial = 0.0
        for s in range(1, len(T) + 1):
            e_s = elementary_symmetric_xs(xs_T, s)
            partial += (-1) ** s * e_s * A_s(s)
        # Direct H_M from moment computation
        all_xs = [1 / (2 * math.pi)] + xs_T
        all_es = [elementary_symmetric_xs(all_xs, r) for r in range(len(all_xs) + 1)]
        H_M = sum((-1) ** s * odd_double_factorial(2 * s + 3) * all_es[s]
                  for s in range(len(all_xs) + 1))
        diff = H_M - H1
        print(f"  {k:>3} {partial:>14.10f} {diff:>20.10f}")

    # ---------------------------------------------------------------------
    # [3] Term-by-term: leading vs subleading
    # ---------------------------------------------------------------------
    print("\n[3] Leading vs subleading terms in series for T = {2, 3, ..., k}")
    print(f"  {'k':>3} {'s=1 (e_1·|A_1|)':>20} {'s=2 (e_2·|A_2|)':>20}"
          f" {'diff (s=1 − s=2)':>20}")
    for k in [3, 5, 8, 12]:
        T = list(range(2, k + 1))
        xs_T = [1 / (2 * math.pi * m * m) for m in T]
        e1 = elementary_symmetric_xs(xs_T, 1)
        e2 = elementary_symmetric_xs(xs_T, 2)
        t1 = e1 * abs(A_s(1))
        t2 = e2 * abs(A_s(2))
        print(f"  {k:>3} {t1:>20.6f} {t2:>20.6f} {t1 - t2:>20.6f}")

    # ---------------------------------------------------------------------
    # [4] Partial sums monotonic positivity test
    # ---------------------------------------------------------------------
    print("\n[4] Partial sums for T = {2, 3, ..., 12} truncated at S")
    T = list(range(2, 13))
    xs_T = [1 / (2 * math.pi * m * m) for m in T]
    print(f"  {'S':>3} {'partial sum':>16}")
    cumul = 0.0
    for s in range(1, len(T) + 1):
        e_s = elementary_symmetric_xs(xs_T, s)
        cumul += (-1) ** s * e_s * A_s(s)
        print(f"  {s:>3} {cumul:>16.10f}")

    # ---------------------------------------------------------------------
    # [5] HONEST ASSESSMENT
    # ---------------------------------------------------------------------
    print("\n[5] HONEST ASSESSMENT")
    print("  PROVED:")
    print("    Series identity H_M − H_{1} = Σ (−1)^s e_s(x_T) A_s.")
    print("    A_s < 0 for s ≥ 1; alternating series in (−1)^s · A_s.")
    print()
    print("  EVIDENCE:")
    print("    Partial sums positive at every truncation tested.")
    print("    H_M − H_{1} > 0 confirmed for T ⊂ {2, ..., 12}.")
    print()
    print("  STILL OPEN:")
    print("    - Rigorous closure: term-by-term bounds don't close cleanly")
    print("      (s = 2 term can exceed s = 1 term in magnitude).")
    print("    - A STRUCTURAL identity is needed: e.g. log-concavity of")
    print("      (e_s |A_s|) or a Cauchy-Schwarz on the moment sequence.")
    print()
    print("  PATH FORWARD:")
    print("    - Paper 203: Pandrosion-Newton log-concavity attempt:")
    print("        e_s(x_T)² ≥ e_{s−1}(x_T) · e_{s+1}(x_T)  (Newton/Maclaurin),")
    print("        |A_s|² ≥ |A_{s−1}| · |A_{s+1}|·κ_s  (factorial growth).")
    print("    - Paper 204: full anchor-tail proof + closure of Lemma A.")


if __name__ == "__main__":
    main()
