/-
  Universitas Pandrosion — §36. **Möbius conjugacy of
  Pandrosion-Steffensen at `p = 2`.**

  Target. Prove the algebraic heart of the real-line McMullen
  statement for `p = 2`: under the Möbius change of coordinate
      `w = φ(s) := (s − α) / (s + α)`,
  the Steffensen accelerator `σ_{2,x}` becomes *exactly* the scaled
  squaring map
      `σ̃(w) = μ² · w²`,   `μ := (1 − α) / (1 + α)`.

  This is the analytic golden key: the Julia set of `w ↦ μ²w²` on
  `ℂP¹` is the circle `|μ²w| = 1`, which meets `ℝ` in exactly two
  points; the complement of this two-point Julia section, minus the
  countable orbit-preimages of the Möbius singularities `s = −1`
  and `s = −α²`, is the full-Lebesgue-measure Fatou set of `σ`,
  foliated into the two basins of `±α`.

  Scope of this module.

    §36.1  `sigma_p2_explicit` — the explicit rational form of the
           `p = 2` Steffensen step, defined directly (independent of
           the generic `steffensen_step` definition) and used as the
           analytic carrier of the conjugacy.

    §36.2  `sigma_p2_minus_alpha_cleared`,
           `sigma_p2_plus_alpha_cleared` — the **cleared-denominator
           factored identities** (pure polynomial algebra, modulo
           `x·α² = 1`):

               `(σ(s) − α) · 2α²(s+1)(xs+1)  =  (s − α)² (1 − α)²`,
               `(σ(s) + α) · 2α²(s+1)(xs+1)  =  (s + α)² (1 + α)²`.

    §36.3  `steffensen_mobius_conjugacy_p2_cleared` — **the headline
           identity**, cleared of denominators for fully polynomial
           content:

               `(σ(s) − α) · (s + α)² · (1 + α)²
                  = (σ(s) + α) · (s − α)² · (1 − α)²`.

           Together with the obvious non-vanishing of `(s ± α)` and
           `(1 ± α)` (given `|s| ≠ α` and `α ≠ 1`), this is
           equivalent to

               `(σ(s) − α)/(σ(s) + α) = μ² · ((s − α)/(s + α))²`,

           the headline Möbius conjugacy in ratio form.

  Roadmap. Discharging `RealMcMullenP2` from §35.4 unconditionally
  requires three further ingredients, outside the scope of this
  module but unlocked by §36.3:

    (A)  Bridge `sigma_p2_explicit x α s = steffensen_step x 2 s`
         on the non-degenerate set (tedious `Sp_p2` unfoldings —
         pure algebra, no new mathematical content).

    (B)  Iterated conjugacy `φ(σⁿ(s)) = μ^{2(2ⁿ−1)} · φ(s)^{2ⁿ}`
         — induction from §36.3.

    (C)  Measure-zero argument: the "bad set" is the countable union
         of orbit-preimages of `{−1, −α², ±α}`, and the Julia
         section `{s : |φ(s)| = 1/μ²}` is a two-point set on `ℝ`.

  (A)–(C) together discharge `RealMcMullenP2 x α` for `x > 1`.
-/

import Pandrosion.Core.SteffensenRealP2

namespace Pandrosion

open Filter Topology

/-! ============================================================
  §36.1  Explicit rational form of `σ_{2,x}`
============================================================ -/

section ExplicitSigma

/-- **Explicit Möbius-form Steffensen step at `p = 2`.**

    On the non-degenerate set `{s : (s + 1)(xs + 1) ≠ 0}`, the
    Pandrosion-Steffensen accelerator admits the closed form
        `σ_{2,x}(s) = s − (s − α)(s + α)(2xs + x + 1) / [2(s+1)(xs+1)]`.
    We define this as a standalone analytic object; the bridge to
    the generic `steffensen_step` (through `Sp_2 = 1 + s` and the
    Möbius square `h²(s) = ((x+1)s + 2)/(2xs + x + 1)`) is a pure
    unfolding exercise handled outside this module. -/
noncomputable def sigma_p2_explicit (x α s : ℝ) : ℝ :=
  s - (s - α) * (s + α) * (2 * x * s + x + 1) / (2 * (s + 1) * (x * s + 1))

end ExplicitSigma

/-! ============================================================
  §36.2  Cleared-denominator factored identities
============================================================ -/

section FactoredFormsCleared

/-- **Cleared-denominator factored form of `σ(s) − α`.**

    At any Pandrosion fixed point (`xα² = 1`), on the non-degenerate
    set `(s + 1)(xs + 1) ≠ 0`, the explicit Steffensen step satisfies
        `(σ(s) − α) · 2α²(s+1)(xs+1)  =  (s − α)² (1 − α)²`.

    Pure polynomial identity modulo `xα² = 1`. -/
theorem sigma_p2_minus_alpha_cleared
    (x α s : ℝ)
    (hs_ne_neg_one : s + 1 ≠ 0)
    (hs_xs1_ne : x * s + 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    (sigma_p2_explicit x α s - α) * (2 * α ^ 2 * (s + 1) * (x * s + 1))
      = (s - α) ^ 2 * (1 - α) ^ 2 := by
  unfold sigma_p2_explicit
  have h_denom_ne : (2 * (s + 1) * (x * s + 1) : ℝ) ≠ 0 :=
    mul_ne_zero (mul_ne_zero two_ne_zero hs_ne_neg_one) hs_xs1_ne
  have h_rewrite :
      (s - (s - α) * (s + α) * (2 * x * s + x + 1) / (2 * (s + 1) * (x * s + 1)) - α)
        * (2 * α ^ 2 * (s + 1) * (x * s + 1))
      = (s - α) * (2 * α ^ 2 * (s + 1) * (x * s + 1))
        - α ^ 2 * ((s - α) * (s + α) * (2 * x * s + x + 1)) := by
    field_simp
    ring
  rw [h_rewrite]
  linear_combination ((1 - 2 * α) * s ^ 2 - 2 * α * (1 - α) * s + α ^ 2) * hxα

/-- **Cleared-denominator factored form of `σ(s) + α`.**

    At any Pandrosion fixed point (`xα² = 1`), on the non-degenerate
    set `(s + 1)(xs + 1) ≠ 0`, the explicit Steffensen step satisfies
        `(σ(s) + α) · 2α²(s+1)(xs+1)  =  (s + α)² (1 + α)²`.

    Pure polynomial identity modulo `xα² = 1`. -/
theorem sigma_p2_plus_alpha_cleared
    (x α s : ℝ)
    (hs_ne_neg_one : s + 1 ≠ 0)
    (hs_xs1_ne : x * s + 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    (sigma_p2_explicit x α s + α) * (2 * α ^ 2 * (s + 1) * (x * s + 1))
      = (s + α) ^ 2 * (1 + α) ^ 2 := by
  unfold sigma_p2_explicit
  have h_denom_ne : (2 * (s + 1) * (x * s + 1) : ℝ) ≠ 0 :=
    mul_ne_zero (mul_ne_zero two_ne_zero hs_ne_neg_one) hs_xs1_ne
  have h_rewrite :
      (s - (s - α) * (s + α) * (2 * x * s + x + 1) / (2 * (s + 1) * (x * s + 1)) + α)
        * (2 * α ^ 2 * (s + 1) * (x * s + 1))
      = (s + α) * (2 * α ^ 2 * (s + 1) * (x * s + 1))
        - α ^ 2 * ((s - α) * (s + α) * (2 * x * s + x + 1)) := by
    field_simp
    ring
  rw [h_rewrite]
  linear_combination ((1 + 2 * α) * s ^ 2 + 2 * α * (1 + α) * s + α ^ 2) * hxα

end FactoredFormsCleared

/-! ============================================================
  §36.3  Headline Möbius conjugacy (cleared form)
============================================================ -/

section MobiusConjugacy

/-- **★★★ Möbius conjugacy of Pandrosion-Steffensen at `p = 2`
    (cleared-denominator form).**

    For `xα² = 1` and `s` in the non-degenerate set
    `(s+1)(xs+1) ≠ 0`, the ratio `φ(s) := (s − α)/(s + α)` satisfies
        `(σ(s) − α)·(s + α)²·(1 + α)²
          = (σ(s) + α)·(s − α)²·(1 − α)²`.

    Equivalently (on the set `s ≠ ±α`, `α ≠ ±1`):
        `φ(σ(s)) = μ² · φ(s)²`,   `μ := (1 − α)/(1 + α)`.

    This is the **analytic golden key** for the real-line McMullen
    a.e. entry at `p = 2`: `σ` is conjugate to the scaled squaring
    map `w ↦ μ²w²`, whose Julia set is a single circle on `ℂP¹` —
    cutting `ℝ` in at most two points, hence Lebesgue-negligible on
    the line. -/
theorem steffensen_mobius_conjugacy_p2_cleared
    (x α s : ℝ)
    (hα_ne_zero : α ≠ 0)
    (hs_ne_neg_one : s + 1 ≠ 0)
    (hs_xs1_ne : x * s + 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    (sigma_p2_explicit x α s - α) * (s + α) ^ 2 * (1 + α) ^ 2
      = (sigma_p2_explicit x α s + α) * (s - α) ^ 2 * (1 - α) ^ 2 := by
  have h_minus := sigma_p2_minus_alpha_cleared x α s hs_ne_neg_one hs_xs1_ne hxα
  have h_plus  := sigma_p2_plus_alpha_cleared  x α s hs_ne_neg_one hs_xs1_ne hxα
  have hα2_ne : α ^ 2 ≠ 0 := pow_ne_zero 2 hα_ne_zero
  have hK_ne : 2 * α ^ 2 * (s + 1) * (x * s + 1) ≠ 0 :=
    mul_ne_zero
      (mul_ne_zero (mul_ne_zero two_ne_zero hα2_ne) hs_ne_neg_one)
      hs_xs1_ne
  -- Multiply the goal by the common denominator K := 2α²(s+1)(xs+1):
  --   (σ − α)(s+α)²(1+α)² · K = (σ + α)(s−α)²(1−α)² · K
  -- Use h_minus and h_plus to collapse both sides.
  have h_key :
      (sigma_p2_explicit x α s - α) * (s + α) ^ 2 * (1 + α) ^ 2
        * (2 * α ^ 2 * (s + 1) * (x * s + 1))
      = (sigma_p2_explicit x α s + α) * (s - α) ^ 2 * (1 - α) ^ 2
        * (2 * α ^ 2 * (s + 1) * (x * s + 1)) := by
    linear_combination (s + α) ^ 2 * (1 + α) ^ 2 * h_minus
                       - (s - α) ^ 2 * (1 - α) ^ 2 * h_plus
  exact mul_right_cancel₀ hK_ne h_key

end MobiusConjugacy

end Pandrosion
