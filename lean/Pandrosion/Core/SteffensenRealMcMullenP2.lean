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

/-! ============================================================
  §36.4  Bridge to the generic `steffensen_step`
============================================================ -/

section BridgeToGeneric

/-- **`p = 2` closed form of `pandrosion_h`.**
    For `x ≠ 0` and `s + 1 ≠ 0`,
        `h_{2,x}(s) = (x·s + 1) / (x·(s + 1))`.
    Normalising form (common denominator `x·(s+1)`) preferred for the
    nested-iterate chain of §36.4. -/
theorem pandrosion_h_p2_closed_form
    (x : ℝ) (hx : x ≠ 0) (s : ℝ) (hs : s + 1 ≠ 0) :
    pandrosion_h x 2 s = (x * s + 1) / (x * (s + 1)) := by
  unfold pandrosion_h
  rw [Sp_p2]
  have hxs1 : x * (s + 1) ≠ 0 := mul_ne_zero hx hs
  have h_1s_eq : (1 + s : ℝ) = s + 1 := by ring
  rw [h_1s_eq]
  field_simp
  ring

/-- **First difference of `h`.** `h(s) − s = (1 − x·s²)/(x·(s+1))`. -/
theorem pandrosion_h_p2_minus_s
    (x : ℝ) (hx : x ≠ 0) (s : ℝ) (hs : s + 1 ≠ 0) :
    pandrosion_h x 2 s - s = (1 - x * s ^ 2) / (x * (s + 1)) := by
  rw [pandrosion_h_p2_closed_form x hx s hs]
  have hxs1 : x * (s + 1) ≠ 0 := mul_ne_zero hx hs
  field_simp
  ring

/-- **`h(s) + 1` in closed form.** Used to certify that the second-level
    iterate `h(h(s))` has non-vanishing denominator. -/
theorem pandrosion_h_p2_plus_one
    (x : ℝ) (hx : x ≠ 0) (s : ℝ) (hs : s + 1 ≠ 0) :
    pandrosion_h x 2 s + 1 = (2 * x * s + x + 1) / (x * (s + 1)) := by
  rw [pandrosion_h_p2_closed_form x hx s hs]
  have hxs1 : x * (s + 1) ≠ 0 := mul_ne_zero hx hs
  field_simp
  ring

/-- **`p = 2` closed form of `h(h(s))`.**
    For `x ≠ 0`, `s + 1 ≠ 0`, and the Möbius-iterate denominator
    `2·x·s + x + 1 ≠ 0`,
        `h(h(s)) = ((x + 1)·s + 2) / (2·x·s + x + 1)`.
    This is the Möbius-matrix square `[[1, 1/x], [1, 1]]^2`. -/
theorem pandrosion_h_p2_hh_closed_form
    (x : ℝ) (hx : x ≠ 0) (s : ℝ) (hs : s + 1 ≠ 0)
    (hq : 2 * x * s + x + 1 ≠ 0) :
    pandrosion_h x 2 (pandrosion_h x 2 s)
      = ((x + 1) * s + 2) / (2 * x * s + x + 1) := by
  have hxs1_ne : x * (s + 1) ≠ 0 := mul_ne_zero hx hs
  have h_plus_one := pandrosion_h_p2_plus_one x hx s hs
  have h_hs_plus_one_ne : pandrosion_h x 2 s + 1 ≠ 0 := by
    rw [h_plus_one]; exact div_ne_zero hq hxs1_ne
  rw [pandrosion_h_p2_closed_form x hx _ h_hs_plus_one_ne]
  rw [h_plus_one, pandrosion_h_p2_closed_form x hx s hs]
  field_simp
  ring

/-- **Closed form of `steffensen_denom` at `p = 2`.**
        `steffensen_denom x 2 s
           = 2·(x·s + 1)·(x·s² − 1) / (x·(s+1)·(2·x·s + x + 1))`.
    Pure algebra from the three Möbius closed forms above. -/
theorem steffensen_denom_p2_closed_form
    (x : ℝ) (hx : x ≠ 0) (s : ℝ) (hs : s + 1 ≠ 0)
    (hq : 2 * x * s + x + 1 ≠ 0) :
    steffensen_denom x 2 s
      = 2 * (x * s + 1) * (x * s ^ 2 - 1)
          / (x * (s + 1) * (2 * x * s + x + 1)) := by
  unfold steffensen_denom
  rw [pandrosion_h_p2_hh_closed_form x hx s hs hq]
  rw [pandrosion_h_p2_closed_form x hx s hs]
  have hxs1_ne : x * (s + 1) ≠ 0 := mul_ne_zero hx hs
  field_simp
  ring

/-- **Non-vanishing of `steffensen_denom` at `p = 2`.**
    On the non-degenerate set — `s + 1 ≠ 0`, `x·s + 1 ≠ 0`,
    `2·x·s + x + 1 ≠ 0`, and `s ≠ ±α` — the Steffensen denominator
    is nonzero, so `steffensen_step` takes its non-identity branch. -/
theorem steffensen_denom_p2_ne_zero
    (x α : ℝ) (hx : x ≠ 0) (hxα : x * α ^ 2 = 1)
    (s : ℝ) (hs : s + 1 ≠ 0)
    (hxs1 : x * s + 1 ≠ 0)
    (hq : 2 * x * s + x + 1 ≠ 0)
    (hsα : s ≠ α) (hs_negα : s ≠ -α) :
    steffensen_denom x 2 s ≠ 0 := by
  rw [steffensen_denom_p2_closed_form x hx s hs hq]
  apply div_ne_zero
  · have hxs2_m1_ne : x * s ^ 2 - 1 ≠ 0 := by
      intro h
      have h_factor : x * (s - α) * (s + α) = 0 := by
        linear_combination h - hxα
      rcases mul_eq_zero.mp h_factor with h1 | h2
      · rcases mul_eq_zero.mp h1 with h3 | h4
        · exact hx h3
        · exact hsα (by linarith)
      · exact hs_negα (by linarith)
    exact mul_ne_zero (mul_ne_zero two_ne_zero hxs1) hxs2_m1_ne
  · exact mul_ne_zero (mul_ne_zero hx hs) hq

/-- **Universal (`α`-free) closed form of `steffensen_step x 2 s`.**
    On the non-degenerate set `s ≠ −1`, `xs + 1 ≠ 0`, `2xs+x+1 ≠ 0`,
    `xs² − 1 ≠ 0`, the generic Steffensen iterator admits the explicit
    rational form
        `steffensen_step x 2 s
           = s − (x·s² − 1)·(2·x·s + x + 1) / (2·x·(s+1)·(x·s + 1))`.
    Pure algebra — no fixed-point relation invoked. -/
theorem steffensen_step_p2_universal
    (x : ℝ) (hx : x ≠ 0)
    (s : ℝ) (hs : s + 1 ≠ 0) (hxs1 : x * s + 1 ≠ 0)
    (hq : 2 * x * s + x + 1 ≠ 0)
    (hxs2_m1_ne : x * s ^ 2 - 1 ≠ 0) :
    steffensen_step x 2 s
      = s - (x * s ^ 2 - 1) * (2 * x * s + x + 1)
              / (2 * x * (s + 1) * (x * s + 1)) := by
  have hxs1_ne : x * (s + 1) ≠ 0 := mul_ne_zero hx hs
  have hD_ne : steffensen_denom x 2 s ≠ 0 := by
    rw [steffensen_denom_p2_closed_form x hx s hs hq]
    exact div_ne_zero
      (mul_ne_zero (mul_ne_zero two_ne_zero hxs1) hxs2_m1_ne)
      (mul_ne_zero hxs1_ne hq)
  unfold steffensen_step
  rw [if_neg hD_ne]
  rw [pandrosion_h_p2_minus_s x hx s hs,
      steffensen_denom_p2_closed_form x hx s hs hq]
  field_simp
  ring

/-- **Explicit `σ_{2,x}` as the universal closed form, under `xα² = 1`.**
    On the same non-degenerate set, the α-parametrised explicit form
    `sigma_p2_explicit x α s` reduces to the α-free universal form,
    via the Pandrosion fixed-point relation `xα² = 1`. -/
theorem sigma_p2_explicit_eq_universal
    (x α : ℝ) (hx : x ≠ 0) (hxα : x * α ^ 2 = 1)
    (s : ℝ) (hs : s + 1 ≠ 0) (hxs1 : x * s + 1 ≠ 0) :
    sigma_p2_explicit x α s
      = s - (x * s ^ 2 - 1) * (2 * x * s + x + 1)
              / (2 * x * (s + 1) * (x * s + 1)) := by
  unfold sigma_p2_explicit
  field_simp
  linear_combination
    2 * (s + 1) * (x * s + 1) * (2 * x * s + x + 1) * hxα

/-- **★★★ Phase A bridge: explicit `σ_{2,x}` coincides with the
    generic Pandrosion-Steffensen iterator at `p = 2`.**

    Under the standard non-degeneracy hypotheses — `x ≠ 0`, fixed-point
    relation `x·α² = 1`, and `s` outside the finite bad set
    `{−1, −1/x, −(x+1)/(2x), ±α}` — the explicit Möbius-closed-form
    iterator `sigma_p2_explicit` of §36.1 **equals** the generic
    `steffensen_step x 2` defined in Foundations §3:

        `steffensen_step x 2 s = sigma_p2_explicit x α s`.

    Together with §36.3, this promotes the cleared Möbius conjugacy
    identity from the explicit formula to the generic Steffensen
    accelerator — removing the last algebraic gap between the §36
    abstract algebraic heart and the downstream `RealMcMullenP2`
    pipeline. -/
theorem steffensen_step_p2_eq_sigma_p2_explicit
    (x α : ℝ) (hx : x ≠ 0) (hxα : x * α ^ 2 = 1)
    (s : ℝ) (hs : s + 1 ≠ 0) (hxs1 : x * s + 1 ≠ 0)
    (hq : 2 * x * s + x + 1 ≠ 0)
    (hsα : s ≠ α) (hs_negα : s ≠ -α) :
    steffensen_step x 2 s = sigma_p2_explicit x α s := by
  have hxs2_m1_ne : x * s ^ 2 - 1 ≠ 0 := by
    intro h
    have h_factor : x * (s - α) * (s + α) = 0 := by
      linear_combination h - hxα
    rcases mul_eq_zero.mp h_factor with h1 | h2
    · rcases mul_eq_zero.mp h1 with h3 | h4
      · exact hx h3
      · exact hsα (by linarith)
    · exact hs_negα (by linarith)
  rw [steffensen_step_p2_universal x hx s hs hxs1 hq hxs2_m1_ne,
      ← sigma_p2_explicit_eq_universal x α hx hxα s hs hxs1]


end BridgeToGeneric

end Pandrosion
