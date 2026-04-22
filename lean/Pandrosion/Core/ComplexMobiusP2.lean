/-
  Universitas Pandrosion — §38. **Möbius conjugacy on ℂ at `p = 2`.**

  Complex-variable port of §36–§37.2: the algebraic core of the
  unconditional `McMullenAEEntry 2 x α` discharge over ℂ.

  Every identity in this module is a **pure polynomial identity**
  modulo `x · α² = 1`, valid verbatim on ℂ — there is no analytic
  obstruction, so the proofs are direct transcripts of their §36–§37
  ℝ counterparts with `ℝ → ℂ` and `|·| → ‖·‖` substitutions.

  Contents.

    §38.1  `sigma_p2_explicit_C` — explicit Möbius form of the
           Steffensen `p = 2` step, defined directly on ℂ.

    §38.2  `sigma_p2_minus_alpha_cleared_C`,
           `sigma_p2_plus_alpha_cleared_C` — cleared-denominator
           factored identities (`x · α² = 1`).

    §38.3  `steffensen_mobius_conjugacy_p2_cleared_C` — the headline
           identity, polynomial form:
               `(σ(z) − α)·(z + α)²·(1 + α)²
                  = (σ(z) + α)·(z − α)²·(1 − α)²`.

    §38.4  `v_p2_C` — Böttcher-normalised Möbius coordinate on ℂ:
           `v_α(z) = μ²·(z − α)/(z + α)`,  `μ := (1 − α)/(1 + α)`.

    §38.5  `v_p2_sq_from_conjugacy_C` — Böttcher functional equation
           `v(σ(z)) = v(z)²`.

    §38.6  `v_p2_iterated_C` — iterated form `v(σⁿ(z)) = v(z)^{2ⁿ}`.

  This module is self-contained and unconditionally proved; the
  measure-theoretic finale (`mcmullen_p2_complex_unconditional`) is
  built on top in §39.
-/

import Pandrosion.Core.Foundations

namespace Pandrosion

open Complex

/-! ============================================================
  §38.1  Explicit complex Möbius form of `σ_{2,x}`
============================================================ -/

section ExplicitSigmaC

/-- **Explicit Möbius-form Steffensen step at `p = 2` on ℂ.**

    On the non-degenerate set `{z : (z + 1)(xz + 1) ≠ 0}`, the
    Pandrosion-Steffensen accelerator admits the closed form
        `σ_{2,x}(z) = z − (z − α)(z + α)(2xz + x + 1) / [2(z+1)(xz+1)]`.

    The bridge to the generic `steffensen_step_C x 2` is established
    through the same `Sp_C 2 = 1 + s` and Möbius-square route used
    on ℝ in §36.4; that bridge is restated in §39 as needed. -/
noncomputable def sigma_p2_explicit_C (x α z : ℂ) : ℂ :=
  z - (z - α) * (z + α) * (2 * x * z + x + 1) / (2 * (z + 1) * (x * z + 1))

end ExplicitSigmaC

/-! ============================================================
  §38.2  Cleared-denominator factored identities
============================================================ -/

section FactoredFormsClearedC

/-- **Cleared-denominator factored form of `σ(z) − α`.**

    At any Pandrosion fixed point (`xα² = 1`), on the non-degenerate
    set `(z + 1)(xz + 1) ≠ 0`,
        `(σ(z) − α) · 2α²(z+1)(xz+1)  =  (z − α)² (1 − α)²`. -/
theorem sigma_p2_minus_alpha_cleared_C
    (x α z : ℂ)
    (hz_ne_neg_one : z + 1 ≠ 0)
    (hz_xz1_ne : x * z + 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    (sigma_p2_explicit_C x α z - α) * (2 * α ^ 2 * (z + 1) * (x * z + 1))
      = (z - α) ^ 2 * (1 - α) ^ 2 := by
  unfold sigma_p2_explicit_C
  have h_denom_ne : (2 * (z + 1) * (x * z + 1) : ℂ) ≠ 0 :=
    mul_ne_zero (mul_ne_zero two_ne_zero hz_ne_neg_one) hz_xz1_ne
  have h_rewrite :
      (z - (z - α) * (z + α) * (2 * x * z + x + 1) / (2 * (z + 1) * (x * z + 1)) - α)
        * (2 * α ^ 2 * (z + 1) * (x * z + 1))
      = (z - α) * (2 * α ^ 2 * (z + 1) * (x * z + 1))
        - α ^ 2 * ((z - α) * (z + α) * (2 * x * z + x + 1)) := by
    field_simp
    ring
  rw [h_rewrite]
  linear_combination ((1 - 2 * α) * z ^ 2 - 2 * α * (1 - α) * z + α ^ 2) * hxα

/-- **Cleared-denominator factored form of `σ(z) + α`.**

    At any Pandrosion fixed point (`xα² = 1`), on the non-degenerate
    set `(z + 1)(xz + 1) ≠ 0`,
        `(σ(z) + α) · 2α²(z+1)(xz+1)  =  (z + α)² (1 + α)²`. -/
theorem sigma_p2_plus_alpha_cleared_C
    (x α z : ℂ)
    (hz_ne_neg_one : z + 1 ≠ 0)
    (hz_xz1_ne : x * z + 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    (sigma_p2_explicit_C x α z + α) * (2 * α ^ 2 * (z + 1) * (x * z + 1))
      = (z + α) ^ 2 * (1 + α) ^ 2 := by
  unfold sigma_p2_explicit_C
  have h_denom_ne : (2 * (z + 1) * (x * z + 1) : ℂ) ≠ 0 :=
    mul_ne_zero (mul_ne_zero two_ne_zero hz_ne_neg_one) hz_xz1_ne
  have h_rewrite :
      (z - (z - α) * (z + α) * (2 * x * z + x + 1) / (2 * (z + 1) * (x * z + 1)) + α)
        * (2 * α ^ 2 * (z + 1) * (x * z + 1))
      = (z + α) * (2 * α ^ 2 * (z + 1) * (x * z + 1))
        - α ^ 2 * ((z - α) * (z + α) * (2 * x * z + x + 1)) := by
    field_simp
    ring
  rw [h_rewrite]
  linear_combination ((1 + 2 * α) * z ^ 2 + 2 * α * (1 + α) * z + α ^ 2) * hxα

end FactoredFormsClearedC

/-! ============================================================
  §38.3  Headline complex Möbius conjugacy (cleared form)
============================================================ -/

section MobiusConjugacyC

/-- **★★★ Möbius conjugacy of Pandrosion-Steffensen on ℂ at `p = 2`.**

    For `xα² = 1`, `α ≠ 0`, and `z` in the non-degenerate set
    `(z + 1)(xz + 1) ≠ 0`,
        `(σ(z) − α) · (z + α)² · (1 + α)²
          = (σ(z) + α) · (z − α)² · (1 − α)²`.

    Equivalently (on the complement of `{z : z = ±α}` and
    `{α = ±1}`):
        `(σ(z) − α)/(σ(z) + α) = μ² · ((z − α)/(z + α))²`,
    `μ := (1 − α)/(1 + α)`.

    This is the **complex analytic golden key** for unconditional
    `McMullenAEEntry 2 x α`: `σ_{2,x}` is conjugate on ℂP¹ to the
    scaled squaring map `w ↦ μ²w²`, whose Julia set is the circle
    `|μ²w| = 1`, Lebesgue-null in 2 dimensions. -/
theorem steffensen_mobius_conjugacy_p2_cleared_C
    (x α z : ℂ)
    (hα_ne_zero : α ≠ 0)
    (hz_ne_neg_one : z + 1 ≠ 0)
    (hz_xz1_ne : x * z + 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    (sigma_p2_explicit_C x α z - α) * (z + α) ^ 2 * (1 + α) ^ 2
      = (sigma_p2_explicit_C x α z + α) * (z - α) ^ 2 * (1 - α) ^ 2 := by
  have h_minus := sigma_p2_minus_alpha_cleared_C x α z hz_ne_neg_one hz_xz1_ne hxα
  have h_plus  := sigma_p2_plus_alpha_cleared_C  x α z hz_ne_neg_one hz_xz1_ne hxα
  have hα2_ne : α ^ 2 ≠ 0 := pow_ne_zero 2 hα_ne_zero
  have hK_ne : (2 * α ^ 2 * (z + 1) * (x * z + 1) : ℂ) ≠ 0 :=
    mul_ne_zero
      (mul_ne_zero (mul_ne_zero two_ne_zero hα2_ne) hz_ne_neg_one)
      hz_xz1_ne
  have h_key :
      (sigma_p2_explicit_C x α z - α) * (z + α) ^ 2 * (1 + α) ^ 2
        * (2 * α ^ 2 * (z + 1) * (x * z + 1))
      = (sigma_p2_explicit_C x α z + α) * (z - α) ^ 2 * (1 - α) ^ 2
        * (2 * α ^ 2 * (z + 1) * (x * z + 1)) := by
    linear_combination (z + α) ^ 2 * (1 + α) ^ 2 * h_minus
                       - (z - α) ^ 2 * (1 - α) ^ 2 * h_plus
  exact mul_right_cancel₀ hK_ne h_key

end MobiusConjugacyC

/-! ============================================================
  §38.4  Böttcher-normalised coordinate `v(z) = μ²·φ(z)` on ℂ
============================================================ -/

section BottcherCoordinateC

/-- **Böttcher-normalised Möbius coordinate on ℂ at `p = 2`.**
    `v_{α}(z) = μ²·(z − α)/(z + α)`, where `μ = (1−α)/(1+α)`.
    Satisfies `v(σ(z)) = v(z)²` (§38.5), the conjugacy `σ ∼ squaring`
    on `ℂP¹`. -/
noncomputable def v_p2_C (α z : ℂ) : ℂ :=
  ((1 - α) / (1 + α)) ^ 2 * ((z - α) / (z + α))

/-- **Böttcher functional equation `v(σ(z)) = v(z)²` on ℂ.**

    For `z` in the non-degenerate set
    `(z + α) ≠ 0, (1 + α) ≠ 0, (z + 1) ≠ 0, (xz + 1) ≠ 0,
     σ(z) + α ≠ 0`, and at any Pandrosion fixed point `xα² = 1`,
    the Möbius conjugacy §38.3 in ratio form gives
        `v(σ(z)) = v(z)²`. -/
theorem v_p2_sq_from_conjugacy_C
    (x α z : ℂ)
    (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0)
    (hz_ne_neg_one : z + 1 ≠ 0)
    (hz_xz1_ne : x * z + 1 ≠ 0)
    (hz_plus_α_ne : z + α ≠ 0)
    (hsigma_plus_α_ne : sigma_p2_explicit_C x α z + α ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    v_p2_C α (sigma_p2_explicit_C x α z) = (v_p2_C α z) ^ 2 := by
  unfold v_p2_C
  have h_conj := steffensen_mobius_conjugacy_p2_cleared_C
    x α z hα_ne_zero hz_ne_neg_one hz_xz1_ne hxα
  have h_ratio : (sigma_p2_explicit_C x α z - α) / (sigma_p2_explicit_C x α z + α)
               = ((1 - α) / (1 + α)) ^ 2 * ((z - α) / (z + α)) ^ 2 := by
    rw [div_eq_iff hsigma_plus_α_ne]
    field_simp
    linear_combination h_conj
  rw [h_ratio]
  ring

end BottcherCoordinateC

/-! ============================================================
  §38.5  Iterated complex Böttcher conjugacy
============================================================ -/

section IteratedConjugacyC

/-- **★★★ Iterated Böttcher conjugacy on ℂ at `p = 2`.**

    For every `n ≥ 0` and every `z` whose forward iterates under
    `sigma_p2_explicit_C x α` stay inside the non-degenerate set
    `{w + 1 ≠ 0, x·w + 1 ≠ 0, w + α ≠ 0}`,
        `v(σⁿ(z)) = v(z)^{2ⁿ}`.

    **Algebraic skeleton of the unconditional `McMullenAEEntry 2 x α`
    on ℂ.** Under the squaring recurrence `vₙ₊₁ = vₙ²`:
      • if `‖v(z)‖ < 1`, then `v(σⁿ(z)) → 0`, so `σⁿ(z) → α`;
      • if `‖v(z)‖ > 1`, then `v(σⁿ(z)) → ∞`, so `σⁿ(z) → −α`;
      • `‖v(z)‖ = 1` is a Lebesgue-null circle on ℂ (§39). -/
theorem v_p2_iterated_C
    (x α : ℂ) (hα_ne_zero : α ≠ 0) (hα_ne_neg_one : 1 + α ≠ 0)
    (hxα : x * α ^ 2 = 1)
    (z : ℂ) (n : ℕ)
    (H : ∀ k ≤ n, (sigma_p2_explicit_C x α)^[k] z + 1 ≠ 0
                ∧ x * (sigma_p2_explicit_C x α)^[k] z + 1 ≠ 0
                ∧ (sigma_p2_explicit_C x α)^[k] z + α ≠ 0) :
    v_p2_C α ((sigma_p2_explicit_C x α)^[n] z) = (v_p2_C α z) ^ (2 ^ n) := by
  induction n with
  | zero => simp [pow_zero, Function.iterate_zero_apply]
  | succ k ih =>
    have H_k : ∀ j ≤ k, (sigma_p2_explicit_C x α)^[j] z + 1 ≠ 0
                      ∧ x * (sigma_p2_explicit_C x α)^[j] z + 1 ≠ 0
                      ∧ (sigma_p2_explicit_C x α)^[j] z + α ≠ 0 := by
      intro j hj; exact H j (Nat.le_succ_of_le hj)
    obtain ⟨h_ne_m1, h_xz1, h_plus_α⟩ := H k (Nat.le_succ k)
    obtain ⟨_, _, h_plus_α'⟩ := H (k + 1) le_rfl
    have h_iter_succ :
        (sigma_p2_explicit_C x α)^[k + 1] z
          = sigma_p2_explicit_C x α ((sigma_p2_explicit_C x α)^[k] z) := by
      rw [Function.iterate_succ_apply']
    rw [h_iter_succ]
    have h_step := v_p2_sq_from_conjugacy_C x α
      ((sigma_p2_explicit_C x α)^[k] z)
      hα_ne_zero hα_ne_neg_one h_ne_m1 h_xz1 h_plus_α
      (by rw [← h_iter_succ]; exact h_plus_α') hxα
    rw [h_step, ih H_k, ← pow_mul, pow_succ]

end IteratedConjugacyC

end Pandrosion
