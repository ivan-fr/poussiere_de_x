/-
  Universitas Pandrosion — §45. **Measure-uniform effective
  complexity for `p = 2` on ℂ.**

  Synthesis of the §40 (pointwise effective rate) and §43
  (slow-set decay in measure) results into a *measure-uniform*
  complexity theorem:

      *For every `δ > 0` and every `r > 0`, there exists an
      explicit iteration count `K(z₀, r)` for every `z₀` outside
      a bad set of Lebesgue measure `< δ`, such that every
      Pandrosion-Steffensen orbit enters the cyclotomic basin
      of radius `r` within `K(z₀, r)` iterations.*

  The bad set is `slowSet α (1/(N+1)) ∪ complex_bad_set ∪
  alphaHitSet`, each component controlled:
    • `slowSet α (1/(N+1)) ∩ B(0, R)`: volume `< δ` for `N` large
      enough (§43.7).
    • `complex_bad_set x α`: Lebesgue-null (§39.5).
    • `alphaHitSet x α`: countable union of finite `σⁿ`-preimages
      of `{α, −α}` — Lebesgue-null.

  Off this combined set, §40.6
  (`mcmullen_p2_complex_effective`) yields basin-entry with the
  explicit bound `max(entryTimeAlpha α r ‖v(z₀)‖,
  entryTimeNegAlpha α r ‖v(z₀)‖)`.

  Contents.

    §45.1  `alphaHitSet` and its Lebesgue-nullity.

    §45.2  `slowSet_iff_norm_v` — equivalence between the
           polynomial §43 `slowSet` and the `|‖v‖ − 1| ≤ ε` form
           (modulo `z + α ≠ 0`).

    §45.3  `mcmullen_p2_complex_measure_uniform` — headline
           theorem: measure bound + per-`z₀` basin-entry with
           closed-form `max(entryTimeAlpha, entryTimeNegAlpha)`.
-/

import Pandrosion.Core.ComplexMcMullenP2EffectiveRate
import Pandrosion.Core.ComplexJuliaSlowSetDecay
import Pandrosion.Core.QuadraticComplexityEffective

namespace Pandrosion

open MeasureTheory Complex Filter Topology

/-! ============================================================
  §45.1  `alphaHitSet` is Lebesgue-null
============================================================ -/

section AlphaHitSetC

/-- **`±α`-hit set.**
    Points `z₀` whose forward orbit under `steffensen_step_C x 2`
    eventually equals `α` or `−α`. Countable union of finite
    `σⁿ`-fibers of `{α, −α}` (§39.4), hence Lebesgue-null. -/
def alphaHitSet (x α : ℂ) : Set ℂ :=
  {z₀ : ℂ | ∃ k : ℕ, (steffensen_step_C x 2)^[k] z₀ = α} ∪
  {z₀ : ℂ | ∃ k : ℕ, (steffensen_step_C x 2)^[k] z₀ = -α}

theorem alphaHitSet_volume_zero
    (x α : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    volume (alphaHitSet x α) = 0 := by
  refine measure_union_null ?_ ?_
  · have h_subset :
        {z₀ : ℂ | ∃ k : ℕ, (steffensen_step_C x 2)^[k] z₀ = α} ⊆
        ⋃ k : ℕ, {z₀ : ℂ | (steffensen_step_C x 2)^[k] z₀ = α} := by
      intro z₀ hz₀
      obtain ⟨k, hk⟩ := hz₀
      rw [Set.mem_iUnion]; exact ⟨k, hk⟩
    refine measure_mono_null h_subset ?_
    have h_each_finite : ∀ k : ℕ,
        Set.Finite {z₀ : ℂ | (steffensen_step_C x 2)^[k] z₀ = α} := fun k =>
      steffensen_step_C_iterate_fiber_finite x α hx_ne hx_ne_one hxα k α
    exact (Set.countable_iUnion (fun k => (h_each_finite k).countable)).measure_zero _
  · have h_subset :
        {z₀ : ℂ | ∃ k : ℕ, (steffensen_step_C x 2)^[k] z₀ = -α} ⊆
        ⋃ k : ℕ, {z₀ : ℂ | (steffensen_step_C x 2)^[k] z₀ = -α} := by
      intro z₀ hz₀
      obtain ⟨k, hk⟩ := hz₀
      rw [Set.mem_iUnion]; exact ⟨k, hk⟩
    refine measure_mono_null h_subset ?_
    have h_each_finite : ∀ k : ℕ,
        Set.Finite {z₀ : ℂ | (steffensen_step_C x 2)^[k] z₀ = -α} := fun k =>
      steffensen_step_C_iterate_fiber_finite x α hx_ne hx_ne_one hxα k (-α)
    exact (Set.countable_iUnion (fun k => (h_each_finite k).countable)).measure_zero _

end AlphaHitSetC

/-! ============================================================
  §45.2  Polynomial ↔ `|‖v‖ − 1|` equivalence for `slowSet`
============================================================ -/

section SlowSetEquivC

/-- **Polynomial `slowSet` ↔ `|‖v‖ − 1| ≤ ε` form.**

    For `0 ≤ ε ≤ 1` and `z + α ≠ 0`, the §43 polynomial
    `slowSet α ε` membership is equivalent to the Böttcher-
    coordinate-ring-form `|‖v_p2_C α z‖ − 1| ≤ ε`. -/
theorem slowSet_iff_norm_v
    (α : ℂ) (ε : ℝ) (hε_le : ε ≤ 1)
    (z : ℂ) (hz_α : z + α ≠ 0) :
    z ∈ slowSet α ε ↔ |‖v_p2_C α z‖ - 1| ≤ ε := by
  constructor
  · intro hz
    obtain ⟨h1, h2⟩ := hz
    -- Divide both by ‖z + α‖² to get (1 - ε)² ≤ ‖v(z)‖² ≤ (1 + ε)².
    have hzα_pos : 0 < ‖z + α‖ ^ 2 := by
      have : 0 < ‖z + α‖ := norm_pos_iff.mpr hz_α
      positivity
    have h_v_sq_eq := v_p2_C_norm_sq_eq α z
    have h_lb : (1 - ε) ^ 2 ≤ ‖v_p2_C α z‖ ^ 2 := by
      have h_div : (1 - ε) ^ 2 ≤
          ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 / ‖z + α‖ ^ 2 := by
        rw [le_div_iff hzα_pos]; exact h1
      rw [← h_v_sq_eq] at h_div; exact h_div
    have h_ub : ‖v_p2_C α z‖ ^ 2 ≤ (1 + ε) ^ 2 := by
      have h_div : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 / ‖z + α‖ ^ 2
                   ≤ (1 + ε) ^ 2 := by
        rw [div_le_iff hzα_pos]; exact h2
      rw [← h_v_sq_eq] at h_div; exact h_div
    -- From (1 - ε)² ≤ ‖v‖² ≤ (1 + ε)² deduce |‖v‖ - 1| ≤ ε.
    have h_v_nn : 0 ≤ ‖v_p2_C α z‖ := norm_nonneg _
    have h_one_minus_nn : 0 ≤ 1 - ε := by linarith
    have h_one_plus_nn : 0 ≤ 1 + ε := by linarith
    -- ‖v‖ ≤ 1 + ε.
    have h_v_le : ‖v_p2_C α z‖ ≤ 1 + ε := by
      have h_sq : ‖v_p2_C α z‖ ^ 2 ≤ (1 + ε) ^ 2 := h_ub
      exact (pow_le_pow_iff_left h_v_nn h_one_plus_nn (by norm_num)).mp h_sq
    -- ‖v‖ ≥ 1 - ε.
    have h_v_ge : 1 - ε ≤ ‖v_p2_C α z‖ := by
      have h_sq : (1 - ε) ^ 2 ≤ ‖v_p2_C α z‖ ^ 2 := h_lb
      exact (pow_le_pow_iff_left h_one_minus_nn h_v_nn (by norm_num)).mp h_sq
    -- Combine to get |‖v‖ - 1| ≤ ε.
    rw [abs_sub_le_iff]
    refine ⟨?_, ?_⟩ <;> linarith
  · intro hz
    have hz_le : |‖v_p2_C α z‖ - 1| ≤ ε := hz
    have hv_range : 1 - ε ≤ ‖v_p2_C α z‖ ∧ ‖v_p2_C α z‖ ≤ 1 + ε := by
      rw [abs_sub_le_iff] at hz_le
      refine ⟨?_, ?_⟩ <;> linarith [hz_le.1, hz_le.2]
    have hzα_pos : 0 < ‖z + α‖ ^ 2 := by
      have : 0 < ‖z + α‖ := norm_pos_iff.mpr hz_α
      positivity
    have h_v_sq_eq := v_p2_C_norm_sq_eq α z
    have h_v_nn : 0 ≤ ‖v_p2_C α z‖ := norm_nonneg _
    have h_one_minus_nn : 0 ≤ 1 - ε := by linarith
    have h_v_sq_le : ‖v_p2_C α z‖ ^ 2 ≤ (1 + ε) ^ 2 :=
      pow_le_pow_left h_v_nn hv_range.2 2
    have h_v_sq_ge : (1 - ε) ^ 2 ≤ ‖v_p2_C α z‖ ^ 2 :=
      pow_le_pow_left h_one_minus_nn hv_range.1 2
    refine ⟨?_, ?_⟩
    · rw [← le_div_iff hzα_pos, ← h_v_sq_eq]; exact h_v_sq_ge
    · rw [← div_le_iff hzα_pos, ← h_v_sq_eq]; exact h_v_sq_le

end SlowSetEquivC

/-! ============================================================
  §45.3  Measure-uniform complexity theorem (★★★★★)
============================================================ -/

section MeasureUniformComplexityC

/-- **★★★★★ Measure-uniform effective complexity for `p = 2` on ℂ.**

    For every `x ∈ ℂ \ {0, 1}` (with `α ≠ 0, ±1` and `α² = 1/x`),
    every `R > 0`, every `δ > 0`, and every `r > 0`, there exist
    `N : ℕ` such that:

      (i)  The "bad set" `slowSet α (1/(N+1)) ∩ B(0,R) ∪
           alphaHitSet ∪ complex_bad_set` has volume `< δ`.

      (ii) For every `z₀` outside the bad set, basin entry into
           `B(cycAnchor α 2 s, r)` holds with the explicit bound
           `K(z₀) := max(entryTimeAlpha α r ‖v(z₀)‖,
                         entryTimeNegAlpha α r ‖v(z₀)‖)`,
           for some `s ∈ Fin 2`.

    *Interpretation.* Basin-entry within a **computable, closed-
    form** iteration count, for Lebesgue-almost-everywhere except
    a prescribed `δ`-mass of `B(0, R)`. Combines §40.6 (pointwise
    effective rate) with §43 (slow-set decay in measure) and
    §39.5 (complex bad set null) and §45.1 (α-hit null). -/
theorem mcmullen_p2_complex_measure_uniform
    (x : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (R : ℝ) (δ : ℝ) (hδ_pos : 0 < δ) (r : ℝ) (hr_pos : 0 < r) :
    ∃ N : ℕ,
      volume ((slowSet α (1 / ((N : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R)
              ∪ alphaHitSet x α
              ∪ complex_bad_set x α) < ENNReal.ofReal δ ∧
      ∀ z₀ : ℂ,
        z₀ ∉ slowSet α (1 / ((N : ℝ) + 1)) →
        z₀ ∉ alphaHitSet x α →
        z₀ ∉ complex_bad_set x α →
        0 < ‖v_p2_C α z₀‖ →
        ∃ s : Fin 2, ∀ k ≥ max (entryTimeAlpha α r ‖v_p2_C α z₀‖)
                             (entryTimeNegAlpha α r ‖v_p2_C α z₀‖),
          ‖(steffensen_step_C x 2)^[k] z₀ - cycAnchor α 2 s‖ < r := by
  have hxα : x * α ^ 2 = 1 := by rw [hα_pow]; field_simp
  obtain ⟨N, hN⟩ := slow_set_decay_qualitative α hα_ne_zero R δ hδ_pos
  refine ⟨N, ?_, ?_⟩
  · -- Measure bound.
    have h_slow := hN N (le_refl N)
    have h_α_hit_zero : volume (alphaHitSet x α) = 0 :=
      alphaHitSet_volume_zero x α hx_ne hx_ne_one hxα
    have h_bad_zero : volume (complex_bad_set x α) = 0 :=
      complex_bad_set_volume_zero x α hx_ne hx_ne_one hα_ne_zero hxα
    calc volume ((slowSet α (1 / ((N : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R)
                  ∪ alphaHitSet x α
                  ∪ complex_bad_set x α)
        ≤ volume ((slowSet α (1 / ((N : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R)
                   ∪ alphaHitSet x α) +
          volume (complex_bad_set x α) :=
          measure_union_le _ _
      _ ≤ volume (slowSet α (1 / ((N : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R) +
          volume (alphaHitSet x α) + volume (complex_bad_set x α) := by
          gcongr; exact measure_union_le _ _
      _ = volume (slowSet α (1 / ((N : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R) := by
          rw [h_α_hit_zero, h_bad_zero, add_zero, add_zero]
      _ < ENNReal.ofReal δ := h_slow
  · -- Per-z₀ basin entry, direct chain through §40.6.
    intro z₀ _h_not_slow h_not_α h_not_bad hv_pos
    have h_no_alpha_hit : ∀ k, (steffensen_step_C x 2)^[k] z₀ ≠ α := by
      intro k h_eq; apply h_not_α; left; exact ⟨k, h_eq⟩
    have h_no_neg_alpha_hit : ∀ k, (steffensen_step_C x 2)^[k] z₀ ≠ -α := by
      intro k h_eq; apply h_not_α; right; exact ⟨k, h_eq⟩
    exact mcmullen_p2_complex_effective x hx_ne α hα_ne_zero hα_ne_neg_one
      hα_ne_one hα_pow z₀ h_not_bad hv_pos h_no_alpha_hit h_no_neg_alpha_hit
      r hr_pos

end MeasureUniformComplexityC

end Pandrosion
