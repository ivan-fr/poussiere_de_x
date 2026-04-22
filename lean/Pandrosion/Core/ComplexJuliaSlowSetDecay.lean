/-
  Universitas Pandrosion — §43. **Slow-set Lebesgue decay around the
  Julia circle (`p = 2` complex).**

  Quantitative refinement in measure of §39: the slow set
      `slowSet α ε := {z ∈ ℂ : (1−ε)²·‖z+α‖² ≤ ‖μ²‖²·‖z−α‖²
                              ≤ (1+ε)²·‖z+α‖²}`
  shrinks (in volume) to the polynomial Julia locus
      `juliaSetP α := {z : ‖μ²‖²·‖z−α‖² = ‖z+α‖²}`
  which has Lebesgue measure 0 (modulo `{−α}`, identified with
  §39.2's `julia_section_complex_volume_zero`).

  **Polynomial formulation rationale.** Defining `slowSet` directly
  in terms of `‖v_p2_C α z‖` would entangle the Lean
  measurability story with `Inv.inv : ℂ → ℂ` (continuous on `ℂˣ`,
  discontinuous at `0` under the `0⁻¹ = 0` convention). The
  polynomial form sidesteps this entirely: every set involved is
  cut out by continuous polynomial inequalities in `(Re z, Im z)`,
  hence trivially measurable.

  The polynomial form coincides with the `‖v(z)‖`-form modulo the
  single point `{−α}` (which has Lebesgue measure 0), so the
  measure-zero conclusion of §39.2 transports verbatim.

  Contents.

    §43.1  `juliaSetP`, `slowSet` polynomial definitions.

    §43.2  Measurability — both sets are level/sublevel sets of
           continuous polynomial-in-norm functions.

    §43.3  `slowSet_antitone`, `juliaSetP_subset_slowSet` — basic
           monotonicity in `ε`.

    §43.4  `slowSet_iInter_eq_juliaSetP` — decreasing intersection.

    §43.5  `juliaSetP_volume_zero` — polynomial Julia locus is
           Lebesgue-null (chained with §39.2 via `{−α}` agreement).

    §43.6  `volume_slowSet_inter_ball_tendsto_zero` — the volume
           of `slowSet α ε ∩ B(0, R)` tends to 0 as `ε ↓ 0`.

    §43.7  `slow_set_decay_qualitative` — paper-citable form.
-/

import Pandrosion.Core.ComplexMcMullenP2Unconditional

namespace Pandrosion

open MeasureTheory Complex Filter Topology

/-! ============================================================
  §43.1  Polynomial Julia locus and slow set
============================================================ -/

section JuliaSlowSetDefsC

/-- **Polynomial Julia locus on ℂ at `p = 2`.**
    `juliaSetP α := {z : ‖μ²‖²·‖z−α‖² = ‖z+α‖²}` where
    `μ² = ((1−α)/(1+α))²`. Coincides with §39.2's
    `julia_section_complex_volume_zero` set on `ℂ \ {−α}`. -/
def juliaSetP (α : ℂ) : Set ℂ :=
  {z : ℂ | ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2}

/-- **Slow set (polynomial form):** points where the Apollonius
    Möbius ratio is within `ε` of `1`.

    `slowSet α ε := {z : (1−ε)²·‖z+α‖² ≤ ‖μ²‖²·‖z−α‖² ≤
                          (1+ε)²·‖z+α‖²}`. -/
def slowSet (α : ℂ) (ε : ℝ) : Set ℂ :=
  {z : ℂ | (1 - ε) ^ 2 * ‖z + α‖ ^ 2 ≤
             ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2
         ∧ ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 ≤
             (1 + ε) ^ 2 * ‖z + α‖ ^ 2}

end JuliaSlowSetDefsC

/-! ============================================================
  §43.2  Measurability via continuous polynomial form
============================================================ -/

section MeasurabilityC

/-- The function `z ↦ ‖z − a‖²` is continuous, hence Borel measurable. -/
private theorem norm_sq_sub_const_continuous (a : ℂ) :
    Continuous (fun z : ℂ => ‖z - a‖ ^ 2) :=
  ((continuous_id.sub continuous_const).norm).pow 2

/-- The function `z ↦ ‖z + a‖²` is continuous. -/
private theorem norm_sq_add_const_continuous (a : ℂ) :
    Continuous (fun z : ℂ => ‖z + a‖ ^ 2) :=
  ((continuous_id.add continuous_const).norm).pow 2

/-- **`juliaSetP` measurability.** -/
theorem juliaSetP_measurable (α : ℂ) : MeasurableSet (juliaSetP α) := by
  unfold juliaSetP
  have h_cont :
      Continuous (fun z : ℂ =>
        ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 - ‖z + α‖ ^ 2) :=
    (continuous_const.mul (norm_sq_sub_const_continuous α)).sub
      (norm_sq_add_const_continuous α)
  have h_eq :
      {z : ℂ | ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2}
        = (fun z : ℂ =>
              ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 - ‖z + α‖ ^ 2) ⁻¹' {0} := by
    ext z; simp [Set.mem_setOf_eq, sub_eq_zero]
  rw [h_eq]
  exact h_cont.measurable (measurableSet_singleton 0)

/-- **`slowSet` measurability.** -/
theorem slowSet_measurable (α : ℂ) (ε : ℝ) : MeasurableSet (slowSet α ε) := by
  unfold slowSet
  have h_cont1 :
      Continuous (fun z : ℂ =>
        ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 - (1 - ε) ^ 2 * ‖z + α‖ ^ 2) :=
    (continuous_const.mul (norm_sq_sub_const_continuous α)).sub
      (continuous_const.mul (norm_sq_add_const_continuous α))
  have h_cont2 :
      Continuous (fun z : ℂ =>
        (1 + ε) ^ 2 * ‖z + α‖ ^ 2 - ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2) :=
    (continuous_const.mul (norm_sq_add_const_continuous α)).sub
      (continuous_const.mul (norm_sq_sub_const_continuous α))
  have h_set_eq :
      {z : ℂ |
          (1 - ε) ^ 2 * ‖z + α‖ ^ 2 ≤ ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2
        ∧ ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 ≤ (1 + ε) ^ 2 * ‖z + α‖ ^ 2}
      = ((fun z : ℂ =>
            ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 - (1 - ε) ^ 2 * ‖z + α‖ ^ 2)
              ⁻¹' Set.Ici 0) ∩
        ((fun z : ℂ =>
            (1 + ε) ^ 2 * ‖z + α‖ ^ 2 - ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2)
              ⁻¹' Set.Ici 0) := by
    ext z
    simp only [Set.mem_setOf_eq, Set.mem_inter_iff, Set.mem_preimage, Set.mem_Ici,
               sub_nonneg]
  rw [h_set_eq]
  exact (h_cont1.measurable measurableSet_Ici).inter
        (h_cont2.measurable measurableSet_Ici)

end MeasurabilityC

/-! ============================================================
  §43.3  Monotonicity of `slowSet` in `ε`
============================================================ -/

section SlowSetMonotoneC

theorem juliaSetP_subset_slowSet (α : ℂ) {ε : ℝ} (hε : 0 ≤ ε) (hε_le : ε ≤ 1) :
    juliaSetP α ⊆ slowSet α ε := by
  intro z hz
  show (1 - ε) ^ 2 * ‖z + α‖ ^ 2 ≤
         ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2
       ∧ ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 ≤
         (1 + ε) ^ 2 * ‖z + α‖ ^ 2
  rw [show ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2 from hz]
  have h_zα_nn : (0 : ℝ) ≤ ‖z + α‖ ^ 2 := by positivity
  have h_one_minus_in : 0 ≤ 1 - ε ∧ 1 - ε ≤ 1 := ⟨by linarith, by linarith⟩
  have h_one_minus : (1 - ε) ^ 2 ≤ 1 := by
    calc (1 - ε) ^ 2 ≤ 1 ^ 2 :=
            pow_le_pow_left h_one_minus_in.1 h_one_minus_in.2 2
      _ = 1 := by ring
  have h_one_plus : (1 : ℝ) ≤ (1 + ε) ^ 2 := by
    calc (1 : ℝ) = 1 ^ 2 := by ring
      _ ≤ (1 + ε) ^ 2 := pow_le_pow_left zero_le_one (by linarith) 2
  refine ⟨?_, ?_⟩
  · calc (1 - ε) ^ 2 * ‖z + α‖ ^ 2 ≤ 1 * ‖z + α‖ ^ 2 :=
            mul_le_mul_of_nonneg_right h_one_minus h_zα_nn
      _ = ‖z + α‖ ^ 2 := by ring
  · calc ‖z + α‖ ^ 2 = 1 * ‖z + α‖ ^ 2 := by ring
      _ ≤ (1 + ε) ^ 2 * ‖z + α‖ ^ 2 :=
            mul_le_mul_of_nonneg_right h_one_plus h_zα_nn

end SlowSetMonotoneC

/-! ============================================================
  §43.4  Decreasing intersection equals `juliaSetP`
============================================================ -/

section SlowSetLimitC

/-- **Antitone family `n ↦ slowSet α (1/(n+1))`.** -/
theorem slowSet_seq_antitone (α : ℂ) :
    Antitone (fun n : ℕ => slowSet α (1 / ((n : ℝ) + 1))) := by
  intro n m hnm z hz
  -- ε_m = 1/(m+1) ≤ 1/(n+1) = ε_n; show z ∈ slowSet α ε_n.
  have h_le_inv : (1 : ℝ) / ((m : ℝ) + 1) ≤ 1 / ((n : ℝ) + 1) := by
    apply div_le_div_of_nonneg_left (by norm_num) (by positivity)
    exact_mod_cast Nat.add_le_add_right hnm 1
  obtain ⟨h1, h2⟩ := hz
  refine ⟨?_, ?_⟩
  · -- (1 - ε_n)² ≤ (1 - ε_m)² when ε_n ≥ ε_m and both ≤ 1.
    have h_ε_n_le_one : 1 / ((n : ℝ) + 1) ≤ 1 := by
      apply div_le_one_of_le
      · linarith [show (0 : ℝ) ≤ (n : ℝ) by positivity]
      · positivity
    have h_step : (1 - 1 / ((n : ℝ) + 1)) ^ 2 ≤ (1 - 1 / ((m : ℝ) + 1)) ^ 2 := by
      have h_a_nn : 0 ≤ 1 - 1 / ((n : ℝ) + 1) := by linarith
      have h_b_ge : 1 - 1 / ((n : ℝ) + 1) ≤ 1 - 1 / ((m : ℝ) + 1) := by linarith
      nlinarith [h_a_nn, h_b_ge]
    have h_zα_nn : (0 : ℝ) ≤ ‖z + α‖ ^ 2 := by positivity
    have h_step_mul :
        (1 - 1 / ((n : ℝ) + 1)) ^ 2 * ‖z + α‖ ^ 2 ≤
          (1 - 1 / ((m : ℝ) + 1)) ^ 2 * ‖z + α‖ ^ 2 :=
      mul_le_mul_of_nonneg_right h_step h_zα_nn
    linarith
  · -- (1 + ε_m)² ≤ (1 + ε_n)² when ε_m ≤ ε_n.
    have h_step : (1 + 1 / ((m : ℝ) + 1)) ^ 2 ≤ (1 + 1 / ((n : ℝ) + 1)) ^ 2 := by
      have h_a_nn : 0 ≤ 1 + 1 / ((m : ℝ) + 1) := by positivity
      have h_b_ge : 1 + 1 / ((m : ℝ) + 1) ≤ 1 + 1 / ((n : ℝ) + 1) := by linarith
      nlinarith [h_a_nn, h_b_ge]
    have h_zα_nn : (0 : ℝ) ≤ ‖z + α‖ ^ 2 := by positivity
    have h_step_mul :
        (1 + 1 / ((m : ℝ) + 1)) ^ 2 * ‖z + α‖ ^ 2 ≤
          (1 + 1 / ((n : ℝ) + 1)) ^ 2 * ‖z + α‖ ^ 2 :=
      mul_le_mul_of_nonneg_right h_step h_zα_nn
    linarith

/-- **★ Intersection equals `juliaSetP`.**

    `⋂ n, slowSet α (1/(n+1)) = juliaSetP α`. -/
theorem slowSet_iInter_eq_juliaSetP (α : ℂ) :
    (⋂ n : ℕ, slowSet α (1 / ((n : ℝ) + 1))) = juliaSetP α := by
  ext z
  simp only [Set.mem_iInter, slowSet, juliaSetP, Set.mem_setOf_eq]
  set A : ℝ := ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 with hA_def
  set B : ℝ := ‖z + α‖ ^ 2 with hB_def
  have hB_nn : 0 ≤ B := by rw [hB_def]; positivity
  constructor
  · intro h
    -- ∀ n, (1 - 1/(n+1))²·B ≤ A ≤ (1 + 1/(n+1))²·B.  As n → ∞, both sides → B.
    by_contra h_ne
    -- Suppose A ≠ B. Two cases: A > B or A < B.
    rcases lt_or_gt_of_ne h_ne with h_lt | h_gt
    · -- A < B: take n large enough that (1 - 1/(n+1))²·B > A.
      have h_diff_pos : 0 < B - A := by linarith
      -- Pick N such that 1/(N+1) < (B - A)/(2B + 1).
      by_cases hB_zero : B = 0
      · -- B = 0 ⟹ z = -α. Then A = ‖μ²‖²·4‖α‖² ≥ 0; A < 0 impossible.
        rw [hB_zero] at h_lt
        have hA_nn : 0 ≤ A := by rw [hA_def]; positivity
        linarith
      · have hB_pos : 0 < B := lt_of_le_of_ne hB_nn (Ne.symm hB_zero)
        -- Choose ε small enough: (1-ε)²·B > A iff (1-ε)² > A/B (with B > 0).
        -- Equivalently (1-ε) > √(A/B) (when 1-ε ≥ 0, A ≥ 0).
        -- Suffices: ε < 1 - √(A/B).  As h_lt: A < B ⟹ A/B < 1 ⟹ √(A/B) < 1.
        have hA_nn : 0 ≤ A := by rw [hA_def]; positivity
        have h_AB_lt : A / B < 1 := (div_lt_one hB_pos).mpr h_lt
        have h_AB_nn : 0 ≤ A / B := div_nonneg hA_nn (le_of_lt hB_pos)
        set s : ℝ := Real.sqrt (A / B) with hs_def
        have hs_nn : 0 ≤ s := Real.sqrt_nonneg _
        have hs_lt_one : s < 1 := by
          rw [hs_def, show (1 : ℝ) = Real.sqrt 1 from (Real.sqrt_one).symm]
          exact Real.sqrt_lt_sqrt h_AB_nn h_AB_lt
        have h_one_sub_s_pos : 0 < 1 - s := by linarith
        -- Pick n with 1/(n+1) < 1 - s.
        obtain ⟨n, hn⟩ := exists_nat_gt (1 / (1 - s))
        have hn_pos : (0 : ℝ) < n + 1 := by
          have : 0 < (n : ℝ) := lt_of_le_of_lt (by positivity) hn
          linarith
        have h_inv_lt : 1 / ((n : ℝ) + 1) < 1 - s := by
          rw [div_lt_iff hn_pos]
          have h_chain : 1 < ((n : ℝ) + 1) * (1 - s) := by
            calc (1 : ℝ) = (1 / (1 - s)) * (1 - s) := by
                    rw [div_mul_cancel₀]; exact ne_of_gt h_one_sub_s_pos
              _ < n * (1 - s) := mul_lt_mul_of_pos_right hn h_one_sub_s_pos
              _ ≤ (n + 1) * (1 - s) := by
                    apply mul_le_mul_of_nonneg_right
                    · linarith
                    · exact le_of_lt h_one_sub_s_pos
          linarith
        -- Then 1 - 1/(n+1) > s, so (1 - 1/(n+1))² > s² = A/B, so (1-1/(n+1))²·B > A.
        have h_step : s < 1 - 1 / ((n : ℝ) + 1) := by linarith
        have h_step_sq : s ^ 2 < (1 - 1 / ((n : ℝ) + 1)) ^ 2 := by
          have h_lhs_nn : 0 ≤ 1 - 1 / ((n : ℝ) + 1) := by linarith
          nlinarith [hs_nn, h_step, h_lhs_nn]
        have h_s_sq_eq : s ^ 2 = A / B := by
          rw [hs_def]; exact Real.sq_sqrt h_AB_nn
        have h_eq_AB : A < (1 - 1 / ((n : ℝ) + 1)) ^ 2 * B := by
          have h_intermediate : A / B < (1 - 1 / ((n : ℝ) + 1)) ^ 2 := by
            rw [← h_s_sq_eq]; exact h_step_sq
          have h_mul : (A / B) * B < (1 - 1 / ((n : ℝ) + 1)) ^ 2 * B :=
            mul_lt_mul_of_pos_right h_intermediate hB_pos
          rwa [div_mul_cancel₀ A (ne_of_gt hB_pos)] at h_mul
        have h_n := (h n).1
        linarith
    · -- A > B: similar with (1 + 1/(n+1))² < A/B for n large.
      have h_diff_pos : 0 < A - B := by linarith
      have hB_pos : 0 < B := by
        rcases hB_nn.lt_or_eq with hB_lt | hB_eq
        · exact hB_lt
        · exfalso
          have hB_zero : B = 0 := hB_eq.symm
          have h_n2 := (h 0).2
          rw [hB_zero, mul_zero] at h_n2
          linarith
      have h_AB_gt : 1 < A / B := (one_lt_div hB_pos).mpr h_gt
      set s : ℝ := Real.sqrt (A / B) with hs_def
      have hs_pos : 0 < s := Real.sqrt_pos.mpr (by linarith)
      have hs_gt_one : 1 < s := by
        rw [hs_def, show (1 : ℝ) = Real.sqrt 1 from (Real.sqrt_one).symm]
        exact Real.sqrt_lt_sqrt (by linarith) h_AB_gt
      have h_s_minus_one_pos : 0 < s - 1 := by linarith
      obtain ⟨n, hn⟩ := exists_nat_gt (1 / (s - 1))
      have hn_pos : (0 : ℝ) < n + 1 := by
        have : 0 < (n : ℝ) := lt_of_le_of_lt (by positivity) hn
        linarith
      have h_inv_lt : 1 / ((n : ℝ) + 1) < s - 1 := by
        rw [div_lt_iff hn_pos]
        have h_chain : 1 < ((n : ℝ) + 1) * (s - 1) := by
          calc (1 : ℝ) = (1 / (s - 1)) * (s - 1) := by
                  rw [div_mul_cancel₀]; exact ne_of_gt h_s_minus_one_pos
            _ < n * (s - 1) := mul_lt_mul_of_pos_right hn h_s_minus_one_pos
            _ ≤ (n + 1) * (s - 1) := by
                  apply mul_le_mul_of_nonneg_right
                  · linarith
                  · exact le_of_lt h_s_minus_one_pos
        linarith
      have h_step : 1 + 1 / ((n : ℝ) + 1) < s := by linarith
      have h_step_sq : (1 + 1 / ((n : ℝ) + 1)) ^ 2 < s ^ 2 := by
        have h_lhs_nn : 0 ≤ 1 + 1 / ((n : ℝ) + 1) := by positivity
        nlinarith [h_lhs_nn, h_step, hs_pos]
      have h_s_sq_eq : s ^ 2 = A / B := by
        rw [hs_def]; exact Real.sq_sqrt (by linarith)
      have h_eq_AB : (1 + 1 / ((n : ℝ) + 1)) ^ 2 * B < A := by
        have h_intermediate : (1 + 1 / ((n : ℝ) + 1)) ^ 2 < A / B := by
          rw [← h_s_sq_eq]; exact h_step_sq
        have h_mul : (1 + 1 / ((n : ℝ) + 1)) ^ 2 * B < (A / B) * B :=
          mul_lt_mul_of_pos_right h_intermediate hB_pos
        rwa [div_mul_cancel₀ A (ne_of_gt hB_pos)] at h_mul
      have h_n := (h n).2
      linarith
  · -- ←: A = B ⟹ slowSet membership for any ε ≥ 0.
    intro hAB n
    rw [hAB]
    have h_one_div_nn : 0 ≤ 1 / ((n : ℝ) + 1) := by positivity
    have h_one_div_le_one : 1 / ((n : ℝ) + 1) ≤ 1 := by
      apply div_le_one_of_le
      · linarith [show (0 : ℝ) ≤ (n : ℝ) by positivity]
      · positivity
    have h_one_minus_in : 0 ≤ 1 - 1 / ((n : ℝ) + 1) ∧ 1 - 1 / ((n : ℝ) + 1) ≤ 1 := by
      refine ⟨by linarith, by linarith⟩
    have h_one_minus : (1 - 1 / ((n : ℝ) + 1)) ^ 2 ≤ 1 := by
      calc (1 - 1 / ((n : ℝ) + 1)) ^ 2
          ≤ 1 ^ 2 :=
            pow_le_pow_left h_one_minus_in.1 h_one_minus_in.2 2
        _ = 1 := by ring
    have h_one_plus : (1 : ℝ) ≤ (1 + 1 / ((n : ℝ) + 1)) ^ 2 := by
      calc (1 : ℝ) = 1 ^ 2 := by ring
        _ ≤ (1 + 1 / ((n : ℝ) + 1)) ^ 2 :=
            pow_le_pow_left zero_le_one (by linarith) 2
    refine ⟨?_, ?_⟩
    · calc (1 - 1 / ((n : ℝ) + 1)) ^ 2 * B ≤ 1 * B :=
              mul_le_mul_of_nonneg_right h_one_minus hB_nn
        _ = B := by ring
    · calc B = 1 * B := by ring
        _ ≤ (1 + 1 / ((n : ℝ) + 1)) ^ 2 * B :=
              mul_le_mul_of_nonneg_right h_one_plus hB_nn

end SlowSetLimitC

/-! ============================================================
  §43.5  `juliaSetP` has Lebesgue measure zero
============================================================ -/

section JuliaSetPVolumeZeroC

/-- **★★ Polynomial Julia locus is Lebesgue-null.**

    `juliaSetP α` and `julia_section_complex_volume_zero α`'s set
    `{z : ‖v_p2_C α z‖ = 1}` agree on `ℂ \ {−α}`, hence have the
    same Lebesgue volume. By §39.2 the latter is null. -/
theorem juliaSetP_volume_zero (α : ℂ) (hα_ne_zero : α ≠ 0) :
    volume (juliaSetP α) = 0 := by
  -- juliaSetP α ⊆ {z : ‖v_p2_C α z‖ = 1} ∪ {-α}.
  have h_subset : juliaSetP α ⊆ {z : ℂ | ‖v_p2_C α z‖ = 1} ∪ {(-α : ℂ)} := by
    intro z hz
    by_cases h_zα : z + α = 0
    · right; show z = -α; linear_combination h_zα
    · left
      simp only [Set.mem_setOf_eq]
      have h_eq : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2 := hz
      have hzα_pos : 0 < ‖z + α‖ ^ 2 := by
        have : 0 < ‖z + α‖ := norm_pos_iff.mpr h_zα
        positivity
      have h_norm_v_sq : ‖v_p2_C α z‖ ^ 2 = 1 := by
        rw [v_p2_C_norm_sq_eq α z]
        rw [div_eq_one_iff_eq (ne_of_gt hzα_pos)]
        exact h_eq
      have h_norm_nn : 0 ≤ ‖v_p2_C α z‖ := norm_nonneg _
      nlinarith [h_norm_v_sq, h_norm_nn, sq_nonneg (‖v_p2_C α z‖ - 1)]
  refine measure_mono_null h_subset ?_
  refine measure_union_null
    (julia_section_complex_volume_zero α hα_ne_zero) ?_
  exact measure_singleton _

end JuliaSetPVolumeZeroC

/-! ============================================================
  §43.6  Volume of `slowSet ∩ B(0, R)` tends to 0
============================================================ -/

section SlowSetVolumeDecayC

/-- **★★★ Continuity of measure: `volume(slowSet α (1/(n+1)) ∩ B(0, R))
    → 0`.** -/
theorem volume_slowSet_inter_ball_tendsto_zero
    (α : ℂ) (hα_ne_zero : α ≠ 0) (R : ℝ) :
    Tendsto (fun n : ℕ => volume (slowSet α (1 / ((n : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R))
            atTop (nhds 0) := by
  set s : ℕ → Set ℂ := fun n => slowSet α (1 / ((n : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R
  have hs_meas : ∀ n, MeasurableSet (s n) := fun n =>
    (slowSet_measurable α _).inter measurableSet_ball
  have hs_anti : Antitone s := by
    intro n m hnm z hz
    refine ⟨slowSet_seq_antitone α hnm hz.1, hz.2⟩
  have hs0_finite : volume (s 0) ≠ (⊤ : ENNReal) := by
    have h_le : s 0 ⊆ Metric.ball (0 : ℂ) R := Set.inter_subset_right _ _
    have h_ball_bd : Bornology.IsBounded (Metric.ball (0 : ℂ) R) := Metric.isBounded_ball
    have h_ball_lt : volume (Metric.ball (0 : ℂ) R) < (⊤ : ENNReal) :=
      h_ball_bd.measure_lt_top
    exact ne_of_lt (lt_of_le_of_lt (measure_mono h_le) h_ball_lt)
  have h_iInter : (⋂ n : ℕ, s n) = juliaSetP α ∩ Metric.ball (0 : ℂ) R := by
    show (⋂ n : ℕ, slowSet α (1 / ((n : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R)
       = juliaSetP α ∩ Metric.ball (0 : ℂ) R
    have h_split : (⋂ n : ℕ, slowSet α (1 / ((n : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R)
                 = (⋂ n : ℕ, slowSet α (1 / ((n : ℝ) + 1))) ∩ Metric.ball (0 : ℂ) R := by
      ext z
      simp only [Set.mem_iInter, Set.mem_inter_iff]
      constructor
      · intro h
        refine ⟨fun i => (h i).1, (h 0).2⟩
      · intro ⟨h1, h2⟩ i
        exact ⟨h1 i, h2⟩
    rw [h_split, slowSet_iInter_eq_juliaSetP α]
  have h_inter_null : volume (juliaSetP α ∩ Metric.ball (0 : ℂ) R) = 0 := by
    refine measure_mono_null (Set.inter_subset_left _ _) ?_
    exact juliaSetP_volume_zero α hα_ne_zero
  have h_tendsto : Tendsto (volume ∘ s) atTop (nhds (volume (⋂ n, s n))) :=
    tendsto_measure_iInter hs_meas hs_anti ⟨0, hs0_finite⟩
  rw [h_iInter, h_inter_null] at h_tendsto
  exact h_tendsto

end SlowSetVolumeDecayC

/-! ============================================================
  §43.7  Paper-citable qualitative slow-set decay
============================================================ -/

section QualitativeDecayC

/-- **★★★★ Slow-set decay (qualitative form).**

    For every `α ≠ 0`, every ball radius `R > 0`, and every target
    `δ > 0`, there exists `N : ℕ` such that for every `n ≥ N`,
        `volume (slowSet α (1/(n+1)) ∩ B(0, R)) < δ`.

    *Interpretation.* Off a set of arbitrarily small Lebesgue
    measure (the slow set), the §40 effective entry-time bound
    `k₀(z₀) ≤ entryTimeAlpha/Neg α r ‖v(z₀)‖` is *bounded by an
    explicit `K(δ)`*. Combined with §39.2 (`julia_section_complex_volume_zero`)
    and §40 (`mcmullen_p2_complex_effective`), this completes the
    quantitative-in-measure picture: arbitrarily fast convergence
    on a set of arbitrarily large Lebesgue measure. -/
theorem slow_set_decay_qualitative
    (α : ℂ) (hα_ne_zero : α ≠ 0) (R : ℝ) (δ : ℝ) (hδ_pos : 0 < δ) :
    ∃ N : ℕ, ∀ n ≥ N,
      volume (slowSet α (1 / ((n : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R) <
        ENNReal.ofReal δ := by
  have h_tendsto := volume_slowSet_inter_ball_tendsto_zero α hα_ne_zero R
  have h_event : ∀ᶠ n : ℕ in atTop,
      volume (slowSet α (1 / ((n : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R) <
        ENNReal.ofReal δ := by
    have h_mem : Set.Iio (ENNReal.ofReal δ) ∈ nhds (0 : ENNReal) := by
      apply IsOpen.mem_nhds isOpen_Iio
      simp [ENNReal.ofReal_pos.mpr hδ_pos]
    exact h_tendsto h_mem
  exact h_event.exists_forall_of_atTop

end QualitativeDecayC

end Pandrosion
