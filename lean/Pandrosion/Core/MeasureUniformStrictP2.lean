/-
  Universitas Pandrosion — §46. **Strict measure-uniform
  complexity for `p = 2` on ℂ.**

  Strengthens §45 by replacing the `z₀`-dependent per-point
  bound `K(z₀) := max(entryTimeAlpha α r ‖v(z₀)‖,
  entryTimeNegAlpha α r ‖v(z₀)‖)` with a **uniform**
  `K*(δ, r) : ℕ` independent of `z₀`:

      *For every `δ, r > 0`, there exists `N : ℕ` and
      `K* : ℕ` such that, off a set of Lebesgue measure `< δ`,
      every orbit enters the cyclotomic basin of radius `r`
      within `K*` iterations.*

  The uniform constant
      `K* := max(entryTimeAlpha α r (1 − 1/(N+1)),
                  entryTimeNegAlpha α r (1 + 1/(N+1)))`
  is obtained via the monotonicity of `entryTimeAux` (in `t`,
  for fixed `δ > 0`) and `entryTimeGrowAux` (antitone in `t`,
  for fixed `M > 0`). Outside the §43 slow set, `‖v(z₀)‖`
  stays bounded away from `1` by `ε_N := 1/(N+1)`, letting
  the bounds evaluate at the boundary values.

  Contents.

    §46.1  `entryTimeAux_mono_t` — monotonicity in `t`.

    §46.2  `entryTimeGrowAux_antitone_t` — dual antitone bound.

    §46.3  `entryTimeAlpha_mono_t`, `entryTimeNegAlpha_antitone_t`.

    §46.4  `mcmullen_p2_complex_measure_uniform_strict` — the
           strict-uniform complexity theorem.
-/

import Pandrosion.Core.MeasureUniformComplexityP2

namespace Pandrosion

open MeasureTheory Complex Filter Topology

/-! ============================================================
  §46.1  `entryTimeAux` monotonicity in `t`
============================================================ -/

section EntryTimeAuxMonoC

/-- **`entryTimeAux` is monotone-increasing in `t` (for `0 < t < 1`,
    `0 < δ`).** Larger `t` ⟹ slower decay ⟹ larger entry time. -/
theorem entryTimeAux_mono_t
    (δ : ℝ) (hδ_pos : 0 < δ) (hδ_lt_one : δ < 1)
    (t₁ t₂ : ℝ) (ht1_pos : 0 < t₁) (ht2_lt : t₂ < 1) (ht_le : t₁ ≤ t₂) :
    entryTimeAux t₁ δ ≤ entryTimeAux t₂ δ := by
  have ht2_pos : 0 < t₂ := lt_of_lt_of_le ht1_pos ht_le
  have ht1_lt : t₁ < 1 := lt_of_le_of_lt ht_le ht2_lt
  unfold entryTimeAux
  have h_cond1 : ¬ (t₁ ≤ 0 ∨ 1 ≤ t₁ ∨ δ ≤ 0) := by
    push_neg; exact ⟨ht1_pos, ht1_lt, hδ_pos⟩
  have h_cond2 : ¬ (t₂ ≤ 0 ∨ 1 ≤ t₂ ∨ δ ≤ 0) := by
    push_neg; exact ⟨ht2_pos, ht2_lt, hδ_pos⟩
  rw [if_neg h_cond1, if_neg h_cond2]
  have hlog_t1_inv_pos : 0 < Real.log (1/t₁) := by
    rw [one_div, Real.log_inv, neg_pos]; exact Real.log_neg ht1_pos ht1_lt
  have hlog_t2_inv_pos : 0 < Real.log (1/t₂) := by
    rw [one_div, Real.log_inv, neg_pos]; exact Real.log_neg ht2_pos ht2_lt
  have h_t2_inv_pos : 0 < 1 / t₂ := one_div_pos.mpr ht2_pos
  have h_t1_inv_pos : 0 < 1 / t₁ := one_div_pos.mpr ht1_pos
  have h_inv_le : (1 : ℝ) / t₂ ≤ 1 / t₁ :=
    one_div_le_one_div_of_le ht1_pos ht_le
  have h_log_inv_le : Real.log (1/t₂) ≤ Real.log (1/t₁) :=
    (Real.log_le_log_iff h_t2_inv_pos h_t1_inv_pos).mpr h_inv_le
  have hlog_δ_inv_pos : 0 < Real.log (1/δ) := by
    rw [one_div, Real.log_inv, neg_pos]; exact Real.log_neg hδ_pos hδ_lt_one
  have h_div_le : Real.log (1/δ) / Real.log (1/t₁)
                ≤ Real.log (1/δ) / Real.log (1/t₂) :=
    div_le_div_of_nonneg_left (le_of_lt hlog_δ_inv_pos) hlog_t2_inv_pos h_log_inv_le
  have h_pos1 : 0 < Real.log (1/δ) / Real.log (1/t₁) :=
    div_pos hlog_δ_inv_pos hlog_t1_inv_pos
  have h_pos2 : 0 < Real.log (1/δ) / Real.log (1/t₂) :=
    div_pos hlog_δ_inv_pos hlog_t2_inv_pos
  have h_logb_le :
      Real.logb 2 (Real.log (1/δ) / Real.log (1/t₁))
        ≤ Real.logb 2 (Real.log (1/δ) / Real.log (1/t₂)) :=
    (Real.logb_le_logb (by norm_num : (1:ℝ) < 2) h_pos1 h_pos2).mpr h_div_le
  have h_ceil_le :
      Nat.ceil (Real.logb 2 (Real.log (1/δ) / Real.log (1/t₁)))
        ≤ Nat.ceil (Real.logb 2 (Real.log (1/δ) / Real.log (1/t₂))) :=
    Nat.ceil_le_ceil _ _ h_logb_le
  exact Nat.add_le_add_right h_ceil_le 1

end EntryTimeAuxMonoC

/-! ============================================================
  §46.2  `entryTimeGrowAux` is antitone in `t`
============================================================ -/

section EntryTimeGrowAntitoneC

/-- **`entryTimeGrowAux` is antitone in `t` (for `1 < t`, `0 < M`).**
    Larger `t` ⟹ faster growth ⟹ smaller entry time. -/
theorem entryTimeGrowAux_antitone_t
    (M : ℝ)
    (t₁ t₂ : ℝ) (ht1_gt : 1 < t₁) (ht_le : t₁ ≤ t₂) :
    entryTimeGrowAux t₂ M ≤ entryTimeGrowAux t₁ M := by
  have ht2_gt : 1 < t₂ := lt_of_lt_of_le ht1_gt ht_le
  unfold entryTimeGrowAux
  by_cases hM_le_one : M ≤ 1
  · have h_cond1 : t₁ ≤ 1 ∨ M ≤ 1 := Or.inr hM_le_one
    have h_cond2 : t₂ ≤ 1 ∨ M ≤ 1 := Or.inr hM_le_one
    rw [if_pos h_cond1, if_pos h_cond2]
  · push_neg at hM_le_one
    have h_cond1 : ¬ (t₁ ≤ 1 ∨ M ≤ 1) := by push_neg; exact ⟨ht1_gt, hM_le_one⟩
    have h_cond2 : ¬ (t₂ ≤ 1 ∨ M ≤ 1) := by push_neg; exact ⟨ht2_gt, hM_le_one⟩
    rw [if_neg h_cond1, if_neg h_cond2]
    have hlog_t1_pos : 0 < Real.log t₁ := Real.log_pos ht1_gt
    have hlog_t2_pos : 0 < Real.log t₂ := Real.log_pos ht2_gt
    have hlog_M_pos : 0 < Real.log M := Real.log_pos hM_le_one
    have h_log_le : Real.log t₁ ≤ Real.log t₂ :=
      (Real.log_le_log_iff (by linarith) (by linarith)).mpr ht_le
    -- log M / log t₁ ≥ log M / log t₂.
    have h_div_le : Real.log M / Real.log t₂ ≤ Real.log M / Real.log t₁ :=
      div_le_div_of_nonneg_left (le_of_lt hlog_M_pos) hlog_t1_pos h_log_le
    have h_pos1 : 0 < Real.log M / Real.log t₁ :=
      div_pos hlog_M_pos hlog_t1_pos
    have h_pos2 : 0 < Real.log M / Real.log t₂ :=
      div_pos hlog_M_pos hlog_t2_pos
    have h_logb_le :
        Real.logb 2 (Real.log M / Real.log t₂)
          ≤ Real.logb 2 (Real.log M / Real.log t₁) :=
      (Real.logb_le_logb (by norm_num : (1:ℝ) < 2) h_pos2 h_pos1).mpr h_div_le
    exact Nat.add_le_add_right (Nat.ceil_le_ceil _ _ h_logb_le) 1

end EntryTimeGrowAntitoneC

/-! ============================================================
  §46.3  Monotonicity of `entryTimeAlpha`/`entryTimeNegAlpha`
============================================================ -/

section EntryTimeAlphaMonoC

/-- **`entryTimeAlpha α r ·` is monotone-increasing in `t`** on
    `(0, 1)`. -/
theorem entryTimeAlpha_mono_t
    (α : ℂ) (r : ℝ) (hr_pos : 0 < r)
    (hα_ne_zero : α ≠ 0) (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hθ_lt_one : min (‖((1 - α) / (1 + α)) ^ 2‖ / 2)
                     (r * ‖((1 - α) / (1 + α)) ^ 2‖ / (4 * ‖α‖)) < 1)
    (t₁ t₂ : ℝ) (ht1_pos : 0 < t₁) (ht2_lt : t₂ < 1) (ht_le : t₁ ≤ t₂) :
    entryTimeAlpha α r t₁ ≤ entryTimeAlpha α r t₂ := by
  unfold entryTimeAlpha
  set θ : ℝ := min (‖((1 - α) / (1 + α)) ^ 2‖ / 2)
                   (r * ‖((1 - α) / (1 + α)) ^ 2‖ / (4 * ‖α‖)) with hθ_def
  have hμ_ne : ((1 - α) / (1 + α)) ^ 2 ≠ 0 :=
    pow_ne_zero _ (div_ne_zero hα_ne_one hα_ne_neg_one)
  have hμ_norm_pos : 0 < ‖((1 - α) / (1 + α)) ^ 2‖ := norm_pos_iff.mpr hμ_ne
  have hα_norm_pos : 0 < ‖α‖ := norm_pos_iff.mpr hα_ne_zero
  have hθ_pos : 0 < θ := by
    rw [hθ_def]
    exact lt_min (half_pos hμ_norm_pos)
      (div_pos (mul_pos hr_pos hμ_norm_pos) (by positivity))
  exact entryTimeAux_mono_t θ hθ_pos hθ_lt_one t₁ t₂ ht1_pos ht2_lt ht_le

/-- **`entryTimeNegAlpha α r ·` is antitone in `t`** on `(1, ∞)`. -/
theorem entryTimeNegAlpha_antitone_t
    (α : ℂ) (r : ℝ) (hr_pos : 0 < r)
    (hα_ne_zero : α ≠ 0) (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (t₁ t₂ : ℝ) (ht1_gt : 1 < t₁) (ht_le : t₁ ≤ t₂) :
    entryTimeNegAlpha α r t₂ ≤ entryTimeNegAlpha α r t₁ := by
  unfold entryTimeNegAlpha
  set M : ℝ := max (2 * ‖((1 - α) / (1 + α)) ^ 2‖)
                   (4 * ‖α‖ * ‖((1 - α) / (1 + α)) ^ 2‖ / r + 1) with hM_def
  have hμ_ne : ((1 - α) / (1 + α)) ^ 2 ≠ 0 :=
    pow_ne_zero _ (div_ne_zero hα_ne_one hα_ne_neg_one)
  have hμ_norm_pos : 0 < ‖((1 - α) / (1 + α)) ^ 2‖ := norm_pos_iff.mpr hμ_ne
  have hα_norm_pos : 0 < ‖α‖ := norm_pos_iff.mpr hα_ne_zero
  have hM_pos : 0 < M := by
    rw [hM_def]
    apply lt_of_lt_of_le _ (le_max_right _ _)
    have : (0 : ℝ) ≤ 4 * ‖α‖ * ‖((1 - α) / (1 + α)) ^ 2‖ / r := by positivity
    linarith
  let _ := hM_pos  -- referenced for documentation; trivially discharged below.
  exact entryTimeGrowAux_antitone_t M t₁ t₂ ht1_gt ht_le

end EntryTimeAlphaMonoC

/-! ============================================================
  §46.4  Strict measure-uniform complexity theorem
============================================================ -/

section StrictMeasureUniformC

/-- **★★★★★★ Strict measure-uniform complexity for `p = 2` on ℂ.**

    For every `x ∈ ℂ \ {0, 1}` (with `α ≠ 0, ±1` and `α² = 1/x`),
    every `R, δ, r > 0`, there exist `N, K* : ℕ` such that:

      (i)  The bad set
           `slowSet α (1/(N+1)) ∩ B(0,R) ∪ alphaHitSet ∪ complex_bad_set`
           has volume `< δ`.

      (ii) For every `z₀` outside the bad set, **some** cyclotomic
           anchor `s ∈ Fin 2` captures `σ^k z₀` within distance
           `r` for all `k ≥ K*` — **`K*` independent of `z₀`**.

    The uniform constant:
        `K* := max(entryTimeAlpha α r (1 − 1/(N+1)),
                     entryTimeNegAlpha α r (1 + 1/(N+1)))`
    via §46.3's monotonicity: for `z₀ ∉ slowSet α (1/(N+1))` with
    `‖v(z₀)‖ < 1`, `‖v(z₀)‖ ≤ 1 − 1/(N+1)`, so
    `entryTimeAlpha α r ‖v(z₀)‖ ≤ entryTimeAlpha α r (1 − 1/(N+1))`;
    symmetrically for `‖v(z₀)‖ > 1`. -/
theorem mcmullen_p2_complex_measure_uniform_strict
    (x : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (R : ℝ) (δ : ℝ) (hδ_pos : 0 < δ) (r : ℝ) (hr_pos : 0 < r)
    (hθ_lt_one : min (‖((1 - α) / (1 + α)) ^ 2‖ / 2)
                     (r * ‖((1 - α) / (1 + α)) ^ 2‖ / (4 * ‖α‖)) < 1) :
    ∃ N : ℕ, ∃ K_star : ℕ,
      volume ((slowSet α (1 / ((N : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R)
              ∪ alphaHitSet x α
              ∪ complex_bad_set x α) < ENNReal.ofReal δ ∧
      ∀ z₀ : ℂ,
        z₀ ∉ slowSet α (1 / ((N : ℝ) + 1)) →
        z₀ ∉ alphaHitSet x α →
        z₀ ∉ complex_bad_set x α →
        0 < ‖v_p2_C α z₀‖ →
        ∃ s : Fin 2, ∀ k ≥ K_star,
          ‖(steffensen_step_C x 2)^[k] z₀ - cycAnchor α 2 s‖ < r := by
  -- Pull the §45 result.
  obtain ⟨N, h_meas, h_entry⟩ :=
    mcmullen_p2_complex_measure_uniform x hx_ne hx_ne_one α hα_ne_zero
      hα_ne_neg_one hα_ne_one hα_pow R δ hδ_pos r hr_pos
  set ε_N : ℝ := 1 / ((N : ℝ) + 1) with hε_N_def
  have hε_N_pos : 0 < ε_N := by rw [hε_N_def]; positivity
  have hε_N_lt_one : ε_N ≤ 1 := by
    rw [hε_N_def]
    apply div_le_one_of_le
    · linarith [show (0 : ℝ) ≤ (N : ℝ) by positivity]
    · positivity
  set K_star : ℕ :=
    max (entryTimeAlpha α r (1 - ε_N)) (entryTimeNegAlpha α r (1 + ε_N))
  refine ⟨N, K_star, h_meas, ?_⟩
  intro z₀ h_not_slow h_not_α h_not_bad hv_pos
  -- z₀ ∉ alphaHitSet ⟹ z + α ≠ 0 (otherwise σ^0 z₀ = z₀ = -α).
  have h_z_α_ne : z₀ + α ≠ 0 := by
    intro h_eq
    apply h_not_α
    right
    refine ⟨0, ?_⟩
    simp; linear_combination h_eq
  -- Equivalence slowSet ↔ |‖v‖ - 1| ≤ ε_N.
  have h_v_far : ε_N < |‖v_p2_C α z₀‖ - 1| := by
    have h_iff := slowSet_iff_norm_v α ε_N hε_N_lt_one z₀ h_z_α_ne
    by_contra h_le
    push_neg at h_le
    exact h_not_slow (h_iff.mpr h_le)
  -- Apply §45 per-z₀ theorem to get s and per-z₀ bound.
  obtain ⟨s, h_s_bound⟩ := h_entry z₀ h_not_slow h_not_α h_not_bad hv_pos
  refine ⟨s, fun k hk => ?_⟩
  -- Show the §46 uniform K_star dominates the per-z₀ max bound.
  have h_max_le :
      max (entryTimeAlpha α r ‖v_p2_C α z₀‖) (entryTimeNegAlpha α r ‖v_p2_C α z₀‖)
        ≤ K_star := by
    -- Dichotomy on ‖v(z₀)‖ relative to 1.
    rcases lt_or_gt_of_ne (abs_pos.mp (lt_trans hε_N_pos h_v_far)) with h_v_neg | h_v_pos'
    · -- ‖v(z₀)‖ < 1: use entryTimeAlpha monotonicity.
      have h_v_lt : ‖v_p2_C α z₀‖ < 1 := by linarith
      have h_v_le_bound : ‖v_p2_C α z₀‖ ≤ 1 - ε_N := by
        have h_eq : |‖v_p2_C α z₀‖ - 1| = 1 - ‖v_p2_C α z₀‖ := by
          rw [abs_sub_comm]
          exact _root_.abs_of_nonneg (by linarith)
        rw [h_eq] at h_v_far; linarith
      -- entryTimeNegAlpha ‖v‖ = 0 (since ‖v‖ ≤ 1, trivial regime).
      have h_neg_zero :
          entryTimeNegAlpha α r ‖v_p2_C α z₀‖ ≤ K_star := by
        unfold entryTimeNegAlpha entryTimeGrowAux
        have h_cond : ‖v_p2_C α z₀‖ ≤ 1 ∨
                      max (2 * ‖((1 - α) / (1 + α)) ^ 2‖)
                        (4 * ‖α‖ * ‖((1 - α) / (1 + α)) ^ 2‖ / r + 1) ≤ 1 :=
          Or.inl (le_of_lt h_v_lt)
        rw [if_pos h_cond]
        exact Nat.zero_le _
      -- entryTimeAlpha ‖v‖ ≤ entryTimeAlpha (1 - ε_N).
      have h_one_sub_lt : 1 - ε_N < 1 := by linarith
      have h_alpha_mono := entryTimeAlpha_mono_t α r hr_pos hα_ne_zero
        hα_ne_neg_one hα_ne_one hθ_lt_one ‖v_p2_C α z₀‖ (1 - ε_N) hv_pos
        h_one_sub_lt h_v_le_bound
      have h_alpha_le_K : entryTimeAlpha α r ‖v_p2_C α z₀‖ ≤ K_star :=
        h_alpha_mono.trans (le_max_left _ _)
      exact max_le h_alpha_le_K h_neg_zero
    · -- ‖v(z₀)‖ > 1: use entryTimeNegAlpha antitone.
      have h_v_gt : 1 < ‖v_p2_C α z₀‖ := by linarith
      have h_v_ge_bound : 1 + ε_N ≤ ‖v_p2_C α z₀‖ := by
        have h_eq : |‖v_p2_C α z₀‖ - 1| = ‖v_p2_C α z₀‖ - 1 :=
          _root_.abs_of_pos (by linarith)
        rw [h_eq] at h_v_far; linarith
      have hε_plus_gt : 1 < 1 + ε_N := by linarith
      have h_neg_antitone := entryTimeNegAlpha_antitone_t α r hr_pos hα_ne_zero
        hα_ne_neg_one hα_ne_one (1 + ε_N) ‖v_p2_C α z₀‖ hε_plus_gt h_v_ge_bound
      have h_neg_le_K : entryTimeNegAlpha α r ‖v_p2_C α z₀‖ ≤ K_star :=
        h_neg_antitone.trans (le_max_right _ _)
      -- entryTimeAlpha ‖v‖ = 0 (since ‖v‖ ≥ 1, trivial regime).
      have h_alpha_zero :
          entryTimeAlpha α r ‖v_p2_C α z₀‖ ≤ K_star := by
        unfold entryTimeAlpha entryTimeAux
        have h_cond : ‖v_p2_C α z₀‖ ≤ 0 ∨ 1 ≤ ‖v_p2_C α z₀‖ ∨
                      min (‖((1 - α) / (1 + α)) ^ 2‖ / 2)
                        (r * ‖((1 - α) / (1 + α)) ^ 2‖ / (4 * ‖α‖)) ≤ 0 :=
          Or.inr (Or.inl (le_of_lt h_v_gt))
        rw [if_pos h_cond]
        exact Nat.zero_le _
      exact max_le h_alpha_zero h_neg_le_K
  exact h_s_bound k (le_trans h_max_le hk)

end StrictMeasureUniformC

end Pandrosion
