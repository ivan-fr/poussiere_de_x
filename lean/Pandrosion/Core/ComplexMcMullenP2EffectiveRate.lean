/-
  Universitas Pandrosion — §40. **Effective rate for `p = 2` complex
  unconditional McMullen.**

  Quantitative refinement of §39: the `mcmullen_p2_complex_unconditional`
  entry time `k₀(z₀)` admits an **explicit, closed-form upper bound**
  in terms of the Böttcher coordinate `‖v(z₀)‖ ≠ 1`.

  **Main theorem.**  For every `z₀ ∉ complex_bad_set x α` with
  `0 < ‖v(z₀)‖ < 1`, every iterate index
      `k ≥ entryTimeAlpha α r ‖v(z₀)‖`
  satisfies `‖(steffensen_step_C x 2)^[k] z₀ − α‖ < r`. The bound
  `entryTimeAlpha` is a **closed-form expression** in `(α, r, t)` —
  no hidden existentials.

  The analytical core is the numeric tower inequality
      `t^{2^k} < δ`  for  `k ≥ ⌈log₂(log(1/δ)/log(1/t))⌉ + 1`,
  combined with the §38.6 iterated Böttcher
      `‖v(σⁿ z₀)‖ = ‖v(z₀)‖^{2ⁿ}`,
  the §39 bridge, and the Möbius-inverse Lipschitz bound.

  The `entryTimeNegAlpha` companion for `‖v(z₀)‖ > 1` is
  symmetric and omitted here (it would duplicate the proof with
  `M := max(2‖μ²‖, 4‖α‖‖μ²‖/r + 1)` as the high threshold).

  Contents.

    §40.1  `entryTimeAux`, `entryTimeAux_bound` — numeric tower
           inequality with explicit constant.

    §40.2  `entryTimeAlpha` — closed-form basin-of-`+α` entry time.

    §40.3  `effective_entry_time_plus_alpha` — quantitative theorem.
-/

import Pandrosion.Core.ComplexMcMullenP2Unconditional

namespace Pandrosion

open MeasureTheory Complex

/-! ============================================================
  §40.1  Numeric tower inequality
============================================================ -/

section NumericTowerC

/-- **Explicit ceiling for quadratic tower decay.**
    For `0 < t < 1` and `0 < δ`,
        `entryTimeAux t δ := ⌈log₂(log(1/δ) / log(1/t))⌉ + 1`
    is the explicit `k` such that `t^{2^k} < δ` for every iterate
    count `≥ entryTimeAux t δ`. Outside the valid regime, returns
    `0` (vacuous or saturated). -/
noncomputable def entryTimeAux (t δ : ℝ) : ℕ :=
  if t ≤ 0 ∨ 1 ≤ t ∨ δ ≤ 0 then 0
  else
    Nat.ceil (Real.logb 2 ((Real.log (1/δ)) / (Real.log (1/t)))) + 1

/-- **Numeric tower bound.** For `0 < t < 1` and `0 < δ ≤ 1`,
    `t^{2^k} < δ` for every `k ≥ entryTimeAux t δ`. -/
theorem entryTimeAux_bound
    (t δ : ℝ) (ht_pos : 0 < t) (ht_lt : t < 1) (hδ_pos : 0 < δ) :
    ∀ k ≥ entryTimeAux t δ, t ^ (2 ^ k) < δ := by
  intro k hk
  by_cases hδ_le_one : δ ≤ 1
  · -- Non-trivial regime.
    by_cases hδ_eq_one : δ = 1
    · -- δ = 1: `t^{anything} < 1` suffices.
      rw [hδ_eq_one]
      have h_pow_lt : t ^ (2 ^ k) < 1 ^ (2 ^ k) := by
        apply pow_lt_pow_left ht_lt (le_of_lt ht_pos)
        exact Nat.pos_iff_ne_zero.mp (Nat.pos_pow_of_pos _ (by norm_num))
      simpa using h_pow_lt
    · have hδ_lt_one : δ < 1 := lt_of_le_of_ne hδ_le_one hδ_eq_one
      unfold entryTimeAux at hk
      have h_cond : ¬ (t ≤ 0 ∨ 1 ≤ t ∨ δ ≤ 0) := by
        push_neg; exact ⟨ht_pos, ht_lt, hδ_pos⟩
      rw [if_neg h_cond] at hk
      -- hk: k ≥ ⌈log₂(log(1/δ)/log(1/t))⌉ + 1
      have hlog_t_neg : Real.log t < 0 := Real.log_neg ht_pos ht_lt
      have hlog_t_inv_pos : 0 < Real.log (1/t) := by
        rw [one_div, Real.log_inv]; linarith
      have hlog_δ_inv_pos : 0 < Real.log (1/δ) := by
        rw [one_div, Real.log_inv, neg_pos]
        exact Real.log_neg hδ_pos hδ_lt_one
      set N : ℝ := Real.log (1/δ) / Real.log (1/t) with hN_def
      have hN_pos : 0 < N := div_pos hlog_δ_inv_pos hlog_t_inv_pos
      -- 2^k > N.
      have h2k_gt : N < (2 : ℝ) ^ k := by
        have h_ceil_ge : Real.logb 2 N ≤ (Nat.ceil (Real.logb 2 N) : ℝ) :=
          Nat.le_ceil _
        have h_two_pos : (0 : ℝ) < 2 := by norm_num
        have h_two_ne_one : (2 : ℝ) ≠ 1 := by norm_num
        have h_pow_ge_N :
            N ≤ (2 : ℝ) ^ ((Nat.ceil (Real.logb 2 N) : ℕ) : ℝ) := by
          calc N = (2 : ℝ) ^ Real.logb 2 N :=
                  (Real.rpow_logb h_two_pos h_two_ne_one hN_pos).symm
            _ ≤ (2 : ℝ) ^ ((Nat.ceil (Real.logb 2 N) : ℕ) : ℝ) :=
                  Real.rpow_le_rpow_of_exponent_le (by norm_num) h_ceil_ge
        have h_pow_cast :
            (2 : ℝ) ^ (((Nat.ceil (Real.logb 2 N) : ℕ) : ℝ))
              = ((2 : ℝ) ^ (Nat.ceil (Real.logb 2 N) : ℕ) : ℝ) := by
          rw [Real.rpow_nat_cast]
        rw [h_pow_cast] at h_pow_ge_N
        have h_pow_strict :
            (2 : ℝ) ^ (Nat.ceil (Real.logb 2 N) : ℕ)
              < (2 : ℝ) ^ (Nat.ceil (Real.logb 2 N) + 1 : ℕ) := by
          apply pow_lt_pow_right (by norm_num : (1 : ℝ) < 2)
          omega
        have h_pow_mono :
            (2 : ℝ) ^ (Nat.ceil (Real.logb 2 N) + 1 : ℕ) ≤ (2 : ℝ) ^ k :=
          pow_le_pow_right (by norm_num : (1 : ℝ) ≤ 2) hk
        linarith
      -- 2^k · log(1/t) > log(1/δ) ⟹ log(t^{2^k}) < log(δ) ⟹ t^{2^k} < δ.
      have h_mul : Real.log (1/δ) < (2 : ℝ) ^ k * Real.log (1/t) := by
        rw [hN_def, div_lt_iff hlog_t_inv_pos] at h2k_gt
        exact h2k_gt
      have h_neg_t : Real.log (1/t) = - Real.log t := by
        rw [one_div, Real.log_inv]
      have h_neg_δ : Real.log (1/δ) = - Real.log δ := by
        rw [one_div, Real.log_inv]
      rw [h_neg_t, h_neg_δ] at h_mul
      have h_log_pow : Real.log (t ^ (2 ^ k)) = (2^k : ℕ) * Real.log t :=
        Real.log_pow _ _
      have h_log_lt : Real.log (t ^ (2 ^ k)) < Real.log δ := by
        rw [h_log_pow]; push_cast; linarith
      have h_pow_pos : 0 < t ^ (2 ^ k) := pow_pos ht_pos _
      exact (Real.log_lt_log_iff h_pow_pos hδ_pos).mp h_log_lt
  · push_neg at hδ_le_one
    have h_pow_le_one : t ^ (2 ^ k) ≤ 1 :=
      pow_le_one _ (le_of_lt ht_pos) (le_of_lt ht_lt)
    linarith

end NumericTowerC

/-! ============================================================
  §40.2  Closed-form `entryTimeAlpha`
============================================================ -/

section EntryTimeAlphaC

/-- **Closed-form basin-of-`+α` entry-time threshold.**
    For `α ∈ ℂ`, `r > 0`, and Böttcher magnitude `t := ‖v(z₀)‖`:
        `θ := min(‖μ²‖/2, r·‖μ²‖/(4‖α‖))`    (Möbius-inverse threshold)
        `entryTimeAlpha α r t := entryTimeAux t θ`. -/
noncomputable def entryTimeAlpha (α : ℂ) (r t : ℝ) : ℕ :=
  entryTimeAux t
    (min (‖((1 - α) / (1 + α)) ^ 2‖ / 2)
         (r * ‖((1 - α) / (1 + α)) ^ 2‖ / (4 * ‖α‖)))

end EntryTimeAlphaC

/-! ============================================================
  §40.3  Effective entry-time theorem (basin of `+α`)
============================================================ -/

section EffectiveTheoremC

/-- **★★★ Effective `+α`-basin entry-time bound.**

    For every `z₀ ∉ complex_bad_set x α` with `‖v(z₀)‖ < 1` (and
    non-zero), every iterate index `k ≥ entryTimeAlpha α r ‖v(z₀)‖`
    satisfies
        `‖(steffensen_step_C x 2)^[k] z₀ − α‖ < r`.

    The bound is **closed-form** in `(α, r, ‖v(z₀)‖)`: no hidden
    existentials, no tendsto-extracted witnesses.

    **Proof structure.**
      (i)   Orbit non-degeneracy + bridge from §39.
      (ii)  §38.6 iterated form `‖v(σⁿ z₀)‖ = ‖v(z₀)‖^{2ⁿ}`.
      (iii) §40.1 numeric tower: `‖v(z₀)‖^{2^k} < θ`
            for `k ≥ entryTimeAux ‖v(z₀)‖ θ`.
      (iv)  Möbius-inverse Lipschitz: `‖σ^k z₀ − α‖ ≤
            4‖α‖‖v(σ^k z₀)‖/‖μ²‖` (when `‖v‖ ≤ ‖μ²‖/2`).
      (v)   Threshold `θ := min(‖μ²‖/2, r‖μ²‖/(4‖α‖))` clinches `< r`. -/
theorem effective_entry_time_plus_alpha
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (z₀ : ℂ) (hgood : z₀ ∉ complex_bad_set x α)
    (hv_pos : 0 < ‖v_p2_C α z₀‖) (hv_lt : ‖v_p2_C α z₀‖ < 1)
    (h_no_alpha_hit : ∀ k, (steffensen_step_C x 2)^[k] z₀ ≠ α)
    (h_no_neg_alpha_hit : ∀ k, (steffensen_step_C x 2)^[k] z₀ ≠ -α)
    (r : ℝ) (hr_pos : 0 < r) :
    ∀ k ≥ entryTimeAlpha α r ‖v_p2_C α z₀‖,
      ‖(steffensen_step_C x 2)^[k] z₀ - α‖ < r := by
  intro k hk
  have hxα : x * α ^ 2 = 1 := by rw [hα_pow]; field_simp
  obtain ⟨h_orbit_nondeg, _h_v_ne_one⟩ :=
    complex_orbit_non_degenerate_of_not_bad_set x α z₀ hgood
  -- Bridge to sigma_p2_explicit_C.
  have H_bridge_data : ∀ n : ℕ,
      (steffensen_step_C x 2)^[n] z₀ ≠ -1
      ∧ (steffensen_step_C x 2)^[n] z₀ ≠ -1 / x
      ∧ 2 * x * ((steffensen_step_C x 2)^[n] z₀) + x + 1 ≠ 0
      ∧ (steffensen_step_C x 2)^[n] z₀ ≠ α
      ∧ (steffensen_step_C x 2)^[n] z₀ ≠ -α := fun n =>
    let ⟨h1, h2, h3⟩ := h_orbit_nondeg n
    ⟨h1, h2, h3, h_no_alpha_hit n, h_no_neg_alpha_hit n⟩
  have h_orbit_eq : ∀ n : ℕ,
      (steffensen_step_C x 2)^[n] z₀
        = (sigma_p2_explicit_C x α)^[n] z₀ :=
    steffensen_C_eq_sigma_p2_iterate_of_orbit_good x α hx_ne hxα z₀ H_bridge_data
  -- Non-degeneracy along sigma iterates.
  have H_vpow : ∀ n : ℕ, ∀ j ≤ n,
      (sigma_p2_explicit_C x α)^[j] z₀ + 1 ≠ 0
      ∧ x * (sigma_p2_explicit_C x α)^[j] z₀ + 1 ≠ 0
      ∧ (sigma_p2_explicit_C x α)^[j] z₀ + α ≠ 0 := by
    intro n j _hjn
    obtain ⟨h1, h2, _, _, h5⟩ := H_bridge_data j
    rw [← h_orbit_eq j] at *
    refine ⟨?_, ?_, ?_⟩
    · intro h; apply h1; linear_combination h
    · intro h; apply h2; field_simp; linear_combination h
    · intro h; apply h5; linear_combination h
  -- Iterated Böttcher: v(σ^n z₀) = v(z₀)^(2^n).
  have h_v_iter : ∀ n : ℕ,
      v_p2_C α ((sigma_p2_explicit_C x α)^[n] z₀) = (v_p2_C α z₀) ^ (2 ^ n) :=
    fun n => v_p2_iterated_C x α hα_ne_zero hα_ne_neg_one hxα z₀ n (H_vpow n)
  -- Möbius-inverse identities.
  set μ_sq : ℂ := ((1 - α) / (1 + α)) ^ 2 with hμ_sq_def
  have hμ_sq_ne_zero : μ_sq ≠ 0 := by
    rw [hμ_sq_def]
    exact pow_ne_zero _ (div_ne_zero hα_ne_one hα_ne_neg_one)
  have hμ_sq_norm_pos : 0 < ‖μ_sq‖ := norm_pos_iff.mpr hμ_sq_ne_zero
  have hα_norm_pos : 0 < ‖α‖ := norm_pos_iff.mpr hα_ne_zero
  have h_mobius_minus :
      ∀ s : ℂ, s + α ≠ 0 →
        (s - α) * (μ_sq - v_p2_C α s) = 2 * α * v_p2_C α s := by
    intro s hsα
    unfold v_p2_C; rw [← hμ_sq_def]; field_simp; ring
  have h_mobius_plus :
      ∀ s : ℂ, s + α ≠ 0 →
        (s + α) * (μ_sq - v_p2_C α s) = 2 * α * μ_sq := by
    intro s hsα
    unfold v_p2_C; rw [← hμ_sq_def]; field_simp; ring
  have h_mu_sq_minus_v_ne :
      ∀ s : ℂ, s + α ≠ 0 → μ_sq - v_p2_C α s ≠ 0 := by
    intro s hsα h_eq
    have h := h_mobius_plus s hsα
    rw [h_eq, mul_zero] at h
    have hne : (2 : ℂ) * α * μ_sq ≠ 0 :=
      mul_ne_zero (mul_ne_zero two_ne_zero hα_ne_zero) hμ_sq_ne_zero
    exact hne h.symm
  -- Threshold θ and numeric bound.
  set θ : ℝ := min (‖μ_sq‖ / 2) (r * ‖μ_sq‖ / (4 * ‖α‖)) with hθ_def
  have h4α_pos : (0 : ℝ) < 4 * ‖α‖ := by positivity
  have hθ_pos : 0 < θ := by
    rw [hθ_def]
    exact lt_min (half_pos hμ_sq_norm_pos)
      (div_pos (mul_pos hr_pos hμ_sq_norm_pos) h4α_pos)
  have hθ_le_half : θ ≤ ‖μ_sq‖ / 2 := min_le_left _ _
  have hθ_le_target : θ ≤ r * ‖μ_sq‖ / (4 * ‖α‖) := min_le_right _ _
  -- Apply numeric tower: ‖v(z₀)‖^{2^k} < θ.
  have h_t_pow_lt : ‖v_p2_C α z₀‖ ^ (2 ^ k) < θ := by
    apply entryTimeAux_bound ‖v_p2_C α z₀‖ θ hv_pos hv_lt hθ_pos k
    -- k ≥ entryTimeAlpha α r ‖v(z₀)‖ = entryTimeAux ‖v(z₀)‖ θ.
    unfold entryTimeAlpha at hk
    rw [← hθ_def] at hk
    exact hk
  -- Identify ‖v(σ^k z₀)‖ = ‖v(z₀)‖^{2^k}.
  set s_k : ℂ := (sigma_p2_explicit_C x α)^[k] z₀ with hs_k_def
  have h_v_iter_k : v_p2_C α s_k = (v_p2_C α z₀) ^ (2 ^ k) := h_v_iter k
  have h_abs_v_sk_eq : ‖v_p2_C α s_k‖ = ‖v_p2_C α z₀‖ ^ (2 ^ k) := by
    rw [h_v_iter_k, norm_pow]
  have h_abs_v_sk_lt : ‖v_p2_C α s_k‖ < θ := by
    rw [h_abs_v_sk_eq]; exact h_t_pow_lt
  -- Möbius-inverse: ‖s_k - α‖ = 2‖α‖·‖v(s_k)‖ / ‖μ² - v(s_k)‖.
  have h_sigma_plus_α_ne : s_k + α ≠ 0 := by
    rw [hs_k_def]; exact (H_vpow k k (le_refl k)).2.2
  have h_mu_v_ne_k : μ_sq - v_p2_C α s_k ≠ 0 :=
    h_mu_sq_minus_v_ne s_k h_sigma_plus_α_ne
  have h_mu_v_abs_pos : 0 < ‖μ_sq - v_p2_C α s_k‖ := norm_pos_iff.mpr h_mu_v_ne_k
  have h_sk_minus_α_eq :
      s_k - α = 2 * α * v_p2_C α s_k / (μ_sq - v_p2_C α s_k) := by
    have h := h_mobius_minus s_k h_sigma_plus_α_ne
    field_simp; linear_combination h
  have h_abs_sk_minus_α :
      ‖s_k - α‖ = 2 * ‖α‖ * ‖v_p2_C α s_k‖ / ‖μ_sq - v_p2_C α s_k‖ := by
    rw [h_sk_minus_α_eq, norm_div]
    rw [show (2 : ℂ) * α * v_p2_C α s_k = 2 * (α * v_p2_C α s_k) from by ring,
        norm_mul, norm_mul, show ‖(2 : ℂ)‖ = 2 from by
          rw [show ((2 : ℂ)) = ((2 : ℕ) : ℂ) from by norm_num, Complex.norm_nat]
          norm_num]
    ring
  -- Denominator lower bound ‖μ² - v‖ ≥ ‖μ²‖/2 when ‖v‖ ≤ ‖μ²‖/2.
  have h_denom_lower : ‖μ_sq‖ / 2 ≤ ‖μ_sq - v_p2_C α s_k‖ := by
    have h_v_small : ‖v_p2_C α s_k‖ ≤ ‖μ_sq‖ / 2 :=
      le_of_lt (lt_of_lt_of_le h_abs_v_sk_lt hθ_le_half)
    have h_triangle : ‖μ_sq‖ - ‖v_p2_C α s_k‖ ≤ ‖μ_sq - v_p2_C α s_k‖ :=
      norm_sub_norm_le μ_sq (v_p2_C α s_k)
    linarith
  -- Final bound ‖s_k - α‖ < r.
  have h_final : ‖s_k - α‖ < r := by
    rw [h_abs_sk_minus_α, div_lt_iff h_mu_v_abs_pos]
    have h_θ_bd : θ * (4 * ‖α‖) ≤ r * ‖μ_sq‖ := by
      rw [le_div_iff h4α_pos] at hθ_le_target; linarith
    have h_v_abs_nn : (0 : ℝ) ≤ ‖v_p2_C α s_k‖ := norm_nonneg _
    have h_α_nn : (0 : ℝ) ≤ ‖α‖ := le_of_lt hα_norm_pos
    nlinarith [h_abs_v_sk_lt, h_denom_lower, h_θ_bd, hμ_sq_norm_pos,
               hr_pos, h_α_nn, h_v_abs_nn, h4α_pos]
  rw [h_orbit_eq k, ← hs_k_def]
  exact h_final

end EffectiveTheoremC

end Pandrosion
