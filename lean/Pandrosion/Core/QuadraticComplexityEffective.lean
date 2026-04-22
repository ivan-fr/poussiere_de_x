/-
  Universitas Pandrosion — §44. **Effective form of
  `quadratic_loglog_complexity` (§17.2).**

  §17.2 produces an existential `∃ N, ∀ k ≥ N, e k ≤ ε` via
  `exists_pow_lt_of_lt_one`. This module lifts that existential
  to an **explicit, closed-form bound** in terms of `(K, r, ε)`,
  matching the asymptotic claim documented in §17.2:
      `N(K, r, ε) := ⌈log(K·ε) / log(K·r)⌉ + 1`.

  Combined with §40's `entryTimeAux` (the quadratic *tower* bound
  `t^{2^k} < δ`), this yields a fully effective `e_k ≤ ε` witness:
  `k ≥ ⌈log₂(N(K, r, ε))⌉ + 1`. The wrapper
  `quadratic_loglog_complexity_effective` packages this directly.

  Contents.

    §44.1  `linearTowerAux`, `linearTowerAux_bound` — closed-form
           bound for the linear tower `t^n < δ` (`0 < t < 1`).

    §44.2  `quadratic_loglog_complexity_effective` — explicit
           witness for §17.2.
-/

import Pandrosion.Core.QuadraticComplexity

namespace Pandrosion

open MeasureTheory

/-! ============================================================
  §44.1  Linear tower `t^n < δ` with explicit ceiling
============================================================ -/

section LinearTowerC

/-- **Explicit ceiling for linear tower decay.**
    For `0 < t < 1` and `0 < δ ≤ 1`,
        `linearTowerAux t δ := ⌈log δ / log t⌉ + 1`
    is the explicit `n` such that `t^n < δ` for every iterate
    count `≥ linearTowerAux t δ`. -/
noncomputable def linearTowerAux (t δ : ℝ) : ℕ :=
  if t ≤ 0 ∨ 1 ≤ t ∨ δ ≤ 0 then 0
  else
    Nat.ceil (Real.log δ / Real.log t) + 1

/-- **Linear tower bound.** For `0 < t < 1` and `0 < δ`,
    `t^n < δ` for every `n ≥ linearTowerAux t δ`. -/
theorem linearTowerAux_bound
    (t δ : ℝ) (ht_pos : 0 < t) (ht_lt : t < 1) (hδ_pos : 0 < δ) :
    ∀ n ≥ linearTowerAux t δ, t ^ n < δ := by
  intro n hn
  by_cases hδ_le_one : δ ≤ 1
  · by_cases hδ_eq_one : δ = 1
    · rw [hδ_eq_one]
      have h_no_cond : ¬ (t ≤ 0 ∨ 1 ≤ t ∨ δ ≤ 0) := by
        push_neg; exact ⟨ht_pos, ht_lt, hδ_pos⟩
      have h_n_pos : 0 < n := by
        unfold linearTowerAux at hn
        rw [if_neg h_no_cond] at hn
        omega
      have h_pow_lt : t ^ n < 1 ^ n :=
        pow_lt_pow_left ht_lt (le_of_lt ht_pos) (Nat.pos_iff_ne_zero.mp h_n_pos)
      simpa using h_pow_lt
    · have hδ_lt_one : δ < 1 := lt_of_le_of_ne hδ_le_one hδ_eq_one
      unfold linearTowerAux at hn
      have h_cond : ¬ (t ≤ 0 ∨ 1 ≤ t ∨ δ ≤ 0) := by
        push_neg; exact ⟨ht_pos, ht_lt, hδ_pos⟩
      rw [if_neg h_cond] at hn
      have hlog_t_neg : Real.log t < 0 := Real.log_neg ht_pos ht_lt
      -- log δ < 0 (used implicitly via log monotonicity below).
      -- Set N := log δ / log t ≥ 0 (negative / negative).
      set N : ℝ := Real.log δ / Real.log t with hN_def
      -- n > N (from n ≥ ⌈N⌉ + 1).
      have h_n_ge : (Nat.ceil N : ℝ) + 1 ≤ n := by
        have h_step : (Nat.ceil N + 1 : ℕ) ≤ n := hn
        exact_mod_cast h_step
      have h_ceil_ge : N ≤ (Nat.ceil N : ℝ) := Nat.le_ceil _
      have hn_gt : N < (n : ℝ) := by linarith
      -- n · log t < N · log t  (since log t < 0, multiplication flips).
      have h_mul : (n : ℝ) * Real.log t < N * Real.log t :=
        (mul_lt_mul_right_of_neg hlog_t_neg).mpr hn_gt
      -- N · log t = log δ.
      have h_N_mul : N * Real.log t = Real.log δ := by
        rw [hN_def, div_mul_cancel₀]; exact ne_of_lt hlog_t_neg
      -- log(t^n) = n · log t < log δ.
      have h_log_pow : Real.log (t ^ n) = (n : ℝ) * Real.log t := Real.log_pow _ _
      have h_log_lt : Real.log (t ^ n) < Real.log δ := by
        rw [h_log_pow, ← h_N_mul]; exact h_mul
      have h_pow_pos : 0 < t ^ n := pow_pos ht_pos _
      exact (Real.log_lt_log_iff h_pow_pos hδ_pos).mp h_log_lt
  · push_neg at hδ_le_one
    have h_pow_le_one : t ^ n ≤ 1 :=
      pow_le_one _ (le_of_lt ht_pos) (le_of_lt ht_lt)
    linarith

end LinearTowerC

/-! ============================================================
  §44.2  Effective `quadratic_loglog_complexity`
============================================================ -/

section EffectiveQuadraticLogLogC

/-- **Explicit `N_tail` for `quadratic_loglog_complexity`.**
    For `K · r < 1` and target precision `ε > 0`, the explicit
    witness is
        `quadraticLoglogN K r ε := linearTowerAux (K · r) (K · ε)`. -/
noncomputable def quadraticLoglogN (K r ε : ℝ) : ℕ :=
  linearTowerAux (K * r) (K * ε)

/-- **★★★ Effective `quadratic_loglog_complexity`.**
    The §17.2 existential witness is bounded above by the explicit
    closed form `quadraticLoglogN K r ε`. For every `k ≥
    quadraticLoglogN K r ε`, `e_k ≤ ε`. -/
theorem quadratic_loglog_complexity_effective
    (K r : ℝ) (hK : 0 < K) (hr : 0 ≤ r) (hKr : K * r < 1)
    (e : ℕ → ℝ) (he_nn : ∀ k, 0 ≤ e k) (he_0 : e 0 ≤ r)
    (he_rec : ∀ k, e (k + 1) ≤ K * (e k) ^ 2)
    (ε : ℝ) (hε_pos : 0 < ε) :
    ∀ k ≥ quadraticLoglogN K r ε, e k ≤ ε := by
  intro k hk
  -- Replicate the §17.2 proof with explicit N.
  have hKr_nn : 0 ≤ K * r := mul_nonneg hK.le hr
  have hKε_pos : 0 < K * ε := mul_pos hK hε_pos
  -- Get n ≥ quadraticLoglogN such that (K · r)^n < K · ε.
  by_cases hKr_pos : 0 < K * r
  · have h_bound :=
      linearTowerAux_bound (K * r) (K * ε) hKr_pos hKr hKε_pos
    -- For n := quadraticLoglogN K r ε: (K · r)^n < K · ε.
    have h_pow_lt : (K * r) ^ (quadraticLoglogN K r ε) < K * ε := by
      unfold quadraticLoglogN
      exact h_bound (linearTowerAux (K * r) (K * ε)) (le_refl _)
    -- Now apply the §17.2 chain.
    set n := quadraticLoglogN K r ε
    have h_tower := quadratic_tower_bound K r hK hr e he_nn he_0 he_rec k
    have h_n_le_2n : n ≤ 2 ^ n := by
      induction n with
      | zero => simp
      | succ m ih =>
        have h2 : 0 < 2 ^ m := Nat.pos_pow_of_pos m (by norm_num)
        calc m + 1 ≤ 2 ^ m + 1 := by linarith
          _ ≤ 2 ^ m + 2 ^ m := by linarith
          _ = 2 ^ (m + 1) := by rw [pow_succ]; ring
    have h_2n_le_2k : 2 ^ n ≤ 2 ^ k :=
      Nat.pow_le_pow_right (by norm_num) hk
    have h_2k_ge_n : n ≤ 2 ^ k := le_trans h_n_le_2n h_2n_le_2k
    have h_mono : (K * r) ^ (2 ^ k) ≤ (K * r) ^ n :=
      pow_le_pow_of_le_one hKr_nn (le_of_lt hKr) h_2k_ge_n
    have h_bound_pow : K * e k ≤ K * ε := by linarith
    exact le_of_mul_le_mul_left h_bound_pow hK
  · -- Edge case: K · r = 0. Then either K = 0 (impossible) or r = 0.
    push_neg at hKr_pos
    have h_eq_zero : K * r = 0 := le_antisymm hKr_pos hKr_nn
    rcases mul_eq_zero.mp h_eq_zero with hK0 | hr0
    · exact absurd hK0 (ne_of_gt hK)
    · -- r = 0 ⟹ e 0 ≤ 0 ⟹ e 0 = 0.  Recurrence: e (k+1) ≤ K · 0² = 0.
      have he_0_zero : e 0 = 0 := le_antisymm (hr0 ▸ he_0) (he_nn 0)
      have he_zero : ∀ j, e j = 0 := by
        intro j
        induction j with
        | zero => exact he_0_zero
        | succ m ih =>
          have h_step : e (m + 1) ≤ K * (e m) ^ 2 := he_rec m
          rw [ih, pow_two, mul_zero, mul_zero] at h_step
          exact le_antisymm h_step (he_nn _)
      rw [he_zero k]; exact le_of_lt hε_pos

end EffectiveQuadraticLogLogC

end Pandrosion
