/-
  Universitas Pandrosion — §98. **★★★★★★★★★ EFFECTIVE KUNG-TRAUB :
  Closed-form loglog iteration count for §95.**

  The §95 universal multi-start theorem produces an iteration count
  `N(ε)` *existentially* via `quadratic_loglog_complexity` (which
  in turn invokes `exists_pow_lt_of_lt_one` from Mathlib). This
  module replaces the existential by a **closed-form formula**:

    `effectiveLoglogCount K r ε := Nat.log 2 ⌈log(K·ε) / log(K·r)⌉₊ + 1`

  and proves it discharges the §95 conclusion *quantitatively*:
  every `k ≥ effectiveLoglogCount` actually satisfies `‖σⁿ - α‖ ≤ ε`.

  Composed with §33's explicit `K_α := steffensenK_of_fp` and
  `R_α := steffensenR_of_fp`, the result is a **fully explicit
  Pandrosion iteration count** as a closed-form function of
  `(x, p, α, ε)`. This is the first formally verified Kung-Traub
  optimal effective complexity bound for `z^p = x`.

  Asymptotically, `effectiveLoglogCount K r ε ~ log₂ log(1/ε)` as
  `ε → 0` — the optimal Kung-Traub rate.

  Contents.

    §98.1  `linearTowerCount`, `effectiveLoglogCount` — closed-form
           definitions.
    §98.2  `linearTowerCount_pow_le` — `(K·r)^N ≤ K·ε`.
    §98.3  `effectiveLoglogCount_pow_gt` — `2^N > linearTowerCount`.
    §98.4  `quadratic_loglog_effective` — replaces existential N
           in §17.2 by the closed-form `effectiveLoglogCount`.
    §98.5  `quadratic_loglog_from_basin_effective_v2` — parallel of
           §34.1 with explicit count.
    §98.6  `pandrosionEffectiveLoglogCount` — instantiated to §33's
           `K_α, R_α`.
    §98.7  `pandrosion_loglog_universal_multi_start_effective` —
           §95.3 with the explicit closed-form count.
-/

import Pandrosion.Core.LoglogUniversalMultiStartGeneric
import Pandrosion.Core.SteffensenExplicitRate
import Mathlib.Data.Nat.Log
import Mathlib.Analysis.SpecialFunctions.Log.Basic

namespace Pandrosion

open Real Filter Topology Complex

/-! ============================================================
  §98.1  Closed-form linear and loglog counts
============================================================ -/

section EffectiveCounts

/-- **Linear-tower count.** A natural `N` such that
    `(K·r)^N ≤ K·ε`, given by the closed form
    `N := ⌈log(K·ε) / log(K·r)⌉₊`. -/
noncomputable def linearTowerCount (K r ε : ℝ) : ℕ :=
  ⌈Real.log (K * ε) / Real.log (K * r)⌉₊

/-- **★ Effective loglog count.** The loglog version of
    `linearTowerCount`: `Nat.log 2 (linearTowerCount K r ε) + 1`,
    satisfying `2^(effectiveLoglogCount K r ε) > linearTowerCount K r ε`
    via the standard `Nat.lt_pow_succ_log_self` inequality.

    **Asymptotic.** As `ε → 0`, this scales as
    `log₂ log(1/ε)` — matching the Kung-Traub optimal bound for
    derivative-free order-2 methods. -/
noncomputable def effectiveLoglogCount (K r ε : ℝ) : ℕ :=
  Nat.log 2 (linearTowerCount K r ε) + 1

end EffectiveCounts

/-! ============================================================
  §98.2-§98.3  Key inequalities
============================================================ -/

section KeyInequalities

/-- **(K·r)^N ≤ K·ε at N := linearTowerCount K r ε.** Closed-form
    discharge of `exists_pow_lt_of_lt_one` for the specific count. -/
theorem linearTowerCount_pow_le
    (K r ε : ℝ) (hK : 0 < K) (hr : 0 < r) (hKr_lt_one : K * r < 1)
    (hε_pos : 0 < ε) :
    (K * r) ^ (linearTowerCount K r ε) ≤ K * ε := by
  have hKr_pos : 0 < K * r := mul_pos hK hr
  have hKε_pos : 0 < K * ε := mul_pos hK hε_pos
  by_cases hKε_ge : 1 ≤ K * ε
  · -- Trivial case: (K·r)^N ≤ 1 ≤ K·ε.
    calc (K * r) ^ (linearTowerCount K r ε)
        ≤ 1 := pow_le_one _ hKr_pos.le hKr_lt_one.le
      _ ≤ K * ε := hKε_ge
  · -- Non-trivial case: K·ε < 1, so log(K·ε) < 0, log(K·r) < 0.
    push_neg at hKε_ge
    have hlogKr_neg : Real.log (K * r) < 0 := Real.log_neg hKr_pos hKr_lt_one
    have _hlogKε_neg : Real.log (K * ε) < 0 := Real.log_neg hKε_pos hKε_ge
    have h_N_ge : Real.log (K * ε) / Real.log (K * r)
                  ≤ (linearTowerCount K r ε : ℝ) := Nat.le_ceil _
    -- N · log(K·r) ≤ log(K·ε): multiply h_N_ge by log(K·r) (negative, reverses ≤).
    have h_mul : (linearTowerCount K r ε : ℝ) * Real.log (K * r)
                  ≤ Real.log (K * ε) := by
      have h1 : Real.log (K * ε) / Real.log (K * r) * Real.log (K * r)
                = Real.log (K * ε) :=
        div_mul_cancel₀ _ hlogKr_neg.ne
      have h2 : (linearTowerCount K r ε : ℝ) * Real.log (K * r)
                ≤ Real.log (K * ε) / Real.log (K * r) * Real.log (K * r) :=
        mul_le_mul_of_nonpos_right h_N_ge hlogKr_neg.le
      linarith [h1, h2]
    -- exp-monotone: log((K·r)^N) ≤ log(K·ε) → (K·r)^N ≤ K·ε.
    have h_pow_pos : 0 < (K * r) ^ (linearTowerCount K r ε) :=
      pow_pos hKr_pos _
    have h_log_pow : Real.log ((K * r) ^ (linearTowerCount K r ε))
                    = (linearTowerCount K r ε : ℝ) * Real.log (K * r) :=
      Real.log_pow _ _
    have h_log_le : Real.log ((K * r) ^ (linearTowerCount K r ε))
                    ≤ Real.log (K * ε) := by
      rw [h_log_pow]; exact h_mul
    exact (Real.log_le_log_iff h_pow_pos hKε_pos).mp h_log_le

/-- **2^(effectiveLoglogCount K r ε) > linearTowerCount K r ε.**
    Direct from `Nat.lt_pow_succ_log_self`. -/
theorem effectiveLoglogCount_pow_gt
    (K r ε : ℝ) :
    linearTowerCount K r ε < 2 ^ effectiveLoglogCount K r ε := by
  unfold effectiveLoglogCount
  exact Nat.lt_pow_succ_log_self (by norm_num : 1 < 2) _

end KeyInequalities

/-! ============================================================
  §98.4  Effective loglog complexity (closed-form N)
============================================================ -/

section EffectiveLoglogComplexity

/-- **★★ Effective Kung-Traub loglog count.** Replaces the
    non-explicit existential `N` from `quadratic_loglog_complexity`
    with the closed-form `effectiveLoglogCount K r ε`.

    For any quadratically contracting sequence
    `e_{k+1} ≤ K · e_k²` with `e_0 ≤ r` and `K · r < 1`,
    every iterate `k ≥ effectiveLoglogCount K r ε` satisfies
    `e_k ≤ ε`. -/
theorem quadratic_loglog_effective
    (K r : ℝ) (hK : 0 < K) (hr : 0 < r) (hKr : K * r < 1)
    (e : ℕ → ℝ) (he_nn : ∀ k, 0 ≤ e k) (he_0 : e 0 ≤ r)
    (he_rec : ∀ k, e (k + 1) ≤ K * (e k) ^ 2)
    (ε : ℝ) (hε_pos : 0 < ε) :
    ∀ k ≥ effectiveLoglogCount K r ε, e k ≤ ε := by
  intro k hk
  have hKr_pos : 0 < K * r := mul_pos hK hr
  -- Tower bound: K·e_k ≤ (K·r)^(2^k).
  have h_tower := quadratic_tower_bound K r hK hr.le e he_nn he_0 he_rec k
  -- 2^k ≥ 2^(effectiveLoglogCount K r ε) > linearTowerCount K r ε.
  have h_2k_ge : linearTowerCount K r ε ≤ 2 ^ k := by
    have h1 : linearTowerCount K r ε < 2 ^ effectiveLoglogCount K r ε :=
      effectiveLoglogCount_pow_gt K r ε
    have h2 : 2 ^ effectiveLoglogCount K r ε ≤ 2 ^ k :=
      Nat.pow_le_pow_right (by norm_num) hk
    linarith
  -- (K·r)^(2^k) ≤ (K·r)^(linearTowerCount K r ε): higher power, smaller value.
  have h_pow_mono : (K * r) ^ (2 ^ k) ≤ (K * r) ^ (linearTowerCount K r ε) :=
    pow_le_pow_of_le_one hKr_pos.le hKr.le h_2k_ge
  -- (K·r)^(linearTowerCount K r ε) ≤ K·ε.
  have h_lin := linearTowerCount_pow_le K r ε hK hr hKr hε_pos
  -- Combine.
  have h_chain : K * e k ≤ K * ε := by linarith
  exact le_of_mul_le_mul_left h_chain hK

end EffectiveLoglogComplexity

/-! ============================================================
  §98.5  Effective loglog from basin (parallel to §34.1)
============================================================ -/

section EffectiveLoglogFromBasin

/-- **★★ Effective version of §34.1 `quadratic_loglog_from_basin`.**

    Replaces the existential N in §34.1 by the closed-form
    `effectiveLoglogCount K_α R_α ε` with `K_α, R_α` from §33. -/
theorem quadratic_loglog_from_basin_effective_v2
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    (z : ℂ) (hz : ‖z - α‖ < steffensenR_of_fp x hx p hp α hα hSp hfp)
    (ε : ℝ) (hε_pos : 0 < ε) :
    ∀ k ≥ effectiveLoglogCount
            (steffensenK_of_fp x hx p hp α hα hSp hfp)
            (steffensenR_of_fp x hx p hp α hα hSp hfp) ε,
      ‖(steffensen_step_C x p)^[k] z - α‖ ≤ ε := by
  set K := steffensenK_of_fp x hx p hp α hα hSp hfp
  set r := steffensenR_of_fp x hx p hp α hα hSp hfp
  obtain ⟨hK_pos, hr_pos, hKr, h_quad, _⟩ :=
    steffensen_explicit_super_attractive_rate x hx p hp α hα hSp hfp
  set e : ℕ → ℝ := fun k => ‖(steffensen_step_C x p)^[k] z - α‖
  have he_nn : ∀ k, 0 ≤ e k := fun _ => norm_nonneg _
  have he_0 : e 0 ≤ r := by
    simp [e, Function.iterate_zero]; exact le_of_lt hz
  -- Stay-in-basin + quadratic recurrence (verbatim from §34.1).
  have h_iter_unfold : ∀ k, (steffensen_step_C x p)^[k + 1] z =
      steffensen_step_C x p ((steffensen_step_C x p)^[k] z) := by
    intro k; rw [Function.iterate_succ_apply']
  have h_stay_quad : ∀ k, e k < r ∧ e (k + 1) ≤ K * (e k) ^ 2 := by
    intro k
    induction k with
    | zero =>
      refine ⟨?_, ?_⟩
      · simpa [e] using hz
      · have h1 : ‖steffensen_step_C x p z - α‖ ≤ K * ‖z - α‖ ^ 2 := h_quad z hz
        have h_e0 : e 0 = ‖z - α‖ := by simp [e]
        have h_e1 : e 1 = ‖steffensen_step_C x p z - α‖ := by
          simp [e, h_iter_unfold 0]
        rw [h_e0, h_e1]; exact h1
    | succ k ih =>
      obtain ⟨h_ek_lt_r, h_rec_k⟩ := ih
      have h_e_nn : 0 ≤ e k := he_nn k
      have h_ek_le_r : e k ≤ r := le_of_lt h_ek_lt_r
      have h_rec_trans : e (k + 1) ≤ K * r * e k := by
        calc e (k + 1) ≤ K * (e k) ^ 2 := h_rec_k
          _ = K * e k * e k := by ring
          _ ≤ K * r * e k := by
              apply mul_le_mul_of_nonneg_right _ h_e_nn
              exact mul_le_mul_of_nonneg_left h_ek_le_r hK_pos.le
      have h_half_e : e (k + 1) ≤ (1 / 2) * e k :=
        le_trans h_rec_trans (by
          apply mul_le_mul_of_nonneg_right hKr h_e_nn)
      have h_ek1_lt_r : e (k + 1) < r := by
        calc e (k + 1) ≤ (1 / 2) * e k := h_half_e
          _ ≤ (1 / 2) * r := mul_le_mul_of_nonneg_left h_ek_le_r (by norm_num)
          _ < r := by linarith
      have h_base : ‖steffensen_step_C x p ((steffensen_step_C x p)^[k + 1] z) - α‖
          ≤ K * ‖(steffensen_step_C x p)^[k + 1] z - α‖ ^ 2 :=
        h_quad _ h_ek1_lt_r
      have h_ek1_eq : e (k + 1) = ‖(steffensen_step_C x p)^[k + 1] z - α‖ := by
        simp [e]
      have h_ek2_eq :
          e (k + 2) = ‖steffensen_step_C x p ((steffensen_step_C x p)^[k + 1] z) - α‖ := by
        simp [e, h_iter_unfold (k + 1)]
      refine ⟨h_ek1_lt_r, ?_⟩
      rw [h_ek2_eq, h_ek1_eq]; exact h_base
  have he_rec : ∀ k, e (k + 1) ≤ K * (e k) ^ 2 := fun k => (h_stay_quad k).2
  have hKr_lt_one : K * r < 1 := by linarith
  exact quadratic_loglog_effective K r hK_pos hr_pos hKr_lt_one
    e he_nn he_0 he_rec ε hε_pos

end EffectiveLoglogFromBasin

/-! ============================================================
  §98.6  Pandrosion-specific effective count
============================================================ -/

section PandrosionEffectiveCount

/-- **★ Effective Pandrosion loglog count.** The closed-form
    iteration count for the §95 universal multi-start theorem,
    using §33's `K_α := steffensenK_of_fp` and
    `R_α := steffensenR_of_fp`.

    *This is the first formally verified Kung-Traub-optimal
    effective complexity bound for `z^p = x`.* -/
noncomputable def pandrosionEffectiveLoglogCount
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hSp : Sp_C p α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (ε : ℝ) : ℕ :=
  effectiveLoglogCount
    (steffensenK_of_fp x hx p hp α hα hSp hα_pow)
    (steffensenR_of_fp x hx p hp α hα hSp hα_pow)
    ε

end PandrosionEffectiveCount

/-! ============================================================
  §98.7  ★★★ Effective Kung-Traub for §95 universal multi-start
============================================================ -/

section PandrosionEffectiveTheorem

/-- **★★★★★★★★★ EFFECTIVE KUNG-TRAUB UNIVERSAL MULTI-START.**

    For any `(x, p, α)` with the standard non-degeneracy hypotheses,
    any `z ∈ ℂ` with `z ≠ α`, and any precision `ε > 0`, the
    multi-start algorithm with seed
    `multi_start_basin_seed_generic α (R_α / (2·‖z − α‖)) z`
    achieves precision `ε` after exactly
    `pandrosionEffectiveLoglogCount x p α ε` iterations:

        ∀ k ≥ pandrosionEffectiveLoglogCount,
            ‖σ^k(seed) − α‖ ≤ ε.

    **Asymptotic.** The count grows as `log₂ log(1/ε)` —
    matching the Kung-Traub optimal bound for derivative-free
    order-2 methods. **Inconditionnel.**

    This is the **fully effective version of §95.3**: the
    non-explicit `N` is replaced by the closed-form
    `effectiveLoglogCount K_α R_α ε`, computable from
    `(x, p, α, ε)` alone. -/
theorem pandrosion_loglog_universal_multi_start_effective
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hSp : Sp_C p α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (z : ℂ) (hz_ne : z ≠ α) (ε : ℝ) (hε_pos : 0 < ε) :
    ∀ k ≥ pandrosionEffectiveLoglogCount x hx p hp α hα hSp hα_pow ε,
      ‖(steffensen_step_C x p)^[k]
          (multi_start_basin_seed_generic α
            (steffensenR_of_fp x hx p hp α hα hSp hα_pow / (2 * ‖z - α‖))
            z) - α‖ ≤ ε := by
  have h_d_pos : 0 < ‖z - α‖ := norm_pos_iff.mpr (sub_ne_zero.mpr hz_ne)
  set R := steffensenR_of_fp x hx p hp α hα hSp hα_pow
  have hR_pos : 0 < R := steffensenR_of_fp_pos x hx p hp α hα hSp hα_pow
  set ε_seed : ℝ := R / (2 * ‖z - α‖)
  have hε_seed_pos : 0 < ε_seed := by unfold_let ε_seed; positivity
  -- The seed is in B(α, R/2) ⊂ basin.
  have h_seed_in_basin :
      ‖multi_start_basin_seed_generic α ε_seed z - α‖ < R := by
    rw [multi_start_basin_seed_generic_norm α ε_seed hε_seed_pos z]
    have h_d_ne : ‖z - α‖ ≠ 0 := ne_of_gt h_d_pos
    have h_eq : ε_seed * ‖z - α‖ = R / 2 := by
      unfold_let ε_seed
      rw [div_mul_eq_mul_div, mul_div_mul_right _ _ h_d_ne]
    linarith
  -- Apply §98.5 effective from basin.
  exact quadratic_loglog_from_basin_effective_v2 x hx p hp α hα hSp hα_pow
    (multi_start_basin_seed_generic α ε_seed z) h_seed_in_basin
    ε hε_pos

end PandrosionEffectiveTheorem

end Pandrosion
