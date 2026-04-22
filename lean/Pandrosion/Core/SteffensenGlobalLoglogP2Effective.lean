/-
  Universitas Pandrosion — §47. **Effective global log-log
  complexity at `p = 2` on ℂ.**

  Lifts the remaining existential in
  `steffensen_global_loglog_p2_complex_unconditional` (§42.1):
  the inner "tail" `N_tail(ε)` from §34.1 `quadratic_loglog_from_basin`
  is replaced by the **closed-form** `quadraticLoglogN` from §44.

  Result:
      `steffensen_p2_solves_loglog_complex_effective` — for
      Lebesgue-a.e. `z₀ ∈ ℂ`, and every `ε > 0`, every iteration
      count `k ≥ k₀(z₀) + quadraticLoglogN(K_α, R_α, ε)` satisfies
      `‖σ^k z₀ - γ_s‖ ≤ ε` for some `s ∈ Fin 2`.

  Here:
    • `k₀(z₀)` is the §32.1 McMullen a.e. entry time, discharged
      unconditionally via §39.8.
    • `K_α, R_α := steffensenK_of_fp, steffensenR_of_fp` are the
      §33 basin-explicit constants.
    • `quadraticLoglogN K R ε := ⌈log(K·ε) / log(K·R)⌉ + 1`
      is the §44 closed-form effective `N_tail`.

  All existentials in the main end-to-end complexity chain are
  now closed-form: the only `∃` left is `k₀(z₀)` (McMullen
  waiting time), which is inherent (z₀-specific).

  Contents.

    §47.1  `quadraticLoglogN_from_basin` — explicit tail bound.

    §47.2  `quadratic_loglog_from_basin_effective` — effective
           version of §34.1.

    §47.3  `steffensen_p2_solves_loglog_complex_effective` —
           effective global a.e. theorem (p = 2, complex).
-/

import Pandrosion.Core.QuadraticComplexityEffective
import Pandrosion.Core.SteffensenGlobalLoglogP2Unconditional
import Pandrosion.Core.ComplexMcMullenP2Unconditional

namespace Pandrosion

open MeasureTheory Complex Filter Topology

/-! ============================================================
  §47.1  Explicit `N_tail` for the basin tail
============================================================ -/

section EffectiveBasinTailC

/-- **Explicit `N_tail` for `quadratic_loglog_from_basin`.**
    Closed-form version of §34.1's existential witness, using the
    §44 `quadraticLoglogN`. -/
noncomputable def quadraticLoglogN_from_basin
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    (ε : ℝ) : ℕ :=
  quadraticLoglogN (steffensenK_of_fp x hx p hp α hα hSp hfp)
                   (steffensenR_of_fp x hx p hp α hα hSp hfp) ε

/-- **★★★ Effective `quadratic_loglog_from_basin`.**
    For every `z` within the §33 explicit basin `B(α, R_α)` and
    every `ε > 0`, every iterate count
        `k ≥ quadraticLoglogN_from_basin(x, p, α, ε)`
    satisfies `‖σ^k z − α‖ ≤ ε`. -/
theorem quadratic_loglog_from_basin_effective
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    (z : ℂ) (hz : ‖z - α‖ < steffensenR_of_fp x hx p hp α hα hSp hfp)
    (ε : ℝ) (hε_pos : 0 < ε) :
    ∀ k ≥ quadraticLoglogN_from_basin x hx p hp α hα hSp hfp ε,
      ‖(steffensen_step_C x p)^[k] z - α‖ ≤ ε := by
  set K := steffensenK_of_fp x hx p hp α hα hSp hfp
  set r := steffensenR_of_fp x hx p hp α hα hSp hfp
  obtain ⟨hK_pos, hr_pos, hKr, h_quad, _⟩ :=
    steffensen_explicit_super_attractive_rate x hx p hp α hα hSp hfp
  set e : ℕ → ℝ := fun k => ‖(steffensen_step_C x p)^[k] z - α‖
  have he_nn : ∀ k, 0 ≤ e k := fun _ => norm_nonneg _
  have he_0 : e 0 ≤ r := by simp [e, Function.iterate_zero]; exact le_of_lt hz
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
  -- Call the effective §44 theorem.
  exact quadratic_loglog_complexity_effective K r hK_pos (le_of_lt hr_pos)
    hKr_lt_one e he_nn he_0 he_rec ε hε_pos

end EffectiveBasinTailC

/-! ============================================================
  §47.2  Effective global log-log (p = 2, ℂ, unconditional)
============================================================ -/

section EffectiveGlobalLogLogC

/-- **★★★★★★ Effective global log-log complexity for `p = 2` on ℂ,
    fully unconditional.**

    For every `x ∈ ℂ \ {0}` (with `α ≠ 0, ±1` and `α² = 1/x`), and
    every `Sp_C 2`-non-degeneracy at the cyclotomic anchors, for
    Lebesgue-almost every `z₀ ∈ ℂ`, there exists a McMullen entry
    time `k₀(z₀)` and an anchor `s ∈ Fin 2` such that for every
    `ε > 0`, every iteration count
        `k ≥ k₀(z₀) + quadraticLoglogN_from_basin(x, 2, γ_s, ε)`
    satisfies `‖σ^k z₀ − γ_s‖ ≤ ε`.

    Here `quadraticLoglogN_from_basin = O(log log 1/ε)` is
    **closed-form** via §47.1. The only residual existential is
    the McMullen waiting time `k₀(z₀)` — inherent and `z₀`-specific,
    existing unconditionally via §39.8. -/
theorem steffensen_p2_solves_loglog_complex_effective
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (hSp : ∀ s : Fin 2, Sp_C 2 (cycAnchor α 2 s) ≠ 0) :
    ∀ᵐ z₀ : ℂ ∂volume, ∃ s : Fin 2, ∃ k₀ : ℕ,
      ∀ ε > 0, ∀ k ≥ k₀ + quadraticLoglogN_from_basin x hx_ne 2 (by norm_num)
                       (cycAnchor α 2 s) (cycAnchor_ne_zero hα_ne_zero 2 s)
                       (hSp s)
                       (by rw [cycAnchor_pow α (by norm_num : 1 ≤ 2) s]; exact hα_pow)
                       ε,
        ‖(steffensen_step_C x 2)^[k] z₀ - cycAnchor α 2 s‖ ≤ ε := by
  -- McMullen discharge via §39.8.
  have hMcM : McMullenAEEntry 2 x α :=
    mcmullen_p2_complex_unconditional x hx_ne α hα_ne_zero hα_ne_neg_one
      hα_ne_one hα_pow
  -- Reuse the §34.2 a.e. entry into explicit basins.
  set r_fn : Fin 2 → ℝ := fun s =>
    steffensenR_of_fp x hx_ne 2 (by norm_num) (cycAnchor α 2 s)
      (cycAnchor_ne_zero hα_ne_zero 2 s) (hSp s)
      (by rw [cycAnchor_pow α (by norm_num : 1 ≤ 2) s]; exact hα_pow)
  have hr_fn_pos : ∀ s, 0 < r_fn s := fun s =>
    steffensenR_of_fp_pos x hx_ne 2 (by norm_num) (cycAnchor α 2 s)
      (cycAnchor_ne_zero hα_ne_zero 2 s) (hSp s) _
  have h_entry := hMcM r_fn hr_fn_pos
  filter_upwards [h_entry] with z₀ hz₀
  obtain ⟨s, k₀, h_enter⟩ := hz₀
  refine ⟨s, k₀, ?_⟩
  intro ε hε_pos k hk
  -- Apply §47.1 at entry point.
  have hfp_s : (cycAnchor α 2 s) ^ 2 = 1 / x := by
    rw [cycAnchor_pow α (by norm_num : 1 ≤ 2) s]; exact hα_pow
  have h_tail :=
    quadratic_loglog_from_basin_effective x hx_ne 2 (by norm_num) (cycAnchor α 2 s)
      (cycAnchor_ne_zero hα_ne_zero 2 s) (hSp s) hfp_s
      ((steffensen_step_C x 2)^[k₀] z₀) h_enter ε hε_pos
  set Ntail : ℕ := quadraticLoglogN_from_basin x hx_ne 2 (by norm_num) (cycAnchor α 2 s)
    (cycAnchor_ne_zero hα_ne_zero 2 s) (hSp s) hfp_s ε
  have h_iter_split :
      (steffensen_step_C x 2)^[k] z₀ =
        (steffensen_step_C x 2)^[k - k₀] ((steffensen_step_C x 2)^[k₀] z₀) := by
    rw [← Function.iterate_add_apply]
    congr 1; omega
  rw [h_iter_split]
  have h_tail_arg : Ntail ≤ k - k₀ := by omega
  exact h_tail _ h_tail_arg

end EffectiveGlobalLogLogC

end Pandrosion
