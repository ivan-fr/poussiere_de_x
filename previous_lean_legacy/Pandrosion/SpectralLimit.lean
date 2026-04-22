/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP XI: SPECTRAL LIMIT D_d → -ln(2)

  Architecture:
  1. Define D_p as the spectral sum
  2. Prove D_p < 0 for p ≥ 2 (sum_lt_sum + cos bounds)
  3. Closed form via product formula + tangent symmetry [classical input]
  4. Prove D_p → -ln(2) from closed form + log(n)/n → 0

  Reference: pandrosion_master.tex, §77 (Spectral Limit)
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.MeasureTheory.Integral.IntervalIntegral
import Mathlib.Tactic
import Pandrosion.FixedPointGeneral

open Finset BigOperators Real MeasureTheory Set Filter

namespace Pandrosion

/-! ## §77. The Spectral Descent Coefficient -/

noncomputable def D (p : ℕ) : ℝ :=
  (1 / (p : ℝ)) * ∑ k in range p, Real.log (Real.cos ((k : ℝ) * π / (2 * (p : ℝ))))

/-! ## §78. Classical Analytic Input: Integral Identity -/

variable (integral_log_cos_eq :
    ∫ t in (0 : ℝ)..1, Real.log (Real.cos (π * t / 2)) = -Real.log 2)

/-! ## §79. Auxiliary: angle bounds and cosine properties -/

theorem angle_nn (p k : ℕ) (_hk : k < p) :
    0 ≤ (k : ℝ) * π / (2 * (p : ℝ)) := by positivity

theorem angle_lt_half_pi (p k : ℕ) (hp : p ≥ 1) (hk : k < p) :
    (k : ℝ) * π / (2 * (p : ℝ)) < π / 2 := by
  rw [div_lt_div_iff (by positivity : (0 : ℝ) < 2 * ↑p) two_pos]
  have : (k : ℝ) < (p : ℝ) := Nat.cast_lt.mpr hk
  nlinarith [Real.pi_pos]

theorem cos_angle_is_pos (p k : ℕ) (hp : p ≥ 1) (hk : k < p) :
    0 < Real.cos ((k : ℝ) * π / (2 * (p : ℝ))) := by
  apply Real.cos_pos_of_mem_Ioo
  constructor
  · have := angle_nn p k hk; linarith [Real.pi_pos]
  · exact angle_lt_half_pi p k hp hk

theorem cos_angle_lt_unit (p k : ℕ) (hp : p ≥ 1) (hk_pos : 0 < k) (hk : k < p) :
    Real.cos ((k : ℝ) * π / (2 * (p : ℝ))) < 1 := by
  have hθ_pos : (0 : ℝ) < (k : ℝ) * π / (2 * (p : ℝ)) := by positivity
  have hθ_le_pi : (k : ℝ) * π / (2 * (p : ℝ)) ≤ π := by
    rw [div_le_iff (by positivity : (0 : ℝ) < 2 * ↑p)]
    have : (k : ℝ) < (p : ℝ) := Nat.cast_lt.mpr hk
    nlinarith [Real.pi_pos]
  calc Real.cos ((k : ℝ) * π / (2 * (p : ℝ)))
      < Real.cos 0 := by
        apply Real.strictAntiOn_cos
        · exact ⟨le_refl _, le_of_lt Real.pi_pos⟩
        · exact ⟨le_of_lt hθ_pos, hθ_le_pi⟩
        · exact hθ_pos
    _ = 1 := Real.cos_zero

/-! ## §80. D_p < 0 for p ≥ 2 -/

theorem D_neg (p : ℕ) (hp : p ≥ 2) : D p < 0 := by
  unfold D
  apply mul_neg_of_pos_of_neg
  · positivity
  · have hle : ∀ k ∈ Finset.range p,
        Real.log (Real.cos ((k : ℝ) * π / (2 * (p : ℝ)))) ≤ 0 := by
      intro k hk
      apply Real.log_nonpos
      · exact le_of_lt (cos_angle_is_pos p k (by omega) (Finset.mem_range.mp hk))
      · exact Real.cos_le_one _
    have hlt : Real.log (Real.cos (((1 : ℕ) : ℝ) * π / (2 * (p : ℝ)))) < 0 := by
      apply Real.log_neg
      · exact cos_angle_is_pos p 1 (by omega) (by omega)
      · exact cos_angle_lt_unit p 1 (by omega) (by omega) (by omega)
    calc ∑ k in Finset.range p, Real.log (Real.cos ((k : ℝ) * π / (2 * (p : ℝ))))
        < ∑ _k in Finset.range p, (0 : ℝ) :=
          Finset.sum_lt_sum hle ⟨1, Finset.mem_range.mpr (by omega),
            by push_cast at hlt ⊢; exact hlt⟩
      _ = 0 := by simp

/-! ## §81. Classical Analytic Input: Closed Form -/

variable (D_eq_closed : ∀ (p : ℕ) (_hp : p ≥ 2),
    D p = (Real.log (2 * ↑p) - (2 * ↑p - 1) * Real.log 2) / (2 * ↑p))

/-! ## §82. D_p → -ln(2) -/

/-- **D_p converges to -ln(2) as p → ∞.** -/
theorem D_tendsto_neg_log2 :
    Tendsto D atTop (nhds (-Real.log 2)) := by
  -- Step 1: For n ≥ 2, D n equals the closed form which we rewrite as
  -- -log 2 + log(2n)/(2n) + log(2)/(2n)
  -- Step 2: Show this → -log 2 + 0 + 0 = -log 2
  -- Use Tendsto.congr': if f =ᶠ g and Tendsto g, then Tendsto f
  have hcf : Tendsto
      (fun n : ℕ => Real.log (2 * (n : ℝ)) / (2 * (n : ℝ)) +
        (- Real.log 2) + Real.log 2 / (2 * (n : ℝ)))
      atTop (nhds (-Real.log 2)) := by
    -- -log 2 = 0 + (-log 2) + 0
    conv_rhs => rw [show (-Real.log 2 : ℝ) = 0 + (-Real.log 2) + 0 from by ring]
    apply Filter.Tendsto.add
    · apply Filter.Tendsto.add
      · -- log(2n)/(2n) → 0
        have hlog : Tendsto (fun x : ℝ => Real.log x / x) atTop (nhds 0) :=
          Real.isLittleO_log_id_atTop.tendsto_div_nhds_zero
        have h2n : Tendsto (fun n : ℕ => (2 : ℝ) * (n : ℝ)) atTop atTop := by
          have := @tendsto_nat_cast_atTop_atTop ℝ _ _
          exact (this.atTop_mul_const (by norm_num : (0 : ℝ) < 2)).congr
            (fun n => by ring)
        exact (hlog.comp h2n).congr (fun _ => rfl)
      · -- const -log 2 → -log 2
        exact tendsto_const_nhds
    · -- log(2)/(2n) → 0
      -- log(2)/(2n) = (log 2 / 2) * (1/n) → 0
      have : Tendsto (fun n : ℕ => (n : ℝ)⁻¹) atTop (nhds 0) :=
        tendsto_inverse_atTop_nhds_zero_nat
      have h0 : Tendsto (fun n : ℕ => Real.log 2 / (2 * (n : ℝ))) atTop (nhds 0) := by
        have h1 := tendsto_const_div_atTop_nhds_zero_nat (Real.log 2 / 2)
        apply h1.congr
        intro n
        by_cases hn : (n : ℝ) = 0
        · simp [hn]
        · field_simp
      exact h0
  -- Now apply congr': D =ᶠ the formula above
  apply hcf.congr'
  filter_upwards [eventually_ge_atTop 2] with n hn
  rw [D_eq_closed n hn]
  have h2n : (2 : ℝ) * ↑n ≠ 0 := by positivity
  field_simp
  ring

/-! ## §83. Explicit Computations -/

theorem D_two_eq : D 2 = (1 / 2) * Real.log (Real.cos (π / 4)) := by
  unfold D
  simp only [Finset.sum_range_succ, Finset.sum_range_zero]
  norm_num

end Pandrosion
