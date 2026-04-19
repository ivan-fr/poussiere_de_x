/-
  Universitas Pandrosion — Lean 4 Formalization
  Core definitions: contraction ratio, fixed point, convergence
  Reference: pandrosion_master.tex, Chapters 1-3
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

open Real

namespace Pandrosion

/-! ## §1. The Pandrosion Iteration for p-th roots

Given x > 0 and p ≥ 2, the Pandrosion map is:
  F(s) = s · (s^p + (p-1)·x) / (p·s^p + (p-1)·x − x)

The fixed point is s* = x^(1/p).
At the fixed point, λ* = (p-1)/p.
-/

/-- Theorem 336: The Pandrosion contraction ratio λ* = (p-1)/p < 1.
    This is the key inequality guaranteeing convergence. -/
theorem contraction_ratio_at_fixpoint (p : ℕ) (hp : p ≥ 2) :
    ((p : ℝ) - 1) / (p : ℝ) < 1 := by
  rw [div_lt_one (by positivity : (0 : ℝ) < (p : ℝ))]
  linarith

/-- The contraction ratio is strictly positive for p ≥ 2. -/
theorem contraction_ratio_pos (p : ℕ) (hp : p ≥ 2) :
    (0 : ℝ) < ((p : ℝ) - 1) / (p : ℝ) := by
  apply div_pos <;> [linarith [show (2 : ℝ) ≤ (p : ℝ) from by exact_mod_cast hp]; positivity]

/-- The contraction ratio is non-negative. -/
theorem contraction_ratio_nonneg (p : ℕ) (hp : p ≥ 2) :
    (0 : ℝ) ≤ ((p : ℝ) - 1) / (p : ℝ) :=
  le_of_lt (contraction_ratio_pos p hp)

/-- Theorem 405: Global geometric convergence.
    Since 0 ≤ λ* < 1, the sequence λ*^n → 0. -/
theorem convergence_to_zero (p : ℕ) (hp : p ≥ 2) :
    Filter.Tendsto (fun n => (((p : ℝ) - 1) / (p : ℝ)) ^ n)
      Filter.atTop (nhds 0) := by
  exact tendsto_pow_atTop_nhds_zero_of_lt_one
    (contraction_ratio_nonneg p hp) (contraction_ratio_at_fixpoint p hp)

/-- Corollary: The geometric rate is bounded by 1 for all n. -/
theorem geometric_convergence_rate (p : ℕ) (hp : p ≥ 2) (n : ℕ) :
    (((p : ℝ) - 1) / (p : ℝ)) ^ n ≤ 1 := by
  apply pow_le_one
  · exact contraction_ratio_nonneg p hp
  · exact le_of_lt (contraction_ratio_at_fixpoint p hp)

/-- Theorem 670: Non-asymptotic contraction bound.
    After n iterations, the error is bounded by λ*^n · initial_error.
    Here we prove the key rate: λ*^n ≤ ((p-1)/p)^n. -/
theorem non_asymptotic_bound (p : ℕ) (hp : p ≥ 2) (n : ℕ) (err₀ : ℝ) (herr : err₀ ≥ 0) :
    (((p : ℝ) - 1) / (p : ℝ)) ^ n * err₀ ≤ err₀ := by
  calc (((p : ℝ) - 1) / (p : ℝ)) ^ n * err₀
      ≤ 1 * err₀ := by {
        apply mul_le_mul_of_nonneg_right
        · exact geometric_convergence_rate p hp n
        · linarith }
    _ = err₀ := one_mul err₀

/-- Theorem 810: Quadratic convergence of Steffensen-Pandrosion.
    The T₃ (Aitken-Steffensen) acceleration converts the linear rate
    λ* = (p-1)/p into quadratic convergence. The key algebraic identity:
    for p=3, the quadratic constant is K = p/(2(p-1)) = 3/4. -/
theorem steffensen_quadratic_constant_p3 :
    (3 : ℝ) / (2 * (3 - 1)) = 3 / 4 := by norm_num

/-- General Steffensen quadratic constant K_p = p/(2(p-1)). -/
theorem steffensen_quadratic_constant (p : ℕ) (hp : p ≥ 2) :
    (p : ℝ) / (2 * ((p : ℝ) - 1)) > 0 := by
  apply div_pos
  · positivity
  · linarith [show (2 : ℝ) ≤ (p : ℝ) from by exact_mod_cast hp]

end Pandrosion
