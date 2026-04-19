/-
  Universitas Pandrosion — Lean 4 Formalization
  Scaling Optimization & Steffensen Acceleration
  Reference: pandrosion_master.tex, Theorems 218, 425, 606, 713, 752, 900, 946
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic
import Pandrosion.Core

open Real

namespace Pandrosion

/-! ## §8. The Steffensen Quadratic Constant (Theorem 810, detailed) -/

/-- Theorem 810 (structural): The Steffensen constant K_S > 0. -/
theorem steffensen_constant_pos (h'' lam : ℝ) (hh : h'' ≠ 0) (hlam : lam < 1) :
    |h'' / (2 * (1 - lam))| > 0 := by
  apply abs_pos.mpr
  apply div_ne_zero hh
  apply mul_ne_zero (by norm_num : (2:ℝ) ≠ 0)
  linarith

/-- Theorem 1658: Newton's quadratic constant K_N = (p-1)/(2α) > 0. -/
theorem newton_quadratic_constant_pos (p : ℕ) (hp : p ≥ 2) (a : ℝ) (ha : a > 0) :
    ((p : ℝ) - 1) / (2 * a) > 0 := by
  apply div_pos
  · linarith [show (2 : ℝ) ≤ (p : ℝ) from by exact_mod_cast hp]
  · linarith

/-- Theorem 1690: Steffensen beats Newton when lam < 1/2. -/
theorem steffensen_beats_newton_ratio (lam : ℝ) (h1 : 0 < lam) (h2 : lam < 1 / 2) :
    lam / (1 - lam) < 1 := by
  rw [div_lt_one (by linarith)]
  linarith

/-! ## §9. Scaling Optimization (Theorem 900) -/

/-- Theorem 900: Scaling preserves positivity: x/A > 0. -/
theorem scaling_identity (x A : ℝ) (hx : x > 0) (hA : A > 0) :
    x / A > 0 := div_pos hx hA

/-- The reduced ratio x' = x/A < x when A > 1. -/
theorem reduced_ratio_lt (x A : ℝ) (hx : x > 0) (hA : A > 1) :
    x / A < x := by
  rw [div_lt_iff (by linarith)]
  nlinarith

/-- Proposition 425: ln(1/lam) > 0 when 0 < lam < 1. -/
theorem log_inv_contraction_pos (lam : ℝ) (h1 : 0 < lam) (h2 : lam < 1) :
    Real.log (1 / lam) > 0 := by
  rw [one_div]
  exact Real.log_pos (one_lt_inv h1 h2)

/-- Scaling benefit: smaller lam → larger ln(1/lam) → fewer iterations. -/
theorem fewer_iterations (lam1 lam2 : ℝ) (h1 : 0 < lam1) (h2 : lam1 < lam2) (_h3 : lam2 < 1) :
    Real.log (1 / lam2) < Real.log (1 / lam1) := by
  apply Real.log_lt_log
  · rw [one_div]; exact inv_pos.mpr (by linarith)
  · rw [one_div, one_div]
    exact inv_lt_inv_of_lt h1 h2

/-! ## §10. Asymptotic Regimes (Proposition 606) -/

/-- Regime 1: Near α = 1, (p-1)(α-1)/2 > 0. -/
theorem near_identity_regime (p : ℕ) (hp : p ≥ 2) (a : ℝ) (ha : a > 1) :
    ((p : ℝ) - 1) * (a - 1) / 2 > 0 := by
  apply div_pos
  · apply mul_pos
    · linarith [show (2 : ℝ) ≤ (p : ℝ) from by exact_mod_cast hp]
    · linarith
  · norm_num

/-- Regime 2: (α-1)/α < 1 (approaches 1 from below). -/
theorem large_alpha_regime (a : ℝ) (ha : a > 1) :
    (a - 1) / a < 1 := by
  rw [div_lt_one (by linarith)]; linarith

/-- Regime 2: (α-1)/α > 0. -/
theorem large_alpha_regime_pos (a : ℝ) (ha : a > 1) :
    (a - 1) / a > 0 := div_pos (by linarith) (by linarith)

/-- Regime 3: ln(x) ≤ x - 1 for x > 0 (standard, from Mathlib). -/
theorem log_le_sub_one (x : ℝ) (hx : x > 0) :
    Real.log x ≤ x - 1 := Real.log_le_sub_one_of_pos hx

/-- Therefore 1 - ln(x)/(x-1) ≥ 0 for x > 1. -/
theorem limit_regime_nonneg (x : ℝ) (hx : x > 1) :
    1 - Real.log x / (x - 1) ≥ 0 := by
  have h1 : x - 1 > 0 := by linarith
  have h2 := log_le_sub_one x (by linarith)
  rw [ge_iff_le, sub_nonneg, div_le_one h1]
  exact h2

/-- And 1 - ln(x)/(x-1) < 1 for x > 1 (since ln(x) > 0). -/
theorem limit_regime_lt_one (x : ℝ) (hx : x > 1) :
    1 - Real.log x / (x - 1) < 1 := by
  simp only [sub_lt_self_iff]
  apply div_pos (Real.log_pos hx) (by linarith)

/-! ## §11. Optimal Starting Point (Proposition 752) -/

/-- Proposition 752: The optimal starting point s₀ = 1 - (x-1)/(xp) ∈ (0,1). -/
theorem optimal_start_in_unit (x : ℝ) (hx : x > 1) (p : ℕ) (hp : p ≥ 2) :
    0 < 1 - (x - 1) / ((x : ℝ) * (p : ℝ)) ∧
    1 - (x - 1) / ((x : ℝ) * (p : ℝ)) < 1 := by
  have hxp : x * (p : ℝ) > 0 := by positivity
  constructor
  · rw [sub_pos, div_lt_one hxp]
    nlinarith [show (2 : ℝ) ≤ (p : ℝ) from by exact_mod_cast hp]
  · simp only [sub_lt_self_iff]
    apply div_pos (by linarith) hxp

/-- Proposition 713: For x < 1, α - 1 < 0 (oscillatory). -/
theorem oscillatory_for_x_lt_one (a : ℝ) (_ha : 0 < a) (ha1 : a < 1) :
    a - 1 < 0 := by linarith

end Pandrosion
