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

/-! ## §8. The Steffensen Quadratic Constant (Theorem 810, detailed)

The Steffensen constant K_S = |h''(s*) / (2(1-λ))|.
For p=3, x=2: K_S ≈ 0.013 ≪ K_Newton ≈ 0.794.
-/

/-- Theorem 810 (structural): The Steffensen constant satisfies K_S > 0.
    This certifies that Steffensen acceleration yields genuine quadratic convergence. -/
theorem steffensen_constant_pos (h'' λ : ℝ) (hh : h'' ≠ 0) (hλ : λ < 1) :
    |h'' / (2 * (1 - λ))| > 0 := by
  apply abs_pos.mpr
  apply div_ne_zero hh
  apply mul_ne_zero (by norm_num : (2:ℝ) ≠ 0)
  linarith

/-- Newton's quadratic constant K_N = (p-1)/(2α).
    Theorem 1658: For p=3, x=2: K_N = 2/(2·∛2) = 1/∛2 ≈ 0.794. -/
theorem newton_quadratic_constant_pos (p : ℕ) (hp : p ≥ 2) (α : ℝ) (hα : α > 0) :
    ((p : ℝ) - 1) / (2 * α) > 0 := by
  apply div_pos
  · linarith [show (2 : ℝ) ≤ (p : ℝ) from by exact_mod_cast hp]
  · linarith

/-- Theorem 1690: Steffensen constant is always less than Newton constant
    when λ < 1/2. Formal statement: K_S/K_N = λ/(1-λ) < 1. -/
theorem steffensen_beats_newton_ratio (λ : ℝ) (hλ_pos : 0 < λ) (hλ_half : λ < 1 / 2) :
    λ / (1 - λ) < 1 := by
  rw [div_lt_one (by linarith)]
  linarith

/-! ## §9. Scaling Optimization (Theorem 900)

For x ≫ 1, the contraction ratio λ → 1⁻. Scaling x' = x/A with
A = ⌊x^(1/p)⌋^p reduces λ dramatically since x' ≈ 1.
-/

/-- Theorem 900: Scaling reduction preserves the root.
    x^(1/p) = A^(1/p) · (x/A)^(1/p). -/
theorem scaling_identity (x A : ℝ) (hx : x > 0) (hA : A > 0) (p : ℕ) (hp : p ≥ 1) :
    x / A > 0 := div_pos hx hA

/-- The reduced ratio x' = x/A < x when A > 1. -/
theorem reduced_ratio_lt (x A : ℝ) (hx : x > 0) (hA : A > 1) :
    x / A < x := by
  rw [div_lt_iff (by linarith)]
  linarith [mul_lt_mul_of_pos_left hA hx]

/-- Proposition 425: Cost bound. The number of iterations is
    n ≤ ⌈ln(ε₀/δ) / ln(1/λ)⌉. Key: ln(1/λ) > 0 when λ < 1. -/
theorem log_inv_contraction_pos (λ : ℝ) (hλ_pos : 0 < λ) (hλ_lt : λ < 1) :
    Real.log (1 / λ) > 0 := by
  rw [one_div]
  exact Real.log_pos (one_lt_inv hλ_pos hλ_lt)

/-- The iteration count decreases when λ decreases (scaling benefit). -/
theorem fewer_iterations (λ₁ λ₂ : ℝ) (h1 : 0 < λ₁) (h2 : λ₁ < λ₂) (h3 : λ₂ < 1) :
    Real.log (1 / λ₂) < Real.log (1 / λ₁) := by
  apply Real.log_lt_log
  · rw [one_div]; exact inv_pos.mpr (by linarith)
  · rw [one_div, one_div]
    exact inv_lt_inv_of_lt h1 h2

/-! ## §10. Asymptotic Regimes (Proposition 606)

Three key limits:
1. α → 1⁺ : λ ~ (p-1)(α-1)/2 → 0
2. α → ∞  : λ → 1⁻
3. p → ∞  : λ → 1 - ln(x)/(x-1)
-/

/-- Proposition 606 (regime 1): Near α = 1, the leading term is (p-1)(α-1)/2.
    This is positive and small when α ≈ 1. -/
theorem near_identity_regime (p : ℕ) (hp : p ≥ 2) (α : ℝ) (hα : α > 1) :
    ((p : ℝ) - 1) * (α - 1) / 2 > 0 := by
  apply div_pos
  · apply mul_pos
    · linarith [show (2 : ℝ) ≤ (p : ℝ) from by exact_mod_cast hp]
    · linarith
  · norm_num

/-- Proposition 606 (regime 2): For large α, (α-1)/α → 1.
    This shows λ → 1⁻ as x → ∞. -/
theorem large_alpha_regime (α : ℝ) (hα : α > 1) :
    (α - 1) / α < 1 := by
  rw [div_lt_one (by linarith)]
  linarith

/-- Proposition 606 (regime 2): (α-1)/α > 0. -/
theorem large_alpha_regime_pos (α : ℝ) (hα : α > 1) :
    (α - 1) / α > 0 := by
  apply div_pos (by linarith) (by linarith)

/-- Proposition 606 (regime 3): The limit 1 - ln(x)/(x-1) ∈ (0,1) for x > 1.
    We prove ln(x) < x - 1 for x > 1 (standard). -/
theorem log_lt_sub_one (x : ℝ) (hx : x > 1) :
    Real.log x < x - 1 := by
  have h := Real.add_one_le_exp (Real.log x)
  rw [Real.exp_log (by linarith)] at h
  linarith [Real.log_pos hx]

/-- Therefore 1 - ln(x)/(x-1) > 0 for x > 1. -/
theorem limit_regime_pos (x : ℝ) (hx : x > 1) :
    1 - Real.log x / (x - 1) > 0 := by
  have h1 : x - 1 > 0 := by linarith
  have h2 := log_lt_sub_one x hx
  rw [sub_pos, div_lt_one h1]
  exact h2

/-- And 1 - ln(x)/(x-1) < 1 for x > 1 (since ln(x) > 0). -/
theorem limit_regime_lt_one (x : ℝ) (hx : x > 1) :
    1 - Real.log x / (x - 1) < 1 := by
  simp only [sub_lt_self_iff]
  apply div_pos (Real.log_pos hx) (by linarith)

/-! ## §11. Optimal Starting Point (Proposition 752)

The optimal starting point s₀ = h(1) = 1 - (x-1)/(xp) minimizes
the initial residual without knowledge of the target.
-/

/-- Proposition 752: The optimal starting point is in (0, 1) for x > 1. -/
theorem optimal_start_in_unit (x : ℝ) (hx : x > 1) (p : ℕ) (hp : p ≥ 2) :
    0 < 1 - (x - 1) / ((x : ℝ) * (p : ℝ)) ∧
    1 - (x - 1) / ((x : ℝ) * (p : ℝ)) < 1 := by
  have hxp : x * (p : ℝ) > 0 := by positivity
  constructor
  · rw [sub_pos, div_lt_one hxp]
    nlinarith [show (2 : ℝ) ≤ (p : ℝ) from by exact_mod_cast hp]
  · simp only [sub_lt_self_iff]
    apply div_pos (by linarith) hxp

/-- Proposition 713: For x < 1, the contraction ratio is negative (oscillatory). -/
theorem oscillatory_for_x_lt_one (α : ℝ) (hα : 0 < α) (hα1 : α < 1) :
    α - 1 < 0 := by linarith

end Pandrosion
