/-
  Universitas Pandrosion — Lean 4 Formalization
  BASE ITERATION COMPLEXITY BOUND

  Provides the formal logarithmic bound for linear contraction
  required to establish the O(log(1/ε) / log(1/c)) step count.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

open Real

namespace Pandrosion

/-! ## §230. Linear Contraction Step Bound -/

/-- **Strict logarithmic bound for linear contraction.**
    If the sequence contracts linearly to zero, the necessary steps to reach ε
    is strictly bounded by the logarithmic ratio. -/
theorem base_complexity_bound (c e₀ ε : ℝ)
    (hc_pos : 0 < c) (hc_lt : c < 1)
    (he_pos : 0 < e₀) (hε_pos : 0 < ε) (n : ℕ)
    (hn : Real.log (e₀ / ε) / Real.log (1 / c) ≤ (n : ℝ)) :
    c ^ n * e₀ ≤ ε := by
  -- Since 0 < c < 1, 1/c > 1, so log(1/c) > 0
  have hc_inv_one : 1 < 1 / c := one_lt_one_div hc_pos hc_lt
  have hlog_pos : 0 < Real.log (1 / c) := Real.log_pos hc_inv_one

  -- Multiply by log(1/c)
  have h1 : Real.log (e₀ / ε) ≤ (n : ℝ) * Real.log (1 / c) := by
    -- hn : Real.log (e₀ / ε) / Real.log (1 / c) ≤ n
    -- We want to show Real.log (e₀ / ε) ≤ n * Real.log (1 / c)
    have step1 : Real.log (e₀ / ε) / Real.log (1 / c) * Real.log (1 / c) ≤ (n : ℝ) * Real.log (1 / c) :=
      mul_le_mul_of_nonneg_right hn (le_of_lt hlog_pos)
    have step2 : Real.log (e₀ / ε) / Real.log (1 / c) * Real.log (1 / c) = Real.log (e₀ / ε) :=
      div_mul_cancel₀ _ (ne_of_gt hlog_pos)
    rwa [step2] at step1

  -- log(1/c) = log(c⁻¹) = -log(c)
  have hlog_inv : Real.log (1 / c) = -Real.log c := by
    rw [one_div, Real.log_inv c]
  rw [hlog_inv] at h1

  -- log(e₀/ε) = log(e₀) - log(ε)
  have hlog_div : Real.log (e₀ / ε) = Real.log e₀ - Real.log ε :=
    Real.log_div he_pos.ne' hε_pos.ne'
  rw [hlog_div] at h1

  -- Rearrange: n * log(c) + log(e₀) ≤ log(ε)
  have h2 : (n : ℝ) * Real.log c + Real.log e₀ ≤ Real.log ε := by
    linarith

  -- Real.exp(a + b) = exp(a) * exp(b)
  have h3 : Real.exp ((n : ℝ) * Real.log c + Real.log e₀) ≤ Real.exp (Real.log ε) :=
    Real.exp_le_exp.mpr h2

  -- exp(log(ε)) = ε
  have hexp_eps : Real.exp (Real.log ε) = ε := Real.exp_log hε_pos
  rw [hexp_eps] at h3

  -- exp(a + b) = exp(a) * exp(b)
  have hexp_add : Real.exp ((n : ℝ) * Real.log c + Real.log e₀) =
                  Real.exp ((n : ℝ) * Real.log c) * Real.exp (Real.log e₀) :=
    Real.exp_add ((n : ℝ) * Real.log c) (Real.log e₀)
  rw [hexp_add] at h3

  -- exp(log(e₀)) = e₀
  have hexp_e0 : Real.exp (Real.log e₀) = e₀ := Real.exp_log he_pos
  rw [hexp_e0] at h3

  -- exp(n * log(c)) = c^n
  have hexp_n_log : Real.exp ((n : ℝ) * Real.log c) = c ^ n := by
    -- Using rpow definition: c ^ (n:ℝ) = exp(log c * n)
    have hc_rpow : c ^ (n : ℝ) = Real.exp (Real.log c * (n : ℝ)) :=
      Real.rpow_def_of_pos hc_pos (n : ℝ)
    have h_mul_comm : Real.log c * (n : ℝ) = (n : ℝ) * Real.log c := mul_comm _ _
    rw [h_mul_comm] at hc_rpow
    rw [← hc_rpow]
    exact Real.rpow_nat_cast c n

  rw [hexp_n_log] at h3
  exact h3

end Pandrosion
