/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP XXII: THE DIFFERENTIAL ATTRACTION THEOREM
  
  This module proves that the continuous mathematical derivative of the
  Generalized Pandrosion map, evaluated at its fixed point r, is
  precisely the predicted snail constant (-1/5). This validates the
  local contraction geometry using Formal Calculus.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Tactic

open Real

namespace Pandrosion

/-! ## §126. Continuous Differential Profile -/

/-- **The Differential Attraction Constant.**
    Evaluates the explicit derivative of the Generalized iteration map
    and verifies that the local linear multiplier is strictly -1/5. -/
theorem pandrosion_contraction_ratio_3 (x r : ℝ) (h_root : r^3 = x) (hr : r ≠ 0) :
    HasDerivAt (fun s => (s^4 + 4 * x * s) / (3 * s^3 + 2 * x)) (-1 / 5) r := by
  have h_den : 3 * r^3 + 2 * x ≠ 0 := by
    rw [← h_root]
    have : 3 * r ^ 3 + 2 * r ^ 3 = 5 * r ^ 3 := by ring
    rw [this]
    exact mul_ne_zero (by norm_num) (pow_ne_zero 3 hr)
    
  have hF : HasDerivAt (fun s => s^4 + 4 * x * s) (4 * r^3 + 4 * x) r := by
    have h1 : HasDerivAt (fun s => s^4) (4 * r^3) r := hasDerivAt_pow 4 r
    have h2 : HasDerivAt (fun s => 4 * x * s) (4 * x * 1) r := HasDerivAt.const_mul (4 * x) (hasDerivAt_id r)
    have h3 := HasDerivAt.add h1 h2
    have h_eq : 4 * r ^ 3 + 4 * x * 1 = 4 * r ^ 3 + 4 * x := by ring
    rw [← h_eq]
    exact h3

  have hH : HasDerivAt (fun s => 3 * s^3 + 2 * x) (9 * r^2) r := by
    have h1 : HasDerivAt (fun s => 3 * s^3) (3 * (3 * r^2)) r := HasDerivAt.const_mul 3 (hasDerivAt_pow 3 r)
    have h2 : HasDerivAt (fun s => 2 * x) 0 r := hasDerivAt_const r (2 * x)
    have h3 := HasDerivAt.add h1 h2
    have h_eq : 3 * (3 * r ^ 2) + 0 = 9 * r ^ 2 := by ring
    rw [← h_eq]
    exact h3

  have h_div := HasDerivAt.div hF hH h_den
  
  have h_eq : (4 * r ^ 3 + 4 * x) * (3 * r ^ 3 + 2 * x) - (r ^ 4 + 4 * x * r) * (9 * r ^ 2) = - (r ^ 6) * 5 := by
    rw [← h_root]
    ring
    
  have h_den_sq : (3 * r ^ 3 + 2 * x) ^ 2 = 25 * r ^ 6 := by
    rw [← h_root]
    ring
    
  have h_final : ((4 * r ^ 3 + 4 * x) * (3 * r ^ 3 + 2 * x) - (r ^ 4 + 4 * x * r) * (9 * r ^ 2)) / (3 * r ^ 3 + 2 * x) ^ 2 = -1 / 5 := by
    rw [h_eq, h_den_sq]
    have h_r6 : r ^ 6 ≠ 0 := pow_ne_zero 6 hr
    calc (- (r ^ 6) * 5) / (25 * r ^ 6)
      _ = (r ^ 6 * -5) / (r ^ 6 * 25) := by ring_nf
      _ = -5 / 25 := by exact mul_div_mul_left (-5) 25 h_r6
      _ = -1 / 5 := by norm_num
      
  rw [← h_final]
  exact h_div

end Pandrosion
