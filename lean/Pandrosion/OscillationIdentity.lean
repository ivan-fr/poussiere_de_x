/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP XX: THE OSCILLATION IDENTITY
  
  This module formalizes the unique signature of the Generalized
  Pandrosion map: its alternating linear convergence pattern (the "snail"),
  caused by a negative ratio identity. We prove this exact polynomial 
  decomposition strictly confirming the existence of the oscillation.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

open Real

namespace Pandrosion

/-! ## §123. The Oscillation Signature for p=3 -/

/-- **The Pandrosion Oscillation Identity.**
    This isolates the exact polynomial factor responsible for the
    alternating negative ratio λ_3 = -1/5 during the approach phase. -/
theorem pandrosion_oscillation_identity (x r s : ℝ) 
    (h_root : r^3 = x) (h_denom : 3 * s^3 + 2 * x ≠ 0) :
    pandrosion_generalized_map 3 x s - r = 
    (s - r) * (s^3 - 2*r*s^2 - 2*r^2*s + 2*r^3) / (3 * s^3 + 2 * x) := by
  have h_mul : (pandrosion_generalized_map 3 x s - r) * (3 * s^3 + 2 * x) = 
               (s - r) * (s^3 - 2*r*s^2 - 2*r^2*s + 2*r^3) := by
    unfold pandrosion_generalized_map
    dsimp
    have h_den_eq : (3 : ℝ) * s ^ 3 + ((3 : ℝ) - 1) * ((3 : ℝ) - 2) * x = 3 * s ^ 3 + 2 * x := by ring
    have h_den_sub : (3 : ℝ) * s ^ 3 + ((3 : ℝ) - 1) * ((3 : ℝ) - 2) * x ≠ 0 := by
      rw [h_den_eq]
      exact h_denom
    rw [if_neg h_den_sub]
    rw [h_den_eq]
    
    have h_div_cancel := div_mul_cancel₀ (s * (s ^ 3 + ((3:ℝ) - 1) ^ 2 * x)) h_denom
    calc (s * (s ^ 3 + ((3:ℝ) - 1) ^ 2 * x) / (3 * s ^ 3 + 2 * x) - r) * (3 * s ^ 3 + 2 * x)
      _ = s * (s ^ 3 + ((3:ℝ) - 1) ^ 2 * x) / (3 * s ^ 3 + 2 * x) * (3 * s ^ 3 + 2 * x) - r * (3 * s ^ 3 + 2 * x) := by ring
      _ = s * (s ^ 3 + ((3:ℝ) - 1) ^ 2 * x) - r * (3 * s ^ 3 + 2 * x) := by rw [h_div_cancel]
      _ = s * (s ^ 3 + 4 * x) - r * (3 * s ^ 3 + 2 * x) := by ring
      _ = (s - r) * (s ^ 3 - 2 * r * s ^ 2 - 2 * r ^ 2 * s + 2 * r ^ 3) := by 
        rw [← h_root]
        ring

  exact eq_div_of_mul_eq h_denom h_mul

end Pandrosion
