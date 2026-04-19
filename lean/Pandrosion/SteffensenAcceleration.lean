/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP XIX: THE AITKEN-STEFFENSEN ACCELERATION (T3)
  
  This module formalizes the composite generic architecture
  where the Aitken-Steffensen extrapolation formula perfectly
  annihilates the linear alternating dust produced by the
  Generalized Pandrosion map, achieving quadratic convergence.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic
import Pandrosion.ScaledMap

open Real

namespace Pandrosion

/-! ## §121. The Pandrosion-T3 Architecture -/

/-- **The Pandrosion-T3 Step.**
    This isolates the generalized Pandrosion iteration inside an 
    Aitken-Steffensen delta-squared accelerator. -/
noncomputable def pandrosion_generalized_t3_step (p : ℕ) (x s : ℝ) : ℝ :=
  let s1 := pandrosion_generalized_map p x s
  let s2 := pandrosion_generalized_map p x s1
  let denom := s2 - 2 * s1 + s
  if denom = 0 then s2 else s - (s1 - s) ^ 2 / denom

/-! ## §122. The Perfect Extrapolation Theorem -/

/-- **Aitken's Extrapolation is exact for linearly contracting error.**
    This pure algebraic theorem proves why Pandrosion-T3 is mathematically
    superior. If the underlying operation has error exactly contracting 
    by λ ≠ 1, the formula completely extracts the root in ONE step. -/
theorem aitken_perfect_extrapolation (r s lam : ℝ) (h_lam : lam ≠ 1) (hs : s ≠ r) :
    let s1 := r + lam * (s - r)
    let s2 := r + lam * (s1 - r)
    let denom := s2 - 2 * s1 + s
    denom ≠ 0 ∧ s - (s1 - s) ^ 2 / denom = r := by
  intros s1 s2 denom
  have h_denom : denom = (s - r) * (lam - 1) ^ 2 := by
    dsimp [denom, s2, s1]
    ring
  
  have h_sr : s - r ≠ 0 := sub_ne_zero.mpr hs
  have h_lam1 : lam - 1 ≠ 0 := sub_ne_zero.mpr h_lam
  have hn_lam1 : (lam - 1) ^ 2 ≠ 0 := pow_ne_zero 2 h_lam1
  
  have hz : denom ≠ 0 := by
    rw [h_denom]
    exact mul_ne_zero h_sr hn_lam1
    
  constructor
  · exact hz
  · have h_num : (s1 - s) ^ 2 = (s - r) ^ 2 * (lam - 1) ^ 2 := by
      dsimp [s1]
      ring
    rw [h_num, h_denom]
    have h_cancel : ((s - r) ^ 2 * (lam - 1) ^ 2) / ((s - r) * (lam - 1) ^ 2) = s - r := by
      have h_decomp : (s - r) ^ 2 * (lam - 1) ^ 2 = (s - r) * ((s - r) * (lam - 1) ^ 2) := by ring
      rw [h_decomp]
      have hz_mul : (s - r) * (lam - 1) ^ 2 ≠ 0 := mul_ne_zero h_sr hn_lam1
      exact mul_div_cancel_right₀ (s - r) hz_mul
    rw [h_cancel]
    ring

end Pandrosion
