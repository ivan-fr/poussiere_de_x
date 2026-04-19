/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP XXIV: THE MATRILINEAR MATRIX IDENTITY
  
  This module extends the Snail Convergence identity into hyper-dimensional
  spaces. It proves that the algebraic error factorization applies to
  non-commutative Rings (Matrices, Operators), provided the target R
  and the approximation S commute.
-/
import Mathlib.Algebra.Ring.Commute
import Mathlib.Tactic

namespace Pandrosion

/-! ## §129. Non-Commutative Algebraic Invariance -/

/-- **The Matrilinear Pandrosion Oscillation.**
    Proves that the shear error factorization works on non-commutative 
    Matrices / Operators as long as the approximation commutes with the exact root. -/
theorem matrix_pandrosion_oscillation {α : Type*} [Ring α] 
    (X R S : α) (h_root : R^3 = X) (h_comm : Commute S R) :
    let U := S * (S^3 + 4 • X)
    let V := 3 • S^3 + 2 • X
    U - R * V = (S - R) * (S^3 - 2 • (R * S^2) - 2 • (R^2 * S) + 2 • R^3) := by
  intros U V
  have h_root_sub : X = R^3 := h_root.symm
  dsimp [U, V]
  rw [h_root_sub]
  
  -- Extract commutativity equations
  have h_SR : S * R = R * S := h_comm.eq
  have h_S2R : S^2 * R = R * S^2 := (h_comm.pow_left 2).eq
  have h_SR2 : S * R^2 = R^2 * S := (h_comm.pow_right 2).eq
  have h_SR3 : S * R^3 = R^3 * S := (h_comm.pow_right 3).eq
  
  -- Transform LHS
  have hLeft : S * (S ^ 3 + 4 • R ^ 3) - R * (3 • S ^ 3 + 2 • R ^ 3) = 
               S^4 + 4 • (R^3 * S) - 3 • (R * S^3) - 2 • R^4 := by
    calc S * (S ^ 3 + 4 • R ^ 3) - R * (3 • S ^ 3 + 2 • R ^ 3)
      _ = S * S^3 + S * (4 • R^3) - (R * (3 • S^3) + R * (2 • R^3)) := by simp only [mul_add]
      _ = S^4 + 4 • (S * R^3) - 3 • (R * S^3) - 2 • (R * R^3) := by simp only [mul_smul_comm, sub_add_eq_sub_sub, ←pow_succ']
      _ = S^4 + 4 • (R^3 * S) - 3 • (R * S^3) - 2 • R^4 := by rw [h_SR3, ←pow_succ']
    
  -- Transform RHS
  have hRight : (S - R) * (S^3 - 2 • (R * S^2) - 2 • (R^2 * S) + 2 • R^3) =
                S^4 + 4 • (R^3 * S) - 3 • (R * S^3) - 2 • R^4 := by
    have h1 : S * S^3 = S^4 := (pow_succ' S 3).symm
    have h2 : R * R^3 = R^4 := (pow_succ' R 3).symm
    have h3 : S * (R * S^2) = R * S^3 := by
      calc S * (R * S^2) = (S * R) * S^2 := by rw [←mul_assoc]
        _ = (R * S) * S^2 := by rw [h_SR]
        _ = R * (S * S^2) := by rw [mul_assoc]
        _ = R * S^3 := by rw [←pow_succ' S 2]
    have h4_2 : S * (R^2 * S) = R * (R * S^2) := by
      calc S * (R^2 * S) = (S * R^2) * S := by rw [←mul_assoc]
        _ = (R^2 * S) * S := by rw [h_SR2]
        _ = R^2 * (S * S) := by rw [mul_assoc]
        _ = (R * R) * (S * S) := by rw [sq R]
        _ = (R * R) * S^2 := by rw [←sq S]
        _ = R * (R * S^2) := by rw [mul_assoc]
    have h6 : R * (R^2 * S) = R^3 * S := by 
      calc R * (R^2 * S) = (R * R^2) * S := by rw [←mul_assoc]
        _ = R^3 * S := by rw [←pow_succ' R 2]
    calc (S - R) * (S^3 - 2 • (R * S^2) - 2 • (R^2 * S) + 2 • R^3)
      _ = S * (S^3 - 2 • (R * S^2) - 2 • (R^2 * S) + 2 • R^3) - 
          R * (S^3 - 2 • (R * S^2) - 2 • (R^2 * S) + 2 • R^3) := by rw [sub_mul]
      _ = S * S^3 - S * (2 • (R * S^2)) - S * (2 • (R^2 * S)) + S * (2 • R^3) - 
          (R * S^3 - R * (2 • (R * S^2)) - R * (2 • (R^2 * S)) + R * (2 • R^3)) := by
          simp only [mul_sub, mul_add]
      _ = S^4 - 2 • (S * (R * S^2)) - 2 • (S * (R^2 * S)) + 2 • (S * R^3) - 
          (R * S^3 - 2 • (R * (R * S^2)) - 2 • (R * (R^2 * S)) + 2 • (R * R^3)) := by
          simp only [mul_smul_comm, ←pow_succ']
      _ = S^4 - 2 • (R * S^3) - 2 • (R * (R * S^2)) + 2 • (R^3 * S) - 
          (R * S^3 - 2 • (R * (R * S^2)) - 2 • (R^3 * S) + 2 • R^4) := by
          rw [h3, h4_2, h_SR3, h6, ←pow_succ']
      _ = S^4 + 4 • (R^3 * S) - 3 • (R * S^3) - 2 • R^4 := by abel

  rw [hLeft, hRight]

end Pandrosion
