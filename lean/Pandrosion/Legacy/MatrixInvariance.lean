/-
  Universitas Pandrosion — Legacy / Matrix invariance.

  Non-commutative algebraic invariance of the cube-root Pandrosion
  numerator `U = S(S³ + 4X)` and denominator `V = 3 S³ + 2 X`,
  lifted from ℂ to `Matrix n n ℂ`:

    §L44  **Commutation propagation.** If `S` commutes with `X`,
          then `S` commutes with `U` and with `V`, and `U` commutes
          with `V`.

    §L45  **Hermitian preservation.** If additionally `S = star S`
          and `X = star X`, then `U, V` are Hermitian, and the
          Pandrosion product `U · V` is Hermitian.

  Ported from `CommutativityPropagation.lean` +
  `HermitianPreservation.lean` on the `main_previous` branch
  (with the stale `Pandrosion.MatrixSpectral` dependency elided
  — nothing from that module is actually used).
-/

import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Complex.Basic
import Mathlib.Algebra.Star.Basic
import Mathlib.Algebra.Star.Pi
import Mathlib.Tactic

namespace Pandrosion

open Matrix

/-! ## §L44. Commutation propagation -/

/-- `S` commutes with the numerator `U = S(S³ + 4X)`. -/
theorem matrix_phase_lock_U {n : Type*} [Fintype n] [DecidableEq n]
    (X S : Matrix n n ℂ) (h_comm : Commute S X) :
    Commute S (S * (S ^ 3 + 4 • X)) := by
  have h_S3 : Commute S (S ^ 3) := (Commute.refl S).pow_right 3
  have h_4X : Commute S (4 • X) := h_comm.smul_right 4
  have h_part : Commute S (S ^ 3 + 4 • X) := h_S3.add_right h_4X
  exact Commute.mul_right (Commute.refl S) h_part

/-- `S` commutes with the denominator `V = 3 S³ + 2X`. -/
theorem matrix_phase_lock_V {n : Type*} [Fintype n] [DecidableEq n]
    (X S : Matrix n n ℂ) (h_comm : Commute S X) :
    Commute S (3 • S ^ 3 + 2 • X) := by
  have h_S3 : Commute S (3 • S ^ 3) :=
    Commute.smul_right ((Commute.refl S).pow_right 3) 3
  have h_2X : Commute S (2 • X) := Commute.smul_right h_comm 2
  exact Commute.add_right h_S3 h_2X

/-- `U` commutes with `V`: the full Pandrosion numerator and
    denominator commute as matrices whenever `S` commutes with `X`. -/
theorem matrix_phase_lock_UV {n : Type*} [Fintype n] [DecidableEq n]
    (X S : Matrix n n ℂ) (h_comm : Commute S X) :
    Commute (S * (S ^ 3 + 4 • X)) (3 • S ^ 3 + 2 • X) := by
  -- S commutes with V
  have h_SV : Commute S (3 • S ^ 3 + 2 • X) := matrix_phase_lock_V X S h_comm
  -- X commutes with V
  have h_XV : Commute X (3 • S ^ 3 + 2 • X) := by
    have h_S3 : Commute X (3 • S ^ 3) :=
      Commute.smul_right (h_comm.symm.pow_right 3) 3
    have h_X : Commute X (2 • X) := Commute.smul_right (Commute.refl X) 2
    exact Commute.add_right h_S3 h_X
  -- S³ commutes with V (via S), 4X commutes with V (via X)
  have h_S3_V : Commute (S ^ 3) (3 • S ^ 3 + 2 • X) := h_SV.pow_left 3
  have h_4X_V : Commute (4 • X) (3 • S ^ 3 + 2 • X) := h_XV.smul_left 4
  -- S³ + 4X commutes with V
  have h_diff_V : Commute (S ^ 3 + 4 • X) (3 • S ^ 3 + 2 • X) :=
    Commute.add_left h_S3_V h_4X_V
  -- Phase-lock: U = S · (S³ + 4X) commutes with V
  exact Commute.mul_left h_SV h_diff_V

/-! ## §L45. Hermitian preservation -/

/-- If `S` and `X` are Hermitian and commute, then the numerator
    `U = S(S³ + 4X)` is Hermitian. -/
theorem matrix_hermitian_U {n : Type*} [Fintype n] [DecidableEq n]
    (X S : Matrix n n ℂ) (h_comm : Commute S X)
    (hX_herm : star X = X) (hS_herm : star S = S) :
    star (S * (S ^ 3 + 4 • X)) = S * (S ^ 3 + 4 • X) := by
  have h_star : star (S * (S ^ 3 + 4 • X))
      = star (S ^ 3 + 4 • X) * star S := StarMul.star_mul S (S ^ 3 + 4 • X)
  have h_star_S3 : star (S ^ 3) = S ^ 3 := by
    calc star (S ^ 3) = (star S) ^ 3 := star_pow S 3
      _ = S ^ 3 := by rw [hS_herm]
  have h_star_4X : star (4 • X) = 4 • X := by
    calc star (4 • X) = 4 • star X := star_nsmul X 4
      _ = 4 • X := by rw [hX_herm]
  have h_star_part : star (S ^ 3 + 4 • X) = S ^ 3 + 4 • X := by
    calc star (S ^ 3 + 4 • X) = star (S ^ 3) + star (4 • X) :=
          star_add (S ^ 3) (4 • X)
      _ = S ^ 3 + 4 • X := by rw [h_star_S3, h_star_4X]
  rw [h_star, h_star_part, hS_herm]
  -- Now prove (S³ + 4X) · S = S · (S³ + 4X), i.e., S commutes with S³ + 4X.
  have h_S3 : Commute S (S ^ 3) := (Commute.refl S).pow_right 3
  have h_4X : Commute S (4 • X) := h_comm.smul_right 4
  exact (h_S3.add_right h_4X).symm.eq

/-- If `S` and `X` are Hermitian, the denominator `V = 3 S³ + 2 X`
    is Hermitian (no commutation hypothesis needed). -/
theorem matrix_hermitian_V {n : Type*} [Fintype n] [DecidableEq n]
    (X S : Matrix n n ℂ)
    (hX_herm : star X = X) (hS_herm : star S = S) :
    star (3 • S ^ 3 + 2 • X) = 3 • S ^ 3 + 2 • X := by
  have h_star_S3 : star (S ^ 3) = S ^ 3 := by
    calc star (S ^ 3) = (star S) ^ 3 := star_pow S 3
      _ = S ^ 3 := by rw [hS_herm]
  have h_star_3S3 : star (3 • S ^ 3) = 3 • S ^ 3 := by
    calc star (3 • S ^ 3) = 3 • star (S ^ 3) := star_nsmul (S ^ 3) 3
      _ = 3 • S ^ 3 := by rw [h_star_S3]
  have h_star_2X : star (2 • X) = 2 • X := by
    calc star (2 • X) = 2 • star X := star_nsmul X 2
      _ = 2 • X := by rw [hX_herm]
  calc star (3 • S ^ 3 + 2 • X)
      = star (3 • S ^ 3) + star (2 • X) := star_add (3 • S ^ 3) (2 • X)
    _ = 3 • S ^ 3 + 2 • X := by rw [h_star_3S3, h_star_2X]

/-- **Pandrosion product is Hermitian.** If `S` and `X` are Hermitian
    and commute, the matrix product `U · V = S(S³+4X) · (3S³+2X)` is
    Hermitian. -/
theorem matrix_hermitian_UV {n : Type*} [Fintype n] [DecidableEq n]
    (X S : Matrix n n ℂ) (h_comm : Commute S X)
    (hX_herm : star X = X) (hS_herm : star S = S) :
    star ((S * (S ^ 3 + 4 • X)) * (3 • S ^ 3 + 2 • X))
      = (S * (S ^ 3 + 4 • X)) * (3 • S ^ 3 + 2 • X) := by
  have h_UV_comm : Commute (S * (S ^ 3 + 4 • X)) (3 • S ^ 3 + 2 • X) :=
    matrix_phase_lock_UV X S h_comm
  have h_star_U : star (S * (S ^ 3 + 4 • X)) = S * (S ^ 3 + 4 • X) :=
    matrix_hermitian_U X S h_comm hX_herm hS_herm
  have h_star_V : star (3 • S ^ 3 + 2 • X) = 3 • S ^ 3 + 2 • X :=
    matrix_hermitian_V X S hX_herm hS_herm
  calc star ((S * (S ^ 3 + 4 • X)) * (3 • S ^ 3 + 2 • X))
      = star (3 • S ^ 3 + 2 • X) * star (S * (S ^ 3 + 4 • X)) :=
        StarMul.star_mul _ _
    _ = (3 • S ^ 3 + 2 • X) * (S * (S ^ 3 + 4 • X)) := by
        rw [h_star_V, h_star_U]
    _ = (S * (S ^ 3 + 4 • X)) * (3 • S ^ 3 + 2 • X) := h_UV_comm.symm.eq

end Pandrosion
