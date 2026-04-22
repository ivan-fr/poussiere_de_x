/-
  Universitas Pandrosion — Legacy / Cubic fixed-point classification
  and residual conservation.

  Two self-contained results about the cube-root Pandrosion map
  `G(s) = s(s³ + 4X) / (3s³ + 2X)`:

    §L37  **Fixed-point classification.**
            G(s) = s   ⟺   s = 0   ∨   s³ = X
          (no parasitic attractors — the fixed-point set is exactly
           {0} together with the three cube roots of X). Proved both
           on ℝ and on ℂ.

    §L38  **Explicit cubic residual conservation.**
          `G(s)³ − X = (s³ − X) · Q(s, X) / (3 s³ + 2 X)³`
          with the explicit polynomial factor
          `Q(s, X) = s⁹ − 14 s⁶ X − 20 s³ X² + 8 X³`.

    §L39  **Complex conjugate symmetry.** G(z̄) = G(z) (real X).

  Ported from `ResidualConservation.lean` + `ComplexFixedPoints.lean`
  on the `main_previous` branch.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Data.Complex.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §L37. Fixed-point classification (real) -/

/-- **No parasitic attractors on ℝ.**
    For the cube-root Pandrosion map `G(s) = s(s³ + 4X)/(3s³ + 2X)`,
    the only fixed points are `s = 0` and the cube roots of `X`. -/
theorem cube_no_parasitic_attractors (X s : ℝ) (hs : 3 * s ^ 3 + 2 * X ≠ 0) :
    s * (s ^ 3 + 4 * X) / (3 * s ^ 3 + 2 * X) = s ↔ (s = 0 ∨ s ^ 3 = X) := by
  have hfrac :
      s * (s ^ 3 + 4 * X) / (3 * s ^ 3 + 2 * X) = s
      ↔ s * (s ^ 3 + 4 * X) = s * (3 * s ^ 3 + 2 * X) :=
    div_eq_iff hs
  rw [hfrac]
  rw [show s * (s ^ 3 + 4 * X) = s * (3 * s ^ 3 + 2 * X)
      ↔ s * (s ^ 3 + 4 * X) - s * (3 * s ^ 3 + 2 * X) = 0 from sub_eq_zero.symm]
  rw [show s * (s ^ 3 + 4 * X) - s * (3 * s ^ 3 + 2 * X)
        = 2 * s * (X - s ^ 3) from by ring]
  rw [show 2 * s * (X - s ^ 3) = 0 ↔ 2 * s = 0 ∨ X - s ^ 3 = 0 from mul_eq_zero]
  constructor
  · rintro (h2 | hX)
    · left; linarith
    · right; linarith
  · rintro (h0 | hcube)
    · left; linarith
    · right; linarith

/-! ## §L38. Explicit residual conservation -/

/-- **Kinematic residual conservation on ℝ.**
    The iterated residual `G(s)³ − X` factors through the input residual
    `s³ − X` with the explicit polynomial cofactor
    `s⁹ − 14 s⁶ X − 20 s³ X² + 8 X³`. -/
theorem cube_kinematic_conservation (X s : ℝ) (hs : 3 * s ^ 3 + 2 * X ≠ 0) :
    (s * (s ^ 3 + 4 * X) / (3 * s ^ 3 + 2 * X)) ^ 3 - X =
    (s ^ 3 - X) * (s ^ 9 - 14 * s ^ 6 * X - 20 * s ^ 3 * X ^ 2 + 8 * X ^ 3) /
      ((3 * s ^ 3 + 2 * X) ^ 3) := by
  set A : ℝ := s * (s ^ 3 + 4 * X) with hA
  set B : ℝ := 3 * s ^ 3 + 2 * X with hB
  have hb3 : B ^ 3 ≠ 0 := pow_ne_zero 3 hs
  calc (A / B) ^ 3 - X
      = A ^ 3 / B ^ 3 - X := by rw [div_pow]
    _ = A ^ 3 / B ^ 3 - (X * B ^ 3) / B ^ 3 := by
        rw [mul_div_cancel_right₀ X hb3]
    _ = (A ^ 3 - X * B ^ 3) / B ^ 3 := by rw [← sub_div]
    _ = (s ^ 3 - X) *
        (s ^ 9 - 14 * s ^ 6 * X - 20 * s ^ 3 * X ^ 2 + 8 * X ^ 3) / B ^ 3 := by
        have halg :
            A ^ 3 - X * B ^ 3 =
              (s ^ 3 - X) *
                (s ^ 9 - 14 * s ^ 6 * X - 20 * s ^ 3 * X ^ 2 + 8 * X ^ 3) := by
          simp only [hA, hB]; ring
        rw [halg]

/-! ## §L39. Complex fixed-point classification and conjugate symmetry -/

/-- **No parasitic attractors on ℂ.** -/
theorem complex_cube_no_parasitic_attractors
    (X z : ℂ) (hz : 3 * z ^ 3 + 2 * X ≠ 0) :
    z * (z ^ 3 + 4 * X) / (3 * z ^ 3 + 2 * X) = z ↔ (z = 0 ∨ z ^ 3 = X) := by
  have hfrac :
      z * (z ^ 3 + 4 * X) / (3 * z ^ 3 + 2 * X) = z
      ↔ z * (z ^ 3 + 4 * X) = z * (3 * z ^ 3 + 2 * X) :=
    div_eq_iff hz
  rw [hfrac]
  rw [show z * (z ^ 3 + 4 * X) = z * (3 * z ^ 3 + 2 * X)
      ↔ z * (z ^ 3 + 4 * X) - z * (3 * z ^ 3 + 2 * X) = 0 from sub_eq_zero.symm]
  rw [show z * (z ^ 3 + 4 * X) - z * (3 * z ^ 3 + 2 * X)
        = 2 * z * (X - z ^ 3) from by ring]
  rw [show 2 * z * (X - z ^ 3) = 0 ↔ 2 * z = 0 ∨ X - z ^ 3 = 0 from mul_eq_zero]
  constructor
  · rintro (h2 | hX)
    · left
      have : (2 : ℂ) ≠ 0 := by norm_num
      exact (mul_eq_zero.mp h2).resolve_left this
    · right; linear_combination -hX
  · rintro (h0 | hcube)
    · left; rw [h0]; ring
    · right; linear_combination -hcube

/-- **Complex conjugate symmetry.** For real `X`, the cube-root Pandrosion
    map commutes with complex conjugation: `G(z̄) = G(z)` where overline
    denotes `star z`. -/
theorem cube_complex_conjugate_symmetry (X : ℝ) (z : ℂ) :
    (fun (w : ℂ) => w * (w ^ 3 + 4 * (X : ℂ)) / (3 * w ^ 3 + 2 * (X : ℂ)))
      (star z) =
    star
      ((fun (w : ℂ) => w * (w ^ 3 + 4 * (X : ℂ)) / (3 * w ^ 3 + 2 * (X : ℂ))) z) := by
  have hX : star (X : ℂ) = (X : ℂ) := by
    apply Complex.ext <;> simp
  calc  star z * (star z ^ 3 + 4 * (X : ℂ)) / (3 * star z ^ 3 + 2 * (X : ℂ))
      = star z * (star z ^ 3 + 4 * star (X : ℂ)) /
          (3 * star z ^ 3 + 2 * star (X : ℂ)) := by rw [hX]
    _ = star (z * (z ^ 3 + 4 * (X : ℂ)) / (3 * z ^ 3 + 2 * (X : ℂ))) := by simp

end Pandrosion
