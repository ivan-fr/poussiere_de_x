/-
  Universitas Pandrosion — §55. **Auxiliary bound `2^{−1/3} ≥ 3/4`.**

  Isolated technical lemma needed by §56 to instantiate the
  p=3, x=2 half-plane contraction theorem at the concrete
  cube root `α = 2^{−1/3}`.

  **Proof strategy:** cube both sides.
    `2^{−1/3} ≥ 3/4`
    ⟺ `(2^{−1/3})³ ≥ (3/4)³`   (both positive, cube is monotone)
    ⟺ `1/2 ≥ 27/64`
    ⟺ `32/64 ≥ 27/64` ✓

  The technical ingredient is `Real.rpow_natCast` and
  `Real.rpow_mul` to compute `(2^{−1/3})³ = 2^{−1} = 1/2`.

  Contents.

    §55.1  `two_rpow_neg_third_pow_three` — algebraic identity
           `(2^{−1/3})³ = 1/2`.

    §55.2  `two_rpow_neg_third_ge_three_fourths` — the bound
           `2^{−1/3} ≥ 3/4`.
-/

import Mathlib.Analysis.SpecialFunctions.Pow.Real

namespace Pandrosion

/-! ============================================================
  §55.1  Algebraic identity `(2^{−1/3})³ = 1/2`
============================================================ -/

section CubeIdentityC

/-- **`(2^{−1/3})³ = 1/2`.**
    Direct from `Real.rpow_nat_cast` and `Real.rpow_mul`. -/
theorem two_rpow_neg_third_pow_three :
    ((2 : ℝ) ^ (-(1 : ℝ) / 3)) ^ 3 = 1 / 2 := by
  have h2_pos : (0 : ℝ) < 2 := by norm_num
  have h_pow_to_rpow : ((2 : ℝ) ^ (-(1 : ℝ) / 3)) ^ 3
                     = (2 : ℝ) ^ ((-(1 : ℝ) / 3) * 3) := by
    rw [← Real.rpow_nat_cast ((2 : ℝ) ^ (-(1 : ℝ) / 3)) 3,
        ← Real.rpow_mul (le_of_lt h2_pos)]
    norm_num
  rw [h_pow_to_rpow]
  have h_exp : (-(1 : ℝ) / 3) * 3 = -1 := by ring
  rw [h_exp, Real.rpow_neg_one]
  norm_num

end CubeIdentityC

/-! ============================================================
  §55.2  Bound `2^{−1/3} ≥ 3/4`
============================================================ -/

section BoundC

/-- **★ `2^{−1/3} ≥ 3/4`.**

    Proved by cubing both sides:
      `(2^{−1/3})³ = 1/2 = 32/64`,
      `(3/4)³ = 27/64`,
    and `32/64 ≥ 27/64`. Since cube is monotone on ℝ₊, the
    inequality transfers. -/
theorem two_rpow_neg_third_ge_three_fourths :
    ((2 : ℝ) ^ (-(1 : ℝ) / 3)) ≥ 3 / 4 := by
  have h2_pos : (0 : ℝ) < 2 := by norm_num
  have h_base_pos : (0 : ℝ) < (2 : ℝ) ^ (-(1 : ℝ) / 3) :=
    Real.rpow_pos_of_pos h2_pos _
  have h_three_four_nn : (0 : ℝ) ≤ 3 / 4 := by norm_num
  have h_cube_eq : ((2 : ℝ) ^ (-(1 : ℝ) / 3)) ^ 3 = 1 / 2 :=
    two_rpow_neg_third_pow_three
  have h_cube_bound : ((3 : ℝ) / 4) ^ 3 ≤ ((2 : ℝ) ^ (-(1 : ℝ) / 3)) ^ 3 := by
    rw [h_cube_eq]
    norm_num
  -- Cube is monotone on ℝ₊. `pow_le_pow_left` with appropriate signs.
  exact (pow_le_pow_iff_left h_three_four_nn (le_of_lt h_base_pos)
    (by norm_num : 3 ≠ 0)).mp h_cube_bound

end BoundC

end Pandrosion
