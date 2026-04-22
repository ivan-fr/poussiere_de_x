/-
  Universitas Pandrosion — Legacy / Block descent.

  Block-form lemmas complementing `Legacy.Descent`:

    §L40  `cos(π/(2d)) ∈ (0, 1)` and `cos(π/(2d)) ≤ cos(π/4)`
          bounds for `d ≥ 2`, used to express the block-descent
          ratio `cos^d(π/(2d))`.
    §L41  Block-descent ratio `cos^d(π/(2d)) ∈ (0, 1)`; equivalence
          with the symmetric-epoch descent `d · log(cos(π/(2d)))`.
    §L42  Product-descent bound `((γ+ε)/(1−ε))^d < 1` under
          `0 < γ < 1`, `0 ≤ ε`, `γ + ε < 1 − ε`.
    §L43  Energy decrease `(ρ/R)^k` is strictly decreasing in `R`.

  Ported from `Descent.lean` + `Descent2.lean` on the
  `main_previous` branch.
-/

import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

namespace Pandrosion

open Real

/-! ## §L40. `cos(π/(2d))` bounds -/

/-- For `d ≥ 2`, the primary descent angle `π/(2d)` is in `(0, π/2)`. -/
theorem primary_angle_in_range (d : ℕ) (hd : d ≥ 2) :
    0 < π / (2 * (d : ℝ)) ∧ π / (2 * (d : ℝ)) < π / 2 := by
  refine ⟨?_, ?_⟩
  · apply div_pos pi_pos; positivity
  · apply div_lt_div_of_pos_left pi_pos (by positivity)
      (by linarith [show (2 : ℝ) ≤ (d : ℝ) from by exact_mod_cast hd])

/-- `cos(π/(2d)) < 1` for `d ≥ 2`. -/
theorem cos_primary_angle_lt_one (d : ℕ) (hd : d ≥ 2) :
    cos (π / (2 * (d : ℝ))) < 1 := by
  have ⟨h_pos, h_lt⟩ := primary_angle_in_range d hd
  have hsin : 0 < sin (π / (2 * (d : ℝ))) :=
    sin_pos_of_pos_of_lt_pi h_pos (by linarith [pi_pos])
  have hsc := sin_sq_add_cos_sq (π / (2 * (d : ℝ)))
  nlinarith [sq_nonneg (sin (π / (2 * (d : ℝ)))),
             sq_nonneg (cos (π / (2 * (d : ℝ))) - 1)]

/-- `0 < cos(π/(2d))` for `d ≥ 2`. -/
theorem cos_primary_angle_pos (d : ℕ) (hd : d ≥ 2) :
    0 < cos (π / (2 * (d : ℝ))) := by
  have ⟨h_pos, h_lt⟩ := primary_angle_in_range d hd
  exact cos_pos_of_mem_Ioo ⟨by linarith, h_lt⟩

/-- `log(cos(π/(2d))) < 0` for `d ≥ 2`. -/
theorem log_cos_primary_neg (d : ℕ) (hd : d ≥ 2) :
    log (cos (π / (2 * (d : ℝ)))) < 0 :=
  log_neg (cos_primary_angle_pos d hd) (cos_primary_angle_lt_one d hd)

/-- Epoch descent `d · log(cos(π/(2d))) < 0` for `d ≥ 2`. -/
theorem epoch_descent_neg (d : ℕ) (hd : d ≥ 2) :
    (d : ℝ) * log (cos (π / (2 * (d : ℝ)))) < 0 := by
  apply mul_neg_of_pos_of_neg
  · exact_mod_cast (show 0 < d by omega)
  · exact log_cos_primary_neg d hd

/-- `π/(2d) ≤ π/4` for `d ≥ 2`. -/
theorem primary_angle_le_pi_div_four (d : ℕ) (hd : d ≥ 2) :
    π / (2 * (d : ℝ)) ≤ π / 4 := by
  have hd_cast : (2 : ℝ) ≤ (d : ℝ) := by exact_mod_cast hd
  have h4 : (4 : ℝ) ≤ 2 * (d : ℝ) := by linarith
  exact div_le_div_of_nonneg_left (le_of_lt pi_pos)
    (by norm_num : (0 : ℝ) < 4) h4

/-- `cos(π/4) ≤ cos(π/(2d))` for `d ≥ 2`. -/
theorem contraction_bounded_below (d : ℕ) (hd : d ≥ 2) :
    cos (π / 4) ≤ cos (π / (2 * (d : ℝ))) := by
  have h1 := primary_angle_in_range d hd
  apply cos_le_cos_of_nonneg_of_le_pi
  · exact le_of_lt h1.1
  · linarith [h1.2, pi_pos]
  · exact primary_angle_le_pi_div_four d hd

/-! ## §L41. Block-descent ratio `cos^d(π/(2d))` -/

/-- `cos^d(π/(2d)) > 0` for `d ≥ 2`. -/
theorem block_descent_ratio_pos (d : ℕ) (hd : d ≥ 2) :
    0 < cos (π / (2 * (d : ℝ))) ^ d :=
  pow_pos (cos_primary_angle_pos d hd) d

/-- `cos^d(π/(2d)) < 1` for `d ≥ 2`. -/
theorem block_descent_ratio_lt_one (d : ℕ) (hd : d ≥ 2) :
    cos (π / (2 * (d : ℝ))) ^ d < 1 :=
  pow_lt_one (le_of_lt (cos_primary_angle_pos d hd))
    (cos_primary_angle_lt_one d hd) (by omega)

/-- `d · log(cos(π/(2d))) = log(cos^d(π/(2d)))`. -/
theorem epoch_descent_equals_log_block (d : ℕ) (_hd : d ≥ 2) :
    (d : ℝ) * log (cos (π / (2 * (d : ℝ)))) =
      log (cos (π / (2 * (d : ℝ))) ^ d) := by
  rw [Real.log_pow]

/-- `log(cos^d(π/(2d))) < 0` for `d ≥ 2`. -/
theorem log_block_descent_neg (d : ℕ) (hd : d ≥ 2) :
    log (cos (π / (2 * (d : ℝ))) ^ d) < 0 :=
  log_neg (block_descent_ratio_pos d hd) (block_descent_ratio_lt_one d hd)

/-! ## §L42. Product-descent bound -/

/-- **Product-descent bound.** Under `0 < γ < 1`, `0 ≤ ε`,
    `γ + ε < 1 − ε`, and `d ≥ 1`, the product rate
    `((γ + ε)/(1 − ε))^d < 1`. -/
theorem product_descent_bound
    (gamma eps : ℝ)
    (hg : 0 < gamma) (_hg1 : gamma < 1)
    (he : 0 ≤ eps) (hge : gamma + eps < 1 - eps)
    (d : ℕ) (hd : d ≥ 1) :
    ((gamma + eps) / (1 - eps)) ^ d < 1 := by
  have h1eps_pos : 0 < 1 - eps := by linarith
  apply pow_lt_one
  · exact div_nonneg (by linarith) (by linarith)
  · rw [div_lt_one h1eps_pos]; exact hge
  · omega

/-- Convenience: for `d ≥ 3` and `R = 3ρ`, the correction mass
    `(1/3)^d < 1/2`. This is the typical instantiation of
    `product_descent_bound`. -/
theorem descent_eps_bound (d : ℕ) (hd : d ≥ 3) :
    ((1 : ℝ) / 3) ^ d < 1 / 2 := by
  calc ((1 : ℝ) / 3) ^ d
      ≤ (1 / 3) ^ 3 := by
        apply pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega)
    _ = 1 / 27 := by norm_num
    _ < 1 / 2 := by norm_num

/-! ## §L43. Energy normalisation -/

/-- `(ρ/R)^k < 1` for every mode `k ≥ 1`, whenever `0 < ρ < R`. -/
theorem energy_normalizes (rho R : ℝ)
    (hrho : 0 < rho) (hR : rho < R)
    (k : ℕ) (hk : k ≥ 1) :
    (rho / R) ^ k < 1 := by
  exact pow_lt_one
    (le_of_lt (div_pos hrho (by linarith)))
    (by rw [div_lt_one (by linarith)]; exact hR)
    (by omega)

/-- Quadratic decay `(ρ/R)² < 1`. -/
theorem energy_excess_quadratic_decay (rho R : ℝ)
    (hrho : 0 < rho) (hR : rho < R) :
    (rho / R) ^ 2 < 1 :=
  energy_normalizes rho R hrho hR 2 (by omega)

/-- Energy is strictly decreasing in the outer radius `R`:
    for `rho < R₁ < R₂` and `k ≥ 1`, `(ρ/R₂)^k < (ρ/R₁)^k`. -/
theorem energy_decreasing (rho R1 R2 : ℝ)
    (hrho : 0 < rho) (hR1 : rho < R1) (hR2 : R1 < R2)
    (k : ℕ) (hk : k ≥ 1) :
    (rho / R2) ^ k < (rho / R1) ^ k := by
  apply pow_lt_pow_left
  · apply div_lt_div_of_pos_left hrho (by linarith) (by linarith)
  · exact le_of_lt (div_pos hrho (by linarith))
  · omega

end Pandrosion
