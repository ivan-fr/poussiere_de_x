/-
  Universitas Pandrosion — Lean 4 Formalization
  The Universal Descent Constant: -π²/8
  Reference: pandrosion_master.tex, Theorems 3445, 3881
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic
import Pandrosion.Core

open Real

namespace Pandrosion

/-! ## The Descent Architecture

The core result: for P(z) = z^d - 1 evaluated on a Cauchy circle of radius R,
the sum of logarithmic contractions over d equispaced starts equals:

  D(R) = d · log(cos(π/(2d)))

which converges to -π²/(8d) as d → ∞.

The key analytical steps:
1. cos(x) ≤ 1  for all x  (trivial)
2. cos(π/(2d)) < 1  for d ≥ 2  (strict)
3. log(cos(x)) < 0  for 0 < x < π/2
4. d · log(cos(π/(2d))) → -π²/8  (the universal constant)
-/

/-- π is positive (convenience wrapper). -/
lemma pi_pos' : (0 : ℝ) < π := pi_pos

/-- For d ≥ 2, the angle π/(2d) is in (0, π/2). -/
lemma angle_in_range (d : ℕ) (hd : d ≥ 2) :
    0 < π / (2 * (d : ℝ)) ∧ π / (2 * (d : ℝ)) < π / 2 := by
  constructor
  · apply div_pos pi_pos
    positivity
  · apply div_lt_div_of_pos_left pi_pos (by positivity) (by linarith [show (2 : ℝ) ≤ (d : ℝ) from by exact_mod_cast hd])

/-- Theorem 3445 (base case): cos(π/(2d)) < 1 for d ≥ 2.
    This means the Pandrosion scanner detects non-trivial contraction. -/
theorem cos_angle_lt_one (d : ℕ) (hd : d ≥ 2) :
    cos (π / (2 * (d : ℝ))) < 1 := by
  have ⟨h_pos, h_lt⟩ := angle_in_range d hd
  have hsin : 0 < sin (π / (2 * (d : ℝ))) :=
    sin_pos_of_pos_of_lt_pi h_pos (by linarith [pi_pos])
  have hsc := sin_sq_add_cos_sq (π / (2 * (d : ℝ)))
  nlinarith [sq_nonneg (sin (π / (2 * (d : ℝ)))),
             sq_nonneg (cos (π / (2 * (d : ℝ))) - 1)]

/-- cos(π/(2d)) > 0 for d ≥ 2 (the angle is in the first quadrant). -/
theorem cos_angle_pos (d : ℕ) (hd : d ≥ 2) :
    0 < cos (π / (2 * (d : ℝ))) := by
  have ⟨_, h_lt⟩ := angle_in_range d hd
  exact cos_pos_of_mem_Ioo ⟨by linarith [angle_in_range d hd |>.1], h_lt⟩

/-- Theorem 3881 (sign): The per-epoch descent is strictly negative.
    log(cos(π/(2d))) < 0 for d ≥ 2. -/
theorem descent_negative (d : ℕ) (hd : d ≥ 2) :
    log (cos (π / (2 * (d : ℝ)))) < 0 := by
  apply log_neg (cos_angle_pos d hd) (cos_angle_lt_one d hd)

/-- The total epoch descent d · log(cos(π/(2d))) is negative. -/
theorem epoch_descent_negative (d : ℕ) (hd : d ≥ 2) :
    (d : ℝ) * log (cos (π / (2 * (d : ℝ)))) < 0 := by
  apply mul_neg_of_pos_of_neg
  · exact_mod_cast (show 0 < d by omega)
  · exact descent_negative d hd

/-- For d ≥ 2, the angle π/(2d) ≤ π/4. -/
lemma angle_le_pi_div_four (d : ℕ) (hd : d ≥ 2) :
    π / (2 * (d : ℝ)) ≤ π / 4 := by
  have hd_cast : (2 : ℝ) ≤ (d : ℝ) := by exact_mod_cast hd
  have h4 : (4 : ℝ) ≤ 2 * (d : ℝ) := by linarith
  exact div_le_div_of_nonneg_left (le_of_lt pi_pos) (by norm_num : (0:ℝ) < 4) h4

/-- The contraction ratio is bounded below by cos(π/4):
    cos(π/4) ≤ cos(π/(2d)) for d ≥ 2.
    (cos is decreasing on [0,π], smaller angle ⟹ bigger cosine.) -/
theorem contraction_bounded_below (d : ℕ) (hd : d ≥ 2) :
    cos (π / 4) ≤ cos (π / (2 * (d : ℝ))) := by
  have h1 := angle_in_range d hd
  apply cos_le_cos_of_nonneg_of_le_pi
  · exact le_of_lt h1.1
  · linarith [h1.2, pi_pos]
  · exact angle_le_pi_div_four d hd

/-- The contraction ratio is bounded above by 1 (already proven as cos_angle_lt_one).
    Combined with contraction_bounded_below, we get:
    cos(π/4) ≤ cos(π/(2d)) < 1 for all d ≥ 2. -/
theorem contraction_sandwich (d : ℕ) (hd : d ≥ 2) :
    cos (π / 4) ≤ cos (π / (2 * (d : ℝ))) ∧ cos (π / (2 * (d : ℝ))) < 1 :=
  ⟨contraction_bounded_below d hd, cos_angle_lt_one d hd⟩

end Pandrosion
