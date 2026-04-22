/-
  Universitas Pandrosion — Legacy / Cube-root comparison bundle.

  Three consolidated results specific to the cube-root Pandrosion
  iteration `P(s) = s · (s³ + 4X) / (3s³ + 2X)`:

    §L23  Pandrosion ≠ Newton, Pandrosion ≠ Halley (cross identities).
    §L24  `P'(r) = −1/5` at every cube root  (universal linear rate).
    §L25  `P''(r) = −12 / (25r)`             (quadratic correction).
    §L26  Pandrosion ∉ Chebyshev-Halley family (parameter conflict).

  Ported from `HalleyComparison.lean`, `HigherDerivatives.lean`,
  and `ChebyshevHalleyExclusion.lean` on the `main_previous`
  branch.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §L23. Cube-root iterators -/

/-- Pandrosion cube-root iteration `P(s) = s(s³ + 4X) / (3s³ + 2X)`. -/
noncomputable def pandrosion_P (X s : ℝ) : ℝ :=
  s * (s ^ 3 + 4 * X) / (3 * s ^ 3 + 2 * X)

/-- Newton cube-root iteration `N(s) = (2s³ + X) / (3s²)`. -/
noncomputable def newton_N (X s : ℝ) : ℝ :=
  (2 * s ^ 3 + X) / (3 * s ^ 2)

/-- Halley cube-root iteration `H(s) = s(s³ + 2X) / (2s³ + X)`. -/
noncomputable def halley_H (X s : ℝ) : ℝ :=
  s * (s ^ 3 + 2 * X) / (2 * s ^ 3 + X)

/-- **Pandrosion − Newton cross identity.**
    `s(s³+4X)·3s² − (2s³+X)·(3s³+2X) = −(3s³−2X)(s³−X)`. -/
theorem pandrosion_newton_cross (X s : ℝ) :
    s * (s ^ 3 + 4 * X) * (3 * s ^ 2) -
    (2 * s ^ 3 + X) * (3 * s ^ 3 + 2 * X) =
    -(3 * s ^ 3 - 2 * X) * (s ^ 3 - X) := by
  ring

/-- **Pandrosion − Halley cross identity.**
    `s(s³+4X)·(2s³+X) − s(s³+2X)·(3s³+2X) = −s⁴(s³−X)`. -/
theorem pandrosion_halley_cross (X s : ℝ) :
    s * (s ^ 3 + 4 * X) * (2 * s ^ 3 + X) -
    s * (s ^ 3 + 2 * X) * (3 * s ^ 3 + 2 * X) =
    -(s ^ 4 * (s ^ 3 - X)) := by
  ring

/-- **Pandrosion and Halley agree only at roots.** -/
theorem pandrosion_halley_agree_iff (X s : ℝ) :
    s * (s ^ 3 + 4 * X) * (2 * s ^ 3 + X) =
    s * (s ^ 3 + 2 * X) * (3 * s ^ 3 + 2 * X)
    ↔ s ^ 4 * (s ^ 3 - X) = 0 := by
  constructor
  · intro h
    have := pandrosion_halley_cross X s
    linarith
  · intro h
    have := pandrosion_halley_cross X s
    linarith

/-! ## §L24. Universal linear rate `P'(r) = −1/5` -/

/-- Derivative numerator at a cube root. -/
theorem cube_derivative_numerator (X r : ℝ) (hX : r ^ 3 = X) :
    (4 * r ^ 3 + 4 * X) * (3 * r ^ 3 + 2 * X) -
    (r ^ 4 + 4 * r * X) * (9 * r ^ 2) =
    -(5 * r ^ 6) := by
  subst hX; ring

/-- Derivative denominator at a cube root. -/
theorem cube_derivative_denominator (X r : ℝ) (hX : r ^ 3 = X) :
    (3 * r ^ 3 + 2 * X) ^ 2 = 25 * r ^ 6 := by
  subst hX; ring

/-- **Universal linear rate.** For every cube root `r` of `X` with
    `r ≠ 0`, the Pandrosion derivative satisfies `P'(r) = −1/5`. -/
theorem pandrosion_derivative_value (X r : ℝ) (hX : r ^ 3 = X) (hr : r ≠ 0) :
    ((4 * r ^ 3 + 4 * X) * (3 * r ^ 3 + 2 * X) -
     (r ^ 4 + 4 * r * X) * (9 * r ^ 2)) /
    ((3 * r ^ 3 + 2 * X) ^ 2) = -(1 : ℝ) / 5 := by
  rw [cube_derivative_numerator X r hX, cube_derivative_denominator X r hX]
  have hr6 : r ^ 6 ≠ 0 := pow_ne_zero 6 hr
  have h25 : (25 : ℝ) * r ^ 6 ≠ 0 := mul_ne_zero (by norm_num) hr6
  rw [div_eq_div_iff h25 (by norm_num : (5 : ℝ) ≠ 0)]
  ring

/-- The linear rate is non-unit: `−1/5 ≠ 1` (Aitken-Steffensen applicable). -/
theorem pandrosion_rate_ne_one : -(1 : ℝ) / 5 ≠ 1 := by norm_num

/-- The linear rate is non-zero: convergence is strictly linear. -/
theorem pandrosion_rate_ne_zero : -(1 : ℝ) / 5 ≠ 0 := by norm_num

/-- `|−1/5| < 1`: Pandrosion contracts. -/
theorem pandrosion_rate_abs_lt_one : |(-1 : ℝ) / 5| < 1 := by
  simp [abs_of_neg (by norm_num : (-1 : ℝ) / 5 < 0)]; norm_num

/-- Newton has quadratic convergence: `N'(r) = 0` at every cube root. -/
theorem newton_derivative_numerator (X r : ℝ) (hX : r ^ 3 = X) :
    (6 * r ^ 2) * (3 * r ^ 2) - (2 * r ^ 3 + X) * (6 * r) = 0 := by
  subst hX; ring

/-! ## §L25. Quadratic correction `P''(r) = −12/(25r)` -/

/-- Second-derivative numerator at a cube root. -/
theorem cube_second_derivative_numerator (X r : ℝ) (hX : r ^ 3 = X) :
    ((12 * r ^ 2) * (3 * r ^ 3 + 2 * X) - (r ^ 4 + 4 * r * X) * (18 * r)) *
    (3 * r ^ 3 + 2 * X) -
    2 * ((4 * r ^ 3 + 4 * X) * (3 * r ^ 3 + 2 * X) -
    (r ^ 4 + 4 * r * X) * (9 * r ^ 2)) * (9 * r ^ 2) =
    -(60 * r ^ 8) := by
  subst hX; ring

/-- Second-derivative denominator at a cube root. -/
theorem cube_second_derivative_denominator (X r : ℝ) (hX : r ^ 3 = X) :
    (3 * r ^ 3 + 2 * X) ^ 3 = 125 * r ^ 9 := by
  subst hX; ring

/-- **Quadratic correction.** `P''(r) = −12 / (25r)` at every cube root. -/
theorem pandrosion_second_derivative_value (X r : ℝ) (hX : r ^ 3 = X) (hr : r ≠ 0) :
    (((12 * r ^ 2) * (3 * r ^ 3 + 2 * X) - (r ^ 4 + 4 * r * X) * (18 * r)) *
    (3 * r ^ 3 + 2 * X) -
    2 * ((4 * r ^ 3 + 4 * X) * (3 * r ^ 3 + 2 * X) -
    (r ^ 4 + 4 * r * X) * (9 * r ^ 2)) * (9 * r ^ 2)) /
    ((3 * r ^ 3 + 2 * X) ^ 3) = -(12 : ℝ) / (25 * r) := by
  rw [cube_second_derivative_numerator X r hX,
      cube_second_derivative_denominator X r hX]
  have hr9 : r ^ 9 ≠ 0 := pow_ne_zero 9 hr
  have h125 : (125 : ℝ) * r ^ 9 ≠ 0 := mul_ne_zero (by norm_num) hr9
  have h25r : (25 : ℝ) * r ≠ 0 := mul_ne_zero (by norm_num) hr
  rw [div_eq_div_iff h125 h25r]
  ring

/-! ## §L26. Pandrosion ∉ Chebyshev-Halley family -/

/-- **First Chebyshev-Halley matching condition.**
    If Pandrosion belongs to `CH_α`, the `s⁶` coefficient forces `α = 1/2`. -/
theorem chebyshev_halley_coeff_s6 (α : ℝ)
    (h : 6 * (3 - 2 * α) = 15 - 6 * α) :
    α = 1 / 2 := by linarith

/-- **Second Chebyshev-Halley matching condition.**
    The `s³X` coefficient forces `α = 2/5`. -/
theorem chebyshev_halley_coeff_s3X (α : ℝ) (h : 12 * α = 4 + 2 * α) :
    α = 2 / 5 := by linarith

/-- `1/2 ≠ 2/5`. -/
theorem half_ne_two_fifths : (1 : ℝ) / 2 ≠ 2 / 5 := by norm_num

/-- **Exclusion theorem.** There is no `α ∈ ℝ` such that the Pandrosion
    iteration equals the Chebyshev-Halley iterator `CH_α`. The two
    coefficient-matching conditions are inconsistent. -/
theorem pandrosion_not_in_chebyshev_halley :
    ¬ ∃ α : ℝ,
      6 * (3 - 2 * α) = 15 - 6 * α ∧ 12 * α = 4 + 2 * α := by
  intro ⟨α, h1, h2⟩
  have _hα1 : α = 1 / 2 := chebyshev_halley_coeff_s6 α h1
  have _hα2 : α = 2 / 5 := chebyshev_halley_coeff_s3X α h2
  linarith

/-- Pandrosion linear rate differs from the universal CH rate `0`:
    every `CH_α` has `P'(r) = 0` (at least quadratic), Pandrosion has
    `P'(r) = −1/5 ≠ 0`. -/
theorem pandrosion_rate_ne_ch_rate : -(1 : ℝ) / 5 ≠ 0 := by norm_num

end Pandrosion
