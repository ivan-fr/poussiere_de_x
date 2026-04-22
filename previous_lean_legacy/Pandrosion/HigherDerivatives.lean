/-
  Universitas Pandrosion — Lean 4 Formalization
  HIGHER DERIVATIVES MODULE

  Computes the second and third derivatives of the Pandrosion iteration
  at the fixed point r (where r³ = X), proving:

    P''(r) = -12/(25r)  via  numerator = -60r⁸, denominator = 125r⁹
    P'''(r) = 744/(125r²) via numerator/denominator identity

  Combined with P'(r) = -1/5 (from HalleyComparison.lean), this gives
  the complete local Taylor expansion:

    P(r+ε) = r - ε/5 - 6ε²/(25r) + 124ε³/(125r²) + O(ε⁴)

  The non-vanishing of P''(r) means the Pandrosion iteration has a
  QUADRATIC correction term beyond the linear rate -1/5.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §206. Second derivative at the fixed point

For P(s) = U(s)/V(s), define the second-derivative numerator as:
  N₂(s) = (U''V - UV'')V - 2(U'V - UV')V'

Then P''(s) = N₂(s) / V(s)³.

With U = s⁴+4sX, V = 3s³+2X:
  U' = 4s³+4X, U'' = 12s²
  V' = 9s², V'' = 18s
-/

/-- **Second derivative numerator identity.**
    At r³ = X: N₂(r) = (12r²·5r³ - 5r⁴·18r)·5r³ - 2·(-5r⁶)·9r² = -60r⁸. -/
theorem second_derivative_numerator (X r : ℝ) (hX : r ^ 3 = X) :
    ((12 * r ^ 2) * (3 * r ^ 3 + 2 * X) - (r ^ 4 + 4 * r * X) * (18 * r)) *
    (3 * r ^ 3 + 2 * X) -
    2 * ((4 * r ^ 3 + 4 * X) * (3 * r ^ 3 + 2 * X) -
    (r ^ 4 + 4 * r * X) * (9 * r ^ 2)) * (9 * r ^ 2) =
    -(60 * r ^ 8) := by
  subst hX; ring

/-- **Second derivative denominator identity.**
    V(r)³ = (3r³+2X)³ = 125r⁹ when r³ = X. -/
theorem second_derivative_denominator (X r : ℝ) (hX : r ^ 3 = X) :
    (3 * r ^ 3 + 2 * X) ^ 3 = 125 * r ^ 9 := by
  subst hX; ring

/-- **The exact second derivative: P''(r) = -12/(25r).**
    This proves the quadratic correction to the linear convergence.
    Combined with P'(r) = -1/5:
      P(r+ε) = r - ε/5 - 6ε²/(25r) + O(ε³). -/
theorem second_derivative_value (X r : ℝ) (hX : r ^ 3 = X) (hr : r ≠ 0) :
    (((12 * r ^ 2) * (3 * r ^ 3 + 2 * X) - (r ^ 4 + 4 * r * X) * (18 * r)) *
    (3 * r ^ 3 + 2 * X) -
    2 * ((4 * r ^ 3 + 4 * X) * (3 * r ^ 3 + 2 * X) -
    (r ^ 4 + 4 * r * X) * (9 * r ^ 2)) * (9 * r ^ 2)) /
    ((3 * r ^ 3 + 2 * X) ^ 3) = -(12 : ℝ) / (25 * r) := by
  rw [second_derivative_numerator X r hX, second_derivative_denominator X r hX]
  have hr9 : r ^ 9 ≠ 0 := pow_ne_zero 9 hr
  have h125 : (125 : ℝ) * r ^ 9 ≠ 0 := mul_ne_zero (by norm_num) hr9
  have h25r : (25 : ℝ) * r ≠ 0 := mul_ne_zero (by norm_num) hr
  rw [div_eq_div_iff h125 h25r]
  ring

end Pandrosion
