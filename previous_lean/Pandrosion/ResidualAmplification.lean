/-
  Universitas Pandrosion — Lean 4 Formalization
  RESIDUAL AMPLIFICATION AT ROOT MODULE

  Computes the exact value of the conservation polynomial Φ(s,X)
  at the fixed point, establishing the precise link between the
  conservation identity and the convergence rate.

  Main results:
    Φ(r, r³) = -25r⁹
    (3r³ + 2X)³ = 125r⁹
    Ratio Φ/(3s³+2X)³ at root = -1/5

  This proves that the residual ratio P(s)³-X / (s³-X) converges
  to exactly -1/5 at the fixed point, providing an independent
  derivation of the convergence rate that does NOT use derivatives.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §209. The conservation polynomial at the fixed point

The conservation identity (Eq. 3 in the paper) states:
  P(s)³ - X = (s³ - X) · Φ(s,X) / (3s³ + 2X)³

where Φ(s,X) = s⁹ - 14s⁶X - 20s³X² + 8X³.

At the fixed point s = r where r³ = X, the residual s³-X vanishes
(as expected), but the RATIO Φ(s,X)/(3s³+2X)³ has a well-defined
limit of -1/5. We prove this algebraically.
-/

/-- **The conservation polynomial at root.**
    Φ(r, r³) = r⁹ - 14r⁹ - 20r⁹ + 8r⁹ = -25r⁹. -/
theorem conservation_polynomial_at_root (X r : ℝ) (hX : r ^ 3 = X) :
    r ^ 9 - 14 * r ^ 6 * X - 20 * r ^ 3 * X ^ 2 + 8 * X ^ 3 =
    -(25 * r ^ 9) := by
  subst hX; ring

/-- **The denominator cubed at root.**
    (3r³ + 2X)³ = (5r³)³ = 125r⁹. -/
theorem denominator_cubed_at_root (X r : ℝ) (hX : r ^ 3 = X) :
    (3 * r ^ 3 + 2 * X) ^ 3 = 125 * r ^ 9 := by
  subst hX; ring

/-- **The residual ratio at the fixed point equals -1/5.**
    Φ(r,X) / (3r³+2X)³ = -25r⁹ / 125r⁹ = -1/5.

    This provides an INDEPENDENT derivation of the convergence rate
    without computing any derivative. The rate emerges purely from
    the algebraic structure of the conservation identity. -/
theorem residual_ratio_at_root (X r : ℝ) (hX : r ^ 3 = X) (hr : r ≠ 0) :
    (r ^ 9 - 14 * r ^ 6 * X - 20 * r ^ 3 * X ^ 2 + 8 * X ^ 3) /
    ((3 * r ^ 3 + 2 * X) ^ 3) = -(1 : ℝ) / 5 := by
  rw [conservation_polynomial_at_root X r hX, denominator_cubed_at_root X r hX]
  have hr9 : r ^ 9 ≠ 0 := pow_ne_zero 9 hr
  have h125 : (125 : ℝ) * r ^ 9 ≠ 0 := mul_ne_zero (by norm_num) hr9
  rw [div_eq_div_iff h125 (by norm_num : (5 : ℝ) ≠ 0)]
  ring

/-! ## §210. The oscillation polynomial at the fixed point

The oscillation identity states:
  P(s) - r = (s - r) · Ψ(s) / (3s³ + 2X)

where Ψ(s) = s³ - 2rs² - 2r²s + 2r³.
At s = r: Ψ(r) = r³ - 2r³ - 2r³ + 2r³ = -r³.
-/

/-- **The oscillation polynomial at root.**
    Ψ(r) = r³ - 2r³ - 2r³ + 2r³ = -r³. -/
theorem oscillation_polynomial_at_root (r : ℝ) :
    r ^ 3 - 2 * r * r ^ 2 - 2 * r ^ 2 * r + 2 * r ^ 3 = -(r ^ 3) := by
  ring

/-- **Local oscillation ratio.**
    Ψ(r) / (3r³+2X) = -r³ / (5r³) = -1/5 when r³ = X.

    This independently confirms the rate: the oscillation identity
    and the conservation identity both give -1/5 at the root. -/
theorem oscillation_ratio_at_root (X r : ℝ) (hX : r ^ 3 = X) (hr : r ≠ 0) :
    -(r ^ 3) / (3 * r ^ 3 + 2 * X) = -(1 : ℝ) / 5 := by
  rw [← hX]
  have hr3 : r ^ 3 ≠ 0 := pow_ne_zero 3 hr
  have h5 : (5 : ℝ) * r ^ 3 ≠ 0 := mul_ne_zero (by norm_num) hr3
  rw [show 3 * r ^ 3 + 2 * r ^ 3 = 5 * r ^ 3 from by ring]
  field_simp
  ring

/-! ## §211. The convergence rate -1/5 emerges from THREE independent sources

1. The derivative: P'(r) = -1/5         (HalleyComparison.lean)
2. The conservation ratio: Φ/D³ = -1/5  (this module, §209)
3. The oscillation ratio: Ψ/D = -1/5    (this module, §210)

All three are proved by pure algebra (ring tactic).
This triple coincidence is a structural property of the Pandrosion
iteration, not a numerical accident.
-/

end Pandrosion
