/-
  Universitas Pandrosion — Lean 4 Formalization
  HALLEY COMPARISON MODULE

  Proves that the Pandrosion cube-root iteration P(s) = s(s³+4X)/(3s³+2X)
  is algebraically DISTINCT from both Newton's method N(s) = (2s³+X)/(3s²)
  and Halley's method H(s) = s(s³+2X)/(2s³+X) for the equation s³ = X.

  Main results:
    1. pandrosion_newton_cross: P and N differ by -(3s³-2X)(s³-X)
    2. pandrosion_halley_cross: P and H differ by -s⁴(s³-X)
    3. derivative_numerator: P'(r) numerator = -5r⁶ when r³ = X
    4. derivative_denominator: P'(r) denominator = 25r⁶ when r³ = X
    5. pandrosion_derivative_value: P'(r) = -1/5 for all cube roots r

  All proofs are pure polynomial identities verified by `ring`.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §200. Three cube-root iterations -/

/-- The Pandrosion iteration for cube roots: P(s) = s(s³+4X)/(3s³+2X). -/
noncomputable def pandrosion_P (X s : ℝ) : ℝ :=
  s * (s ^ 3 + 4 * X) / (3 * s ^ 3 + 2 * X)

/-- Newton's iteration for cube roots: N(s) = (2s³+X)/(3s²). -/
noncomputable def newton_N (X s : ℝ) : ℝ :=
  (2 * s ^ 3 + X) / (3 * s ^ 2)

/-- Halley's iteration for cube roots: H(s) = s(s³+2X)/(2s³+X). -/
noncomputable def halley_H (X s : ℝ) : ℝ :=
  s * (s ^ 3 + 2 * X) / (2 * s ^ 3 + X)

/-! ## §201. Pandrosion ≠ Newton

The cross-multiplication identity shows that P(s) and N(s) differ by
an exact polynomial factor proportional to (s³ - X).
When s³ = X (at the root), the difference vanishes, confirming
they agree at the fixed point. Away from the root, they differ. -/

/-- **Pandrosion − Newton: cross-multiplication form.**
    s(s³+4X)·(3s²) − (2s³+X)·(3s³+2X) = −(3s³−2X)(s³−X). -/
theorem pandrosion_newton_cross (X s : ℝ) :
    s * (s ^ 3 + 4 * X) * (3 * s ^ 2) -
    (2 * s ^ 3 + X) * (3 * s ^ 3 + 2 * X) =
    -(3 * s ^ 3 - 2 * X) * (s ^ 3 - X) := by
  ring

/-! ## §202. Pandrosion ≠ Halley

This is the key algebraic refutation: the Pandrosion iteration
is NOT Halley's method applied to s³ − X = 0. The explicit
polynomial difference is −s⁴(s³ − X), which is nonzero whenever
s ≠ 0 and s³ ≠ X. -/

/-- **Pandrosion − Halley: cross-multiplication form.**
    s(s³+4X)·(2s³+X) − s(s³+2X)·(3s³+2X) = −s⁴(s³−X). -/
theorem pandrosion_halley_cross (X s : ℝ) :
    s * (s ^ 3 + 4 * X) * (2 * s ^ 3 + X) -
    s * (s ^ 3 + 2 * X) * (3 * s ^ 3 + 2 * X) =
    -(s ^ 4 * (s ^ 3 - X)) := by
  ring

/-- **Pandrosion and Halley agree ONLY at roots.**
    P(s) · denom_H = H(s) · denom_P iff s⁴(s³−X) = 0. -/
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

/-! ## §203. The Universal Derivative P'(r) = −1/5

For P(s) = U(s)/V(s) with U(s) = s⁴+4sX, V(s) = 3s³+2X:
  P'(s) = (U'(s)V(s) − U(s)V'(s)) / V(s)²

where U'(s) = 4s³+4X, V'(s) = 9s².

At the fixed point r satisfying r³ = X:
  U'(r)V(r) − U(r)V'(r) = −5r⁶
  V(r)² = 25r⁶
  P'(r) = −5r⁶ / 25r⁶ = −1/5

This is the UNIVERSAL linear convergence rate:
the Pandrosion iteration contracts errors by exactly 1/5
at every cube root, regardless of the target X.
The negative sign produces the characteristic alternating
approach pattern (overshooting/undershooting). -/

/-- **Derivative numerator at root.** U'(r)V(r) − U(r)V'(r) = −5r⁶. -/
theorem derivative_numerator (X r : ℝ) (hX : r ^ 3 = X) :
    (4 * r ^ 3 + 4 * X) * (3 * r ^ 3 + 2 * X) -
    (r ^ 4 + 4 * r * X) * (9 * r ^ 2) =
    -(5 * r ^ 6) := by
  subst hX; ring

/-- **Derivative denominator at root.** V(r)² = 25r⁶. -/
theorem derivative_denominator (X r : ℝ) (hX : r ^ 3 = X) :
    (3 * r ^ 3 + 2 * X) ^ 2 = 25 * r ^ 6 := by
  subst hX; ring

/-- **The universal convergence rate: P'(r) = −1/5.**
    For ANY cube root r (with r³ = X and r ≠ 0),
    the derivative of the Pandrosion iteration at r equals −1/5.

    This is independent of X: whether computing ∛2, ∛1000, or ∛π,
    the local contraction rate is always exactly 1/5. -/
theorem pandrosion_derivative_value (X r : ℝ) (hX : r ^ 3 = X) (hr : r ≠ 0) :
    ((4 * r ^ 3 + 4 * X) * (3 * r ^ 3 + 2 * X) -
     (r ^ 4 + 4 * r * X) * (9 * r ^ 2)) /
    ((3 * r ^ 3 + 2 * X) ^ 2) = -(1 : ℝ) / 5 := by
  rw [derivative_numerator X r hX, derivative_denominator X r hX]
  have hr6 : r ^ 6 ≠ 0 := pow_ne_zero 6 hr
  have h25 : (25 : ℝ) * r ^ 6 ≠ 0 := mul_ne_zero (by norm_num) hr6
  rw [div_eq_div_iff h25 (by norm_num : (5 : ℝ) ≠ 0)]
  ring

/-! ## §204. Corollaries

Since P'(r) = −1/5 ≠ 0, the convergence is strictly LINEAR.
This is SLOWER than Newton (quadratic, P'(r) = 0) and Halley (cubic).

Since P'(r) = −1/5 ≠ 1, Aitken's Δ² extrapolation (Steffensen's
method) can be applied to accelerate convergence from linear to
at least quadratic. See SteffensenAcceleration.lean. -/

/-- The absolute contraction rate |P'(r)| = 1/5 < 1, confirming convergence. -/
theorem rate_abs_lt_one : |(-1 : ℝ) / 5| < 1 := by
  simp [abs_of_neg (by norm_num : (-1 : ℝ) / 5 < 0)]; norm_num

/-! ## §205. Newton's derivative vanishes at the root (for comparison)

This shows Newton has QUADRATIC convergence (P' = 0 at root),
while Pandrosion has only LINEAR convergence (P' = −1/5 ≠ 0). -/

/-- **Newton's derivative numerator at root vanishes.**
    N'(r) numerator = (6r²)(3r²) − (2r³+X)(6r) = 0 when r³ = X. -/
theorem newton_derivative_numerator (X r : ℝ) (hX : r ^ 3 = X) :
    (6 * r ^ 2) * (3 * r ^ 2) - (2 * r ^ 3 + X) * (6 * r) = 0 := by
  subst hX; ring

end Pandrosion
