/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP V: THE DERIVATIVE h'(s) and asymptotic convergence rate

  For p = 2: h(s) = 1 - (x-1)/(x(1+s))
  h'(s) = (x-1)/(x(1+s)²)
  At the fixed point s*: h'(s*) = (x-1)/(x(1+s*)²)

  This uses Mathlib's HasDerivAt calculus API.
  Reference: pandrosion_master.tex, Theorems 336, 405
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Tactic
import Pandrosion.Deep
import Pandrosion.Deep2

open Real

namespace Pandrosion

/-! ## §60. The Derivative of h for p = 2

The central computation:
h(s) = 1 - (x-1)/(x(1+s)) = 1 - (x-1)/x · 1/(1+s)

h'(s) = (x-1)/x · 1/(1+s)²

This is the ASYMPTOTIC CONTRACTION RATE.
-/

/-- **The derivative of g(s) = 1/(1+s) is -1/(1+s)².** -/
theorem hasDerivAt_inv_one_add (s : ℝ) (hs : 1 + s ≠ 0) :
    HasDerivAt (fun t => (1 + t)⁻¹) (-1 / (1 + s) ^ 2) s := by
  have h1 : HasDerivAt (fun t => 1 + t) 1 s := by
    have := (hasDerivAt_const s (1 : ℝ)).add (hasDerivAt_id s)
    simp at this; exact this
  exact h1.inv hs

/-- **The derivative of h for p = 2.**
    h(s) = 1 - (x-1)/(x(1+s))
    h'(s) = (x-1) / (x · (1+s)²)

    Proof: Write h(s) = 1 - c/(1+s) where c = (x-1)/x.
    Then h'(s) = c/(1+s)² = (x-1)/(x(1+s)²). -/
theorem h_deriv_p2 (x s : ℝ) (hx : x ≠ 0) (hs : (1 : ℝ) + s ≠ 0) :
    HasDerivAt (pandrosion_h x 2) ((x - 1) / (x * (1 + s) ^ 2)) s := by
  -- h(s) = 1 - (x-1)/(x·(1+s)) = 1 - ((x-1)/x) · (1+s)⁻¹
  -- Derivative of constant 1 is 0
  have hd1 : HasDerivAt (fun _ : ℝ => (1 : ℝ)) 0 s := hasDerivAt_const s 1
  -- Derivative of (x-1)/x · (1+t)⁻¹
  have hd2 : HasDerivAt (fun t => (x - 1) / x * (1 + t)⁻¹)
      ((x - 1) / x * (-1 / (1 + s) ^ 2)) s := by
    have := (hasDerivAt_const s ((x - 1) / x)).mul (hasDerivAt_inv_one_add s hs)
    simp at this; exact this
  -- h'(s) = 0 - [(x-1)/x · (-1/(1+s)²)]
  have hd : HasDerivAt (fun t => 1 - (x - 1) / x * (1 + t)⁻¹)
      (0 - (x - 1) / x * (-1 / (1 + s) ^ 2)) s := hd1.sub hd2
  -- Now convert: pandrosion_h x 2 = fun t => 1 - (x-1)/x * (1+t)⁻¹
  have heq : pandrosion_h x 2 = fun t => 1 - (x - 1) / x * (1 + t)⁻¹ := by
    ext t; simp [pandrosion_h, Sp2_eq]; field_simp
  rw [heq]
  have hval : (x - 1) / (x * (1 + s) ^ 2) = 0 - (x - 1) / x * (-1 / (1 + s) ^ 2) := by
    have hs2 : (1 + s) ^ 2 ≠ 0 := pow_ne_zero 2 hs
    field_simp
  rw [hval]
  exact hd

/-- **The derivative at the fixed point s* gives the asymptotic rate.**
    h'(s*) = (x-1)/(x(1+s*)²).

    For x = 4 (square root of 4):
    s* = 1/2, h'(1/2) = 3/(4·(3/2)²) = 3/9 = 1/3.

    This is strictly less than 1, proving linear convergence. -/
theorem h_deriv_at_fixpoint_p2 (x sstar : ℝ) (_hx_pos : x > 1) (hx : x ≠ 0)
    (hss_pos : sstar > 0) (_hss_lt : sstar < 1) :
    HasDerivAt (pandrosion_h x 2) ((x - 1) / (x * (1 + sstar) ^ 2)) sstar :=
  h_deriv_p2 x sstar hx (by linarith)

/-- **The asymptotic rate h'(s*) is in (0,1).**
    h'(s*) = (x-1)/(x(1+s*)²). Since s* > 0 and x > 1:
    x(1+s*)² > x > x-1, so h'(s*) < 1.
    h'(s*) > 0 since x-1 > 0 and x(1+s*)² > 0. -/
theorem asymptotic_rate_in_unit (x sstar : ℝ) (hx : x > 1) (hss : sstar > 0) :
    0 < (x - 1) / (x * (1 + sstar) ^ 2) ∧
    (x - 1) / (x * (1 + sstar) ^ 2) < 1 := by
  have hx_pos : x > 0 := by linarith
  have h1s : (1 : ℝ) + sstar > 1 := by linarith
  have h1s2 : (1 + sstar) ^ 2 > 1 := by nlinarith
  have hD : x * (1 + sstar) ^ 2 > 0 := by positivity
  constructor
  · exact div_pos (by linarith) hD
  · rw [div_lt_one hD]
    -- x - 1 < x · (1+s*)²
    -- x · (1+s*)² > x · 1 = x > x - 1
    calc x - 1 < x := by linarith
      _ = x * 1 := (mul_one x).symm
      _ < x * (1 + sstar) ^ 2 := by
          exact mul_lt_mul_of_pos_left h1s2 hx_pos

/-! ## §61. Contraction for p = 3 (Cube Root Case) -/

/-- **S₃(s) = 1 + s + s² for p = 3.** -/
theorem Sp3_eq (s : ℝ) : Sp 3 s = 1 + s + s ^ 2 := by
  unfold Sp; simp [Finset.sum_range_succ]

/-- **h for p = 3: h(s) = 1 - (x-1)/(x(1+s+s²)).** -/
theorem h_p3 (x s : ℝ) (_hx : x ≠ 0) (_hs : 1 + s + s ^ 2 ≠ 0) :
    pandrosion_h x 3 s = 1 - (x - 1) / (x * (1 + s + s ^ 2)) := by
  unfold pandrosion_h
  have : Sp 3 s = 1 + s + s ^ 2 := Sp3_eq s
  rw [this]

/-- **Fixed point for cube root: h(s) = s ⟺ s³ = 1/x.** -/
theorem cube_root_fixpoint (x : ℝ) (hx : x > 0) (s : ℝ) (hs : 0 ≤ s) (hs1 : s ≠ 1) :
    pandrosion_h x 3 s = s ↔ s ^ 3 = 1 / x :=
  fixed_point_iff x hx 3 (by omega) s hs hs1

end Pandrosion
