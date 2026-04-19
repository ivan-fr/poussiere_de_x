/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP X: DERIVATIVE OF Sp — HasDerivAt for the geometric partial sum

  Sp(s) = Σ_{k<p} s^k  ⟹  Sp'(s) = Σ_{k=1}^{p-1} k · s^{k-1}

  Using Mathlib's HasDerivAt.sum + hasDerivAt_pow.

  Reference: pandrosion_master.tex, §74 (Differential Structure)
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Tactic
import Pandrosion.Deep

open Finset BigOperators

namespace Pandrosion

/-! ## §74. Derivative of Sp -/

/-- The formal derivative of Sp: Sp'(s) = Σ_{k<p} k · s^{k-1}. -/
noncomputable def Sp' (p : ℕ) (s : ℝ) : ℝ :=
  ∑ k in range p, (k : ℝ) * s ^ (k - 1)

/-- **HasDerivAt for Sp.** The geometric partial sum Σ s^k is differentiable
    with derivative Σ k · s^{k-1}.

    Proof: apply HasDerivAt.sum to the sum of hasDerivAt_pow terms. -/
theorem Sp_hasDerivAt (p : ℕ) (s : ℝ) :
    HasDerivAt (Sp p) (Sp' p s) s := by
  unfold Sp Sp'
  exact HasDerivAt.sum (fun k _ => hasDerivAt_pow k s)

/-- **Sp is differentiable at every point.** -/
theorem Sp_differentiableAt (p : ℕ) (s : ℝ) :
    DifferentiableAt ℝ (Sp p) s :=
  (Sp_hasDerivAt p s).differentiableAt

/-- **The derivative of Sp equals Sp'.** -/
theorem Sp_deriv (p : ℕ) (s : ℝ) :
    deriv (Sp p) s = Sp' p s :=
  (Sp_hasDerivAt p s).deriv

/-! ## §75. Properties of Sp' -/

/-- **Sp' is non-negative for s ≥ 0.** -/
theorem Sp'_nonneg (p : ℕ) (s : ℝ) (hs : s ≥ 0) :
    Sp' p s ≥ 0 := by
  unfold Sp'
  apply Finset.sum_nonneg
  intro k _
  exact mul_nonneg (Nat.cast_nonneg' k) (pow_nonneg hs _)

/-- **Sp' is strictly positive for s ≥ 0 and p ≥ 2.**
    The k=1 term contributes 1 · s^0 = 1 > 0. -/
theorem Sp'_pos (p : ℕ) (hp : p ≥ 2) (s : ℝ) (hs : s ≥ 0) :
    Sp' p s > 0 := by
  unfold Sp'
  -- The k=1 term is 1 * s^0 = 1
  have h1mem : (1 : ℕ) ∈ range p := Finset.mem_range.mpr (by omega)
  have hle : (1 : ℝ) ≤ ∑ k in range p, (k : ℝ) * s ^ (k - 1) := by
    have := Finset.single_le_sum
      (fun (k : ℕ) (_ : k ∈ range p) => mul_nonneg (Nat.cast_nonneg' k) (pow_nonneg hs (k-1)))
      h1mem
    simp at this; exact this
  linarith

/-! ## §76. Derivative of h(s) = 1 - (x-1)/(x · Sp(s)) -/

/-- **HasDerivAt for the Pandrosion iteration h.**
    h'(s) = (x-1) · Sp'(s) / (x · Sp(s)²) -/
theorem h_hasDerivAt (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1) (s : ℝ) (hs : s ≥ 0) :
    HasDerivAt (pandrosion_h x p) ((x - 1) * Sp' p s / (x * (Sp p s) ^ 2)) s := by
  unfold pandrosion_h
  have hSp_ne : Sp p s ≠ 0 := ne_of_gt (Sp_pos p (by omega) s hs)
  have hx_ne : (x : ℝ) ≠ 0 := ne_of_gt hx
  have hxSp_ne : x * Sp p s ≠ 0 := mul_ne_zero hx_ne hSp_ne
  have hSp : HasDerivAt (Sp p) (Sp' p s) s := Sp_hasDerivAt p s
  -- d/ds [x * Sp(s)] = x * Sp'(s)
  have hxSp : HasDerivAt (fun s => x * Sp p s) (0 * Sp p s + x * Sp' p s) s :=
    (hasDerivAt_const s x).mul hSp
  -- d/ds [(x-1)/(x*Sp(s))]
  have hdiv : HasDerivAt (fun s => (x - 1) / (x * Sp p s))
      ((0 * (x * Sp p s) - (x - 1) * (0 * Sp p s + x * Sp' p s)) / (x * Sp p s) ^ 2) s :=
    (hasDerivAt_const s (x - 1)).div hxSp hxSp_ne
  -- h(s) = 1 - (x-1)/(x·Sp(s))
  have hfinal : HasDerivAt (fun s => 1 - (x - 1) / (x * Sp p s))
      (0 - (0 * (x * Sp p s) - (x - 1) * (0 * Sp p s + x * Sp' p s)) / (x * Sp p s) ^ 2) s :=
    (hasDerivAt_const s 1).sub hdiv
  convert hfinal using 1
  field_simp
  ring

end Pandrosion
