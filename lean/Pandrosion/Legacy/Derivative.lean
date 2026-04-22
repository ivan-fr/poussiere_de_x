/-
  Universitas Pandrosion — Legacy / Derivative.

  `HasDerivAt` facts for `Sp` and the Pandrosion iterator `h x p`.

    • `Sp_hasDerivAt` — `Sp` is everywhere differentiable with
      derivative `Sp'(s) = Σ k·s^{k-1}`.
    • `h_hasDerivAt` — generic derivative at any `s ≥ 0`.
    • `h_deriv_p2`, `h_deriv_at_fixpoint_p2`, `asymptotic_rate_in_unit`
      — the explicit `p = 2` derivative formula and the asymptotic
      convergence rate `h'(s*) ∈ (0, 1)`.
    • `Sp3_eq`, `h_p3`, `cube_root_fixpoint` — the cube-root shape.

  Ported from `DerivativeSp.lean` + `DerivativeRate.lean` on the
  `main_previous` branch.
-/

import Pandrosion.Legacy.Basic
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.MeanValue

namespace Pandrosion

open Finset BigOperators

/-! ## §L14. Derivative of `Sp` -/

/-- Formal derivative `Sp'(s) = Σ_{k<p} k · s^{k-1}`. -/
noncomputable def Sp' (p : ℕ) (s : ℝ) : ℝ :=
  ∑ k in range p, (k : ℝ) * s ^ (k - 1)

/-- `HasDerivAt (Sp p) (Sp' p s) s`. -/
theorem Sp_hasDerivAt (p : ℕ) (s : ℝ) :
    HasDerivAt (Sp p) (Sp' p s) s := by
  unfold Sp Sp'
  exact HasDerivAt.sum (fun k _ => hasDerivAt_pow k s)

/-- `Sp` is everywhere differentiable. -/
theorem Sp_differentiableAt (p : ℕ) (s : ℝ) :
    DifferentiableAt ℝ (Sp p) s :=
  (Sp_hasDerivAt p s).differentiableAt

/-- `deriv (Sp p) s = Sp' p s`. -/
theorem Sp_deriv (p : ℕ) (s : ℝ) :
    deriv (Sp p) s = Sp' p s :=
  (Sp_hasDerivAt p s).deriv

/-- `Sp'(s) ≥ 0` for `s ≥ 0`. -/
theorem Sp'_nonneg (p : ℕ) (s : ℝ) (hs : s ≥ 0) :
    Sp' p s ≥ 0 := by
  unfold Sp'
  apply Finset.sum_nonneg
  intro k _
  exact mul_nonneg (Nat.cast_nonneg' k) (pow_nonneg hs _)

/-- `Sp'(s) > 0` for `s ≥ 0` and `p ≥ 2` (the `k = 1` term is `1`). -/
theorem Sp'_pos (p : ℕ) (hp : p ≥ 2) (s : ℝ) (hs : s ≥ 0) :
    Sp' p s > 0 := by
  unfold Sp'
  have h1mem : (1 : ℕ) ∈ range p := Finset.mem_range.mpr (by omega)
  have hle : (1 : ℝ) ≤ ∑ k in range p, (k : ℝ) * s ^ (k - 1) := by
    have := Finset.single_le_sum
      (fun (k : ℕ) (_ : k ∈ range p) =>
        mul_nonneg (Nat.cast_nonneg' k) (pow_nonneg hs (k - 1)))
      h1mem
    simp at this; exact this
  linarith

/-! ## §L15. Derivative of the Pandrosion map `h` -/

/-- **Generic derivative of `h`.** For `x > 0`, `p ≥ 1`, `s ≥ 0`:
    `h'(s) = (x − 1) · Sp'(s) / (x · Sp(s)²)`. -/
theorem h_hasDerivAt (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s : ℝ) (hs : s ≥ 0) :
    HasDerivAt (pandrosion_h x p)
      ((x - 1) * Sp' p s / (x * (Sp p s) ^ 2)) s := by
  unfold pandrosion_h
  have hSp_ne : Sp p s ≠ 0 := ne_of_gt (Sp_pos p hp s hs)
  have hx_ne : (x : ℝ) ≠ 0 := ne_of_gt hx
  have hxSp_ne : x * Sp p s ≠ 0 := mul_ne_zero hx_ne hSp_ne
  have hSp : HasDerivAt (Sp p) (Sp' p s) s := Sp_hasDerivAt p s
  have hxSp : HasDerivAt (fun s => x * Sp p s)
      (0 * Sp p s + x * Sp' p s) s :=
    (hasDerivAt_const s x).mul hSp
  have hdiv : HasDerivAt (fun s => (x - 1) / (x * Sp p s))
      ((0 * (x * Sp p s) - (x - 1) *
        (0 * Sp p s + x * Sp' p s)) / (x * Sp p s) ^ 2) s :=
    (hasDerivAt_const s (x - 1)).div hxSp hxSp_ne
  have hfinal : HasDerivAt (fun s => 1 - (x - 1) / (x * Sp p s))
      (0 - (0 * (x * Sp p s) - (x - 1) *
        (0 * Sp p s + x * Sp' p s)) / (x * Sp p s) ^ 2) s :=
    (hasDerivAt_const s 1).sub hdiv
  convert hfinal using 1
  field_simp
  ring

/-! ## §L16. The `p = 2` derivative and asymptotic rate -/

/-- `HasDerivAt (fun t => (1 + t)⁻¹) (−1/(1+s)²) s` when `1 + s ≠ 0`. -/
theorem hasDerivAt_inv_one_add (s : ℝ) (hs : 1 + s ≠ 0) :
    HasDerivAt (fun t => (1 + t)⁻¹) (-1 / (1 + s) ^ 2) s := by
  have h1 : HasDerivAt (fun t => 1 + t) 1 s := by
    have := (hasDerivAt_const s (1 : ℝ)).add (hasDerivAt_id s)
    simp at this; exact this
  exact h1.inv hs

/-- **`p = 2` derivative.** `h'(s) = (x − 1) / (x · (1 + s)²)`. -/
theorem h_deriv_p2 (x s : ℝ) (hx : x ≠ 0) (hs : (1 : ℝ) + s ≠ 0) :
    HasDerivAt (pandrosion_h x 2)
      ((x - 1) / (x * (1 + s) ^ 2)) s := by
  have hd1 : HasDerivAt (fun _ : ℝ => (1 : ℝ)) 0 s := hasDerivAt_const s 1
  have hd2 : HasDerivAt (fun t => (x - 1) / x * (1 + t)⁻¹)
      ((x - 1) / x * (-1 / (1 + s) ^ 2)) s := by
    have := (hasDerivAt_const s ((x - 1) / x)).mul
      (hasDerivAt_inv_one_add s hs)
    simp at this; exact this
  have hd : HasDerivAt (fun t => 1 - (x - 1) / x * (1 + t)⁻¹)
      (0 - (x - 1) / x * (-1 / (1 + s) ^ 2)) s := hd1.sub hd2
  have heq : pandrosion_h x 2 = fun t => 1 - (x - 1) / x * (1 + t)⁻¹ := by
    ext t
    have hSp : Sp 2 t = 1 + t := by
      unfold Sp; simp [Finset.sum_range_succ]
    simp [pandrosion_h, hSp]
    field_simp
  rw [heq]
  have hval : (x - 1) / (x * (1 + s) ^ 2)
      = 0 - (x - 1) / x * (-1 / (1 + s) ^ 2) := by
    have hs2 : (1 + s) ^ 2 ≠ 0 := pow_ne_zero 2 hs
    field_simp
  rw [hval]
  exact hd

/-- Derivative at a positive fixed point `sstar` (`p = 2`). -/
theorem h_deriv_at_fixpoint_p2 (x sstar : ℝ) (_hx_pos : x > 1) (hx : x ≠ 0)
    (hss_pos : sstar > 0) (_hss_lt : sstar < 1) :
    HasDerivAt (pandrosion_h x 2)
      ((x - 1) / (x * (1 + sstar) ^ 2)) sstar :=
  h_deriv_p2 x sstar hx (by linarith)

/-- **Asymptotic rate is in `(0, 1)`.**
    `0 < (x − 1) / (x · (1 + sstar)²) < 1` for `x > 1`, `sstar > 0`. -/
theorem asymptotic_rate_in_unit (x sstar : ℝ) (hx : x > 1) (hss : sstar > 0) :
    0 < (x - 1) / (x * (1 + sstar) ^ 2) ∧
    (x - 1) / (x * (1 + sstar) ^ 2) < 1 := by
  have hx_pos : x > 0 := by linarith
  have h1s2 : (1 + sstar) ^ 2 > 1 := by nlinarith
  have hD : x * (1 + sstar) ^ 2 > 0 := by positivity
  refine ⟨div_pos (by linarith) hD, ?_⟩
  rw [div_lt_one hD]
  calc x - 1 < x := by linarith
    _ = x * 1 := (mul_one x).symm
    _ < x * (1 + sstar) ^ 2 := by
        exact mul_lt_mul_of_pos_left h1s2 hx_pos

/-! ## §L17. Cube-root shape -/

/-- `S_3(s) = 1 + s + s²`. -/
theorem Sp3_eq (s : ℝ) : Sp 3 s = 1 + s + s ^ 2 := by
  unfold Sp; simp [Finset.sum_range_succ]

/-- `h x 3 s = 1 − (x − 1)/(x(1 + s + s²))`. -/
theorem h_p3 (x s : ℝ) (_hx : x ≠ 0) (_hs : 1 + s + s ^ 2 ≠ 0) :
    pandrosion_h x 3 s = 1 - (x - 1) / (x * (1 + s + s ^ 2)) := by
  unfold pandrosion_h
  rw [Sp3_eq s]

/-- **Cube-root fixed-point characterisation.** -/
theorem cube_root_fixpoint (x : ℝ) (hx : x > 0) (s : ℝ) (hs : 0 ≤ s) :
    pandrosion_h x 3 s = s ↔ s ^ 3 = 1 / x :=
  pandrosion_fixed_point_iff x hx 3 (by omega) s hs

end Pandrosion
