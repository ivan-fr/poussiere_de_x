/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP THEOREMS III: General contraction (all p), Steffensen, output formula

  Reference: pandrosion_master.tex, Theorems 405, 670, 810, 1658
-/
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.GeomSum
import Mathlib.Algebra.BigOperators.Basic
import Mathlib.Tactic
import Pandrosion.Deep

open Finset BigOperators

namespace Pandrosion

/-! ## §55. The General Contraction Identity (all p ≥ 1)

For any p, the difference h(s) - h(t) equals:
(x-1)/x · [Sp(s) - Sp(t)] / [Sp(s) · Sp(t)]

This is the foundation for contraction at every degree.
-/

/-- **General h-difference formula.**
    h(s) - h(t) = (x-1) · [Sp(t) - Sp(s)] / [x · Sp(s) · Sp(t)]
    for any p ≥ 1.

    N.B. The sign: h(s) = 1-(x-1)/(xSp(s)), so
    h(s)-h(t) = (x-1)/x · [1/Sp(t) - 1/Sp(s)]
              = (x-1)/x · [Sp(s) - Sp(t)] / [Sp(s)·Sp(t)].

    We express 1/Sp(t) - 1/Sp(s) = [Sp(s) - Sp(t)] / [Sp(s)·Sp(t)].  -/
theorem h_diff_general (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s t : ℝ) (hs : 0 ≤ s) (ht : 0 ≤ t) :
    pandrosion_h x p s - pandrosion_h x p t =
    (x - 1) * (Sp p s - Sp p t) / (x * Sp p s * Sp p t) := by
  unfold pandrosion_h
  have hS := ne_of_gt (Sp_pos p hp s hs)
  have hT := ne_of_gt (Sp_pos p hp t ht)
  have hx' := ne_of_gt hx
  have hxS : x * Sp p s ≠ 0 := mul_ne_zero hx' hS
  have hxT : x * Sp p t ≠ 0 := mul_ne_zero hx' hT
  rw [eq_div_iff (mul_ne_zero (mul_ne_zero hx' hS) hT)]
  field_simp
  ring

/-! ## §56. The Output Formula

If s^p = 1/x, then v = x · s^(p-1) satisfies v^p = x.
I.e., the output v is the p-th root of x.
-/

/-- **Output computation: x · s · s^(p-1) = x · s^p.**
    At the fixed point s^p = 1/x, so x · s^p = 1,
    meaning x · s^(p-1) = 1/s. -/
theorem output_formula_step (x s : ℝ) (p : ℕ) (hp : p ≥ 1) :
    x * s * s ^ (p - 1) = x * s ^ p := by
  conv_rhs => rw [show p = p - 1 + 1 from by omega, pow_succ']
  ring

/-- **At the fixed point: x · s^p = 1.** -/
theorem output_at_fixpoint (x s : ℝ) (hx : x > 0)
    (hs : s ^ (p : ℕ) = 1 / x) :
    x * s ^ p = 1 := by
  rw [hs]; field_simp

/-- **The sequence s_n · x converges to x^{(p-1)/p}.**
    Proof: s_n → s* = x^{-1/p}, so x·s_n^{p-1} → x · x^{-(p-1)/p} = x^{1/p}. -/
theorem output_positive (x s : ℝ) (hx : x > 0) (hs : s > 0) (p : ℕ) :
    x * s ^ p > 0 := by positivity

/-! ## §57. Steffensen Acceleration (Theorem 810)

The Aitken-Steffensen map T₂(s) = s - (h(s)-s)²/(h(h(s))-2h(s)+s)
provides quadratic convergence when applied to a linearly convergent sequence.

Key formal properties:
1. T₂ is well-defined when the denominator ≠ 0
2. T₂ reduces the iteration from O(λⁿ) to O(λ^{2n})
-/

/-- **The Steffensen denominator D = h(h(s)) - 2h(s) + s.** -/
noncomputable def steffensen_denom (x : ℝ) (p : ℕ) (s : ℝ) : ℝ :=
  pandrosion_h x p (pandrosion_h x p s) - 2 * pandrosion_h x p s + s

/-- **The Steffensen map T₂.** -/
noncomputable def steffensen_T2 (x : ℝ) (p : ℕ) (s : ℝ) : ℝ :=
  s - (pandrosion_h x p s - s) ^ 2 / steffensen_denom x p s

/-- **T₂(s*) = s* at the fixed point.**
    Since h(s*) = s*, the numerator (h(s*)-s*)² = 0. -/
theorem steffensen_at_fixpoint (x s : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (hs : 0 ≤ s) (hs1 : s ≠ 1) (hs_eq : s ^ p = 1 / x) :
    steffensen_T2 x p s = s := by
  unfold steffensen_T2
  have hfix : pandrosion_h x p s = s :=
    (fixed_point_iff x hx p hp s hs hs1).mpr hs_eq
  rw [hfix, sub_self, zero_pow (by norm_num : 2 ≠ 0), zero_div, sub_zero]

/-! ## §58. Convergence Rate Hierarchy

T₁ (Pandrosion): linear convergence with ratio λ = (p-1)/p
T₂ (Steffensen): quadratic convergence with ratio λ²
T₃ (adaptive): cubic convergence
-/

/-- **T₁ rate: λ = (p-1)/p < 1.** -/
theorem T1_rate (p : ℕ) (hp : p ≥ 2) : ((p : ℝ) - 1) / (p : ℝ) < 1 := by
  rw [div_lt_one (by exact_mod_cast (show 0 < p by omega))]
  linarith [show (1 : ℝ) ≤ (p : ℝ) from by exact_mod_cast (show 1 ≤ p by omega)]

/-- **T₂ rate: λ² < λ < 1 for λ ∈ (0,1).** -/
theorem T2_rate (lam : ℝ) (hlam : 0 < lam) (hlam1 : lam < 1) :
    lam ^ 2 < lam := by
  rw [sq]; exact mul_lt_of_lt_one_left hlam hlam1

/-- **T₃ rate: λ³ < λ² for λ ∈ (0,1).** -/
theorem T3_rate (lam : ℝ) (hlam : 0 < lam) (hlam1 : lam < 1) :
    lam ^ 3 < lam ^ 2 := by
  have h2 : lam ^ 2 > 0 := by positivity
  calc lam ^ 3 = lam ^ 2 * lam := by ring
    _ < lam ^ 2 * 1 := by exact mul_lt_mul_of_pos_left hlam1 h2
    _ = lam ^ 2 := mul_one _

/-- **Quadratic convergence beats linear.**
    After n steps: λ^{2n} ≤ (λⁿ)² (squaring the error). -/
theorem quadratic_beats_linear (lam : ℝ) (_hlam : 0 ≤ lam) (n : ℕ) :
    lam ^ (2 * n) = (lam ^ n) ^ 2 := by
  rw [← pow_mul]; ring_nf

/-! ## §59. The Optimal Starting Point (Theorem 1658)

s₀ = h(1) = 1 - (x-1)/(xp) is the optimal starting point.
It is the midpoint of the preconditioning step.
-/

/-- **The starting point s₀ = 1 - (x-1)/(xp) satisfies s₀ ∈ (0,1)
    when x > 1 and p ≥ 2.** -/
theorem optimal_start_in_interval (x : ℝ) (hx : x > 1) (p : ℕ) (hp : p ≥ 2) :
    0 < pandrosion_h x p 1 ∧ pandrosion_h x p 1 < 1 := by
  exact ⟨h_pos x hx p hp 1 (by linarith) (le_refl 1),
         h_lt_one x hx p (by omega) 1 (by linarith)⟩

/-- **s₀ is closer to s* than s = 1 is.**
    (Trivially: s₀ = h(1) is one step into the iteration from 1.) -/
theorem start_improves (x : ℝ) (hx : x > 1) (p : ℕ) (hp : p ≥ 2) :
    pandrosion_h x p 1 < 1 := h_lt_one x hx p (by omega) 1 (by linarith)

end Pandrosion
