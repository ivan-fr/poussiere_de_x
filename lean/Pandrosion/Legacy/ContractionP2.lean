/-
  Universitas Pandrosion — Legacy / Contraction (p = 2).

  The square-root case, with the explicit contraction identity
    h(s) − h(t) = (x−1)(s − t) / [x(1+s)(1+t)]
  that drives strict distance-decrease and monotone convergence.

  Ported from `ContractionIdentity.lean` + `MonotoneConvergence.lean`
  on the `main_previous` branch.
-/

import Pandrosion.Legacy.Basic

namespace Pandrosion

open Finset BigOperators

/-! ## §L5. `Sp 2 s = 1 + s`, `h x 2 s` shape -/

/-- `S_2(s) = 1 + s`. -/
theorem Sp2_eq (s : ℝ) : Sp 2 s = 1 + s := by
  unfold Sp; simp [Finset.sum_range_succ]

/-- `h x 2 s = 1 − (x−1)/(x(1+s))`. -/
theorem h_p2 (x s : ℝ) : pandrosion_h x 2 s = 1 - (x - 1) / (x * (1 + s)) := by
  unfold pandrosion_h; rw [Sp2_eq]

/-! ## §L6. The contraction identity (p = 2) -/

/-- **Contraction identity (p = 2).**
    `h(s) − h(t) = (x−1)(s−t) / [x(1+s)(1+t)]`. -/
theorem contraction_identity_p2 (x s t : ℝ) (hs1 : (1 : ℝ) + s ≠ 0)
    (ht1 : (1 : ℝ) + t ≠ 0) (hx : x ≠ 0) :
    pandrosion_h x 2 s - pandrosion_h x 2 t =
    (x - 1) * (s - t) / (x * (1 + s) * (1 + t)) := by
  simp only [pandrosion_h, Sp2_eq]
  have h1 : x * (1 + s) ≠ 0 := mul_ne_zero hx hs1
  have h3 : x * (1 + s) * (1 + t) ≠ 0 := mul_ne_zero h1 ht1
  rw [eq_div_iff h3]
  field_simp
  ring

/-! ## §L7. Strict distance decrease and monotonicity (p = 2) -/

/-- `h(s) − s* = (x−1)(s − s*) / [x(1+s)(1+s*)]`. -/
theorem convergence_p2 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 2 = 1 / x) (s : ℝ) (_hs : s ≥ 0) :
    pandrosion_h x 2 s - sstar =
    (x - 1) * (s - sstar) / (x * (1 + s) * (1 + sstar)) := by
  have h_fix : pandrosion_h x 2 sstar = sstar :=
    (pandrosion_fixed_point_iff x (by linarith) 2 (by omega) sstar
      (le_of_lt hss_pos)).mpr hss_eq
  calc pandrosion_h x 2 s - sstar
      = pandrosion_h x 2 s - pandrosion_h x 2 sstar := by rw [h_fix]
    _ = (x - 1) * (s - sstar) / (x * (1 + s) * (1 + sstar)) :=
        contraction_identity_p2 x s sstar (by linarith) (by linarith) (by linarith)

/-- The contraction factor `(x−1)/[x(1+s)(1+t)] < 1`. -/
theorem contraction_factor_lt_one (x : ℝ) (hx : x > 1)
    (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    (x - 1) / (x * (1 + s) * (1 + t)) < 1 := by
  have hD : x * (1 + s) * (1 + t) > 0 :=
    mul_pos (mul_pos (by linarith) (by linarith)) (by linarith)
  rw [div_lt_one hD]
  have : x * (1 + s) ≥ x := by nlinarith
  nlinarith

/-- The contraction factor is positive. -/
theorem contraction_factor_pos (x : ℝ) (hx : x > 1)
    (s t : ℝ) (_hs : s ≥ 0) (_ht : t ≥ 0) :
    (x - 1) / (x * (1 + s) * (1 + t)) > 0 := by
  apply div_pos (by linarith)
  exact mul_pos (mul_pos (by linarith) (by linarith)) (by linarith)

/-- **Strict distance decrease (p = 2).**
    `|h(s) − s*| < |s − s*|` for `s ≥ 0`, `s ≠ s*`. -/
theorem distance_decreases_p2 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 2 = 1 / x) (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x 2 s - sstar| < |s - sstar| := by
  rw [convergence_p2 x sstar hx hss_pos hss_lt hss_eq s hs]
  set D := x * (1 + s) * (1 + sstar)
  have hD_pos : D > 0 :=
    mul_pos (mul_pos (by linarith) (by linarith)) (by linarith)
  rw [abs_div, abs_mul, abs_of_pos hD_pos]
  rw [div_lt_iff hD_pos]
  have hne : s - sstar ≠ 0 := sub_ne_zero.mpr hs_ne
  rw [mul_comm]
  apply mul_lt_mul_of_pos_left _ (abs_pos.mpr hne)
  rw [abs_of_pos (by linarith : x - 1 > 0)]
  show x - 1 < D
  have h1 : (1 : ℝ) + s ≥ 1 := by linarith
  have h2 : (1 : ℝ) + sstar > 1 := by linarith
  have h3 : x * (1 + s) ≥ x * 1 := mul_le_mul_of_nonneg_left h1 (by linarith)
  calc x - 1 < x := by linarith
    _ = x * 1 := (mul_one x).symm
    _ ≤ x * (1 + s) := h3
    _ = x * (1 + s) * 1 := (mul_one _).symm
    _ < x * (1 + s) * (1 + sstar) := by
        apply mul_lt_mul_of_pos_left h2
        exact mul_pos (by linarith) (by linarith)

/-- `h` is strictly increasing for `p = 2`, `x > 1`. -/
theorem h_strict_mono_p2 (x : ℝ) (hx : x > 1) (s t : ℝ)
    (hs : s ≥ 0) (_ht : t ≥ 0) (hst : s < t) :
    pandrosion_h x 2 s < pandrosion_h x 2 t := by
  simp only [pandrosion_h, Sp2_eq]
  have hxm1 : x - 1 > 0 := by linarith
  have hs1 : x * (1 + s) > 0 := by positivity
  have hlt : x * (1 + s) < x * (1 + t) := by nlinarith
  linarith [div_lt_div_of_pos_left hxm1 hs1 hlt]

/-- **Monotone from above (p = 2).** `sstar < s ≤ 1` ⇒ `sstar < h(s) < s`. -/
theorem monotone_from_above_p2 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 2 = 1 / x)
    (s : ℝ) (hs : s > sstar) (hs1 : s ≤ 1) :
    sstar < pandrosion_h x 2 s ∧ pandrosion_h x 2 s < s := by
  have hfix : pandrosion_h x 2 sstar = sstar :=
    (pandrosion_fixed_point_iff x (by linarith) 2 (by omega) sstar
      (le_of_lt hss_pos)).mpr hss_eq
  refine ⟨?_, ?_⟩
  · rw [← hfix]
    exact h_strict_mono_p2 x hx sstar s (le_of_lt hss_pos)
      (le_of_lt (lt_trans hss_pos hs)) hs
  · have hdist := distance_decreases_p2 x sstar hx hss_pos hss_lt hss_eq s
      (le_of_lt (lt_trans hss_pos hs)) (ne_of_gt hs)
    have hhs_gt : pandrosion_h x 2 s > sstar := by
      rw [← hfix]
      exact h_strict_mono_p2 x hx sstar s (le_of_lt hss_pos)
        (le_of_lt (lt_trans hss_pos hs)) hs
    rw [abs_of_pos (by linarith : s - sstar > 0)] at hdist
    rw [abs_of_pos (by linarith : pandrosion_h x 2 s - sstar > 0)] at hdist
    linarith

/-- **Monotone from below (p = 2).** `0 < s < sstar` ⇒ `s < h(s) < sstar`. -/
theorem monotone_from_below_p2 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 2 = 1 / x)
    (s : ℝ) (hs_pos : s > 0) (hs : s < sstar) :
    s < pandrosion_h x 2 s ∧ pandrosion_h x 2 s < sstar := by
  have hfix : pandrosion_h x 2 sstar = sstar :=
    (pandrosion_fixed_point_iff x (by linarith) 2 (by omega) sstar
      (le_of_lt hss_pos)).mpr hss_eq
  refine ⟨?_, ?_⟩
  · have hdist := distance_decreases_p2 x sstar hx hss_pos hss_lt hss_eq s
      (le_of_lt hs_pos) (ne_of_lt hs)
    have hhs_lt : pandrosion_h x 2 s < sstar := by
      rw [← hfix]
      exact h_strict_mono_p2 x hx s sstar (le_of_lt hs_pos)
        (le_of_lt hss_pos) hs
    rw [abs_of_neg (by linarith : s - sstar < 0)] at hdist
    rw [abs_of_neg (by linarith : pandrosion_h x 2 s - sstar < 0)] at hdist
    linarith
  · rw [← hfix]
    exact h_strict_mono_p2 x hx s sstar (le_of_lt hs_pos) (le_of_lt hss_pos) hs

/-- **Uniqueness of the positive fixed point (p = 2).**
    Two positive solutions of `s^2 = 1/x` must coincide. -/
theorem fixed_point_unique_p2 (x : ℝ) (_hx : x > 1)
    (s₁ s₂ : ℝ) (hs1_pos : s₁ > 0) (_hs1_lt : s₁ < 1)
    (hs2_pos : s₂ > 0) (_hs2_lt : s₂ < 1)
    (hs1_eq : s₁ ^ 2 = 1 / x) (hs2_eq : s₂ ^ 2 = 1 / x) :
    s₁ = s₂ := by
  have h2 : (s₁ - s₂) * (s₁ + s₂) = 0 := by nlinarith
  rcases mul_eq_zero.mp h2 with h3 | h3
  · linarith
  · linarith

end Pandrosion
