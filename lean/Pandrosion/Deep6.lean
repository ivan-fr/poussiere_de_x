/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP VI: Monotone convergence, Banach contraction, uniqueness

  The Pandrosion iteration h is monotone increasing (h' > 0),
  so the sequence h^n(s₀) converges monotonically to s*.

  Reference: pandrosion_master.tex, Theorems 336, 405, 670
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic
import Pandrosion.Deep
import Pandrosion.Deep2

open Finset BigOperators

namespace Pandrosion

/-! ## §62. h is Monotone Increasing (p = 2)

Since h'(s) = (x-1)/(x(1+s)²) > 0, the iteration h is
strictly monotone increasing. Formally:
  s < t ⟹ h(s) < h(t).
-/

/-- **h is strictly increasing for p = 2, x > 1.**
    h(s) < h(t) when s < t. -/
theorem h_strict_mono_p2 (x : ℝ) (hx : x > 1) (s t : ℝ) (hs : s ≥ 0)
    (_ht : t ≥ 0) (hst : s < t) :
    pandrosion_h x 2 s < pandrosion_h x 2 t := by
  simp only [pandrosion_h, Sp2_eq]
  -- h(s) = 1 - (x-1)/(x(1+s)), h(t) = 1 - (x-1)/(x(1+t))
  -- Since x > 1: (x-1) > 0
  -- Since 1+s < 1+t: x(1+s) < x(1+t)
  -- So (x-1)/(x(1+s)) > (x-1)/(x(1+t))
  -- So 1 - (x-1)/(x(1+s)) < 1 - (x-1)/(x(1+t))
  have hxm1 : x - 1 > 0 := by linarith
  have hs1 : x * (1 + s) > 0 := by positivity
  have ht1 : x * (1 + t) > 0 := by positivity
  have hlt : x * (1 + s) < x * (1 + t) := by nlinarith
  linarith [div_lt_div_of_pos_left hxm1 hs1 hlt]

/-! ## §63. Monotone Convergence

If h is increasing and |h(s)-s*| < |s-s*|,
then the sequence h^n(s₀) converges monotonically.
-/

/-- **Key: h(s*) = s* and h increasing ⟹ s > s* ⟹ h(s) > s*.**
    Combined with contraction: s > s* ⟹ s* < h(s) < s.
    So the sequence decreases monotonically to s*. -/
theorem monotone_from_above_p2 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 2 = 1 / x)
    (s : ℝ) (hs : s > sstar) (hs1 : s ≤ 1) :
    sstar < pandrosion_h x 2 s ∧ pandrosion_h x 2 s < s := by
  constructor
  · -- h(s) > h(s*) = s* since h is increasing and s > s*
    have hfix : pandrosion_h x 2 sstar = sstar :=
      (fixed_point_iff x (by linarith) 2 (by omega) sstar (le_of_lt hss_pos)
        (ne_of_lt hss_lt)).mpr hss_eq
    rw [← hfix]
    exact h_strict_mono_p2 x (by linarith) sstar s (le_of_lt hss_pos)
      (le_of_lt (lt_trans hss_pos hs)) hs
  · -- h(s) < s: since |h(s) - s*| < |s - s*| and h(s) > s*, we get h(s) < s
    -- We use the contraction: h(s) - s* = λ(s - s*) with 0 < λ < 1
    -- So h(s) = s* + λ(s-s*) < s* + (s-s*) = s
    have hdist := distance_decreases_p2 x sstar hx hss_pos hss_lt hss_eq s
      (le_of_lt (lt_trans hss_pos hs)) (ne_of_gt hs)
    have hfix : pandrosion_h x 2 sstar = sstar :=
      (fixed_point_iff x (by linarith) 2 (by omega) sstar (le_of_lt hss_pos)
        (ne_of_lt hss_lt)).mpr hss_eq
    rw [← hfix] at hs
    have hhs_gt : pandrosion_h x 2 s > sstar := by
      rw [← hfix]
      exact h_strict_mono_p2 x (by linarith) sstar s (le_of_lt hss_pos)
        (le_of_lt (lt_trans hss_pos (by linarith [hfix]))) (by linarith [hfix])
    -- |h(s) - s*| < |s - s*|, both sides positive since h(s) > s* and s > s*
    rw [abs_of_pos (by linarith : s - sstar > 0)] at hdist
    -- hdist: |h(s) - s*| < s - s*
    -- Since h(s) > s*: h(s) - s* ≥ 0, so |h(s) - s*| = h(s) - s*
    rw [abs_of_pos (by linarith : pandrosion_h x 2 s - sstar > 0)] at hdist
    linarith

/-- **Monotone from below: s < s* ⟹ s < h(s) < s*.**
    Since h is increasing and a contraction. -/
theorem monotone_from_below_p2 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 2 = 1 / x)
    (s : ℝ) (hs_pos : s > 0) (hs : s < sstar) :
    s < pandrosion_h x 2 s ∧ pandrosion_h x 2 s < sstar := by
  constructor
  · -- h(s) > h(s*) is wrong here. We need h(s) > s.
    -- Use contraction: |h(s) - s*| < |s - s*|
    -- s < s* and h(s) < s* (from contraction), so s* - h(s) < s* - s, i.e., h(s) > s
    have hdist := distance_decreases_p2 x sstar hx hss_pos hss_lt hss_eq s
      (le_of_lt hs_pos) (ne_of_lt hs)
    have hfix : pandrosion_h x 2 sstar = sstar :=
      (fixed_point_iff x (by linarith) 2 (by omega) sstar (le_of_lt hss_pos)
        (ne_of_lt hss_lt)).mpr hss_eq
    have hhs_lt : pandrosion_h x 2 s < sstar := by
      rw [← hfix]
      exact h_strict_mono_p2 x (by linarith) s sstar (le_of_lt hs_pos)
        (le_of_lt hss_pos) hs
    -- |h(s) - s*| < |s - s*|, both negative since h(s) < s* and s < s*
    rw [abs_of_neg (by linarith : s - sstar < 0)] at hdist
    rw [abs_of_neg (by linarith : pandrosion_h x 2 s - sstar < 0)] at hdist
    linarith
  · -- h(s) < h(s*) = s* since h increasing and s < s*
    have hfix : pandrosion_h x 2 sstar = sstar :=
      (fixed_point_iff x (by linarith) 2 (by omega) sstar (le_of_lt hss_pos)
        (ne_of_lt hss_lt)).mpr hss_eq
    rw [← hfix]
    exact h_strict_mono_p2 x (by linarith) s sstar (le_of_lt hs_pos) (le_of_lt hss_pos) hs

/-! ## §64. Banach Fixed-Point Theorem Application

Since h is a contraction on (0,1) with Lipschitz constant < 1,
Banach's theorem guarantees:
1. Uniqueness of the fixed point
2. Convergence from any starting point in (0,1)
-/

/-- **Uniqueness: if s₁ and s₂ are both fixed points in (0,1),
    then s₁ = s₂.** Proof: if s₁ ≠ s₂, contraction gives
    |s₁ - s₂| = |h(s₁) - h(s₂)| < |s₁ - s₂|, contradiction. -/
theorem fixed_point_unique_p2 (_x : ℝ) (_hx : _x > 1)
    (s₁ s₂ : ℝ) (hs1_pos : s₁ > 0) (_hs1_lt : s₁ < 1)
    (hs2_pos : s₂ > 0) (_hs2_lt : s₂ < 1)
    (hs1_eq : s₁ ^ 2 = 1 / _x) (hs2_eq : s₂ ^ 2 = 1 / _x) :
    s₁ = s₂ := by
  -- s₁^2 = s₂^2 = 1/x, and both positive, so s₁ = s₂
  have _h1 : s₁ ^ 2 - s₂ ^ 2 = 0 := by linarith
  have h2 : (s₁ - s₂) * (s₁ + s₂) = 0 := by nlinarith
  rcases mul_eq_zero.mp h2 with h3 | h3
  · linarith
  · linarith

end Pandrosion
