/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP VIII: CONTRACTION for p = 4 (fourth root) and p = 5 (fifth root)

  Same structure as Deep7 (p=3):
  1. Factor Sp(s)-Sp(t) = (s-t)·Qp(s,t)
  2. Bound Qp ≤ Sp(s)·Sp(t)
  3. Conclude |h(s)-s*| ≤ (x-1)/x · |s-s*| < |s-s*|

  Reference: pandrosion_master.tex, Theorem 1729 (extended)
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic
import Pandrosion.Deep

open Finset BigOperators

namespace Pandrosion

/-! ## §68. Contraction for p = 4 (Fourth Root)

Sp₄(s) = 1 + s + s² + s³
Sp₄(s) - Sp₄(t) = (s-t)(1 + s + t + s² + st + t²)
-/

/-- **Sp₄ factoring.** -/
theorem Sp4_sub (s t : ℝ) :
    Sp 4 s - Sp 4 t = (s - t) * (1 + s + t + s ^ 2 + s * t + t ^ 2) := by
  simp [Sp, Finset.sum_range_succ]; ring

/-- **Key bound for p=4: Q₄ ≤ Sp₄(s)·Sp₄(t) for s,t ≥ 0.** -/
theorem Sp4_prod_bound (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    1 + s + t + s ^ 2 + s * t + t ^ 2 ≤
    (1 + s + s ^ 2 + s ^ 3) * (1 + t + t ^ 2 + t ^ 3) := by
  -- RHS - LHS = s³ + t³ + s³t + st³ + s²t + st² + s²t² + s³t² + s²t³ + s³t³ ≥ 0
  nlinarith [sq_nonneg s, sq_nonneg t, sq_nonneg (s*t),
             mul_nonneg hs ht,
             mul_nonneg (mul_nonneg hs hs) hs,
             mul_nonneg (mul_nonneg ht ht) ht,
             mul_nonneg (mul_nonneg hs hs) ht,
             mul_nonneg hs (mul_nonneg ht ht),
             mul_nonneg (mul_nonneg hs hs) (mul_nonneg ht ht)]

/-- **Contraction for p = 4.** -/
theorem contraction_p4 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 4 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (_hs_ne : s ≠ sstar) :
    |pandrosion_h x 4 s - sstar| ≤ (x - 1) / x * |s - sstar| := by
  have hx_pos : x > 0 := by linarith
  have hSs : Sp 4 s > 0 := Sp_pos 4 (by omega) s hs
  have hSt : Sp 4 sstar > 0 := Sp_pos 4 (by omega) sstar (le_of_lt hss_pos)
  have hfix : pandrosion_h x 4 sstar = sstar :=
    (fixed_point_iff x hx_pos 4 (by omega) sstar (le_of_lt hss_pos) (ne_of_lt hss_lt)).mpr hss_eq
  have hD : x * Sp 4 s * Sp 4 sstar > 0 := by positivity
  have hD_ne : x * Sp 4 s * Sp 4 sstar ≠ 0 := ne_of_gt hD
  have hdiff : pandrosion_h x 4 s - pandrosion_h x 4 sstar =
      (x - 1) * (Sp 4 s - Sp 4 sstar) / (x * Sp 4 s * Sp 4 sstar) := by
    unfold pandrosion_h; rw [eq_div_iff hD_ne]; field_simp; ring
  rw [← hfix, hdiff, Sp4_sub]
  -- Factor and take abs
  have hQ : 1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2 > 0 := by positivity
  rw [show (x - 1) * ((s - sstar) * (1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2)) =
      (x - 1) * (1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2) * (s - sstar) from by ring]
  rw [abs_div, abs_mul, abs_mul,
      abs_of_pos (by linarith : x - 1 > 0),
      abs_of_pos hQ, abs_of_pos hD, hfix]
  -- Cross-multiply
  have hbound := Sp4_prod_bound s sstar hs (le_of_lt hss_pos)
  have hxm1 : x - 1 > 0 := by linarith
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, add_zero] at hbound
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, add_zero]
  rw [div_le_iff (by positivity)]
  have hsimp : (x - 1) / x * |s - sstar| *
      (x * (0 + 1 + s ^ 1 + s ^ 2 + s ^ 3) * (0 + 1 + sstar ^ 1 + sstar ^ 2 + sstar ^ 3))
    = (x - 1) * |s - sstar| *
      ((0 + 1 + s ^ 1 + s ^ 2 + s ^ 3) * (0 + 1 + sstar ^ 1 + sstar ^ 2 + sstar ^ 3)) := by
    field_simp; ring
  rw [hsimp]
  have hnn : (x - 1) * |s - sstar| ≥ 0 := mul_nonneg (le_of_lt hxm1) (abs_nonneg _)
  calc (x - 1) * (1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2) * |s - sstar|
      = (x - 1) * |s - sstar| * (1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2) := by ring
    _ ≤ (x - 1) * |s - sstar| *
        ((0 + 1 + s ^ 1 + s ^ 2 + s ^ 3) * (0 + 1 + sstar ^ 1 + sstar ^ 2 + sstar ^ 3)) := by
        apply mul_le_mul_of_nonneg_left _ hnn
        simp only [zero_add, pow_one]; exact hbound
    _ = _ := rfl

/-- **Distance decrease for p = 4.** -/
theorem distance_decreases_p4 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 4 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x 4 s - sstar| < |s - sstar| := by
  have hbound := contraction_p4 x sstar hx hss_pos hss_lt hss_eq s hs hs_ne
  have hfrac : (x - 1) / x < 1 := by rw [div_lt_one (by linarith)]; linarith
  linarith [mul_lt_mul_of_pos_right hfrac (abs_pos.mpr (sub_ne_zero.mpr hs_ne))]

/-! ## §69. Contraction for p = 5 (Fifth Root)

Sp₅(s) = 1 + s + s² + s³ + s⁴
Sp₅(s) - Sp₅(t) = (s-t)(1+s+t+s²+st+t²+s³+s²t+st²+t³)
-/

/-- **Sp₅ factoring.** -/
theorem Sp5_sub (s t : ℝ) :
    Sp 5 s - Sp 5 t = (s - t) *
    (1 + s + t + s^2 + s*t + t^2 + s^3 + s^2*t + s*t^2 + t^3) := by
  simp [Sp, Finset.sum_range_succ]; ring

/-- **Key bound for p=5: Q₅ ≤ Sp₅(s)·Sp₅(t).** -/
theorem Sp5_prod_bound (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    1 + s + t + s^2 + s*t + t^2 + s^3 + s^2*t + s*t^2 + t^3 ≤
    (1 + s + s^2 + s^3 + s^4) * (1 + t + t^2 + t^3 + t^4) := by
  nlinarith [sq_nonneg s, sq_nonneg t, sq_nonneg (s*t),
             mul_nonneg hs ht,
             mul_nonneg (mul_nonneg hs hs) hs,
             mul_nonneg (mul_nonneg ht ht) ht,
             mul_nonneg (mul_nonneg hs hs) ht,
             mul_nonneg hs (mul_nonneg ht ht),
             mul_nonneg (mul_nonneg hs hs) (mul_nonneg ht ht),
             mul_nonneg (mul_nonneg (mul_nonneg hs hs) hs) ht,
             mul_nonneg hs (mul_nonneg (mul_nonneg ht ht) ht),
             mul_nonneg (mul_nonneg (mul_nonneg hs hs) hs) (mul_nonneg ht ht),
             mul_nonneg (mul_nonneg hs hs) (mul_nonneg (mul_nonneg ht ht) ht),
             mul_nonneg (mul_nonneg (mul_nonneg hs hs) hs) (mul_nonneg (mul_nonneg ht ht) ht)]

/-- **Contraction for p = 5.** -/
theorem contraction_p5 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 5 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (_hs_ne : s ≠ sstar) :
    |pandrosion_h x 5 s - sstar| ≤ (x - 1) / x * |s - sstar| := by
  have hx_pos : x > 0 := by linarith
  have hSs : Sp 5 s > 0 := Sp_pos 5 (by omega) s hs
  have hSt : Sp 5 sstar > 0 := Sp_pos 5 (by omega) sstar (le_of_lt hss_pos)
  have hfix : pandrosion_h x 5 sstar = sstar :=
    (fixed_point_iff x hx_pos 5 (by omega) sstar (le_of_lt hss_pos) (ne_of_lt hss_lt)).mpr hss_eq
  have hD : x * Sp 5 s * Sp 5 sstar > 0 := by positivity
  have hD_ne : x * Sp 5 s * Sp 5 sstar ≠ 0 := ne_of_gt hD
  have hdiff : pandrosion_h x 5 s - pandrosion_h x 5 sstar =
      (x - 1) * (Sp 5 s - Sp 5 sstar) / (x * Sp 5 s * Sp 5 sstar) := by
    unfold pandrosion_h; rw [eq_div_iff hD_ne]; field_simp; ring
  rw [← hfix, hdiff, Sp5_sub]
  have hQ : 1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3 > 0 := by positivity
  rw [show (x-1) * ((s-sstar) * (1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3))
      = (x-1) * (1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3) * (s-sstar) from by ring]
  rw [abs_div, abs_mul, abs_mul,
      abs_of_pos (by linarith : x - 1 > 0),
      abs_of_pos hQ, abs_of_pos hD, hfix]
  have hbound := Sp5_prod_bound s sstar hs (le_of_lt hss_pos)
  have hxm1 : x - 1 > 0 := by linarith
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, add_zero] at hbound
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, add_zero]
  rw [div_le_iff (by positivity)]
  have hsimp : (x-1)/x * |s-sstar| *
      (x * (0+1+s^1+s^2+s^3+s^4) * (0+1+sstar^1+sstar^2+sstar^3+sstar^4))
    = (x-1) * |s-sstar| *
      ((0+1+s^1+s^2+s^3+s^4) * (0+1+sstar^1+sstar^2+sstar^3+sstar^4)) := by
    field_simp; ring
  rw [hsimp]
  have hnn : (x-1) * |s-sstar| ≥ 0 := mul_nonneg (le_of_lt hxm1) (abs_nonneg _)
  calc (x-1) * (1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3) * |s-sstar|
      = (x-1) * |s-sstar| * (1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3) := by ring
    _ ≤ (x-1) * |s-sstar| *
        ((0+1+s^1+s^2+s^3+s^4) * (0+1+sstar^1+sstar^2+sstar^3+sstar^4)) := by
        apply mul_le_mul_of_nonneg_left _ hnn
        simp only [zero_add, pow_one]; exact hbound
    _ = _ := rfl

/-- **Distance decrease for p = 5.** -/
theorem distance_decreases_p5 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 5 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x 5 s - sstar| < |s - sstar| := by
  have hbound := contraction_p5 x sstar hx hss_pos hss_lt hss_eq s hs hs_ne
  have hfrac : (x - 1) / x < 1 := by rw [div_lt_one (by linarith)]; linarith
  linarith [mul_lt_mul_of_pos_right hfrac (abs_pos.mpr (sub_ne_zero.mpr hs_ne))]

end Pandrosion
