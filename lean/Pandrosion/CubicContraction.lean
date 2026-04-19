/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP VII: CONTRACTION for p = 3 (cube root case)

  Key algebraic identity:
    Sp₃(s) - Sp₃(t) = (s-t)(1+s+t)

  Key bound for contraction:
    (1+s+t) ≤ (1+s+s²)(1+t+t²)  for s,t ≥ 0

  Therefore the contraction ratio is at most (x-1)/x < 1.

  Reference: pandrosion_master.tex, Theorem 1729
-/
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.GeomSum
import Mathlib.Tactic
import Pandrosion.FixedPointGeneral

open Finset BigOperators

namespace Pandrosion

/-! ## §65. Explicit Factoring for p = 3

For p = 3: Sp(s) = 1 + s + s²
Sp(s) - Sp(t) = (s-t) + (s²-t²) = (s-t)(1+s+t)
-/

/-- **Sp₃ factoring: Sp₃(s) - Sp₃(t) = (s-t)(1+s+t).** -/
theorem Sp3_sub (s t : ℝ) :
    Sp 3 s - Sp 3 t = (s - t) * (1 + s + t) := by
  simp [Sp, Finset.sum_range_succ]; ring

/-- **Key bound: (1+s+t) ≤ (1+s+s²)(1+t+t²) for s,t ≥ 0.**
    Proof: expanding, the RHS minus LHS equals
    t² + st + st² + s² + s²t + s²t² ≥ 0. -/
theorem Sp3_prod_bound (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    1 + s + t ≤ (1 + s + s ^ 2) * (1 + t + t ^ 2) := by
  nlinarith [sq_nonneg s, sq_nonneg t, sq_nonneg (s * t),
             mul_nonneg hs ht, mul_nonneg (mul_nonneg hs hs) ht,
             mul_nonneg hs (mul_nonneg ht ht)]

/-- **THE CONTRACTION FOR p = 3.**
    |h(s) - s*| ≤ (x-1)/x · |s - s*|.

    Proof outline:
    1. h(s) - s* = (x-1)(Sp(s)-Sp(s*))/(x·Sp(s)·Sp(s*))
    2. Sp(s)-Sp(s*) = (s-s*)(1+s+s*)
    3. (1+s+s*) ≤ Sp(s)·Sp(s*)  [the key bound]
    4. Therefore |h(s)-s*| ≤ (x-1)/x · |s-s*| -/
theorem contraction_p3 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 3 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (_hs_ne : s ≠ sstar) :
    |pandrosion_h x 3 s - sstar| ≤ (x - 1) / x * |s - sstar| := by
  have hx_pos : x > 0 := by linarith
  have hSs : Sp 3 s > 0 := Sp_pos 3 (by omega) s hs
  have hSt : Sp 3 sstar > 0 := Sp_pos 3 (by omega) sstar (le_of_lt hss_pos)
  -- h(s*) = s*
  have hfix : pandrosion_h x 3 sstar = sstar :=
    (fixed_point_iff x hx_pos 3 (by omega) sstar (le_of_lt hss_pos) (ne_of_lt hss_lt)).mpr hss_eq
  -- h(s) - s* = h(s) - h(s*) = (x-1)(Sp s - Sp s*)/(x·Sp s·Sp s*)
  -- Using the definition directly:
  -- h(s) = 1 - (x-1)/(x·Sp s)
  -- h(s) - s* = 1-(x-1)/(x·Sp s) - s* = (1-s*) - (x-1)/(x·Sp s)
  -- Instead, use algebra on the h_diff directly.
  -- pandrosion_h x 3 s - pandrosion_h x 3 sstar
  --   = (x-1)/x · [1/Sp(sstar) - 1/Sp(s)]
  --   = (x-1) · [Sp(s) - Sp(sstar)] / (x · Sp(s) · Sp(sstar))
  --   = (x-1)(s-sstar)(1+s+sstar) / (x · Sp(s) · Sp(sstar))
  have hD : x * Sp 3 s * Sp 3 sstar > 0 := by positivity
  have hD_ne : x * Sp 3 s * Sp 3 sstar ≠ 0 := ne_of_gt hD
  -- Compute the difference via the definition
  have hdiff : pandrosion_h x 3 s - pandrosion_h x 3 sstar =
      (x - 1) * (Sp 3 s - Sp 3 sstar) / (x * Sp 3 s * Sp 3 sstar) := by
    unfold pandrosion_h
    rw [eq_div_iff hD_ne]; field_simp; ring
  rw [← hfix, hdiff, Sp3_sub]
  -- Now we have: |(x-1)·(s-s*)·(1+s+s*) / (x·Sp s·Sp s*)| ≤ (x-1)/x · |s-s*|
  rw [mul_comm (s - sstar) (1 + s + sstar)]
  rw [show (x - 1) * ((1 + s + sstar) * (s - sstar)) = (x - 1) * (1 + s + sstar) * (s - sstar) from by ring]
  rw [abs_div, abs_mul, abs_mul,
      abs_of_pos (by linarith : x - 1 > 0),
      abs_of_pos (by linarith : 1 + s + sstar > 0),
      abs_of_pos hD]
  -- Goal: (x-1) * (1+s+s*) * |s-s*| / (x * Sp s * Sp s*) ≤ (x-1)/x * |s-s*|
  rw [hfix]
  -- Reduce to showing (1+s+sstar) ≤ Sp 3 s * Sp 3 sstar
  have hbound := Sp3_prod_bound s sstar hs (le_of_lt hss_pos)
  have hxm1 : x - 1 > 0 := by linarith
  -- habs_nn used implicitly
  -- Factor: (x-1)*(1+s+s*)*|s-s*|/(x*Sp*Sp*) ≤ (x-1)*|s-s*|/x
  -- ⟺ (1+s+s*)/(Sp*Sp*) ≤ 1
  -- ⟺ (1+s+s*) ≤ Sp*Sp*
  -- Which is exactly hbound (after unfolding Sp 3)
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, add_zero] at hbound
  -- Now hbound : 1+s+sstar ≤ (1+s+s^2)*(1+sstar+sstar^2)
  -- And we need to show the div inequality
  -- (x-1)*(1+s+sstar)*|s-sstar| / (x*(1+s+s^2)*(1+sstar+sstar^2)) ≤ (x-1)/x*|s-sstar|
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, add_zero]
  -- Cross-multiply everything
  -- hprod, hprod2 used via positivity in the div_le_iff call
  -- (x-1)*(1+s+s*)*|s-s*| / (x*(Sp s)*(Sp s*)) ≤ (x-1)/x * |s-s*|
  -- ⟺ (x-1)*(1+s+s*)*|s-s*| ≤ (x-1)/x * |s-s*| * x * (Sp s) * (Sp s*)
  -- ⟺ (1+s+s*) ≤ (Sp s) * (Sp s*)  [cancel (x-1)*|s-s*|*... ]
  -- Use: a/b ≤ c ⟺ a ≤ c * b for b > 0
  rw [div_le_iff (by positivity : x * (0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2) > 0)]
  -- Simplify (x-1)/x * |s-s*| * (x * Sp * Sp*)
  -- = (x-1) * |s-s*| * Sp * Sp*  [since (x-1)/x * x = x-1]
  have hsimp : (x - 1) / x * |s - sstar| * (x * (0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2))
    = (x - 1) * |s - sstar| * ((0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2)) := by
    field_simp; ring
  rw [hsimp]
  -- Now: (x-1)*(1+s+s*)*|s-s*| ≤ (x-1)*|s-s*|*((1+s+s²)*(1+s*+s*²))
  -- This is (x-1)*|s-s*| · (1+s+s*) ≤ (x-1)*|s-s*| · (Sp*Sp*)
  -- Follows from hbound by multiplying by (x-1)*|s-s*| ≥ 0
  have hnn : (x - 1) * |s - sstar| ≥ 0 := mul_nonneg (le_of_lt hxm1) (abs_nonneg _)
  calc (x - 1) * (1 + s + sstar) * |s - sstar|
      = (x - 1) * |s - sstar| * (1 + s + sstar) := by ring
    _ ≤ (x - 1) * |s - sstar| * ((0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2)) := by
        apply mul_le_mul_of_nonneg_left _ hnn
        simp only [zero_add, pow_one]; exact hbound
    _ = (x - 1) * |s - sstar| * ((0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2)) := rfl

/-- **Corollary: strict distance decrease for p = 3.** -/
theorem distance_decreases_p3 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 3 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x 3 s - sstar| < |s - sstar| := by
  have hbound := contraction_p3 x sstar hx hss_pos hss_lt hss_eq s hs hs_ne
  have hfrac : (x - 1) / x < 1 := by rw [div_lt_one (by linarith)]; linarith
  have habs_pos : |s - sstar| > 0 := abs_pos.mpr (sub_ne_zero.mpr hs_ne)
  linarith [mul_lt_mul_of_pos_right hfrac habs_pos]

/-- **The universal factoring identity: s^n - t^n = (s-t) · ∑ s^j · t^{n-1-j}.**
    Direct from Mathlib's `geom_sum₂_mul`. -/
theorem pow_sub_factor (s t : ℝ) (n : ℕ) :
    s ^ n - t ^ n = (s - t) * ∑ j in range n, s ^ j * t ^ (n - 1 - j) := by
  have := geom_sum₂_mul s t n
  linarith [mul_comm (∑ i in range n, s ^ i * t ^ (n - 1 - i)) (s - t)]

/-- **The contraction ratio (x-1)/x < 1 for x > 1.** -/
theorem contraction_ratio_lt_one (x : ℝ) (hx : x > 1) :
    (x - 1) / x < 1 := by
  rw [div_lt_one (by linarith)]; linarith

end Pandrosion
