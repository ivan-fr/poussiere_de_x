/-
  Universitas Pandrosion — Legacy / Universal contraction.

  The uniform bound `|h(s) − s*| ≤ (x−1)/x · |s − s*|` for every
  `p ≥ 2` and every `x > 1`. Three layers:

    • concrete factorisations for `p = 3, 4, 5`
      (`Sp3_sub`, `Sp4_sub`, `Sp5_sub` and the accompanying
      product bounds),
    • the universal divided difference `Qp` with the general
      factoring identity `Sp(s) − Sp(t) = (s − t)·Qp(s,t)`,
    • the general contraction theorem `contraction_general`.

  Ported from `CubicContraction.lean`, `FourthFifthRoot.lean`,
  and `GeneralContractionAll.lean` on the `main_previous` branch.
-/

import Pandrosion.Legacy.Basic
import Mathlib.Algebra.GeomSum

namespace Pandrosion

open Finset BigOperators

/-! ## §L8. The universal factoring identity `s^n − t^n = (s−t)·Σ s^j t^{n-1-j}` -/

/-- `s^n − t^n = (s − t) · Σ_{j<n} s^j · t^{n-1-j}`.
    Direct from Mathlib's `geom_sum₂_mul`. -/
theorem pow_sub_factor (s t : ℝ) (n : ℕ) :
    s ^ n - t ^ n = (s - t) * ∑ j in range n, s ^ j * t ^ (n - 1 - j) := by
  have := geom_sum₂_mul s t n
  linarith [mul_comm (∑ i in range n, s ^ i * t ^ (n - 1 - i)) (s - t)]

/-- `(x−1)/x < 1` for `x > 1`. -/
theorem contraction_ratio_lt_one (x : ℝ) (hx : x > 1) :
    (x - 1) / x < 1 := by
  rw [div_lt_one (by linarith)]; linarith

/-! ## §L9. Contraction for `p = 3` -/

/-- `S_3(s) − S_3(t) = (s − t)(1 + s + t)`. -/
theorem Sp3_sub (s t : ℝ) :
    Sp 3 s - Sp 3 t = (s - t) * (1 + s + t) := by
  simp [Sp, Finset.sum_range_succ]; ring

/-- `1 + s + t ≤ (1 + s + s²)(1 + t + t²)` for `s, t ≥ 0`. -/
theorem Sp3_prod_bound (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    1 + s + t ≤ (1 + s + s ^ 2) * (1 + t + t ^ 2) := by
  nlinarith [sq_nonneg s, sq_nonneg t, sq_nonneg (s * t),
             mul_nonneg hs ht, mul_nonneg (mul_nonneg hs hs) ht,
             mul_nonneg hs (mul_nonneg ht ht)]

/-- Contraction `|h(s) − s*| ≤ (x−1)/x · |s − s*|` for `p = 3`. -/
theorem contraction_p3 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 3 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (_hs_ne : s ≠ sstar) :
    |pandrosion_h x 3 s - sstar| ≤ (x - 1) / x * |s - sstar| := by
  have hx_pos : x > 0 := by linarith
  have hSs : Sp 3 s > 0 := Sp_pos 3 (by omega) s hs
  have hSt : Sp 3 sstar > 0 := Sp_pos 3 (by omega) sstar (le_of_lt hss_pos)
  have hfix : pandrosion_h x 3 sstar = sstar :=
    (pandrosion_fixed_point_iff x hx_pos 3 (by omega) sstar
      (le_of_lt hss_pos)).mpr hss_eq
  have hD : x * Sp 3 s * Sp 3 sstar > 0 := by positivity
  have hD_ne : x * Sp 3 s * Sp 3 sstar ≠ 0 := ne_of_gt hD
  have hdiff : pandrosion_h x 3 s - pandrosion_h x 3 sstar =
      (x - 1) * (Sp 3 s - Sp 3 sstar) / (x * Sp 3 s * Sp 3 sstar) := by
    unfold pandrosion_h
    rw [eq_div_iff hD_ne]; field_simp; ring
  rw [← hfix, hdiff, Sp3_sub]
  rw [mul_comm (s - sstar) (1 + s + sstar)]
  rw [show (x - 1) * ((1 + s + sstar) * (s - sstar))
      = (x - 1) * (1 + s + sstar) * (s - sstar) from by ring]
  rw [abs_div, abs_mul, abs_mul,
      abs_of_pos (by linarith : x - 1 > 0),
      abs_of_pos (by linarith : 1 + s + sstar > 0),
      abs_of_pos hD]
  rw [hfix]
  have hbound := Sp3_prod_bound s sstar hs (le_of_lt hss_pos)
  have hxm1 : x - 1 > 0 := by linarith
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero,
             pow_zero, add_zero] at hbound
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero,
             pow_zero, add_zero]
  rw [div_le_iff (by positivity : x * (0 + 1 + s ^ 1 + s ^ 2) *
    (0 + 1 + sstar ^ 1 + sstar ^ 2) > 0)]
  have hsimp : (x - 1) / x * |s - sstar| *
      (x * (0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2))
    = (x - 1) * |s - sstar| *
      ((0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2)) := by
    field_simp; ring
  rw [hsimp]
  have hnn : (x - 1) * |s - sstar| ≥ 0 :=
    mul_nonneg (le_of_lt hxm1) (abs_nonneg _)
  calc (x - 1) * (1 + s + sstar) * |s - sstar|
      = (x - 1) * |s - sstar| * (1 + s + sstar) := by ring
    _ ≤ (x - 1) * |s - sstar| *
        ((0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2)) := by
        apply mul_le_mul_of_nonneg_left _ hnn
        simp only [zero_add, pow_one]; exact hbound
    _ = (x - 1) * |s - sstar| *
        ((0 + 1 + s ^ 1 + s ^ 2) * (0 + 1 + sstar ^ 1 + sstar ^ 2)) := rfl

/-- Strict distance decrease for `p = 3`. -/
theorem distance_decreases_p3 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 3 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x 3 s - sstar| < |s - sstar| := by
  have hbound := contraction_p3 x sstar hx hss_pos hss_lt hss_eq s hs hs_ne
  have hfrac : (x - 1) / x < 1 := contraction_ratio_lt_one x hx
  linarith [mul_lt_mul_of_pos_right hfrac
    (abs_pos.mpr (sub_ne_zero.mpr hs_ne))]

/-! ## §L10. Contraction for `p = 4` -/

/-- `S_4(s) − S_4(t) = (s − t)(1 + s + t + s² + st + t²)`. -/
theorem Sp4_sub (s t : ℝ) :
    Sp 4 s - Sp 4 t =
    (s - t) * (1 + s + t + s ^ 2 + s * t + t ^ 2) := by
  simp [Sp, Finset.sum_range_succ]; ring

/-- `Q₄(s,t) ≤ S_4(s) · S_4(t)` for `s, t ≥ 0`. -/
theorem Sp4_prod_bound (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    1 + s + t + s ^ 2 + s * t + t ^ 2 ≤
    (1 + s + s ^ 2 + s ^ 3) * (1 + t + t ^ 2 + t ^ 3) := by
  nlinarith [sq_nonneg s, sq_nonneg t, sq_nonneg (s*t),
             mul_nonneg hs ht,
             mul_nonneg (mul_nonneg hs hs) hs,
             mul_nonneg (mul_nonneg ht ht) ht,
             mul_nonneg (mul_nonneg hs hs) ht,
             mul_nonneg hs (mul_nonneg ht ht),
             mul_nonneg (mul_nonneg hs hs) (mul_nonneg ht ht)]

/-- Contraction for `p = 4`. -/
theorem contraction_p4 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 4 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (_hs_ne : s ≠ sstar) :
    |pandrosion_h x 4 s - sstar| ≤ (x - 1) / x * |s - sstar| := by
  have hx_pos : x > 0 := by linarith
  have hSs : Sp 4 s > 0 := Sp_pos 4 (by omega) s hs
  have hSt : Sp 4 sstar > 0 := Sp_pos 4 (by omega) sstar (le_of_lt hss_pos)
  have hfix : pandrosion_h x 4 sstar = sstar :=
    (pandrosion_fixed_point_iff x hx_pos 4 (by omega) sstar
      (le_of_lt hss_pos)).mpr hss_eq
  have hD : x * Sp 4 s * Sp 4 sstar > 0 := by positivity
  have hD_ne : x * Sp 4 s * Sp 4 sstar ≠ 0 := ne_of_gt hD
  have hdiff : pandrosion_h x 4 s - pandrosion_h x 4 sstar =
      (x - 1) * (Sp 4 s - Sp 4 sstar) / (x * Sp 4 s * Sp 4 sstar) := by
    unfold pandrosion_h; rw [eq_div_iff hD_ne]; field_simp; ring
  rw [← hfix, hdiff, Sp4_sub]
  have hQ : 1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2 > 0 := by positivity
  rw [show (x - 1) * ((s - sstar) *
        (1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2))
      = (x - 1) * (1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2) *
        (s - sstar) from by ring]
  rw [abs_div, abs_mul, abs_mul,
      abs_of_pos (by linarith : x - 1 > 0),
      abs_of_pos hQ, abs_of_pos hD, hfix]
  have hbound := Sp4_prod_bound s sstar hs (le_of_lt hss_pos)
  have hxm1 : x - 1 > 0 := by linarith
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero,
             pow_zero, add_zero] at hbound
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero,
             pow_zero, add_zero]
  rw [div_le_iff (by positivity)]
  have hsimp : (x - 1) / x * |s - sstar| *
      (x * (0 + 1 + s ^ 1 + s ^ 2 + s ^ 3) *
        (0 + 1 + sstar ^ 1 + sstar ^ 2 + sstar ^ 3))
    = (x - 1) * |s - sstar| *
      ((0 + 1 + s ^ 1 + s ^ 2 + s ^ 3) *
        (0 + 1 + sstar ^ 1 + sstar ^ 2 + sstar ^ 3)) := by
    field_simp; ring
  rw [hsimp]
  have hnn : (x - 1) * |s - sstar| ≥ 0 :=
    mul_nonneg (le_of_lt hxm1) (abs_nonneg _)
  calc (x - 1) * (1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2) *
         |s - sstar|
      = (x - 1) * |s - sstar| *
         (1 + s + sstar + s ^ 2 + s * sstar + sstar ^ 2) := by ring
    _ ≤ (x - 1) * |s - sstar| *
         ((0 + 1 + s ^ 1 + s ^ 2 + s ^ 3) *
           (0 + 1 + sstar ^ 1 + sstar ^ 2 + sstar ^ 3)) := by
        apply mul_le_mul_of_nonneg_left _ hnn
        simp only [zero_add, pow_one]; exact hbound
    _ = (x - 1) * |s - sstar| *
         ((0 + 1 + s ^ 1 + s ^ 2 + s ^ 3) *
           (0 + 1 + sstar ^ 1 + sstar ^ 2 + sstar ^ 3)) := rfl

/-- Strict distance decrease for `p = 4`. -/
theorem distance_decreases_p4 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 4 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x 4 s - sstar| < |s - sstar| := by
  have hbound := contraction_p4 x sstar hx hss_pos hss_lt hss_eq s hs hs_ne
  have hfrac : (x - 1) / x < 1 := contraction_ratio_lt_one x hx
  linarith [mul_lt_mul_of_pos_right hfrac
    (abs_pos.mpr (sub_ne_zero.mpr hs_ne))]

/-! ## §L11. Contraction for `p = 5` -/

/-- `S_5(s) − S_5(t) = (s − t) · (1 + s + t + s² + st + t² + s³ + s²t + st² + t³)`. -/
theorem Sp5_sub (s t : ℝ) :
    Sp 5 s - Sp 5 t = (s - t) *
    (1 + s + t + s^2 + s*t + t^2 + s^3 + s^2*t + s*t^2 + t^3) := by
  simp [Sp, Finset.sum_range_succ]; ring

/-- `Q₅(s,t) ≤ S_5(s) · S_5(t)` for `s, t ≥ 0`. -/
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
             mul_nonneg (mul_nonneg (mul_nonneg hs hs) hs)
               (mul_nonneg ht ht),
             mul_nonneg (mul_nonneg hs hs)
               (mul_nonneg (mul_nonneg ht ht) ht),
             mul_nonneg (mul_nonneg (mul_nonneg hs hs) hs)
               (mul_nonneg (mul_nonneg ht ht) ht)]

/-- Contraction for `p = 5`. -/
theorem contraction_p5 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 5 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (_hs_ne : s ≠ sstar) :
    |pandrosion_h x 5 s - sstar| ≤ (x - 1) / x * |s - sstar| := by
  have hx_pos : x > 0 := by linarith
  have hSs : Sp 5 s > 0 := Sp_pos 5 (by omega) s hs
  have hSt : Sp 5 sstar > 0 := Sp_pos 5 (by omega) sstar (le_of_lt hss_pos)
  have hfix : pandrosion_h x 5 sstar = sstar :=
    (pandrosion_fixed_point_iff x hx_pos 5 (by omega) sstar
      (le_of_lt hss_pos)).mpr hss_eq
  have hD : x * Sp 5 s * Sp 5 sstar > 0 := by positivity
  have hD_ne : x * Sp 5 s * Sp 5 sstar ≠ 0 := ne_of_gt hD
  have hdiff : pandrosion_h x 5 s - pandrosion_h x 5 sstar =
      (x - 1) * (Sp 5 s - Sp 5 sstar) / (x * Sp 5 s * Sp 5 sstar) := by
    unfold pandrosion_h; rw [eq_div_iff hD_ne]; field_simp; ring
  rw [← hfix, hdiff, Sp5_sub]
  have hQ : 1 + s + sstar + s^2 + s*sstar + sstar^2 + s^3 + s^2*sstar +
            s*sstar^2 + sstar^3 > 0 := by positivity
  rw [show (x-1) * ((s-sstar) *
        (1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3))
      = (x-1) *
        (1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3) *
        (s-sstar) from by ring]
  rw [abs_div, abs_mul, abs_mul,
      abs_of_pos (by linarith : x - 1 > 0),
      abs_of_pos hQ, abs_of_pos hD, hfix]
  have hbound := Sp5_prod_bound s sstar hs (le_of_lt hss_pos)
  have hxm1 : x - 1 > 0 := by linarith
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero,
             pow_zero, add_zero] at hbound
  simp only [Sp, Finset.sum_range_succ, Finset.sum_range_zero,
             pow_zero, add_zero]
  rw [div_le_iff (by positivity)]
  have hsimp : (x-1)/x * |s-sstar| *
      (x * (0+1+s^1+s^2+s^3+s^4) *
        (0+1+sstar^1+sstar^2+sstar^3+sstar^4))
    = (x-1) * |s-sstar| *
      ((0+1+s^1+s^2+s^3+s^4) *
        (0+1+sstar^1+sstar^2+sstar^3+sstar^4)) := by
    field_simp; ring
  rw [hsimp]
  have hnn : (x-1) * |s-sstar| ≥ 0 :=
    mul_nonneg (le_of_lt hxm1) (abs_nonneg _)
  calc (x-1) *
         (1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3) *
         |s-sstar|
      = (x-1) * |s-sstar| *
         (1+s+sstar+s^2+s*sstar+sstar^2+s^3+s^2*sstar+s*sstar^2+sstar^3) := by
         ring
    _ ≤ (x-1) * |s-sstar| *
         ((0+1+s^1+s^2+s^3+s^4) *
           (0+1+sstar^1+sstar^2+sstar^3+sstar^4)) := by
        apply mul_le_mul_of_nonneg_left _ hnn
        simp only [zero_add, pow_one]; exact hbound
    _ = (x-1) * |s-sstar| *
         ((0+1+s^1+s^2+s^3+s^4) *
           (0+1+sstar^1+sstar^2+sstar^3+sstar^4)) := rfl

/-- Strict distance decrease for `p = 5`. -/
theorem distance_decreases_p5 (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ 5 = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x 5 s - sstar| < |s - sstar| := by
  have hbound := contraction_p5 x sstar hx hss_pos hss_lt hss_eq s hs hs_ne
  have hfrac : (x - 1) / x < 1 := contraction_ratio_lt_one x hx
  linarith [mul_lt_mul_of_pos_right hfrac
    (abs_pos.mpr (sub_ne_zero.mpr hs_ne))]

/-! ## §L12. The universal divided difference `Qp` -/

/-- Divided difference `Qp(s,t) = Σ_{k<p} Σ_{j<k} s^j · t^{k-1-j}`. -/
noncomputable def Qp (p : ℕ) (s t : ℝ) : ℝ :=
  ∑ k in range p, ∑ j in range k, s ^ j * t ^ (k - 1 - j)

/-- **Factoring: `S_p(s) − S_p(t) = (s − t)·Qp(s,t)`.** -/
theorem Sp_sub_factor (p : ℕ) (s t : ℝ) :
    Sp p s - Sp p t = (s - t) * Qp p s t := by
  unfold Sp Qp
  rw [← Finset.sum_sub_distrib]
  rw [Finset.mul_sum]
  congr 1; ext k
  exact pow_sub_factor s t k

/-- `Qp(s, t) ≥ 0` for `s, t ≥ 0`. -/
theorem Qp_nonneg (p : ℕ) (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    Qp p s t ≥ 0 := by
  unfold Qp
  apply Finset.sum_nonneg; intro k _
  apply Finset.sum_nonneg; intro j _
  exact mul_nonneg (pow_nonneg hs _) (pow_nonneg ht _)

/-- `Qp(s,t) ≤ S_p(s) · S_p(t)` for `s, t ≥ 0`. The `p(p-1)/2`
    monomials of `Qp` sit inside the `p²` monomials of the product. -/
theorem Qp_le_Sp_mul (p : ℕ) (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    Qp p s t ≤ Sp p s * Sp p t := by
  unfold Sp Qp
  rw [Finset.sum_mul_sum]
  calc ∑ k in range p, ∑ j in range k, s ^ j * t ^ (k - 1 - j)
      = ∑ ab in (range p ×ˢ range p).filter
          (fun ab : ℕ × ℕ => ab.1 + ab.2 + 2 ≤ p),
          s ^ ab.1 * t ^ ab.2 := by
        induction p with
        | zero =>
          simp only [Finset.range_zero, Finset.sum_empty]
          rw [Finset.empty_product]
          simp
        | succ n ihn =>
          rw [Finset.sum_range_succ]
          conv_rhs =>
            rw [show (range (n + 1) ×ˢ range (n + 1)).filter
                  (fun ab : ℕ × ℕ => ab.1 + ab.2 + 2 ≤ n + 1)
              = (range n ×ˢ range n).filter
                  (fun ab : ℕ × ℕ => ab.1 + ab.2 + 2 ≤ n)
              ∪ (range n).map
                  ⟨fun j => (j, n - 1 - j),
                   fun a b h => by simp [Prod.ext_iff] at h; omega⟩
              from by
                ext ⟨a, b⟩
                simp only [Finset.mem_filter, Finset.mem_product,
                           Finset.mem_range, Finset.mem_union,
                           Finset.mem_map, Function.Embedding.coeFn_mk]
                constructor
                · intro ⟨⟨ha, hb⟩, hab⟩
                  by_cases h : a + b + 2 ≤ n
                  · left; exact ⟨⟨by omega, by omega⟩, h⟩
                  · right
                    refine ⟨a, by omega, ?_⟩
                    exact Prod.ext rfl (by omega)
                · intro h
                  rcases h with ⟨⟨ha, hb⟩, hab⟩ | ⟨j, hj, hprod⟩
                  · exact ⟨⟨by omega, by omega⟩, by omega⟩
                  · simp only [Prod.mk.injEq] at hprod
                    exact ⟨⟨by omega, by omega⟩, by omega⟩]
          rw [Finset.sum_union]
          · conv_rhs => rw [Finset.sum_map]
            rw [← ihn]
            congr 1
          · rw [Finset.disjoint_left]
            intro ⟨a, b⟩ hmem1 hmem2
            simp only [Finset.mem_filter, Finset.mem_product,
                       Finset.mem_range] at hmem1
            simp only [Finset.mem_map,
                       Function.Embedding.coeFn_mk] at hmem2
            obtain ⟨j, _, hprod⟩ := hmem2
            simp only [Prod.mk.injEq] at hprod
            omega
    _ ≤ ∑ ab in (range p ×ˢ range p), s ^ ab.1 * t ^ ab.2 := by
        apply Finset.sum_le_sum_of_subset_of_nonneg
          (Finset.filter_subset _ _)
        intro i _ _
        exact mul_nonneg (pow_nonneg hs _) (pow_nonneg ht _)
    _ = _ := by rw [Finset.sum_product]

/-! ## §L13. The general contraction theorem -/

/-- **General contraction theorem.** For every `p ≥ 2`, every `x > 1`,
    and every `s ≥ 0`: `|h(s) − s*| ≤ (x−1)/x · |s − s*|`. -/
theorem contraction_general (p : ℕ) (hp : p ≥ 2) (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ p = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (_hs_ne : s ≠ sstar) :
    |pandrosion_h x p s - sstar| ≤ (x - 1) / x * |s - sstar| := by
  have hx_pos : x > 0 := by linarith
  have hSs : Sp p s > 0 := Sp_pos p (by omega) s hs
  have hSt : Sp p sstar > 0 := Sp_pos p (by omega) sstar (le_of_lt hss_pos)
  have hfix : pandrosion_h x p sstar = sstar :=
    (pandrosion_fixed_point_iff x hx_pos p (by omega) sstar
      (le_of_lt hss_pos)).mpr hss_eq
  have hD_ne : x * Sp p s * Sp p sstar ≠ 0 := by positivity
  have hdiff : pandrosion_h x p s - pandrosion_h x p sstar =
      (x - 1) * (Sp p s - Sp p sstar) / (x * Sp p s * Sp p sstar) := by
    unfold pandrosion_h; rw [eq_div_iff hD_ne]; field_simp; ring
  rw [← hfix, hdiff, Sp_sub_factor]
  have hQnn : Qp p s sstar ≥ 0 := Qp_nonneg p s sstar hs (le_of_lt hss_pos)
  have hQbound : Qp p s sstar ≤ Sp p s * Sp p sstar :=
    Qp_le_Sp_mul p s sstar hs (le_of_lt hss_pos)
  rw [show (x - 1) * ((s - sstar) * Qp p s sstar) =
      (x - 1) * Qp p s sstar * (s - sstar) from by ring]
  rw [abs_div, abs_mul, abs_mul,
      abs_of_pos (by linarith : x - 1 > 0),
      abs_of_nonneg hQnn,
      abs_of_pos (by positivity : x * Sp p s * Sp p sstar > 0)]
  rw [hfix]
  rw [div_le_iff (by positivity : x * Sp p s * Sp p sstar > 0)]
  have hsimp : (x - 1) / x * |s - sstar| * (x * Sp p s * Sp p sstar)
      = (x - 1) * |s - sstar| * (Sp p s * Sp p sstar) := by
    field_simp; ring
  rw [hsimp]
  have hnn : (x - 1) * |s - sstar| ≥ 0 :=
    mul_nonneg (by linarith) (abs_nonneg _)
  calc (x - 1) * Qp p s sstar * |s - sstar|
      = (x - 1) * |s - sstar| * Qp p s sstar := by ring
    _ ≤ (x - 1) * |s - sstar| * (Sp p s * Sp p sstar) :=
        mul_le_mul_of_nonneg_left hQbound hnn

/-- **Strict distance decrease for every `p ≥ 2`.** -/
theorem distance_decreases_general (p : ℕ) (hp : p ≥ 2) (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ p = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x p s - sstar| < |s - sstar| := by
  have hbound := contraction_general p hp x sstar hx hss_pos hss_lt hss_eq s hs hs_ne
  have hfrac : (x - 1) / x < 1 := contraction_ratio_lt_one x hx
  linarith [mul_lt_mul_of_pos_right hfrac
    (abs_pos.mpr (sub_ne_zero.mpr hs_ne))]

end Pandrosion
