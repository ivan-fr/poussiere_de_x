/-
  Universitas Pandrosion — §84. **Invariance demi-plan pour
  `pandrosion_h_C x 3` à `x` réel ∈ (1, 8].**

  Généralisation paramétrée de §67.3 (qui était fixé à `x = 2`).
  Pour `x` variable, `Re h(z) ≥ 1/2 ⟺ (x-1)/(x·‖Sp(z)‖) ≤ 1/2`.
  Sur le demi-plan `Re z ≥ 1/2`, `‖Sp(z)‖ ≥ 7/4`, donc cette borne
  se ramène à `(x-1)/(x·7/4) = 4(x-1)/(7x) ≤ 1/2 ⟺ x ≤ 8`.

  **Importance** : démontre que la chaîne §67 → §57 (Banach abstrait)
  s'applique uniformément sur la plage `x ∈ (1, 8]` (qui correspond
  à `α(x) ∈ [1/2, 1)` via §83).

  Contents.

    §84.1  `pandrosion_h_C_p3_half_plane_invariant_generic_x` —
           `Re h(x, 3, z) ≥ 1/2` pour `x ∈ (1, 8], Re z ≥ 1/2`.
-/

import Pandrosion.Core.HalfPlaneTendstoX2
import Pandrosion.Core.AlphaGenericX

namespace Pandrosion

open Complex

/-- **Invariance demi-plan paramétrée** : pour tout `x ∈ ℝ` avec
    `1 < x ≤ 8` et tout `z ∈ ℂ` avec `Re z ≥ 1/2`,

        `Re (pandrosion_h_C x 3 z) ≥ 1/2`.

    Preuve directe : majoration `(x-1)/(x·‖Sp‖) ≤ 1/2` via
    `‖Sp‖ ≥ 7/4` et `x ≤ 8`. -/
theorem pandrosion_h_C_p3_half_plane_invariant_generic_x
    (x : ℝ) (hx_lo : 1 < x) (hx_hi : x ≤ 8) (z : ℂ) (hz : z.re ≥ 1/2) :
    (pandrosion_h_C ((x : ℂ)) 3 z).re ≥ 1/2 := by
  -- ‖Sp(z)‖ ≥ 7/4 sur demi-plan.
  have h_Sz_lb_sq : ‖Sp_C 3 z‖^2 ≥ (1 + z.re + z.re^2)^2 :=
    Sp_C_p3_norm_sq_lower_bound z hz
  have h_rhs_nn : 0 ≤ 1 + z.re + z.re^2 := by
    have := sq_nonneg z.re; linarith
  have h_Sz_lb : ‖Sp_C 3 z‖ ≥ 1 + z.re + z.re^2 := by
    have hS := Real.sqrt_le_sqrt h_Sz_lb_sq
    rwa [Real.sqrt_sq h_rhs_nn, Real.sqrt_sq (norm_nonneg _)] at hS
  have h_74 : 1 + z.re + z.re^2 ≥ 7/4 := by
    linarith [sq_nonneg (z.re - 1/2),
      show 1 + z.re + z.re^2 - 7/4 = (z.re - 1/2)^2 + 2*(z.re - 1/2) from by ring]
  have h_Sz_74 : ‖Sp_C 3 z‖ ≥ 7/4 := by linarith
  -- Normes de x et x-1 (réels positifs).
  have h_x_pos : (0 : ℝ) < x := by linarith
  have h_xN : ‖((x : ℂ))‖ = x := by
    rw [show ((x : ℂ)) = (((x : ℝ)) : ℂ) from rfl,
        Complex.norm_real, Real.norm_eq_abs]
    exact abs_of_pos h_x_pos
  have h_x_sub_one_pos : (0 : ℝ) < x - 1 := by linarith
  have h_x_sub_one_N : ‖((x : ℂ)) - 1‖ = x - 1 := by
    have h_eq : ((x : ℂ)) - 1 = (((x - 1 : ℝ)) : ℂ) := by push_cast; ring
    rw [h_eq, Complex.norm_real, Real.norm_eq_abs]
    exact abs_of_pos h_x_sub_one_pos
  -- Re h(z) = 1 - Re((x-1)/(x·Sp(z))).
  have h_re_decomp : (pandrosion_h_C ((x : ℂ)) 3 z).re
                   = 1 - (((x : ℂ) - 1) / ((x : ℂ) * Sp_C 3 z)).re := by
    unfold pandrosion_h_C
    simp [Complex.sub_re, Complex.one_re]
  rw [h_re_decomp]
  -- Norme : ‖(x-1)/(x·Sp(z))‖ = (x-1)/(x·‖Sp(z)‖).
  have h_xSz_pos : 0 < x * ‖Sp_C 3 z‖ := mul_pos h_x_pos (by linarith)
  have h_norm_eq : ‖((x : ℂ) - 1) / ((x : ℂ) * Sp_C 3 z)‖
                 = (x - 1) / (x * ‖Sp_C 3 z‖) := by
    rw [norm_div, norm_mul, h_xN, h_x_sub_one_N]
  -- Re ≤ |Re| ≤ ‖·‖ = (x-1)/(x·‖Sp‖).
  have h_re_le_norm : (((x : ℂ) - 1) / ((x : ℂ) * Sp_C 3 z)).re
                  ≤ (x - 1) / (x * ‖Sp_C 3 z‖) := by
    rw [← h_norm_eq]
    exact le_trans (le_abs_self _) (Complex.abs_re_le_abs _)
  -- (x-1)/(x·‖Sp‖) ≤ 1/2 sur x ≤ 8 et ‖Sp‖ ≥ 7/4.
  have h_div_le_half : (x - 1) / (x * ‖Sp_C 3 z‖) ≤ 1/2 := by
    rw [div_le_iff h_xSz_pos]
    -- Want : x - 1 ≤ (1/2)·(x·‖Sp(z)‖). On a `x·‖Sp‖ ≥ 7x/4`.
    have h_xSz_ge_raw : x * (7/4) ≤ x * ‖Sp_C 3 z‖ :=
      mul_le_mul_of_nonneg_left h_Sz_74 (le_of_lt h_x_pos)
    have h_xSz_ge : x * ‖Sp_C 3 z‖ ≥ 7*x/4 := by linarith
    -- 2(x-1) ≤ 7x/4 ⟺ 8(x-1) ≤ 7x ⟺ x ≤ 8. ✓
    linarith
  linarith

end Pandrosion
