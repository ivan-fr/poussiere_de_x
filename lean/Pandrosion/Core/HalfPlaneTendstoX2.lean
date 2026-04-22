/-
  Universitas Pandrosion — §67. **h-contraction sur le demi-plan droit
  à `x = 2`.**

  Extension de §64 (disque `B(α₀, 1/4)`) au demi-plan `Re z ≥ 1/2`.
  Couverture inconditionnelle de ~50 % du plan complexe.

  Contents.
    §67.1  `half_plane_quadratic_bound_at_x2` — inégalité clef.
    §67.2  `pandrosion_h_C_p3_half_plane_contraction_x2` — K ≤ 1/2.
    §67.3  `pandrosion_h_C_p3_half_plane_invariant_x2` — Re h(z) ≥ 1/2.
    §67.4  `pandrosion_h_C_p3_half_plane_orbit_bound_x2` — (1/2)^n.
    §67.5  `pandrosion_h_C_p3_tendsto_half_plane_x2` — Tendsto.
-/

import Pandrosion.Core.BanachX2Concrete

namespace Pandrosion

open Complex Filter Topology

/-! ======================== §67.0 Helpers ======================== -/

private theorem alphaX2_C_ne_zero_loc : ((alphaX2 : ℝ) : ℂ) ≠ 0 := by
  rw [Ne, Complex.ofReal_eq_zero]; exact ne_of_gt alphaX2_pos

private theorem Sp_C_p3_alphaX2_ne_zero_loc :
    Sp_C 3 ((alphaX2 : ℝ) : ℂ) ≠ 0 := by
  rw [Sp_C_p3]
  have h_real : (1 : ℂ) + ((alphaX2 : ℝ) : ℂ) + ((alphaX2 : ℝ) : ℂ) ^ 2
              = (((1 + alphaX2 + alphaX2 ^ 2) : ℝ) : ℂ) := by push_cast; ring
  rw [h_real, Ne, Complex.ofReal_eq_zero]
  have := sq_nonneg alphaX2; linarith [alphaX2_pos]

private theorem alphaX2_cube_C_loc : ((alphaX2 : ℝ) : ℂ) ^ 3 = 1 / (2 : ℂ) := by
  have h_cast : ((alphaX2 : ℝ) : ℂ) ^ 3 = (((alphaX2 ^ 3) : ℝ) : ℂ) := by
    push_cast; ring
  rw [h_cast, alphaX2_cube]; push_cast; ring

/-! ======================== §67.1 Inégalité clef ======================== -/

/-- `‖1 + z + α₀‖² ≤ (1+α₀+α₀²)² · ‖Sp_C 3 z‖²` sur `Re z ≥ 1/2`. -/
theorem half_plane_quadratic_bound_at_x2 (z : ℂ) (hz : z.re ≥ 1/2) :
    ‖1 + z + ((alphaX2 : ℝ) : ℂ)‖ ^ 2 ≤
      (1 + alphaX2 + alphaX2 ^ 2) ^ 2 * ‖Sp_C 3 z‖ ^ 2 := by
  have hα34 : alphaX2 ≥ 3/4 := alphaX2_ge_three_fourths
  have hα1  : alphaX2 ≤ 1 := alphaX2_le_one
  have hαp  : (0 : ℝ) < alphaX2 := alphaX2_pos
  have hId := Sp_C_p3_norm_sq_identity z
  have hNum :
      ‖1 + z + ((alphaX2 : ℝ) : ℂ)‖ ^ 2
        = (1 + z.re + alphaX2) ^ 2 + z.im ^ 2 := by
    rw [Complex.norm_eq_abs, Complex.sq_abs, Complex.normSq_apply]
    simp [Complex.add_re, Complex.add_im, Complex.one_re, Complex.one_im,
          Complex.ofReal_re, Complex.ofReal_im]; ring
  rw [hNum]
  set a := z.re
  set b := z.im
  set α := alphaX2
  have hSpSq : ‖Sp_C 3 z‖ ^ 2
             = (1 + a + a^2)^2 + b^2 * (2*a^2 + 2*a - 1 + b^2) := by linarith
  rw [hSpSq]
  -- (A) 1+a+α ≤ (1+α+α²)(1+a+a²).
  have hα_nn : 0 ≤ α := hαp.le
  have ha_nn : 0 ≤ a := by linarith
  have hA : (1 + a + α) ≤ (1 + α + α^2) * (1 + a + a^2) := by
    have h_mid : 0 ≤ α * a * (1 + α) :=
      mul_nonneg (mul_nonneg hα_nn ha_nn) (by linarith)
    have h_sum_nn : 0 ≤ 1 + α + α^2 := by
      have := sq_nonneg α; linarith
    have h_last : 0 ≤ a^2 * (1 + α + α^2) :=
      mul_nonneg (sq_nonneg _) h_sum_nn
    linarith [sq_nonneg α,
      show (1 + α + α^2) * (1 + a + a^2) - (1 + a + α)
         = α^2 + α * a * (1 + α) + a^2 * (1 + α + α^2) from by ring]
  have hA_nn : 0 ≤ 1 + a + α := by linarith
  have hA_sq : (1 + a + α)^2 ≤ (1 + α + α^2)^2 * (1 + a + a^2)^2 := by
    have h_sq_prod : (1 + a + α)^2 ≤ ((1 + α + α^2) * (1 + a + a^2))^2 :=
      pow_le_pow_left hA_nn hA 2
    linarith [h_sq_prod,
      show ((1 + α + α^2) * (1 + a + a^2))^2
         = (1 + α + α^2)^2 * (1 + a + a^2)^2 from by ring]
  -- (B) b² ≤ (1+α+α²)²·b²·(2a²+2a-1+b²).
  have hB_inner : (1/2 : ℝ) ≤ 2*a^2 + 2*a - 1 + b^2 := by
    have h1 : 0 ≤ 2*(a - 1/2)^2 := by
      have : 0 ≤ (a - 1/2)^2 := sq_nonneg _; linarith
    linarith [sq_nonneg b,
      show 2*a^2 + 2*a - 1 + b^2 - 1/2
         = 2*(a - 1/2)^2 + 4*(a - 1/2) + b^2 from by ring]
  have hB_P : (37/16 : ℝ) ≤ 1 + α + α^2 := by
    linarith [sq_nonneg (α - 3/4),
      show 1 + α + α^2 - 37/16 = (α - 3/4)^2 + (5/2)*(α - 3/4) from by ring]
  have hB_Psq : (37/16 : ℝ)^2 ≤ (1 + α + α^2)^2 :=
    pow_le_pow_left (by norm_num) hB_P 2
  have hB : (1 : ℝ) ≤ (1 + α + α^2)^2 * (2*a^2 + 2*a - 1 + b^2) := by
    have h_prod :
        (37/16 : ℝ)^2 * (1/2) ≤ (1 + α + α^2)^2 * (2*a^2 + 2*a - 1 + b^2) :=
      mul_le_mul hB_Psq hB_inner (by norm_num) (sq_nonneg _)
    linarith [show (37/16 : ℝ)^2 * (1/2) = 1369/512 from by norm_num]
  have hB_b : b^2 ≤ (1 + α + α^2)^2 * (b^2 * (2*a^2 + 2*a - 1 + b^2)) := by
    have h_mul := mul_le_mul_of_nonneg_left hB (sq_nonneg b)
    linarith [h_mul,
      show b^2 * ((1 + α + α^2)^2 * (2*a^2 + 2*a - 1 + b^2))
         = (1 + α + α^2)^2 * (b^2 * (2*a^2 + 2*a - 1 + b^2)) from by ring]
  -- Combine.
  linarith [hA_sq, hB_b,
    show (1 + α + α^2)^2 * ((1 + a + a^2)^2 + b^2 * (2*a^2 + 2*a - 1 + b^2))
       = (1 + α + α^2)^2 * (1 + a + a^2)^2
       + (1 + α + α^2)^2 * (b^2 * (2*a^2 + 2*a - 1 + b^2)) from by ring]

/-! ======================== §67.2 Contraction ≤ 1/2 ======================== -/

theorem pandrosion_h_C_p3_half_plane_contraction_x2
    (z : ℂ) (hz : z.re ≥ 1/2) :
    ‖pandrosion_h_C (2 : ℂ) 3 z - ((alphaX2 : ℝ) : ℂ)‖
      ≤ (1/2) * ‖z - ((alphaX2 : ℝ) : ℂ)‖ := by
  have h_Sα_ne : Sp_C 3 ((alphaX2 : ℝ) : ℂ) ≠ 0 := Sp_C_p3_alphaX2_ne_zero_loc
  have h_Sz_ne : Sp_C 3 z ≠ 0 := Sp_C_p3_ne_zero_of_re_ge_half z hz
  have h_2_ne : ((2 : ℂ)) ≠ 0 := by norm_num
  have h_α_pow : ((alphaX2 : ℝ) : ℂ) ^ 3 = 1 / (2 : ℂ) := alphaX2_cube_C_loc
  have hDiff := pandrosion_h_C_p3_diff_identity
    (2 : ℂ) ((alphaX2 : ℝ) : ℂ) z h_2_ne h_Sz_ne h_Sα_ne h_α_pow
  have h_2_norm : ‖(2 : ℂ)‖ = 2 := by
    rw [show ((2 : ℂ)) = (((2 : ℝ)) : ℂ) from by norm_num,
        Complex.norm_real, Real.norm_eq_abs]; norm_num
  have hα_pos : 0 < 1 + alphaX2 + alphaX2 ^ 2 := by
    have := sq_nonneg alphaX2; linarith [alphaX2_pos]
  have h_Sα_norm : ‖Sp_C 3 ((alphaX2 : ℝ) : ℂ)‖ = 1 + alphaX2 + alphaX2 ^ 2 := by
    rw [Sp_C_p3]
    have h_eq : (1 : ℂ) + ((alphaX2 : ℝ) : ℂ) + ((alphaX2 : ℝ) : ℂ) ^ 2
              = (((1 + alphaX2 + alphaX2 ^ 2) : ℝ) : ℂ) := by push_cast; ring
    rw [h_eq, Complex.norm_real, Real.norm_eq_abs]; exact abs_of_pos hα_pos
  have hSzp : 0 < ‖Sp_C 3 z‖ := by
    rcases (norm_nonneg (Sp_C 3 z)).lt_or_eq with hp | hp
    · exact hp
    · exfalso; exact h_Sz_ne (norm_eq_zero.mp hp.symm)
  -- Take norm of hDiff.
  have h_norm_hz :
      ‖pandrosion_h_C (2 : ℂ) 3 z - ((alphaX2 : ℝ) : ℂ)‖
        = ‖1 + z + ((alphaX2 : ℝ) : ℂ)‖ * ‖z - ((alphaX2 : ℝ) : ℂ)‖
            / (2 * ‖Sp_C 3 z‖ * (1 + alphaX2 + alphaX2 ^ 2)) := by
    rw [hDiff]
    rw [norm_div]
    have h_num : ‖((2 : ℂ) - 1) * (1 + z + ((alphaX2 : ℝ) : ℂ)) *
                    (z - ((alphaX2 : ℝ) : ℂ))‖
               = ‖1 + z + ((alphaX2 : ℝ) : ℂ)‖ * ‖z - ((alphaX2 : ℝ) : ℂ)‖ := by
      rw [show ((2 : ℂ)) - 1 = 1 from by norm_num]
      rw [norm_mul, norm_mul, norm_one, one_mul]
    rw [h_num]
    have h_den : ‖(2 : ℂ) * Sp_C 3 z * Sp_C 3 ((alphaX2 : ℝ) : ℂ)‖
               = 2 * ‖Sp_C 3 z‖ * (1 + alphaX2 + alphaX2 ^ 2) := by
      rw [norm_mul, norm_mul, h_2_norm, h_Sα_norm]
    rw [h_den]
  rw [h_norm_hz]
  have h_den_pos : 0 < 2 * ‖Sp_C 3 z‖ * (1 + alphaX2 + alphaX2 ^ 2) :=
    mul_pos (mul_pos (by norm_num) hSzp) hα_pos
  rw [div_le_iff h_den_pos]
  -- Now: ‖1+z+α‖·‖z-α‖ ≤ (1/2)·‖z-α‖·(2·‖Sp‖·(1+α+α²))
  --                   = ‖Sp(z)‖·(1+α+α²)·‖z-α‖
  have h_rhs_rw :
      (1/2) * ‖z - ((alphaX2 : ℝ) : ℂ)‖
        * (2 * ‖Sp_C 3 z‖ * (1 + alphaX2 + alphaX2 ^ 2))
      = ((1 + alphaX2 + alphaX2 ^ 2) * ‖Sp_C 3 z‖)
          * ‖z - ((alphaX2 : ℝ) : ℂ)‖ := by ring
  rw [h_rhs_rw]
  -- §67.1 ⟹ ‖1+z+α‖ ≤ (1+α+α²)·‖Sp(z)‖ via sqrt.
  have hKeySq := half_plane_quadratic_bound_at_x2 z hz
  have hAN : 0 ≤ (1 + alphaX2 + alphaX2 ^ 2) * ‖Sp_C 3 z‖ :=
    mul_nonneg hα_pos.le (norm_nonneg _)
  have hASqEq : ((1 + alphaX2 + alphaX2 ^ 2) * ‖Sp_C 3 z‖) ^ 2
              = (1 + alphaX2 + alphaX2 ^ 2) ^ 2 * ‖Sp_C 3 z‖ ^ 2 := by ring
  have hKey : ‖1 + z + ((alphaX2 : ℝ) : ℂ)‖
              ≤ (1 + alphaX2 + alphaX2 ^ 2) * ‖Sp_C 3 z‖ := by
    have hS := Real.sqrt_le_sqrt
               (show ‖1 + z + ((alphaX2 : ℝ) : ℂ)‖ ^ 2
                 ≤ ((1 + alphaX2 + alphaX2 ^ 2) * ‖Sp_C 3 z‖) ^ 2 by
                 rw [hASqEq]; exact hKeySq)
    rwa [Real.sqrt_sq (norm_nonneg _), Real.sqrt_sq hAN] at hS
  exact mul_le_mul_of_nonneg_right hKey (norm_nonneg _)

/-! ======================== §67.3 Invariance ======================== -/

theorem pandrosion_h_C_p3_half_plane_invariant_x2
    (z : ℂ) (hz : z.re ≥ 1/2) :
    (pandrosion_h_C (2 : ℂ) 3 z).re ≥ 1/2 := by
  have h_Sz_lb_sq : ‖Sp_C 3 z‖ ^ 2 ≥ (1 + z.re + z.re ^ 2) ^ 2 :=
    Sp_C_p3_norm_sq_lower_bound z hz
  have h_rhs_nn : 0 ≤ 1 + z.re + z.re ^ 2 := by
    have := sq_nonneg z.re; linarith
  have h_Sz_nn : 0 ≤ ‖Sp_C 3 z‖ := norm_nonneg _
  have h_Sz_lb : ‖Sp_C 3 z‖ ≥ 1 + z.re + z.re ^ 2 := by
    have hS := Real.sqrt_le_sqrt h_Sz_lb_sq
    rwa [Real.sqrt_sq h_rhs_nn, Real.sqrt_sq h_Sz_nn] at hS
  have h_74 : 1 + z.re + z.re ^ 2 ≥ 7/4 := by
    linarith [sq_nonneg (z.re - 1/2),
      show 1 + z.re + z.re ^ 2 - 7/4
         = (z.re - 1/2)^2 + 2 * (z.re - 1/2) from by ring]
  have h2N : ‖(2 : ℂ)‖ = 2 := by
    rw [show ((2 : ℂ)) = (((2 : ℝ)) : ℂ) from by norm_num,
        Complex.norm_real, Real.norm_eq_abs]; norm_num
  have h_inv_bound : ‖(1 : ℂ) / ((2 : ℂ) * Sp_C 3 z)‖ ≤ 2/7 := by
    rw [norm_div, norm_mul, norm_one, h2N,
        div_le_iff (by linarith : (0:ℝ) < 2 * ‖Sp_C 3 z‖)]
    linarith
  have h_re_decomp : (pandrosion_h_C (2 : ℂ) 3 z).re
                   = 1 - ((1 : ℂ) / ((2 : ℂ) * Sp_C 3 z)).re := by
    unfold pandrosion_h_C
    simp [Complex.sub_re, Complex.one_re,
          show ((2 : ℂ)) - 1 = 1 from by norm_num]
  rw [h_re_decomp]
  have h_abs_le : |((1 : ℂ) / ((2 : ℂ) * Sp_C 3 z)).re|
                 ≤ ‖(1 : ℂ) / ((2 : ℂ) * Sp_C 3 z)‖ :=
    Complex.abs_re_le_abs _
  linarith [le_abs_self (((1 : ℂ) / ((2 : ℂ) * Sp_C 3 z)).re),
            h_abs_le, h_inv_bound]

/-! ======================== §67.4 Orbit bound ======================== -/

theorem pandrosion_h_C_p3_half_plane_orbit_bound_x2
    (z : ℂ) (hz : z.re ≥ 1/2) (n : ℕ) :
    ‖(pandrosion_h_C (2 : ℂ) 3)^[n] z - ((alphaX2 : ℝ) : ℂ)‖
      ≤ (1/2) ^ n * ‖z - ((alphaX2 : ℝ) : ℂ)‖ := by
  induction n with
  | zero => simp
  | succ k ih =>
    have h_inv : ∀ m, ((pandrosion_h_C (2 : ℂ) 3)^[m] z).re ≥ 1/2 := by
      intro m
      induction m with
      | zero => simpa using hz
      | succ j ihj =>
        rw [Function.iterate_succ', Function.comp_apply]
        exact pandrosion_h_C_p3_half_plane_invariant_x2 _ ihj
    have h_step := pandrosion_h_C_p3_half_plane_contraction_x2 _ (h_inv k)
    rw [Function.iterate_succ', Function.comp_apply]
    have h_chain :
        (1/2) * ‖(pandrosion_h_C (2 : ℂ) 3)^[k] z - ((alphaX2 : ℝ) : ℂ)‖
        ≤ (1/2) * ((1/2) ^ k * ‖z - ((alphaX2 : ℝ) : ℂ)‖) :=
      mul_le_mul_of_nonneg_left ih (by norm_num)
    have h_pow_succ :
        (1/2 : ℝ) * ((1/2) ^ k * ‖z - ((alphaX2 : ℝ) : ℂ)‖)
        = (1/2) ^ (k + 1) * ‖z - ((alphaX2 : ℝ) : ℂ)‖ := by
      rw [pow_succ]; ring
    show ‖pandrosion_h_C (2 : ℂ) 3 ((pandrosion_h_C (2 : ℂ) 3)^[k] z)
            - ((alphaX2 : ℝ) : ℂ)‖
          ≤ (1/2) ^ (k + 1) * ‖z - ((alphaX2 : ℝ) : ℂ)‖
    linarith [h_step, h_chain, h_pow_succ.le, h_pow_succ.ge]

/-! ======================== §67.5 Tendsto ======================== -/

/-- **★★★★★★ h-convergence sur le demi-plan droit entier à x = 2.** -/
theorem pandrosion_h_C_p3_tendsto_half_plane_x2
    (z : ℂ) (hz : z.re ≥ 1/2) :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  rw [tendsto_iff_norm_sub_tendsto_zero]
  have h_orbit := pandrosion_h_C_p3_half_plane_orbit_bound_x2 z hz
  have h_pow :
      Tendsto (fun n : ℕ => (1/2 : ℝ) ^ n) atTop (𝓝 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)
  exact squeeze_zero (fun _ => norm_nonneg _) h_orbit
    (by simpa using h_pow.mul_const ‖z - ((alphaX2 : ℝ) : ℂ)‖)

end Pandrosion
