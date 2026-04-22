/-
  Universitas Pandrosion — §57. **Quantitative contraction bound
  for `pandrosion_h_C x 3` on a disk.**

  Chaîne §54 + §55 + §56 → borne quantitative multiplicative :
  sur le disque `B(α₀, r₀)` avec `α₀ ≥ 3/4`, `r₀ ≤ α₀ − 1/2`,
  `x > 1`, `α₀³ = 1/x`:

      `‖h(z) − α₀‖ · [x · D(α₀, r₀) · Sp₀]
         ≤ (x − 1) · N(α₀, r₀) · ‖z − α₀‖`

  où `D(α₀, r₀) := 1 + (α₀ − r₀) + (α₀ − r₀)²`, `Sp₀ := 1 + α₀ + α₀²`,
  `N(α₀, r₀) := 1 + 2α₀ + r₀`. En divisant par le facteur positif
  gauche, cela donne le K < 1 explicite.

  Module **court et focalisé** : petits lemmes isolés, preuves
  par `ring_nf` + `nlinarith` bornés, sans `calc` imbriqués.

  Contents.

    §57.1  `sp_lb_on_disk` — `‖Sp_C 3 z‖ ≥ 1 + (α₀-r₀) + (α₀-r₀)²`.

    §57.2  `q_ub_on_disk` — `‖1 + z + α₀‖ ≤ 1 + 2α₀ + r₀`.

    §57.3  `pandrosion_h_C_p3_quantitative_bound` — borne
           multiplicative finale.
-/

import Pandrosion.Core.HalfPlaneContractionP3Complete

namespace Pandrosion

open Complex

/-! ============================================================
  Helpers (small isolated lemmas)
============================================================ -/

section HelpersC

/-- Re z ≥ α₀ − r₀ sur le disque `B(α₀, r₀)`. -/
private theorem re_ge_of_disk (α₀ r₀ : ℝ) (z : ℂ)
    (hz : ‖z - ((α₀ : ℂ))‖ ≤ r₀) : z.re ≥ α₀ - r₀ := by
  have h_re_sub : (z - ((α₀ : ℂ))).re = z.re - α₀ := by simp [Complex.sub_re]
  have h_abs_re : |z.re - α₀| ≤ ‖z - ((α₀ : ℂ))‖ := by
    rw [← h_re_sub]; exact Complex.abs_re_le_abs _
  have := (abs_le.mp (le_trans h_abs_re hz)).1
  linarith

/-- ‖z‖ ≤ α₀ + r₀ sur le disque. -/
private theorem norm_ub_of_disk (α₀ r₀ : ℝ) (hα : 0 < α₀) (z : ℂ)
    (hz : ‖z - ((α₀ : ℂ))‖ ≤ r₀) : ‖z‖ ≤ α₀ + r₀ := by
  have h_α : ‖((α₀ : ℂ))‖ = α₀ := by
    rw [show ((α₀ : ℂ)) = ((α₀ : ℝ) : ℂ) from rfl, Complex.norm_real,
        Real.norm_eq_abs]
    exact abs_of_pos hα
  have h_tri : ‖z‖ ≤ ‖z - ((α₀ : ℂ))‖ + ‖((α₀ : ℂ))‖ := by
    have h_eq : z = (z - ((α₀ : ℂ))) + ((α₀ : ℂ)) := by ring
    calc ‖z‖ = ‖(z - ((α₀ : ℂ))) + ((α₀ : ℂ))‖ := by rw [← h_eq]
      _ ≤ ‖z - ((α₀ : ℂ))‖ + ‖((α₀ : ℂ))‖ := norm_add_le _ _
  rw [h_α] at h_tri
  linarith

/-- `‖Sp_C 3 α₀‖ = 1 + α₀ + α₀²` (réel positif). -/
private theorem sp_α_norm (α₀ : ℝ) (hα : 0 < α₀) :
    ‖Sp_C 3 ((α₀ : ℂ))‖ = 1 + α₀ + α₀ ^ 2 := by
  have h_eq : Sp_C 3 ((α₀ : ℂ)) = (((1 + α₀ + α₀ ^ 2) : ℝ) : ℂ) := by
    rw [Sp_C_p3]; push_cast; ring
  rw [h_eq, Complex.norm_real, Real.norm_eq_abs]
  have : 0 < 1 + α₀ + α₀ ^ 2 := by positivity
  exact abs_of_pos this

/-- `‖Q_3(z, α₀)‖ = ‖1 + z + α₀‖ ≤ 1 + 2α₀ + r₀` sur le disque. -/
theorem q_ub_on_disk (α₀ r₀ : ℝ) (hα : 0 < α₀) (z : ℂ)
    (hz : ‖z - ((α₀ : ℂ))‖ ≤ r₀) :
    ‖1 + z + ((α₀ : ℂ))‖ ≤ 1 + 2 * α₀ + r₀ := by
  have h_α : ‖((α₀ : ℂ))‖ = α₀ := by
    rw [show ((α₀ : ℂ)) = ((α₀ : ℝ) : ℂ) from rfl, Complex.norm_real,
        Real.norm_eq_abs]
    exact abs_of_pos hα
  have h_norm_z := norm_ub_of_disk α₀ r₀ hα z hz
  have h1 : ‖(1 : ℂ) + z + ((α₀ : ℂ))‖ ≤ ‖(1 : ℂ) + z‖ + ‖((α₀ : ℂ))‖ :=
    norm_add_le _ _
  have h2 : ‖(1 : ℂ) + z‖ ≤ 1 + ‖z‖ := by
    have := norm_add_le (1 : ℂ) z
    rw [norm_one] at this
    exact this
  rw [h_α] at h1
  linarith

end HelpersC

/-! ============================================================
  §57.1  Lower bound `‖Sp_C 3 z‖ ≥ 1 + (α₀-r₀) + (α₀-r₀)²`
============================================================ -/

section SpLBC

/-- Monotonie polynomiale : `(α₀-r₀) ≤ z.re ∧ z.re ≥ 1/2 ⟹
    1 + (α₀-r₀) + (α₀-r₀)² ≤ 1 + z.re + z.re²`. -/
private theorem sp_coord_mono (α₀ r₀ : ℝ) (hα_r : α₀ - r₀ ≥ 1 / 2) (z : ℂ)
    (hz_re : z.re ≥ α₀ - r₀) :
    1 + (α₀ - r₀) + (α₀ - r₀) ^ 2 ≤ 1 + z.re + z.re ^ 2 := by
  have h_diff : 0 ≤ z.re - (α₀ - r₀) := by linarith
  have h_sum : 0 ≤ z.re + (α₀ - r₀) := by linarith
  have h_sq : (α₀ - r₀) ^ 2 ≤ z.re ^ 2 := by
    have h_prod : 0 ≤ (z.re - (α₀ - r₀)) * (z.re + (α₀ - r₀)) :=
      mul_nonneg h_diff h_sum
    nlinarith [h_prod]
  linarith

/-- **`‖Sp_C 3 z‖ ≥ 1 + (α₀ − r₀) + (α₀ − r₀)²`** sur `B(α₀, r₀)`
    avec `α₀ − r₀ ≥ 1/2`. -/
theorem sp_lb_on_disk (α₀ r₀ : ℝ) (hα_r : α₀ - r₀ ≥ 1 / 2) (z : ℂ)
    (hz : ‖z - ((α₀ : ℂ))‖ ≤ r₀) :
    ‖Sp_C 3 z‖ ≥ 1 + (α₀ - r₀) + (α₀ - r₀) ^ 2 := by
  have hz_re : z.re ≥ α₀ - r₀ := re_ge_of_disk α₀ r₀ z hz
  have hz_re_half : z.re ≥ 1 / 2 := by linarith
  have h_lb_sq := Sp_C_p3_norm_sq_lower_bound z hz_re_half
  have h_mono := sp_coord_mono α₀ r₀ hα_r z hz_re
  have h_target_nn : 0 ≤ 1 + (α₀ - r₀) + (α₀ - r₀) ^ 2 := by
    have : 0 ≤ (α₀ - r₀) ^ 2 := sq_nonneg _
    linarith
  have h_sq_mono : (1 + (α₀ - r₀) + (α₀ - r₀) ^ 2) ^ 2
                 ≤ (1 + z.re + z.re ^ 2) ^ 2 :=
    pow_le_pow_left h_target_nn h_mono 2
  have h_combined : (1 + (α₀ - r₀) + (α₀ - r₀) ^ 2) ^ 2 ≤ ‖Sp_C 3 z‖ ^ 2 :=
    le_trans h_sq_mono h_lb_sq
  have h_norm_nn : 0 ≤ ‖Sp_C 3 z‖ := norm_nonneg _
  -- Pass to square roots: a² ≤ b² with a, b ≥ 0 ⟹ a ≤ b.
  have h_sqrt := Real.sqrt_le_sqrt h_combined
  rw [Real.sqrt_sq h_target_nn, Real.sqrt_sq h_norm_nn] at h_sqrt
  exact h_sqrt

end SpLBC

/-! ============================================================
  §57.3  Quantitative contraction bound
============================================================ -/

section QuantitativeC

/-- **★★★ Borne quantitative de contraction pour `pandrosion_h_C 3`
    sur un disque complexe.**

    Sur `B(α₀, r₀)` ⊂ ℂ avec paramètres satisfaisant :
      - `α₀ ≥ 3/4, 0 < r₀ ≤ α₀ − 1/2`,
      - `x > 1`, `α₀³ = 1/x`,
      - `‖z − α₀‖ ≤ r₀`,
    on a la borne multiplicative :

        `‖h(z) − α₀‖ · [x · (1 + (α₀-r₀) + (α₀-r₀)²) · (1 + α₀ + α₀²)]
         ≤ (x − 1) · (1 + 2α₀ + r₀) · ‖z − α₀‖`. -/
theorem pandrosion_h_C_p3_quantitative_bound
    (α₀ : ℝ) (hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_pos : 0 < r₀) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (x : ℝ) (hx : x > 1) (hα_fp : (α₀ : ℝ) ^ 3 = 1 / x)
    (z : ℂ) (hz_bound : ‖z - (α₀ : ℂ)‖ ≤ r₀) :
    ‖pandrosion_h_C (x : ℂ) 3 z - (α₀ : ℂ)‖ *
      (x * (1 + (α₀ - r₀) + (α₀ - r₀) ^ 2) * (1 + α₀ + α₀ ^ 2))
    ≤ (x - 1) * (1 + 2 * α₀ + r₀) * ‖z - (α₀ : ℂ)‖ := by
  -- Constantes de base.
  have hα_pos : (0 : ℝ) < α₀ := by linarith
  have hx_pos : (0 : ℝ) < x := by linarith
  have hα_r_diff : α₀ - r₀ ≥ 1 / 2 := by linarith
  have hz_re_half : z.re ≥ 1 / 2 := by
    have := re_ge_of_disk α₀ r₀ z hz_bound
    linarith
  -- Non-annulations.
  have hSα_pos : 0 < 1 + α₀ + α₀ ^ 2 := by positivity
  have h_Sα_norm : ‖Sp_C 3 ((α₀ : ℂ))‖ = 1 + α₀ + α₀ ^ 2 :=
    sp_α_norm α₀ hα_pos
  have h_Sα_ne : Sp_C 3 ((α₀ : ℂ)) ≠ 0 := by
    intro h; rw [h, norm_zero] at h_Sα_norm; linarith
  have h_Sz_ne : Sp_C 3 z ≠ 0 :=
    Sp_C_p3_ne_zero_of_re_ge_half z hz_re_half
  have hx_ne : ((x : ℂ)) ≠ 0 := by exact_mod_cast ne_of_gt hx_pos
  have hα_fp_C : ((α₀ : ℂ)) ^ 3 = 1 / ((x : ℂ)) := by
    have := hα_fp; push_cast; exact_mod_cast this
  -- Identité §56.3.
  have h_diff := pandrosion_h_C_p3_diff_identity ((x : ℂ)) ((α₀ : ℂ)) z
    hx_ne h_Sz_ne h_Sα_ne hα_fp_C
  -- Norme de la différence.
  have h_x_sub_norm : ‖((x : ℂ)) - 1‖ = x - 1 := by
    have h_eq : ((x : ℂ)) - 1 = (((x - 1 : ℝ)) : ℂ) := by push_cast; ring
    rw [h_eq, Complex.norm_real, Real.norm_eq_abs]
    exact abs_of_pos (by linarith : (0 : ℝ) < x - 1)
  have h_x_norm : ‖((x : ℂ))‖ = x := by
    rw [show ((x : ℂ)) = ((x : ℝ) : ℂ) from rfl, Complex.norm_real, Real.norm_eq_abs]
    exact abs_of_pos hx_pos
  have h_D_pos : 0 < ‖((x : ℂ))‖ * ‖Sp_C 3 z‖ * ‖Sp_C 3 ((α₀ : ℂ))‖ := by
    rw [h_x_norm, h_Sα_norm]
    have h_Sz_pos : 0 < ‖Sp_C 3 z‖ := by
      rcases (norm_nonneg (Sp_C 3 z)).lt_or_eq with h | h
      · exact h
      · exfalso; exact h_Sz_ne (norm_eq_zero.mp h.symm)
    positivity
  have h_D_ne : ‖((x : ℂ))‖ * ‖Sp_C 3 z‖ * ‖Sp_C 3 ((α₀ : ℂ))‖ ≠ 0 :=
    ne_of_gt h_D_pos
  have h_norm_diff :
      ‖pandrosion_h_C ((x : ℂ)) 3 z - ((α₀ : ℂ))‖ *
        (‖((x : ℂ))‖ * ‖Sp_C 3 z‖ * ‖Sp_C 3 ((α₀ : ℂ))‖) =
      ‖((x : ℂ)) - 1‖ * ‖1 + z + ((α₀ : ℂ))‖ * ‖z - ((α₀ : ℂ))‖ := by
    have h_diff_norm :
        ‖pandrosion_h_C ((x : ℂ)) 3 z - ((α₀ : ℂ))‖ =
        ‖((x : ℂ)) - 1‖ * ‖1 + z + ((α₀ : ℂ))‖ * ‖z - ((α₀ : ℂ))‖ /
          (‖((x : ℂ))‖ * ‖Sp_C 3 z‖ * ‖Sp_C 3 ((α₀ : ℂ))‖) := by
      rw [h_diff, norm_div]
      congr 1
      · rw [norm_mul, norm_mul]
      · rw [norm_mul, norm_mul]
    rw [h_diff_norm, div_mul_cancel₀ _ h_D_ne]
  -- Bornes sur les termes.
  have h_Q_ub : ‖1 + z + ((α₀ : ℂ))‖ ≤ 1 + 2 * α₀ + r₀ :=
    q_ub_on_disk α₀ r₀ hα_pos z hz_bound
  have h_Sz_lb : ‖Sp_C 3 z‖ ≥ 1 + (α₀ - r₀) + (α₀ - r₀) ^ 2 :=
    sp_lb_on_disk α₀ r₀ hα_r_diff z hz_bound
  have h_zm_nn : 0 ≤ ‖z - ((α₀ : ℂ))‖ := norm_nonneg _
  have h_hz_diff_nn : 0 ≤ ‖pandrosion_h_C ((x : ℂ)) 3 z - ((α₀ : ℂ))‖ :=
    norm_nonneg _
  -- Réécriture du LHS via h_norm_diff.
  -- On veut : LHS · (1+(α₀-r₀)+(α₀-r₀)²) ≤ ...
  -- Pas (norm_diff * bigprod). On compare directement.
  have h_lhs_via_diff :
      ‖pandrosion_h_C ((x : ℂ)) 3 z - ((α₀ : ℂ))‖ *
        (x * (1 + (α₀ - r₀) + (α₀ - r₀) ^ 2) * (1 + α₀ + α₀ ^ 2))
      ≤ ‖pandrosion_h_C ((x : ℂ)) 3 z - ((α₀ : ℂ))‖ *
        (‖((x : ℂ))‖ * ‖Sp_C 3 z‖ * ‖Sp_C 3 ((α₀ : ℂ))‖) := by
    apply mul_le_mul_of_nonneg_left _ h_hz_diff_nn
    rw [h_x_norm, h_Sα_norm]
    have h_mid : x * (1 + (α₀ - r₀) + (α₀ - r₀) ^ 2) ≤ x * ‖Sp_C 3 z‖ := by
      apply mul_le_mul_of_nonneg_left h_Sz_lb (le_of_lt hx_pos)
    have h_Sα_nn : 0 ≤ 1 + α₀ + α₀ ^ 2 := le_of_lt hSα_pos
    exact mul_le_mul_of_nonneg_right h_mid h_Sα_nn
  -- LHS ≤ ‖norm_diff expr via h_lhs_via_diff‖, puis = RHS_norm via h_norm_diff,
  -- puis borné par (x-1)·(1+2α₀+r₀)·‖z-α₀‖.
  have h_equality :
      ‖pandrosion_h_C ((x : ℂ)) 3 z - ((α₀ : ℂ))‖ *
        (‖((x : ℂ))‖ * ‖Sp_C 3 z‖ * ‖Sp_C 3 ((α₀ : ℂ))‖) ≤
      (x - 1) * (1 + 2 * α₀ + r₀) * ‖z - ((α₀ : ℂ))‖ := by
    rw [h_norm_diff, h_x_sub_norm]
    have h_num_bd :
        (x - 1) * ‖1 + z + ((α₀ : ℂ))‖ * ‖z - ((α₀ : ℂ))‖ ≤
        (x - 1) * (1 + 2 * α₀ + r₀) * ‖z - ((α₀ : ℂ))‖ := by
      apply mul_le_mul_of_nonneg_right _ h_zm_nn
      exact mul_le_mul_of_nonneg_left h_Q_ub (by linarith)
    linarith [h_num_bd]
  linarith [h_lhs_via_diff, h_equality]

end QuantitativeC

end Pandrosion
