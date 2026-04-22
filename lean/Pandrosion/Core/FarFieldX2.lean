/-
  Universitas Pandrosion — §77. **Far-field convergence h à x = 2** :
  `|z| ≥ 2 ⟹ h^n(z) → α₀`.

  Extension de §67 (demi-plan droit `Re z ≥ 1/2`) au **complément de
  la boule `|z| < 2`**. Argument : pour `|z| ≥ 2`,

      `|Sp(z)| = |1 + z + z²| ≥ |z|² − |z| − 1 ≥ 4 − 2 − 1 = 1`

  donc `|1/(2·Sp(z))| ≤ 1/2`, d'où `h(z) ∈ B(1, 1/2) ⊂ {Re z ≥ 1/2}`.
  Une itération h amène dans le demi-plan, puis §67 conclut.

  **Couverture résultante** : `{|z| ≥ 2} ∪ {Re z ≥ 1/2}` est dans
  le h-principal-basin. Le complément `{|z| < 2 ∧ Re z < 1/2}` est
  un *compact borné* — significativement plus petit que la bande
  infinie initiale.

  Contents.

    §77.1  `Sp_norm_lower_bound_on_far_field` — `|Sp(z)| ≥ 1` pour `|z| ≥ 2`.
    §77.2  `h_maps_far_field_to_half_plane` — `|z| ≥ 2 ⟹ Re h(z) ≥ 1/2`.
    §77.3  `pandrosion_h_C_p3_tendsto_far_field_x2` — Tendsto sur far-field.
-/

import Pandrosion.Core.HalfPlaneTendstoX2

namespace Pandrosion

open Complex Filter Topology

/-! §77.1  `|Sp(z)| ≥ 1` sur `|z| ≥ 2` -/

/-- **`|Sp_C 3 z| ≥ 1` pour `|z| ≥ 2`.**

    Inégalité triangulaire inverse : `|1+z+z²| ≥ |z²| − |z| − 1`. -/
theorem Sp_norm_lower_bound_on_far_field (z : ℂ) (hz : ‖z‖ ≥ 2) :
    ‖Sp_C 3 z‖ ≥ 1 := by
  rw [Sp_C_p3]
  -- Triangle inégalité inverse, via preuve sur forme décomposée.
  have h_tri_dec :
      ‖(1 + z + z^2) - z - 1‖ ≤ ‖1 + z + z^2‖ + ‖z‖ + 1 := by
    calc ‖(1 + z + z^2) - z - 1‖
        ≤ ‖(1 + z + z^2) - z‖ + ‖(1 : ℂ)‖ := norm_sub_le _ _
      _ ≤ (‖1 + z + z^2‖ + ‖z‖) + ‖(1 : ℂ)‖ :=
          add_le_add_right (norm_sub_le _ _) _
      _ = ‖1 + z + z^2‖ + ‖z‖ + 1 := by rw [norm_one]
  have h_eq : ‖(z^2 : ℂ)‖ = ‖(1 + z + z^2) - z - 1‖ := by
    congr 1; ring
  have h_tri : ‖z^2‖ ≤ ‖1 + z + z^2‖ + ‖z‖ + 1 := by linarith
  have h_z_sq : ‖z^2‖ = ‖z‖^2 := norm_pow _ _
  -- ‖z‖² − ‖z‖ − 1 ≥ 1 à ‖z‖ = 2, via factorisation `(‖z‖-2)(‖z‖+1) ≥ 0`.
  have h_nn1 : 0 ≤ ‖z‖ - 2 := by linarith
  have h_nn2 : 0 ≤ ‖z‖ + 1 := by linarith [norm_nonneg z]
  have h_prod : 0 ≤ (‖z‖ - 2) * (‖z‖ + 1) := mul_nonneg h_nn1 h_nn2
  linarith [show ‖z‖^2 - ‖z‖ - 1 - 1 = (‖z‖ - 2) * (‖z‖ + 1) from by ring]

/-! §77.2  `|z| ≥ 2 ⟹ Re h(z) ≥ 1/2` -/

/-- **h envoie `|z| ≥ 2` dans le demi-plan droit.** -/
theorem h_maps_far_field_to_half_plane (z : ℂ) (hz : ‖z‖ ≥ 2) :
    (pandrosion_h_C (2 : ℂ) 3 z).re ≥ 1/2 := by
  have h_Sp_lb : ‖Sp_C 3 z‖ ≥ 1 := Sp_norm_lower_bound_on_far_field z hz
  have h_Sp_pos : 0 < ‖Sp_C 3 z‖ := by linarith
  have h2N : ‖(2 : ℂ)‖ = 2 := by
    rw [show ((2 : ℂ)) = (((2 : ℝ)) : ℂ) from by norm_num,
        Complex.norm_real, Real.norm_eq_abs]; norm_num
  -- `|1/(2·Sp(z))| ≤ 1/2`
  have h_inv_bound : ‖(1 : ℂ) / ((2 : ℂ) * Sp_C 3 z)‖ ≤ 1/2 := by
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

/-! §77.3  Far-field Tendsto -/

/-- **★★★★★ Far-field h-convergence à `x = 2, p = 3`.**

    Pour tout `z ∈ ℂ` avec `|z| ≥ 2`, `h^n(z) → α₀`.

    *Argument* : une itération h amène dans le demi-plan droit
    (§77.2), puis §67 conclut. -/
theorem pandrosion_h_C_p3_tendsto_far_field_x2
    (z : ℂ) (hz : ‖z‖ ≥ 2) :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  -- h(z) ∈ {Re ≥ 1/2}, donc par §67 Tendsto.
  have h_image_in_half_plane : (pandrosion_h_C (2 : ℂ) 3 z).re ≥ 1/2 :=
    h_maps_far_field_to_half_plane z hz
  have h_tail :=
    pandrosion_h_C_p3_tendsto_half_plane_x2
      (pandrosion_h_C (2 : ℂ) 3 z) h_image_in_half_plane
  -- Tendsto de σⁿ(σ z) = σ^(n+1)(z) ⟹ Tendsto de σⁿ(z).
  have h_comp :
      (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] (pandrosion_h_C (2 : ℂ) 3 z))
    = fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n + 1] z := by
    funext n; rw [Function.iterate_succ, Function.comp_apply]
  rw [h_comp] at h_tail
  exact (Filter.tendsto_add_atTop_iff_nat 1).mp h_tail

end Pandrosion
