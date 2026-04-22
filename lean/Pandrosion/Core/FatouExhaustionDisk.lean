/-
  Universitas Pandrosion — §63. **Fatou-exhaustion inconditionnelle
  sur le disque §58.**

  Premier théorème Fatou-exhaustion inconditionnel du corpus :
  pour tout `z₀ ∈ B(α₀, r₀)` (disque Banach de §58, avec `K < 1`),
  l'itération de `pandrosion_h_C x 3` converge vers `α₀`.

      `∀ z ∈ B(α₀, r₀), Tendsto (fun n ↦ h^n z) atTop (𝓝 α₀).`

  **Pourquoi "Fatou-exhaustion"** : un point `z` est dans le bassin
  de Fatou si l'itération converge ; ce théorème prouve que *tout*
  le disque §58 est dans le bassin de `α₀`. C'est la première
  instance **entièrement inconditionnelle** de cette propriété dans
  le corpus (par contraste avec `McMullenAEEntry`, toujours
  conditionnelle hors de `p = 2`).

  Chaîne : §58.4 `orbit_bound` + `tendsto_pow_atTop_nhds_zero` + squeeze.

  Contents.

    §63.1  `pandrosion_h_C_p3_tendsto_on_disk` — convergence sur le
           disque §58.

    §63.2  `pandrosion_h_C_p3_disk_in_half_plane` — le disque §58
           est contenu dans `{Re z ≥ 1/2}`, donc disjoint des zéros
           de `Sp_C 3`.
-/

import Pandrosion.Core.HalfPlaneContractionP3Ratio
import Pandrosion.Core.JuliaNullX2P3Preimages

namespace Pandrosion

open Complex MeasureTheory Filter Topology

/-! ============================================================
  §63.1  Tendsto on disk
============================================================ -/

/-- **★★★ Fatou-exhaustion inconditionnelle sur le disque §58.**

    Hypothèses : paramètres valides + `K < 1` (strict). Sous ces
    hypothèses, pour tout `z ∈ B(α₀, r₀)`,

        `Tendsto (fun n ↦ (pandrosion_h_C x 3)^[n] z) atTop (𝓝 α₀).`

    *Preuve.* §58.4 donne `‖h^n z − α₀‖ ≤ K^n · ‖z − α₀‖`.
    `K^n → 0` (standard), puis squeeze. -/
theorem pandrosion_h_C_p3_tendsto_on_disk
    (α₀ : ℝ) (hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_pos : 0 < r₀) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (x : ℝ) (hx : x > 1) (hα_fp : (α₀ : ℝ) ^ 3 = 1 / x)
    (hK_lt_one : pandrosion_p3_rate α₀ r₀ x < 1)
    (z : ℂ) (hz_bound : ‖z - (α₀ : ℂ)‖ ≤ r₀) :
    Tendsto (fun n : ℕ => (pandrosion_h_C (x : ℂ) 3)^[n] z)
            atTop (𝓝 ((α₀ : ℂ))) := by
  rw [tendsto_iff_norm_sub_tendsto_zero]
  have h_rate_nn : 0 ≤ pandrosion_p3_rate α₀ r₀ x :=
    pandrosion_p3_rate_nonneg α₀ hα₀_ge r₀ hr₀_pos hr₀_le x hx
  have h_K_le_one : pandrosion_p3_rate α₀ r₀ x ≤ 1 := le_of_lt hK_lt_one
  have h_orbit := pandrosion_h_C_p3_orbit_bound α₀ hα₀_ge r₀ hr₀_pos hr₀_le
    x hx hα_fp h_K_le_one z hz_bound
  have h_pow_tendsto :
      Tendsto (fun n : ℕ => pandrosion_p3_rate α₀ r₀ x ^ n) atTop (𝓝 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one h_rate_nn hK_lt_one
  have h_bound_tendsto :
      Tendsto (fun n : ℕ => pandrosion_p3_rate α₀ r₀ x ^ n * ‖z - (α₀ : ℂ)‖)
              atTop (𝓝 0) := by
    have h_mul := h_pow_tendsto.mul_const ‖z - (α₀ : ℂ)‖
    simpa using h_mul
  exact squeeze_zero (fun n => norm_nonneg _) h_orbit h_bound_tendsto

/-! ============================================================
  §63.2  Disk ⊂ {Re z ≥ 1/2}, disjoint from Sp-zeros
============================================================ -/

/-- **Disque §58 ⊂ `{Re z ≥ 1/2}`.** -/
theorem pandrosion_h_C_p3_disk_in_half_plane
    (α₀ : ℝ) (_hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (z : ℂ) (hz_bound : ‖z - (α₀ : ℂ)‖ ≤ r₀) :
    z.re ≥ 1 / 2 := by
  have h_re_sub : (z - ((α₀ : ℂ))).re = z.re - α₀ := by simp [Complex.sub_re]
  have h_abs_re : |z.re - α₀| ≤ ‖z - ((α₀ : ℂ))‖ := by
    rw [← h_re_sub]; exact Complex.abs_re_le_abs _
  have h_bound_re : |z.re - α₀| ≤ r₀ := le_trans h_abs_re hz_bound
  have h_lb : z.re - α₀ ≥ -r₀ := (abs_le.mp h_bound_re).1
  linarith

/-- **★★ Disque §58 disjoint des zéros de `Sp_C 3`.**

    Conséquence de §63.2 + §56.2 (`Sp_C_p3_ne_zero_of_re_ge_half`). -/
theorem pandrosion_h_C_p3_disk_disjoint_sp_zeros
    (α₀ : ℝ) (hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (z : ℂ) (hz_bound : ‖z - (α₀ : ℂ)‖ ≤ r₀) :
    Sp_C 3 z ≠ 0 :=
  Sp_C_p3_ne_zero_of_re_ge_half z
    (pandrosion_h_C_p3_disk_in_half_plane α₀ hα₀_ge r₀ hr₀_le z hz_bound)

end Pandrosion
