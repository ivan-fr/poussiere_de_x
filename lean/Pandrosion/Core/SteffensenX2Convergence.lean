/-
  Universitas Pandrosion — §65. **σ-convergence inconditionnelle
  sur le disque Steffensen à `x = 2`.**

  Passerelle décisive : le corpus §58–§64 établit la convergence de
  `pandrosion_h_C 2 3` (itération linéaire). Ici on la transfère à
  **`steffensen_step_C 2 3`** (itération σ), qui est la vraie cible
  du papier McMullen-style.

  Combinaison :
    §33.2 `steffensen_explicit_super_attractive_rate` → taux quadratique
    local pour σ au voisinage d'un point fixe Pandrosion ;
    §34.1 `quadratic_loglog_from_basin` → convergence loglog à ε-près
    depuis l'intérieur du bassin explicite.

  Instanciation à `x = 2, α₀ = 2^{-1/3}` :
    • `α₀` est non nul ;
    • `Sp_C 3 α₀ = 1 + α₀ + α₀²` est non nul (réel positif) ;
    • `α₀³ = 1/2 = 1/x` (§55).

  Conséquence : `Tendsto (σ^[n] z) → α₀` pour tout `z` dans le disque
  Steffensen `B(α₀, R_σ)` où `R_σ := steffensenR_of_fp 2 3 α₀` est
  positif et concret (skolemisé).

  **Premier théorème σ-convergence inconditionnel du corpus.**

  Contents.

    §65.1  `alphaX2_C_*` — coercions en ℂ de l'ancre `α₀`.

    §65.2  `steffensenR_at_x2` — rayon Steffensen instancié.

    §65.3  `steffensen_step_C_p3_tendsto_at_x2` — Tendsto σ^n z → α₀.
-/

import Pandrosion.Core.BanachX2Concrete
import Pandrosion.Core.SteffensenGlobalLoglog

namespace Pandrosion

open Complex Filter Topology

/-! ============================================================
  §65.1  Complex coercions
============================================================ -/

/-- `α = 2^{-1/3} ≠ 0` en ℂ. -/
theorem alphaX2_C_ne_zero : ((alphaX2 : ℝ) : ℂ) ≠ 0 := by
  rw [Ne, Complex.ofReal_eq_zero]
  exact ne_of_gt alphaX2_pos

/-- `Sp_C 3 α ≠ 0` at `α = 2^{-1/3}` (car `1 + α + α² > 0` pour α réel > 0). -/
theorem Sp_C_p3_alphaX2_ne_zero : Sp_C 3 ((alphaX2 : ℝ) : ℂ) ≠ 0 := by
  rw [Sp_C_p3]
  have h_real :
      (1 : ℂ) + ((alphaX2 : ℝ) : ℂ) + ((alphaX2 : ℝ) : ℂ) ^ 2
        = (((1 + alphaX2 + alphaX2 ^ 2) : ℝ) : ℂ) := by
    push_cast; ring
  rw [h_real, Ne, Complex.ofReal_eq_zero]
  have h_pos : 0 < 1 + alphaX2 + alphaX2 ^ 2 := by
    have := sq_nonneg alphaX2
    linarith [alphaX2_pos]
  linarith

/-- `α³ = 1/2` en ℂ. -/
theorem alphaX2_cube_C : ((alphaX2 : ℝ) : ℂ) ^ 3 = 1 / (2 : ℂ) := by
  have h := alphaX2_cube
  have h_cast : ((alphaX2 : ℝ) : ℂ) ^ 3 = (((alphaX2 ^ 3) : ℝ) : ℂ) := by
    push_cast; ring
  rw [h_cast, h]
  push_cast; ring

/-! ============================================================
  §65.2  Instantiated Steffensen radius
============================================================ -/

/-- Rayon Steffensen positif à `x = 2`, `α₀ = 2^{-1/3}`. -/
noncomputable def steffensenR_at_x2 : ℝ :=
  steffensenR_of_fp (2 : ℂ) (by norm_num) 3 (by norm_num)
    ((alphaX2 : ℝ) : ℂ) alphaX2_C_ne_zero
    Sp_C_p3_alphaX2_ne_zero alphaX2_cube_C

/-- **`steffensenR_at_x2 > 0`**. -/
theorem steffensenR_at_x2_pos : 0 < steffensenR_at_x2 :=
  steffensenR_of_fp_pos (2 : ℂ) (by norm_num) 3 (by norm_num)
    ((alphaX2 : ℝ) : ℂ) alphaX2_C_ne_zero
    Sp_C_p3_alphaX2_ne_zero alphaX2_cube_C

/-! ============================================================
  §65.3  Main theorem: σ-Tendsto at x=2
============================================================ -/

/-- **★★★★★ σ-convergence inconditionnelle à `x = 2` sur le disque
    Steffensen.**

    Pour tout `z ∈ B(2^{-1/3}, R_σ)` où `R_σ := steffensenR_at_x2 > 0`,

        `(steffensen_step_C 2 3)^[n] z  →  2^{-1/3}`

    quand `n → ∞`. **Aucune hypothèse McMullen, aucune conjecture
    Julia.** Premier théorème σ-convergence inconditionnel du corpus.

    *Preuve.* §34.1 `quadratic_loglog_from_basin` donne la convergence
    à ε-près depuis le disque. Déballage en Tendsto via
    `Metric.tendsto_atTop`. -/
theorem steffensen_step_C_p3_tendsto_at_x2
    (z : ℂ) (hz : ‖z - ((alphaX2 : ℝ) : ℂ)‖ < steffensenR_at_x2) :
    Tendsto (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  rw [Metric.tendsto_atTop]
  intro ε hε_pos
  have h_half : (0 : ℝ) < ε / 2 := by linarith
  obtain ⟨N, hN⟩ :=
    quadratic_loglog_from_basin (2 : ℂ) (by norm_num) 3 (by norm_num)
      ((alphaX2 : ℝ) : ℂ) alphaX2_C_ne_zero Sp_C_p3_alphaX2_ne_zero
      alphaX2_cube_C z hz (ε / 2) h_half
  refine ⟨N, fun n hn => ?_⟩
  rw [dist_eq_norm]
  have h_bound := hN n hn
  linarith

end Pandrosion
