/-
  Universitas Pandrosion — §76. **σ-Tendsto à `z = 1` (chaînage
  conditionnel via §74 + §75 + §65).**

  Cible Niveau 1 originale (chaînage σ-Tendsto sur tout demi-plan)
  reste hors budget. Cette session ferme le chaînage **à un point
  concret** : `z = 1` est dans le demi-plan, on a déjà :

    • §75 : `steffensen_step_C 2 3 (1) = sigma_x2_closed (1) = 353/444`.
    • §74 : `|sigma_x2_closed (1) - α₀| < 1/4`.

  Si en plus `R_σ ≥ 1/4` (conjecture §70), alors `σ(1) ∈ B(α₀, R_σ)`,
  d'où σⁿ⁺¹(1) → α₀ par §65, donc σⁿ(1) → α₀ par shift.

  **Validation empirique** : Python confirme σⁿ(1) → α₀ rapidement
  (entry into B(α₀, 1/4) à n = 1).

  Contents.

    §76.1  `sigma_tendsto_at_one_conditional` — chaînage conditionnel
           sur `R_σ ≥ |σ_closed(1) - α₀|`.
-/

import Pandrosion.Core.SigmaEquivPointsX2
import Pandrosion.Core.SteffensenX2Convergence

namespace Pandrosion

open Complex Filter Topology

/-! §76.1  Chaînage conditionnel à `z = 1` -/

/-- **Chaînage conditionnel `σⁿ(1) → α₀` à `x = 2, p = 3`.**

    Hypothèse résiduelle : `‖σ(1) - α₀‖ < R_σ` (la position de σ(1)
    relativement au rayon Steffensen). Sous cette hypothèse, §65
    `steffensen_step_C_p3_tendsto_at_x2` s'applique depuis σ(1),
    puis shift d'itération conclut.

    *Note* : la valeur `‖σ(1) - α₀‖ = ‖353/444 - α₀‖` est connue
    (§74), bornée par `1/4`. Donc l'hypothèse `R_σ ≥ 1/4` (§70
    `SteffensenRadiusAtLeast14`) suffit. -/
theorem sigma_tendsto_at_one_conditional
    (h_R : ‖steffensen_step_C (2 : ℂ) 3 (1 : ℂ)
              - ((alphaX2 : ℝ) : ℂ)‖ < steffensenR_at_x2) :
    Tendsto (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n] (1 : ℂ))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  -- σⁿ(σ(1)) → α₀ par §65.
  have h_tendsto_shifted :=
    steffensen_step_C_p3_tendsto_at_x2
      (steffensen_step_C (2 : ℂ) 3 (1 : ℂ)) h_R
  -- σⁿ(σ(1)) = σⁿ⁺¹(1).
  have h_comp :
      (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n]
                       (steffensen_step_C (2 : ℂ) 3 (1 : ℂ)))
    = fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n + 1] (1 : ℂ) := by
    funext n; rw [Function.iterate_succ, Function.comp_apply]
  rw [h_comp] at h_tendsto_shifted
  exact (Filter.tendsto_add_atTop_iff_nat 1).mp h_tendsto_shifted

/-- **Chaînage conditionnel via `R_σ ≥ 1/4`.** Combine §70
    `SteffensenRadiusAtLeast14` avec §75 + §74 pour obtenir σⁿ(1) → α₀. -/
theorem sigma_tendsto_at_one_from_radius_at_least_quarter
    (h_R : steffensenR_at_x2 ≥ 1/4) :
    Tendsto (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n] (1 : ℂ))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply sigma_tendsto_at_one_conditional
  rw [steffensen_step_C_p3_at_x2_eq_at_one]
  -- ‖σ_closed(1) - α₀‖ < 1/4 ≤ R_σ.
  exact lt_of_lt_of_le sigma_x2_closed_in_basin_at_one h_R

end Pandrosion
