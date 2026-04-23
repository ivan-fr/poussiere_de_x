/-
  Universitas Pandrosion — §74. **Bornes concrètes `|σ_closed(z) − α₀| < 1/4`
  à des points spécifiques du demi-plan.**

  Cible Niveau 1 : prouver `SigmaClosedStepIntoBasinX2` (§73). Le full
  bound polynomial sur tout `Re z ≥ 1/2` reste hors budget. Cette
  session ajoute des **vérifications concrètes** à des points nommés :

    • `z = 1` : σ_closed(1) = 353/444 ≈ 0.795.
    • `z = 2` : σ_closed(2) = 5571/6888 ≈ 0.809.

  À chaque point, l'écart à `α₀ = 2^{-1/3} ∈ [3/4, 1]` est borné par
  `1/4`. Démontre que C1 tient sur ces points concrets — building
  blocks pour une éventuelle preuve full.

  Contents.

    §74.1  `sigma_x2_closed_at_one` — `σ_closed(1) = 353/444`.
    §74.2  `sigma_x2_closed_in_basin_at_one` — `|σ_closed(1) - α₀| < 1/4`.
    §74.3  `sigma_x2_closed_at_two` — `σ_closed(2) = 5571/6888`.
    §74.4  `sigma_x2_closed_in_basin_at_two` — `|σ_closed(2) - α₀| < 1/4`.
-/

import Pandrosion.Core.SigmaStepBoundX2

namespace Pandrosion

open Complex

/-! §74.1  Valeur exacte à `z = 1` -/

theorem sigma_x2_closed_at_one :
    sigma_x2_closed (1 : ℂ) = 353 / 444 := by
  unfold sigma_x2_closed
  norm_num

/-! §74.2  Borne C1 à `z = 1` -/

/-- **`|σ_closed(1) − α₀| < 1/4`.**

    Argument : `σ_closed(1) = 353/444 ≈ 0.795` est dans `B(α₀, 1/4)`
    car `α₀ ∈ [3/4, 1]` et `353/444 ∈ (3/4 − 1/4, 1 + 1/4) =
    (1/2, 5/4)`. -/
theorem sigma_x2_closed_in_basin_at_one :
    ‖sigma_x2_closed (1 : ℂ) - ((alphaX2 : ℝ) : ℂ)‖ < 1/4 := by
  rw [sigma_x2_closed_at_one]
  -- Réduire à norme réelle.
  have h_sub_real :
      ((353 : ℂ) / 444) - ((alphaX2 : ℝ) : ℂ)
      = (((353 / 444 - alphaX2) : ℝ) : ℂ) := by
    push_cast; ring
  rw [h_sub_real, Complex.norm_real, Real.norm_eq_abs]
  -- |353/444 - alphaX2| < 1/4 via α₀ ∈ [3/4, 1].
  have h_α_lb : alphaX2 ≥ 3/4 := alphaX2_ge_three_fourths
  have h_α_ub : alphaX2 ≤ 1 := alphaX2_le_one
  -- 353/444 - 1 = -91/444 ≈ -0.205, 353/444 - 3/4 = 20/444 ≈ 0.045
  -- |·| ≤ max(0.205, 0.045) = 0.205 < 0.25.
  rw [abs_lt]
  refine ⟨?_, ?_⟩
  · -- -1/4 < 353/444 - alphaX2 ⟺ alphaX2 < 353/444 + 1/4 = 464/444
    linarith
  · -- 353/444 - alphaX2 < 1/4 ⟺ alphaX2 > 353/444 - 1/4 = 242/444
    linarith

/-! §74.3  Valeur exacte à `z = 2` -/

theorem sigma_x2_closed_at_two :
    sigma_x2_closed (2 : ℂ) = 5571 / 6888 := by
  unfold sigma_x2_closed
  norm_num

/-! §74.4  Borne C1 à `z = 2` -/

theorem sigma_x2_closed_in_basin_at_two :
    ‖sigma_x2_closed (2 : ℂ) - ((alphaX2 : ℝ) : ℂ)‖ < 1/4 := by
  rw [sigma_x2_closed_at_two]
  have h_sub_real :
      ((5571 : ℂ) / 6888) - ((alphaX2 : ℝ) : ℂ)
      = (((5571 / 6888 - alphaX2) : ℝ) : ℂ) := by
    push_cast; ring
  rw [h_sub_real, Complex.norm_real, Real.norm_eq_abs]
  have h_α_lb : alphaX2 ≥ 3/4 := alphaX2_ge_three_fourths
  have h_α_ub : alphaX2 ≤ 1 := alphaX2_le_one
  -- 5571/6888 - 1 = -1317/6888, 5571/6888 - 3/4 = 405/6888
  -- max abs ≤ 1317/6888 ≈ 0.191 < 0.25
  rw [abs_lt]
  refine ⟨?_, ?_⟩
  · linarith
  · linarith

end Pandrosion
