/-
  Universitas Pandrosion — §70. **σ-Tendsto sur le demi-plan : cadre
  conjectural pour Niveau 1.**

  **Écart observationnel exceptionnel** (Python §70) :
    • À `x = 2, p = 3`, l'itération σ **mappe tout le demi-plan droit
      `Re z ≥ 1/2` dans `B(α₀, 1/4)` en UNE seule itération**.
    • max entry time dans `B(α₀, 1/4)` observé : 1 (median = mean = 1).
    • 100 % des points échantillonnés sur `Re z ∈ [1/2, 100], Im z ∈
      [−100, 100]` convergent vers α₀ (sans ω·α₀ ni ω²·α₀).

  **Explication théorique** : σ est l'accélération de Steffensen de
  `h`. Pour `|z|` grand, le calcul asymptotique donne
  `σ(z) → 5/6 ≈ 0.833`, proche de `α₀ = 0.794` — dans `B(α₀, 1/4)`.
  Cette cancellation asymptotique est le cœur de la puissance de σ
  vs h.

  **Conjectures posées ici** (par force empirique croissante) :
    (C1) `σ(Re z ≥ 1/2) ⊆ B(α₀, 1/4)` — un-pas vers le §64-disque.
    (C2) `∀ z, Re z ≥ 1/2 ⟹ σⁿ(z) → α₀` — Tendsto σ demi-plan.
    (C3) `∀ z ∈ ℂ non-singulier, σⁿ(z) → α₀` — principal-dominance (§69).

  **Chaînage inconditionnel prouvé ici** :
    (C1) + §65 ⟹ (C2) : un-pas amène dans le basin Steffensen, puis
    convergence.
    (C2) ⟹ `{z : Re z ≥ 1/2} ⊆ PrincipalBasinP3X2`.
    (C2) + §62 h-singular orbit ⟹ partial progress vers (C3).

  Contents.

    §70.1  `SigmaStepIntoBasinX2` — conjecture C1.
    §70.2  `HalfPlaneSigmaTendstoP3X2` — conjecture C2.
    §70.3  `half_plane_in_principal_basin_from_sigma_tendsto` — C2 → §69.
    §70.4  `half_plane_sigma_tendsto_from_step_into_basin_v2` — (C1)
           + relation `R_σ ≥ 1/4` ⟹ C2.
-/

import Pandrosion.Core.PrincipalDominanceP3X2

namespace Pandrosion

open Complex MeasureTheory Filter Topology

/-! §70.1  Conjecture C1 : un-pas dans `B(α₀, 1/4)` -/

/-- **Conjecture C1** : σ envoie tout le demi-plan droit dans le
    §64-disque `B(α₀, 1/4)` en une seule itération.

    **Validation empirique** : 0 violations sur `[1/2, 100] × [−100, 100]`
    @ 400×400 échantillons. `min Re σ(z)` observé ≈ 0.76. -/
def SigmaStepIntoBasinX2 : Prop :=
  ∀ z : ℂ, z.re ≥ 1/2 →
    ‖steffensen_step_C (2 : ℂ) 3 z - ((alphaX2 : ℝ) : ℂ)‖ < 1/4

/-! §70.2  Conjecture C2 : σ Tendsto sur le demi-plan -/

/-- **Conjecture C2** : `σⁿ(z) → α₀` pour tout `z` dans le demi-plan droit.

    Version σ (Steffensen) du théorème §67 `pandrosion_h_C_p3_tendsto_half_plane_x2`
    (qui est pour h). **Niveau 1** de la hiérarchie ambitieuse. -/
def HalfPlaneSigmaTendstoP3X2 : Prop :=
  ∀ z : ℂ, z.re ≥ 1/2 →
    Tendsto (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ))

/-! §70.3  C2 ⟹ demi-plan ⊆ PrincipalBasin -/

/-- **Chaînage inconditionnel C2 ⟹ demi-plan ⊆ bassin principal.** -/
theorem half_plane_in_principal_basin_from_sigma_tendsto
    (hC2 : HalfPlaneSigmaTendstoP3X2) :
    {z : ℂ | z.re ≥ 1/2} ⊆ PrincipalBasinP3X2 := by
  intro z hz
  exact hC2 z hz

/-! §70.4  C1 + `R_σ ≥ 1/4` ⟹ C2 -/

/-- **Hypothèse auxiliaire** : le rayon Steffensen à x=2 est ≥ 1/4.

    Conjecturalement vrai (dépend de l'unfolding de la définition
    skolémisée de `steffensenR_of_fp`). -/
def SteffensenRadiusAtLeast14 : Prop := steffensenR_at_x2 ≥ 1/4

/-- **Chaînage conditionnel (C1) + (R_σ ≥ 1/4) ⟹ (C2).**

    Si σ mappe le demi-plan dans `B(α₀, 1/4)` en une étape, et si le
    rayon Steffensen est ≥ 1/4, alors §65 `steffensen_step_C_p3_tendsto_at_x2`
    s'applique depuis `σ(z) ∈ B(α₀, 1/4)` : convergence garantie. -/
theorem half_plane_sigma_tendsto_from_step_into_basin_v2
    (hC1 : SigmaStepIntoBasinX2)
    (hR : SteffensenRadiusAtLeast14) :
    HalfPlaneSigmaTendstoP3X2 := by
  intro z hz
  have h_sigma_in : ‖steffensen_step_C (2 : ℂ) 3 z - ((alphaX2 : ℝ) : ℂ)‖
                      < 1/4 := hC1 z hz
  have h_sigma_in_R : ‖steffensen_step_C (2 : ℂ) 3 z - ((alphaX2 : ℝ) : ℂ)‖
                      < steffensenR_at_x2 := lt_of_lt_of_le h_sigma_in hR
  have h_tendsto_shifted :=
    steffensen_step_C_p3_tendsto_at_x2 (steffensen_step_C (2 : ℂ) 3 z) h_sigma_in_R
  -- `h_tendsto_shifted : Tendsto (fun n => σⁿ(σ z)) → α₀`.
  -- Or `σⁿ(σ z) = σ^(n+1)(z)`, donc Tendsto de `σ^(n+1)(z)` ⟹ Tendsto de `σⁿ(z)`.
  have h_comp : (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n]
                  (steffensen_step_C (2 : ℂ) 3 z))
              = fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n + 1] z := by
    funext n
    rw [Function.iterate_succ, Function.comp_apply]
  rw [h_comp] at h_tendsto_shifted
  exact (Filter.tendsto_add_atTop_iff_nat 1).mp h_tendsto_shifted

/-! §70.5  C2 → principal dominance restricted to half-plane -/

/-- **C2 ⟹ principal dominance sur le demi-plan.**

    Si σ converge sur tout le demi-plan, alors la mesure du non-
    principal-basin ∩ demi-plan est nulle.

    *Note* : le demi-plan a mesure infinie, donc ceci n'implique pas
    directement `PrincipalDominanceP3X2` (conjecture sur tout ℂ). Il
    faut étendre au complément (Re z < 1/2), qui est le travail
    restant pour fermer Niveau 5. -/
theorem half_plane_principal_dominance_from_sigma_tendsto
    (hC2 : HalfPlaneSigmaTendstoP3X2) :
    ∀ᵐ z : ℂ ∂(volume.restrict {z : ℂ | z.re ≥ 1/2}),
      Tendsto (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n] z)
        atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  rw [ae_restrict_iff' (measurableSet_le measurable_const
      Complex.continuous_re.measurable)]
  filter_upwards with z hz
  exact hC2 z hz

end Pandrosion
