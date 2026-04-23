/-
  Universitas Pandrosion — §69. **Bassins de Fatou de σ à `p = 3, x = 2`.**

  ⚠️ **HISTORIQUE — refondu après découverte empirique** ⚠️

  La version originale de ce module définissait `PrincipalDominanceP3X2`
  comme `∀ᵐ z, σⁿ(z) → α₀` (Niveau 5 strict). Les analyses Python
  (sigma_structural_analysis.py) ont **réfuté** cette conjecture :

    • Bassin de `ω·α₀` a mesure positive (~π·0.025² ≈ 0.002).
    • Bassin de `ω²·α₀` a mesure positive (par conjugaison complexe).
    • Donc `∀ᵐ z, σⁿ(z) → α₀` est **strictement faux**.

  La conjecture **correcte** est :
    • **Niveau 2** (`McMullenAEEntry 3 2 α₀`) : `∀ᵐ z, ∃ s, σⁿ(z) → γ_s`.
    • **Borne quantitative** (§78 reformulé) : `NonPrincipalBasin ⊆
      petits disques explicites autour de ω·α₀, ω²·α₀`.

  Ce module garde donc :
    • Définitions des bassins (`PrincipalBasinP3X2`, `CyclotomicBasinP3X2`).
    • Inclusions inconditionnelles (Steffensen disk ⊆ principal basin).
    • Lien avec McMullenAEEntry via §65 + §32.

  Contents.

    §69.1  `PrincipalBasinP3X2`, `CyclotomicBasinP3X2 s` — bassins.

    §69.2  `principal_basin_contains_steffensen_disk` — Steffensen
           disk ⊆ bassin principal (inconditionnel).

    §69.3  `h_principal_basin_contains_half_plane` — h-version pour
           demi-plan.
-/

import Pandrosion.Core.HalfPlaneTendstoX2
import Pandrosion.Core.McMullenX2Refinement

namespace Pandrosion

open Complex MeasureTheory Filter Topology

/-! §69.1  Bassins de Fatou par racine -/

/-- **Bassin principal de σ à x = 2** : l'ensemble des `z ∈ ℂ` dont
    l'itération σ converge vers la racine principale `α₀ = 2^{-1/3}`. -/
def PrincipalBasinP3X2 : Set ℂ :=
  {z : ℂ |
    Tendsto (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ))}

/-- **Bassin cyclotomique de σ à x = 2** : l'ensemble des `z ∈ ℂ`
    dont l'itération σ converge vers la racine cyclotomique `γ_s`. -/
def CyclotomicBasinP3X2 (s : Fin 3) : Set ℂ :=
  {z : ℂ |
    Tendsto (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s))}

/-- **`PrincipalBasinP3X2 = CyclotomicBasinP3X2 0`.** Le bassin principal
    est le bassin cyclotomique `s = 0` (par `cycAnchor α 3 0 = α`). -/
theorem principal_basin_eq_cyclotomic_zero :
    PrincipalBasinP3X2 = CyclotomicBasinP3X2 (0 : Fin 3) := by
  unfold PrincipalBasinP3X2 CyclotomicBasinP3X2
  ext z
  simp only [Set.mem_setOf_eq]
  rw [cycAnchor_p3_zero]

/-! §69.2  Steffensen disk ⊆ bassin principal (inconditionnel) -/

/-- **★★ Steffensen disk ⊆ bassin principal.**

    Conséquence directe de §65 `steffensen_step_C_p3_tendsto_at_x2`. -/
theorem principal_basin_contains_steffensen_disk :
    {z : ℂ | ‖z - ((alphaX2 : ℝ) : ℂ)‖ < steffensenR_at_x2}
      ⊆ PrincipalBasinP3X2 := by
  intro z hz
  exact steffensen_step_C_p3_tendsto_at_x2 z hz

/-! §69.3  Half-plane ⊆ bassin principal pour h -/

/-- **Le demi-plan `Re z ≥ 1/2` est dans le bassin principal pour `h`.**

    Version h-Tendsto de §67. *Note* : ce n'est PAS directement le
    bassin principal de σ (iteration différente), mais un lower bound
    conjecturel via le chaînage h → σ. -/
theorem h_principal_basin_contains_half_plane (z : ℂ) (hz : z.re ≥ 1/2) :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) :=
  pandrosion_h_C_p3_tendsto_half_plane_x2 z hz

end Pandrosion
