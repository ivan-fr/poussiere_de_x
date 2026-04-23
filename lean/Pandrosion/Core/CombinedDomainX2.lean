/-
  Universitas Pandrosion — §90. **Tendsto h^n sur l'union demi-plan ∪ far-field
  à `x = 2, p = 3`.**

  Combine §67 (h-Tendsto sur `Re z ≥ 1/2`) + §77 (h-Tendsto sur
  `|z| ≥ 2`) en un seul énoncé : pour tout `z ∈ ℂ` dans
  `{Re z ≥ 1/2} ∪ {|z| ≥ 2}`, l'itération h converge vers `α₀`.

  **Couverture complémentaire** : ℂ \ (demi-plan ∪ far-field) = la
  région bornée gap `{Re z < 1/2 ∧ |z| < 2}`. Cette gap a mesure
  ≤ 4·(1/2 + 2)·... = bounded compact. Le complément a mesure infinie
  (le demi-plan domine).

  Contents.

    §90.1  `pandrosion_h_C_p3_tendsto_combined_x2` — Tendsto sur l'union.
    §90.2  `pandrosion_h_C_p3_combined_set_def` — set explicite.
-/

import Pandrosion.Core.FarFieldX2

namespace Pandrosion

open Complex Filter Topology

/-! §90.1  Tendsto sur l'union -/

/-- **★★ h-convergence sur l'union demi-plan ∪ far-field.**

    Pour tout `z ∈ ℂ` avec `Re z ≥ 1/2` OU `|z| ≥ 2`, on a
    `h^n(z) → α₀`. -/
theorem pandrosion_h_C_p3_tendsto_combined_x2
    (z : ℂ) (hz : z.re ≥ 1/2 ∨ ‖z‖ ≥ 2) :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  rcases hz with h_half | h_far
  · exact pandrosion_h_C_p3_tendsto_half_plane_x2 z h_half
  · exact pandrosion_h_C_p3_tendsto_far_field_x2 z h_far

/-! §90.2  Le set combiné -/

/-- **L'ensemble couvert** : `{z ∈ ℂ : Re z ≥ 1/2 ∨ |z| ≥ 2}`. -/
def CombinedDomainP3X2 : Set ℂ :=
  {z : ℂ | z.re ≥ 1/2 ∨ ‖z‖ ≥ 2}

/-- **Tout point du combined-domain est dans le h-bassin principal**
    (au sens : `h^n(z) → α₀`). -/
theorem combined_domain_in_h_principal_basin (z : ℂ) (hz : z ∈ CombinedDomainP3X2) :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) :=
  pandrosion_h_C_p3_tendsto_combined_x2 z hz

end Pandrosion
