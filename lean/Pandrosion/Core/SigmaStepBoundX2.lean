/-
  Universitas Pandrosion — §73. **Borne `|σ_closed(z) − α₀| < 1/4`
  sur le demi-plan à x=2, p=3 — framework + cas trivial.**

  Cible Niveau 1 : prouver `SigmaStepIntoBasinX2` (§70 conjecture C1),
  réformulée sur `sigma_x2_closed` (§72).

  **Version atteinte dans cette session** : définition du Prop,
  vérification au point fixe `α₀`, et structure pour attaque future
  via identité pivot §72.

  **Argument visé (non-formalisé)** : `σ_closed(z) − α₀ = (N − α₀·4·BigPoly)
  / (4·BigPoly)`. Sur demi-plan `Re z ≥ 1/2`, `|BigPoly(z)|` est borné
  inférieurement (via §54 pour `|1+z+z²|` + analyse pour `|SmallPoly|`).
  Le numérateur `N − α₀·4·BigPoly` est un polynôme de degré 6 avec
  coefficients dépendant de `α₀` : majoration par triangle.

  **Validation empirique** (§70) : 100 % OK sur `[1/2, 100] × [-100, 100]`.

  **Formalisation complète différée** : le raisonnement polynomial
  asymptotique + borne finie sur `[1/2, 2]` demande un découpage
  fin dépassant le budget session. Cible pour session future.

  Contents.

    §73.1  `SigmaClosedStepIntoBasinX2` — conjecture sur `sigma_x2_closed`.
    §73.2  `sigma_closed_in_basin_at_alpha` — cas trivial `z = α₀`.
    §73.3  `SigmaClosedStepIntoBasinX2_implies_SigmaStepIntoBasinX2` —
           chaînage conditionnel vers §70 si §75 prouve l'équivalence.
-/

import Pandrosion.Core.SigmaClosedFormX2
import Pandrosion.Core.SigmaHalfPlaneP3X2

namespace Pandrosion

open Complex

/-! §73.1  Conjecture `sigma_x2_closed` sur demi-plan -/

/-- **Conjecture C1 reformulée sur `sigma_x2_closed`** : pour
    `z ∈ ℂ` avec `Re z ≥ 1/2`, `|σ_closed(z) − α₀| < 1/4`.

    Forme alternative de §70 `SigmaStepIntoBasinX2`, utile tant que
    §75 (équivalence `steffensen_step_C ↔ sigma_x2_closed`) reste
    ouverte. -/
def SigmaClosedStepIntoBasinX2 : Prop :=
  ∀ z : ℂ, z.re ≥ 1/2 →
    ‖sigma_x2_closed z - ((alphaX2 : ℝ) : ℂ)‖ < 1/4

/-! §73.2  Cas trivial : `z = α₀` (point fixe) -/

/-- **`|σ_closed(α₀) − α₀| = 0 < 1/4`.**

    Conséquence directe de §72 `sigma_x2_closed_at_alpha`. -/
theorem sigma_closed_in_basin_at_alpha :
    ‖sigma_x2_closed ((alphaX2 : ℝ) : ℂ) - ((alphaX2 : ℝ) : ℂ)‖ < 1/4 := by
  rw [sigma_x2_closed_at_alpha]
  rw [sub_self, norm_zero]
  norm_num

/-! §73.3  Chaînage conditionnel vers §70 -/

/-- **Équivalence conditionnelle `SigmaClosedStepIntoBasinX2 ↔
    SigmaStepIntoBasinX2`** si §75 (équivalence rationnelle) est
    établie.

    Actuellement déclarée comme implication conditionnelle sur
    l'hypothèse `h_eq : ∀ z admissible, steffensen_step_C (2:ℂ) 3 z
    = sigma_x2_closed z`. Quand §75 sera prouvé, cette implication
    devient automatique. -/
theorem sigma_closed_step_implies_step_into_basin
    (h_eq : ∀ z : ℂ, z.re ≥ 1/2 →
      steffensen_step_C (2 : ℂ) 3 z = sigma_x2_closed z)
    (h_closed : SigmaClosedStepIntoBasinX2) :
    SigmaStepIntoBasinX2 := by
  intro z hz
  rw [h_eq z hz]
  exact h_closed z hz

end Pandrosion
