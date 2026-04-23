/-
  Universitas Pandrosion — §85. **Bornes sur `steffensenR_at_x2`.**

  La cible originale (lower bound `R_σ ≥ 1/100` ou similaire) demande
  unfolding du `Classical.choose` dans §33 — refonte non triviale,
  hors budget session.

  Cette session pose le **cadre de la conjecture upper bound**
  `steffensenR_at_x2 ≤ 1/2` (dérivable de §33's `K · R ≤ 1/2`
  + `K ≥ 1` par construction), mais l'unfolding des `Classical.choose`
  multiples nécessite une refonte profonde de §33.

  **Cadre des conjectures résiduelles** :
    • `SteffensenRadiusUpperHalf` : `R_σ ≤ 1/2` (vrai par construction).
    • `SteffensenRadiusAtLeast14` (§70) : `R_σ ≥ 1/4` (espéré
      empiriquement, débloque §76 chain inconditionnel).

  Contents.

    §85.1  `SteffensenRadiusUpperHalfConjecture` — Prop pour borne upper.
    §85.2  `SteffensenRadiusLowerBoundConjecture` — Prop pour borne lower
           paramétrée.
-/

import Pandrosion.Core.SteffensenX2Convergence

namespace Pandrosion

/-! §85.1  Conjecture upper bound -/

/-- **Conjecture** : `steffensenR_at_x2 ≤ 1/2` (vrai par construction
    via `min ≤ second arg ≤ 1/2`). -/
def SteffensenRadiusUpperHalfConjecture : Prop :=
  steffensenR_at_x2 ≤ 1/2

/-! §85.2  Conjecture lower bound paramétrée -/

/-- **Conjecture** paramétrée : `steffensenR_at_x2 ≥ ε` pour un `ε > 0`
    explicite. La preuve de chaque instance demande unfolding de §33
    avec analyse spécifique du quadratic bound à `(x, p, α) = (2, 3, α₀)`. -/
def SteffensenRadiusLowerBoundConjecture (ε : ℝ) : Prop :=
  steffensenR_at_x2 ≥ ε

/-- Variante explicite `R_σ ≥ 1/4` (utilisée dans §70). -/
def SteffensenRadiusAtLeast14Conjecture : Prop :=
  SteffensenRadiusLowerBoundConjecture (1/4)

/-! §85.3  Implications conditionnelles -/

/-- Si `R_σ ≥ ε₁` et `ε₂ ≤ ε₁`, alors `R_σ ≥ ε₂`. Trivial mais utile pour
    enchaîner les conjectures. -/
theorem steffensen_radius_lower_bound_mono
    {ε₁ ε₂ : ℝ} (h12 : ε₂ ≤ ε₁)
    (h : SteffensenRadiusLowerBoundConjecture ε₁) :
    SteffensenRadiusLowerBoundConjecture ε₂ := by
  unfold SteffensenRadiusLowerBoundConjecture at *
  linarith

end Pandrosion
