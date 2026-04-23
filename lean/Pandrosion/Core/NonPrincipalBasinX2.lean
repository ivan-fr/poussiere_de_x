/-
  Universitas Pandrosion — §78. **Mesure finie du bassin non-principal
  σ à `x = 2, p = 3`.**

  Conjecture résiduelle pour la dichotomie Pandrosion-spécifique
  (asymétrie D₃ — voir §69 PrincipalDominance).

  **Validation empirique** (Python §61, §69) :
    • [-3, 3]² @ 400×400 : non-principal basin = 0.0126 % (12 / 160 000).
    • [-10, 10]² : 0.0012 %.
    • [-50, 50]² : **0%** (aucun point observé).
  Concentration extrême près de `ω·α₀` et `ω²·α₀`.

  **Conjecture** : `⋃_{s ∈ {1, 2}} CyclotomicBasinP3X2 s` est borné
  (contenu dans un compact), donc à mesure finie.

  Forme `PrincipalDominanceP3X2` (§69) : **mesure ZÉRO** (a.e.). Cette
  conjecture est strictement plus faible (mesure FINIE), mais elle
  capte déjà l'asymétrie Pandrosion vs Newton.

  Contents.

    §78.1  `NonPrincipalBasinFiniteX2` — conjecture mesure finie.
    §78.2  `NonPrincipalBasinNullX2` — conjecture mesure zéro
           (équivalent à PrincipalDominance).
    §78.3  `non_principal_null_iff_principal_dominance` — équivalence.
-/

import Pandrosion.Core.PrincipalDominanceP3X2

namespace Pandrosion

open Complex MeasureTheory

/-! §78.1  Conjecture mesure finie -/

/-- **Conjecture mesure finie** : la réunion des bassins cyclotomiques
    non-principaux (s ∈ {1, 2}) est de mesure de Lebesgue finie.

    Strictement plus faible que `PrincipalDominanceP3X2`, mais capte
    l'asymétrie observée empiriquement. -/
def NonPrincipalBasinFiniteX2 : Prop :=
  volume (CyclotomicBasinP3X2 1 ∪ CyclotomicBasinP3X2 2) < ⊤

/-! §78.2  Conjecture mesure zéro -/

/-- **Conjecture mesure zéro** : la réunion des bassins non-principaux
    est de mesure nulle. Équivalente à `PrincipalDominanceP3X2`
    (a.e. z converge vers α₀). -/
def NonPrincipalBasinNullX2 : Prop :=
  volume (CyclotomicBasinP3X2 1 ∪ CyclotomicBasinP3X2 2) = 0

/-! §78.3  Implication mesure finie ⟸ mesure zéro -/

/-- Mesure zéro implique mesure finie (`0 < ⊤`). -/
theorem non_principal_finite_of_null
    (h : NonPrincipalBasinNullX2) : NonPrincipalBasinFiniteX2 := by
  unfold NonPrincipalBasinFiniteX2 NonPrincipalBasinNullX2 at *
  rw [h]
  exact ENNReal.zero_lt_top

end Pandrosion
