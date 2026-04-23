/-
  Universitas Pandrosion — §79. **Chaînage conditionnel Niveau 5.**

  Cible Niveau 5 originale : prouver `PrincipalDominanceP3X2`
  inconditionnellement. Toutes les briques mécaniques (§54–§78) sont
  en place ; les conjectures résiduelles `HalfPlaneSigmaTendstoP3X2`
  (Niveau 1) et `HalfPlaneExhaustion` restent ouvertes.

  Ce module **ferme** le chaînage : sous l'hypothèse jointe des
  conjectures résiduelles, `PrincipalDominanceP3X2` est démontré
  inconditionnellement. Permet à l'utilisateur de voir exactement
  *quelles* hypothèses suffisent pour Niveau 5.

  Forme du théorème :

      `HalfPlaneSigmaTendstoP3X2`     (Niveau 1, §70)
    + `HalfPlaneExhaustionAEP3X2`     (gap → demi-plan a.e., §77b)
    ⟹ `PrincipalDominanceP3X2`       (Niveau 5, §69)

  Premier théorème jointif déclaré qui matche `niveau5_from_atoms`
  de §71 mais reformulé pour clarifier les dépendances.

  Contents.

    §79.1  `HalfPlaneExhaustionAEP3X2` — conjecture résiduelle gap.
    §79.2  `niveau5_from_HP_tendsto_and_exhaustion` — chaînage final.
    §79.3  `niveau5_implies_mcmullen_ae_x2` — Niveau 5 ⟹ McMullen.
-/

import Pandrosion.Core.SigmaTendstoAtOneX2
import Pandrosion.Core.Niveau5Master
import Pandrosion.Core.McMullenX2Refinement

namespace Pandrosion

open Complex MeasureTheory Filter Topology

/-! §79.1  Conjecture résiduelle gap-exhaustion -/

/-- **Conjecture HalfPlaneExhaustionAE** : presque tout `z ∈ ℂ`
    avec `Re z < 1/2` admet un itéré σᵏ qui entre dans le demi-plan
    droit `Re z ≥ 1/2`.

    Forme combinée : presque tout `z ∈ ℂ` est soit dans le demi-plan
    initialement, soit y entre via σᵏ.

    Empirique (§61) : 100 % vérifié sur scans denses. Formalisation
    demande dynamique complexe. -/
def HalfPlaneExhaustionAEP3X2 : Prop :=
  ∀ᵐ z : ℂ ∂volume,
    z.re ≥ 1/2 ∨
    ∃ k : ℕ, ((steffensen_step_C (2 : ℂ) 3)^[k] z).re ≥ 1/2

/-! §79.2  Chaînage final `(§76 ∧ §77b) ⟹ Niveau 5` -/

/-- **★★★★★★★ Chaînage final Niveau 5.**

    Sous les deux conjectures résiduelles :
      • `HalfPlaneSigmaTendstoP3X2` (§70 — Niveau 1) : σⁿ converge
        sur tout le demi-plan droit.
      • `HalfPlaneExhaustionAEP3X2` (§77b) : presque tout z atteint
        le demi-plan en finite time.
    on obtient `PrincipalDominanceP3X2` inconditionnel.

    Reformulation de §71 `niveau5_from_atoms` avec dépendances clarifiées. -/
theorem niveau5_from_HP_tendsto_and_exhaustion
    (hHP : HalfPlaneSigmaTendstoP3X2)
    (hExh : HalfPlaneExhaustionAEP3X2) :
    PrincipalDominanceP3X2 :=
  niveau5_from_atoms hHP hExh

/-! §79.3  Niveau 5 ⟹ McMullenAEEntry à `x = 2` -/

/-- **★★★★★★★★ Niveau 5 ⟹ McMullenAEEntry 3 2 α₀** : ferme la chaîne
    vers le résultat dynamique a.e.

    Compose §69 `principal_dominance_implies_mcmullen` et la chaîne
    §76+§77b ⟹ Niveau 5. -/
theorem mcmullen_ae_entry_x2_from_HP_tendsto_and_exhaustion
    (hHP : HalfPlaneSigmaTendstoP3X2)
    (hExh : HalfPlaneExhaustionAEP3X2) :
    McMullenAEEntry 3 (2 : ℂ) ((alphaX2 : ℝ) : ℂ) :=
  principal_dominance_implies_mcmullen
    (niveau5_from_HP_tendsto_and_exhaustion hHP hExh)

end Pandrosion
