/-
  Universitas Pandrosion — §75. **Équivalence `steffensen_step_C 2 3
  ↔ sigma_x2_closed` à des points spécifiques.**

  Cible Niveau 1 : prouver l'équivalence formelle à `x = 2, p = 3`.
  Le full equivalence demande `field_simp` sur expressions de degré 18+
  (voir §72 — déféré). Cette session ajoute des **vérifications
  ponctuelles** : prouve `steffensen_step_C 2 3 z = sigma_x2_closed z`
  à `z = 1` par calcul numérique direct.

  Building block pour la preuve générale future, ou pour utilisations
  isolées (orbit constant à un point fixe spécifique).

  Contents.

    §75.1  `steffensen_step_C_p3_at_x2_eq_at_one` — `σ(1) = σ_closed(1)`.
-/

import Pandrosion.Core.SigmaConcretePointsX2

namespace Pandrosion

open Complex

/-! §75.1  Équivalence ponctuelle à `z = 1` -/

/-- **`steffensen_step_C 2 3 (1) = sigma_x2_closed (1) = 353/444`.**

    Calcul direct : `h(1) = 5/6`, `h(h(1)) = 73/91`,
    `denom = h(h(1)) − 2·h(1) + 1 = 37/273`,
    `σ(1) = 1 − (h(1) − 1)² / denom = 1 − (1/36) / (37/273) = 353/444`.

    Match avec la forme rationnelle `σ_closed(1) = 353/444` (§74.1). -/
theorem steffensen_step_C_p3_at_x2_eq_at_one :
    steffensen_step_C (2 : ℂ) 3 (1 : ℂ) = sigma_x2_closed (1 : ℂ) := by
  rw [sigma_x2_closed_at_one]
  unfold steffensen_step_C steffensen_denom_C pandrosion_h_C
  simp only [Sp_C_p3]
  -- All numeric, norm_num closes.
  norm_num

end Pandrosion
