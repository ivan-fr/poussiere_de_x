/-
  Universitas Pandrosion — §91. **Compte d'itérations concret pour
  `ε`-précision à `x = 2, p = 3`.**

  Instance numérique de §86.2 `pandrosion_h_C_p3_iter_count_at_x2_geometric_bound` :

    `(104/259)^4 · (1/4) ≤ 1/100`

  donc **4 itérations h suffisent** pour atteindre `ε = 1/100` depuis
  tout `z ∈ B(α₀, 1/4)`.

  Calcul vérification : `(104/259)^4 = 116985856/4493920401 ≈ 0.026`.
  `× 1/4 ≈ 0.0065 < 0.01` ✓.

  Contents.

    §91.1  `four_iter_suffices_for_one_hundredth` — 4 itérations
           ≤ 1/100.
    §91.2  `pandrosion_h_C_p3_four_iter_at_x2_within_one_hundredth` —
           h^4(z) à 1/100 près de α₀.
-/

import Pandrosion.Core.IterCountX2

namespace Pandrosion

open Complex

/-! §91.1  Inégalité numérique pivotale -/

/-- **`(104/259)^4 · (1/4) ≤ 1/100`** : 4 itérations suffisent pour
    précision `1/100` selon §86.2. -/
theorem four_iter_suffices_for_one_hundredth :
    (104/259 : ℝ)^4 * (1/4) ≤ 1/100 := by
  norm_num

/-! §91.2  4 itérations h donnent `1/100`-précision -/

/-- **★★ 4 itérations de `h` à `x = 2, p = 3` donnent une précision
    `1/100`** depuis tout `z ∈ B(α₀, 1/4)`. -/
theorem pandrosion_h_C_p3_four_iter_at_x2_within_one_hundredth
    (z : ℂ) (hz : ‖z - ((alphaX2 : ℝ) : ℂ)‖ ≤ (1/4 : ℝ)) :
    ‖(pandrosion_h_C (2 : ℂ) 3)^[4] z - ((alphaX2 : ℝ) : ℂ)‖ ≤ 1/100 :=
  pandrosion_h_C_p3_iter_count_at_x2_geometric_bound (1/100 : ℝ) z hz 4
    four_iter_suffices_for_one_hundredth 4 (le_refl 4)

end Pandrosion
