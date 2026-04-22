/-
  Universitas Pandrosion — §81. **Ancre `α₀ := 2^{-1/4}` à `x = 2, p = 4`.**

  Analogue de §55+§64.1 (qui établissent `α := 2^{-1/3}` pour `p = 3`)
  pour le cas `p = 4`. Données arithmétiques de base :

    • `α₀^4 = 1/2` (forme normalisée Pandrosion `α^p = 1/x`).
    • `α₀ ≥ 3/4` (car `(3/4)^4 = 81/256 ≤ 128/256 = 1/2 = α₀^4`).
    • `α₀ ≤ 1`.

  **Portée du module** : *seulement* l'arithmétique de base de `α₀`,
  pas le théorème Banach complet (analogue §64.3). Le chaînage
  Banach pour p=4 demanderait :
    – Identité §56.3 généralisée (factor `Q_4 = 1+z+α+z²+zα+α²`).
    – Borne polynomiale §80 (déjà acquise pour Re z ≥ 1, pas 1/2).
    – Vérification `K(p=4) < 1` numérique.

  Contents.

    §81.1  `alphaX2P4` — définition `2^{-1/4}`.
    §81.2  `alphaX2P4_pos`, `alphaX2P4_pow_four_eq_half`.
    §81.3  `alphaX2P4_ge_three_fourths` — borne 3/4 ≤ α₀.
    §81.4  `alphaX2P4_le_one` — α₀ ≤ 1.
-/

import Pandrosion.Core.HalfPlaneContractionP4
import Mathlib.Analysis.SpecialFunctions.Pow.Real

namespace Pandrosion

/-! §81.1  Définition de l'ancre p=4 -/

/-- **`α₀ := 2^{-1/4}`**, racine quatrième positive de `1/2`. -/
noncomputable def alphaX2P4 : ℝ := (2 : ℝ) ^ (-(1 : ℝ) / 4)

/-! §81.2  Identité `α₀^4 = 1/2` et positivité -/

theorem alphaX2P4_pos : 0 < alphaX2P4 :=
  Real.rpow_pos_of_pos (by norm_num) _

theorem alphaX2P4_pow_four_eq_half :
    alphaX2P4 ^ 4 = 1 / 2 := by
  unfold alphaX2P4
  have h2_pos : (0 : ℝ) < 2 := by norm_num
  rw [← Real.rpow_nat_cast ((2 : ℝ) ^ (-(1 : ℝ) / 4)) 4,
      ← Real.rpow_mul (le_of_lt h2_pos)]
  have h_exp : (-(1 : ℝ) / 4) * (4 : ℕ) = -1 := by push_cast; ring
  rw [h_exp, Real.rpow_neg_one]
  norm_num

/-! §81.3  Borne `α₀ ≥ 3/4` -/

/-- **`2^{-1/4} ≥ 3/4`.**

    *Preuve* : élever à la puissance 4.
    `(2^{-1/4})^4 = 1/2 = 128/256`,
    `(3/4)^4 = 81/256`,
    et `128/256 ≥ 81/256`. Monotonie de `^4` sur ℝ₊ transfère. -/
theorem alphaX2P4_ge_three_fourths : alphaX2P4 ≥ 3 / 4 := by
  have h_three_four_nn : (0 : ℝ) ≤ 3 / 4 := by norm_num
  have h_base_nn : (0 : ℝ) ≤ alphaX2P4 := le_of_lt alphaX2P4_pos
  have h_pow_eq : alphaX2P4 ^ 4 = 1 / 2 := alphaX2P4_pow_four_eq_half
  have h_pow_bound : ((3 : ℝ) / 4) ^ 4 ≤ alphaX2P4 ^ 4 := by
    rw [h_pow_eq]; norm_num
  exact (pow_le_pow_iff_left h_three_four_nn h_base_nn
    (by norm_num : 4 ≠ 0)).mp h_pow_bound

/-! §81.4  Borne `α₀ ≤ 1` -/

theorem alphaX2P4_le_one : alphaX2P4 ≤ 1 := by
  by_contra h
  push_neg at h
  have h_pow_gt : (1 : ℝ) < alphaX2P4 ^ 4 := by
    calc (1 : ℝ) = 1^4 := by norm_num
      _ < alphaX2P4 ^ 4 := pow_lt_pow_left h (by norm_num) (by norm_num)
  rw [alphaX2P4_pow_four_eq_half] at h_pow_gt
  linarith

end Pandrosion
