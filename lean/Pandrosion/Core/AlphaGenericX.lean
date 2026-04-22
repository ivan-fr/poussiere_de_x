/-
  Universitas Pandrosion — §83. **Ancre `α(x) := x^{-1/3}` pour
  `x > 1` arbitraire à `p = 3`.**

  Généralisation de §64.1 (`alphaX2 := 2^{-1/3}`) à un paramètre
  `x` variable. Données arithmétiques de base de la famille des
  ancres réelles positives :

    • `α(x) > 0` pour `x > 0`.
    • `α(x)^3 = 1/x`.
    • `α(x) < 1` pour `x > 1`.
    • `α(x) ≥ 1/2` pour `x ≤ 8` (car `(1/2)^3 = 1/8 = 1/x` à `x = 8`).
    • `α(x) > 3/4` pour `x ≤ (4/3)^3 = 64/27 ≈ 2.37`.

  **Portée du module** : *seulement* l'arithmétique de base de la famille
  `α(x)`, pas le théorème Banach généralisé. Le chaînage Banach §64
  full pour x variable demanderait paramétrer toutes les bornes
  polynomiales sur `x`, hors budget session.

  Contents.

    §83.1  `alphaX` — définition `x^{-1/3}` pour `x ∈ ℝ`.
    §83.2  `alphaX_pos`, `alphaX_pow_three_eq_inv`.
    §83.3  `alphaX_lt_one`, `alphaX_le_one_of_one_le`.
    §83.4  `alphaX_ge_half_of_le_eight`.
-/

import Pandrosion.Core.AlphaX2P4
import Mathlib.Analysis.SpecialFunctions.Pow.Real

namespace Pandrosion

/-! §83.1  Définition générique -/

/-- **Ancre réelle positive `α(x) := x^{-1/3}`** pour `x > 0`. -/
noncomputable def alphaX (x : ℝ) : ℝ := x ^ (-(1 : ℝ) / 3)

/-! §83.2  Identités basiques -/

theorem alphaX_pos {x : ℝ} (hx : 0 < x) : 0 < alphaX x :=
  Real.rpow_pos_of_pos hx _

theorem alphaX_pow_three_eq_inv {x : ℝ} (hx : 0 < x) :
    alphaX x ^ 3 = 1 / x := by
  unfold alphaX
  rw [← Real.rpow_nat_cast (x ^ (-(1 : ℝ) / 3)) 3,
      ← Real.rpow_mul (le_of_lt hx)]
  have h_exp : (-(1 : ℝ) / 3) * (3 : ℕ) = -1 := by push_cast; ring
  rw [h_exp, Real.rpow_neg_one, inv_eq_one_div]

/-! §83.3  Bornes supérieures -/

/-- `α(x) ≤ 1` pour `x ≥ 1`. -/
theorem alphaX_le_one_of_one_le {x : ℝ} (hx : 1 ≤ x) : alphaX x ≤ 1 := by
  by_contra h
  push_neg at h
  have h_pow_gt : (1 : ℝ) < alphaX x ^ 3 := by
    calc (1 : ℝ) = 1^3 := by norm_num
      _ < alphaX x ^ 3 := pow_lt_pow_left h (by norm_num) (by norm_num)
  rw [alphaX_pow_three_eq_inv (by linarith)] at h_pow_gt
  -- 1 < 1/x avec x ≥ 1 : contradiction (1/x ≤ 1).
  have h_inv_le_one : 1/x ≤ 1 := by
    rcases lt_or_eq_of_le hx with hx_lt | hx_eq
    · exact le_of_lt ((div_lt_one (by linarith : (0:ℝ) < x)).mpr hx_lt)
    · rw [← hx_eq]; norm_num
  linarith

/-- `α(x) < 1` pour `x > 1`. -/
theorem alphaX_lt_one {x : ℝ} (hx : 1 < x) : alphaX x < 1 := by
  rcases lt_or_eq_of_le (alphaX_le_one_of_one_le (le_of_lt hx)) with h_lt | h_eq
  · exact h_lt
  · exfalso
    have h_pow := alphaX_pow_three_eq_inv (by linarith : (0:ℝ) < x)
    rw [h_eq] at h_pow
    -- 1^3 = 1 ≠ 1/x quand x > 1.
    have h_eq_inv : (1 : ℝ) = 1/x := by
      linarith [show (1 : ℝ)^3 = 1 from by norm_num]
    have h_x_eq : x = 1 := by
      field_simp at h_eq_inv
      linarith
    linarith

/-! §83.4  Borne inférieure pour `x ≤ 8` -/

/-- `α(x) ≥ 1/2` pour `x ≤ 8`.

    *Preuve* : élever à la puissance 3.
    `α(x)^3 = 1/x ≥ 1/8 = (1/2)^3` ⟺ `α(x) ≥ 1/2`. -/
theorem alphaX_ge_half_of_le_eight {x : ℝ} (hx_pos : 0 < x) (hx_le : x ≤ 8) :
    alphaX x ≥ 1 / 2 := by
  have h_half_nn : (0 : ℝ) ≤ 1/2 := by norm_num
  have h_base_nn : (0 : ℝ) ≤ alphaX x := le_of_lt (alphaX_pos hx_pos)
  have h_pow_eq : alphaX x ^ 3 = 1 / x := alphaX_pow_three_eq_inv hx_pos
  have h_pow_bound : ((1 : ℝ) / 2) ^ 3 ≤ alphaX x ^ 3 := by
    rw [h_pow_eq]
    have h_inv_ge : (1 : ℝ) / x ≥ 1 / 8 :=
      one_div_le_one_div_of_le hx_pos hx_le
    linarith [show ((1 : ℝ) / 2) ^ 3 = 1/8 from by norm_num]
  exact (pow_le_pow_iff_left h_half_nn h_base_nn
    (by norm_num : 3 ≠ 0)).mp h_pow_bound

end Pandrosion
