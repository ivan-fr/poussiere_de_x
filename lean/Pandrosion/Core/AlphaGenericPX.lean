/-
  Universitas Pandrosion — §82. **Ancre générique `α(p, x) := x^{-1/p}`**.

  Unification de §81 (`alphaX2P4 = 2^{-1/4}`) et §83 (`alphaX x = x^{-1/3}`)
  en une définition unique paramétrée par `(p, x)`.

  Identités algébriques fondamentales :
    • `α(p, x) > 0` pour `x > 0`.
    • `α(p, x)^p = 1/x` (forme normalisée Pandrosion).
    • `α(p, 2) ≥ (1/2)^(1/p)` quand `p ≥ 1` (positif).

  Spécialisations :
    • `α(3, 2) = alphaX2` (§64).
    • `α(4, 2) = alphaX2P4` (§81).
    • `α(3, x) = alphaX x` (§83).

  **Portée** : seulement la définition unifiée + identités basiques.
  Les chaînages Banach paramétrés sur `(p, x)` resteraient à
  formaliser (demande analyse polynomiale paramétrée).

  Contents.

    §82.1  `alphaXp` — définition `x^{-1/p}` pour `p ≥ 1, x > 0`.
    §82.2  `alphaXp_pos`, `alphaXp_pow_eq_inv`.
    §82.3  Spécialisations vers `alphaX2`, `alphaX2P4`, `alphaX`.
-/

import Pandrosion.Core.AlphaGenericX
import Pandrosion.Core.BanachX2Concrete

namespace Pandrosion

/-! §82.1  Définition générique -/

/-- **Ancre générique `α(p, x) := x^{-1/p}`** pour `p : ℕ, x : ℝ`. -/
noncomputable def alphaXp (p : ℕ) (x : ℝ) : ℝ := x ^ (-(1 : ℝ) / p)

/-! §82.2  Propriétés basiques -/

theorem alphaXp_pos {p : ℕ} {x : ℝ} (hx : 0 < x) : 0 < alphaXp p x :=
  Real.rpow_pos_of_pos hx _

/-- **`α(p, x)^p = 1/x`** pour `x > 0, p ≥ 1`. -/
theorem alphaXp_pow_eq_inv {p : ℕ} (hp : 1 ≤ p) {x : ℝ} (hx : 0 < x) :
    alphaXp p x ^ p = 1 / x := by
  unfold alphaXp
  rw [← Real.rpow_nat_cast (x ^ (-(1 : ℝ) / p)) p,
      ← Real.rpow_mul (le_of_lt hx)]
  have h_p_pos : (0 : ℝ) < p := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hp
  have h_p_ne : (p : ℝ) ≠ 0 := ne_of_gt h_p_pos
  have h_exp : (-(1 : ℝ) / p) * (p : ℕ) = -1 := by
    push_cast
    field_simp
  rw [h_exp, Real.rpow_neg_one, inv_eq_one_div]

/-! §82.3  Spécialisations -/

/-- **`α(3, 2) = alphaX2`** (§64). -/
theorem alphaXp_three_two_eq_alphaX2 :
    alphaXp 3 2 = alphaX2 := by
  unfold alphaXp alphaX2
  norm_num

/-- **`α(4, 2) = alphaX2P4`** (§81). -/
theorem alphaXp_four_two_eq_alphaX2P4 :
    alphaXp 4 2 = alphaX2P4 := by
  unfold alphaXp alphaX2P4
  norm_num

/-- **`α(3, x) = alphaX x`** (§83). -/
theorem alphaXp_three_eq_alphaX (x : ℝ) :
    alphaXp 3 x = alphaX x := by
  unfold alphaXp alphaX
  norm_num

end Pandrosion
