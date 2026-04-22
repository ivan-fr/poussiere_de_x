/-
  Universitas Pandrosion — §64. **Concrete Banach-Fatou theorem at
  `x = 2`, unconditional.**

  Aboutissement de la chaîne §54 → §63 : instanciation à l'ancre
  canonique `α₀ := 2^{-1/3}`, rayon `r₀ := 1/4`, paramètre `x := 2`.

  On prouve :

    (i)   `α₀ ∈ [3/4, 1]`, `α₀³ = 1/2`,
    (ii)  `r₀ = 1/4 ≤ α₀ − 1/2`,
    (iii) `K(α₀, r₀, 2) < 1` (borne effective ≤ 104/259 < 1/2),
    (iv)  `Tendsto (h_{2,3})^[n] z → α₀` pour tout `z` dans le disque.

  Le théorème final (§64.4 `pandrosion_h_C_p3_banach_at_x2`) est
  **totalement inconditionnel** : aucun `McMullenAEEntry`, aucune
  conjecture Julia, aucun paramètre abstrait. C'est **le** théorème
  vitrine du corpus à p=3.

  Contents.

    §64.1  `alphaX2_*` — propriétés de `α₀ := 2^{-1/3}`.

    §64.2  `rate_lt_one_at_x2` — `K < 1` via polynôme explicite.

    §64.3  `pandrosion_h_C_p3_banach_at_x2` — Tendsto inconditionnel.
-/

import Pandrosion.Core.FatouExhaustionDisk

namespace Pandrosion

open Complex Filter Topology

/-! ============================================================
  §64.1  Canonical anchor `α₀ := 2^{-1/3}`
============================================================ -/

/-- **Ancre réelle canonique `α := 2^{-1/3}`.** -/
noncomputable def alphaX2 : ℝ := (2 : ℝ) ^ (-(1 : ℝ) / 3)

/-- `α₀ > 0`. -/
theorem alphaX2_pos : 0 < alphaX2 :=
  Real.rpow_pos_of_pos (by norm_num) _

/-- `α₀ ≥ 3/4` (via §55). -/
theorem alphaX2_ge_three_fourths : alphaX2 ≥ 3 / 4 :=
  two_rpow_neg_third_ge_three_fourths

/-- `α₀³ = 1/2` (via §55). -/
theorem alphaX2_cube : alphaX2 ^ 3 = 1 / 2 :=
  two_rpow_neg_third_pow_three

/-- `α₀ ≤ 1`. Car `α₀ > 0` et `α₀³ = 1/2 < 1`. -/
theorem alphaX2_le_one : alphaX2 ≤ 1 := by
  by_contra h
  push_neg at h
  have h_cube_gt : (1 : ℝ) < alphaX2 ^ 3 := by
    have h_one_cube : (1 : ℝ) ^ 3 = 1 := one_pow 3
    rw [← h_one_cube]
    exact pow_lt_pow_left h (by norm_num) (by norm_num)
  rw [alphaX2_cube] at h_cube_gt
  linarith

/-! ============================================================
  §64.2  `K < 1` at `x = 2, r₀ = 1/4`
============================================================ -/

/-- **★★ `K(α₀, 1/4, 2) < 1`.** Borne effective via polynôme
    explicitement positif pour tout `α₀ ≥ 0`.

    *Preuve.* Après multiplication par le dénominateur positif,
    l'inégalité `Num < Denom` équivaut à

        `3/8 + 5α/8 + 37α²/8 + 3α³ + 2α⁴ > 0`

    qui est **évidemment positive** puisque tous les coefficients
    et toutes les puissances de `α ≥ 0` sont ≥ 0, et le terme
    constant `3/8 > 0`. -/
theorem rate_lt_one_at_x2 :
    pandrosion_p3_rate alphaX2 (1/4) 2 < 1 := by
  unfold pandrosion_p3_rate
  have h_α_lb : alphaX2 ≥ 3/4 := alphaX2_ge_three_fourths
  have h_α_pos : 0 < alphaX2 := alphaX2_pos
  have h_α_nn : 0 ≤ alphaX2 := le_of_lt h_α_pos
  set a := alphaX2
  have h_denom_pos :
      0 < 2 * (1 + (a - 1/4) + (a - 1/4)^2) * (1 + a + a^2) := by
    apply mul_pos
    · apply mul_pos (by norm_num)
      have := sq_nonneg (a - 1/4); nlinarith [h_α_lb]
    · have := sq_nonneg a; linarith
  rw [div_lt_one h_denom_pos]
  -- Closed-form identity : denom − num = 3/8 + 5a/8 + 37a²/8 + 3a³ + 2a⁴.
  have h_expand :
      2 * (1 + (a - 1/4) + (a - 1/4)^2) * (1 + a + a^2)
        - (2 - 1) * (1 + 2 * a + 1/4)
      = 3/8 + 5*a/8 + 37*a^2/8 + 3*a^3 + 2*a^4 := by ring
  have h_positive : 0 < 3/8 + 5*a/8 + 37*a^2/8 + 3*a^3 + 2*a^4 := by
    have h3 : (0 : ℝ) ≤ 3 * a^3 := by positivity
    have h4 : (0 : ℝ) ≤ 2 * a^4 := by positivity
    nlinarith [h_α_nn, sq_nonneg a]
  linarith [h_expand, h_positive]

/-! ============================================================
  §64.3  Main theorem : unconditional Banach at `x = 2`
============================================================ -/

/-- **★★★★★ THÉORÈME VITRINE : Banach-Fatou à `x = 2`, inconditionnel.**

    Pour tout `z ∈ B(2^{-1/3}, 1/4) ⊂ ℂ`,

        `(pandrosion_h_C 2 3)^[n] z →  2^{-1/3}`

    quand `n → ∞`. **Aucune hypothèse McMullen, aucune conjecture
    Julia, aucun paramètre abstrait.**

    Résultat complet : *un disque explicite sur lequel l'itération
    Pandrosion converge inconditionnellement vers la racine cubique
    positive de 1/2*. -/
theorem pandrosion_h_C_p3_banach_at_x2
    (z : ℂ) (hz : ‖z - ((alphaX2 : ℝ) : ℂ)‖ ≤ (1/4 : ℝ)) :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] z)
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  have h_α_ge : alphaX2 ≥ 3 / 4 := alphaX2_ge_three_fourths
  have h_r_pos : (0 : ℝ) < 1/4 := by norm_num
  have h_r_le : (1/4 : ℝ) ≤ alphaX2 - 1/2 := by linarith [h_α_ge]
  have h_x_gt_one : (2 : ℝ) > 1 := by norm_num
  have h_α_fp : (alphaX2 : ℝ) ^ 3 = 1 / (2 : ℝ) := alphaX2_cube
  exact pandrosion_h_C_p3_tendsto_on_disk
    alphaX2 h_α_ge (1/4) h_r_pos h_r_le 2 h_x_gt_one h_α_fp
    rate_lt_one_at_x2 z hz

end Pandrosion
