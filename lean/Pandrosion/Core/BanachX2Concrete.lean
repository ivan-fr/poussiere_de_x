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

  Le théorème final (§64.3 `pandrosion_h_C_p3_banach_at_x2`) est
  **totalement inconditionnel** : aucun `McMullenAEEntry`, aucune
  conjecture Julia, aucun paramètre abstrait. C'est **le** théorème
  vitrine du corpus à p=3.

  **§64.4–§64.5 : complexité certifiée.** On renforce la Banach
  abstraite en un taux numérique explicite `K ≤ 104/259 < 1/2`,
  donnant la borne géométrique effective

      `‖h^n z − α₀‖ ≤ (104/259)^n · (1/4)`,

  soit **convergence linéaire certifiée** avec compte d'itérations
  calculable : `n ≥ log(1/(4ε)) / log(259/104)` suffit pour
  `ε`-précision.

  Contents.

    §64.1  `alphaX2_*` — propriétés de `α₀ := 2^{-1/3}`.

    §64.2  `rate_lt_one_at_x2` — `K < 1` via polynôme explicite.

    §64.3  `pandrosion_h_C_p3_banach_at_x2` — Tendsto inconditionnel.

    §64.4  `rate_le_104_over_259_at_x2` — borne effective `K ≤ 104/259`.

    §64.5  `pandrosion_h_C_p3_rate_bound_at_x2` — borne géométrique
           explicite `(104/259)^n · (1/4)` par itération.
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

/-! ============================================================
  §64.4  Effective bound `K ≤ 104/259`
============================================================ -/

/-- **★★★★ Borne effective `K(α₀, 1/4, 2) ≤ 104/259 < 1/2`.**

    *Preuve.* Avec `a = α₀ ∈ [3/4, 1]` :
      • `Num = 5/4 + 2a ≤ 13/4` (car `a ≤ 1`).
      • `Denom = 2·(1 + (a−1/4) + (a−1/4)²)·(1 + a + a²)` :
          - `1 + (a − 1/4) + (a − 1/4)² ≥ 1 + 1/2 + 1/4 = 7/4`,
          - `1 + a + a² ≥ 1 + 3/4 + 9/16 = 37/16`,
          - donc `Denom ≥ 2 · 7/4 · 37/16 = 259/32`.
      • `K = Num/Denom ≤ (13/4) / (259/32) = 104/259`. -/
theorem rate_le_104_over_259_at_x2 :
    pandrosion_p3_rate alphaX2 (1/4) 2 ≤ 104/259 := by
  unfold pandrosion_p3_rate
  have h_α_lb : alphaX2 ≥ 3/4 := alphaX2_ge_three_fourths
  have h_α_ub : alphaX2 ≤ 1 := alphaX2_le_one
  have h_α_pos : 0 < alphaX2 := alphaX2_pos
  set a := alphaX2
  -- Majoration du numérateur.
  have h_num_le : (2 - 1) * (1 + 2 * a + 1/4) ≤ 13/4 := by linarith
  -- Minorations des deux facteurs du dénominateur.
  have h_a_sq_lb : a^2 ≥ 9/16 := by nlinarith [sq_nonneg (a - 3/4)]
  have h_sub_sq_lb : (a - 1/4)^2 ≥ 1/4 := by nlinarith [sq_nonneg (a - 3/4)]
  have h_f1_ge : 1 + (a - 1/4) + (a - 1/4)^2 ≥ 7/4 := by linarith
  have h_f2_ge : 1 + a + a^2 ≥ 37/16 := by linarith
  have h_f1_nn : 0 ≤ 1 + (a - 1/4) + (a - 1/4)^2 := by linarith
  -- Produit des deux minorations.
  have h_prod_ge :
      (7/4 : ℝ) * (37/16) ≤ (1 + (a - 1/4) + (a - 1/4)^2) * (1 + a + a^2) :=
    mul_le_mul h_f1_ge h_f2_ge (by norm_num) h_f1_nn
  have h_denom_ge :
      2 * (1 + (a - 1/4) + (a - 1/4)^2) * (1 + a + a^2) ≥ 259/32 := by
    have h_eq : (7/4 : ℝ) * (37/16) = 259/64 := by norm_num
    nlinarith [h_prod_ge, h_eq]
  have h_denom_pos :
      0 < 2 * (1 + (a - 1/4) + (a - 1/4)^2) * (1 + a + a^2) := by
    have h_f2_pos : 0 < 1 + a + a^2 := by linarith
    have h_f1_pos : 0 < 1 + (a - 1/4) + (a - 1/4)^2 := by linarith
    exact mul_pos (mul_pos (by norm_num) h_f1_pos) h_f2_pos
  rw [div_le_iff h_denom_pos]
  calc (2 - 1) * (1 + 2 * a + 1/4)
      ≤ 13/4 := h_num_le
    _ = (104/259) * (259/32) := by norm_num
    _ ≤ (104/259) * (2 * (1 + (a - 1/4) + (a - 1/4)^2) * (1 + a + a^2)) :=
        mul_le_mul_of_nonneg_left h_denom_ge (by norm_num)

/-! ============================================================
  §64.5  Explicit geometric rate per iteration
============================================================ -/

/-- **★★★★★ Borne géométrique effective par itération.**

    Pour tout `z ∈ B(2^{-1/3}, 1/4)` et tout `n ∈ ℕ`,

        `‖(pandrosion_h_C 2 3)^[n] z − 2^{-1/3}‖ ≤ (104/259)^n · (1/4)`.

    **Convergence linéaire certifiée** à taux `104/259 ≈ 0.401 < 1/2`.
    Pour `ε`-précision, `n ≥ log(1/(4ε)) / log(259/104)` suffit. -/
theorem pandrosion_h_C_p3_rate_bound_at_x2
    (z : ℂ) (hz : ‖z - ((alphaX2 : ℝ) : ℂ)‖ ≤ (1/4 : ℝ)) (n : ℕ) :
    ‖(pandrosion_h_C (2 : ℂ) 3)^[n] z - ((alphaX2 : ℝ) : ℂ)‖
      ≤ (104/259 : ℝ)^n * (1/4) := by
  have h_α_ge : alphaX2 ≥ 3/4 := alphaX2_ge_three_fourths
  have h_r_pos : (0 : ℝ) < 1/4 := by norm_num
  have h_r_le : (1/4 : ℝ) ≤ alphaX2 - 1/2 := by linarith
  have h_rate_le : pandrosion_p3_rate alphaX2 (1/4) 2 ≤ 104/259 :=
    rate_le_104_over_259_at_x2
  have h_rate_nn : 0 ≤ pandrosion_p3_rate alphaX2 (1/4) 2 :=
    pandrosion_p3_rate_nonneg alphaX2 h_α_ge (1/4) h_r_pos h_r_le 2 (by norm_num)
  have h_rate_le_one : pandrosion_p3_rate alphaX2 (1/4) 2 ≤ 1 :=
    le_of_lt rate_lt_one_at_x2
  have h_orbit := pandrosion_h_C_p3_orbit_bound alphaX2 h_α_ge
    (1/4) h_r_pos h_r_le 2 (by norm_num) alphaX2_cube h_rate_le_one z hz n
  have h_104_nn : (0 : ℝ) ≤ 104/259 := by norm_num
  have h_pow_le : pandrosion_p3_rate alphaX2 (1/4) 2 ^ n ≤ (104/259 : ℝ)^n :=
    pow_le_pow_left h_rate_nn h_rate_le n
  have h_104_pow_nn : (0 : ℝ) ≤ (104/259 : ℝ)^n := pow_nonneg h_104_nn n
  calc ‖(pandrosion_h_C (2 : ℂ) 3)^[n] z - ((alphaX2 : ℝ) : ℂ)‖
      ≤ pandrosion_p3_rate alphaX2 (1/4) 2 ^ n * ‖z - ((alphaX2 : ℝ) : ℂ)‖ :=
        h_orbit
    _ ≤ (104/259 : ℝ)^n * ‖z - ((alphaX2 : ℝ) : ℂ)‖ :=
        mul_le_mul_of_nonneg_right h_pow_le (norm_nonneg _)
    _ ≤ (104/259 : ℝ)^n * (1/4) :=
        mul_le_mul_of_nonneg_left hz h_104_pow_nn

end Pandrosion
