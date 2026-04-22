/-
  Universitas Pandrosion — §58. **Banach-style ratio form for
  `pandrosion_h_C x 3` on a disk.**

  Déroule §57 en une chaîne Banach-Cacciopoli :

      ratio → invariance du disque → borne d'orbite `K^n`.

  Constante explicite :

      `K(α₀, r₀, x) := ((x − 1) · (1 + 2α₀ + r₀)) /
                         (x · (1 + (α₀ − r₀) + (α₀ − r₀)²) · (1 + α₀ + α₀²))`.

  Sous une hypothèse `K < 1`, `h` contracte strictement le disque
  complexe `B(α₀, r₀)` et l'itération converge géométriquement.

  Contents.

    §58.1  `pandrosion_h_C_p3_ratio_bound` — `‖h(z)−α₀‖ ≤ K·‖z−α₀‖`.

    §58.2  `pandrosion_h_C_p3_banach_contraction` — existentiel
           empaqueté `∃ K, 0 ≤ K ∧ K < 1 ∧ …`.

    §58.3  `pandrosion_h_C_p3_disk_invariance` — `h(B(α₀,r₀)) ⊆ B(α₀,r₀)`
           sous `K ≤ 1`.

    §58.4  `pandrosion_h_C_p3_orbit_bound` — `‖hⁿ(z)−α₀‖ ≤ Kⁿ·‖z−α₀‖`.
-/

import Pandrosion.Core.HalfPlaneContractionP3Final

namespace Pandrosion

open Complex

/-- Constante `K` de contraction, notation abrégée. -/
noncomputable def pandrosion_p3_rate (α₀ r₀ x : ℝ) : ℝ :=
  ((x - 1) * (1 + 2 * α₀ + r₀)) /
    (x * (1 + (α₀ - r₀) + (α₀ - r₀) ^ 2) * (1 + α₀ + α₀ ^ 2))

/-! ============================================================
  §58.1  Ratio form
============================================================ -/

/-- **★★ Ratio form of the p=3 contraction on a disk.** -/
theorem pandrosion_h_C_p3_ratio_bound
    (α₀ : ℝ) (hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_pos : 0 < r₀) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (x : ℝ) (hx : x > 1) (hα_fp : (α₀ : ℝ) ^ 3 = 1 / x)
    (z : ℂ) (hz_bound : ‖z - (α₀ : ℂ)‖ ≤ r₀) :
    ‖pandrosion_h_C (x : ℂ) 3 z - (α₀ : ℂ)‖
      ≤ pandrosion_p3_rate α₀ r₀ x * ‖z - (α₀ : ℂ)‖ := by
  unfold pandrosion_p3_rate
  have h_B_pos : 0 < x * (1 + (α₀ - r₀) + (α₀ - r₀) ^ 2) * (1 + α₀ + α₀ ^ 2) := by
    have h1 : 0 < x := by linarith
    have h2 : 0 < 1 + (α₀ - r₀) + (α₀ - r₀) ^ 2 := by
      have := sq_nonneg (α₀ - r₀); linarith
    have h3 : 0 < 1 + α₀ + α₀ ^ 2 := by
      have := sq_nonneg α₀; linarith
    exact mul_pos (mul_pos h1 h2) h3
  have h57 := pandrosion_h_C_p3_quantitative_bound
    α₀ hα₀_ge r₀ hr₀_pos hr₀_le x hx hα_fp z hz_bound
  rw [div_mul_eq_mul_div, le_div_iff h_B_pos]
  linarith [h57]

/-! ============================================================
  §58.2  Banach existentiel
============================================================ -/

/-- `K ≥ 0` est automatique dès que les paramètres sont valides. -/
theorem pandrosion_p3_rate_nonneg
    (α₀ : ℝ) (hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_pos : 0 < r₀) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (x : ℝ) (hx : x > 1) :
    0 ≤ pandrosion_p3_rate α₀ r₀ x := by
  unfold pandrosion_p3_rate
  have h_num_nn : 0 ≤ (x - 1) * (1 + 2 * α₀ + r₀) :=
    mul_nonneg (by linarith) (by linarith)
  have h_B_pos : 0 < x * (1 + (α₀ - r₀) + (α₀ - r₀) ^ 2) * (1 + α₀ + α₀ ^ 2) := by
    have h1 : 0 < x := by linarith
    have h2 : 0 < 1 + (α₀ - r₀) + (α₀ - r₀) ^ 2 := by
      have := sq_nonneg (α₀ - r₀); linarith
    have h3 : 0 < 1 + α₀ + α₀ ^ 2 := by
      have := sq_nonneg α₀; linarith
    exact mul_pos (mul_pos h1 h2) h3
  exact div_nonneg h_num_nn (le_of_lt h_B_pos)

/-- **★★★ Banach-style contraction for p=3 on a disk (existential).** -/
theorem pandrosion_h_C_p3_banach_contraction
    (α₀ : ℝ) (hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_pos : 0 < r₀) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (x : ℝ) (hx : x > 1) (hα_fp : (α₀ : ℝ) ^ 3 = 1 / x)
    (hK_lt_one : pandrosion_p3_rate α₀ r₀ x < 1)
    (z : ℂ) (hz_bound : ‖z - (α₀ : ℂ)‖ ≤ r₀) :
    ∃ K : ℝ, 0 ≤ K ∧ K < 1 ∧
      ‖pandrosion_h_C (x : ℂ) 3 z - (α₀ : ℂ)‖ ≤ K * ‖z - (α₀ : ℂ)‖ :=
  ⟨pandrosion_p3_rate α₀ r₀ x,
   pandrosion_p3_rate_nonneg α₀ hα₀_ge r₀ hr₀_pos hr₀_le x hx,
   hK_lt_one,
   pandrosion_h_C_p3_ratio_bound α₀ hα₀_ge r₀ hr₀_pos hr₀_le x hx hα_fp z hz_bound⟩

/-! ============================================================
  §58.3  Disk invariance
============================================================ -/

/-- **Disk invariance** : sous `K ≤ 1`, `h` envoie `B(α₀, r₀)` dans
    lui-même. -/
theorem pandrosion_h_C_p3_disk_invariance
    (α₀ : ℝ) (hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_pos : 0 < r₀) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (x : ℝ) (hx : x > 1) (hα_fp : (α₀ : ℝ) ^ 3 = 1 / x)
    (hK_le_one : pandrosion_p3_rate α₀ r₀ x ≤ 1)
    (z : ℂ) (hz_bound : ‖z - (α₀ : ℂ)‖ ≤ r₀) :
    ‖pandrosion_h_C (x : ℂ) 3 z - (α₀ : ℂ)‖ ≤ r₀ := by
  have h_rate_nn : 0 ≤ pandrosion_p3_rate α₀ r₀ x :=
    pandrosion_p3_rate_nonneg α₀ hα₀_ge r₀ hr₀_pos hr₀_le x hx
  have h_ratio := pandrosion_h_C_p3_ratio_bound
    α₀ hα₀_ge r₀ hr₀_pos hr₀_le x hx hα_fp z hz_bound
  have h_chain : pandrosion_p3_rate α₀ r₀ x * ‖z - (α₀ : ℂ)‖
               ≤ 1 * r₀ := by
    have h1 : pandrosion_p3_rate α₀ r₀ x * ‖z - (α₀ : ℂ)‖
            ≤ pandrosion_p3_rate α₀ r₀ x * r₀ :=
      mul_le_mul_of_nonneg_left hz_bound h_rate_nn
    have h2 : pandrosion_p3_rate α₀ r₀ x * r₀ ≤ 1 * r₀ :=
      mul_le_mul_of_nonneg_right hK_le_one (le_of_lt hr₀_pos)
    linarith
  linarith [h_ratio, h_chain]

/-! ============================================================
  §58.4  Orbit bound `‖hⁿ(z) − α₀‖ ≤ Kⁿ · ‖z − α₀‖`
============================================================ -/

/-- **★★★★ Borne d'orbite de Banach-Cacciopoli pour p=3.**

    Sous une hypothèse `K ≤ 1` (contraction non-expansive sur le
    disque), l'itération de `h` satisfait la borne géométrique

        `‖hⁿ(z) − α₀‖ ≤ Kⁿ · ‖z − α₀‖`.

    Si de plus `K < 1`, cela donne convergence linéaire vers `α₀`. -/
theorem pandrosion_h_C_p3_orbit_bound
    (α₀ : ℝ) (hα₀_ge : α₀ ≥ 3 / 4)
    (r₀ : ℝ) (hr₀_pos : 0 < r₀) (hr₀_le : r₀ ≤ α₀ - 1 / 2)
    (x : ℝ) (hx : x > 1) (hα_fp : (α₀ : ℝ) ^ 3 = 1 / x)
    (hK_le_one : pandrosion_p3_rate α₀ r₀ x ≤ 1)
    (z : ℂ) (hz_bound : ‖z - (α₀ : ℂ)‖ ≤ r₀) :
    ∀ n : ℕ,
      ‖(pandrosion_h_C (x : ℂ) 3)^[n] z - (α₀ : ℂ)‖
        ≤ (pandrosion_p3_rate α₀ r₀ x) ^ n * ‖z - (α₀ : ℂ)‖ := by
  intro n
  induction n with
  | zero =>
    simp
  | succ k ih =>
    -- Invariance: the k-th iterate stays in B(α₀, r₀).
    have h_inv_step : ∀ m,
        ‖(pandrosion_h_C (x : ℂ) 3)^[m] z - (α₀ : ℂ)‖ ≤ r₀ := by
      intro m
      induction m with
      | zero => simpa using hz_bound
      | succ j ih_j =>
        rw [Function.iterate_succ', Function.comp_apply]
        exact pandrosion_h_C_p3_disk_invariance α₀ hα₀_ge r₀ hr₀_pos hr₀_le
          x hx hα_fp hK_le_one _ ih_j
    -- Apply ratio bound at the k-th iterate.
    have h_k_in := h_inv_step k
    have h_step := pandrosion_h_C_p3_ratio_bound α₀ hα₀_ge r₀ hr₀_pos hr₀_le
      x hx hα_fp _ h_k_in
    have h_rate_nn : 0 ≤ pandrosion_p3_rate α₀ r₀ x :=
      pandrosion_p3_rate_nonneg α₀ hα₀_ge r₀ hr₀_pos hr₀_le x hx
    -- Chain: ratio · ih.
    rw [Function.iterate_succ', Function.comp_apply]
    have h_chain : pandrosion_p3_rate α₀ r₀ x *
          ‖(pandrosion_h_C (x : ℂ) 3)^[k] z - (α₀ : ℂ)‖
        ≤ pandrosion_p3_rate α₀ r₀ x *
          (pandrosion_p3_rate α₀ r₀ x ^ k * ‖z - (α₀ : ℂ)‖) :=
      mul_le_mul_of_nonneg_left ih h_rate_nn
    have h_goal_form : pandrosion_p3_rate α₀ r₀ x *
          (pandrosion_p3_rate α₀ r₀ x ^ k * ‖z - (α₀ : ℂ)‖)
        = pandrosion_p3_rate α₀ r₀ x ^ (k + 1) * ‖z - (α₀ : ℂ)‖ := by
      rw [pow_succ]; ring
    show ‖pandrosion_h_C (x : ℂ) 3 ((pandrosion_h_C (x : ℂ) 3)^[k] z) -
            (α₀ : ℂ)‖
          ≤ pandrosion_p3_rate α₀ r₀ x ^ (k + 1) * ‖z - (α₀ : ℂ)‖
    linarith [h_step, h_chain, h_goal_form.ge, h_goal_form.le]

end Pandrosion
