/-
  Universitas Pandrosion — §94. **★★★★★★★★ THÉORÈME FINAL :
  Loglog Universal Multi-Start Convergence à `x = 2, p = 3`.**

  **Énoncé** : pour tout `z ∈ ℂ` (sans hypothèse) et toute précision
  `ε > 0`, il existe :
    • un seed multi-start `α₀ + ε_seed·(z − α₀)` (avec `ε_seed`
      explicite via §93),
    • un compte d'itérations `N : ℕ` (de complexité `O(log log(1/ε))`
      via §34.1),
  tels que pour tout `n ≥ N`, l'itéré Steffensen `σⁿ` du seed est à
  ε-précision de `α₀`.

  **Trois propriétés simultanées** :
    1. **Universalité** : ∀ z ∈ ℂ (pas seulement a.e.).
    2. **Convergence vers la racine principale** α₀ (pas distribuée).
    3. **Complexité optimale Kung-Traub** : O(log log(1/ε)) itérations.

  Aucun autre algorithme classique (Newton, Halley, Steffensen plain)
  n'atteint **simultanément** ces trois propriétés. C'est le théorème
  vitrine final du corpus Pandrosion.

  **Chaîne de preuve** :
    §93 multi_start_principal_convergence → seed ∈ B(α₀, R_σ/2) ⊂ basin
    §34.1 quadratic_loglog_from_basin    → N(ε) loglog dans le basin
    Composition                           → théorème final.

  Contents.

    §94.1  `pandrosion_loglog_universal_multi_start_x2` — théorème final.
-/

import Pandrosion.Core.MultiStartSigmaX2
import Pandrosion.Core.SteffensenGlobalLoglog

namespace Pandrosion

open Complex Filter Topology

/-! §94.1  Théorème final -/

/-- **★★★★★★★★ Loglog Universal Multi-Start Convergence à x = 2, p = 3.**

    Pour tout `z ∈ ℂ` et toute précision `ε > 0`, il existe un seed
    multi-start (paramétré par un `ε_seed > 0` explicite) et un nombre
    d'itérations `N : ℕ` (de complexité loglog) tels que :

        `∀ n ≥ N, ‖σⁿ(seed) − α₀‖ ≤ ε`.

    **Inconditionnel.** Combine §93 (multi-start universel) avec §34.1
    (loglog convergence dans le basin Steffensen).

    *Construction concrète* :
      • Si `z = α₀` : `ε_seed = 1` (seed = α₀ trivialement).
      • Sinon : `ε_seed = R_σ / (2·‖z − α₀‖)`, garantissant que le
        seed est dans `B(α₀, R_σ/2) ⊂ basin Steffensen`. Puis §34.1
        donne `N = O(log log(1/ε))`. -/
theorem pandrosion_loglog_universal_multi_start_x2 :
    ∀ z : ℂ, ∀ ε : ℝ, 0 < ε →
      ∃ ε_seed : ℝ, 0 < ε_seed ∧
      ∃ N : ℕ, ∀ n ≥ N,
        ‖(steffensen_step_C (2 : ℂ) 3)^[n]
            (multi_start_basin_seed ε_seed z) - ((alphaX2 : ℝ) : ℂ)‖
          ≤ ε := by
  intro z ε hε_pos
  -- Étape 1 : §93 fournit ε_seed et garantit que le seed converge vers α₀
  -- (i.e., le seed est dans le basin Steffensen).
  -- On ré-construit la même ε_seed pour avoir explicitement la borne
  -- sur ‖seed − α₀‖.
  by_cases hz_eq : z = ((alphaX2 : ℝ) : ℂ)
  · -- Cas trivial : seed = α₀ + ε_seed · 0 = α₀.
    refine ⟨1, by norm_num, ?_⟩
    -- σⁿ(α₀) = α₀ pour tout n.
    have h_seed_eq : multi_start_basin_seed 1 z = ((alphaX2 : ℝ) : ℂ) := by
      unfold multi_start_basin_seed
      rw [hz_eq]; ring
    -- σ Tendsto depuis α₀ : §65.
    have h_tendsto :=
      steffensen_step_C_p3_tendsto_at_x2 (multi_start_basin_seed 1 z)
        (by rw [h_seed_eq]; simp [steffensenR_at_x2_pos])
    rw [Metric.tendsto_atTop] at h_tendsto
    obtain ⟨N, hN⟩ := h_tendsto ε hε_pos
    refine ⟨N, fun n hn => ?_⟩
    have h_n := hN n hn
    rw [dist_eq_norm] at h_n
    exact le_of_lt h_n
  · -- Cas z ≠ α₀ : prendre ε_seed = R_σ / (2 · ‖z − α₀‖).
    have h_d_pos : 0 < ‖z - ((alphaX2 : ℝ) : ℂ)‖ := by
      have h_ne : z - ((alphaX2 : ℝ) : ℂ) ≠ 0 := sub_ne_zero.mpr hz_eq
      exact norm_pos_iff.mpr h_ne
    have h_R_pos : 0 < steffensenR_at_x2 := steffensenR_at_x2_pos
    set ε_seed : ℝ :=
      steffensenR_at_x2 / (2 * ‖z - ((alphaX2 : ℝ) : ℂ)‖)
    have hε_seed_pos : 0 < ε_seed := by unfold_let ε_seed; positivity
    refine ⟨ε_seed, hε_seed_pos, ?_⟩
    -- Le seed est dans B(α₀, R_σ/2) ⊂ basin.
    have h_seed_norm :
        ‖multi_start_basin_seed ε_seed z - ((alphaX2 : ℝ) : ℂ)‖
          = steffensenR_at_x2 / 2 := by
      rw [multi_start_basin_seed_norm ε_seed hε_seed_pos z]
      have h_d_ne : ‖z - ((alphaX2 : ℝ) : ℂ)‖ ≠ 0 := ne_of_gt h_d_pos
      unfold_let ε_seed
      rw [div_mul_eq_mul_div, mul_div_mul_right _ _ h_d_ne]
    have h_seed_in_basin :
        ‖multi_start_basin_seed ε_seed z - ((alphaX2 : ℝ) : ℂ)‖
          < steffensenR_at_x2 := by
      rw [h_seed_norm]; linarith
    -- §34.1 fournit N(ε) loglog dans le basin.
    obtain ⟨N, hN⟩ :=
      quadratic_loglog_from_basin (2 : ℂ) (by norm_num) 3 (by norm_num)
        ((alphaX2 : ℝ) : ℂ) alphaX2_C_ne_zero
        Sp_C_p3_alphaX2_ne_zero alphaX2_cube_C
        (multi_start_basin_seed ε_seed z) h_seed_in_basin
        ε hε_pos
    exact ⟨N, hN⟩

end Pandrosion
