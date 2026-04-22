/-
  Universitas Pandrosion — §86. **Compte d'itérations pour `ε`-précision
  à `x = 2, p = 3`.**

  Formalise la complexité numérique : pour tout `ε > 0` et tout
  `z ∈ B(α₀, 1/4)`, il existe un nombre d'itérations `N` tel que
  `‖h^n(z) − α₀‖ ≤ ε` pour tout `n ≥ N`.

  **Forme existentielle** (cette session) : déballage direct du
  Tendsto §64.3 via `Metric.tendsto_atTop`.

  **Forme explicite future** (§86b) : `N(ε) = ⌈log(1/(4ε)) / log(259/104)⌉`,
  exigeant manipulation `Real.log` + chaînage avec §64.5
  `pandrosion_h_C_p3_rate_bound_at_x2`.

  Contents.

    §86.1  `pandrosion_h_C_p3_iter_count_at_x2_existential` — existence
           du `N`.
    §86.2  `pandrosion_h_C_p3_iter_count_at_x2_geometric_bound` — borne
           explicite `(104/259)^N · (1/4) ≤ ε` ⟹ `N` suffit.
-/

import Pandrosion.Core.BanachX2Concrete

namespace Pandrosion

open Complex Filter Topology

/-! §86.1  Forme existentielle -/

/-- **Existence d'un compte d'itérations pour `ε`-précision.**

    Pour tout `ε > 0` et `z ∈ B(α₀, 1/4)`, il existe `N : ℕ` tel que
    pour tout `n ≥ N`, `‖h^n(z) − α₀‖ ≤ ε`. -/
theorem pandrosion_h_C_p3_iter_count_at_x2_existential
    (ε : ℝ) (hε : 0 < ε) (z : ℂ)
    (hz : ‖z - ((alphaX2 : ℝ) : ℂ)‖ ≤ (1/4 : ℝ)) :
    ∃ N : ℕ, ∀ n ≥ N,
      ‖(pandrosion_h_C (2 : ℂ) 3)^[n] z - ((alphaX2 : ℝ) : ℂ)‖ ≤ ε := by
  have h_tendsto := pandrosion_h_C_p3_banach_at_x2 z hz
  rw [Metric.tendsto_atTop] at h_tendsto
  obtain ⟨N, hN⟩ := h_tendsto ε hε
  refine ⟨N, fun n hn => ?_⟩
  have h := hN n hn
  rw [dist_eq_norm] at h
  exact le_of_lt h

/-! §86.2  Borne explicite `(104/259)^N · (1/4) ≤ ε` ⟹ N suffit -/

/-- **Borne explicite `(104/259)^N · (1/4) ≤ ε` ⟹ N itérations
    suffisent.**

    Permet d'extraire un `N` concret en résolvant la borne géométrique
    (sans passer par `Real.log` formel). -/
theorem pandrosion_h_C_p3_iter_count_at_x2_geometric_bound
    (ε : ℝ) (z : ℂ)
    (hz : ‖z - ((alphaX2 : ℝ) : ℂ)‖ ≤ (1/4 : ℝ))
    (N : ℕ) (hN : (104/259 : ℝ)^N * (1/4) ≤ ε)
    (n : ℕ) (hn : n ≥ N) :
    ‖(pandrosion_h_C (2 : ℂ) 3)^[n] z - ((alphaX2 : ℝ) : ℂ)‖ ≤ ε := by
  have h_rate := pandrosion_h_C_p3_rate_bound_at_x2 z hz n
  -- h_rate : ‖h^n z - α₀‖ ≤ (104/259)^n · (1/4)
  -- want : ≤ ε. Suffit (104/259)^n · (1/4) ≤ ε.
  have h_104_nn : (0 : ℝ) ≤ 104/259 := by norm_num
  have h_pow_mono : (104/259 : ℝ)^n ≤ (104/259 : ℝ)^N :=
    pow_le_pow_of_le_one h_104_nn (by norm_num : (104/259 : ℝ) ≤ 1) hn
  have h_zm_nn : (0 : ℝ) ≤ (1 : ℝ)/4 := by norm_num
  have h_chain : (104/259 : ℝ)^n * (1/4) ≤ (104/259 : ℝ)^N * (1/4) :=
    mul_le_mul_of_nonneg_right h_pow_mono h_zm_nn
  linarith

end Pandrosion
