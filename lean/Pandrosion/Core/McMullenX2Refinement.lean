/-
  Universitas Pandrosion — §66. **Refinement of `McMullenAEEntry 3 2 α`
  à `x = 2` : réduction à un conjecture Fatou-entering explicite.**

  **But.** Repackager la conjecture ouverte `ComplexJuliaP3MeasureZero 2 α`
  (≡ `McMullenAEEntry 3 2 α`) en une forme plus concrète et plus
  petite : **Fatou-entering hypothesis**, qui ne demande que l'entrée
  dans les disques Steffensen explicites §33 (pas un contrôle fin
  du Julia).

  Chaîne démontrée inconditionnellement ici :

      `FatouEnteringP3X2`
          ⟹ `McMullenAEEntry 3 (2 : ℂ) α`     (§66.2 — démontré)
          ⟹ convergence a.e. vers racines cyclotomiques  (§41.5 — chaîné)

  La **conjecture résiduelle** `FatouEnteringP3X2` est :

      *"Pour Lebesgue-presque tout `z₀ ∈ ℂ`, il existe `s : Fin 3` et
      `k : ℕ` tels que `σᵏ z₀` entre dans le disque Steffensen autour
      de la racine cyclotomique `γ_s = cycAnchor α 3 s`."*

  **Validation empirique forte** (Python §61) : à `x = 2`, sur
  360 000 points et rings ε = 10⁻⁹ autour des singularités,
  **100 % des points non-singuliers** convergent (donc entrent). La
  conjecture est validée numériquement sans ambiguïté.

  **Status** : la partie mesure-théorique (§62, §65) est acquise
  inconditionnellement. La conjecture résiduelle est purement
  dynamique, et plus concrète que le `McMullenAEEntry` original car
  elle n'exige que l'entrée dans des disques *explicites* (pas tout
  voisinage arbitrairement petit).

  Contents.

    §66.1  `Sp_C_p3_cycAnchor_ne_zero_at_x2` — non-dégénérescence
           aux trois anchors.

    §66.2  `steffensenR_at_x2_s` — rayon Steffensen par anchor.

    §66.3  `FatouEnteringP3X2` — conjecture résiduelle explicite.

    §66.4  `mcmullen_ae_entry_3_2_alphaX2_from_fatou` — chaînage
           inconditionnel Fatou ⟹ McMullen.

    §66.5  `steffensen_p3_solves_at_x2_from_fatou` — chaînage
           final vers la convergence a.e.
-/

import Pandrosion.Core.SteffensenX2Convergence
import Pandrosion.Core.ComplexMcMullenP3Conditional

namespace Pandrosion

open Complex MeasureTheory Filter Topology

/-! ============================================================
  §66.1  Non-dégénérescence des anchors à `x = 2`
============================================================ -/

/-- **`Sp_C 3 γ_s ≠ 0`** à `α = 2^{-1/3}` pour chaque `s : Fin 3`.
    Les trois racines cyclotomiques (`α, ω·α, ω²·α`) sont toutes
    distinctes des `Sp`-zéros `{ω, ω²}` car `α ≠ 1, ω, ω²`. -/
theorem Sp_C_p3_cycAnchor_ne_zero_at_x2 (s : Fin 3) :
    Sp_C 3 (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s) ≠ 0 := by
  -- `(cycAnchor α 3 s)³ = α³ = 1/2`, donc cycAnchor ≠ {ω, ω²}
  -- (car ω³ = ω²³ = 1 ≠ 1/2).
  intro h_zero
  rw [Sp_C_p3_zero_iff_cube_root] at h_zero
  obtain ⟨h_cube, _⟩ := h_zero
  have h_anchor_cube : (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s) ^ 3
                     = ((alphaX2 : ℝ) : ℂ) ^ 3 :=
    cycAnchor_pow _ (by norm_num : 1 ≤ 3) s
  rw [alphaX2_cube_C] at h_anchor_cube
  rw [h_anchor_cube] at h_cube
  -- h_cube : (1 : ℂ)/2 = 1, contradiction via norm_num.
  norm_num at h_cube

/-! ============================================================
  §66.2  Per-anchor Steffensen radius at `x = 2`
============================================================ -/

/-- Rayon Steffensen autour de l'anchor `γ_s` à `x = 2`. -/
noncomputable def steffensenR_at_x2_s (s : Fin 3) : ℝ :=
  steffensenR_of_fp (2 : ℂ) (by norm_num) 3 (by norm_num)
    (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s)
    (cycAnchor_ne_zero alphaX2_C_ne_zero 3 s)
    (Sp_C_p3_cycAnchor_ne_zero_at_x2 s)
    (by rw [cycAnchor_pow _ (by norm_num : 1 ≤ 3) s]; exact alphaX2_cube_C)

/-- **`steffensenR_at_x2_s s > 0`** pour tout `s`. -/
theorem steffensenR_at_x2_s_pos (s : Fin 3) :
    0 < steffensenR_at_x2_s s :=
  steffensenR_of_fp_pos _ _ _ _ _ _ _ _

/-! ============================================================
  §66.3  Fatou-entering hypothesis (résiduelle)
============================================================ -/

/-- **Conjecture Fatou-entering à `x = 2`.**

    Pour Lebesgue-presque tout `z₀ ∈ ℂ`, l'itération σ entre en
    temps fini dans l'un des trois disques Steffensen explicites
    autour des racines cyclotomiques cubiques de `1/2`. -/
def FatouEnteringP3X2 : Prop :=
  ∀ᵐ z₀ : ℂ ∂volume,
    ∃ s : Fin 3, ∃ k : ℕ,
      ‖(steffensen_step_C (2 : ℂ) 3)^[k] z₀
          - cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s‖
        < steffensenR_at_x2_s s

/-! ============================================================
  §66.4  Chaînage inconditionnel Fatou ⟹ McMullen
============================================================ -/

/-- **★★★★★ `FatouEnteringP3X2 ⟹ McMullenAEEntry 3 2 α`.**

    Chaînage **inconditionnel** : une fois entré dans un disque
    Steffensen explicite, §34.1 `quadratic_loglog_from_basin` donne
    l'entrée dans tout disque plus petit (convergence quadratique
    dans le basin). -/
theorem mcmullen_ae_entry_3_2_alphaX2_from_fatou
    (hFatou : FatouEnteringP3X2) :
    McMullenAEEntry 3 (2 : ℂ) ((alphaX2 : ℝ) : ℂ) := by
  intro r hr
  filter_upwards [hFatou] with z₀ hz₀
  obtain ⟨s, k₀, hk₀⟩ := hz₀
  -- Once inside the Steffensen disk, quadratic_loglog_from_basin
  -- gives entry into any smaller r s.
  have hα_pow_s : (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s) ^ 3 = 1 / (2 : ℂ) := by
    rw [cycAnchor_pow _ (by norm_num : 1 ≤ 3) s]; exact alphaX2_cube_C
  obtain ⟨N, hN⟩ :=
    quadratic_loglog_from_basin (2 : ℂ) (by norm_num) 3 (by norm_num)
      (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s)
      (cycAnchor_ne_zero alphaX2_C_ne_zero 3 s)
      (Sp_C_p3_cycAnchor_ne_zero_at_x2 s)
      hα_pow_s
      ((steffensen_step_C (2 : ℂ) 3)^[k₀] z₀) hk₀
      (r s / 2) (by linarith [hr s])
  refine ⟨s, N + k₀, ?_⟩
  rw [Function.iterate_add_apply]
  have h_le : ‖(steffensen_step_C (2 : ℂ) 3)^[N]
                ((steffensen_step_C (2 : ℂ) 3)^[k₀] z₀)
              - cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s‖
            ≤ r s / 2 := hN N (le_refl N)
  linarith [hr s, h_le]

/-! ============================================================
  §66.5  Final chain: a.e. solve at `x = 2`
============================================================ -/

/-- **★★★★★★ `FatouEnteringP3X2 ⟹ Pandrosion résout `z³ = 2` a.e.`** -/
theorem steffensen_p3_solves_at_x2_from_fatou
    (hFatou : FatouEnteringP3X2) :
    ∀ᵐ z₀ : ℂ ∂volume,
      ∃ s : Fin 3,
        Tendsto (fun k => (steffensen_step_C (2 : ℂ) 3)^[k] z₀) atTop
          (𝓝 (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s)) :=
  steffensen_solves_ae_mod_mcmullen 3 (by norm_num) (2 : ℂ) (by norm_num)
    ((alphaX2 : ℝ) : ℂ) alphaX2_C_ne_zero alphaX2_cube_C
    (fun s => Sp_C_p3_cycAnchor_ne_zero_at_x2 s)
    (mcmullen_ae_entry_3_2_alphaX2_from_fatou hFatou)

end Pandrosion
