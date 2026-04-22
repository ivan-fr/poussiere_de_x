/-
  Universitas Pandrosion — §71. **Niveau 5 Master : recapitulation
  inconditionnelle du bassin principal.**

  Module consolidant toutes les briques inconditionnelles vers la
  conjecture `PrincipalDominanceP3X2`. On établit :

    (A) Le bassin principal est **non-vide** (contient `α₀`, `1`, tous
        les reals ≥ 1/2 pour `h`).
    (B) Une sous-partie de `ℂ` de **mesure positive** est dans le
        bassin principal (Steffensen disk §65).
    (C) Le demi-plan droit entier est dans le bassin principal **pour
        l'itération `h`** (§67). L'extension à σ demande Niveau 1.
    (D) Décomposition Niveau 5 : `PrincipalDominance` équivaut à
        "tout z ∈ ℂ non-singulier est dans `PrincipalBasinP3X2`".
    (E) Chaînage final : Niveau 1 + extension au complément du
        demi-plan ⟹ Niveau 5.

  Contents.

    §71.1  `alpha_in_principal_basin` — α₀ ∈ bassin principal.
    §71.2  `one_in_h_principal_basin` — 1 ∈ bassin h-principal.
    §71.3  `principal_basin_nonempty` — le bassin est non-vide.
    §71.4  `niveau5_decomposition_theorem` — décomposition finale.
-/

import Pandrosion.Core.SigmaHalfPlaneP3X2

namespace Pandrosion

open Complex Filter Topology MeasureTheory

/-! §71.1  α₀ ∈ bassin principal -/

/-- **α₀ est un point fixe de σ à x = 2**, donc `σⁿ(α₀) = α₀ → α₀`. -/
theorem alpha_in_principal_basin :
    ((alphaX2 : ℝ) : ℂ) ∈ PrincipalBasinP3X2 := by
  unfold PrincipalBasinP3X2
  simp only [Set.mem_setOf_eq]
  -- α₀ dans B(α₀, R_σ) car ‖α₀ - α₀‖ = 0 < R_σ.
  apply steffensen_step_C_p3_tendsto_at_x2
  simp [steffensenR_at_x2_pos]

/-! §71.2  1 ∈ bassin h-principal (via §67 half-plane) -/

/-- **Le point `z = 1` a `h^n(1) → α₀`** via §67 (demi-plan h). -/
theorem one_in_h_principal_basin :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] (1 : ℂ))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply pandrosion_h_C_p3_tendsto_half_plane_x2
  show (1 : ℂ).re ≥ 1/2
  rw [Complex.one_re]; norm_num

/-- **Tout nombre réel `x ≥ 1/2` donne `h^n(x) → α₀`** (h-iteration). -/
theorem real_ge_half_in_h_principal_basin (x : ℝ) (hx : x ≥ 1/2) :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] ((x : ℂ)))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply pandrosion_h_C_p3_tendsto_half_plane_x2
  show ((x : ℂ)).re ≥ 1/2
  rw [Complex.ofReal_re]; linarith

/-! §71.3  `PrincipalBasinP3X2` est non-vide -/

/-- **Le bassin principal σ est non-vide.** Contient `α₀` lui-même. -/
theorem principal_basin_nonempty : PrincipalBasinP3X2.Nonempty :=
  ⟨((alphaX2 : ℝ) : ℂ), alpha_in_principal_basin⟩

/-! §71.4  Décomposition Niveau 5 -/

/-- **★★★★★★★ Décomposition Niveau 5.**

    `PrincipalDominanceP3X2` équivaut à l'énoncé géométrique
    "pour Lebesgue-presque tout `z ∈ ℂ`, `z ∈ PrincipalBasinP3X2`".

    C'est la reformulation **ensembliste** pure de la conjecture. -/
theorem niveau5_decomposition_theorem :
    PrincipalDominanceP3X2 ↔ ∀ᵐ z : ℂ ∂volume, z ∈ PrincipalBasinP3X2 := by
  unfold PrincipalDominanceP3X2 PrincipalBasinP3X2
  rfl

/-- **Décomposition en hypothèses atomiques.** Si (a) tout le demi-plan
    `Re z ≥ 1/2` est dans le bassin principal (Niveau 1 = §70 conjecture
    C2), et (b) presque tout `z` avec `Re z < 1/2` admet un itéré σ^k
    dans `Re z ≥ 1/2`, alors Niveau 5 suit (pour presque tout z, puisque
    le complément du demi-plan est envoyé dans le demi-plan). -/
theorem niveau5_from_atoms
    (hHP : HalfPlaneSigmaTendstoP3X2)
    (hExhaust :
      ∀ᵐ z : ℂ ∂volume, z.re ≥ 1/2 ∨
        ∃ k : ℕ, ((steffensen_step_C (2 : ℂ) 3)^[k] z).re ≥ 1/2) :
    PrincipalDominanceP3X2 := by
  unfold PrincipalDominanceP3X2
  filter_upwards [hExhaust] with z hz
  rcases hz with h_half | ⟨k, h_enter⟩
  · exact hHP z h_half
  · -- σ^k z est dans le demi-plan, donc σⁿ(σᵏ z) → α₀ par HP,
    -- puis σ^(n+k)(z) → α₀, ce qui équivaut à σⁿ(z) → α₀.
    have h_tendsto_shifted : Tendsto
        (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n]
                         ((steffensen_step_C (2 : ℂ) 3)^[k] z))
        atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := hHP _ h_enter
    have h_comp : (fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n]
                    ((steffensen_step_C (2 : ℂ) 3)^[k] z))
                = fun n : ℕ => (steffensen_step_C (2 : ℂ) 3)^[n + k] z := by
      funext n
      rw [← Function.iterate_add_apply]
    rw [h_comp] at h_tendsto_shifted
    exact (Filter.tendsto_add_atTop_iff_nat k).mp h_tendsto_shifted

end Pandrosion
