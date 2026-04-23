/-
  Universitas Pandrosion — §79. **Chaînage vers `McMullenAEEntry 3 2 α₀`
  (Niveau 2).**

  ⚠️ **REFONDU après réfutation Niveau 5 strict** ⚠️

  La version originale chaînait vers `PrincipalDominanceP3X2`
  (`∀ᵐ z, σⁿ(z) → α₀`). Cette conjecture est **strictement fausse**
  (voir §69 header).

  Cette version chaîne vers la conjecture **correcte** :
  `McMullenAEEntry 3 2 α₀` (Niveau 2 — `∀ᵐ z, ∃ s, σⁿ(z) → γ_s`).

  Forme du théorème :

      `MultiBasinAEConvergence`     (⟺ McMullenAEEntry 3 2 α₀)
    + `HalfPlaneSigmaTendstoP3X2`   (Niveau 1, §70 — auxiliary)
    ⟹ `McMullenAEEntry 3 2 α₀`     (chaîne directe)

  Note : la conjecture `HalfPlaneExhaustionAEP3X2` originale est
  également suspecte — les points dans le bassin de ω·α₀ ne sortent
  PAS vers le demi-plan (ils restent près de ω·α₀). Donc elle n'est
  pas universellement vraie ; elle vaut seulement pour z dans le
  bassin principal.

  Contents.

    §79.1  `MultiBasinAEConvergence` — la conjecture cible (= Niveau 2).
    §79.2  `mcmullen_ae_entry_x2_iff_multibasin` — équivalence avec
           McMullenAEEntry.
-/

import Pandrosion.Core.SigmaTendstoAtOneX2
import Pandrosion.Core.Niveau5Master
import Pandrosion.Core.McMullenX2Refinement

namespace Pandrosion

open Complex MeasureTheory Filter Topology

/-! §79.1  Conjecture cible (Niveau 2, version σ) -/

/-- **Conjecture multi-bassin AE** : presque tout `z ∈ ℂ` converge
    vers **une** des trois racines cyclotomiques. C'est la forme
    correcte (vraie empiriquement) de la conjecture Niveau 2/5.

    Différence avec la fausse `PrincipalDominanceP3X2` :
      • PrincipalDominance : `∀ᵐ z, σⁿ(z) → α₀` (FAUX).
      • MultiBasinAE : `∀ᵐ z, ∃ s, σⁿ(z) → γ_s` (vrai). -/
def MultiBasinAEConvergence : Prop :=
  ∀ᵐ z : ℂ ∂volume,
    ∃ s : Fin 3,
      Tendsto (fun k : ℕ => (steffensen_step_C (2 : ℂ) 3)^[k] z) atTop
        (𝓝 (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s))

/-! §79.2  Équivalence avec McMullenAEEntry -/

/-- **`MultiBasinAEConvergence ⟹ McMullenAEEntry 3 2 α₀`.**

    Direct par déballage du Tendsto. -/
theorem mcmullen_ae_entry_x2_from_multibasin
    (hMB : MultiBasinAEConvergence) :
    McMullenAEEntry 3 (2 : ℂ) ((alphaX2 : ℝ) : ℂ) := by
  intro r hr
  filter_upwards [hMB] with z₀ hz₀
  obtain ⟨s, h_tendsto⟩ := hz₀
  refine ⟨s, ?_⟩
  have h_pos : 0 < r s := hr s
  have h_eventually :
      ∀ᶠ n in atTop,
        ‖(steffensen_step_C (2 : ℂ) 3)^[n] z₀
            - cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s‖ < r s := by
    have := Metric.tendsto_atTop.mp h_tendsto (r s) h_pos
    obtain ⟨N, hN⟩ := this
    filter_upwards [Filter.eventually_atTop.mpr ⟨N, fun n hn => hn⟩] with n hn
    have := hN n hn
    rwa [dist_eq_norm] at this
  exact h_eventually.exists

/-! §79.3  Niveau 1 (HalfPlaneSigmaTendsto) ⟹ MultiBasin sur demi-plan -/

/-- **Niveau 1 ⟹ MultiBasin restreint au demi-plan**.

    Sous `HalfPlaneSigmaTendstoP3X2`, le demi-plan est entièrement dans
    `PrincipalBasinP3X2` (= `CyclotomicBasinP3X2 0`), donc satisfait
    a fortiori la convergence vers une racine cyclotomique. -/
theorem multibasin_on_half_plane_from_HP_tendsto
    (hHP : HalfPlaneSigmaTendstoP3X2) (z : ℂ) (hz : z.re ≥ 1/2) :
    ∃ s : Fin 3,
      Tendsto (fun k : ℕ => (steffensen_step_C (2 : ℂ) 3)^[k] z) atTop
        (𝓝 (cycAnchor ((alphaX2 : ℝ) : ℂ) 3 s)) :=
  ⟨0, by rw [cycAnchor_p3_zero]; exact hHP z hz⟩

end Pandrosion
