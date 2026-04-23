/-
  Universitas Pandrosion — §71. **Bassin principal Master : briques
  inconditionnelles.**

  ⚠️ **REFONDU après réfutation Niveau 5 strict** ⚠️

  La version originale prouvait `PrincipalDominanceP3X2 ⟺ ∀ᵐ z ∈
  PrincipalBasinP3X2`. Cette conjecture est **strictement fausse**
  (bassins de ω·α₀ ont mesure positive — voir §69 header).

  Ce module garde donc seulement les briques inconditionnelles sur le
  bassin principal :
    (A) `α₀ ∈ PrincipalBasinP3X2` (point fixe).
    (B) Points concrets `z = 1, real ≥ 1/2` dans le h-bassin.
    (C) Le bassin principal est non-vide.

  Le chaînage vers Niveau 5 est supprimé. Pour atteindre la
  conjecture **correcte** (Niveau 2 = `McMullenAEEntry 3 2 α₀`),
  voir §79 refondu.

  Contents.

    §71.1  `alpha_in_principal_basin` — α₀ ∈ bassin principal.
    §71.2  `one_in_h_principal_basin`, `real_ge_half_in_h_principal_basin`.
    §71.3  `principal_basin_nonempty`.
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
  apply steffensen_step_C_p3_tendsto_at_x2
  simp [steffensenR_at_x2_pos]

/-! §71.2  Points concrets dans h-bassin (via §67 half-plane) -/

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

end Pandrosion
