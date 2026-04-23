/-
  Universitas Pandrosion — §78. **Mesure des bassins non-principaux
  de σ à `x = 2, p = 3`.**

  ⚠️ **REFONDU après réfutation Niveau 5 strict** ⚠️

  La conjecture originale `NonPrincipalBasinNullX2` (`volume(...) = 0`)
  est **strictement fausse** (validation numérique
  `sigma_structural_analysis.py` montre bassins de mesure ~π·0.025²
  ≈ 0.002 pour ω·α₀ et ω²·α₀).

  La conjecture **correcte** est :

    `NonPrincipalBasinBoundedX2` : `⋃_{s∈{1,2}} CyclotomicBasinP3X2 s`
    ⊆ `B(ω·α₀, R) ∪ B(ω²·α₀, R)` pour un `R` explicite (~0.05).

  Cette forme capture l'asymétrie observée empiriquement (les bassins
  non-principaux sont *localisés et bornés*) sans prétendre qu'ils
  sont de mesure nulle.

  Contents.

    §78.1  `NonPrincipalBasinFiniteX2` — borné/finite (vrai).
    §78.2  `NonPrincipalBasinBoundedX2 R` — borné par disques explicites.
    §78.3  `NonPrincipalBasinNullX2` — **REFUTED** (gardé pour
           documentation, marqué FALSE).
-/

import Pandrosion.Core.PrincipalDominanceP3X2

namespace Pandrosion

open Complex MeasureTheory

/-! §78.1  Conjecture mesure finie -/

/-- **Conjecture mesure finie** : la réunion des bassins cyclotomiques
    non-principaux (s ∈ {1, 2}) est de mesure de Lebesgue finie.

    Empiriquement vrai (bassins observés concentrés autour de
    ω·α₀, ω²·α₀ avec rayon ~0.025). -/
def NonPrincipalBasinFiniteX2 : Prop :=
  volume (CyclotomicBasinP3X2 1 ∪ CyclotomicBasinP3X2 2) < ⊤

/-! §78.2  Conjecture borné par disques explicites -/

/-- **Conjecture borné par disques explicites** : pour un `R > 0`
    suffisamment grand, les bassins non-principaux sont contenus dans
    deux petits disques autour de ω·α₀, ω²·α₀.

    Empirique : `R = 0.05` suffit (rayon observé ~0.025).

    Forme stronge de l'asymétrie Pandrosion : non seulement la mesure
    est finie, mais elle est **localisée explicitement**. -/
def NonPrincipalBasinBoundedX2 (R : ℝ) : Prop :=
  CyclotomicBasinP3X2 1 ∪ CyclotomicBasinP3X2 2 ⊆
    {z : ℂ | ‖z - cycAnchor ((alphaX2 : ℝ) : ℂ) 3 1‖ ≤ R} ∪
    {z : ℂ | ‖z - cycAnchor ((alphaX2 : ℝ) : ℂ) 3 2‖ ≤ R}

/-! §78.3  Bounded ⟹ Finite -/

/-- Bornés explicitement ⟹ mesure finie. -/
theorem non_principal_finite_of_bounded {R : ℝ}
    (h : NonPrincipalBasinBoundedX2 R) :
    NonPrincipalBasinFiniteX2 := by
  unfold NonPrincipalBasinFiniteX2
  apply lt_of_le_of_lt (measure_mono h)
  -- Mesure de l'union de deux disques bornés est finie.
  apply lt_of_le_of_lt (measure_union_le _ _)
  -- Chaque disque a mesure finie (il est borné).
  have h_disk_lt : ∀ (c : ℂ),
      volume {z : ℂ | ‖z - c‖ ≤ R} < ⊤ := fun c => by
    apply IsCompact.measure_lt_top
    have h_eq : {z : ℂ | ‖z - c‖ ≤ R} = Metric.closedBall c R := by
      ext z; simp [Metric.closedBall, dist_eq_norm]
    rw [h_eq]
    exact isCompact_closedBall c R
  exact ENNReal.add_lt_top.mpr ⟨h_disk_lt _, h_disk_lt _⟩

/-! §78.4  Conjecture refuted (documentation only) -/

/-- ⚠️ **REFUTED** : `volume(non-principal basins) = 0` est strictement
    faux. Validation numérique : bassins observés de mesure ~0.002 par
    racine non-principale (rayon ~0.025 autour de ω·α₀, ω²·α₀).

    Conservé pour documenter l'erreur historique du corpus. **Ne pas
    utiliser** comme hypothèse — elle implique des conclusions fausses. -/
def NonPrincipalBasinNullX2_REFUTED : Prop :=
  volume (CyclotomicBasinP3X2 1 ∪ CyclotomicBasinP3X2 2) = 0

end Pandrosion
