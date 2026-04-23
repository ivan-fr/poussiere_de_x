/-
  Universitas Pandrosion — §105. **★★★★★★★★★★ Continuité a.e. de
  la fonction solution multi-start Pandrosion.**

  Complète la trilogie analytique §102 (mesurabilité) + §103 (Julia-
  null) + §105 (continuité a.e.) pour caractériser complètement la
  fonction solution :

    **(I) Mesurable partout** (§102).
    **(II) Julia-null set** (§103).
    **(III) Continuë hors du Julia set** (§105 — le présent module).

  Résultat principal :

    `∀ᵐ z : ℂ ∂volume, ContinuousAt (niveau5_multi_start_target_borel α p hp) z`

  C'est le théorème de type **Luzin** pour Pandrosion — la fonction
  solution est mesurable partout et continue presque-partout.

  Stratégie de preuve :

    1. Sur le Fatou set, il y a une racine cyclotomique *strictement*
       plus proche que toutes les autres.
    2. Par continuité de `z ↦ ‖z − γ_t‖`, cette inégalité stricte se
       propage dans un voisinage ouvert.
    3. Donc `nearest_anchor_borel` est **localement constante** sur
       le Fatou set.
    4. Localement constante ⇒ continue.
    5. `niveau5_multi_start_target_borel = cycAnchor ∘ nearest_anchor_borel`
       est continue sur le Fatou set.
    6. Fatou set a mesure pleine (§103) ⇒ continue a.e.

  Contents.

    §105.1  `fatou_strict_min` — sur Fatou, min strict parmi anchors.
    §105.2  `nearest_anchor_borel_eventually_const` — locally constant.
    §105.3  `nearest_anchor_borel_continuousAt_fatou` — ContinuousAt.
    §105.4  `niveau5_multi_start_target_borel_continuousAt_fatou` —
             ContinuousAt pour la composition.
    §105.5  ★★★ `niveau5_multi_start_target_borel_ae_continuous` —
             a.e. ContinuousAt (théorème principal).
-/

import Pandrosion.Core.PandrosionJuliaNull

namespace Pandrosion

open Complex Filter Topology MeasureTheory

/-! ============================================================
  §105.1  Min strict sur le Fatou set
============================================================ -/

section FatouStrictMin

/-- **Min strict sur Fatou.** Pour `z₀ ∈ PandrosionFatouSet α p`, la
    racine cyclotomique `nearest_anchor_borel α p hp z₀` est
    *strictement* plus proche de `z₀` que toutes les autres. -/
theorem fatou_strict_min
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z₀ : ℂ)
    (hz₀ : z₀ ∈ PandrosionFatouSet α p) :
    ∀ t : Fin p, t ≠ nearest_anchor_borel α p hp z₀ →
      ‖z₀ - cycAnchor α p (nearest_anchor_borel α p hp z₀)‖
        < ‖z₀ - cycAnchor α p t‖ := by
  intro t ht_ne
  have h_le := nearest_anchor_borel_spec α p hp z₀ t
  rcases lt_or_eq_of_le h_le with h_lt | h_eq
  · exact h_lt
  · exfalso
    unfold PandrosionFatouSet PandrosionJuliaSet at hz₀
    rw [Set.mem_compl_iff] at hz₀
    apply hz₀
    refine ⟨nearest_anchor_borel α p hp z₀, t, fun h => ht_ne h.symm, ?_, ?_⟩
    · rw [norm_sub_rev (cycAnchor α p (nearest_anchor_borel α p hp z₀)) z₀,
          norm_sub_rev (cycAnchor α p t) z₀]
      exact h_eq
    · intro u
      rw [norm_sub_rev (cycAnchor α p (nearest_anchor_borel α p hp z₀)) z₀,
          norm_sub_rev (cycAnchor α p u) z₀]
      exact nearest_anchor_borel_spec α p hp z₀ u

end FatouStrictMin

/-! ============================================================
  §105.2  `nearest_anchor_borel` localement constante sur Fatou
============================================================ -/

section LocallyConstant

/-- **★★ Locally constant au point Fatou.** Pour `z₀ ∈ Fatou`, il
    existe un voisinage ouvert `U ∋ z₀` où `nearest_anchor_borel z = s`
    pour tout `z ∈ U`, avec `s = nearest_anchor_borel z₀`.

    Preuve : par continuité de `z ↦ ‖z − γ_t‖ − ‖z − γ_s‖`, chaque
    inégalité stricte (min strict sur Fatou, §105.1) se propage dans
    un voisinage. Intersection finie de voisinages = voisinage où
    l'argmin est `{s}`. -/
theorem nearest_anchor_borel_eventually_const
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z₀ : ℂ)
    (hz₀ : z₀ ∈ PandrosionFatouSet α p) :
    ∀ᶠ z in 𝓝 z₀,
      nearest_anchor_borel α p hp z = nearest_anchor_borel α p hp z₀ := by
  set s := nearest_anchor_borel α p hp z₀
  have h_strict := fatou_strict_min α p hp z₀ hz₀
  -- Continuité de z ↦ ‖z − γ_t‖ − ‖z − γ_s‖.
  have h_cont : ∀ t : Fin p, Continuous (fun z : ℂ =>
      ‖z - cycAnchor α p t‖ - ‖z - cycAnchor α p s‖) := fun t =>
    (continuous_norm.comp (continuous_id.sub continuous_const)).sub
    (continuous_norm.comp (continuous_id.sub continuous_const))
  -- Pour chaque t, eventually `t = s ∨ ‖z − γ_s‖ < ‖z − γ_t‖`.
  have h_per_t : ∀ t : Fin p, ∀ᶠ z in 𝓝 z₀,
      t = s ∨ ‖z - cycAnchor α p s‖ < ‖z - cycAnchor α p t‖ := by
    intro t
    by_cases ht : t = s
    · exact Filter.eventually_of_forall (fun _ => Or.inl ht)
    · -- Au point z₀ : gap strict.
      have h_pos : 0 < ‖z₀ - cycAnchor α p t‖ - ‖z₀ - cycAnchor α p s‖ := by
        have h := h_strict t ht
        linarith
      -- Ouvert par continuité.
      have h_open : IsOpen {z : ℂ |
          0 < ‖z - cycAnchor α p t‖ - ‖z - cycAnchor α p s‖} :=
        isOpen_lt continuous_const (h_cont t)
      have h_mem : z₀ ∈ {z : ℂ |
          0 < ‖z - cycAnchor α p t‖ - ‖z - cycAnchor α p s‖} := h_pos
      filter_upwards [h_open.mem_nhds h_mem] with z hz
      exact Or.inr (by linarith [hz])
  -- Finite intersection via Filter.eventually_all.
  have h_all : ∀ᶠ z in 𝓝 z₀, ∀ t : Fin p,
      t = s ∨ ‖z - cycAnchor α p s‖ < ‖z - cycAnchor α p t‖ := by
    exact (Filter.eventually_all).mpr h_per_t
  -- Conclure : sur le voisinage, argmin = {s}, donc nearest = s.
  filter_upwards [h_all] with z hz
  unfold nearest_anchor_borel
  have h_argmin_singleton : argminAnchorFinset α p z = {s} := by
    ext k
    unfold argminAnchorFinset
    simp only [Finset.mem_filter, Finset.mem_singleton, Finset.mem_univ, true_and]
    constructor
    · intro h_k_min
      by_contra h_ne_s
      have h_gt : ‖z - cycAnchor α p s‖ < ‖z - cycAnchor α p k‖ := by
        rcases hz k with h_eq | h_lt
        · exact absurd h_eq h_ne_s
        · exact h_lt
      have h_le := h_k_min s
      linarith
    · intro h_eq
      subst h_eq
      intro u
      rcases hz u with h_u_eq | h_u_lt
      · rw [h_u_eq]
      · exact le_of_lt h_u_lt
  -- min' ∈ argminAnchorFinset = {s}, donc min' = s.
  set m := (argminAnchorFinset α p z).min' (argminAnchorFinset_nonempty α p hp z)
  have hm_mem : m ∈ argminAnchorFinset α p z :=
    Finset.min'_mem _ _
  have hm_in_sing : m ∈ ({s} : Finset (Fin p)) := by
    rw [← h_argmin_singleton]; exact hm_mem
  exact Finset.mem_singleton.mp hm_in_sing

end LocallyConstant

/-! ============================================================
  §105.3  `nearest_anchor_borel` continuousAt sur Fatou
============================================================ -/

section ContinuousAtFatou

/-- **★★ `nearest_anchor_borel` est `ContinuousAt` sur Fatou.** Suit
    de la locally constant par §105.2, puisque localement constant
    implique continuité (via `Tendsto.congr'`). -/
theorem nearest_anchor_borel_continuousAt_fatou
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z₀ : ℂ)
    (hz₀ : z₀ ∈ PandrosionFatouSet α p) :
    ContinuousAt (nearest_anchor_borel α p hp) z₀ := by
  have h_const := nearest_anchor_borel_eventually_const α p hp z₀ hz₀
  -- Convert `∀ᶠ z, f z = c` to `EventuallyEq`.
  have h_ee : (nearest_anchor_borel α p hp) =ᶠ[𝓝 z₀]
              (fun _ => nearest_anchor_borel α p hp z₀) := h_const
  unfold ContinuousAt
  have h_tendsto_const : Tendsto (fun _ : ℂ => nearest_anchor_borel α p hp z₀)
      (𝓝 z₀) (𝓝 (nearest_anchor_borel α p hp z₀)) := tendsto_const_nhds
  exact h_tendsto_const.congr' h_ee.symm

end ContinuousAtFatou

/-! ============================================================
  §105.4  `niveau5_multi_start_target_borel` continuousAt sur Fatou
============================================================ -/

section Niveau5BorelContinuousAt

/-- **★★ `niveau5_multi_start_target_borel` est `ContinuousAt` sur
    Fatou.** Composition de `cycAnchor` (continue, domaine discret)
    et `nearest_anchor_borel` (ContinuousAt sur Fatou). -/
theorem niveau5_multi_start_target_borel_continuousAt_fatou
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z₀ : ℂ)
    (hz₀ : z₀ ∈ PandrosionFatouSet α p) :
    ContinuousAt (niveau5_multi_start_target_borel α p hp) z₀ := by
  unfold niveau5_multi_start_target_borel
  -- Fin p a topologie discrète ⇒ toute fonction de Fin p est continue.
  have h_cycAnchor_cont : Continuous (cycAnchor α p) :=
    continuous_of_discreteTopology
  exact h_cycAnchor_cont.continuousAt.comp
    (nearest_anchor_borel_continuousAt_fatou α p hp z₀ hz₀)

end Niveau5BorelContinuousAt

/-! ============================================================
  §105.5  ★★★★★★★★★★ Théorème principal : a.e. continuous
============================================================ -/

section AEContinuous

/-- **★★★★★★★★★★ CONTINUITÉ A.E. DE LA FONCTION SOLUTION PANDROSION.**

    Pour tout `(α, p)` avec `α ≠ 0`, la fonction solution multi-start
    `niveau5_multi_start_target_borel α p hp : ℂ → ℂ` est
    `ContinuousAt` à *Lebesgue-presque tout* `z ∈ ℂ`.

    **Caractérisation complète (Luzin-style).** Cette fonction est :

      • **Mesurable partout** (§102 `niveau5_multi_start_target_borel_measurable`).
      • **Continue à Lebesgue-presque tout `z`** (présent théorème).
      • **Constante localement** sur le Fatou set (`∀ z ∈ Fatou`,
        `nearest_anchor_borel` est localement constante).

    **Conséquence analytique.** La décomposition `ℂ = Fatou ∪ Julia` :
      • Sur Fatou (ouvert, de mesure pleine) : la fonction est
        constante par morceaux, continue partout.
      • Sur Julia (fermé, de mesure zéro) : la fonction peut être
        discontinue, mais c'est un ensemble Lebesgue-null (§103). -/
theorem niveau5_multi_start_target_borel_ae_continuous
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (hα : α ≠ 0) :
    ∀ᵐ z : ℂ ∂volume,
      ContinuousAt (niveau5_multi_start_target_borel α p hp) z := by
  have h_fatou_ae := pandrosion_fatou_ae_full_measure α p hα hp
  filter_upwards [h_fatou_ae] with z hz
  exact niveau5_multi_start_target_borel_continuousAt_fatou α p hp z hz

end AEContinuous

end Pandrosion
