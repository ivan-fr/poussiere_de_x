/-
  Universitas Pandrosion — §102. **★★★★★★★★★ Mesurabilité de la
  fonction solution multi-start Pandrosion.**

  §99 définit `niveau5_multi_start_target : ℂ → ℂ` via `Classical.choose`
  (dans `nearest_anchor_selector`), donc non-mesurable a priori.

  §102 construit une version *déterministe* `niveau5_multi_start_target_borel`
  via `Finset.min'` (choix de l'indice minimal parmi les équidistants),
  et prouve sa Borel-mesurabilité :

    `Measurable (niveau5_multi_start_target_borel α p hp)`

  Ouvre la porte aux théorèmes d'intégration (E[‖σⁿ(z) − α‖], etc.)
  et aux probabilités sur les trajectoires Pandrosion.

  Contents.

    §102.1  `argminAnchorFinset`, `nearest_anchor_borel` — sélecteur
            déterministe via `Finset.min'`.
    §102.2  `halfplane_measurableSet` — les half-planes Voronoï sont
            Borel.
    §102.3  `argmin_mem_set_measurable` + `nearest_anchor_borel_preimage_measurable`
            — les level-sets sont Borel.
    §102.4  `nearest_anchor_borel_measurable` — measurability.
    §102.5  `niveau5_multi_start_target_borel_measurable` — théorème
            principal.
-/

import Pandrosion.Core.Niveau5StrictMultiStart

namespace Pandrosion

open Complex MeasureTheory Filter Topology

/-! ============================================================
  §102.1  Sélecteur déterministe via `Finset.min'`
============================================================ -/

section DeterministicSelector

/-- **Argmin Finset.** L'ensemble des indices `k ∈ Fin p` minimisant
    `‖z − γ_k‖`. -/
noncomputable def argminAnchorFinset (α : ℂ) (p : ℕ) (z : ℂ) : Finset (Fin p) :=
  (Finset.univ : Finset (Fin p)).filter (fun k : Fin p =>
    ∀ t : Fin p, ‖z - cycAnchor α p k‖ ≤ ‖z - cycAnchor α p t‖)

/-- **L'argmin finset est non-vide** (par `Finset.exists_min_image`). -/
theorem argminAnchorFinset_nonempty (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z : ℂ) :
    (argminAnchorFinset α p z).Nonempty := by
  unfold argminAnchorFinset
  obtain ⟨k₀, _, hk₀_min⟩ := (Finset.univ : Finset (Fin p)).exists_min_image
      (fun k : Fin p => ‖z - cycAnchor α p k‖)
      (Finset.univ_nonempty_iff.mpr ⟨⟨0, hp⟩⟩)
  refine ⟨k₀, ?_⟩
  rw [Finset.mem_filter]
  exact ⟨Finset.mem_univ _, fun t => hk₀_min t (Finset.mem_univ _)⟩

/-- **Sélecteur Borel-mesurable** : retourne le plus petit indice
    parmi les argmin. Déterministe. -/
noncomputable def nearest_anchor_borel (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z : ℂ) : Fin p :=
  (argminAnchorFinset α p z).min' (argminAnchorFinset_nonempty α p hp z)

/-- **Version Borel de `niveau5_multi_start_target`.** -/
noncomputable def niveau5_multi_start_target_borel
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z : ℂ) : ℂ :=
  cycAnchor α p (nearest_anchor_borel α p hp z)

/-- **Spec du sélecteur Borel** : il minimise toujours `‖z − γ_k‖`. -/
theorem nearest_anchor_borel_spec (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z : ℂ) :
    ∀ t : Fin p,
      ‖z - cycAnchor α p (nearest_anchor_borel α p hp z)‖
        ≤ ‖z - cycAnchor α p t‖ := by
  intro t
  have h_mem : nearest_anchor_borel α p hp z ∈ argminAnchorFinset α p z :=
    Finset.min'_mem _ (argminAnchorFinset_nonempty α p hp z)
  unfold argminAnchorFinset at h_mem
  rw [Finset.mem_filter] at h_mem
  exact h_mem.2 t

end DeterministicSelector

/-! ============================================================
  §102.2  Half-planes sont Borel
============================================================ -/

section HalfplanesMeasurable

/-- **Half-plane `‖z − γ_s‖ ≤ ‖z − γ_t‖` est fermé (et donc Borel).** -/
theorem halfplane_measurableSet (α : ℂ) (p : ℕ) (s t : Fin p) :
    MeasurableSet {z : ℂ | ‖z - cycAnchor α p s‖ ≤ ‖z - cycAnchor α p t‖} := by
  have hf : Continuous (fun z : ℂ => ‖z - cycAnchor α p s‖) :=
    continuous_norm.comp (continuous_id.sub continuous_const)
  have hg : Continuous (fun z : ℂ => ‖z - cycAnchor α p t‖) :=
    continuous_norm.comp (continuous_id.sub continuous_const)
  exact (isClosed_le hf hg).measurableSet

end HalfplanesMeasurable

/-! ============================================================
  §102.3  Level-sets Borel
============================================================ -/

section LevelSetsBorel

/-- **Ensemble `{z | s ∈ argmin(z)}` est Borel** (intersection finie
    de half-planes fermés). -/
theorem argmin_mem_set_measurable (α : ℂ) (p : ℕ) (s : Fin p) :
    MeasurableSet {z : ℂ | s ∈ argminAnchorFinset α p z} := by
  have h_char : {z : ℂ | s ∈ argminAnchorFinset α p z} =
                ⋂ t : Fin p, {z : ℂ | ‖z - cycAnchor α p s‖ ≤ ‖z - cycAnchor α p t‖} := by
    ext z
    unfold argminAnchorFinset
    simp [Finset.mem_filter]
  rw [h_char]
  exact MeasurableSet.iInter (fun t => halfplane_measurableSet α p s t)

/-- **Level-set `{z | nearest_anchor_borel z = s}` est Borel**
    (intersection de `{s ∈ argmin}` et finite intersection de
    complements `{k ∉ argmin}` pour `k < s`). -/
theorem nearest_anchor_borel_preimage_measurable
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (s : Fin p) :
    MeasurableSet {z : ℂ | nearest_anchor_borel α p hp z = s} := by
  have h_char : {z : ℂ | nearest_anchor_borel α p hp z = s} =
                {z : ℂ | s ∈ argminAnchorFinset α p z} ∩
                ⋂ k : Fin p, ⋂ (_ : k < s),
                  {z : ℂ | k ∉ argminAnchorFinset α p z} := by
    ext z
    simp only [Set.mem_setOf_eq, Set.mem_inter_iff, Set.mem_iInter]
    constructor
    · intro h_min
      refine ⟨?_, ?_⟩
      · rw [← h_min]
        exact Finset.min'_mem _ _
      · intro k hks h_mem
        have h_le : (argminAnchorFinset α p z).min' _ ≤ k :=
          Finset.min'_le _ k h_mem
        unfold nearest_anchor_borel at h_min
        rw [h_min] at h_le
        exact absurd h_le (not_le.mpr hks)
    · intro ⟨h_s_in, h_smaller⟩
      apply le_antisymm
      · unfold nearest_anchor_borel
        exact Finset.min'_le _ _ h_s_in
      · by_contra h_lt
        push_neg at h_lt
        unfold nearest_anchor_borel at h_lt
        have h_min_mem := Finset.min'_mem
          (argminAnchorFinset α p z) (argminAnchorFinset_nonempty α p hp z)
        exact h_smaller _ h_lt h_min_mem
  rw [h_char]
  apply MeasurableSet.inter (argmin_mem_set_measurable α p s)
  apply MeasurableSet.iInter
  intro k
  apply MeasurableSet.iInter
  intro _
  exact (argmin_mem_set_measurable α p k).compl

end LevelSetsBorel

/-! ============================================================
  §102.4  `nearest_anchor_borel` est Borel-mesurable
============================================================ -/

section NearestAnchorMeasurable

/-- **★★ Le sélecteur `nearest_anchor_borel α p hp : ℂ → Fin p` est
    Borel-mesurable.** -/
theorem nearest_anchor_borel_measurable (α : ℂ) (p : ℕ) (hp : 1 ≤ p) :
    Measurable (nearest_anchor_borel α p hp) := by
  intro T _hT
  have h_decomp : nearest_anchor_borel α p hp ⁻¹' T =
                  ⋃ s ∈ T, {z : ℂ | nearest_anchor_borel α p hp z = s} := by
    ext z
    simp only [Set.mem_preimage, Set.mem_iUnion, Set.mem_setOf_eq]
    exact ⟨fun hz => ⟨_, hz, rfl⟩, fun ⟨_, hs_T, hs_eq⟩ => hs_eq ▸ hs_T⟩
  rw [h_decomp]
  apply MeasurableSet.biUnion T.toFinite.countable
  intro s _
  exact nearest_anchor_borel_preimage_measurable α p hp s

end NearestAnchorMeasurable

/-! ============================================================
  §102.5  ★★★ Théorème principal : `niveau5_multi_start_target_borel`
              est Borel-mesurable
============================================================ -/

section MainTheorem

/-- **Mesurabilité de `cycAnchor α p : Fin p → ℂ`** : triviale car
    `Fin p` est fini discret. -/
theorem cycAnchor_measurable (α : ℂ) (p : ℕ) :
    Measurable (cycAnchor α p) := by
  intro s _hs
  exact (Set.toFinite _).measurableSet

/-- **★★★★★★★★★ Mesurabilité de la fonction solution Pandrosion.**

    `niveau5_multi_start_target_borel α p hp : ℂ → ℂ` est Borel-
    mesurable. C'est le théorème qui permet de parler de
    probabilités, d'intégrales, d'espérances sur les trajectoires
    Pandrosion. -/
theorem niveau5_multi_start_target_borel_measurable
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) :
    Measurable (niveau5_multi_start_target_borel α p hp) := by
  unfold niveau5_multi_start_target_borel
  exact (cycAnchor_measurable α p).comp (nearest_anchor_borel_measurable α p hp)

end MainTheorem

end Pandrosion
