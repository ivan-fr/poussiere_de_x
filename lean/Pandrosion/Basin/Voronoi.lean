/-
  Universitas Pandrosion — Basin / Voronoi
  Split from Advanced2.lean: §III, IX
  (Strict Voronoï openness, finite intersection)
-/
import Mathlib.Topology.MetricSpace.Basic
import Mathlib.Topology.Instances.Real
import Mathlib.Topology.Order.Basic
import Mathlib.Tactic

namespace Pandrosion

section Voronoi

/-! ## §III. The Strict Voronoï Interior is Open -/

/-- **The squared-distance difference is continuous in z.** -/
theorem squared_distance_diff_continuous (r1x r1y r2x r2y : ℝ) :
    Continuous (fun p : ℝ × ℝ =>
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 -
      ((p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2)) := by
  have c1 : Continuous (fun p : ℝ × ℝ => p.1) := continuous_fst
  have c2 : Continuous (fun p : ℝ × ℝ => p.2) := continuous_snd
  exact ((c1.sub continuous_const).pow 2).add ((c2.sub continuous_const).pow 2)
        |>.sub (((c1.sub continuous_const).pow 2).add ((c2.sub continuous_const).pow 2))

/-- **The strict Voronoï half-plane is open in ℝ².** -/
theorem voronoi_strict_open (r1x r1y r2x r2y : ℝ) :
    IsOpen {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2} := by
  let D : ℝ × ℝ → ℝ := fun p =>
    (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 -
    ((p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2)
  have hcont : Continuous D := squared_distance_diff_continuous r1x r1y r2x r2y
  have h_eq : {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2} = D ⁻¹' Set.Iio 0 := by
    ext p
    simp only [D, Set.mem_setOf_eq, Set.mem_preimage, Set.mem_Iio]
    constructor
    · intro h; linarith
    · intro h; linarith
  rw [h_eq]
  exact isOpen_Iio.preimage hcont

/-- **Explicit δ-neighborhood witness.** -/
theorem voronoi_strict_delta_exists (r1x r1y r2x r2y z1 z2 : ℝ)
    (h : (z1 - r1x) ^ 2 + (z2 - r1y) ^ 2 <
         (z1 - r2x) ^ 2 + (z2 - r2y) ^ 2) :
    ∃ δ > 0, ∀ w1 w2 : ℝ, |w1 - z1| < δ → |w2 - z2| < δ →
      (w1 - r1x) ^ 2 + (w2 - r1y) ^ 2 <
      (w1 - r2x) ^ 2 + (w2 - r2y) ^ 2 := by
  have hopen := voronoi_strict_open r1x r1y r2x r2y
  have hmem : (z1, z2) ∈ {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2} := by
    show (z1 - r1x) ^ 2 + (z2 - r1y) ^ 2 < (z1 - r2x) ^ 2 + (z2 - r2y) ^ 2
    exact h
  obtain ⟨δ, hδ_pos, hδ⟩ := Metric.isOpen_iff.mp hopen (z1, z2) hmem
  refine ⟨δ, hδ_pos, ?_⟩
  intro w1 w2 hw1 hw2
  have hball : (w1, w2) ∈ Metric.ball (z1, z2) δ := by
    rw [Metric.mem_ball, Prod.dist_eq, max_lt_iff]
    refine ⟨?_, ?_⟩
    · simpa [Real.dist_eq] using hw1
    · simpa [Real.dist_eq] using hw2
  exact hδ hball

/-! ## §IX. Finite Voronoï Intersection -/

/-- **Finite Voronoï cell is open in ℝ².** -/
theorem voronoi_cell_finite_open (r0 : ℝ × ℝ) (anchors : Finset (ℝ × ℝ)) :
    IsOpen {p : ℝ × ℝ | ∀ a ∈ anchors,
      (p.1 - r0.1) ^ 2 + (p.2 - r0.2) ^ 2 <
      (p.1 - a.1) ^ 2 + (p.2 - a.2) ^ 2} := by
  have h_eq : {p : ℝ × ℝ | ∀ a ∈ anchors,
      (p.1 - r0.1) ^ 2 + (p.2 - r0.2) ^ 2 <
      (p.1 - a.1) ^ 2 + (p.2 - a.2) ^ 2} =
      ⋂ a ∈ (anchors : Set (ℝ × ℝ)), {p : ℝ × ℝ |
        (p.1 - r0.1) ^ 2 + (p.2 - r0.2) ^ 2 <
        (p.1 - a.1) ^ 2 + (p.2 - a.2) ^ 2} := by
    ext p
    simp only [Set.mem_setOf_eq, Set.mem_iInter, Finset.mem_coe]
  rw [h_eq]
  refine Set.Finite.isOpen_biInter (Finset.finite_toSet _) ?_
  intro a _
  exact voronoi_strict_open r0.1 r0.2 a.1 a.2

end Voronoi

end Pandrosion
