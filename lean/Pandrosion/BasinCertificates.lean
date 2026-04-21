/-
  Universitas Pandrosion — Lean 4 Formalization
  ADVANCED² — Non-Trivial Extensions (Aggregator)

  Post-split: the substantive lemmas of §I–§XI, §XIII have been
  relocated to `Pandrosion.Basin.BabylonianContraction`,
  `Pandrosion.Basin.PellCone`, and `Pandrosion.Basin.Voronoi`.

  This file re-exports those modules and keeps the two grand
  certificates (§V `advanced2_grand_certificate`, §XII
  `advanced2_super_certificate`) which reference results across
  all three.
-/
import Pandrosion.Basin.BabylonianContraction
import Pandrosion.Basin.PellCone
import Pandrosion.Basin.Voronoi

namespace Pandrosion

section Advanced2

/-! ## §V. Unified Advanced² Certificate -/

/-- **Unified Advanced² certificate.**
    Four unconditional facts, each non-trivial, assembled from
    §I, §II, §III, §IV. -/
theorem advanced2_grand_certificate
    (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ)
    (A B X : ℤ) (hX : 0 ≤ X) (h_cone : A ^ 2 ≥ X * B ^ 2)
    (r1x r1y r2x r2y : ℝ) :
    -- (I) Unconditional convergence:
    (error (iterate 2 x s₀ n) r ≤ (1 / 2) ^ n * error s₀ r) ∧
    -- (II) Pell-cone preservation:
    ((pell_iterate X (A, B) n).1 ^ 2 ≥ X * (pell_iterate X (A, B) n).2 ^ 2) ∧
    -- (III) Voronoï openness:
    (IsOpen {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2}) ∧
    -- (IV) Lipschitz-½ contraction on the basin:
    (error (pandrosion_map 2 x s₀) r ≤ (1 / 2) * error s₀ r) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · exact global_convergence_p2 x r s₀ hx hr h_s₀ h_root n
  · exact pell_cone_invariant X hX (A, B) h_cone n
  · exact voronoi_strict_open r1x r1y r2x r2y
  · exact pandrosion_map_lipschitz_half_p2 x r s₀ hx hr h_root h_s₀

/-! ## §XII. Super-Certificate — Ten Unconditional Facts -/

/-- **Super-certificate.** Packages every non-trivial fact in
    Advanced2 into one compound theorem. -/
theorem advanced2_super_certificate
    (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ)
    (A B X : ℤ) (hX : 1 ≤ X) (hA : 1 ≤ A) (hB : 1 ≤ B)
    (h_cone_strict : A ^ 2 > X * B ^ 2)
    (r1x r1y r2x r2y : ℝ) (anchors : Finset (ℝ × ℝ))
    (r0 : ℝ × ℝ) (s t : ℝ) (hs : s ≥ r) (ht : t ≥ r) :
    -- (I) Unconditional convergence:
    (error (iterate 2 x s₀ n) r ≤ (1 / 2) ^ n * error s₀ r) ∧
    -- (II) Pell-cone preservation (weak):
    ((pell_iterate X (A, B) n).1 ^ 2 ≥ X * (pell_iterate X (A, B) n).2 ^ 2) ∧
    -- (III) Voronoï openness (pairwise):
    (IsOpen {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2}) ∧
    -- (IV) Lipschitz-½ to fixed point:
    (error (pandrosion_map 2 x s₀) r ≤ (1 / 2) * error s₀ r) ∧
    -- (VI) Monotone descent:
    (iterate 2 x s₀ (n + 1) ≤ iterate 2 x s₀ n) ∧
    -- (VII) Quadratic convergence:
    (|pandrosion_map 2 x s₀ - r| ≤ (s₀ - r) ^ 2 / (2 * r)) ∧
    -- (VIII) Strict Pell-cone preservation:
    ((pell_iterate X (A, B) n).1 ^ 2 > X * (pell_iterate X (A, B) n).2 ^ 2) ∧
    -- (IX) Finite Voronoï openness:
    (IsOpen {p : ℝ × ℝ | ∀ a ∈ anchors,
      (p.1 - r0.1) ^ 2 + (p.2 - r0.2) ^ 2 <
      (p.1 - a.1) ^ 2 + (p.2 - a.2) ^ 2}) ∧
    -- (X) Pell orbit strict growth:
    ((pell_iterate X (A, B) n).1 < (pell_iterate X (A, B) (n + 1)).1) ∧
    -- (XI) Pairwise Banach contraction:
    (|pandrosion_map 2 x s - pandrosion_map 2 x t| ≤ (1 / 2) * |s - t|) := by
  have hX0 : 0 ≤ X := by linarith
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact global_convergence_p2 x r s₀ hx hr h_s₀ h_root n
  · exact pell_cone_invariant X hX0 (A, B) (by
      show A ^ 2 ≥ X * B ^ 2; linarith) n
  · exact voronoi_strict_open r1x r1y r2x r2y
  · exact pandrosion_map_lipschitz_half_p2 x r s₀ hx hr h_root h_s₀
  · exact iterate_monotone_descent_p2 x r s₀ hr h_s₀ h_root n
  · exact pandrosion_quadratic_convergence_p2 x r s₀ hr (by linarith) h_root h_s₀
  · exact pell_cone_invariant_strict X hX0 (A, B) h_cone_strict n
  · exact voronoi_cell_finite_open r0 anchors
  · exact (pell_iterate_strict_growth X hX (A, B) hA hB n).1
  · exact pandrosion_map_lipschitz_pairwise_p2 x r s t hr h_root hs ht

end Advanced2

end Pandrosion
