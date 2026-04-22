/-
  Universitas Pandrosion — Legacy / Voronoï affine form.

  Two legacy facts complementing `Pandrosion.Core.MultiStartBasins`:

    • `voronoi_halfplane_affine` — the Voronoï bisector inequality
      `|z − r₁|² ≤ |z − r₂|²` is *affine* in `z`. This is the
      coordinate-level statement behind the convexity lemma
      `halfPlane_convex` in `MultiStartBasins`.
    • `multistart_distinct_roots_guarantee` — on `Fin d`, any
      surjection `σ : Fin d → Fin d` is bijective. This complements
      `fin_injective_iff_bijective` in `MultiStartBasins`, giving the
      other pigeonhole direction (surjective ⇒ bijective).

  Ported from `MultiStart.lean` + `VoronoiInvariance.lean` on the
  `main_previous` branch.
-/

import Mathlib.Data.Complex.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §L18. The Voronoï bisector as an affine inequality -/

/-- **Voronoï bisector = affine half-space.**
    `|z − r₁|² ≤ |z − r₂|²` is equivalent to the affine inequality
    `2⟨z, r₂ − r₁⟩ ≤ |r₂|² − |r₁|²` (coordinate form). This is the
    algebraic identity underlying `MultiStartBasins.halfPlane_convex`. -/
theorem voronoi_halfplane_affine (r1 r2 z : ℂ) :
    Complex.normSq (z - r1) ≤ Complex.normSq (z - r2) ↔
    2 * (z.re * (r2.re - r1.re) + z.im * (r2.im - r1.im)) ≤
      Complex.normSq r2 - Complex.normSq r1 := by
  change (z.re - r1.re) * (z.re - r1.re) + (z.im - r1.im) * (z.im - r1.im) ≤
         (z.re - r2.re) * (z.re - r2.re) + (z.im - r2.im) * (z.im - r2.im) ↔
    2 * (z.re * (r2.re - r1.re) + z.im * (r2.im - r1.im)) ≤
      (r2.re * r2.re + r2.im * r2.im) - (r1.re * r1.re + r1.im * r1.im)
  constructor
  · intro h; ring_nf at h ⊢; linarith
  · intro h; ring_nf at h ⊢; linarith

/-! ## §L19. Surjective ⇒ bijective on `Fin d` -/

/-- **Multi-start coverage guarantee (surjective direction).**
    Under a basin-assignment `σ : Fin d → Fin d` that covers every
    root (surjective), the assignment is automatically bijective.
    Combined with `MultiStartBasins.fin_injective_iff_bijective` this
    gives pigeonhole in both directions. -/
theorem multistart_distinct_roots_guarantee (d : ℕ)
    (basin_assignment : Fin d → Fin d)
    (h_coverage : Function.Surjective basin_assignment) :
    Function.Bijective basin_assignment :=
  (Fintype.bijective_iff_surjective_and_card basin_assignment).mpr
    ⟨h_coverage, rfl⟩

end Pandrosion
