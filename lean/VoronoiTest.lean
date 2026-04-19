import Mathlib.Data.Complex.Basic
import Mathlib.Tactic

open Complex

theorem voronoi_halfplane_affine (r1 r2 z : ℂ) :
    Complex.normSq (z - r1) ≤ Complex.normSq (z - r2) ↔
    2 * (z.re * (r2.re - r1.re) + z.im * (r2.im - r1.im)) ≤
      Complex.normSq r2 - Complex.normSq r1 := by
  have h1 : Complex.normSq (z - r1) = (z.re - r1.re)^2 + (z.im - r1.im)^2 := rfl
  have h2 : Complex.normSq (z - r2) = (z.re - r2.re)^2 + (z.im - r2.im)^2 := rfl
  have h3 : Complex.normSq r1 = r1.re^2 + r1.im^2 := rfl
  have h4 : Complex.normSq r2 = r2.re^2 + r2.im^2 := rfl
  rw [h1, h2, h3, h4]
  constructor
  · intro h; ring_nf at h ⊢; linarith
  · intro h; ring_nf at h ⊢; linarith
