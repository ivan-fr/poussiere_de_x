/-
  Universitas Pandrosion — §18. Lebesgue-null Voronoï boundary.

  The `constructive_mcmullen` theorem (Foundations §13) establishes
  uniqueness of the multi-start selector *off* the Voronoï boundary
  `VoronoiBoundary γ`. This file promotes that statement from "off
  the boundary" to "almost everywhere" by proving that the Voronoï
  boundary has 2-dimensional Lebesgue measure zero.

  The strategy:

    §18.1  Each perpendicular bisector between two *distinct* anchors
           `a ≠ b` is a strict affine subspace of `ℂ` (an affine line),
           hence `volume (perpBisector a b) = 0`.

    §18.2  `VoronoiBoundary γ` is contained in the finite union of all
           such bisectors (Foundations §13 lemma
           `VoronoiBoundary_subset_bisectors`). A finite union of
           measure-zero sets is measure-zero.

    ★ §18.3  `voronoiBoundary_volume_zero`
           For any injective family of anchors `γ : Fin p → ℂ`,
           `volume (VoronoiBoundary γ) = 0`.

  Combined with `constructive_mcmullen`, this yields the "almost-
  everywhere" global-convergence corollary
  `constructive_mcmullen_ae`.
-/

import Pandrosion.Core.MultiStartBasins
import Mathlib.MeasureTheory.Measure.Lebesgue.EqHaar
import Mathlib.MeasureTheory.Measure.Lebesgue.Complex

namespace Pandrosion

open Complex MeasureTheory

/-! ============================================================
  §18.1  Each perpendicular bisector has Lebesgue measure zero
============================================================ -/

section PerpBisectorMeasure

/-- **Real-linear form** whose kernel is the translate of
    `perpBisector a b`. Explicitly,
    `L z = (b − a).re · z.re + (b − a).im · z.im`,
    i.e. twice the real inner product of `b − a` and `z`. -/
noncomputable def perpBisectorForm (a b : ℂ) : ℂ →ₗ[ℝ] ℝ where
  toFun z := (b - a).re * z.re + (b - a).im * z.im
  map_add' x y := by
    simp only [Complex.add_re, Complex.add_im]; ring
  map_smul' r z := by
    show (b - a).re * (r • z).re + (b - a).im * (r • z).im
       = r • ((b - a).re * z.re + (b - a).im * z.im)
    simp only [Complex.real_smul, Complex.mul_re, Complex.mul_im,
               Complex.ofReal_re, Complex.ofReal_im, smul_eq_mul]
    ring

@[simp]
theorem perpBisectorForm_apply (a b z : ℂ) :
    perpBisectorForm a b z = (b - a).re * z.re + (b - a).im * z.im := rfl

/-- `a ≠ b` implies the real-linear form `perpBisectorForm a b`
    is non-zero (so its kernel is a strict submodule). -/
theorem perpBisectorForm_ne_zero_of_ne {a b : ℂ} (hab : a ≠ b) :
    perpBisectorForm a b ≠ 0 := by
  intro h
  apply hab
  have hre : (b - a).re = 0 := by
    have h1 : perpBisectorForm a b 1 = 0 := by rw [h]; rfl
    simpa [perpBisectorForm_apply, Complex.one_re, Complex.one_im] using h1
  have him : (b - a).im = 0 := by
    have h2 : perpBisectorForm a b Complex.I = 0 := by rw [h]; rfl
    simpa [perpBisectorForm_apply, Complex.I_re, Complex.I_im] using h2
  have hba : b - a = 0 := Complex.ext hre him
  exact (sub_eq_zero.mp hba).symm

/-- **Perpendicular bisector measure zero.** For two distinct
    anchors `a ≠ b ∈ ℂ`, the set `{z : ‖a − z‖ = ‖b − z‖}` is a
    Lebesgue-null subset of `ℂ`.

    **Proof.** Squaring the defining identity `‖a − z‖ = ‖b − z‖`
    and expanding in real coordinates yields the affine equation
    `L z = (‖b‖² − ‖a‖²)/2`, where `L = perpBisectorForm a b` is
    a non-zero `ℝ`-linear form (§18.1). Hence the bisector is
    contained in the affine subspace `(a+b)/2 + ker L`, which is a
    strict affine subspace of `ℂ` and therefore has Haar-measure
    zero (`MeasureTheory.Measure.addHaar_affineSubspace`). -/
theorem perpBisector_volume_zero (a b : ℂ) (hab : a ≠ b) :
    volume (perpBisector a b) = 0 := by
  -- Midpoint of a and b lies on the bisector (and is the anchor of our affine model).
  let z₀ : ℂ := (a + b) / 2
  have hz₀_re : z₀.re = (a.re + b.re) / 2 := by
    show ((a + b) / 2).re = _
    simp [Complex.add_re, Complex.div_re, Complex.ofReal_re, Complex.ofReal_im,
          Complex.normSq, Complex.mul_re]
  have hz₀_im : z₀.im = (a.im + b.im) / 2 := by
    show ((a + b) / 2).im = _
    simp [Complex.add_im, Complex.div_im, Complex.ofReal_re, Complex.ofReal_im,
          Complex.normSq, Complex.mul_im]
  -- Non-zero linear form L and the affine subspace `z₀ + ker L`.
  let L : ℂ →ₗ[ℝ] ℝ := perpBisectorForm a b
  have hL_ne : L ≠ 0 := perpBisectorForm_ne_zero_of_ne hab
  have hker_ne_top : (LinearMap.ker L) ≠ ⊤ := by
    rw [Ne, LinearMap.ker_eq_top]
    exact hL_ne
  -- Step 1: perpBisector a b ⊆ (mk' z₀ (ker L) : Set ℂ).
  have h_subset :
      perpBisector a b ⊆
        ((AffineSubspace.mk' z₀ (LinearMap.ker L) : AffineSubspace ℝ ℂ) : Set ℂ) := by
    intro z hz
    have hz_eq : ‖a - z‖ = ‖b - z‖ := hz
    have hz_sq : ‖a - z‖ ^ 2 = ‖b - z‖ ^ 2 := by rw [hz_eq]
    have h_nm : ∀ w : ℂ, ‖w‖ ^ 2 = w.re ^ 2 + w.im ^ 2 := fun w => by
      rw [Complex.norm_eq_abs, Complex.sq_abs, Complex.normSq_apply]; ring
    rw [h_nm, h_nm] at hz_sq
    simp only [Complex.sub_re, Complex.sub_im] at hz_sq
    show z ∈ AffineSubspace.mk' z₀ (LinearMap.ker L)
    rw [AffineSubspace.mem_mk'_iff_vsub_mem, LinearMap.mem_ker]
    show (b - a).re * (z -ᵥ z₀).re + (b - a).im * (z -ᵥ z₀).im = 0
    show (b - a).re * (z - z₀).re + (b - a).im * (z - z₀).im = 0
    simp only [Complex.sub_re, Complex.sub_im, hz₀_re, hz₀_im]
    linear_combination (1/2 : ℝ) * hz_sq
  -- Step 2: mk' z₀ (ker L) ≠ ⊤.
  have hS_ne_top :
      (AffineSubspace.mk' z₀ (LinearMap.ker L) : AffineSubspace ℝ ℂ) ≠ ⊤ := by
    intro h_top
    apply hker_ne_top
    have hdir :
        (AffineSubspace.mk' z₀ (LinearMap.ker L) : AffineSubspace ℝ ℂ).direction = ⊤ := by
      rw [h_top]
      exact (AffineSubspace.direction_eq_top_iff_of_nonempty
              (s := (⊤ : AffineSubspace ℝ ℂ)) ⟨z₀, trivial⟩).mpr rfl
    rw [AffineSubspace.direction_mk'] at hdir
    exact hdir
  -- Step 3: monotonicity + addHaar_affineSubspace.
  have h_measure_S :
      volume ((AffineSubspace.mk' z₀ (LinearMap.ker L) : AffineSubspace ℝ ℂ) : Set ℂ) = 0 :=
    MeasureTheory.Measure.addHaar_affineSubspace volume _ hS_ne_top
  exact measure_mono_null h_subset h_measure_S

end PerpBisectorMeasure

/-! ============================================================
  §18.2  Voronoï boundary has Lebesgue measure zero
============================================================ -/

section VoronoiBoundaryMeasure

/-- **★ Voronoï boundary is Lebesgue-null.**
    For any injective family of anchors `γ : Fin p → ℂ`, the
    ambiguous-nearest-anchor set `VoronoiBoundary γ` has
    2-dimensional Lebesgue measure zero.

    *Proof.* The boundary lies inside the finite union of
    perpendicular bisectors between distinct anchors
    (`VoronoiBoundary_subset_bisectors`, §13). Each bisector is
    an affine line in `ℂ` and has measure zero (§18.1). A finite
    union of measure-zero sets is measure-zero. -/
theorem voronoiBoundary_volume_zero {p : ℕ}
    (γ : Fin p → ℂ) (h_inj : Function.Injective γ) :
    volume (VoronoiBoundary γ) = 0 := by
  refine measure_mono_null (VoronoiBoundary_subset_bisectors γ) ?_
  apply measure_iUnion_null
  intro s
  apply measure_iUnion_null
  intro t
  apply measure_iUnion_null
  intro hst
  exact perpBisector_volume_zero (γ s) (γ t) (fun h => hst (h_inj h))

end VoronoiBoundaryMeasure

/-! ============================================================
  §18.3  Almost-everywhere constructive McMullen
============================================================ -/

section ConstructiveMcMullenAE

/-- **★ Almost-everywhere constructive McMullen.**
    The `constructive_mcmullen` conclusion `(B)` (uniqueness of the
    nearest-anchor selector) holds *almost everywhere*: the set of
    query points where it can fail is a Lebesgue-null subset of `ℂ`.

    This upgrades `constructive_mcmullen` from "off a specific
    boundary set" to the full "almost everywhere" statement, which
    is the complex-analytic content of the result. -/
theorem constructive_mcmullen_ae
    {p : ℕ} (hp : 0 < p) (x A : ℂ) (hx : x ≠ 0) (hA : A ≠ 0)
    (γ : Fin p → ℂ) (h_inj : Function.Injective γ)
    (hSp : ∀ s, Sp_C p (γ s) ≠ 0)
    (hfp : ∀ s, (γ s) ^ p = A / x) :
    (∀ k : Fin p, γ k ∈ interior (basin γ k)) ∧
    (⋃ k : Fin p, basin γ k) = Set.univ ∧
    volume (VoronoiBoundary γ) = 0 ∧
    (∀ z₀ ∉ VoronoiBoundary γ,
        ∃! s : Fin p, ∀ t : Fin p, ‖γ s - z₀‖ ≤ ‖γ t - z₀‖) ∧
    (∀ s : Fin p, steffensen_step_C (x / A) p (γ s) = γ s) := by
  obtain ⟨hA_int, hcover, hunique, hfix⟩ :=
    constructive_mcmullen hp x A hx hA γ h_inj hSp hfp
  exact ⟨hA_int, hcover, voronoiBoundary_volume_zero γ h_inj, hunique, hfix⟩

end ConstructiveMcMullenAE

end Pandrosion
