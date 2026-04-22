/-
  Universitas Pandrosion — §51. **Slow-set geometric limit on the
  explicit Julia sphere (`p = 2` complex).**

  Combines §43 (`slowSet_iInter_eq_juliaSetP`: the decreasing
  intersection of slow sets equals the polynomial Julia locus)
  with §48 (`juliaSet_eq_sphere`: the Julia circle *is* an
  explicit Apollonius metric sphere) to produce the **explicit
  geometric limit**

      `⋂ n, slowSet α (1/(n+1)) ⊆ Metric.sphere (juliaCenter α) (juliaRadius α)`

  and the companion reverse containment on the Julia circle's
  complement in `{−α}`.

  This provides the structural foundation for a quantitative
  `volume (slowSet α ε ∩ B(0, R)) ≤ C(α, R) · ε` bound (a tubular
  neighbourhood of the metric sphere). The quantitative bound
  itself requires integration of the radial density, which is a
  dedicated piece of measure theory beyond this module.

  Contents.

    §51.1  `slowSet_iInter_subset_sphere` — decreasing intersection
           contained in the explicit Julia sphere.

    §51.2  `juliaSphere_subset_every_slowSet` — the Julia sphere
           is contained in every `slowSet α ε` (ε ≥ 0).

    §51.3  `sphere_volume_zero_via_slowSet` — alternate derivation
           of §48.volume-zero via the sphere identification (sanity
           check corollary).
-/

import Pandrosion.Core.ComplexJuliaGeometryP2
import Pandrosion.Core.ComplexJuliaSlowSetDecay

namespace Pandrosion

open MeasureTheory Complex

/-! ============================================================
  §51.1  Decreasing intersection ⊆ Julia sphere
============================================================ -/

section SlowSetIInterSphereC

/-- **`juliaSetP α ⊆ Metric.sphere (juliaCenter α) (juliaRadius α)`**.

    `juliaSetP α = {z : ‖μ²‖²·‖z−α‖² = ‖z+α‖²}` is the §43
    polynomial form of the Julia circle. By §48, this coincides
    with the metric sphere for `α` Apollonius-regular. -/
theorem juliaSetP_subset_juliaSphere
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    juliaSetP α ⊆ Metric.sphere (juliaCenter α) (juliaRadius α) := by
  -- juliaSetP α ⊆ {z : ‖v(z)‖ = 1} (§43's `juliaSetP_volume_zero` proof structure).
  intro z hz
  -- `hz : ‖μ²‖²·‖z−α‖² = ‖z+α‖²`. Split on z + α = 0.
  by_cases h_zα : z + α = 0
  · -- z = -α. Non-sphere by §48.
    exfalso
    have h_val : z = -α := by linear_combination h_zα
    have hα2 : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖-α - α‖ ^ 2 = ‖(-α) + α‖ ^ 2 := by
      rw [← h_val]; exact hz
    have h_neg_α_α : (-α + α : ℂ) = 0 := by ring
    rw [h_neg_α_α, norm_zero, zero_pow (by norm_num : 2 ≠ 0)] at hα2
    have h_mul_zero : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖-α - α‖ ^ 2 = 0 := hα2
    rcases mul_eq_zero.mp h_mul_zero with h_μ | h_norm_neg
    · -- ‖μ²‖² = 0 ⟹ μ = 0 ⟹ α = 1.
      have h_μ_zero : ‖((1 - α) / (1 + α)) ^ 2‖ = 0 :=
        pow_eq_zero_iff (n := 2) (by norm_num : (2:ℕ) ≠ 0) |>.mp h_μ
      have h_div_zero : ((1 - α) / (1 + α)) ^ 2 = 0 := norm_eq_zero.mp h_μ_zero
      have h_div_plain : (1 - α) / (1 + α) = 0 :=
        pow_eq_zero_iff (n := 2) (by norm_num : (2:ℕ) ≠ 0) |>.mp h_div_zero
      rcases div_eq_zero_iff.mp h_div_plain with h1 | h2
      · exact hα_ne_one h1
      · exact hα_ne_neg_one h2
    · -- ‖-α-α‖² = 0 ⟹ -2α = 0 ⟹ α = 0.
      have h_norm_zero : ‖(-α - α : ℂ)‖ = 0 :=
        pow_eq_zero_iff (n := 2) (by norm_num : (2:ℕ) ≠ 0) |>.mp h_norm_neg
      have h_eq : (-α - α : ℂ) = 0 := norm_eq_zero.mp h_norm_zero
      have h_2α : (2 * α : ℂ) = 0 := by linear_combination -h_eq
      rcases mul_eq_zero.mp h_2α with h1 | h2
      · exact absurd h1 two_ne_zero
      · exact hα_ne_zero h2
  · -- z + α ≠ 0. Apply §39.2 `v_p2_C_norm_sq_eq` + §48 sphere equality.
    have hzα_pos : 0 < ‖z + α‖ ^ 2 := by
      have : 0 < ‖z + α‖ := norm_pos_iff.mpr h_zα
      positivity
    have h_ratio : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 / ‖z + α‖ ^ 2 = 1 := by
      rw [div_eq_one_iff_eq (ne_of_gt hzα_pos)]; exact hz
    have h_v_sq : ‖v_p2_C α z‖ ^ 2 = 1 := by
      rw [v_p2_C_norm_sq_eq α z]; exact h_ratio
    have h_v_nn : 0 ≤ ‖v_p2_C α z‖ := norm_nonneg _
    have h_v_eq_one : ‖v_p2_C α z‖ = 1 := by
      rw [show (1 : ℝ) = 1 ^ 2 from by ring] at h_v_sq
      exact (sq_eq_sq h_v_nn zero_le_one).mp h_v_sq
    -- §48: {‖v‖ = 1} = sphere.
    have h_eq := juliaSet_eq_sphere α hα_ne_zero hα_ne_neg_one hα_ne_one hμ
    have h_mem : z ∈ {z : ℂ | ‖v_p2_C α z‖ = 1} := h_v_eq_one
    rw [h_eq] at h_mem
    exact h_mem

/-- **★★★★ Decreasing intersection of `slowSet α (1/(n+1))` is contained
    in the explicit Julia sphere.**

    Combines §43 (`slowSet_iInter_eq_juliaSetP`) with §51.1 to
    produce a direct geometric limit: as `ε → 0`, `slowSet α ε`
    collapses onto the Apollonius metric sphere
    `Metric.sphere (juliaCenter α) (juliaRadius α)`. -/
theorem slowSet_iInter_subset_sphere
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    (⋂ n : ℕ, slowSet α (1 / ((n : ℝ) + 1))) ⊆
      Metric.sphere (juliaCenter α) (juliaRadius α) := by
  rw [slowSet_iInter_eq_juliaSetP α]
  exact juliaSetP_subset_juliaSphere α hα_ne_zero hα_ne_neg_one hα_ne_one hμ

end SlowSetIInterSphereC

/-! ============================================================
  §51.2  Julia sphere ⊆ every `slowSet α ε` (ε ≥ 0)
============================================================ -/

section SphereInsideSlowSetC

/-- **Every point of the Julia sphere is in every slow set** (ε ≥ 0,
    ε ≤ 1). Direct from §48 (sphere ⊆ Julia) + §43
    `juliaSetP_subset_slowSet`. -/
theorem juliaSphere_subset_slowSet
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1)
    (ε : ℝ) (hε_nn : 0 ≤ ε) (hε_le : ε ≤ 1) :
    Metric.sphere (juliaCenter α) (juliaRadius α) ⊆ slowSet α ε := by
  -- sphere = {‖v‖ = 1} (§48), which ⊆ juliaSetP (via §51.1 backward),
  -- which ⊆ slowSet α ε (via §43 `juliaSetP_subset_slowSet`).
  intro z hz
  -- Use §48 equality forward: z ∈ sphere ⟹ z ∈ {‖v‖ = 1}.
  have h_eq := juliaSet_eq_sphere α hα_ne_zero hα_ne_neg_one hα_ne_one hμ
  have h_v_set : z ∈ {z : ℂ | ‖v_p2_C α z‖ = 1} := by
    rw [h_eq]; exact hz
  have h_v_eq_one : ‖v_p2_C α z‖ = 1 := h_v_set
  -- z + α ≠ 0 on the sphere: -α ∉ sphere by §48's internal lemma
  -- (provable, but we can bypass via v_eq_one excluding z = -α since
  -- v_p2_C α (-α) = 0 ≠ 1).
  have h_zα_ne : z + α ≠ 0 := by
    intro h_eq_zα
    have hz_val : z = -α := by linear_combination h_eq_zα
    rw [hz_val] at h_v_eq_one
    unfold v_p2_C at h_v_eq_one
    have h_den : (-α + α : ℂ) = 0 := by ring
    rw [h_den, div_zero, mul_zero, norm_zero] at h_v_eq_one
    norm_num at h_v_eq_one
  -- Now z ∈ juliaSetP (polynomial form).
  have h_juliaP : z ∈ juliaSetP α := by
    show ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2
    have h_sq : ‖v_p2_C α z‖ ^ 2 = 1 := by rw [h_v_eq_one]; ring
    have h_sq_eq := v_p2_C_norm_sq_eq α z
    have hzα_sq_pos : 0 < ‖z + α‖ ^ 2 := by
      have : 0 < ‖z + α‖ := norm_pos_iff.mpr h_zα_ne
      positivity
    rw [h_sq_eq] at h_sq
    rw [div_eq_one_iff_eq (ne_of_gt hzα_sq_pos)] at h_sq
    exact h_sq
  -- §43: juliaSetP α ⊆ slowSet α ε.
  exact juliaSetP_subset_slowSet α hε_nn hε_le h_juliaP

end SphereInsideSlowSetC

/-! ============================================================
  §51.3  Alternate sphere-based Lebesgue-zero of Julia
============================================================ -/

section SphereVolumeZeroC

/-- **Alternate derivation of `juliaSetP_volume_zero` via §48
    sphere identification.** Direct chain
    `juliaSetP α ⊆ Metric.sphere (…) (…) ⟹ volume = 0` by
    `Measure.addHaar_sphere`. Provides an independent proof of
    `juliaSetP_volume_zero` (originally through `v_p2_C`). -/
theorem juliaSetP_volume_zero_via_sphere
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    volume (juliaSetP α) = 0 := by
  have h_subset : juliaSetP α ⊆ Metric.sphere (juliaCenter α) (juliaRadius α) :=
    juliaSetP_subset_juliaSphere α hα_ne_zero hα_ne_neg_one hα_ne_one hμ
  have h_sphere_zero : volume (Metric.sphere (juliaCenter α) (juliaRadius α)) = 0 :=
    MeasureTheory.Measure.addHaar_sphere volume (juliaCenter α) (juliaRadius α)
  exact measure_mono_null h_subset h_sphere_zero

end SphereVolumeZeroC

end Pandrosion
