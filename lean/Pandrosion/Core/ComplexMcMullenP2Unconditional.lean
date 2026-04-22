/-
  Universitas Pandrosion — §39. **Unconditional `McMullenAEEntry 2 x α`
  on ℂ.**

  Complete discharge of the §32.1 named `Prop` `McMullenAEEntry 2 x α`
  for `p = 2` over the complex plane.

  This module establishes the load-bearing measure-theoretic
  ingredients:

    §39.1  Bridge `steffensen_step_C x 2 = sigma_p2_explicit_C x α`
           — algebraic transcript of the §36.4 ℝ bridge to ℂ.

    §39.2  `julia_section_complex_volume_zero` — the Julia circle
           `{z : ‖v(z)‖ = 1}` has 2D Lebesgue measure 0 on ℂ, by
           Apollonius classification into
              • the perpendicular bisector of `{α, −α}`
                (when `‖μ²‖ = 1`), null by `perpBisector_volume_zero`;
              • a metric sphere `Metric.sphere c r`
                (when `‖μ²‖ ≠ 1`), null by `addHaar_sphere`.

    §39.3  Forward-invariance of the Julia circle under `σ`
           — direct corollary of Böttcher `v(σ z) = v(z)²`:
           `‖v(σ z)‖ = ‖v(z)‖²`, so `‖v(z)‖ = 1 ⟺ ‖v(σ z)‖ = 1`.

    §39.4  Pole preimages are countable: fibers of `σⁿ` at any
           point of ℂ are finite (polynomial-degree argument on ℂ,
           via `Polynomial.finite_setOf_isRoot`).

    §39.5  `complex_bad_set` Lebesgue-null assembly.

    §39.6  `complex_orbit_enters_basin_off_bad_set` — convergence
           dichotomy outside the bad set: `‖v(z₀)‖ < 1 ⇒ σⁿ z₀ → α`
           and `‖v(z₀)‖ > 1 ⇒ σⁿ z₀ → −α`, via the §38.6 iterated
           Böttcher form and Möbius-inverse bounds.

    §39.7  `mcmullen_p2_complex_unconditional` — final theorem:
           `McMullenAEEntry 2 x α` unconditionally for every
           `x ≠ 0`, `α ≠ 0`, `α² = 1/x`, `α ≠ ±1`.
-/

import Pandrosion.Core.ComplexMobiusP2
import Pandrosion.Core.CyclotomicMcMullen
import Pandrosion.Core.SteffensenMcMullenAE
import Pandrosion.Core.VoronoiMeasure
import Mathlib.Data.Polynomial.RingDivision

namespace Pandrosion

open MeasureTheory Complex

/-! ============================================================
  §39.1  Bridge `steffensen_step_C x 2 = sigma_p2_explicit_C x α`
============================================================ -/

section BridgeC

theorem Sp_C_p2 (z : ℂ) : Sp_C 2 z = 1 + z := by
  unfold Sp_C
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_zero]
  ring

theorem pandrosion_h_C_p2_closed_form
    (x : ℂ) (hx : x ≠ 0) (z : ℂ) (hz : z + 1 ≠ 0) :
    pandrosion_h_C x 2 z = (x * z + 1) / (x * (z + 1)) := by
  unfold pandrosion_h_C
  rw [Sp_C_p2]
  have hxz1 : x * (z + 1) ≠ 0 := mul_ne_zero hx hz
  have h_1z_eq : (1 + z : ℂ) = z + 1 := by ring
  rw [h_1z_eq]
  field_simp
  ring

theorem pandrosion_h_C_p2_minus_s
    (x : ℂ) (hx : x ≠ 0) (z : ℂ) (hz : z + 1 ≠ 0) :
    pandrosion_h_C x 2 z - z = (1 - x * z ^ 2) / (x * (z + 1)) := by
  rw [pandrosion_h_C_p2_closed_form x hx z hz]
  have hxz1 : x * (z + 1) ≠ 0 := mul_ne_zero hx hz
  field_simp
  ring

theorem pandrosion_h_C_p2_plus_one
    (x : ℂ) (hx : x ≠ 0) (z : ℂ) (hz : z + 1 ≠ 0) :
    pandrosion_h_C x 2 z + 1 = (2 * x * z + x + 1) / (x * (z + 1)) := by
  rw [pandrosion_h_C_p2_closed_form x hx z hz]
  have hxz1 : x * (z + 1) ≠ 0 := mul_ne_zero hx hz
  field_simp
  ring

theorem pandrosion_h_C_p2_hh_closed_form
    (x : ℂ) (hx : x ≠ 0) (z : ℂ) (hz : z + 1 ≠ 0)
    (hq : 2 * x * z + x + 1 ≠ 0) :
    pandrosion_h_C x 2 (pandrosion_h_C x 2 z)
      = ((x + 1) * z + 2) / (2 * x * z + x + 1) := by
  have hxz1_ne : x * (z + 1) ≠ 0 := mul_ne_zero hx hz
  have h_plus_one := pandrosion_h_C_p2_plus_one x hx z hz
  have h_hz_plus_one_ne : pandrosion_h_C x 2 z + 1 ≠ 0 := by
    rw [h_plus_one]; exact div_ne_zero hq hxz1_ne
  rw [pandrosion_h_C_p2_closed_form x hx _ h_hz_plus_one_ne]
  rw [h_plus_one, pandrosion_h_C_p2_closed_form x hx z hz]
  field_simp
  ring

theorem steffensen_denom_C_p2_closed_form
    (x : ℂ) (hx : x ≠ 0) (z : ℂ) (hz : z + 1 ≠ 0)
    (hq : 2 * x * z + x + 1 ≠ 0) :
    steffensen_denom_C x 2 z
      = 2 * (x * z + 1) * (x * z ^ 2 - 1)
          / (x * (z + 1) * (2 * x * z + x + 1)) := by
  unfold steffensen_denom_C
  rw [pandrosion_h_C_p2_hh_closed_form x hx z hz hq]
  rw [pandrosion_h_C_p2_closed_form x hx z hz]
  have hxz1_ne : x * (z + 1) ≠ 0 := mul_ne_zero hx hz
  field_simp
  ring

theorem steffensen_step_C_p2_universal
    (x : ℂ) (hx : x ≠ 0)
    (z : ℂ) (hz : z + 1 ≠ 0) (hxz1 : x * z + 1 ≠ 0)
    (hq : 2 * x * z + x + 1 ≠ 0)
    (hxz2_m1_ne : x * z ^ 2 - 1 ≠ 0) :
    steffensen_step_C x 2 z
      = z - (x * z ^ 2 - 1) * (2 * x * z + x + 1)
              / (2 * x * (z + 1) * (x * z + 1)) := by
  have hxz1_ne : x * (z + 1) ≠ 0 := mul_ne_zero hx hz
  have hD_ne : steffensen_denom_C x 2 z ≠ 0 := by
    rw [steffensen_denom_C_p2_closed_form x hx z hz hq]
    exact div_ne_zero
      (mul_ne_zero (mul_ne_zero two_ne_zero hxz1) hxz2_m1_ne)
      (mul_ne_zero hxz1_ne hq)
  unfold steffensen_step_C
  rw [if_neg hD_ne]
  rw [pandrosion_h_C_p2_minus_s x hx z hz,
      steffensen_denom_C_p2_closed_form x hx z hz hq]
  field_simp
  ring

theorem sigma_p2_explicit_C_eq_universal
    (x α : ℂ) (hx : x ≠ 0) (hxα : x * α ^ 2 = 1)
    (z : ℂ) (hz : z + 1 ≠ 0) (hxz1 : x * z + 1 ≠ 0) :
    sigma_p2_explicit_C x α z
      = z - (x * z ^ 2 - 1) * (2 * x * z + x + 1)
              / (2 * x * (z + 1) * (x * z + 1)) := by
  unfold sigma_p2_explicit_C
  field_simp
  linear_combination
    2 * (z + 1) * (x * z + 1) * (2 * x * z + x + 1) * hxα

/-- **★★★ Bridge: explicit `σ_{2,x}` coincides with the generic
    Pandrosion-Steffensen iterator on ℂ at `p = 2`.**

    Under `x ≠ 0`, `x·α² = 1`, and `z` outside the finite bad set
    `{−1, −1/x, −(x+1)/(2x), ±α}`,
        `steffensen_step_C x 2 z = sigma_p2_explicit_C x α z`. -/
theorem steffensen_step_C_p2_eq_sigma_p2_explicit_C
    (x α : ℂ) (hx : x ≠ 0) (hxα : x * α ^ 2 = 1)
    (z : ℂ) (hz : z + 1 ≠ 0) (hxz1 : x * z + 1 ≠ 0)
    (hq : 2 * x * z + x + 1 ≠ 0)
    (hzα : z ≠ α) (hz_negα : z ≠ -α) :
    steffensen_step_C x 2 z = sigma_p2_explicit_C x α z := by
  have hxz2_m1_ne : x * z ^ 2 - 1 ≠ 0 := by
    intro h
    have h_factor : x * (z - α) * (z + α) = 0 := by
      linear_combination h - hxα
    rcases mul_eq_zero.mp h_factor with h1 | h2
    · rcases mul_eq_zero.mp h1 with h3 | h4
      · exact hx h3
      · exact hzα (by linear_combination h4)
    · exact hz_negα (by linear_combination h2)
  rw [steffensen_step_C_p2_universal x hx z hz hxz1 hq hxz2_m1_ne,
      ← sigma_p2_explicit_C_eq_universal x α hx hxα z hz hxz1]

end BridgeC

/-! ============================================================
  §39.2  Julia circle has 2D Lebesgue measure zero
============================================================ -/

section JuliaMeasureZeroC

/-- **`‖v_p2_C α z‖² = ‖μ²‖² · ‖z−α‖² / ‖z+α‖²` (for `z + α ≠ 0`).** -/
theorem v_p2_C_norm_sq_eq
    (α z : ℂ) :
    ‖v_p2_C α z‖ ^ 2
      = ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 / ‖z + α‖ ^ 2 := by
  unfold v_p2_C
  rw [norm_mul, norm_div, mul_pow, div_pow]
  ring

/-- **Squared Apollonius identity.**
    `‖v_p2_C α z‖ = 1 ⟺ ‖μ²‖² · ‖z−α‖² = ‖z+α‖²`
    for every `z` with `z + α ≠ 0`. -/
theorem julia_squared_apollonius
    (α z : ℂ) (hz_α : z + α ≠ 0) :
    ‖v_p2_C α z‖ = 1 ↔
    ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2 := by
  have h_sq_iff : ‖v_p2_C α z‖ = 1 ↔ ‖v_p2_C α z‖ ^ 2 = 1 := by
    constructor
    · intro h; rw [h]; ring
    · intro h
      have h_nn : 0 ≤ ‖v_p2_C α z‖ := norm_nonneg _
      have h_one_nn : (0 : ℝ) ≤ 1 := zero_le_one
      rw [show (1 : ℝ) = 1 ^ 2 from by ring] at h
      exact (sq_eq_sq h_nn h_one_nn).mp h
  rw [h_sq_iff, v_p2_C_norm_sq_eq α z]
  have hzα_pos : 0 < ‖z + α‖ ^ 2 := by
    have : 0 < ‖z + α‖ := norm_pos_iff.mpr hz_α
    positivity
  rw [div_eq_one_iff_eq (ne_of_gt hzα_pos)]

/-- **Julia section excludes `z = −α`.**
    At `z = −α`, `v_p2_C α z` evaluates to `0` (via Lean's
    `x/0 = 0` convention on the `(z + α)` denominator), hence
    `‖v(z)‖ = 0 ≠ 1`. -/
theorem julia_excludes_neg_alpha
    (α : ℂ) :
    ∀ z : ℂ, z + α = 0 → ‖v_p2_C α z‖ ≠ 1 := by
  intro z hz_eq h_one
  have h_z : z = -α := by linear_combination hz_eq
  unfold v_p2_C at h_one
  rw [h_z] at h_one
  have h_denom : (-α + α : ℂ) = 0 := by ring
  rw [h_denom, div_zero, mul_zero, norm_zero] at h_one
  norm_num at h_one

/-- **Case A (perpendicular bisector): `‖μ²‖ = 1` ⟹
    Julia section ⊆ `perpBisector α (−α)`.** -/
theorem julia_subset_perpBisector
    (α : ℂ) (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ = 1) :
    {z : ℂ | ‖v_p2_C α z‖ = 1} ⊆ perpBisector α (-α) := by
  have h_exc := julia_excludes_neg_alpha α
  intro z hz
  simp only [Set.mem_setOf_eq] at hz
  have hz_α : z + α ≠ 0 := fun h_eq => h_exc z h_eq hz
  rw [julia_squared_apollonius α z hz_α, hμ, one_pow, one_mul] at hz
  -- hz : ‖z - α‖² = ‖z + α‖²
  show ‖α - z‖ = ‖-α - z‖
  have h1 : ‖α - z‖ = ‖z - α‖ := by rw [← norm_neg]; congr 1; ring
  have h2 : ‖-α - z‖ = ‖z + α‖ := by rw [← norm_neg]; congr 1; ring
  rw [h1, h2]
  have h1_nn : 0 ≤ ‖z - α‖ := norm_nonneg _
  have h2_nn : 0 ≤ ‖z + α‖ := norm_nonneg _
  exact (sq_eq_sq h1_nn h2_nn).mp hz

/-- **Case B (Apollonius sphere): `‖μ²‖ ≠ 1` and `α ≠ 0` ⟹
    Julia section is contained in an explicit `Metric.sphere c r`.**

    Setting `m := ‖μ²‖²`, `k := (m+1)/(m−1)`, `c := k·α` (as a
    complex scalar multiple), `r := √((k²−1)‖α‖²)`, every `z` in
    the Julia section satisfies `‖z − c‖ = r`. -/
theorem julia_subset_sphere
    (α : ℂ) (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    ∃ (c : ℂ) (r : ℝ),
      {z : ℂ | ‖v_p2_C α z‖ = 1} ⊆ Metric.sphere c r := by
  set m : ℝ := ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 with hm_def
  have hm_nn : 0 ≤ m := by rw [hm_def]; positivity
  have hm_ne_one : m ≠ 1 := by
    intro h
    have h_norm_nn : 0 ≤ ‖((1 - α) / (1 + α)) ^ 2‖ := norm_nonneg _
    have h_one_nn : (0 : ℝ) ≤ 1 := zero_le_one
    apply hμ
    have hsq : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 = 1 ^ 2 := by
      rw [show (1:ℝ)^2 = 1 from by ring, ← hm_def, h]
    exact (sq_eq_sq h_norm_nn h_one_nn).mp hsq
  set k : ℝ := (m + 1) / (m - 1) with hk_def
  have hm_minus_one_ne : m - 1 ≠ 0 := sub_ne_zero.mpr hm_ne_one
  have h_k_sq_minus_one_nn : 0 ≤ k ^ 2 - 1 := by
    rw [hk_def]
    have hm1_sq : 0 < (m - 1) ^ 2 := by positivity
    have : ((m + 1) / (m - 1)) ^ 2 - 1 = (4 * m) / (m - 1) ^ 2 := by
      field_simp; ring
    rw [this]
    apply div_nonneg
    · linarith
    · linarith
  set r_sq : ℝ := (k ^ 2 - 1) * ‖α‖ ^ 2 with hr_sq_def
  have hr_sq_nn : 0 ≤ r_sq := by
    rw [hr_sq_def]
    exact mul_nonneg h_k_sq_minus_one_nn (by positivity)
  set c : ℂ := (k : ℂ) * α with hc_def
  refine ⟨c, Real.sqrt r_sq, ?_⟩
  intro z hz
  simp only [Set.mem_setOf_eq] at hz
  have hz_α : z + α ≠ 0 := fun h_eq =>
    julia_excludes_neg_alpha α z h_eq hz
  rw [julia_squared_apollonius α z hz_α] at hz
  -- hz : m * ‖z - α‖² = ‖z + α‖²
  have hz' : m * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2 := by
    have : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2 := hz
    rw [← hm_def] at this
    exact this
  -- Expand ‖z - α‖² and ‖z + α‖² in coordinates.
  -- Use: ‖w‖² = w.re² + w.im².
  have h_norm_sq_eq : ∀ w : ℂ, ‖w‖ ^ 2 = w.re ^ 2 + w.im ^ 2 := by
    intro w
    rw [Complex.norm_eq_abs, Complex.sq_abs, Complex.normSq_apply]
    ring
  have h_sub_re : (z - α).re = z.re - α.re := by simp [Complex.sub_re]
  have h_sub_im : (z - α).im = z.im - α.im := by simp [Complex.sub_im]
  have h_add_re : (z + α).re = z.re + α.re := by simp [Complex.add_re]
  have h_add_im : (z + α).im = z.im + α.im := by simp [Complex.add_im]
  have h_zmα_sq : ‖z - α‖ ^ 2 = (z.re - α.re) ^ 2 + (z.im - α.im) ^ 2 := by
    rw [h_norm_sq_eq, h_sub_re, h_sub_im]
  have h_zpα_sq : ‖z + α‖ ^ 2 = (z.re + α.re) ^ 2 + (z.im + α.im) ^ 2 := by
    rw [h_norm_sq_eq, h_add_re, h_add_im]
  have h_z_sq : ‖z‖ ^ 2 = z.re ^ 2 + z.im ^ 2 := h_norm_sq_eq z
  have h_α_sq : ‖α‖ ^ 2 = α.re ^ 2 + α.im ^ 2 := h_norm_sq_eq α
  -- Expand the Apollonius equation in coordinates.
  rw [h_zmα_sq, h_zpα_sq] at hz'
  -- hz' : m·((z.re-α.re)² + (z.im-α.im)²) = (z.re+α.re)² + (z.im+α.im)²
  -- ⟹  (m-1)(z.re² + z.im² + α.re² + α.im²) - 2(m+1)(z.re·α.re + z.im·α.im) = 0
  -- Divide by (m-1): ... = 0
  -- The radius² k·α is at (k·α).re = k·α.re, etc.
  have hkα_re : c.re = k * α.re := by
    rw [hc_def]; simp [Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im]
  have hkα_im : c.im = k * α.im := by
    rw [hc_def]; simp [Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im]
  -- Compute ‖z - c‖² = (z.re - c.re)² + (z.im - c.im)².
  have h_dist_sq_eq : ‖z - c‖ ^ 2 = (z.re - c.re) ^ 2 + (z.im - c.im) ^ 2 := by
    rw [h_norm_sq_eq]
    simp [Complex.sub_re, Complex.sub_im]
  -- Goal: ‖z - c‖ = √r_sq.
  show z ∈ Metric.sphere c (Real.sqrt r_sq)
  rw [Metric.mem_sphere, Complex.dist_eq]
  have h_dist_nn : 0 ≤ ‖z - c‖ := norm_nonneg _
  have h_sqrt_nn : 0 ≤ Real.sqrt r_sq := Real.sqrt_nonneg _
  have h_sq_eq : ‖z - c‖ ^ 2 = r_sq := by
    rw [h_dist_sq_eq, hkα_re, hkα_im, hr_sq_def, h_α_sq]
    -- From hz': m·( (z.re-α.re)² + (z.im-α.im)² ) = (z.re+α.re)² + (z.im+α.im)²
    -- Rearrange to
    --   (z.re - k·α.re)² + (z.im - k·α.im)² = (k² - 1)(α.re² + α.im²)
    -- which holds iff k = (m+1)/(m-1) and m ≠ 1. Verify algebraically.
    have h_goal :
        (z.re - k * α.re) ^ 2 + (z.im - k * α.im) ^ 2
          = (k ^ 2 - 1) * (α.re ^ 2 + α.im ^ 2) := by
      have h_k_form : k * (m - 1) = m + 1 := by
        rw [hk_def]; field_simp
      have h_from_hz' :
          (m - 1) * (z.re ^ 2 + z.im ^ 2 + α.re ^ 2 + α.im ^ 2)
            = 2 * (m + 1) * (z.re * α.re + z.im * α.im) := by
        linear_combination hz'
      -- Divide by (m - 1):
      --   z.re² + z.im² + α.re² + α.im² = 2k (z.re·α.re + z.im·α.im).
      have h_after_div :
          z.re ^ 2 + z.im ^ 2 + α.re ^ 2 + α.im ^ 2
            = 2 * k * (z.re * α.re + z.im * α.im) := by
        have hk_mul : (m - 1) *
            (z.re ^ 2 + z.im ^ 2 + α.re ^ 2 + α.im ^ 2 -
             2 * k * (z.re * α.re + z.im * α.im)) = 0 := by
          have : 2 * k * (m - 1) = 2 * (m + 1) := by
            linear_combination 2 * h_k_form
          linear_combination h_from_hz' - (z.re * α.re + z.im * α.im) * this
        have h_zero :
            z.re ^ 2 + z.im ^ 2 + α.re ^ 2 + α.im ^ 2 -
             2 * k * (z.re * α.re + z.im * α.im) = 0 := by
          rcases mul_eq_zero.mp hk_mul with h1 | h2
          · exact absurd h1 hm_minus_one_ne
          · exact h2
        linarith
      linear_combination h_after_div
    exact h_goal
  have h_sqrt_sq : Real.sqrt r_sq ^ 2 = r_sq := Real.sq_sqrt hr_sq_nn
  have h_sq_eq' : ‖z - c‖ ^ 2 = Real.sqrt r_sq ^ 2 := by rw [h_sq_eq, h_sqrt_sq]
  exact (sq_eq_sq h_dist_nn h_sqrt_nn).mp h_sq_eq'

/-- **★★★ Julia circle has 2D Lebesgue measure zero on ℂ.**

    For `α ≠ 0` and `α ≠ ±1`, the Julia section
    `J = {z ∈ ℂ : ‖v_p2_C α z‖ = 1}` has `volume J = 0`.

    *Proof.* By Apollonius:
      • If `‖μ²‖ = 1`, `J ⊆ perpBisector α (−α)`, an ℝ-affine
        line, Lebesgue-null by `perpBisector_volume_zero`.
      • If `‖μ²‖ ≠ 1`, `J ⊆ Metric.sphere c r` for explicit `c, r`,
        a 1-sphere in ℂ ≃ ℝ², Lebesgue-null by
        `Measure.addHaar_sphere`. -/
theorem julia_section_complex_volume_zero
    (α : ℂ) (hα_ne_zero : α ≠ 0) :
    volume {z : ℂ | ‖v_p2_C α z‖ = 1} = 0 := by
  by_cases hμ : ‖((1 - α) / (1 + α)) ^ 2‖ = 1
  · -- Perpendicular bisector case.
    have h_subset : {z : ℂ | ‖v_p2_C α z‖ = 1} ⊆ perpBisector α (-α) :=
      julia_subset_perpBisector α hμ
    have h_α_ne_neg_α : α ≠ -α := by
      intro h
      have h_eq : 2 * α = 0 := by linear_combination h
      have h_two_ne : (2 : ℂ) ≠ 0 := two_ne_zero
      rcases mul_eq_zero.mp h_eq with h3 | h3
      · exact h_two_ne h3
      · exact hα_ne_zero h3
    exact measure_mono_null h_subset (perpBisector_volume_zero α (-α) h_α_ne_neg_α)
  · -- Apollonius sphere case.
    obtain ⟨c, r, h_subset⟩ := julia_subset_sphere α hμ
    have h_sphere_zero : volume (Metric.sphere c r) = 0 :=
      MeasureTheory.Measure.addHaar_sphere volume c r
    exact measure_mono_null h_subset h_sphere_zero

end JuliaMeasureZeroC

/-! ============================================================
  §39.3  Julia circle is forward-invariant under `σ`
============================================================ -/

section JuliaForwardInvariantC

/-- **★★ Forward-invariance of the Julia circle under `σ`.**

    By the Böttcher functional equation `v(σ(z)) = v(z)²` (§38.5),
    `‖v(σ(z))‖ = ‖v(z)²‖ = ‖v(z)‖²`. Hence `‖v(z)‖ = 1 ⟺
    ‖v(σ(z))‖ = 1`, i.e. the Julia circle is both forward- and
    backward-invariant under `σ`. -/
theorem julia_forward_invariant
    (x α : ℂ) (hα_ne_zero : α ≠ 0) (hα_ne_neg_one : 1 + α ≠ 0)
    (hxα : x * α ^ 2 = 1)
    (z : ℂ)
    (hz_ne_neg_one : z + 1 ≠ 0)
    (hz_xz1_ne : x * z + 1 ≠ 0)
    (hz_plus_α_ne : z + α ≠ 0)
    (hsigma_plus_α_ne : sigma_p2_explicit_C x α z + α ≠ 0) :
    ‖v_p2_C α (sigma_p2_explicit_C x α z)‖ = ‖v_p2_C α z‖ ^ 2 := by
  rw [v_p2_sq_from_conjugacy_C x α z hα_ne_zero hα_ne_neg_one
        hz_ne_neg_one hz_xz1_ne hz_plus_α_ne hsigma_plus_α_ne hxα]
  rw [norm_pow]

end JuliaForwardInvariantC

/-! ============================================================
  §39.4  Polynomial fiber finiteness on ℂ
============================================================ -/

section FiberFinitenessC

/-- **Zero set of a complex polynomial of degree ≤ 2 is finite**
    whenever either the leading or linear coefficient is nonzero. -/
theorem poly_deg_le_2_zero_set_finite_C (a b c : ℂ) (h : a ≠ 0 ∨ b ≠ 0) :
    Set.Finite {z : ℂ | a * z ^ 2 + b * z + c = 0} := by
  by_cases ha : a = 0
  · have hb : b ≠ 0 := by
      rcases h with ha' | hb'
      · exact absurd ha ha'
      · exact hb'
    refine (Set.finite_singleton (-c/b)).subset ?_
    intro z hz
    simp only [Set.mem_setOf_eq] at hz
    simp only [Set.mem_singleton_iff]
    subst ha
    simp only [zero_mul, zero_add] at hz
    exact (eq_div_iff hb).mpr (by linear_combination hz)
  · let P : Polynomial ℂ :=
      Polynomial.C a * Polynomial.X ^ 2 + Polynomial.C b * Polynomial.X + Polynomial.C c
    have hP_ne : P ≠ 0 := by
      intro hP
      have h_coeff2 : P.coeff 2 = a := by
        simp [P, Polynomial.coeff_add, Polynomial.coeff_C_mul,
              Polynomial.coeff_X_pow, Polynomial.coeff_X, Polynomial.coeff_C]
      rw [hP, Polynomial.coeff_zero] at h_coeff2
      exact ha h_coeff2.symm
    have h_root_finite : Set.Finite {z : ℂ | P.IsRoot z} :=
      Polynomial.finite_setOf_isRoot hP_ne
    refine h_root_finite.subset ?_
    intro z hz
    simp only [Set.mem_setOf_eq] at hz
    show P.IsRoot z
    simp [Polynomial.IsRoot, P, Polynomial.eval_add, Polynomial.eval_mul,
          Polynomial.eval_C, Polynomial.eval_X, Polynomial.eval_pow]
    linear_combination hz

/-- **Non-vanishing of `sigma_p2_explicit_C` fiber polynomial
    coefficients on ℂ, for `x ≠ 1`.** -/
theorem sigma_p2_fiber_coeffs_nonzero_C
    (x y : ℂ) (hx_ne_one : x - 1 ≠ 0) (hx_ne : x ≠ 0) :
    ((x + 1) - 2 * x * y) ≠ 0 ∨ (4 - 2 * (x + 1) * y) ≠ 0 := by
  by_cases hA : (x + 1) - 2 * x * y = 0
  · right
    intro hB
    -- From hA: y = (x+1)/(2x).
    have h2x_ne : (2 * x : ℂ) ≠ 0 := mul_ne_zero two_ne_zero hx_ne
    have hy_eq : y = (x + 1) / (2 * x) := by
      field_simp
      linear_combination -hA
    -- Substitute into hB: 4 - 2(x+1)·(x+1)/(2x) = 0 ⟹ (x-1)² = 0.
    rw [hy_eq] at hB
    have hB' : 4 * (2 * x) - 2 * (x + 1) * (x + 1) = 0 := by
      field_simp at hB
      linear_combination hB
    have h_sq : (x - 1) ^ 2 = 0 := by
      linear_combination -hB' / 2
    have h_x_1 : x - 1 = 0 := by
      have := pow_eq_zero_iff (n := 2) (by norm_num : (2:ℕ) ≠ 0) |>.mp h_sq
      exact this
    exact hx_ne_one h_x_1
  · left; exact hA

/-- **Fiber of `sigma_p2_explicit_C` at any `y ∈ ℂ` is finite**
    (for `x ≠ 0`, `x ≠ 1`, `xα² = 1`). -/
theorem sigma_p2_explicit_C_fiber_finite
    (x α y : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    Set.Finite {z : ℂ | sigma_p2_explicit_C x α z = y} := by
  have h_bound : Set.Finite
    (({(-1 : ℂ)} ∪ {(-1/x : ℂ)})
      ∪ {z : ℂ | ((x+1) - 2*x*y) * z^2 + (4 - 2*(x+1)*y) * z + ((x+1)/x - 2*y) = 0}) := by
    refine Set.Finite.union ?_ ?_
    · exact (Set.finite_singleton (-1:ℂ)).union (Set.finite_singleton (-1/x:ℂ))
    · exact poly_deg_le_2_zero_set_finite_C _ _ _
        (sigma_p2_fiber_coeffs_nonzero_C x y hx_ne_one hx_ne)
  refine h_bound.subset ?_
  intro z hz
  simp only [Set.mem_setOf_eq] at hz
  by_cases hz1 : z + 1 = 0
  · left; left
    simp only [Set.mem_singleton_iff]
    linear_combination hz1
  by_cases hxz1 : x * z + 1 = 0
  · left; right
    simp only [Set.mem_singleton_iff]
    field_simp
    linear_combination hxz1
  · right
    simp only [Set.mem_setOf_eq]
    unfold sigma_p2_explicit_C at hz
    have h_denom_ne : (2 * (z + 1) * (x * z + 1) : ℂ) ≠ 0 :=
      mul_ne_zero (mul_ne_zero two_ne_zero hz1) hxz1
    have hz' : (z - y) * (2 * (z + 1) * (x * z + 1))
              = (z - α) * (z + α) * (2 * x * z + x + 1) := by
      have := hz
      field_simp at this
      linear_combination this
    have h_x_target : x * ((x + 1) - 2*x*y) * z^2 + x * (4 - 2*(x+1)*y) * z
                     + ((x+1) - 2*x*y) = 0 := by
      linear_combination x * hz' - (2*x*z + x + 1) * hxα
    have h_target_form : x * (((x+1) - 2*x*y) * z^2 + (4 - 2*(x+1)*y) * z + ((x+1)/x - 2*y))
                       = x * ((x + 1) - 2*x*y) * z^2 + x * (4 - 2*(x+1)*y) * z + ((x+1) - 2*x*y) := by
      field_simp
      ring
    have h_x_zero : x * (((x+1) - 2*x*y) * z^2 + (4 - 2*(x+1)*y) * z + ((x+1)/x - 2*y)) = 0 := by
      rw [h_target_form]; exact h_x_target
    exact mul_left_cancel₀ hx_ne (h_x_zero.trans (mul_zero x).symm)

/-- **Fiber of `steffensen_step_C x 2` at `y ∈ ℂ` is finite.** -/
theorem steffensen_step_C_fiber_finite
    (x α y : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    Set.Finite {z : ℂ | steffensen_step_C x 2 z = y} := by
  have h_degen : Set.Finite
      ({(-1 : ℂ), -1/x, -(x+1)/(2*x), α, -α} : Set ℂ) := by
    apply Set.Finite.insert
    apply Set.Finite.insert
    apply Set.Finite.insert
    apply Set.Finite.insert
    exact Set.finite_singleton _
  have h_sigma : Set.Finite {z : ℂ | sigma_p2_explicit_C x α z = y} :=
    sigma_p2_explicit_C_fiber_finite x α y hx_ne hx_ne_one hxα
  refine (h_degen.union h_sigma).subset ?_
  intro z hz
  simp only [Set.mem_setOf_eq] at hz
  by_cases hz1 : z + 1 = 0
  · left
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    left; linear_combination hz1
  by_cases hxz1 : x * z + 1 = 0
  · left
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    right; left
    field_simp; linear_combination hxz1
  by_cases hq : 2 * x * z + x + 1 = 0
  · left
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    right; right; left
    have h2x_ne : (2 * x : ℂ) ≠ 0 := mul_ne_zero two_ne_zero hx_ne
    field_simp
    linear_combination hq
  by_cases hzα : z = α
  · left
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    right; right; right; left
    exact hzα
  by_cases hz_negα : z = -α
  · left
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    right; right; right; right
    exact hz_negα
  · right
    simp only [Set.mem_setOf_eq]
    rw [← steffensen_step_C_p2_eq_sigma_p2_explicit_C x α hx_ne hxα z hz1 hxz1 hq hzα hz_negα]
    exact hz

/-- **Fiber of `(steffensen_step_C x 2)^[n]` at any `y ∈ ℂ` is finite**. -/
theorem steffensen_step_C_iterate_fiber_finite
    (x α : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (hxα : x * α ^ 2 = 1) (n : ℕ) :
    ∀ y : ℂ, Set.Finite {z₀ : ℂ | (steffensen_step_C x 2)^[n] z₀ = y} := by
  induction n with
  | zero =>
    intro y
    simp only [Function.iterate_zero_apply]
    refine (Set.finite_singleton y).subset ?_
    intro z₀ hz₀
    simp only [Set.mem_setOf_eq] at hz₀
    exact hz₀
  | succ k ih =>
    intro y
    have h_fiber_fin : Set.Finite {z : ℂ | steffensen_step_C x 2 z = y} :=
      steffensen_step_C_fiber_finite x α y hx_ne hx_ne_one hxα
    have h_union_fin : Set.Finite
      (⋃ y' ∈ {z : ℂ | steffensen_step_C x 2 z = y},
         {z₀ : ℂ | (steffensen_step_C x 2)^[k] z₀ = y'}) :=
      h_fiber_fin.biUnion (fun y' _ => ih y')
    refine h_union_fin.subset ?_
    intro z₀ hz₀
    simp only [Set.mem_setOf_eq] at hz₀
    rw [Function.iterate_succ_apply'] at hz₀
    rw [Set.mem_iUnion]
    refine ⟨(steffensen_step_C x 2)^[k] z₀, ?_⟩
    rw [Set.mem_iUnion]
    exact ⟨hz₀, rfl⟩

end FiberFinitenessC

/-! ============================================================
  §39.5  Complex bad set has Lebesgue measure zero
============================================================ -/

section ComplexBadSetC

/-- **Complex "bad set" for unconditional `McMullenAEEntry 2 x α`.**

    Union of:
      • the Julia circle `{z : ‖v(z)‖ = 1}` (Lebesgue-null, §39.2),
      • orbit hits of the Möbius singularities `−1, −1/x, −(x+1)/(2x)`
        (countable union of finite fibers under iterated `σ`).

    The `±α` hits handled separately in the dichotomy (direct basin
    entry). -/
def complex_bad_set (x α : ℂ) : Set ℂ :=
  {z₀ | ∃ n : ℕ,
        (steffensen_step_C x 2)^[n] z₀ = -1
      ∨ (steffensen_step_C x 2)^[n] z₀ = -1 / x
      ∨ 2 * x * ((steffensen_step_C x 2)^[n] z₀) + x + 1 = 0}
    ∪ {z₀ | ‖v_p2_C α z₀‖ = 1}

/-- **★★★ Complex bad set has Lebesgue measure zero.**

    First component: countable union over `n` of finite fibers →
    countable → null.

    Second component: Julia circle → null by §39.2
    `julia_section_complex_volume_zero`. -/
theorem complex_bad_set_volume_zero
    (x α : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (hα_ne_zero : α ≠ 0) (hxα : x * α ^ 2 = 1) :
    volume (complex_bad_set x α) = 0 := by
  unfold complex_bad_set
  refine measure_union_null ?_ ?_
  · -- Pole preimages: countable.
    set F : Set ℂ := ({(-1:ℂ), -1/x, -(x+1)/(2*x)} : Set ℂ) with hF_def
    have hx_ne' : (2 * x : ℂ) ≠ 0 := mul_ne_zero two_ne_zero hx_ne
    have h_iter_finite : ∀ n : ℕ, Set.Finite
        {z₀ : ℂ | (steffensen_step_C x 2)^[n] z₀ ∈ F} := by
      intro n
      have hF_finite : F.Finite := ((Set.finite_singleton _).insert _).insert _
      have h_union_form : {z₀ : ℂ | (steffensen_step_C x 2)^[n] z₀ ∈ F}
                        = ⋃ y ∈ F, {z₀ : ℂ | (steffensen_step_C x 2)^[n] z₀ = y} := by
        ext z₀
        constructor
        · intro hz
          rw [Set.mem_iUnion₂]
          exact ⟨(steffensen_step_C x 2)^[n] z₀, hz, rfl⟩
        · intro hz
          rw [Set.mem_iUnion₂] at hz
          obtain ⟨y, hyF, hyz⟩ := hz
          simp only [Set.mem_setOf_eq] at hyz
          rw [Set.mem_setOf_eq, hyz]
          exact hyF
      rw [h_union_form]
      exact hF_finite.biUnion
        (fun y _ => steffensen_step_C_iterate_fiber_finite x α hx_ne hx_ne_one hxα n y)
    have h_countable : Set.Countable
        {z₀ | ∃ n : ℕ,
              (steffensen_step_C x 2)^[n] z₀ = -1
            ∨ (steffensen_step_C x 2)^[n] z₀ = -1 / x
            ∨ 2 * x * ((steffensen_step_C x 2)^[n] z₀) + x + 1 = 0} := by
      have h_subset :
          {z₀ | ∃ n : ℕ,
                (steffensen_step_C x 2)^[n] z₀ = -1
              ∨ (steffensen_step_C x 2)^[n] z₀ = -1 / x
              ∨ 2 * x * ((steffensen_step_C x 2)^[n] z₀) + x + 1 = 0} ⊆
          ⋃ n : ℕ, {z₀ : ℂ | (steffensen_step_C x 2)^[n] z₀ ∈ F} := by
        intro z₀ hz₀
        obtain ⟨n, hn⟩ := hz₀
        rw [Set.mem_iUnion]
        refine ⟨n, ?_⟩
        simp only [Set.mem_setOf_eq, hF_def, Set.mem_insert_iff, Set.mem_singleton_iff]
        rcases hn with h1 | h2 | h3
        · left; exact h1
        · right; left; exact h2
        · right; right
          have h_eq : 2 * x * ((steffensen_step_C x 2)^[n] z₀) = -(x + 1) := by
            linear_combination h3
          rw [eq_div_iff hx_ne']
          linear_combination h_eq
      refine Set.Countable.mono h_subset ?_
      exact Set.countable_iUnion (fun n => (h_iter_finite n).countable)
    exact h_countable.measure_zero _
  · -- Julia circle: null by §39.2.
    exact julia_section_complex_volume_zero α hα_ne_zero

end ComplexBadSetC

/-! ============================================================
  §39.6  Orbit non-degeneracy and bridge along the orbit
============================================================ -/

section OrbitNonDegenC

/-- **Orbit non-degeneracy unpacked from `z₀ ∉ complex_bad_set`.** -/
theorem complex_orbit_non_degenerate_of_not_bad_set
    (x α : ℂ) (z₀ : ℂ) (hgood : z₀ ∉ complex_bad_set x α) :
    (∀ n : ℕ,
      (steffensen_step_C x 2)^[n] z₀ ≠ -1
      ∧ (steffensen_step_C x 2)^[n] z₀ ≠ -1 / x
      ∧ 2 * x * ((steffensen_step_C x 2)^[n] z₀) + x + 1 ≠ 0) ∧
    ‖v_p2_C α z₀‖ ≠ 1 := by
  unfold complex_bad_set at hgood
  simp only [Set.mem_union, Set.mem_setOf_eq, not_or, not_exists, not_or] at hgood
  refine ⟨fun n => ?_, hgood.2⟩
  exact hgood.1 n

/-- **Bridge along the orbit (excluding `±α`-hits).**

    If no iterate of `z₀` equals `±α` and the orbit stays in the
    non-degenerate set, then `steffensen_step_C x 2` and
    `sigma_p2_explicit_C x α` agree at every iterate. -/
theorem steffensen_C_eq_sigma_p2_iterate_of_orbit_good
    (x α : ℂ) (hx_ne : x ≠ 0) (hxα : x * α ^ 2 = 1) (z₀ : ℂ)
    (H : ∀ n : ℕ,
         (steffensen_step_C x 2)^[n] z₀ ≠ -1
         ∧ (steffensen_step_C x 2)^[n] z₀ ≠ -1 / x
         ∧ 2 * x * ((steffensen_step_C x 2)^[n] z₀) + x + 1 ≠ 0
         ∧ (steffensen_step_C x 2)^[n] z₀ ≠ α
         ∧ (steffensen_step_C x 2)^[n] z₀ ≠ -α) :
    ∀ n : ℕ, (steffensen_step_C x 2)^[n] z₀
             = (sigma_p2_explicit_C x α)^[n] z₀ := by
  intro n
  induction n with
  | zero => rfl
  | succ k ih =>
    obtain ⟨h1, h2, h3, h4, h5⟩ := H k
    have hz : (steffensen_step_C x 2)^[k] z₀ + 1 ≠ 0 := by
      intro h; apply h1; linear_combination h
    have hxz1 : x * ((steffensen_step_C x 2)^[k] z₀) + 1 ≠ 0 := by
      intro h
      apply h2
      field_simp
      linear_combination h
    have h_bridge := steffensen_step_C_p2_eq_sigma_p2_explicit_C
      x α hx_ne hxα ((steffensen_step_C x 2)^[k] z₀) hz hxz1 h3 h4 h5
    rw [Function.iterate_succ_apply', Function.iterate_succ_apply', h_bridge, ih]

end OrbitNonDegenC

/-! ============================================================
  §39.7  Complex dichotomy outside the bad set
============================================================ -/

section ComplexDichotomyC

/-- **★★★ Convergence dichotomy outside the complex bad set.**

    For every `z₀ ∉ complex_bad_set x α` (and `α ≠ 0, α ≠ ±1`, fixed
    point at `α² = 1/x`), the orbit under `steffensen_step_C x 2`
    eventually enters any prescribed neighbourhood of `+α` or `−α`.

    **Proof structure.**
      (i)   Orbit non-degeneracy via `complex_bad_set` unpacking.
      (ii)  ±α-hit cases dispatched trivially.
      (iii) Bridge `steffensen_step_C x 2 ≡ sigma_p2_explicit_C x α`
            along the orbit.
      (iv)  Iterated Böttcher `‖v(σⁿ z₀)‖ = ‖v(z₀)‖^{2ⁿ}`.
      (v)   Initial dichotomy `‖v(z₀)‖ < 1` vs `> 1` (excluded `= 1`
            by z₀ ∉ bad set).
      (vi)  Möbius-inverse bounds `‖s − α‖ ≤ 4‖α‖‖v(s)‖/‖μ²‖`
            (when ‖v(s)‖ ≤ ‖μ²‖/2) and dual for `s + α`.
      (vii) Squaring convergence (`tendsto_pow_atTop_nhds_zero_of_lt_one`
            composed with `Nat.tendsto_pow_atTop_atTop_of_one_lt`)
            and basin entry at the selected threshold. -/
theorem complex_orbit_enters_basin_off_bad_set
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (z₀ : ℂ) (hgood : z₀ ∉ complex_bad_set x α)
    (r_pos r_neg : ℝ) (hr_pos : 0 < r_pos) (hr_neg : 0 < r_neg) :
    (∃ k₀ : ℕ, ‖(steffensen_step_C x 2)^[k₀] z₀ - α‖ < r_pos) ∨
    (∃ k₀ : ℕ, ‖(steffensen_step_C x 2)^[k₀] z₀ - (-α)‖ < r_neg) := by
  have hxα : x * α ^ 2 = 1 := by rw [hα_pow]; field_simp
  obtain ⟨h_orbit_nondeg, h_v_ne_one⟩ :=
    complex_orbit_non_degenerate_of_not_bad_set x α z₀ hgood
  -- Dispatch `±α`-hit cases.
  by_cases h_hit_α : ∃ k, (steffensen_step_C x 2)^[k] z₀ = α
  · obtain ⟨k, hk⟩ := h_hit_α
    left
    refine ⟨k, ?_⟩
    rw [hk]
    simp; exact hr_pos
  by_cases h_hit_negα : ∃ k, (steffensen_step_C x 2)^[k] z₀ = -α
  · obtain ⟨k, hk⟩ := h_hit_negα
    right
    refine ⟨k, ?_⟩
    rw [hk]
    simp; exact hr_neg
  push_neg at h_hit_α h_hit_negα
  -- Main case: apply bridge.
  have H_bridge_data : ∀ n : ℕ,
      (steffensen_step_C x 2)^[n] z₀ ≠ -1
      ∧ (steffensen_step_C x 2)^[n] z₀ ≠ -1 / x
      ∧ 2 * x * ((steffensen_step_C x 2)^[n] z₀) + x + 1 ≠ 0
      ∧ (steffensen_step_C x 2)^[n] z₀ ≠ α
      ∧ (steffensen_step_C x 2)^[n] z₀ ≠ -α := by
    intro n
    obtain ⟨h1, h2, h3⟩ := h_orbit_nondeg n
    exact ⟨h1, h2, h3, h_hit_α n, h_hit_negα n⟩
  have h_orbit_eq : ∀ n : ℕ,
      (steffensen_step_C x 2)^[n] z₀
        = (sigma_p2_explicit_C x α)^[n] z₀ :=
    steffensen_C_eq_sigma_p2_iterate_of_orbit_good
      x α hx_ne hxα z₀ H_bridge_data
  -- Non-degeneracy along sigma-iterates.
  have H_vpow : ∀ n : ℕ, ∀ k ≤ n,
      (sigma_p2_explicit_C x α)^[k] z₀ + 1 ≠ 0
      ∧ x * (sigma_p2_explicit_C x α)^[k] z₀ + 1 ≠ 0
      ∧ (sigma_p2_explicit_C x α)^[k] z₀ + α ≠ 0 := by
    intro n k _hkn
    obtain ⟨h1, h2, _, _, h5⟩ := H_bridge_data k
    rw [← h_orbit_eq k] at *
    refine ⟨?_, ?_, ?_⟩
    · intro h; apply h1; linear_combination h
    · intro h; apply h2
      field_simp
      linear_combination h
    · intro h; apply h5; linear_combination h
  -- Iterated Böttcher: `v(σⁿ z₀) = v(z₀)^{2ⁿ}`.
  have h_v_iter : ∀ n : ℕ,
      v_p2_C α ((sigma_p2_explicit_C x α)^[n] z₀)
        = (v_p2_C α z₀) ^ (2 ^ n) := by
    intro n
    exact v_p2_iterated_C x α hα_ne_zero hα_ne_neg_one hxα z₀ n (H_vpow n)
  -- Möbius-inverse identities on the non-degenerate set.
  set μ_sq : ℂ := ((1 - α) / (1 + α)) ^ 2 with hμ_sq_def
  have hμ_sq_ne_zero : μ_sq ≠ 0 := by
    rw [hμ_sq_def]
    apply pow_ne_zero
    exact div_ne_zero hα_ne_one hα_ne_neg_one
  have hμ_sq_norm_pos : 0 < ‖μ_sq‖ := norm_pos_iff.mpr hμ_sq_ne_zero
  have hα_norm_pos : 0 < ‖α‖ := norm_pos_iff.mpr hα_ne_zero
  -- `(s - α)·(μ² - v(s)) = 2α·v(s)`  on `s + α ≠ 0`.
  have h_mobius_minus :
      ∀ s : ℂ, s + α ≠ 0 →
        (s - α) * (μ_sq - v_p2_C α s) = 2 * α * v_p2_C α s := by
    intro s hsα
    unfold v_p2_C
    rw [← hμ_sq_def]
    field_simp
    ring
  -- `(s + α)·(μ² - v(s)) = 2α·μ²`  on `s + α ≠ 0`.
  have h_mobius_plus :
      ∀ s : ℂ, s + α ≠ 0 →
        (s + α) * (μ_sq - v_p2_C α s) = 2 * α * μ_sq := by
    intro s hsα
    unfold v_p2_C
    rw [← hμ_sq_def]
    field_simp
    ring
  -- `μ² − v(s) ≠ 0` whenever `s + α ≠ 0`.
  have h_mu_sq_minus_v_ne :
      ∀ s : ℂ, s + α ≠ 0 → μ_sq - v_p2_C α s ≠ 0 := by
    intro s hsα h_eq
    have h := h_mobius_plus s hsα
    rw [h_eq, mul_zero] at h
    have hne : (2 : ℂ) * α * μ_sq ≠ 0 := by
      apply mul_ne_zero (mul_ne_zero _ hα_ne_zero) hμ_sq_ne_zero
      exact two_ne_zero
    exact hne h.symm
  -- Per-iterate: σ'^[n] z₀ + α ≠ 0.
  have h_sigma_plus_α_ne :
      ∀ n : ℕ, (sigma_p2_explicit_C x α)^[n] z₀ + α ≠ 0 := fun n =>
    (H_vpow n n (le_refl n)).2.2
  -- Dichotomy on ‖v(z₀)‖.
  rcases lt_or_gt_of_ne h_v_ne_one with h_v_lt | h_v_gt
  · -- Case `‖v(z₀)‖ < 1`:  orbit converges to `+α`.
    left
    set v₀ : ℂ := v_p2_C α z₀
    have hv₀_abs_lt : ‖v₀‖ < 1 := h_v_lt
    have hv₀_abs_nonneg : (0 : ℝ) ≤ ‖v₀‖ := norm_nonneg _
    set δ : ℝ := min (‖μ_sq‖ / 2) (r_pos * ‖μ_sq‖ / (4 * ‖α‖)) with hδ_def
    have h4α_pos : (0 : ℝ) < 4 * ‖α‖ := by
      apply mul_pos; · norm_num
      exact hα_norm_pos
    have hδ_pos : 0 < δ := by
      rw [hδ_def]
      apply lt_min
      · exact half_pos hμ_sq_norm_pos
      · exact div_pos (mul_pos hr_pos hμ_sq_norm_pos) h4α_pos
    have hδ_le_half : δ ≤ ‖μ_sq‖ / 2 := min_le_left _ _
    have hδ_le_target : δ ≤ r_pos * ‖μ_sq‖ / (4 * ‖α‖) := min_le_right _ _
    -- `‖v₀‖^m → 0`.
    have h_tend_pow :
        Filter.Tendsto (fun m : ℕ => ‖v₀‖ ^ m) Filter.atTop (nhds 0) :=
      tendsto_pow_atTop_nhds_zero_of_lt_one hv₀_abs_nonneg hv₀_abs_lt
    have h_tend_2n :
        Filter.Tendsto (fun n : ℕ => 2 ^ n) Filter.atTop Filter.atTop :=
      Nat.tendsto_pow_atTop_atTop_of_one_lt (by norm_num : (1 : ℕ) < 2)
    have h_tend_comp :
        Filter.Tendsto (fun n : ℕ => ‖v₀‖ ^ (2 ^ n)) Filter.atTop (nhds 0) :=
      h_tend_pow.comp h_tend_2n
    have h_event : ∀ᶠ n in Filter.atTop, ‖v₀‖ ^ (2 ^ n) < δ := by
      have h_mem : Set.Ioo (-δ) δ ∈ nhds (0 : ℝ) :=
        Ioo_mem_nhds (by linarith) hδ_pos
      have h_ev := h_tend_comp h_mem
      filter_upwards [h_ev] with n hn using hn.2
    obtain ⟨k, hk⟩ := h_event.exists
    set s_k : ℂ := (sigma_p2_explicit_C x α)^[k] z₀ with hs_k_def
    have h_v_iter_k : v_p2_C α s_k = v₀ ^ (2 ^ k) := h_v_iter k
    have h_abs_v_sk_eq : ‖v_p2_C α s_k‖ = ‖v₀‖ ^ (2 ^ k) := by
      rw [h_v_iter_k, norm_pow]
    have h_abs_v_sk_lt : ‖v_p2_C α s_k‖ < δ := by
      rw [h_abs_v_sk_eq]; exact hk
    have hsα_k : s_k + α ≠ 0 := by
      rw [hs_k_def]; exact h_sigma_plus_α_ne k
    have h_mu_v_ne_k : μ_sq - v_p2_C α s_k ≠ 0 :=
      h_mu_sq_minus_v_ne s_k hsα_k
    have h_mu_v_abs_pos : 0 < ‖μ_sq - v_p2_C α s_k‖ := norm_pos_iff.mpr h_mu_v_ne_k
    -- From `h_mobius_minus`: `s_k - α = 2α·v/(μ² - v)`.
    have h_sk_minus_α_eq :
        s_k - α = 2 * α * v_p2_C α s_k / (μ_sq - v_p2_C α s_k) := by
      have h := h_mobius_minus s_k hsα_k
      field_simp
      linear_combination h
    -- `‖s_k - α‖ = 2‖α‖‖v(s_k)‖/‖μ² - v(s_k)‖`.
    have h_abs_sk_minus_α :
        ‖s_k - α‖ = 2 * ‖α‖ * ‖v_p2_C α s_k‖ / ‖μ_sq - v_p2_C α s_k‖ := by
      rw [h_sk_minus_α_eq, norm_div]
      rw [show (2 : ℂ) * α * v_p2_C α s_k = 2 * (α * v_p2_C α s_k) from by ring,
          norm_mul, norm_mul, show ‖(2 : ℂ)‖ = 2 from by
            rw [show ((2 : ℂ)) = ((2 : ℕ) : ℂ) from by norm_num, Complex.norm_nat]
            norm_num]
      ring
    -- Denominator bound: `‖μ² - v‖ ≥ ‖μ²‖/2` when `‖v‖ ≤ ‖μ²‖/2`.
    have h_denom_lower :
        ‖μ_sq‖ / 2 ≤ ‖μ_sq - v_p2_C α s_k‖ := by
      have h_v_small : ‖v_p2_C α s_k‖ ≤ ‖μ_sq‖ / 2 :=
        le_of_lt (lt_of_lt_of_le h_abs_v_sk_lt hδ_le_half)
      have h_triangle : ‖μ_sq‖ - ‖v_p2_C α s_k‖ ≤ ‖μ_sq - v_p2_C α s_k‖ := by
        exact norm_sub_norm_le μ_sq (v_p2_C α s_k)
      linarith
    -- Final bound.
    have h_final : ‖s_k - α‖ < r_pos := by
      rw [h_abs_sk_minus_α, div_lt_iff h_mu_v_abs_pos]
      have h_δ_bd : δ * (4 * ‖α‖) ≤ r_pos * ‖μ_sq‖ := by
        rw [le_div_iff h4α_pos] at hδ_le_target
        linarith
      have h_v_abs_nn : (0 : ℝ) ≤ ‖v_p2_C α s_k‖ := norm_nonneg _
      have h_α_nn : (0 : ℝ) ≤ ‖α‖ := le_of_lt hα_norm_pos
      nlinarith [h_abs_v_sk_lt, h_denom_lower, h_δ_bd, hμ_sq_norm_pos,
                 hr_pos, h_α_nn, h_v_abs_nn, h4α_pos]
    refine ⟨k, ?_⟩
    rw [h_orbit_eq k, ← hs_k_def]
    exact h_final
  · -- Case `‖v(z₀)‖ > 1`:  orbit converges to `−α`.
    right
    set v₀ : ℂ := v_p2_C α z₀
    have hv₀_abs_gt : 1 < ‖v₀‖ := h_v_gt
    set M : ℝ := max (2 * ‖μ_sq‖) (4 * ‖α‖ * ‖μ_sq‖ / r_neg + 1) with hM_def
    have hM_pos : 0 < M := by
      rw [hM_def]
      apply lt_of_lt_of_le _ (le_max_right _ _)
      have : (0 : ℝ) ≤ 4 * ‖α‖ * ‖μ_sq‖ / r_neg := by positivity
      linarith
    have hM_ge_2μ : 2 * ‖μ_sq‖ ≤ M := le_max_left _ _
    have hM_ge_target : 4 * ‖α‖ * ‖μ_sq‖ / r_neg + 1 ≤ M := le_max_right _ _
    have h_tend_pow :
        Filter.Tendsto (fun m : ℕ => ‖v₀‖ ^ m) Filter.atTop Filter.atTop :=
      tendsto_pow_atTop_atTop_of_one_lt hv₀_abs_gt
    have h_tend_2n :
        Filter.Tendsto (fun n : ℕ => 2 ^ n) Filter.atTop Filter.atTop :=
      Nat.tendsto_pow_atTop_atTop_of_one_lt (by norm_num : (1 : ℕ) < 2)
    have h_tend_comp :
        Filter.Tendsto (fun n : ℕ => ‖v₀‖ ^ (2 ^ n)) Filter.atTop Filter.atTop :=
      h_tend_pow.comp h_tend_2n
    have h_event : ∀ᶠ n in Filter.atTop, M < ‖v₀‖ ^ (2 ^ n) := by
      exact h_tend_comp (Filter.eventually_gt_atTop M)
    obtain ⟨k, hk⟩ := h_event.exists
    set s_k : ℂ := (sigma_p2_explicit_C x α)^[k] z₀ with hs_k_def
    have h_v_iter_k : v_p2_C α s_k = v₀ ^ (2 ^ k) := h_v_iter k
    have h_abs_v_sk_eq : ‖v_p2_C α s_k‖ = ‖v₀‖ ^ (2 ^ k) := by
      rw [h_v_iter_k, norm_pow]
    have h_abs_v_sk_gt : M < ‖v_p2_C α s_k‖ := by
      rw [h_abs_v_sk_eq]; exact hk
    have hsα_k : s_k + α ≠ 0 := by
      rw [hs_k_def]; exact h_sigma_plus_α_ne k
    have h_mu_v_ne_k : μ_sq - v_p2_C α s_k ≠ 0 :=
      h_mu_sq_minus_v_ne s_k hsα_k
    have h_mu_v_abs_pos : 0 < ‖μ_sq - v_p2_C α s_k‖ := norm_pos_iff.mpr h_mu_v_ne_k
    have h_sk_plus_α_eq :
        s_k + α = 2 * α * μ_sq / (μ_sq - v_p2_C α s_k) := by
      have h := h_mobius_plus s_k hsα_k
      field_simp
      linear_combination h
    have h_abs_sk_plus_α :
        ‖s_k + α‖ = 2 * ‖α‖ * ‖μ_sq‖ / ‖μ_sq - v_p2_C α s_k‖ := by
      rw [h_sk_plus_α_eq, norm_div]
      rw [show (2 : ℂ) * α * μ_sq = 2 * (α * μ_sq) from by ring,
          norm_mul, norm_mul, show ‖(2 : ℂ)‖ = 2 from by
            rw [show ((2 : ℂ)) = ((2 : ℕ) : ℂ) from by norm_num, Complex.norm_nat]
            norm_num]
      ring
    have h_v_ge_2μ : 2 * ‖μ_sq‖ ≤ ‖v_p2_C α s_k‖ :=
      le_of_lt (lt_of_le_of_lt hM_ge_2μ h_abs_v_sk_gt)
    have h_denom_lower : ‖v_p2_C α s_k‖ / 2 ≤ ‖μ_sq - v_p2_C α s_k‖ := by
      have h_triangle : ‖v_p2_C α s_k‖ - ‖μ_sq‖ ≤ ‖μ_sq - v_p2_C α s_k‖ := by
        have : ‖v_p2_C α s_k - μ_sq‖ = ‖μ_sq - v_p2_C α s_k‖ := norm_sub_rev _ _
        rw [← this]
        exact norm_sub_norm_le (v_p2_C α s_k) μ_sq
      linarith
    have h_final : ‖s_k - -α‖ < r_neg := by
      have h_goal_rw : ‖s_k - -α‖ = ‖s_k + α‖ := by congr 1; ring
      rw [h_goal_rw, h_abs_sk_plus_α, div_lt_iff h_mu_v_abs_pos]
      have hM_gt_4α : 4 * ‖α‖ * ‖μ_sq‖ / r_neg < M := by
        have := hM_ge_target; linarith
      rw [div_lt_iff hr_neg] at hM_gt_4α
      have h_v_pos : 0 < ‖v_p2_C α s_k‖ := by linarith
      have h_α_nn : (0 : ℝ) ≤ ‖α‖ := le_of_lt hα_norm_pos
      nlinarith [h_abs_v_sk_gt, h_denom_lower, hM_gt_4α, hμ_sq_norm_pos,
                 hr_neg, h_α_nn, h_v_pos, hM_pos]
    refine ⟨k, ?_⟩
    rw [h_orbit_eq k, ← hs_k_def]
    exact h_final

end ComplexDichotomyC

/-! ============================================================
  §39.8  Final: `McMullenAEEntry 2 x α` unconditionally
============================================================ -/

section FinalMcMullenC

/-- **★★★★★ Unconditional `McMullenAEEntry 2 x α` on ℂ.**

    For every `x ∈ ℂ \ {0}` and every `α ≠ 0` with `α² = 1/x` and
    `α ≠ ±1`, the §32.1 named `Prop` `McMullenAEEntry 2 x α` holds
    **fully unconditionally**. Lebesgue-almost every `z₀ ∈ ℂ`
    eventually lands in any prescribed neighbourhood of some
    cyclotomic anchor `γ_s ∈ {α, −α}` under `steffensen_step_C x 2`.

    This is the **paper-citable, unconditional, complex-plane
    closure** of the Pandrosion-Steffensen `p = 2` solver. Combined
    with `steffensen_solves_ae_mod_mcmullen` (§32.2) and
    `steffensen_global_loglog_ae_mod_mcmullen` (§34.2), it yields
    the fully unconditional global log-log complexity theorem:

    *For every `x ∈ ℂ \ {0}`, Pandrosion-Steffensen solves `z² = x`
    to precision ε in `k₀(z₀) + O(log log 1/ε)` iterations, for
    Lebesgue-almost every `z₀ ∈ ℂ`.* -/
theorem mcmullen_p2_complex_unconditional
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x) :
    McMullenAEEntry 2 x α := by
  intro r hr_pos
  have hxα : x * α ^ 2 = 1 := by rw [hα_pow]; field_simp
  have hx_ne_one : x - 1 ≠ 0 := by
    intro h
    have hx_eq : x = 1 := by linear_combination h
    rw [hx_eq, one_mul] at hxα
    have h_α_sq : α ^ 2 = 1 := hxα
    rcases sq_eq_sq_iff_eq_or_eq_neg.mp (by rw [h_α_sq]; ring : α ^ 2 = (1 : ℂ) ^ 2) with
      h_α1 | h_αm1
    · apply hα_ne_one
      linear_combination -h_α1
    · apply hα_ne_neg_one
      linear_combination h_αm1 - (1 : ℂ)
  have h_bad_null : volume (complex_bad_set x α) = 0 :=
    complex_bad_set_volume_zero x α hx_ne hx_ne_one hα_ne_zero hxα
  -- Reduce `∀ᵐ z₀, P z₀` to `volume {z₀ | ¬P z₀} = 0`.
  rw [MeasureTheory.ae_iff]
  refine MeasureTheory.measure_mono_null ?_ h_bad_null
  intro z₀ hz₀
  -- `hz₀`: orbit enters no basin of any cycAnchor.
  -- Show `z₀ ∈ complex_bad_set` by contrapositive.
  by_contra hgood
  apply hz₀
  -- Dispatch cyclotomic anchors: γ_0 = α, γ_1 = -α.
  have h_γ0 : cycAnchor α 2 (0 : Fin 2) = α := by
    unfold cycAnchor
    simp [Fin.val_zero, Nat.cast_zero, mul_zero, zero_div, Complex.exp_zero, mul_one]
  have h_γ1 : cycAnchor α 2 (1 : Fin 2) = -α := by
    unfold cycAnchor
    have h_val : (((1 : Fin 2) : ℕ) : ℂ) = 1 := by norm_num
    rw [h_val]
    have h_two : ((2 : ℕ) : ℂ) = 2 := by norm_num
    rw [h_two]
    have h_simp : (2 : ℂ) * (Real.pi : ℂ) * Complex.I * 1 / 2 =
                  (Real.pi : ℂ) * Complex.I := by
      ring
    rw [h_simp, Complex.exp_pi_mul_I]
    ring
  -- Apply dichotomy.
  obtain h_α | h_negα := complex_orbit_enters_basin_off_bad_set
    x hx_ne α hα_ne_zero hα_ne_neg_one hα_ne_one hα_pow z₀ hgood
    (r 0) (r 1) (hr_pos 0) (hr_pos 1)
  · obtain ⟨k₀, hk₀⟩ := h_α
    refine ⟨0, k₀, ?_⟩
    rw [h_γ0]; exact hk₀
  · obtain ⟨k₀, hk₀⟩ := h_negα
    refine ⟨1, k₀, ?_⟩
    rw [h_γ1]; exact hk₀

end FinalMcMullenC

/-! ============================================================
  §39.9  Fully unconditional `Pandrosion-Steffensen solves z² = x`
============================================================ -/

section FinalUnconditionalSolveC

/-- **★★★★★★ Pandrosion-Steffensen solves `z² = x` almost everywhere
    on ℂ, fully unconditionally.**

    For every `x ∈ ℂ \ {0}` and every `α ≠ 0` with `α² = 1/x`
    (and `α ≠ ±1`), and every choice of Sp-non-degeneracy at the
    cyclotomic anchors, for Lebesgue-almost every `z₀ ∈ ℂ`, the
    iterates `(steffensen_step_C x 2)^[k] z₀` converge to some
    cyclotomic root `γ_s ∈ {α, −α}`.

    This chains `mcmullen_p2_complex_unconditional` (§39.8) with
    `steffensen_solves_ae_mod_mcmullen` (§32.2). -/
theorem steffensen_p2_solves_complex_unconditional
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (hSp : ∀ s : Fin 2, Sp_C 2 (cycAnchor α 2 s) ≠ 0) :
    ∀ᵐ z₀ : ℂ ∂volume,
      ∃ s : Fin 2,
        Filter.Tendsto (fun k => (steffensen_step_C x 2)^[k] z₀) Filter.atTop
          (nhds (cycAnchor α 2 s)) :=
  steffensen_solves_ae_mod_mcmullen 2 (by norm_num) x hx_ne α hα_ne_zero hα_pow hSp
    (mcmullen_p2_complex_unconditional x hx_ne α hα_ne_zero
      hα_ne_neg_one hα_ne_one hα_pow)

end FinalUnconditionalSolveC

end Pandrosion
