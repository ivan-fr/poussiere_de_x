/-
  Universitas Pandrosion — §48. **Explicit geometric form of the
  `p = 2` complex Julia circle.**

  Refines §39.2 by exposing the Julia circle
      `J = {z ∈ ℂ : ‖v_p2_C α z‖ = 1}`
  as an **explicit metric sphere** with closed-form centre and
  radius (in the non-degenerate Apollonius regime `‖μ²‖ ≠ 1`).

  Setting
      `m := ‖μ²‖²`    (μ² := ((1−α)/(1+α))²),
      `k := (m + 1) / (m − 1)`,
      `juliaCenter α  := k · α`,
      `juliaRadius α  := √((k² − 1) · ‖α‖²)`,
  we obtain
      `juliaSet α = Metric.sphere (juliaCenter α) (juliaRadius α)`.

  Contents.

    §48.1  `juliaCenter`, `juliaRadius` — closed-form definitions.

    §48.2  `julia_apollonius_bridge` — polynomial pivot lemma:
           the Apollonius condition in ℝ² coordinates is equivalent
           to `‖z − k·α‖² = (k² − 1)·‖α‖²`. Shared by the two
           inclusion proofs.

    §48.3  `juliaSet_eq_sphere` — full equality.
-/

import Pandrosion.Core.ComplexMcMullenP2Unconditional

namespace Pandrosion

open MeasureTheory Complex

/-! ============================================================
  §48.1  Closed-form centre and radius
============================================================ -/

section JuliaGeometryDefsC

/-- **Apollonius ratio** `k := (‖μ²‖² + 1) / (‖μ²‖² − 1)`. -/
noncomputable def juliaApollonius_k (α : ℂ) : ℝ :=
  (‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 + 1) /
  (‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 - 1)

/-- **Julia-circle centre** `juliaCenter α := k · α`. -/
noncomputable def juliaCenter (α : ℂ) : ℂ :=
  (juliaApollonius_k α : ℝ) * α

/-- **Julia-circle radius** `juliaRadius α := √((k² − 1) · ‖α‖²)`. -/
noncomputable def juliaRadius (α : ℂ) : ℝ :=
  Real.sqrt (((juliaApollonius_k α) ^ 2 - 1) * ‖α‖ ^ 2)

end JuliaGeometryDefsC

/-! ============================================================
  §48.2  Apollonius polynomial pivot
============================================================ -/

section ApollonPivotC

/-- Abbreviation for `‖w‖² = w.re² + w.im²`. -/
private theorem norm_sq_re_im (w : ℂ) : ‖w‖ ^ 2 = w.re ^ 2 + w.im ^ 2 := by
  rw [Complex.norm_eq_abs, Complex.sq_abs, Complex.normSq_apply]; ring

/-- **`k² − 1 ≥ 0` in the Apollonius regime.** -/
private theorem k_sq_minus_one_nn
    (α : ℂ) (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    0 ≤ (juliaApollonius_k α) ^ 2 - 1 := by
  unfold juliaApollonius_k
  set m : ℝ := ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 with hm_def
  have hm_ne_one : m ≠ 1 := by
    intro h
    apply hμ
    have h_norm_nn : 0 ≤ ‖((1 - α) / (1 + α)) ^ 2‖ := norm_nonneg _
    have hsq : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 = 1 ^ 2 := by
      rw [show (1:ℝ)^2 = 1 from by ring, ← hm_def, h]
    exact (sq_eq_sq h_norm_nn zero_le_one).mp hsq
  have hm_nn : 0 ≤ m := by rw [hm_def]; positivity
  have hm_minus_one_ne : m - 1 ≠ 0 := sub_ne_zero.mpr hm_ne_one
  have h_eq : ((m + 1) / (m - 1)) ^ 2 - 1 = (4 * m) / (m - 1) ^ 2 := by
    field_simp; ring
  rw [h_eq]
  exact div_nonneg (by linarith) (by positivity)

/-- **`k · (m − 1) = m + 1`.** -/
private theorem k_times_m_minus_one (α : ℂ)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    juliaApollonius_k α * (‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 - 1) =
      ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 + 1 := by
  unfold juliaApollonius_k
  set m : ℝ := ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 with hm_def
  have hm_ne_one : m ≠ 1 := by
    intro h
    apply hμ
    have h_norm_nn : 0 ≤ ‖((1 - α) / (1 + α)) ^ 2‖ := norm_nonneg _
    have hsq : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 = 1 ^ 2 := by
      rw [show (1:ℝ)^2 = 1 from by ring, ← hm_def, h]
    exact (sq_eq_sq h_norm_nn zero_le_one).mp hsq
  have hm_minus_one_ne : m - 1 ≠ 0 := sub_ne_zero.mpr hm_ne_one
  field_simp

/-- **Centre and radius in real coordinates.** -/
private theorem juliaCenter_re (α : ℂ) :
    (juliaCenter α).re = juliaApollonius_k α * α.re := by
  unfold juliaCenter
  simp [Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im]

private theorem juliaCenter_im (α : ℂ) :
    (juliaCenter α).im = juliaApollonius_k α * α.im := by
  unfold juliaCenter
  simp [Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im]

private theorem juliaRadius_sq (α : ℂ)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    juliaRadius α ^ 2 = ((juliaApollonius_k α) ^ 2 - 1) * ‖α‖ ^ 2 := by
  unfold juliaRadius
  exact Real.sq_sqrt (mul_nonneg (k_sq_minus_one_nn α hμ) (by positivity))

/-- **★ Apollonius polynomial equivalence.**

    Setting `m := ‖μ²‖²`, `k := juliaApollonius_k α`, the two
    polynomial-in-coordinates identities
      (A) `m · ‖z − α‖² = ‖z + α‖²`   (Apollonius form)
      (B) `‖z − k·α‖² = (k² − 1) · ‖α‖²`   (sphere form)
    are equivalent (for `‖μ²‖ ≠ 1`). -/
private theorem julia_apollonius_bridge (α z : ℂ)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2 ↔
    ‖z - juliaCenter α‖ ^ 2 = ((juliaApollonius_k α) ^ 2 - 1) * ‖α‖ ^ 2 := by
  set m : ℝ := ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 with hm_def
  set k : ℝ := juliaApollonius_k α with hk_def
  have hm_ne_one : m ≠ 1 := by
    intro h
    apply hμ
    have h_norm_nn : 0 ≤ ‖((1 - α) / (1 + α)) ^ 2‖ := norm_nonneg _
    have hsq : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 = 1 ^ 2 := by
      rw [show (1:ℝ)^2 = 1 from by ring, ← hm_def, h]
    exact (sq_eq_sq h_norm_nn zero_le_one).mp hsq
  have hm_minus_one_ne : m - 1 ≠ 0 := sub_ne_zero.mpr hm_ne_one
  have h_k_form : k * (m - 1) = m + 1 := by
    rw [hk_def, hm_def]
    exact k_times_m_minus_one α hμ
  -- Expand both sides in coordinates.
  have h_A_coord : m * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2 ↔
      m * ((z.re - α.re) ^ 2 + (z.im - α.im) ^ 2) =
        (z.re + α.re) ^ 2 + (z.im + α.im) ^ 2 := by
    rw [norm_sq_re_im (z - α), norm_sq_re_im (z + α)]
    simp [Complex.sub_re, Complex.sub_im, Complex.add_re, Complex.add_im]
  have h_B_coord : ‖z - juliaCenter α‖ ^ 2 = (k ^ 2 - 1) * ‖α‖ ^ 2 ↔
      (z.re - k * α.re) ^ 2 + (z.im - k * α.im) ^ 2 =
        (k ^ 2 - 1) * (α.re ^ 2 + α.im ^ 2) := by
    rw [norm_sq_re_im (z - juliaCenter α), norm_sq_re_im α]
    simp [Complex.sub_re, Complex.sub_im, juliaCenter_re, juliaCenter_im]
  rw [h_A_coord, h_B_coord]
  -- Direct linear-combination equivalence.
  constructor
  · intro h
    -- (A)_coord ⟹ (B)_coord via k·(m−1) = m+1.
    -- Multiplying (B)_coord_difference by (m−1) yields
    -- (m−1)·T₁ − 2k·(m−1)·T₂, which by h_k_form = (m−1)·T₁ − 2(m+1)·T₂
    -- = (A)_coord_difference. So linear combination works with
    -- coefficient 1 for h and coefficient −2·T₂ for h_k_form.
    have h_lin :
        (m - 1) *
          ((z.re - k * α.re) ^ 2 + (z.im - k * α.im) ^ 2 -
           (k ^ 2 - 1) * (α.re ^ 2 + α.im ^ 2)) = 0 := by
      linear_combination h - 2 * (z.re * α.re + z.im * α.im) * h_k_form
    have h_zero :
        (z.re - k * α.re) ^ 2 + (z.im - k * α.im) ^ 2 -
           (k ^ 2 - 1) * (α.re ^ 2 + α.im ^ 2) = 0 := by
      rcases mul_eq_zero.mp h_lin with h1 | h2
      · exact absurd h1 hm_minus_one_ne
      · exact h2
    linarith
  · intro h
    -- (B)_coord ⟹ (A)_coord: multiply (A)_coord_difference by 1,
    -- express via (m−1)·(B)_coord + 2T₂·h_k_form.
    -- (A_diff) = (m−1)·T₁ − 2(m+1)·T₂
    --         = (m−1)·[(B_diff) + 2k·T₂] − 2(m+1)·T₂  (since B_diff = T₁ − 2k·T₂)
    --         = (m−1)·(B_diff) + 2T₂·(k(m−1) − (m+1))
    --         = (m−1)·(B_diff) + 0 via h_k_form.
    linear_combination (m - 1) * h + 2 * (z.re * α.re + z.im * α.im) * h_k_form

end ApollonPivotC

/-! ============================================================
  §48.3  Full equality `juliaSet α = Metric.sphere …`
============================================================ -/

section JuliaSetEqSphereC

/-- **`−α` is not on the sphere** (for `α ≠ 0, ±1`, `‖μ²‖ ≠ 1`).
    If it were, `‖(k+1)α‖² = (k²−1)‖α‖²`, i.e., `(k+1)² = (k²−1)`
    (since `‖α‖² ≠ 0`), yielding `k = −1`, i.e., `m = 0`,
    i.e., `α = 1` — excluded. -/
private theorem neg_alpha_not_on_sphere
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    -α ∉ Metric.sphere (juliaCenter α) (juliaRadius α) := by
  intro h
  have h_dist : dist (-α) (juliaCenter α) = juliaRadius α := h
  set m : ℝ := ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 with hm_def
  set k : ℝ := juliaApollonius_k α with hk_def
  have hm_ne_one : m ≠ 1 := by
    intro h_eq
    apply hμ
    have h_norm_nn : 0 ≤ ‖((1 - α) / (1 + α)) ^ 2‖ := norm_nonneg _
    have hsq : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 = 1 ^ 2 := by
      rw [show (1:ℝ)^2 = 1 from by ring, ← hm_def, h_eq]
    exact (sq_eq_sq h_norm_nn zero_le_one).mp hsq
  have hm_minus_one_ne : m - 1 ≠ 0 := sub_ne_zero.mpr hm_ne_one
  have h_r_sq : juliaRadius α ^ 2 = (k ^ 2 - 1) * ‖α‖ ^ 2 := by
    rw [hk_def]; exact juliaRadius_sq α hμ
  have h_dist_sq : ‖-α - juliaCenter α‖ ^ 2 = (k ^ 2 - 1) * ‖α‖ ^ 2 := by
    rw [dist_eq_norm] at h_dist
    have : ‖-α - juliaCenter α‖ ^ 2 = juliaRadius α ^ 2 := by rw [h_dist]
    rw [this, h_r_sq]
  -- ‖-α - juliaCenter α‖² = (k+1)² · ‖α‖².
  have h_lhs : ‖-α - juliaCenter α‖ ^ 2 = (k + 1) ^ 2 * ‖α‖ ^ 2 := by
    rw [norm_sq_re_im (-α - juliaCenter α), norm_sq_re_im α]
    simp [Complex.sub_re, Complex.sub_im, Complex.neg_re, Complex.neg_im,
          juliaCenter_re, juliaCenter_im]
    ring
  rw [h_lhs] at h_dist_sq
  -- (k+1)²·‖α‖² = (k²-1)·‖α‖² ⟹ (2k+2)·‖α‖² = 0.
  have h_α_sq_pos : 0 < ‖α‖ ^ 2 := by
    have : 0 < ‖α‖ := norm_pos_iff.mpr hα_ne_zero
    positivity
  have h_factor : (2 * k + 2) * ‖α‖ ^ 2 = 0 := by linarith [h_dist_sq]
  have h_k_val : k = -1 := by
    rcases mul_eq_zero.mp h_factor with h1 | h2
    · linarith
    · linarith
  rw [hk_def] at h_k_val
  unfold juliaApollonius_k at h_k_val
  rw [← hm_def] at h_k_val
  have h_eq : (m + 1) = -1 * (m - 1) := by
    have h_mul : (m + 1) / (m - 1) * (m - 1) = -1 * (m - 1) := by
      rw [h_k_val]
    rw [div_mul_cancel₀ _ hm_minus_one_ne] at h_mul
    exact h_mul
  have h_m_zero : m = 0 := by linarith
  have h_μ_sq_zero : ‖((1 - α) / (1 + α)) ^ 2‖ = 0 := by
    have : ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 = 0 := by rw [← hm_def]; exact h_m_zero
    exact pow_eq_zero_iff (n := 2) (by norm_num : (2:ℕ) ≠ 0) |>.mp this
  have h_μ_zero : ((1 - α) / (1 + α)) ^ 2 = 0 := norm_eq_zero.mp h_μ_sq_zero
  have h_div_zero : (1 - α) / (1 + α) = 0 :=
    pow_eq_zero_iff (n := 2) (by norm_num : (2:ℕ) ≠ 0) |>.mp h_μ_zero
  rcases div_eq_zero_iff.mp h_div_zero with h_num | h_denom
  · exact hα_ne_one h_num
  · exact hα_ne_neg_one h_denom

/-- **★★★★ `juliaSet α = Metric.sphere (juliaCenter α) (juliaRadius α)`.**

    The Julia circle on ℂ is **exactly** the Apollonius metric
    sphere with the closed-form centre `k·α` and radius
    `√((k² − 1)‖α‖²)` for `k = (‖μ²‖² + 1)/(‖μ²‖² − 1)`. -/
theorem juliaSet_eq_sphere
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hμ : ‖((1 - α) / (1 + α)) ^ 2‖ ≠ 1) :
    {z : ℂ | ‖v_p2_C α z‖ = 1} =
      Metric.sphere (juliaCenter α) (juliaRadius α) := by
  ext z
  constructor
  · -- Forward: ‖v(z)‖ = 1 ⟹ z ∈ sphere.
    intro hz
    simp only [Set.mem_setOf_eq] at hz
    have hz_α : z + α ≠ 0 := fun h_eq => julia_excludes_neg_alpha α z h_eq hz
    rw [julia_squared_apollonius α z hz_α] at hz
    -- hz: ‖μ²‖² · ‖z-α‖² = ‖z+α‖². Apply bridge.
    have h_sphere := (julia_apollonius_bridge α z hμ).mp hz
    -- h_sphere: ‖z - juliaCenter α‖² = (k² - 1) · ‖α‖².
    show z ∈ Metric.sphere (juliaCenter α) (juliaRadius α)
    rw [Metric.mem_sphere, dist_eq_norm]
    have h_zc_nn : 0 ≤ ‖z - juliaCenter α‖ := norm_nonneg _
    have h_r_nn : 0 ≤ juliaRadius α := Real.sqrt_nonneg _
    have h_r_sq := juliaRadius_sq α hμ
    have h_sq_match : ‖z - juliaCenter α‖ ^ 2 = juliaRadius α ^ 2 := by
      rw [h_sphere, h_r_sq]
    exact (sq_eq_sq h_zc_nn h_r_nn).mp h_sq_match
  · -- Reverse: z ∈ sphere ⟹ ‖v(z)‖ = 1.
    intro hz
    simp only [Set.mem_setOf_eq]
    have h_dist : dist z (juliaCenter α) = juliaRadius α := hz
    have h_r_sq := juliaRadius_sq α hμ
    have h_zc_sq : ‖z - juliaCenter α‖ ^ 2 =
        ((juliaApollonius_k α) ^ 2 - 1) * ‖α‖ ^ 2 := by
      rw [dist_eq_norm] at h_dist
      have h_eq_sq : ‖z - juliaCenter α‖ ^ 2 = juliaRadius α ^ 2 := by rw [h_dist]
      rw [h_eq_sq, h_r_sq]
    have h_apollonius :
        ‖((1 - α) / (1 + α)) ^ 2‖ ^ 2 * ‖z - α‖ ^ 2 = ‖z + α‖ ^ 2 :=
      (julia_apollonius_bridge α z hμ).mpr h_zc_sq
    -- z + α ≠ 0 (else z = -α, contradicting neg_alpha_not_on_sphere).
    have h_zα_ne : z + α ≠ 0 := by
      intro h_eq
      have hz_eq : z = -α := by linear_combination h_eq
      apply neg_alpha_not_on_sphere α hα_ne_zero hα_ne_neg_one hα_ne_one hμ
      rw [← hz_eq]; exact hz
    have hzα_sq_pos : 0 < ‖z + α‖ ^ 2 := by
      have : 0 < ‖z + α‖ := norm_pos_iff.mpr h_zα_ne
      positivity
    have h_v_sq : ‖v_p2_C α z‖ ^ 2 = 1 := by
      rw [v_p2_C_norm_sq_eq α z]
      rw [div_eq_one_iff_eq (ne_of_gt hzα_sq_pos)]
      exact h_apollonius
    have h_v_nn : 0 ≤ ‖v_p2_C α z‖ := norm_nonneg _
    rw [show (1 : ℝ) = 1 ^ 2 from by ring] at h_v_sq
    exact (sq_eq_sq h_v_nn zero_le_one).mp h_v_sq

end JuliaSetEqSphereC

end Pandrosion
