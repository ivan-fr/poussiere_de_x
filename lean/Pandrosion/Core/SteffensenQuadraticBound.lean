/-
  Universitas Pandrosion — §31. **Unconditional quadratic bound for
  `steffensen_step_C`.**

  Load-bearing theorem: at any Pandrosion fixed point α (with p ≥ 1,
  α ≠ 0, S_p(α) ≠ 0, α^p = 1/x), there exist K > 0 and r > 0 such
  that for every z ∈ B(α, r),
      `‖steffensen_step_C x p z − α‖ ≤ K · ‖z − α‖²`.

  This discharges the `h_quadratic` hypothesis of §27
  `steffensen_local_attraction` and §28
  `steffensen_dynamical_convergence_ae`, upgrading both to fully
  unconditional statements (§31.2-§31.3).

  Proof outline. Denote `h := pandrosion_h_C x p` and
  `lam := lambdaClosedC p α`. Writing `R(z) := h(z) − α − lam(z − α)`,
  §30 gives `|R(z)| ≤ M·|z − α|²` on `B(α, r₀)`. From this:

    (i)   `h(z) − z = (lam − 1)(z − α) + R(z)`;

    (ii)  `h(h(z)) − α = lam²(z − α) + lamR(z) + R(h(z))`;

    (iii) Subtracting, the Steffensen denominator has the form
              `D(z) = (lam − 1)²(z − α) + (lam − 2)R(z) + R(h(z))`;

    (iv)  Since `lam ≠ 1` (§29.1), for `|z − α|` small enough,
              `|D(z)| ≥ (|lam − 1|²/2)·|z − α|`;

    (v)   `(h(z) − z)² = (lam − 1)²(z − α)² + 2(lam − 1)(z − α)R(z) + R(z)²`;

    (vi)  Exact rearrangement:
              `S(z) − α = [(z − α)·D(z) − (h(z) − z)²] / D(z)`,
          whose numerator is `O(|z − α|³)`, hence
              `|S(z) − α| ≤ K·|z − α|²`.
-/

import Pandrosion.Core.QuadraticBound
import Pandrosion.Core.LocalAttraction
import Pandrosion.Core.DynamicalConvergence

namespace Pandrosion

open Complex Filter Topology MeasureTheory

/-! ============================================================
  §31.1  Main theorem — `steffensen_step_C_quadratic_bound`
============================================================ -/

/-- **★★★ Unconditional quadratic bound for `steffensen_step_C`.**

  At any Pandrosion fixed point α with p ≥ 1, α ≠ 0, S_p(α) ≠ 0,
  α^p = 1/x, there exist constants K > 0 and r > 0 such that for
  every z ∈ B(α, r),
      `‖steffensen_step_C x p z − α‖ ≤ K · ‖z − α‖²`.

  This is the super-attracting content of Steffensen's method on
  the Pandrosion iterator h, proved unconditionally from §30's
  Pandrosion Taylor bound and §29.1's `lam − 1 ≠ 0`. -/
theorem steffensen_step_C_quadratic_bound
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    ∃ K > 0, ∃ r > 0, ∀ z : ℂ, ‖z - α‖ < r →
      ‖steffensen_step_C x p z - α‖ ≤ K * ‖z - α‖ ^ 2 := by
  -- Pandrosion §30 Taylor bound for h.
  obtain ⟨M, hM_pos, r₀, hr₀_pos, h_taylor⟩ :=
    pandrosion_h_C_quadratic_bound x hx p hp α hSp hfp
  -- Fixed point value of h.
  have h_h_fp : pandrosion_h_C x p α = α :=
    pandrosion_h_C_fixed_point_reverse x hx p α hSp hfp
  -- Non-vanishing of lam − 1 (§29.1).
  have h_lam_ne : lambdaClosedC p α ≠ 1 :=
    lambdaClosedC_ne_one_at_fp hp hα hSp
  set lam : ℂ := lambdaClosedC p α with hlam_def
  have h_lam_sub_ne : lam - 1 ≠ 0 := sub_ne_zero.mpr h_lam_ne
  have h_lam_abs_pos : 0 < ‖lam - 1‖ := norm_pos_iff.mpr h_lam_sub_ne
  have h_lam_sq_pos : 0 < ‖lam - 1‖ ^ 2 := pow_pos h_lam_abs_pos 2
  -- Abbreviations.
  set A : ℝ := ‖lam‖ with hA_def
  have hA_nn : 0 ≤ A := norm_nonneg _
  -- Constants B, C, Mq aggregating the coefficient budgets.
  -- B := A + 1 bounds |lam| + M·r₁ when M·r₁ ≤ 1.
  -- Numerator residue budget  M_n := 2·‖lam − 1‖·M + M²
  -- Denominator residue budget M_d := ‖lam − 2‖·M + M·(A+1)²
  set Md : ℝ := ‖lam - 2‖ * M + M * (A + 1) ^ 2 with hMd_def
  set Mn : ℝ := 2 * ‖lam - 1‖ * M + M ^ 2 with hMn_def
  have hMd_nn : 0 ≤ Md := by
    rw [hMd_def]; positivity
  have hMn_nn : 0 ≤ Mn := by
    rw [hMn_def]; positivity
  -- Choose r₁ small enough for all the estimates to work.
  -- Constraints:
  --   (a) r₁ ≤ r₀               (so §30 applies to z)
  --   (b) M · r₁ ≤ 1            (controls residue amplification)
  --   (c) (A + 1) · r₁ ≤ r₀     (so h(z) ∈ B(α, r₀), enabling §30 on h(z))
  --   (d) Md · r₁ ≤ ‖lam − 1‖²/2 (lower bound on Steffensen denom)
  --   (e) r₁ ≤ 1                (so |z − α|³ ≤ |z − α|² on the ball)
  have hA2_pos : 0 < A + 2 := by linarith
  have hM1_pos : 0 < M + 1 := by linarith
  have hMd1_pos : 0 < Md + 1 := by linarith
  have h_half_pos : 0 < ‖lam - 1‖ ^ 2 / (2 * (Md + 1)) := by positivity
  have h_inv_M1_pos : 0 < 1 / (M + 1) := by positivity
  have h_r₀_A2_pos : 0 < r₀ / (A + 2) := by positivity
  -- Build r₁ from sequential mins so each constraint is a single min_le.
  set r_a : ℝ := min r₀ (1 / (M + 1)) with hra_def
  set r_b : ℝ := min r_a (r₀ / (A + 2)) with hrb_def
  set r_c : ℝ := min r_b (‖lam - 1‖ ^ 2 / (2 * (Md + 1))) with hrc_def
  set r₁ : ℝ := min r_c 1 with hr₁_def
  have hra_pos : 0 < r_a := lt_min hr₀_pos h_inv_M1_pos
  have hrb_pos : 0 < r_b := lt_min hra_pos h_r₀_A2_pos
  have hrc_pos : 0 < r_c := lt_min hrb_pos h_half_pos
  have hr₁_pos : 0 < r₁ := lt_min hrc_pos one_pos
  have hr₁_le_rc : r₁ ≤ r_c := min_le_left _ _
  have hrc_le_rb : r_c ≤ r_b := min_le_left _ _
  have hrb_le_ra : r_b ≤ r_a := min_le_left _ _
  have hra_le_r₀ : r_a ≤ r₀ := min_le_left _ _
  have hra_le_invM1 : r_a ≤ 1 / (M + 1) := min_le_right _ _
  have hrb_le_r₀A2 : r_b ≤ r₀ / (A + 2) := min_le_right _ _
  have hrc_le_half : r_c ≤ ‖lam - 1‖ ^ 2 / (2 * (Md + 1)) := min_le_right _ _
  have hr₁_le_1 : r₁ ≤ 1 := min_le_right _ _
  have hr₁_le_r₀ : r₁ ≤ r₀ :=
    (hr₁_le_rc.trans hrc_le_rb).trans (hrb_le_ra.trans hra_le_r₀)
  have hr₁_le_invM1 : r₁ ≤ 1 / (M + 1) :=
    (hr₁_le_rc.trans hrc_le_rb).trans (hrb_le_ra.trans hra_le_invM1)
  have hr₁_le_r₀A2 : r₁ ≤ r₀ / (A + 2) :=
    (hr₁_le_rc.trans hrc_le_rb).trans hrb_le_r₀A2
  have hr₁_le_half : r₁ ≤ ‖lam - 1‖ ^ 2 / (2 * (Md + 1)) :=
    hr₁_le_rc.trans hrc_le_half
  have hr₁_M_le_1 : M * r₁ ≤ 1 := by
    calc M * r₁ ≤ M * (1 / (M + 1)) :=
          mul_le_mul_of_nonneg_left hr₁_le_invM1 (le_of_lt hM_pos)
      _ ≤ (M + 1) * (1 / (M + 1)) :=
          mul_le_mul_of_nonneg_right (by linarith) (le_of_lt h_inv_M1_pos)
      _ = 1 := by field_simp
  have hr₁_A_le_r₀ : (A + 1) * r₁ ≤ r₀ := by
    calc (A + 1) * r₁ ≤ (A + 1) * (r₀ / (A + 2)) :=
          mul_le_mul_of_nonneg_left hr₁_le_r₀A2 (by linarith)
      _ ≤ (A + 2) * (r₀ / (A + 2)) :=
          mul_le_mul_of_nonneg_right (by linarith) (le_of_lt h_r₀_A2_pos)
      _ = r₀ := by field_simp; ring
  have hr₁_Md_le : Md * r₁ ≤ ‖lam - 1‖ ^ 2 / 2 := by
    calc Md * r₁ ≤ Md * (‖lam - 1‖ ^ 2 / (2 * (Md + 1))) :=
          mul_le_mul_of_nonneg_left hr₁_le_half hMd_nn
      _ ≤ (Md + 1) * (‖lam - 1‖ ^ 2 / (2 * (Md + 1))) :=
          mul_le_mul_of_nonneg_right (by linarith) (le_of_lt h_half_pos)
      _ = ‖lam - 1‖ ^ 2 / 2 := by field_simp; ring
  -- Final constants.
  set K : ℝ := 2 * (Md + Mn) / ‖lam - 1‖ ^ 2 with hK_def
  have hK_pos : 0 < K := by
    rw [hK_def]
    apply div_pos
    · have : 0 < Md + Mn := by positivity
      linarith
    · exact h_lam_sq_pos
  refine ⟨K, hK_pos, r₁, hr₁_pos, ?_⟩
  intro z hz
  -- Case split on z = α vs z ≠ α.
  by_cases hz_eq : z = α
  · -- z = α: everything is α, zero on both sides.
    have h_step_α : steffensen_step_C x p α = α := by
      unfold steffensen_step_C steffensen_denom_C
      rw [h_h_fp]
      simp
    rw [hz_eq, h_step_α]
    simp
  -- Now z ≠ α, so ‖z − α‖ > 0.
  have hz_ne : z - α ≠ 0 := sub_ne_zero.mpr hz_eq
  have hz_norm_pos : 0 < ‖z - α‖ := norm_pos_iff.mpr hz_ne
  -- §30 bound on h(z).
  have hz_r₀ : ‖z - α‖ < r₀ := lt_of_lt_of_le hz hr₁_le_r₀
  have h_Rz_bd : ‖pandrosion_h_C x p z - α - lam * (z - α)‖ ≤ M * ‖z - α‖ ^ 2 :=
    h_taylor z hz_r₀
  -- |h(z) − α| ≤ (A + M·‖z − α‖)·‖z − α‖ ≤ (A + 1)·‖z − α‖.
  have h_h_norm : ‖pandrosion_h_C x p z - α‖ ≤ (A + 1) * ‖z - α‖ := by
    have hz_triangle :
        ‖pandrosion_h_C x p z - α‖
          = ‖(pandrosion_h_C x p z - α - lam * (z - α)) + lam * (z - α)‖ := by
      congr 1; ring
    rw [hz_triangle]
    calc ‖(pandrosion_h_C x p z - α - lam * (z - α)) + lam * (z - α)‖
        ≤ ‖pandrosion_h_C x p z - α - lam * (z - α)‖ + ‖lam * (z - α)‖ :=
          norm_add_le _ _
      _ ≤ M * ‖z - α‖ ^ 2 + A * ‖z - α‖ := by
          apply add_le_add h_Rz_bd
          rw [norm_mul]
      _ = (A + M * ‖z - α‖) * ‖z - α‖ := by ring
      _ ≤ (A + 1) * ‖z - α‖ := by
          apply mul_le_mul_of_nonneg_right _ (norm_nonneg _)
          have h1 : M * ‖z - α‖ ≤ M * r₁ :=
            mul_le_mul_of_nonneg_left (le_of_lt hz) (le_of_lt hM_pos)
          linarith
  -- h(z) is within r₀ of α, so §30 applies to h(z) too.
  have h_hz_r₀ : ‖pandrosion_h_C x p z - α‖ < r₀ := by
    calc ‖pandrosion_h_C x p z - α‖ ≤ (A + 1) * ‖z - α‖ := h_h_norm
      _ < (A + 1) * r₁ := by
          apply mul_lt_mul_of_pos_left hz
          linarith
      _ ≤ r₀ := hr₁_A_le_r₀
  -- §30 applied to h(z): Taylor bound for h(h(z)) around α.
  have h_Rhz_bd :
      ‖pandrosion_h_C x p (pandrosion_h_C x p z) - α - lam * (pandrosion_h_C x p z - α)‖
        ≤ M * ‖pandrosion_h_C x p z - α‖ ^ 2 := by
    have := h_taylor (pandrosion_h_C x p z) h_hz_r₀
    convert this using 2
  -- Bound ‖pandrosion_h_C x p z - α‖² ≤ (A + 1)² · ‖z − α‖².
  have h_hz_sq : ‖pandrosion_h_C x p z - α‖ ^ 2 ≤ (A + 1) ^ 2 * ‖z - α‖ ^ 2 := by
    have h_A1_nn : 0 ≤ A + 1 := by linarith
    have h_norm_nn : 0 ≤ ‖pandrosion_h_C x p z - α‖ := norm_nonneg _
    calc ‖pandrosion_h_C x p z - α‖ ^ 2
        ≤ ((A + 1) * ‖z - α‖) ^ 2 := by
          apply sq_le_sq' _ h_h_norm
          have := mul_nonneg h_A1_nn (norm_nonneg (z - α))
          linarith
      _ = (A + 1) ^ 2 * ‖z - α‖ ^ 2 := by ring
  -- Steffensen denominator: D(z) = (lam − 1)²(z − α) + Rd(z),
  -- where Rd(z) := (lam − 2)·R(z) + R(h(z)) (with R := h − id at α via §30).
  -- Setting R_z := pandrosion_h_C x p z − α − lam(z − α) and
  --         R_hz := h(h(z)) − α − lam(h(z) − α),
  -- we have:
  --   steffensen_denom_C = (lam − 1)²(z − α) + (lam − 2)·R_z + R_hz.
  set Rz : ℂ := pandrosion_h_C x p z - α - lam * (z - α) with hRz_def
  set Rhz : ℂ :=
    pandrosion_h_C x p (pandrosion_h_C x p z) - α - lam * (pandrosion_h_C x p z - α)
    with hRhz_def
  have hRz_bd : ‖Rz‖ ≤ M * ‖z - α‖ ^ 2 := h_Rz_bd
  have hRhz_bd : ‖Rhz‖ ≤ M * (A + 1) ^ 2 * ‖z - α‖ ^ 2 := by
    calc ‖Rhz‖ ≤ M * ‖pandrosion_h_C x p z - α‖ ^ 2 := h_Rhz_bd
      _ ≤ M * ((A + 1) ^ 2 * ‖z - α‖ ^ 2) :=
          mul_le_mul_of_nonneg_left h_hz_sq (le_of_lt hM_pos)
      _ = M * (A + 1) ^ 2 * ‖z - α‖ ^ 2 := by ring
  -- Denominator algebraic identity:
  --   D(z) = (lam − 1)²·(z − α) + (lam − 2)·R_z + R_hz.
  have h_denom_id :
      steffensen_denom_C x p z = (lam - 1) ^ 2 * (z - α) + ((lam - 2) * Rz + Rhz) := by
    unfold steffensen_denom_C
    rw [hRz_def, hRhz_def]
    ring
  -- Denominator residue bound.
  have h_denom_residue_bd :
      ‖(lam - 2) * Rz + Rhz‖ ≤ Md * ‖z - α‖ ^ 2 := by
    calc ‖(lam - 2) * Rz + Rhz‖
        ≤ ‖(lam - 2) * Rz‖ + ‖Rhz‖ := norm_add_le _ _
      _ = ‖lam - 2‖ * ‖Rz‖ + ‖Rhz‖ := by rw [norm_mul]
      _ ≤ ‖lam - 2‖ * (M * ‖z - α‖ ^ 2) + M * (A + 1) ^ 2 * ‖z - α‖ ^ 2 := by
          apply add_le_add _ hRhz_bd
          apply mul_le_mul_of_nonneg_left hRz_bd (norm_nonneg _)
      _ = (‖lam - 2‖ * M + M * (A + 1) ^ 2) * ‖z - α‖ ^ 2 := by ring
      _ = Md * ‖z - α‖ ^ 2 := by rw [hMd_def]
  -- Denominator lower bound: ‖D(z)‖ ≥ (‖lam − 1‖²/2)·‖z − α‖.
  have h_denom_lb :
      (‖lam - 1‖ ^ 2 / 2) * ‖z - α‖ ≤ ‖steffensen_denom_C x p z‖ := by
    rw [h_denom_id]
    have h_tri :
        ‖(lam - 1) ^ 2 * (z - α)‖
          ≤ ‖(lam - 1) ^ 2 * (z - α) + ((lam - 2) * Rz + Rhz)‖ + ‖(lam - 2) * Rz + Rhz‖ := by
      have := norm_sub_le
        ((lam - 1) ^ 2 * (z - α) + ((lam - 2) * Rz + Rhz))
        ((lam - 2) * Rz + Rhz)
      have h_eq :
          (lam - 1) ^ 2 * (z - α) + ((lam - 2) * Rz + Rhz) - ((lam - 2) * Rz + Rhz)
            = (lam - 1) ^ 2 * (z - α) := by ring
      rw [h_eq] at this
      linarith
    have h_main_norm :
        ‖(lam - 1) ^ 2 * (z - α)‖ = ‖lam - 1‖ ^ 2 * ‖z - α‖ := by
      rw [norm_mul, norm_pow]
    rw [h_main_norm] at h_tri
    -- So ‖lam − 1‖²·‖z − α‖ ≤ ‖D‖ + ‖residue‖, i.e.,
    --   ‖D‖ ≥ ‖lam − 1‖²·‖z − α‖ − ‖residue‖
    --       ≥ ‖lam − 1‖²·‖z − α‖ − Md·‖z − α‖²
    --       = (‖lam − 1‖² − Md·‖z − α‖)·‖z − α‖
    --       ≥ (‖lam − 1‖²/2)·‖z − α‖
    -- using Md·‖z − α‖ ≤ Md·r₁ ≤ ‖lam − 1‖²/2.
    have h_Md_le : Md * ‖z - α‖ ≤ ‖lam - 1‖ ^ 2 / 2 := by
      calc Md * ‖z - α‖ ≤ Md * r₁ :=
            mul_le_mul_of_nonneg_left (le_of_lt hz) hMd_nn
        _ ≤ ‖lam - 1‖ ^ 2 / 2 := hr₁_Md_le
    have h_residue_bd' :
        ‖(lam - 2) * Rz + Rhz‖ ≤ Md * ‖z - α‖ * ‖z - α‖ := by
      have : ‖(lam - 2) * Rz + Rhz‖ ≤ Md * ‖z - α‖ ^ 2 := h_denom_residue_bd
      have h_eq : ‖z - α‖ ^ 2 = ‖z - α‖ * ‖z - α‖ := by ring
      rw [h_eq] at this
      linarith
    calc (‖lam - 1‖ ^ 2 / 2) * ‖z - α‖
        = (‖lam - 1‖ ^ 2 - ‖lam - 1‖ ^ 2 / 2) * ‖z - α‖ := by ring
      _ ≤ (‖lam - 1‖ ^ 2 - Md * ‖z - α‖) * ‖z - α‖ := by
          apply mul_le_mul_of_nonneg_right _ (norm_nonneg _)
          linarith
      _ = ‖lam - 1‖ ^ 2 * ‖z - α‖ - Md * ‖z - α‖ * ‖z - α‖ := by ring
      _ ≤ ‖(lam - 1) ^ 2 * (z - α) + ((lam - 2) * Rz + Rhz)‖ := by
          linarith [h_tri, h_residue_bd']
  -- So the denominator is non-zero when z ≠ α.
  have h_denom_ne :
      steffensen_denom_C x p z ≠ 0 := by
    intro h_zero
    rw [h_zero, norm_zero] at h_denom_lb
    have : (‖lam - 1‖ ^ 2 / 2) * ‖z - α‖ > 0 := by
      apply mul_pos _ hz_norm_pos
      positivity
    linarith
  -- Unfold `steffensen_step_C` using the non-zero denominator.
  unfold steffensen_step_C
  rw [if_neg h_denom_ne]
  -- Goal: ‖z − (h(z) − z)² / D(z) − α‖ ≤ K · ‖z − α‖²
  -- Rearrange: z − (h(z) − z)²/D(z) − α = [(z − α)·D(z) − (h(z) − z)²]/D(z)
  have h_numer_id :
      (pandrosion_h_C x p z - z) ^ 2
        = (lam - 1) ^ 2 * (z - α) ^ 2 + (2 * (lam - 1) * (z - α) * Rz + Rz ^ 2) := by
    have h_g : pandrosion_h_C x p z - z = (lam - 1) * (z - α) + Rz := by
      rw [hRz_def]; ring
    rw [h_g]
    ring
  -- Numerator residue bound.
  have h_Rz_sq_bd : ‖Rz‖ ^ 2 ≤ M ^ 2 * ‖z - α‖ ^ 4 := by
    have h_Rz_nn : 0 ≤ ‖Rz‖ := norm_nonneg _
    have h_rhs_nn : 0 ≤ M * ‖z - α‖ ^ 2 := by positivity
    calc ‖Rz‖ ^ 2 ≤ (M * ‖z - α‖ ^ 2) ^ 2 := by
          apply sq_le_sq' _ hRz_bd
          linarith
      _ = M ^ 2 * ‖z - α‖ ^ 4 := by ring
  have h_numer_residue_bd :
      ‖2 * (lam - 1) * (z - α) * Rz + Rz ^ 2‖ ≤ Mn * ‖z - α‖ ^ 3 := by
    calc ‖2 * (lam - 1) * (z - α) * Rz + Rz ^ 2‖
        ≤ ‖2 * (lam - 1) * (z - α) * Rz‖ + ‖Rz ^ 2‖ := norm_add_le _ _
      _ = 2 * ‖lam - 1‖ * ‖z - α‖ * ‖Rz‖ + ‖Rz‖ ^ 2 := by
          rw [show (2 : ℂ) * (lam - 1) * (z - α) * Rz
                = ((2 : ℂ) * (lam - 1) * (z - α)) * Rz by ring]
          rw [norm_mul, norm_mul, norm_mul, sq, norm_mul]
          have h2 : ‖(2 : ℂ)‖ = 2 := by norm_num
          rw [h2]
          ring
      _ ≤ 2 * ‖lam - 1‖ * ‖z - α‖ * (M * ‖z - α‖ ^ 2) + M ^ 2 * ‖z - α‖ ^ 4 := by
          apply add_le_add _ h_Rz_sq_bd
          apply mul_le_mul_of_nonneg_left hRz_bd
          positivity
      _ = (2 * ‖lam - 1‖ * M) * ‖z - α‖ ^ 3 + M ^ 2 * ‖z - α‖ ^ 4 := by ring
      _ ≤ (2 * ‖lam - 1‖ * M) * ‖z - α‖ ^ 3 + M ^ 2 * ‖z - α‖ ^ 3 := by
          have h_sq_bd : ‖z - α‖ ^ 4 ≤ ‖z - α‖ ^ 3 := by
            have h4 : ‖z - α‖ ^ 4 = ‖z - α‖ ^ 3 * ‖z - α‖ := by ring
            have h3_nn : 0 ≤ ‖z - α‖ ^ 3 := by positivity
            rw [h4]
            calc ‖z - α‖ ^ 3 * ‖z - α‖ ≤ ‖z - α‖ ^ 3 * 1 := by
                  apply mul_le_mul_of_nonneg_left _ h3_nn
                  have : ‖z - α‖ ≤ r₁ := le_of_lt hz
                  linarith
              _ = ‖z - α‖ ^ 3 := by ring
          have h_M2_nn : 0 ≤ M ^ 2 := by positivity
          have : M ^ 2 * ‖z - α‖ ^ 4 ≤ M ^ 2 * ‖z - α‖ ^ 3 :=
            mul_le_mul_of_nonneg_left h_sq_bd h_M2_nn
          linarith
      _ = (2 * ‖lam - 1‖ * M + M ^ 2) * ‖z - α‖ ^ 3 := by ring
      _ = Mn * ‖z - α‖ ^ 3 := by rw [hMn_def]
  -- Final algebraic identity:
  --   S(z) − α = z − (h(z) − z)²/D(z) − α
  --           = ((z − α)·D(z) − (h(z) − z)²) / D(z)
  --           = [(z − α)·((lam−1)²(z−α) + Rdenom) − ((lam−1)²(z−α)² + Rnum)] / D(z)
  --           = [(z−α)·Rdenom − Rnum] / D(z)
  have h_step_id :
      z - (pandrosion_h_C x p z - z) ^ 2 / steffensen_denom_C x p z - α
        = ((z - α) * ((lam - 2) * Rz + Rhz) - (2 * (lam - 1) * (z - α) * Rz + Rz ^ 2))
            / steffensen_denom_C x p z := by
    rw [eq_div_iff h_denom_ne]
    have h_numer_exp := h_numer_id
    have h_denom_exp := h_denom_id
    have h_lhs_expand :
        (z - (pandrosion_h_C x p z - z) ^ 2 / steffensen_denom_C x p z - α)
            * steffensen_denom_C x p z
          = (z - α) * steffensen_denom_C x p z - (pandrosion_h_C x p z - z) ^ 2 := by
      field_simp
      ring
    rw [h_lhs_expand, h_denom_exp, h_numer_exp]
    ring
  rw [h_step_id]
  -- Norm bound on the numerator.
  have h_step_numer_norm :
      ‖(z - α) * ((lam - 2) * Rz + Rhz) - (2 * (lam - 1) * (z - α) * Rz + Rz ^ 2)‖
        ≤ (Md + Mn) * ‖z - α‖ ^ 3 := by
    calc ‖(z - α) * ((lam - 2) * Rz + Rhz) - (2 * (lam - 1) * (z - α) * Rz + Rz ^ 2)‖
        ≤ ‖(z - α) * ((lam - 2) * Rz + Rhz)‖
            + ‖2 * (lam - 1) * (z - α) * Rz + Rz ^ 2‖ := by
          rw [sub_eq_add_neg]
          calc ‖_ + _‖ ≤ _ + ‖_‖ := norm_add_le _ _
            _ = _ := by rw [norm_neg]
      _ = ‖z - α‖ * ‖(lam - 2) * Rz + Rhz‖
            + ‖2 * (lam - 1) * (z - α) * Rz + Rz ^ 2‖ := by
          rw [norm_mul]
      _ ≤ ‖z - α‖ * (Md * ‖z - α‖ ^ 2) + Mn * ‖z - α‖ ^ 3 := by
          apply add_le_add _ h_numer_residue_bd
          exact mul_le_mul_of_nonneg_left h_denom_residue_bd (norm_nonneg _)
      _ = (Md + Mn) * ‖z - α‖ ^ 3 := by ring
  -- Division norm bound.
  rw [norm_div]
  have h_denom_norm_pos : 0 < ‖steffensen_denom_C x p z‖ :=
    norm_pos_iff.mpr h_denom_ne
  rw [div_le_iff h_denom_norm_pos]
  -- Goal: ‖numer‖ ≤ K · ‖z − α‖² · ‖D(z)‖
  -- We have: ‖numer‖ ≤ (Md + Mn)·‖z − α‖³ and
  --          ‖D(z)‖ ≥ (‖lam − 1‖²/2)·‖z − α‖
  -- So K · ‖z − α‖² · ‖D(z)‖ ≥ K · ‖z − α‖² · (‖lam − 1‖²/2) · ‖z − α‖
  --                         = K · (‖lam − 1‖²/2) · ‖z − α‖³
  --                         = (Md + Mn) · ‖z − α‖³  [by choice of K]
  --                         ≥ ‖numer‖.
  have hK_eq : K * (‖lam - 1‖ ^ 2 / 2) = Md + Mn := by
    have h_ne : ‖lam - 1‖ ^ 2 ≠ 0 := ne_of_gt h_lam_sq_pos
    have h_ne' : (Complex.abs (lam - 1)) ^ 2 ≠ 0 := by
      rw [← Complex.norm_eq_abs]; exact h_ne
    rw [hK_def]
    field_simp
    ring
  calc ‖(z - α) * ((lam - 2) * Rz + Rhz) - (2 * (lam - 1) * (z - α) * Rz + Rz ^ 2)‖
      ≤ (Md + Mn) * ‖z - α‖ ^ 3 := h_step_numer_norm
    _ = K * (‖lam - 1‖ ^ 2 / 2) * ‖z - α‖ ^ 3 := by rw [hK_eq]
    _ = K * ‖z - α‖ ^ 2 * ((‖lam - 1‖ ^ 2 / 2) * ‖z - α‖) := by ring
    _ ≤ K * ‖z - α‖ ^ 2 * ‖steffensen_denom_C x p z‖ := by
        apply mul_le_mul_of_nonneg_left h_denom_lb
        positivity

/-! ============================================================
  §31.2  Unconditional `steffensen_local_attraction`
============================================================ -/

/-- **★★ Unconditional Steffensen local attraction.**

  The hypothesis `h_quadratic` of §27.3 `steffensen_local_attraction`
  is discharged by §31.1: at any Pandrosion fixed point α, the
  Steffensen iterator has an explicit super-attracting basin around
  α with no analytic hypothesis required. -/
theorem steffensen_local_attraction_unconditional
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    ∃ r > 0, ∀ z₀ : ℂ, ‖z₀ - α‖ < r →
      Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop (𝓝 α) := by
  obtain ⟨K, hK_pos, r₀, hr₀_pos, h_quad⟩ :=
    steffensen_step_C_quadratic_bound x hx p hp α hα hSp hfp
  -- Shrink r₀ if needed so that K · r ≤ 1/2.
  set r : ℝ := min r₀ (1 / (2 * K))
  have h_halfK_pos : 0 < 1 / (2 * K) := by positivity
  have hr_pos : 0 < r := lt_min hr₀_pos h_halfK_pos
  have h_Kr_le : K * r ≤ 1 / 2 := by
    have : K * r ≤ K * (1 / (2 * K)) :=
      mul_le_mul_of_nonneg_left (min_le_right _ _) (le_of_lt hK_pos)
    have h_eq : K * (1 / (2 * K)) = 1 / 2 := by
      field_simp
      ring
    linarith
  have h_quad_r : ∀ z : ℂ, ‖z - α‖ < r →
      ‖steffensen_step_C x p z - α‖ ≤ K * ‖z - α‖ ^ 2 := by
    intro z hz
    exact h_quad z (lt_of_lt_of_le hz (min_le_left _ _))
  refine ⟨r, hr_pos, ?_⟩
  exact steffensen_local_attraction x hx p α hSp hfp hr_pos hK_pos h_Kr_le h_quad_r

/-! ============================================================
  §31.3  Unconditional `steffensen_dynamical_convergence_ae`
============================================================ -/

/-- **★★★ Unconditional Steffensen dynamical convergence a.e.**

  The per-anchor `h_quadratic` hypothesis of §28.3
  `steffensen_dynamical_convergence_ae` is discharged by §31.1.
  Only the a.e. trajectory-entry hypothesis `h_entry` remains,
  which is the genuine global dynamical content (Fatou/McMullen
  at the measure level).

  Conclusion: modulo the a.e. trajectory-entry hypothesis, the
  Steffensen orbit converges in ℂ to a cyclotomic anchor for
  almost every `z₀`. -/
theorem steffensen_dynamical_convergence_ae_unconditional
    (p : ℕ) (hp : 1 ≤ p)
    (x : ℂ) (hx : x ≠ 0)
    (α : ℂ) (hα : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ s : Fin p, Sp_C p (cycAnchor α p s) ≠ 0) :
    ∃ r : Fin p → ℝ, (∀ s, 0 < r s) ∧
      ((∀ᵐ z₀ : ℂ ∂volume, ∃ s : Fin p, ∃ k₀ : ℕ,
          ‖(steffensen_step_C x p)^[k₀] z₀ - cycAnchor α p s‖ < r s) →
        ∀ᵐ z₀ : ℂ ∂volume, ∃ s : Fin p,
          Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop
            (𝓝 (cycAnchor α p s))) := by
  -- For each anchor, extract the local quadratic bound via §31.1.
  classical
  -- Choose per-anchor constants using choice.
  have h_per_s : ∀ s : Fin p, ∃ K r, 0 < K ∧ 0 < r ∧ K * r ≤ 1 / 2 ∧
      ∀ z : ℂ, ‖z - cycAnchor α p s‖ < r →
        ‖steffensen_step_C x p z - cycAnchor α p s‖ ≤ K * ‖z - cycAnchor α p s‖ ^ 2 := by
    intro s
    have hα_s : cycAnchor α p s ≠ 0 := cycAnchor_ne_zero hα p s
    have hfp_s : (cycAnchor α p s) ^ p = 1 / x := by
      rw [cycAnchor_pow α hp s]; exact hα_pow
    obtain ⟨K, hK_pos, r₀, hr₀_pos, h_quad⟩ :=
      steffensen_step_C_quadratic_bound x hx p hp (cycAnchor α p s) hα_s (hSp s) hfp_s
    set r : ℝ := min r₀ (1 / (2 * K))
    have h_halfK_pos : 0 < 1 / (2 * K) := by positivity
    have hr_pos : 0 < r := lt_min hr₀_pos h_halfK_pos
    have h_Kr_le : K * r ≤ 1 / 2 := by
      have : K * r ≤ K * (1 / (2 * K)) :=
        mul_le_mul_of_nonneg_left (min_le_right _ _) (le_of_lt hK_pos)
      have h_eq : K * (1 / (2 * K)) = 1 / 2 := by
        field_simp; ring
      linarith
    refine ⟨K, r, hK_pos, hr_pos, h_Kr_le, ?_⟩
    intro z hz
    exact h_quad z (lt_of_lt_of_le hz (min_le_left _ _))
  choose K_fn r_fn hK_pos hr_pos h_Kr_le h_quad using h_per_s
  refine ⟨r_fn, hr_pos, ?_⟩
  intro h_entry
  exact steffensen_dynamical_convergence_ae p hp x hx α hα_pow hSp r_fn K_fn
    hr_pos hK_pos h_Kr_le h_quad h_entry

end Pandrosion
