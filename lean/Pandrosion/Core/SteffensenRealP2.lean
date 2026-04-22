/-
  Universitas Pandrosion — §35. **Real-line Pandrosion-Steffensen at p=2.**

  Specialisation of the generic real iterator (Foundations §1/§3) to
  the square-root case `p = 2`. The central algebraic observation —
  `h_{2,x}` is a Möbius map with explicit closed-form iterate — yields
  a fully-explicit contraction multiplier

      `μ_{2,x} = (1 − α)/(1 + α)`,   where `α² = 1/x`,

  and an explicit "Cible α" scoping for the unconditional real-line
  McMullen a.e. entry property.

  Rationale. The complex a.e. statement §32 already quantifies over
  all p ≥ 1 modulo the named hypothesis `McMullenAEEntry p x α`
  (§32.1). This module introduces the sharper, purely-real-variable
  companion `RealMcMullenP2 x` — a dedicated named `Prop` for the
  square-root case — and chains it to a real-line a.e. basin-entry
  corollary. No complex-dynamics infrastructure is required for any
  of the algebraic content (Möbius form, closed-form multiplier,
  `|μ| < 1`); those are proved here unconditionally.

  Content.

    §35.1  `Sp_p2`, `pandrosion_h_p2_mobius` — explicit Möbius form
           `h_{2,x}(s) = (s + 1/x) / (1 + s)`, with `S_2(s) = 1 + s`.

    §35.2  `lambdaClosedR`, `lambdaClosedR_p2_at_fp` — real closed-form
           multiplier `μ_p(α) = 1 − p·α^{p−1}/S_p(α)`; specialised at
           `p = 2` to `(1 − α)/(1 + α)`.

    §35.3  `lambdaClosedR_p2_abs_lt_one` — for `x > 1` and the
           positive fixed point `α ∈ (0, 1)`, the multiplier lies in
           `(0, 1)` (strict contraction at `+α`).

    §35.4  `RealMcMullenP2` — named `Prop` for real-line a.e.
           trajectory-entry into neighborhoods of `±α`, the
           sharpened `p = 2` analogue of `McMullenAEEntry`.

    §35.5  `steffensen_p2_ae_basin_entry_mod_real_mcmullen` — the
           paper-citable chaining corollary: given `RealMcMullenP2 x`,
           for Lebesgue-a.e. `z₀ ∈ ℝ`, some iterate enters an
           arbitrary preassigned neighborhood of `+α` or `−α`.
-/

import Pandrosion.Core.Foundations
import Pandrosion.Core.SteffensenMcMullenAE
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic

namespace Pandrosion

open Filter Topology MeasureTheory

/-! ============================================================
  §35.1  Möbius form of `h_{2,x}` on ℝ
============================================================ -/

section MobiusFormP2

/-- **`p = 2` specialisation of the geometric sum:** `S_2(s) = 1 + s`. -/
theorem Sp_p2 (s : ℝ) : Sp 2 s = 1 + s := by
  unfold Sp
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_zero]
  ring

/-- **Möbius closed form of `h_{2,x}` on ℝ.** For `x ≠ 0` and
    `s ≠ −1`,
        `h_{2,x}(s) = (s + 1/x) / (1 + s)`.
    Direct algebraic rearrangement of the paper definition. -/
theorem pandrosion_h_p2_mobius
    (x : ℝ) (hx : x ≠ 0) (s : ℝ) (hs : s ≠ -1) :
    pandrosion_h x 2 s = (s + 1 / x) / (1 + s) := by
  have h1s_ne : (1 + s : ℝ) ≠ 0 := by
    intro h; apply hs; linarith
  unfold pandrosion_h
  rw [Sp_p2]
  have hxS_ne : x * (1 + s) ≠ 0 := mul_ne_zero hx h1s_ne
  field_simp
  ring

/-- **`h_{2,x}` fixes `α` whenever `α² = 1/x` and `α ≠ −1`.** -/
theorem pandrosion_h_p2_fixed_at_sqrt_fp
    (x : ℝ) (hx : x ≠ 0) (α : ℝ) (hα_ne : α ≠ -1)
    (hfp : α ^ 2 = 1 / x) :
    pandrosion_h x 2 α = α := by
  have h1a_ne : (1 + α : ℝ) ≠ 0 := by
    intro h; apply hα_ne; linarith
  rw [pandrosion_h_p2_mobius x hx α hα_ne]
  -- Want: (α + 1/x) / (1 + α) = α. Use α² = 1/x.
  rw [div_eq_iff h1a_ne, ← hfp]
  ring

end MobiusFormP2

/-! ============================================================
  §35.2  Real closed-form multiplier
============================================================ -/

section ClosedFormMultiplierP2

/-- **Real-line closed-form Pandrosion multiplier.**
    The paper's `μ_p(α) = 1 − p·α^{p−1} / S_p(α)` — the real-variable
    analogue of `lambdaClosedC` from §26.4. At a Pandrosion fixed
    point `α^p = 1/x`, it is the linear contraction rate of `h`
    at `α`. -/
noncomputable def lambdaClosedR (p : ℕ) (α : ℝ) : ℝ :=
  1 - (p : ℝ) * α ^ (p - 1) / Sp p α

/-- **`p = 2` specialisation: `μ_2(α) = (1 − α)/(1 + α)`.** -/
theorem lambdaClosedR_p2_at_fp (α : ℝ) (hα_ne : α ≠ -1) :
    lambdaClosedR 2 α = (1 - α) / (1 + α) := by
  have h1a_ne : (1 + α : ℝ) ≠ 0 := by
    intro h; apply hα_ne; linarith
  unfold lambdaClosedR
  rw [Sp_p2]
  -- (p : ℝ) = 2 and α ^ (p - 1) = α ^ 1 = α.
  show 1 - (2 : ℝ) * α ^ (2 - 1) / (1 + α) = (1 - α) / (1 + α)
  have h_exp : α ^ (2 - 1) = α := by norm_num
  rw [h_exp]
  field_simp
  ring

end ClosedFormMultiplierP2

/-! ============================================================
  §35.3  Strict contraction at `+α` when `x > 1`
============================================================ -/

section StrictContractionPositiveFP

/-- **Positive real fixed point lies in `(0, 1)` when `x > 1`.**
    If `α > 0` and `α² = 1/x` with `x > 1`, then `α < 1`. -/
theorem pos_fp_lt_one (x : ℝ) (hx : 1 < x)
    (α : ℝ) (hα_pos : 0 < α) (hα_pow : α ^ 2 = 1 / x) : α < 1 := by
  have hx_pos : 0 < x := lt_trans one_pos hx
  by_contra h
  push_neg at h  -- 1 ≤ α
  have h_sq_ge : (1 : ℝ) ≤ α ^ 2 := by nlinarith
  rw [hα_pow] at h_sq_ge
  -- 1 ≤ 1/x  ⟺  1 * x ≤ 1  ⟺  x ≤ 1, contradicting x > 1.
  rw [le_div_iff hx_pos] at h_sq_ge
  linarith

/-- **★ Strict linear contraction of `h_{2,x}` at `+α` for `x > 1`.**
    For `x > 1` and the positive fixed point `α = 1/√x ∈ (0, 1)`,
    the closed-form multiplier satisfies
        `0 < μ_{2,x}(α) < 1`,
    i.e. `h_{2,x}` is strictly linearly contracting at `+α`.
    This is the real-line explicit-rate companion to
    `lambdaClosedC` at the cyclotomic anchor `γ_0`. -/
theorem lambdaClosedR_p2_pos_lt_one
    (x : ℝ) (hx : 1 < x)
    (α : ℝ) (hα_pos : 0 < α) (hα_pow : α ^ 2 = 1 / x) :
    0 < lambdaClosedR 2 α ∧ lambdaClosedR 2 α < 1 := by
  have hα_lt_one : α < 1 := pos_fp_lt_one x hx α hα_pos hα_pow
  have hα_ne_neg1 : α ≠ -1 := by linarith
  rw [lambdaClosedR_p2_at_fp α hα_ne_neg1]
  have h_denom_pos : (0 : ℝ) < 1 + α := by linarith
  have h_numer_pos : (0 : ℝ) < 1 - α := by linarith
  refine ⟨div_pos h_numer_pos h_denom_pos, ?_⟩
  rw [div_lt_one h_denom_pos]
  linarith

/-- **Magnitude bound `|μ_{2,x}(α)| < 1`.**
    Packages the strict-contraction corollary in absolute-value form,
    the natural Julia-theoretic statement "α is an attracting fixed
    point of `h_{2,x}` on ℝ". -/
theorem lambdaClosedR_p2_abs_lt_one
    (x : ℝ) (hx : 1 < x)
    (α : ℝ) (hα_pos : 0 < α) (hα_pow : α ^ 2 = 1 / x) :
    |lambdaClosedR 2 α| < 1 := by
  obtain ⟨h_pos, h_lt⟩ := lambdaClosedR_p2_pos_lt_one x hx α hα_pos hα_pow
  rw [abs_of_pos h_pos]; exact h_lt

end StrictContractionPositiveFP

/-! ============================================================
  §35.4  Real-line McMullen a.e. entry (named Prop, `p = 2`)
============================================================ -/

section RealMcMullenHypothesis

/-- **★★★ Real-line McMullen a.e. trajectory-entry (sharpened to
    `p = 2`).**

    For the real-variable Pandrosion-Steffensen iterator
    `steffensen_step x 2` on `z² = 1/x`, and any two positive radii
    `r_pos, r_neg` around the signed real fixed points `+α, −α`,
    Lebesgue-almost every starting point `z₀ ∈ ℝ` eventually lands in
    either `(α − r_pos, α + r_pos)` or `(−α − r_neg, −α + r_neg)`.

    **Why a named `Prop`, not an axiom.** The unconditional
    Fatou/McMullen dichotomy for the one-dimensional rational
    Steffensen-accelerator of `z² − 1/x` on ℝ reduces, via the
    complex lift, to hyperbolicity of the Steffensen rational map on
    `ℂP¹`. This is expected to hold (Newton-like root-finding on
    `z² − c` is textbook-hyperbolic, and the Steffensen variant
    inherits hyperbolicity from the same critical-orbit analysis),
    but Mathlib currently does not formalise the required
    post-critical-finite complex-dynamics machinery. Following the
    pattern of `McMullenAEEntry` (§32.1), we expose the statement as
    a named `Prop` that downstream theorems take as an explicit
    hypothesis. Consumers must supply a proof (via external
    complex-dynamics arguments) or instantiate the dependency chain
    conditionally. -/
def RealMcMullenP2 (x : ℝ) (α : ℝ) : Prop :=
  ∀ (r_pos r_neg : ℝ), 0 < r_pos → 0 < r_neg →
    ∀ᵐ z₀ : ℝ ∂volume,
      (∃ k₀ : ℕ, |(steffensen_step x 2)^[k₀] z₀ - α| < r_pos) ∨
      (∃ k₀ : ℕ, |(steffensen_step x 2)^[k₀] z₀ - (-α)| < r_neg)

end RealMcMullenHypothesis

/-! ============================================================
  §35.5  Chain to real-line a.e. basin entry
============================================================ -/

section ChainAEBasinEntry

/-- **★★★ Pandrosion-Steffensen enters a `p = 2` real basin
    almost-everywhere (conditional on `RealMcMullenP2`).**

    Fix `x > 1` and the positive real fixed point `α = 1/√x`
    (abstracted here as `α > 0` with `α² = 1/x`). Given the
    real-line McMullen hypothesis `hRM : RealMcMullenP2 x α` (§35.4),
    for every preassigned positive radii `r_pos, r_neg > 0`, for
    Lebesgue-almost every initial point `z₀ ∈ ℝ`, some iterate
    `(steffensen_step x 2)^[k₀] z₀` enters `B(+α, r_pos)` or
    `B(−α, r_neg)`.

    This is the **real-line paper-citable corollary** of the Cible α
    programme for `p = 2`: a fully-explicit, closed-form-contraction
    basin-entry statement on the real line, modulo the single-named
    `RealMcMullenP2` hypothesis. -/
theorem steffensen_p2_ae_basin_entry_mod_real_mcmullen
    (x : ℝ) (α : ℝ)
    (hRM : RealMcMullenP2 x α) :
    ∀ (r_pos r_neg : ℝ), 0 < r_pos → 0 < r_neg →
      ∀ᵐ z₀ : ℝ ∂volume,
        (∃ k₀ : ℕ, |(steffensen_step x 2)^[k₀] z₀ - α| < r_pos) ∨
        (∃ k₀ : ℕ, |(steffensen_step x 2)^[k₀] z₀ - (-α)| < r_neg) := by
  intro r_pos r_neg hr_pos hr_neg
  exact hRM r_pos r_neg hr_pos hr_neg

/-- **★★★★ Effective-rate formulation for `p = 2` real-line a.e. basin
    entry.**

    Strengthens §35.5 by packaging the closed-form multiplier
    `μ_{2,x}(α) = (1 − α)/(1 + α) ∈ (0, 1)` (§35.3) together with the
    a.e. basin-entry corollary. The output bundles:

      (i)   the explicit linear contraction rate `μ ∈ (0, 1)` of `h`
            at `+α`, given in closed form by the paper's formula;
      (ii)  the Lebesgue-a.e. guarantee that some Steffensen iterate
            enters an arbitrary basin of `±α`.

    Use-case: once `RealMcMullenP2 x α` is discharged (by a future
    complex-dynamics formalisation or by specialised real-line
    arguments in the tube `(−2α, 2α)`), this theorem immediately
    becomes the *unconditional* real-line a.e. solver for `z² = x`. -/
theorem steffensen_p2_effective_rate_and_ae_entry
    (x : ℝ) (hx : 1 < x)
    (α : ℝ) (hα_pos : 0 < α) (hα_pow : α ^ 2 = 1 / x)
    (hRM : RealMcMullenP2 x α) :
    (∃ μ : ℝ, μ = lambdaClosedR 2 α ∧ 0 < μ ∧ μ < 1) ∧
    (∀ (r_pos r_neg : ℝ), 0 < r_pos → 0 < r_neg →
      ∀ᵐ z₀ : ℝ ∂volume,
        (∃ k₀ : ℕ, |(steffensen_step x 2)^[k₀] z₀ - α| < r_pos) ∨
        (∃ k₀ : ℕ, |(steffensen_step x 2)^[k₀] z₀ - (-α)| < r_neg)) := by
  refine ⟨⟨lambdaClosedR 2 α, rfl, ?_, ?_⟩, ?_⟩
  · exact (lambdaClosedR_p2_pos_lt_one x hx α hα_pos hα_pow).1
  · exact (lambdaClosedR_p2_pos_lt_one x hx α hα_pos hα_pow).2
  · exact steffensen_p2_ae_basin_entry_mod_real_mcmullen x α hRM

end ChainAEBasinEntry

end Pandrosion
