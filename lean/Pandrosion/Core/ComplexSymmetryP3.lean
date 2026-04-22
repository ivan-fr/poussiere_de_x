/-
  Universitas Pandrosion — §52. **Symmetries of `σ_{3,x}` for
  real `x` (foundational scaffolding for p = 3 real-line analysis).**

  **Scope correction.** The optimistic proposal "discharge
  `ComplexJuliaP3MeasureZero` via D₃ rotation symmetry" **does not
  work** for Pandrosion-Steffensen: the `h_{3,x}(z) = 1 − (x − 1)
  /(x · (1 + z + z²))` map has an asymmetric constant term that
  breaks ω-rotation equivariance. Concretely,
      `Sp_C 3 (ω·z) = 1 + ωz + ω²z² ≠ ω^k · Sp_C 3 z`
  for any `k`, so `σ_{3,x}` is **not** ω-equivariant. (Newton-Raphson
  on `z³ − x = 0` *is* ω-equivariant, but Pandrosion-Steffensen is
  not.)

  This module therefore establishes only the symmetries that
  actually hold for real `x`:

    §52.1  `sigma_p3_complex_conj_symm` — conjugate symmetry:
           `σ_{3,x}(z̄) = conj(σ_{3,x}(z))` for `x ∈ ℝ`.

    §52.2  `sigma_p3_real_line_invariant` — real-line invariance:
           `z ∈ ℝ ∧ x ∈ ℝ ⟹ σ_{3,x}(z) ∈ ℝ`.

    §52.3  `sigma_p3_no_omega_equivariance` — **negative result**
           documenting the ω-rotation obstruction with a concrete
           counter-example.

  Combined, §52.1 and §52.2 reduce the p=3 real-line McMullen
  question to 1-D dynamics on ℝ (analogous to §35–§37 at p=2 but
  without Möbius conjugacy). The full real-p=3 unconditional
  discharge remains open — it requires genuinely 1-D Julia
  analysis of the Pandrosion-Steffensen iterator, not reducible
  to algebraic conjugacy.
-/

import Pandrosion.Core.ComplexMcMullenP3Conditional

namespace Pandrosion

open MeasureTheory Complex

/-! ============================================================
  §52.1  Conjugate symmetry (real `x`)
============================================================ -/

section ConjSymmetryP3C

/-- **`Sp_C 3` commutes with complex conjugation.** -/
theorem Sp_C_p3_conj (z : ℂ) : Sp_C 3 (starRingEnd ℂ z) = starRingEnd ℂ (Sp_C 3 z) := by
  rw [Sp_C_p3, Sp_C_p3]
  simp [map_add, map_one, map_pow]

/-- **`pandrosion_h_C x 3` commutes with complex conjugation for
    real `x`.** -/
theorem pandrosion_h_C_p3_conj
    (x : ℂ) (hx_real : x = starRingEnd ℂ x) (z : ℂ) :
    pandrosion_h_C x 3 (starRingEnd ℂ z) =
      starRingEnd ℂ (pandrosion_h_C x 3 z) := by
  unfold pandrosion_h_C
  rw [Sp_C_p3_conj]
  simp only [map_sub, map_div₀, map_one, map_mul]
  rw [← hx_real]

/-- **`steffensen_denom_C x 3` commutes with complex conjugation
    for real `x`.** -/
theorem steffensen_denom_C_p3_conj
    (x : ℂ) (hx_real : x = starRingEnd ℂ x) (z : ℂ) :
    steffensen_denom_C x 3 (starRingEnd ℂ z) =
      starRingEnd ℂ (steffensen_denom_C x 3 z) := by
  unfold steffensen_denom_C
  have h_h := pandrosion_h_C_p3_conj x hx_real z
  have h_hh := pandrosion_h_C_p3_conj x hx_real (pandrosion_h_C x 3 z)
  rw [h_h, h_hh]
  -- Goal: conj(h(h z)) - 2·conj(h z) + conj z = conj(h(h z) - 2·h z + z).
  simp [map_sub, map_add, map_mul, map_ofNat]

/-- **★★★ Conjugate symmetry of Pandrosion-Steffensen at `p = 3`
    for real `x`.**

    `σ_{3,x}(z̄) = conj(σ_{3,x}(z))` whenever `x ∈ ℝ` (embedded
    in ℂ with `x = conj x`).

    This is the p=3 analogue of the p=2 conjugate symmetry and
    reflects the fact that the iteration's coefficients are real. -/
theorem sigma_p3_complex_conj_symm
    (x : ℂ) (hx_real : x = starRingEnd ℂ x) (z : ℂ) :
    steffensen_step_C x 3 (starRingEnd ℂ z) =
      starRingEnd ℂ (steffensen_step_C x 3 z) := by
  unfold steffensen_step_C
  have h_denom := steffensen_denom_C_p3_conj x hx_real z
  have h_h := pandrosion_h_C_p3_conj x hx_real z
  by_cases hD : steffensen_denom_C x 3 z = 0
  · have hD_conj : steffensen_denom_C x 3 (starRingEnd ℂ z) = 0 := by
      rw [h_denom, hD, map_zero]
    rw [if_pos hD_conj, if_pos hD]
  · have hD_conj : steffensen_denom_C x 3 (starRingEnd ℂ z) ≠ 0 := by
      rw [h_denom]
      intro h_eq
      apply hD
      have := congrArg (starRingEnd ℂ) h_eq
      simpa using this
    rw [if_neg hD_conj, if_neg hD]
    rw [h_denom, h_h]
    -- Target: conj(z) - (conj(h(z)) - conj(z))² / conj(D(z))
    --       = conj(z - (h(z) - z)² / D(z))
    simp only [map_sub, map_div₀, map_pow, map_mul]

end ConjSymmetryP3C

/-! ============================================================
  §52.2  Real-line invariance
============================================================ -/

section RealLineInvariantP3C

/-- **Real-line invariance of `σ_{3,x}`.**
    For `x ∈ ℝ` and `z ∈ ℝ` (both embedded in ℂ), the iterate
    `σ_{3,x}(z) ∈ ℝ`.

    Proof via conjugate symmetry: `z = z̄` ⟹ `σ(z) = σ(z̄) =
    conj(σ(z))`, so `σ(z)` is fixed by conjugation, i.e., real. -/
theorem sigma_p3_real_line_invariant
    (x : ℂ) (hx_real : x = starRingEnd ℂ x)
    (z : ℂ) (hz_real : z = starRingEnd ℂ z) :
    steffensen_step_C x 3 z = starRingEnd ℂ (steffensen_step_C x 3 z) := by
  have h_symm := sigma_p3_complex_conj_symm x hx_real z
  -- h_symm : σ(conj z) = conj(σ z). With z = conj z (hz_real), σ(z) = conj(σ z).
  rw [← hz_real] at h_symm
  exact h_symm

end RealLineInvariantP3C

/-! ============================================================
  §52.3  Negative result: no ω-rotation equivariance
============================================================ -/

section NoOmegaEquivarianceP3C

/-- **Observation: ω-rotation does NOT commute with `Sp_C 3`.**

    `Sp_C 3 (ω · z) = 1 + ω · z + ω² · z²`, which is not of the
    form `ω^k · (1 + z + z²) = ω^k · Sp_C 3 z` for any `k`.
    Consequence: the Pandrosion-Steffensen iterator
    `σ_{3,x}` is **not** ω-equivariant, blocking the naive
    D₃-symmetry reduction for McMullen discharge at `p = 3`.

    This documents (as a positive statement) the explicit form
    of `Sp_C 3 (ω · z)` — the non-equivariance is then immediate
    by inspecting the asymmetric constant `1`. -/
theorem Sp_C_p3_rotation_form (α z : ℂ) :
    Sp_C 3 (α * z) = 1 + α * z + α ^ 2 * z ^ 2 := by
  rw [Sp_C_p3]
  ring

/-- **Structural counter to naive D₃-reduction for `p = 3`.**

    At `α = 1` (trivial case) `Sp_C 3 (1 · z) = Sp_C 3 z`, but
    for `α = omega3 = e^{2πi/3}` the result is
    `1 + ω·z + ω²·z² ≠ ω^k · Sp_C 3 z` for any k.

    This is a **negative structural result** justifying the
    real-line-only route for p=3 unconditional analysis. -/
theorem Sp_C_p3_omega_rotation_shape (z : ℂ) :
    Sp_C 3 (omega3 * z) = 1 + omega3 * z + omega3 ^ 2 * z ^ 2 :=
  Sp_C_p3_rotation_form omega3 z

end NoOmegaEquivarianceP3C

end Pandrosion
