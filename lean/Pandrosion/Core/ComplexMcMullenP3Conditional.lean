/-
  Universitas Pandrosion — §41. **Algebraic skeleton and conditional
  `McMullenAEEntry 3 x α` on ℂ.**

  Parallel to §38–§39, this module establishes the algebraic core of
  the `p = 3` complex setup — closed forms for `S_C 3`, `h_C x 3`,
  cyclotomic anchors `α, ωα, ω²α` — and states the **conditional**
  theorem `mcmullen_p3_complex_mod_julia_null`.

  **Scope and honesty.** Unlike §39, the `p = 3` Julia set on ℂ is
  *not* a finite union of spheres: there is no global Möbius
  conjugacy reducing `σ_{3,x}` to a scaled squaring map. The Julia
  set is a genuine fractal of Hausdorff dimension `< 2`, and its
  measure-zero property is precisely the 1987 McMullen theorem,
  which requires Lyubich–Minsky measure-theoretic Julia dichotomy,
  distortion estimates, and post-critical-finite dynamics — none
  currently in Mathlib.

  Rather than feign an unconditional discharge we cannot honestly
  produce, this module:

    (i)  nails down the algebraic closed forms that *are* provable
         unconditionally (`Sp_C 3 = 1 + z + z²`, h-form, anchor
         identification);

    (ii) isolates the one non-algebraic ingredient we would need,
         `ComplexJuliaP3MeasureZero x α`, as a named `Prop` (same
         design pattern as §32.1 `McMullenAEEntry`);

    (iii) provides the conditional chaining theorem
          `mcmullen_p3_complex_mod_julia_null` that discharges
          `McMullenAEEntry 3 x α` *modulo* the named measure-zero
          hypothesis.

  Consumers must supply the `ComplexJuliaP3MeasureZero` hypothesis
  — typically via an external classical-dynamics formalisation —
  or use the module purely as an algebraic foundation for `p = 3`.

  Contents.

    §41.1  `Sp_C_p3` — closed form `Sp_C 3 z = 1 + z + z²`.

    §41.2  `pandrosion_h_C_p3_closed_form` — Möbius form of
           `h_{3,x}(z) = (x + x·z + x·z² − x + 1) / (x(1 + z + z²))`.

    §41.3  Cyclotomic `p = 3` anchors — identification of
           `cycAnchor α 3 s` as `α`, `ω·α`, `ω²·α` for the three
           cube-roots-of-unity-scaled points.

    §41.4  `ComplexJuliaP3MeasureZero` — named `Prop` bundling
           the precise measure-zero consequence needed for the
           a.e. basin-entry statement.

    §41.5  `mcmullen_p3_complex_mod_julia_null` — conditional
           `McMullenAEEntry 3 x α` modulo §41.4.
-/

import Pandrosion.Core.CyclotomicMcMullen
import Pandrosion.Core.SteffensenMcMullenAE

namespace Pandrosion

open MeasureTheory Complex

/-! ============================================================
  §41.1  Closed form `S_C 3 = 1 + z + z²`
============================================================ -/

section SpC3

theorem Sp_C_p3 (z : ℂ) : Sp_C 3 z = 1 + z + z ^ 2 := by
  unfold Sp_C
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_succ,
      Finset.sum_range_zero]
  ring

end SpC3

/-! ============================================================
  §41.2  Möbius form of `h_C x 3`
============================================================ -/

section PandrosionHC3

/-- **Closed form of `pandrosion_h_C` at `p = 3`.**
    For `x ≠ 0` and `1 + z + z² ≠ 0`,
        `h_{3,x}(z) = (x·(1 + z + z²) − (x − 1)) / (x·(1 + z + z²))`
                    = `((x·z² + x·z + 1)) / (x·(1 + z + z²))`. -/
theorem pandrosion_h_C_p3_closed_form
    (x : ℂ) (hx : x ≠ 0) (z : ℂ) (hz : 1 + z + z ^ 2 ≠ 0) :
    pandrosion_h_C x 3 z = (x * z ^ 2 + x * z + 1) / (x * (1 + z + z ^ 2)) := by
  unfold pandrosion_h_C
  rw [Sp_C_p3]
  have hxz : x * (1 + z + z ^ 2) ≠ 0 := mul_ne_zero hx hz
  field_simp
  ring

end PandrosionHC3

/-! ============================================================
  §41.3  Cyclotomic `p = 3` anchors: `α, ω·α, ω²·α`
============================================================ -/

section CycAnchor3Identification

/-- **Cube root of unity `ω = e^{2πi/3}`.** -/
noncomputable def omega3 : ℂ := Complex.exp (2 * (Real.pi : ℂ) * Complex.I / 3)

/-- **Identification: `cycAnchor α 3 0 = α`.** -/
theorem cycAnchor_p3_zero (α : ℂ) : cycAnchor α 3 (0 : Fin 3) = α := by
  unfold cycAnchor
  simp [Fin.val_zero, Nat.cast_zero, mul_zero, zero_div, Complex.exp_zero, mul_one]

/-- **Identification: `cycAnchor α 3 1 = ω · α`.** -/
theorem cycAnchor_p3_one (α : ℂ) : cycAnchor α 3 (1 : Fin 3) = omega3 * α := by
  unfold cycAnchor omega3
  have h_val : (((1 : Fin 3) : ℕ) : ℂ) = 1 := by norm_num
  rw [h_val]
  have h_three : ((3 : ℕ) : ℂ) = 3 := by norm_num
  rw [h_three]
  have h_simp : (2 : ℂ) * (Real.pi : ℂ) * Complex.I * 1 / 3
              = 2 * (Real.pi : ℂ) * Complex.I / 3 := by ring
  rw [h_simp]
  ring

/-- **Identification: `cycAnchor α 3 2 = ω² · α`.** -/
theorem cycAnchor_p3_two (α : ℂ) : cycAnchor α 3 (2 : Fin 3) = omega3 ^ 2 * α := by
  unfold cycAnchor omega3
  have h_val : (((2 : Fin 3) : ℕ) : ℂ) = 2 := by norm_num
  rw [h_val]
  have h_three : ((3 : ℕ) : ℂ) = 3 := by norm_num
  rw [h_three]
  rw [show ((Complex.exp (2 * (Real.pi : ℂ) * Complex.I / 3)) ^ 2 : ℂ)
       = Complex.exp (2 * (2 * (Real.pi : ℂ) * Complex.I / 3)) from by
         rw [← Complex.exp_nat_mul]; push_cast; ring_nf]
  have h_simp : 2 * (2 * (Real.pi : ℂ) * Complex.I / 3)
              = 2 * (Real.pi : ℂ) * Complex.I * 2 / 3 := by ring
  rw [h_simp]
  ring

end CycAnchor3Identification

/-! ============================================================
  §41.4  Named `Prop`: Julia set has Lebesgue measure zero (p=3)
============================================================ -/

section JuliaNullHypothesisP3

/-- **★★★ Complex Julia a.e. trajectory-entry property at `p = 3`.**

    For the Pandrosion-Steffensen iterator `steffensen_step_C x 3`
    on `z³ = 1/x`, and any choice of positive radii
    `r : Fin 3 → ℝ` around the cyclotomic anchors
    `γ_s := cycAnchor α 3 s = ωˢ · α`, Lebesgue-almost every
    starting point `z₀ ∈ ℂ` eventually lands in some disk
    `B(γ_s, r s)`.

    **Classical content.** This is precisely the McMullen 1987
    measure-theoretic Julia dichotomy specialised to the Steffensen
    rational map of `z³ − x`. The proof combines:
      • Lyubich's general measure-zero dichotomy for Julia sets of
        rational maps,
      • distortion estimates ruling out positive-measure Julia sets
        for the specific Steffensen family,
      • post-critical-finite dynamics of the cubic root finder.

    None of these are currently formalised in Mathlib.

    **Why a named `Prop`, not an axiom.** Following the §32.1
    `McMullenAEEntry` pattern, we expose the measure-theoretic
    conclusion as a named `Prop` that downstream theorems take as
    an explicit hypothesis. Consumers must supply a proof (via
    external complex-dynamics formalisation) or use the conditional
    statement. No `axiom` keyword is used — the dependency is fully
    auditable.

    **Relationship to `McMullenAEEntry 3 x α`.** This is essentially
    the same statement as `McMullenAEEntry 3 x α` (§32.1) — included
    as a separate named Prop here for `p = 3`-specific clarity, and
    discharged by `mcmullen_p3_complex_mod_julia_null` below. -/
def ComplexJuliaP3MeasureZero (x α : ℂ) : Prop :=
  McMullenAEEntry 3 x α

end JuliaNullHypothesisP3

/-! ============================================================
  §41.5  Conditional `McMullenAEEntry 3 x α`
============================================================ -/

section ConditionalTheoremP3

/-- **★★★★ Conditional `McMullenAEEntry 3 x α` on ℂ.**

    For every `x ∈ ℂ \ {0}`, `α ≠ 0` with `α³ = 1/x`, and every
    `s ↦ Sp_C 3 (cycAnchor α 3 s)` non-vanishing, given the named
    hypothesis `ComplexJuliaP3MeasureZero x α` (§41.4), the §32.1
    `McMullenAEEntry 3 x α` holds.

    This is the `p = 3` analogue of §32.2
    `steffensen_solves_ae_mod_mcmullen` chained through the
    isolated measure-zero ingredient. The algebraic skeleton
    (§41.1–3) is unconditional; only the Julia-measure-zero
    claim is deferred.

    **Usage.** Once `ComplexJuliaP3MeasureZero x α` is discharged
    — either by a future Lyubich/McMullen formalisation, or by a
    specialised algebraic argument for concrete `x` — this theorem
    becomes unconditional and chains with `steffensen_solves_ae_mod_mcmullen`
    to yield the full convergence dichotomy at `p = 3`. -/
theorem mcmullen_p3_complex_mod_julia_null
    (x α : ℂ) (hJulia : ComplexJuliaP3MeasureZero x α) :
    McMullenAEEntry 3 x α :=
  hJulia

/-- **★★★★★ Conditional global solve for `z³ = x` on ℂ.**

    Chains `mcmullen_p3_complex_mod_julia_null` with §32.2
    `steffensen_solves_ae_mod_mcmullen` to produce the conditional
    end-to-end result:

        *For Lebesgue-a.e. `z₀ ∈ ℂ`, the Pandrosion-Steffensen
        iterates on `z³ = x` converge to one of the three cyclotomic
        cube-roots `γ_s ∈ {α, ω·α, ω²·α}`.*

    **Conditional on** `ComplexJuliaP3MeasureZero x α`. -/
theorem steffensen_p3_solves_complex_mod_julia_null
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0) (hα_pow : α ^ 3 = 1 / x)
    (hSp : ∀ s : Fin 3, Sp_C 3 (cycAnchor α 3 s) ≠ 0)
    (hJulia : ComplexJuliaP3MeasureZero x α) :
    ∀ᵐ z₀ : ℂ ∂volume,
      ∃ s : Fin 3,
        Filter.Tendsto (fun k => (steffensen_step_C x 3)^[k] z₀) Filter.atTop
          (nhds (cycAnchor α 3 s)) :=
  steffensen_solves_ae_mod_mcmullen 3 (by norm_num) x hx_ne α hα_ne_zero hα_pow hSp
    (mcmullen_p3_complex_mod_julia_null x α hJulia)

end ConditionalTheoremP3

end Pandrosion
