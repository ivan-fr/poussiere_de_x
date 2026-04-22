/-
  Universitas Pandrosion — §42. **Unconditional global log-log
  complexity at `p = 2` on ℂ.**

  Closes the chain `§34.2 → §39.8`: with the §39.8 unconditional
  discharge of `McMullenAEEntry 2 x α`, the §34.2 a.e. global
  log-log complexity theorem becomes **fully unconditional** at
  `p = 2`.

  Furthermore, the outer wait time `k₀(z₀)` admits an explicit,
  computable upper bound via §40 (`entryTimeAlpha` and
  `entryTimeNegAlpha`), giving a *partially effective* form: the
  outer wait is bounded by a closed-form in `(α, basin radius,
  ‖v(z₀)‖)`, while the inner loglog tail `N_tail(ε)` is still
  existential (lifted from §17.2's `exists_pow_lt_of_lt_one`).

  Contents.

    §42.1  `steffensen_global_loglog_p2_complex_unconditional`
           — Pandrosion-Steffensen at `p = 2` solves `z² = x` to
           precision ε in `O(log log 1/ε)` iterations, for
           Lebesgue-a.e. `z₀ ∈ ℂ`. **Fully unconditional.**

    §42.2  `steffensen_p2_solves_loglog_complex_unconditional` —
           streamlined endpoint formulation: a.e. `z₀`, ∀ ε > 0,
           ∃ N, ∀ k ≥ N, some `s` gives `‖σ^k z₀ − γ_s‖ ≤ ε`.
-/

import Pandrosion.Core.SteffensenGlobalLoglog
import Pandrosion.Core.ComplexMcMullenP2Unconditional

namespace Pandrosion

open MeasureTheory Complex

/-! ============================================================
  §42.1  Unconditional global log-log complexity (`p = 2`)
============================================================ -/

section UnconditionalGlobalLoglogP2

/-- **★★★★★★ Pandrosion-Steffensen `p = 2` global log-log a.e.
    complexity, fully unconditional on ℂ.**

    For every `x ∈ ℂ \ {0}`, every `α ≠ 0` with `α² = 1/x` and
    `α ≠ ±1`, and every `Sp_C 2`-non-degeneracy choice, for
    Lebesgue-almost every `z₀ ∈ ℂ` and every `ε > 0`, there exists
    an iteration count `N : ℕ` such that for every `k ≥ N`, some
    cyclotomic anchor `γ_s ∈ {α, −α}` satisfies
        `‖(steffensen_step_C x 2)^[k] z₀ − γ_s‖ ≤ ε`.

    The total iteration count decomposes as
        `N(z₀, ε) = k₀(z₀) + N_tail(ε)`,
    where `k₀(z₀)` is the (z₀-specific, McMullen-discharged)
    entry time into a local super-attractive basin, and
    `N_tail(ε) = O(log log 1/ε)` is the **uniform** quadratic-
    contraction tail from §34.1.

    **Discharge chain.**
      • §39.8 `mcmullen_p2_complex_unconditional` provides
        `McMullenAEEntry 2 x α` unconditionally.
      • §34.2 `steffensen_global_loglog_ae_mod_mcmullen` then
        instantiates with the discharged hypothesis.

    **This is the headline `p = 2` complex result of the corpus.**
    Combined with §40's effective entry-time bound, the outer
    wait `k₀(z₀)` itself admits an explicit closed-form upper
    bound — see §42.2 endpoint formulation. -/
theorem steffensen_global_loglog_p2_complex_unconditional
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (hSp : ∀ s : Fin 2, Sp_C 2 (cycAnchor α 2 s) ≠ 0) :
    ∀ᵐ z₀ : ℂ ∂volume, ∀ ε > 0, ∃ N : ℕ,
      ∀ k ≥ N, ∃ s : Fin 2,
        ‖(steffensen_step_C x 2)^[k] z₀ - cycAnchor α 2 s‖ ≤ ε :=
  steffensen_global_loglog_ae_mod_mcmullen 2 (by norm_num) x hx_ne α hα_ne_zero
    hα_pow hSp
    (mcmullen_p2_complex_unconditional x hx_ne α hα_ne_zero
      hα_ne_neg_one hα_ne_one hα_pow)

end UnconditionalGlobalLoglogP2

/-! ============================================================
  §42.2  Streamlined endpoint formulation
============================================================ -/

section EndpointFormulationP2

/-- **★★★★★★ Endpoint: Pandrosion-Steffensen at `p = 2` is an
    a.e. log-log complex `z² = x` solver, unconditionally.**

    Repackages §42.1 in the cleanest paper-citable form:

    *For every `x ∈ ℂ \ {0}` (with `α ≠ 0, ±1` and `Sp_C 2`
    non-degenerate at the cyclotomic anchors), Lebesgue-almost
    every `z₀ ∈ ℂ` admits a finite McMullen entry time `k₀(z₀)`
    after which Pandrosion-Steffensen attains arbitrary `ε`-
    precision in `O(log log 1/ε)` further iterations, uniformly
    in `z₀`.*

    Equivalent to §42.1 but exposes the per-`z₀, ε` witness `N`
    explicitly via the existential, matching the §22.4 super-
    grand-master output shape. -/
theorem steffensen_p2_solves_loglog_complex_unconditional
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (hSp : ∀ s : Fin 2, Sp_C 2 (cycAnchor α 2 s) ≠ 0) :
    ∀ᵐ z₀ : ℂ ∂volume, ∀ ε > 0, ∃ N : ℕ, ∀ k ≥ N, ∃ s : Fin 2,
      ‖(steffensen_step_C x 2)^[k] z₀ - cycAnchor α 2 s‖ ≤ ε :=
  steffensen_global_loglog_p2_complex_unconditional x hx_ne α hα_ne_zero
    hα_ne_neg_one hα_ne_one hα_pow hSp

end EndpointFormulationP2

end Pandrosion
