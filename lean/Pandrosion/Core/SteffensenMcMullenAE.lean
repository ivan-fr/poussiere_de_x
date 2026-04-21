/-
  Universitas Pandrosion — §32. **Steffensen a.e. trajectory-entry
  (Fatou/McMullen).**

  The final missing ingredient for a fully unconditional
  "Pandrosion-Steffensen solves `z^p = x` almost everywhere"
  statement.

  Classical content (McMullen, *Families of rational maps and
  iterative root-finding algorithms*, Annals of Mathematics 125
  (1987), pp. 467–493). For the root-finding iteration associated
  with `z^p − x`, the Julia set has Lebesgue measure zero; hence,
  for Lebesgue-almost every initial point `z₀ ∈ ℂ`, the Steffensen
  orbit enters a neighborhood of one of the `p` cyclotomic roots.

  Status in this corpus. Proving McMullen's theorem unconditionally
  requires substantial complex-analytic machinery (Fatou/Julia
  theory, distortion estimates, post-critical-finite dynamics,
  Lyubich-Minsky laminations, …) not yet formalised in Mathlib.
  We therefore isolate the specific consequence we need as a
  named `Prop` and thread it as an explicit hypothesis through
  every theorem that depends on it. No `axiom` keyword is used:
  downstream consumers must supply the McMullen property at the
  call site, making the dependency fully auditable.

  Content.

    §32.1  `McMullenAEEntry` — named `Prop` bundling the a.e.
           trajectory-entry statement for the canonical cyclotomic
           anchor family of `z^p = 1/x`.

    §32.2  `steffensen_solves_ae_mod_mcmullen`
           — the final end-to-end theorem: given the McMullen
           hypothesis, for Lebesgue-almost every `z₀`, the
           Pandrosion-Steffensen iterates converge in ℂ to one of
           the cyclotomic roots `γ_s`. Chains §31.3 with the
           McMullen hypothesis to discharge the last remaining
           premise of `steffensen_dynamical_convergence_ae`.
-/

import Pandrosion.Core.SteffensenQuadraticBound

namespace Pandrosion

open Filter Topology MeasureTheory Complex

/-! ============================================================
  §32.1  McMullen a.e. trajectory-entry (named hypothesis)
============================================================ -/

section McMullenHypothesis

/-- **★★★ McMullen almost-everywhere trajectory-entry property.**

    For the Pandrosion-Steffensen iterator `steffensen_step_C x p`
    on `z^p = 1/x`, and any choice of positive radii
    `r : Fin p → ℝ` around the cyclotomic anchors
    `γ_s := cycAnchor α p s`, Lebesgue-almost every starting
    point `z₀ ∈ ℂ` eventually lands in some disk `B(γ_s, r s)`.

    **Classical content.** This is the Fatou/McMullen measure-
    theoretic dichotomy specialised to the family `z ↦ steffensen
    step for z^p − 1/x`. McMullen (1987) proved the general
    statement for iterative root-finding on polynomials of fixed
    degree: the Julia set of the iteration has Lebesgue measure
    zero, hence a.e. orbit accumulates at the post-critical set
    (here, the `p` roots).

    **Why a named `Prop`, not an axiom.** McMullen's proof
    combines (i) Lyubich's general measure-zero dichotomy for
    Julia sets of rational maps, (ii) distortion estimates that
    rule out positive-measure Julia sets for root-finding
    iterations of `z^p − x`, (iii) post-critical-finite dynamics
    — none of which Mathlib currently formalises. Instead of
    sealing this behind an `axiom` keyword, we expose the precise
    statement as a named `Prop` that downstream theorems take as
    an explicit hypothesis. Consumers are obligated to supply a
    proof (via external machinery) or to instantiate the
    dependency-chain conditionally. -/
def McMullenAEEntry
    (p : ℕ) (x : ℂ) (α : ℂ) : Prop :=
  ∀ (r : Fin p → ℝ), (∀ s, 0 < r s) →
    ∀ᵐ z₀ : ℂ ∂volume,
      ∃ s : Fin p, ∃ k₀ : ℕ,
        ‖(steffensen_step_C x p)^[k₀] z₀ - cycAnchor α p s‖ < r s

end McMullenHypothesis

/-! ============================================================
  §32.2  Final theorem — Pandrosion-Steffensen solves a.e.
============================================================ -/

section SolvesAE

/-- **★★★★ Pandrosion-Steffensen solves `z^p = x` almost everywhere
    (conditional on McMullen).**

    Fix `p ≥ 1`, `x ≠ 0`, and any `α ≠ 0` with `α^p = 1/x`. Assume
    `S_p(γ_s) ≠ 0` for every cyclotomic anchor `γ_s := cycAnchor α
    p s` (generic non-degeneracy — fails only on a proper algebraic
    subvariety of `x`).

    Given the McMullen hypothesis `hMcM : McMullenAEEntry p x α`
    (§32.1), for Lebesgue-almost every initial point `z₀ ∈ ℂ`,
    there exists an anchor index `s : Fin p` such that the
    Pandrosion-Steffensen iterates converge in ℂ to the cyclotomic
    root `γ_s`:
        `(steffensen_step_C x p)^[k] z₀  →  γ_s  as k → ∞`.

    This is the **paper-citable global result**: combining
    §31.1's unconditional local quadratic bound with McMullen's
    classical measure-theoretic theorem, the Pandrosion-Steffensen
    algorithm solves `z^p = x` dynamically (convergence, not
    merely basin occupancy) on a set of full Lebesgue measure.

    All other hypotheses (local quadratic bounds, `K_s · r_s ≤
    1/2`) are discharged unconditionally via §31. -/
theorem steffensen_solves_ae_mod_mcmullen
    (p : ℕ) (hp : 1 ≤ p)
    (x : ℂ) (hx : x ≠ 0)
    (α : ℂ) (hα : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ s : Fin p, Sp_C p (cycAnchor α p s) ≠ 0)
    (hMcM : McMullenAEEntry p x α) :
    ∀ᵐ z₀ : ℂ ∂volume,
      ∃ s : Fin p,
        Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop
          (𝓝 (cycAnchor α p s)) := by
  -- Extract the per-anchor unconditional quadratic bounds from §31.3.
  obtain ⟨r_fn, hr_pos, h_chain⟩ :=
    steffensen_dynamical_convergence_ae_unconditional
      p hp x hx α hα hα_pow hSp
  -- Discharge the a.e. trajectory-entry hypothesis via McMullen.
  exact h_chain (hMcM r_fn hr_pos)

end SolvesAE

end Pandrosion
