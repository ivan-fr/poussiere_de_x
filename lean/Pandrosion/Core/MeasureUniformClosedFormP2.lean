/-
  Universitas Pandrosion — §50. **Closed-form computable
  `K_star(α, r, N)` for `p = 2` complex uniform complexity.**

  Exposes §46/§49's uniform iteration count `K_star` as a **named
  function** (`mcmullen_K_star`) of `(α, r, N)`, instead of an
  existential witness buried inside `∃ K*, …`.

  This lets downstream code evaluate `K_star` concretely on
  specific `(α, r)` pairs — e.g. for numerical examples like
  `z² = 2` with `α = 1/√2`.

  Contents.

    §50.1  `mcmullen_K_star` — closed-form uniform iteration count.

    §50.2  `mcmullen_p2_complex_closed_form` — strict-uniform
           theorem with `K_star` exposed as the named function.
-/

import Pandrosion.Core.MeasureUniformStrictP2
import Pandrosion.Core.MeasureUniformEffectiveP2

namespace Pandrosion

open MeasureTheory Complex Filter Topology

/-! ============================================================
  §50.1  Strict-uniform theorem with named `mcmullen_K_star`
============================================================ -/

section ClosedFormTheoremC

/-- **★★★★★★★ Strict-uniform complexity with closed-form
    `mcmullen_K_star` on ℂ at p = 2.**

    For every `x ∈ ℂ \ {0, 1}` (with `α ≠ 0, ±1`, `α² = 1/x`),
    every `R, δ, r > 0`, and the §46 non-degeneracy `θ < 1`:
    there exists `N : ℕ` such that:

      (i)  The bad set has Lebesgue volume `< δ`.

      (ii) For every `z₀` outside the bad set,
           `∃ s : Fin 2, ∀ k ≥ mcmullen_K_star α r N,
              ‖σ^k z₀ − γ_s‖ < r`.

    Unlike §49's `∃ K_star, …` form, here `K_star` is the **named
    computable function** `mcmullen_K_star α r N`. Given numerical
    `(α, r, N)`, it evaluates to a concrete `ℕ`. -/
theorem mcmullen_p2_complex_closed_form
    (x : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (R : ℝ) (δ : ℝ) (hδ_pos : 0 < δ) (r : ℝ) (hr_pos : 0 < r)
    (hθ_lt_one : min (‖((1 - α) / (1 + α)) ^ 2‖ / 2)
                     (r * ‖((1 - α) / (1 + α)) ^ 2‖ / (4 * ‖α‖)) < 1) :
    ∃ N : ℕ,
      volume ((slowSet α (1 / ((N : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R)
              ∪ alphaHitSet x α
              ∪ complex_bad_set x α) < ENNReal.ofReal δ ∧
      ∀ z₀ : ℂ,
        z₀ ∉ slowSet α (1 / ((N : ℝ) + 1)) →
        z₀ ∉ alphaHitSet x α →
        z₀ ∉ complex_bad_set x α →
        0 < ‖v_p2_C α z₀‖ →
        ∃ s : Fin 2, ∀ k ≥ mcmullen_K_star α r N,
          ‖(steffensen_step_C x 2)^[k] z₀ - cycAnchor α 2 s‖ < r :=
  mcmullen_p2_complex_measure_uniform_strict x hx_ne hx_ne_one α hα_ne_zero
    hα_ne_neg_one hα_ne_one hα_pow R δ hδ_pos r hr_pos hθ_lt_one

end ClosedFormTheoremC

end Pandrosion
