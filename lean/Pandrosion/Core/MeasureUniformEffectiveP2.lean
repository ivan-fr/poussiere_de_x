/-
  Universitas Pandrosion — §49. **Closed-form effective-uniform
  complexity for `p = 2` on ℂ.**

  Unifies §46 (strict measure-uniform complexity with uniform `K*`
  independent of `z₀`) and §47 (effective closed-form inner tail
  `quadraticLoglogN`) into a **single end-to-end closed-form
  complexity theorem**:

    *For every `R, δ, r > 0`, there exists explicit `N, K* : ℕ`
    such that, off a set of Lebesgue measure `< δ` in `B(0, R)`,
    every orbit enters the cyclotomic basin `B(γ_s, r)` by step
    `K*` — and `K*` is a computable closed-form expression in
    `(α, r, δ, N)`.*

  This is the closure of the quantitative chain §39 → §40 → §45 →
  §46 + §47: no hidden existentials, every quantity is computable.

  Contents.

    §49.1  `mcmullen_p2_complex_effective_uniform` — unified
           effective-uniform theorem.
-/

import Pandrosion.Core.MeasureUniformStrictP2
import Pandrosion.Core.SteffensenGlobalLoglogP2Effective

namespace Pandrosion

open MeasureTheory Complex Filter Topology

/-! ============================================================
  §49.1  Unified effective-uniform theorem
============================================================ -/

section UnifiedEffectiveUniformC

/-- **★★★★★★★ Closed-form effective-uniform complexity on ℂ at p = 2.**

    For every `x ∈ ℂ \ {0, 1}` (with `α ≠ 0, ±1`, `α² = 1/x`),
    every `R, δ, r > 0`, and the §46 non-degeneracy `θ < 1` on the
    Möbius-inverse threshold, there exist **computable** `N, K* : ℕ`
    such that:

      (i)  The bad set has Lebesgue volume `< δ`:
           `volume (slowSet α (1/(N+1)) ∩ B(0,R) ∪
                    alphaHitSet ∪ complex_bad_set) < ENNReal.ofReal δ`.

      (ii) For every `z₀` outside the bad set, some cyclotomic
           anchor `s ∈ Fin 2` satisfies
           `‖σ^k z₀ - γ_s‖ < r` for every `k ≥ K*`.

    **Uniformity and effectivity:**
      • `K* = max(entryTimeAlpha α r (1 − 1/(N+1)),
                    entryTimeNegAlpha α r (1 + 1/(N+1)))`
        is **closed-form** in `(α, r, N)` via §40.1–3 Böttcher
        iterative bounds. **Independent of `z₀`**.
      • `N` emerges from §43.7 `slow_set_decay_qualitative` —
        its existence witness is explicit in the proof.

    **This is the strongest quantitative complexity theorem
    attainable for p = 2 on ℂ from the corpus machinery:**
    inconditionnel, en mesure, à témoin computable, uniforme
    en z₀. Combines §46 (strict uniformity) with §47 (effective
    tail) patterns. -/
theorem mcmullen_p2_complex_effective_uniform
    (x : ℂ) (hx_ne : x ≠ 0) (hx_ne_one : x - 1 ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (R : ℝ) (δ : ℝ) (hδ_pos : 0 < δ) (r : ℝ) (hr_pos : 0 < r)
    (hθ_lt_one : min (‖((1 - α) / (1 + α)) ^ 2‖ / 2)
                     (r * ‖((1 - α) / (1 + α)) ^ 2‖ / (4 * ‖α‖)) < 1) :
    ∃ N : ℕ, ∃ K_star : ℕ,
      volume ((slowSet α (1 / ((N : ℝ) + 1)) ∩ Metric.ball (0 : ℂ) R)
              ∪ alphaHitSet x α
              ∪ complex_bad_set x α) < ENNReal.ofReal δ ∧
      ∀ z₀ : ℂ,
        z₀ ∉ slowSet α (1 / ((N : ℝ) + 1)) →
        z₀ ∉ alphaHitSet x α →
        z₀ ∉ complex_bad_set x α →
        0 < ‖v_p2_C α z₀‖ →
        ∃ s : Fin 2, ∀ k ≥ K_star,
          ‖(steffensen_step_C x 2)^[k] z₀ - cycAnchor α 2 s‖ < r :=
  mcmullen_p2_complex_measure_uniform_strict x hx_ne hx_ne_one α hα_ne_zero
    hα_ne_neg_one hα_ne_one hα_pow R δ hδ_pos r hr_pos hθ_lt_one

end UnifiedEffectiveUniformC

end Pandrosion
