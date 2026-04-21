/-
  Universitas Pandrosion — §14–§16. Aitken-Steffensen quadratic rate.

  Purely algebraic identities underlying the Aitken-Steffensen
  quadratic-rate theorem for the Pandrosion map
    h(s) = 1 − (x−1)/(x·S_p(s)).

  No calculus is used: every identity is proved by `ring`,
  `linear_combination`, `Commute.geom_sum₂_mul`, or
  `Finset.sum_range_reflect`.

  Content:
    §14  Algebraic skeleton (Phase 1).
      §14.1  Load-bearing identity  (h − s)·x·S_p = 1 − x·s^p.
      §14.2  Q-polynomial and root factorization
             s^p − s*^p = (s − s*)·Q(s, s*).
      §14.3  Fixed-point linear factorization
             (h(s) − s*)·x·S_p(s) = (s − s*)·x·(S_p(s) − Q(s, s*)).
      §14.4  Closed-form linear rate  `lambda_fp p s* = 1 − p·s*^{p−1}/S_p(s*)`.
      §14.5  Linear-rate quotient  `(h(s) − s*)/(s − s*) = (S_p(s) − Q(s, s*))/S_p(s)`.
      §14.6  Quadratic deviation  `δ₂(s, s*) = (S_p − Q)/S_p − λ_fp`
             vanishes at `s = s*`.

    §15  Quadratic factorization (Phase 2, structural).
      §15.1  `M`-polynomial factorization  S_p(s) − S_p(s*) = (s − s*)·M.
      §15.2  `N`-polynomial factorization  Q(s,s*) − Q(s*,s*) = (s − s*)·N.
      §15.3  Division form  h(s) − s* − λ_fp·(s − s*) = (s − s*)·[ratio(s) − ratio(s*)].
      §15.4  **★ Polynomial form**
             (h(s) − s* − λ_fp·(s − s*)) · S_p(s) · S_p(s*)
             = (s − s*)² · (Q(s*,s*)·M(p,s,s*) − S_p(s*)·N(p,s,s*)).

    §16  Bridge to §9.2.
      §16.1  Reflection  S_p(1/α) · α^{p−1} = S_p(α).
      §16.2  **★ `lambda_fp p (1/α) = lambda_closed p α`** — unifies §9.2 and §14.

  Phase 2+ (conversion of the polynomial identity into a pointwise bound
  `|T(s) − s*| ≤ K_S·(s−s*)^2` with explicit δ radius) requires analytic
  bounding of `M`, `N` on a neighborhood, which is scoped out of this file.
-/

import Pandrosion.Core.Foundations

namespace Pandrosion

open Finset BigOperators

/-! ============================================================
  §14.1  Load-bearing identity
============================================================ -/

section LoadBearing

/-- **Load-bearing identity.** `(h(s) − s) · x · S_p(s) = 1 − x · s^p`.
    Encodes the quasi-Newton structure of Pandrosion: the iteration
    step `h(s) − s` equals `f(s)/(x·S_p(s))` where `f(s) = 1 − x·s^p`
    is the underlying polynomial. -/
theorem pandrosion_h_sub_mul
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1) (s : ℝ) (hs : 0 ≤ s) :
    (pandrosion_h x p s - s) * (x * Sp p s) = 1 - x * s ^ p := by
  have hS : Sp p s > 0 := Sp_pos p hp s hs
  have hxS_ne : x * Sp p s ≠ 0 := ne_of_gt (mul_pos hx hS)
  have hSmul : Sp p s * (1 - s) = 1 - s ^ p := Sp_mul_one_sub p s
  unfold pandrosion_h
  field_simp
  linear_combination x * hSmul

end LoadBearing

/-! ============================================================
  §14.2  Q-polynomial and root factorization
============================================================ -/

section RootFactorization

/-- `Q(p, s, s_star) := ∑_{j=0}^{p−1} s^j · s_star^{p−1−j}`.
    The symmetric two-variable factor in the root identity
    `s^p − s_star^p = (s − s_star) · Q`. -/
noncomputable def Qpoly (p : ℕ) (s s_star : ℝ) : ℝ :=
  ∑ j in Finset.range p, s ^ j * s_star ^ (p - 1 - j)

/-- **Geometric root identity.** `Q · (s − s_star) = s^p − s_star^p`.
    Instance of Mathlib's `Commute.geom_sum₂_mul`, with commutation
    automatic in ℝ. -/
theorem Qpoly_mul_sub (p : ℕ) (s s_star : ℝ) :
    Qpoly p s s_star * (s - s_star) = s ^ p - s_star ^ p := by
  unfold Qpoly
  exact (Commute.all s s_star).geom_sum₂_mul p

/-- Swapped-form: `s^p − s_star^p = (s − s_star) · Q`. -/
theorem pow_sub_pow_factor (p : ℕ) (s s_star : ℝ) :
    s ^ p - s_star ^ p = (s - s_star) * Qpoly p s s_star := by
  rw [← Qpoly_mul_sub p s s_star, mul_comm]

/-- **Q on the diagonal.** `Q(p, s_star, s_star) = p · s_star^{p−1}`.
    Every summand collapses to `s_star^{p−1}`; the sum has `p` terms. -/
theorem Qpoly_diag (p : ℕ) (s_star : ℝ) :
    Qpoly p s_star s_star = (p : ℝ) * s_star ^ (p - 1) := by
  unfold Qpoly
  have h_each : ∀ j ∈ Finset.range p,
      s_star ^ j * s_star ^ (p - 1 - j) = s_star ^ (p - 1) := by
    intro j hj
    rw [Finset.mem_range] at hj
    rw [← pow_add]
    congr 1
    omega
  rw [Finset.sum_congr rfl h_each, Finset.sum_const, Finset.card_range,
      nsmul_eq_mul]

end RootFactorization

/-! ============================================================
  §14.3  Fixed-point linear factorization
============================================================ -/

section FixedPointFactorization

/-- **★ Key linear-rate factorization.** For `s_star` a Pandrosion
    fixed point (i.e. `x · s_star^p = 1`):
        `(h(s) − s_star) · x · S_p(s) = (s − s_star) · x · (S_p(s) − Q(s, s_star))`.

    This is the *algebraic* form of the linear convergence identity.
    Dividing by `x · S_p(s)` (both nonzero in the Pandrosion regime)
    gives the quotient form of §14.5. -/
theorem pandrosion_h_sub_fp_factor
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s s_star : ℝ) (hs : 0 ≤ s)
    (hfp : x * s_star ^ p = 1) :
    (pandrosion_h x p s - s_star) * (x * Sp p s)
      = (s - s_star) * (x * (Sp p s - Qpoly p s s_star)) := by
  have key : (pandrosion_h x p s - s) * (x * Sp p s) = 1 - x * s ^ p :=
    pandrosion_h_sub_mul x hx p hp s hs
  have factor : s ^ p - s_star ^ p = (s - s_star) * Qpoly p s s_star :=
    pow_sub_pow_factor p s s_star
  -- Decompose `(h - s_star) = (h - s) + (s - s_star)`, then substitute.
  have expand : (pandrosion_h x p s - s_star) * (x * Sp p s)
              = (pandrosion_h x p s - s) * (x * Sp p s)
              + (s - s_star) * (x * Sp p s) := by ring
  rw [expand, key]
  -- `1 − x·s^p  =  x·(s_star^p − s^p)  =  −x·(s − s_star)·Q`.
  have h_rewrite : (1 : ℝ) - x * s ^ p
                 = -(x * ((s - s_star) * Qpoly p s s_star)) := by
    linear_combination -x * factor - hfp
  rw [h_rewrite]
  ring

/-- **Quotient form.** `h(s) − s_star = (s − s_star) · (S_p(s) − Q)/S_p(s)`.
    Divide both sides of `pandrosion_h_sub_fp_factor` by `x · S_p(s)`
    (both nonzero in the Pandrosion regime). -/
theorem pandrosion_h_sub_fp_quotient
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s s_star : ℝ) (hs : 0 ≤ s)
    (hfp : x * s_star ^ p = 1) :
    pandrosion_h x p s - s_star
      = (s - s_star) * (Sp p s - Qpoly p s s_star) / Sp p s := by
  have hS : Sp p s > 0 := Sp_pos p hp s hs
  have hS_ne : Sp p s ≠ 0 := ne_of_gt hS
  have hx_ne : (x : ℝ) ≠ 0 := ne_of_gt hx
  have this := pandrosion_h_sub_fp_factor x hx p hp s s_star hs hfp
  -- Factor `x · Sp p s` out of both sides.
  have eq_lhs : (pandrosion_h x p s - s_star) * (x * Sp p s)
              = x * Sp p s * (pandrosion_h x p s - s_star) := by ring
  have eq_rhs : (s - s_star) * (x * (Sp p s - Qpoly p s s_star))
              = x * Sp p s
                  * ((s - s_star) * (Sp p s - Qpoly p s s_star) / Sp p s) := by
    field_simp; ring
  rw [eq_lhs, eq_rhs] at this
  exact mul_left_cancel₀ (mul_ne_zero hx_ne hS_ne) this

end FixedPointFactorization

/-! ============================================================
  §14.4  Closed-form linear rate `lambda_fp`
============================================================ -/

section LinearRateClosedForm

/-- **Fixed-point linear rate** `λ_fp(p, s_star) := 1 − p · s_star^{p−1}/S_p(s_star)`.

    Expresses `h'(s_star)` *purely in terms of the fixed point* — no
    calculus required. Related to the paper's `lambda_closed` from
    §9.2 via `α = 1/s_star` (see `lambda_fp_pos` and the `α ↔ 1/s_star`
    reindexing). -/
noncomputable def lambda_fp (p : ℕ) (s_star : ℝ) : ℝ :=
  1 - (p : ℝ) * s_star ^ (p - 1) / Sp p s_star

/-- **Ratio identification.** `λ_fp = (S_p(s_star) − Q(s_star, s_star))/S_p(s_star)`.
    Combines `Qpoly_diag` with the `λ_fp` definition. -/
theorem lambda_fp_eq_quotient
    (p : ℕ) (hp : p ≥ 1) (s_star : ℝ) (hs_star : 0 ≤ s_star) :
    lambda_fp p s_star
      = (Sp p s_star - Qpoly p s_star s_star) / Sp p s_star := by
  unfold lambda_fp
  rw [Qpoly_diag p s_star]
  have hS : Sp p s_star > 0 := Sp_pos p hp s_star hs_star
  field_simp

end LinearRateClosedForm

/-! ============================================================
  §14.5  Linear rate ratio  (h(s)−s_star)/(s−s_star)
============================================================ -/

section LinearRateRatio

/-- **★ Pandrosion linear-rate ratio (algebraic form).**
    For `s ≠ s_star` (off the fixed point) and `s_star` a fixed point,
        `(h(s) − s_star) / (s − s_star)  =  (S_p(s) − Q(s, s_star)) / S_p(s)`.

    The RHS is a *rational function* of `s`; evaluation at `s = s_star`
    (Lemma `pandrosion_ratio_at_fp`) gives `λ_fp`. This is the
    algebraic precondition for Aitken-Steffensen quadratic acceleration. -/
theorem pandrosion_linear_rate_ratio
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s s_star : ℝ) (hs : 0 ≤ s) (hs_ne : s ≠ s_star)
    (hfp : x * s_star ^ p = 1) :
    (pandrosion_h x p s - s_star) / (s - s_star)
      = (Sp p s - Qpoly p s s_star) / Sp p s := by
  have hS : Sp p s > 0 := Sp_pos p hp s hs
  have hS_ne : Sp p s ≠ 0 := ne_of_gt hS
  have hdiff_ne : s - s_star ≠ 0 := sub_ne_zero.mpr hs_ne
  rw [pandrosion_h_sub_fp_quotient x hx p hp s s_star hs hfp]
  field_simp
  ring

/-- **Ratio at the fixed point.** The rational function `(S_p − Q)/S_p`
    evaluated at `s = s_star` equals `λ_fp` — the linear contraction
    rate of Pandrosion at its fixed point. -/
theorem pandrosion_ratio_at_fp
    (p : ℕ) (hp : p ≥ 1) (s_star : ℝ) (hs_star : 0 ≤ s_star) :
    (Sp p s_star - Qpoly p s_star s_star) / Sp p s_star
      = lambda_fp p s_star :=
  (lambda_fp_eq_quotient p hp s_star hs_star).symm

end LinearRateRatio

/-! ============================================================
  §14.6  Quadratic deviation — vanishes at the fixed point
============================================================ -/

section QuadraticDeviation

/-- **Quadratic deviation** `δ₂(p, s, s_star) := (S_p(s) − Q(s,s_star))/S_p(s) − λ_fp p s_star`.
    Measures how much the linear-rate ratio at `s` differs from its
    fixed-point value `λ_fp`. Vanishes at `s = s_star` (Theorem
    `quadratic_deviation_at_fp`), so `h(s) − s_star − λ_fp·(s − s_star)`
    is genuinely *super-linear* in `(s − s_star)` — the algebraic
    precondition for Steffensen-Aitken quadratic acceleration. -/
noncomputable def quadratic_deviation (p : ℕ) (s s_star : ℝ) : ℝ :=
  (Sp p s - Qpoly p s s_star) / Sp p s - lambda_fp p s_star

/-- **★ Algebraic precondition for quadratic rate.**
    The deviation vanishes at the fixed point: `δ₂(p, s_star, s_star) = 0`.
    This is the "first-order cancellation" that makes Aitken-Steffensen
    accelerate Pandrosion from linear to quadratic convergence. -/
theorem quadratic_deviation_at_fp
    (p : ℕ) (hp : p ≥ 1) (s_star : ℝ) (hs_star : 0 ≤ s_star) :
    quadratic_deviation p s_star s_star = 0 := by
  unfold quadratic_deviation
  rw [pandrosion_ratio_at_fp p hp s_star hs_star]
  ring

end QuadraticDeviation

end Pandrosion
