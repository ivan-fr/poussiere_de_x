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

  The pointwise bound `|T(s) − s*| ≤ K_S·(s−s*)²` with explicit δ radius
  reduces to §15.4 (polynomial form) plus neighborhood bounds on `M`, `N`;
  the abstract log-log complexity corollary lives in `QuadraticComplexity.lean`.
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

/-! ============================================================
  §15. Quadratic factorization (Phase 2, structural)
============================================================ -/

/-! ## §15.1  M-polynomial factorization  S_p(s) − S_p(s*) = (s − s*) · M -/

section MFactorization

/-- **M-polynomial** `M(p, s, s_star) := Σ_{k=0}^{p−1} Q(k, s, s_star)`.
    Coefficient polynomial in the telescoping factorization of
    `S_p(s) − S_p(s_star)`. Since `Q(0, s, s_star) = 0` (empty sum),
    the first nonzero summand occurs at `k = 1`. -/
noncomputable def Mpoly (p : ℕ) (s s_star : ℝ) : ℝ :=
  ∑ k in Finset.range p, Qpoly k s s_star

/-- **M-factorization identity.**
    `S_p(s) − S_p(s_star) = (s − s_star) · M(p, s, s_star)`.

    Proved by telescoping the elementary factorization
    `s^k − s_star^k = (s − s_star) · Q(k, s, s_star)` (`pow_sub_pow_factor`)
    across `k ∈ [0, p)`. -/
theorem Sp_sub_Sp_factor (p : ℕ) (s s_star : ℝ) :
    Sp p s - Sp p s_star = (s - s_star) * Mpoly p s s_star := by
  unfold Sp Mpoly
  induction p with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ (fun k => s ^ k),
        Finset.sum_range_succ (fun k => s_star ^ k),
        Finset.sum_range_succ (fun k => Qpoly k s s_star)]
    have hn := pow_sub_pow_factor n s s_star
    linear_combination ih + hn

end MFactorization

/-! ## §15.2  N-polynomial factorization
    `Q(p, s, s*) − Q(p, s*, s*) = (s − s*) · N` -/

section NFactorization

/-- **N-polynomial**
    `N(p, s, s_star) := Σ_{j=0}^{p−1} Q(j, s, s_star) · s_star^{p−1−j}`.
    Coefficient polynomial in the factorization
    `Q(p, s, s_star) − Q(p, s_star, s_star) = (s − s_star) · N`. -/
noncomputable def Npoly (p : ℕ) (s s_star : ℝ) : ℝ :=
  ∑ j in Finset.range p, Qpoly j s s_star * s_star ^ (p - 1 - j)

/-- **N-factorization identity.**
    `Q(p, s, s_star) − Q(p, s_star, s_star) = (s − s_star) · N(p, s, s_star)`.

    Factor `s^j − s_star^j = (s − s_star) · Q(j, s, s_star)` inside each
    summand of `Qpoly p s s_star − Qpoly p s_star s_star`. -/
theorem Qpoly_sub_Qpoly_factor (p : ℕ) (s s_star : ℝ) :
    Qpoly p s s_star - Qpoly p s_star s_star = (s - s_star) * Npoly p s s_star := by
  unfold Qpoly Npoly
  rw [← Finset.sum_sub_distrib, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  have hj := pow_sub_pow_factor j s s_star
  linear_combination s_star ^ (p - 1 - j) * hj

end NFactorization

/-! ## §15.3  Lambda product form + h factorization without `x` -/

section LambdaAndFactor

/-- **λ_fp product form.** `λ_fp(p, s_star) · S_p(s_star) = S_p(s_star) − Q(s_star, s_star)`.
    Multiplies the quotient form of §14.4 through by the denominator. -/
theorem lambda_fp_mul_Sp (p : ℕ) (hp : p ≥ 1) (s_star : ℝ) (hs_star : 0 ≤ s_star) :
    lambda_fp p s_star * Sp p s_star = Sp p s_star - Qpoly p s_star s_star := by
  have hS : Sp p s_star > 0 := Sp_pos p hp s_star hs_star
  unfold lambda_fp
  rw [Qpoly_diag p s_star]
  field_simp

/-- **Pandrosion fixed-point factorization, x-free form.**
    `(h(s) − s_star) · S_p(s) = (s − s_star) · (S_p(s) − Q(s, s_star))`.

    Divides `pandrosion_h_sub_fp_factor` through by `x ≠ 0`. -/
theorem pandrosion_h_sub_fp_factor_no_x
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s s_star : ℝ) (hs : 0 ≤ s)
    (hfp : x * s_star ^ p = 1) :
    (pandrosion_h x p s - s_star) * Sp p s
      = (s - s_star) * (Sp p s - Qpoly p s s_star) := by
  have hx_ne : x ≠ 0 := ne_of_gt hx
  have h := pandrosion_h_sub_fp_factor x hx p hp s s_star hs hfp
  have hx_mul :
      x * ((pandrosion_h x p s - s_star) * Sp p s)
        = x * ((s - s_star) * (Sp p s - Qpoly p s s_star)) := by
    linear_combination h
  exact mul_left_cancel₀ hx_ne hx_mul

end LambdaAndFactor

/-! ## §15.4  ★ Polynomial form of the quadratic identity -/

section QuadraticPolynomial

/-- **★ Polynomial form of the quadratic residue identity.**
    For `s_star` a Pandrosion fixed point (`x · s_star^p = 1`):

    `(h(s) − s_star − λ_fp · (s − s_star)) · S_p(s) · S_p(s_star)`
    `= (s − s_star)² · (Q(s_star, s_star) · M(p, s, s_star) − S_p(s_star) · N(p, s, s_star))`.

    The LHS is the super-linear residue `h(s) − s_star − λ_fp·(s − s_star)`
    multiplied by the positive quantity `S_p(s) · S_p(s_star)`. The RHS
    factors as `(s − s_star)²` times a polynomial in `s, s_star`. This is
    the *algebraic proof* that `h(s) − s_star − λ_fp·(s − s_star) = O((s − s_star)²)`
    — the first-order cancellation underlying Aitken–Steffensen quadratic
    acceleration. -/
theorem pandrosion_h_quadratic_polynomial
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s s_star : ℝ) (hs : 0 ≤ s) (hs_star : 0 ≤ s_star)
    (hfp : x * s_star ^ p = 1) :
    (pandrosion_h x p s - s_star - lambda_fp p s_star * (s - s_star))
      * Sp p s * Sp p s_star
      = (s - s_star) ^ 2 *
        (Qpoly p s_star s_star * Mpoly p s s_star
         - Sp p s_star * Npoly p s s_star) := by
  have h_factor' :=
    pandrosion_h_sub_fp_factor_no_x x hx p hp s s_star hs hfp
  have h_lambda := lambda_fp_mul_Sp p hp s_star hs_star
  have hM := Sp_sub_Sp_factor p s s_star
  have hN := Qpoly_sub_Qpoly_factor p s s_star
  linear_combination
    Sp p s_star * h_factor'
      - (s - s_star) * Sp p s * h_lambda
      + (s - s_star) * Qpoly p s_star s_star * hM
      - (s - s_star) * Sp p s_star * hN

/-- **Quotient form of the quadratic residue.**
    `h(s) − s_star − λ_fp·(s − s_star) = (s − s_star)² · R(p, s, s_star)`
    where `R := (Q(s_star,s_star)·M − S_p(s_star)·N) / (S_p(s)·S_p(s_star))`
    is a rational function of `s, s_star`. Dividing the polynomial form
    through by the strictly positive product `S_p(s) · S_p(s_star)`. -/
theorem pandrosion_h_quadratic_quotient
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s s_star : ℝ) (hs : 0 ≤ s) (hs_star : 0 ≤ s_star)
    (hfp : x * s_star ^ p = 1) :
    pandrosion_h x p s - s_star - lambda_fp p s_star * (s - s_star)
      = (s - s_star) ^ 2 *
        ((Qpoly p s_star s_star * Mpoly p s s_star
          - Sp p s_star * Npoly p s s_star) / (Sp p s * Sp p s_star)) := by
  have hS := Sp_pos p hp s hs
  have hS_star := Sp_pos p hp s_star hs_star
  have h_poly :=
    pandrosion_h_quadratic_polynomial x hx p hp s s_star hs hs_star hfp
  have hS_ne : Sp p s ≠ 0 := ne_of_gt hS
  have hS_star_ne : Sp p s_star ≠ 0 := ne_of_gt hS_star
  rw [← mul_div_assoc, eq_div_iff (mul_ne_zero hS_ne hS_star_ne)]
  linear_combination h_poly

/-- **Pointwise quadratic bound for `h`** (conditional).
    Given any upper bound `C` on `|Q(s_star,s_star)·M − S_p(s_star)·N| / (S_p(s)·S_p(s_star))`
    over a neighborhood, the super-linear residue is quadratically bounded
    there. This is the load-bearing step from the polynomial identity to
    the analytic rate; the constant `C` is obtained by explicit bounding
    of the polynomial sums `M`, `N` on `[s_star − δ, s_star + δ]`. -/
theorem pandrosion_h_quadratic_bound_of
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s s_star : ℝ) (hs : 0 ≤ s) (hs_star : 0 ≤ s_star)
    (hfp : x * s_star ^ p = 1)
    (C : ℝ)
    (hC : |(Qpoly p s_star s_star * Mpoly p s s_star
            - Sp p s_star * Npoly p s s_star) / (Sp p s * Sp p s_star)| ≤ C) :
    |pandrosion_h x p s - s_star - lambda_fp p s_star * (s - s_star)|
      ≤ C * (s - s_star) ^ 2 := by
  rw [pandrosion_h_quadratic_quotient x hx p hp s s_star hs hs_star hfp]
  rw [abs_mul, abs_of_nonneg (sq_nonneg (s - s_star))]
  have h_nn : 0 ≤ (s - s_star) ^ 2 := sq_nonneg _
  calc
    (s - s_star) ^ 2 * |((Qpoly p s_star s_star * Mpoly p s s_star
                          - Sp p s_star * Npoly p s s_star) / (Sp p s * Sp p s_star))|
        ≤ (s - s_star) ^ 2 * C :=
          mul_le_mul_of_nonneg_left hC h_nn
    _ = C * (s - s_star) ^ 2 := by ring

end QuadraticPolynomial

/-! ============================================================
  §16. Bridge to §9.2 — unifies `lambda_fp` with `lambda_closed`
============================================================ -/

/-! ## §16.1  Reflection identity  S_p(1/α) · α^{p−1} = S_p(α) -/

section Reflection

/-- **Reflection identity.** `S_p(1/α) · α^{p−1} = S_p(α)` for `α ≠ 0`.
    Equivalent to `Σ_k α^{−k} · α^{p−1} = Σ_k α^{p−1−k}` followed by the
    reindexing `j = p − 1 − k` (`Finset.sum_range_reflect`). -/
theorem Sp_reflect (p : ℕ) (α : ℝ) (hα : α ≠ 0) :
    Sp p (1 / α) * α ^ (p - 1) = Sp p α := by
  unfold Sp
  rw [Finset.sum_mul, ← Finset.sum_range_reflect (fun j => α ^ j) p]
  apply Finset.sum_congr rfl
  intro k hk
  rw [Finset.mem_range] at hk
  have hk' : k ≤ p - 1 := by omega
  rw [one_div, inv_pow, mul_comm ((α ^ k)⁻¹) (α ^ (p - 1))]
  exact (pow_sub₀ α hα hk').symm

end Reflection

/-! ## §16.2  ★ `lambda_fp p (1/α) = lambda_closed p α` -/

section LambdaBridge

/-- **★ Bridge theorem.** `λ_fp(p, 1/α) = λ_closed(p, α)` — the fixed-point
    rate evaluated at `1/α` equals the paper's closed-form rate
    `1 − p/S_p(α)`. Unifies the algebraic (§14.4) and analytic (§9.2)
    formulations of the Pandrosion contraction rate. -/
theorem lambda_fp_eq_lambda_closed
    (p : ℕ) (hp : p ≥ 1) (α : ℝ) (hα_pos : α > 0) :
    lambda_fp p (1 / α) = lambda_closed p α := by
  have hα_ne : α ≠ 0 := ne_of_gt hα_pos
  have h_inv_nonneg : (0 : ℝ) ≤ 1 / α := le_of_lt (one_div_pos.mpr hα_pos)
  have h_refl := Sp_reflect p α hα_ne
  have hSp_α : Sp p α > 0 := Sp_pos p hp α (le_of_lt hα_pos)
  have hSp_inv : Sp p (1 / α) > 0 := Sp_pos p hp (1 / α) h_inv_nonneg
  have hα_pow : α ^ (p - 1) > 0 := pow_pos hα_pos _
  unfold lambda_fp lambda_closed
  have hSp_inv_ne : Sp p (1 / α) ≠ 0 := ne_of_gt hSp_inv
  have hSp_α_ne : Sp p α ≠ 0 := ne_of_gt hSp_α
  have hα_pow_ne : α ^ (p - 1) ≠ 0 := ne_of_gt hα_pow
  have h_inv : (1 / α) ^ (p - 1) = (α ^ (p - 1))⁻¹ := by
    rw [one_div, inv_pow]
  have key : (1 / α) ^ (p - 1) * Sp p α = Sp p (1 / α) := by
    rw [h_inv, ← h_refl]
    rw [mul_comm (Sp p (1 / α)) (α ^ (p - 1)), ← mul_assoc]
    rw [inv_mul_cancel hα_pow_ne, one_mul]
  have h_eq :
      (p : ℝ) * (1 / α) ^ (p - 1) / Sp p (1 / α) = (p : ℝ) / Sp p α := by
    rw [div_eq_div_iff hSp_inv_ne hSp_α_ne]
    linear_combination (p : ℝ) * key
  rw [h_eq]

end LambdaBridge

end Pandrosion
