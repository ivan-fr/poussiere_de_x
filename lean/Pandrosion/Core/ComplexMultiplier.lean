/-
  Universitas Pandrosion — §26. Complex multiplier of the Pandrosion
  map `h` at a fixed point (Julia-dynamics data).

  The *multiplier* of a holomorphic map at a fixed point is its
  complex derivative there. For the Pandrosion iterator
      `h(z) = 1 − (x − 1) / (x · S_p(z))`
  this multiplier governs whether α is attracting (`|h′(α)| < 1`),
  repelling, parabolic, or super-attracting (`h′(α) = 0`).

  This file bridges the algebraic `lambda_fp` from §14.4 to
  complex calculus: we exhibit `h′(α)` via `HasDerivAt` and
  identify it with the closed-form paper multiplier
      `λ_C(p, α) = 1 − p · α^{p−1} / S_p(α)`
  at any Pandrosion fixed point α (i.e. `α^p = 1/x, S_p(α) ≠ 0`).

  **Scope.** This is the *algebraic heart* of Option M in the
  §24 roadmap — the rigorous complex-calculus multiplier of the
  linear iterator `h`. The Steffensen-accelerator derivative-zero
  claim (the full Julia super-attracting property of
  `steffensen_step_C`) is a *separate* theorem requiring handling
  the removable 0/0 in `(h − id)² / D` at α; it is not proved here.

  Content.

    §26.1  `SpDerivC`, `Sp_C_hasDerivAt` — formal derivative of
           the geometric sum `S_p`, with `HasDerivAt` statement.

    §26.2  `SpDerivC_mul_one_sub` — Leibniz-algebraic identity
           `S_p′(α) · (1 − α) = S_p(α) − p · α^{p−1}`, proved by
           differentiating `S_p · (1 − ·) = 1 − (·)^p`.

    §26.3  `hMultC`, `pandrosion_h_C_hasDerivAt` — complex
           derivative of `h` at any point where `x · S_p(α) ≠ 0`.

    §26.4  `lambdaClosedC`, `hMultC_eq_lambdaClosedC_at_fp` — at
           a Pandrosion fixed point, `h′(α)` equals the paper's
           closed-form multiplier.

    ★ §26.5  `pandrosion_h_C_hasDerivAt_at_fp` — combined crown:
           `HasDerivAt h = λ_C` at any Pandrosion fixed point.
-/

import Pandrosion.Core.Foundations
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Inv

namespace Pandrosion

open Complex BigOperators Finset

/-! ============================================================
  §26.1  Complex derivative of `Sp_C`
============================================================ -/

section SpDerivC

/-- **Formal complex derivative of `S_p(z) = ∑_{k<p} z^k`** at α.
    Termwise differentiation: `∑_{k<p} k · α^{k−1}`. -/
noncomputable def SpDerivC (p : ℕ) (α : ℂ) : ℂ :=
  ∑ k in Finset.range p, (k : ℂ) * α ^ (k - 1)

/-- **Complex derivative of `Sp_C`** via `HasDerivAt`. -/
theorem Sp_C_hasDerivAt (p : ℕ) (α : ℂ) :
    HasDerivAt (Sp_C p) (SpDerivC p α) α := by
  show HasDerivAt (fun z : ℂ => ∑ k in Finset.range p, z ^ k)
                  (∑ k in Finset.range p, (k : ℂ) * α ^ (k - 1)) α
  exact HasDerivAt.sum (fun k _ => hasDerivAt_pow k α)

end SpDerivC

/-! ============================================================
  §26.2  Leibniz-algebraic identity for `SpDerivC`
============================================================ -/

section SpDerivIdentity

/-- **★ Leibniz identity for `S_p`.**
    Follows by differentiating the polynomial identity
    `S_p(z) · (1 − z) = 1 − z^p` at α: term-by-term rules give
        `S_p′(α) · (1 − α) − S_p(α) = −p · α^{p−1}`,
    rearranged as
        `S_p′(α) · (1 − α) = S_p(α) − p · α^{p−1}`.

    This identity is the algebraic bridge from `hMultC` (the raw
    multiplier, a ratio of derivative expressions) to the paper's
    closed-form `λ_C = 1 − p · α^{p−1} / S_p(α)`. -/
theorem SpDerivC_mul_one_sub (p : ℕ) (α : ℂ) :
    SpDerivC p α * (1 - α) = Sp_C p α - (p : ℂ) * α ^ (p - 1) := by
  -- Derivative of `z ↦ 1 − z` at α.
  have h_one_sub : HasDerivAt (fun z : ℂ => (1 : ℂ) - z) (-1 : ℂ) α := by
    have := (hasDerivAt_id α).const_sub (1 : ℂ)
    simpa using this
  -- Derivative of `z ↦ S_p(z) · (1 − z)` at α via product rule.
  have h_mul := (Sp_C_hasDerivAt p α).mul h_one_sub
  -- Derivative of `z ↦ 1 − z^p` at α.
  have h_pow : HasDerivAt (fun z : ℂ => z ^ p) ((p : ℂ) * α ^ (p - 1)) α :=
    hasDerivAt_pow p α
  have h_rhs := h_pow.const_sub (1 : ℂ)
  -- Functional equality `S_p · (1 − ·) = 1 − (·)^p` as a literal equality of functions.
  have h_eq : (fun z : ℂ => Sp_C p z * (1 - z)) = (fun z : ℂ => 1 - z ^ p) := by
    funext z; exact Sp_C_mul_one_sub p z
  rw [h_eq] at h_mul
  -- Uniqueness of derivative at α.
  have h_unique := h_mul.unique h_rhs
  -- `h_unique : SpDerivC p α · (1 − α) + Sp_C p α · (−1) = −(p · α^{p−1})`
  linear_combination h_unique

end SpDerivIdentity

/-! ============================================================
  §26.3  General-point complex derivative of `pandrosion_h_C`
============================================================ -/

section HMultiplier

/-- **Raw complex multiplier** of the Pandrosion iterator at
    a general point α where `x · S_p(α) ≠ 0`:
        `h′(α) = (x − 1) · x · S_p′(α) / (x · S_p(α))²`.
    Derived from the quotient rule applied to
    `h(z) = 1 − (x − 1) / (x · S_p(z))`. -/
noncomputable def hMultC (x : ℂ) (p : ℕ) (α : ℂ) : ℂ :=
  (x - 1) * x * SpDerivC p α / (x * Sp_C p α) ^ 2

/-- **Complex derivative of the Pandrosion iterator at a general
    non-critical point.** At any α with `x ≠ 0` and `S_p(α) ≠ 0`,
    the iterator `h(z) = 1 − (x − 1) / (x · S_p(z))` has complex
    derivative `hMultC x p α`. -/
theorem pandrosion_h_C_hasDerivAt
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (α : ℂ) (hSp : Sp_C p α ≠ 0) :
    HasDerivAt (pandrosion_h_C x p) (hMultC x p α) α := by
  unfold pandrosion_h_C
  have hxSp : x * Sp_C p α ≠ 0 := mul_ne_zero hx hSp
  have h_xSp : HasDerivAt (fun z : ℂ => x * Sp_C p z) (x * SpDerivC p α) α :=
    (Sp_C_hasDerivAt p α).const_mul x
  have h_const : HasDerivAt (fun _ : ℂ => x - 1) (0 : ℂ) α :=
    hasDerivAt_const α (x - 1)
  have h_div := h_const.div h_xSp hxSp
  have h_final := h_div.const_sub (1 : ℂ)
  have h_deriv_eq : -((0 * (x * Sp_C p α) - (x - 1) * (x * SpDerivC p α))
                      / (x * Sp_C p α) ^ 2) = hMultC x p α := by
    unfold hMultC; field_simp; ring
  rw [← h_deriv_eq]
  exact h_final

end HMultiplier

/-! ============================================================
  §26.4  Closed-form multiplier at a fixed point
============================================================ -/

section ClosedFormAtFP

/-- **Closed-form complex multiplier** at a Pandrosion fixed point.
    The paper's `λ_C(p, α) = 1 − p · α^{p−1} / S_p(α)` — direct
    complex analogue of `lambda_fp` from §14.4 and `lambda_closed`
    from §9.2. -/
noncomputable def lambdaClosedC (p : ℕ) (α : ℂ) : ℂ :=
  1 - (p : ℂ) * α ^ (p - 1) / Sp_C p α

/-- **★ Raw multiplier collapses to closed form at a fixed point.**
    At a Pandrosion fixed point α (satisfying `α^p = 1/x` and
    `S_p(α) ≠ 0`), the quotient-rule multiplier
    `(x − 1) · x · S_p′(α) / (x · S_p(α))²` equals the paper's
    closed-form `1 − p · α^{p−1} / S_p(α)`.

    Proof is pure polynomial algebra, combining three identities:
    the Leibniz derivative identity (§26.2), the geometric-sum
    factorization `S_p(α) · (1 − α) = 1 − α^p`, and the fixed-point
    condition `x · α^p = 1`. -/
theorem hMultC_eq_lambdaClosedC_at_fp
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (α : ℂ)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    hMultC x p α = lambdaClosedC p α := by
  have h_xαp : x * α ^ p = 1 := by rw [hfp]; field_simp
  have h_Leibniz := SpDerivC_mul_one_sub p α
  have h_geom := Sp_C_mul_one_sub p α
  have hxSp : x * Sp_C p α ≠ 0 := mul_ne_zero hx hSp
  have hxSp_sq : (x * Sp_C p α) ^ 2 ≠ 0 := pow_ne_zero _ hxSp
  -- Cleared-denominator form of the target identity.
  have h_key : (x - 1) * x * SpDerivC p α * Sp_C p α
             = (Sp_C p α - (p : ℂ) * α ^ (p - 1)) * (x * Sp_C p α) ^ 2 := by
    linear_combination x ^ 2 * Sp_C p α ^ 2 * h_Leibniz
                     - x ^ 2 * SpDerivC p α * Sp_C p α * h_geom
                     + x * SpDerivC p α * Sp_C p α * h_xαp
  unfold hMultC lambdaClosedC
  -- Rewrite `1 − p·α^{p−1}/S_p(α)` as a single fraction.
  have h_rhs_eq : (1 : ℂ) - (p : ℂ) * α ^ (p - 1) / Sp_C p α
                = (Sp_C p α - (p : ℂ) * α ^ (p - 1)) / Sp_C p α := by
    field_simp
  rw [h_rhs_eq]
  -- Cross-multiplication reduces to `h_key`.
  rw [div_eq_div_iff hxSp_sq hSp]
  linear_combination h_key

end ClosedFormAtFP

/-! ============================================================
  §26.5  ★ Combined: `HasDerivAt h = λ_C` at a fixed point
============================================================ -/

section Crown

/-- **★★ Pandrosion multiplier formula at a fixed point.** At a
    Pandrosion fixed point α (where `α^p = 1/x` and `S_p(α) ≠ 0`),
    the Pandrosion iterator `h(z) = 1 − (x−1)/(x·S_p(z))` has
    complex derivative
        `h′(α) = λ_C(p, α) = 1 − p · α^{p−1} / S_p(α)`.

    This is the rigorous complex-calculus realization of the
    algebraic `lambda_fp` from §14.4 — the Julia-dynamics
    multiplier of the Pandrosion iterator at a fixed point,
    given in the paper's closed form. The magnitude `|λ_C|`
    governs local dynamics: `|λ_C| < 1` ⇒ linear contraction,
    `λ_C = 0` ⇒ super-attracting (the Steffensen-accelerated
    regime, treated separately). -/
theorem pandrosion_h_C_hasDerivAt_at_fp
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (α : ℂ)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    HasDerivAt (pandrosion_h_C x p) (lambdaClosedC p α) α := by
  rw [← hMultC_eq_lambdaClosedC_at_fp x hx p α hSp hfp]
  exact pandrosion_h_C_hasDerivAt x hx p α hSp

end Crown

end Pandrosion
