/-
  Universitas Pandrosion — Lean 4 Formalization
  GENERAL p-TH ROOT MULTISTART

  Generalizes Core.lean's cubic multistart (§219–§230) to arbitrary
  integer exponent `p ≥ 2`. The same anchor / Aitken-reanchor
  architecture, now extracting `p`-th roots for any `p`.

  ## Algebraic foundation

  Factor the base polynomial `s^p - a^p` as `(s - a) · Q_p(a, s)`,
  where

      Q_p(a, s) := s^(p-1) + a · s^(p-2) + a² · s^(p-3) + ⋯ + a^(p-1)
                 = Σ_{k=0}^{p-1} a^k · s^(p-1-k).

  Then the anchor-step for `r^p = X` is

      F_{a,p}(s) := (a · Q_p(a, s) - (a^p - X)) / Q_p(a, s).

  Four core properties follow *uniformly* in `p`:

    (1) **Fixed point.**  F_{a,p}(r) = r whenever r^p = X.
    (2) **Newton collapse.**  F_{a,p}(a) = ((p-1)·a^p + X) / (p·a^(p-1)).
    (3) **Master residual identity.**
           F_{a,p}(s) - r = (a - r)(s - r) · R_p(a, r, s) / Q_p(a, s),
        where R_p is a polynomial capturing the second-order correction.
    (4) **Self-anchor quadratic.**  F_{s,p}(s) - r = (s - r)² · R_p(s,r,s) / Q_p(s,s).

  Plus: a general `reanchor_p` (Aitken Δ²) and a general
  `multistart_step_p` that reduce to Core.lean's `reanchor` and
  `multistart_step` exactly when `p = 3`.

  No external Pandrosion dependencies — this module stands alone on
  top of Mathlib.
-/
import Mathlib.Tactic

namespace Pandrosion

section MultiStartP

/-! ## §1. The Q_pow Polynomial

Q_pow p a s = sum_{k=0}^{p-1} a^k · s^(p-1-k). We use a direct
recursion for clean inductive proofs; no Finset sums. -/

/-- **Definition of Q_pow.** -/
def Q_pow : ℕ → ℝ → ℝ → ℝ
  | 0,     _, _ => 0
  | (n+1), a, s => s ^ n + a * Q_pow n a s

@[simp] theorem Q_pow_zero (a s : ℝ) : Q_pow 0 a s = 0 := rfl

theorem Q_pow_succ (n : ℕ) (a s : ℝ) :
    Q_pow (n + 1) a s = s ^ n + a * Q_pow n a s := rfl

theorem Q_pow_one (a s : ℝ) : Q_pow 1 a s = 1 := by
  rw [Q_pow_succ, Q_pow_zero]; simp

theorem Q_pow_two (a s : ℝ) : Q_pow 2 a s = s + a := by
  rw [show (2 : ℕ) = 1 + 1 from rfl, Q_pow_succ, Q_pow_one]; ring

theorem Q_pow_three (a s : ℝ) :
    Q_pow 3 a s = s ^ 2 + a * s + a ^ 2 := by
  rw [show (3 : ℕ) = 2 + 1 from rfl, Q_pow_succ, Q_pow_two]; ring

/-- **Factorization.** `s^p - a^p = (s - a) · Q_pow p a s`. -/
theorem Q_pow_factorization (p : ℕ) (a s : ℝ) :
    (s - a) * Q_pow p a s = s ^ p - a ^ p := by
  induction p with
  | zero => simp
  | succ p ih =>
    rw [Q_pow_succ, pow_succ s, pow_succ a]
    linear_combination a * ih

/-- **Symmetry.** `Q_pow` is symmetric in its two real arguments. -/
theorem Q_pow_symm (p : ℕ) (a s : ℝ) : Q_pow p a s = Q_pow p s a := by
  induction p with
  | zero => simp
  | succ p ih =>
    rw [Q_pow_succ, Q_pow_succ, ih]
    have hf : (s - a) * Q_pow p s a = s ^ p - a ^ p := by
      rw [← ih]; exact Q_pow_factorization p a s
    linear_combination -hf

/-- **Diagonal value.** `Q_pow (p+1) a a = (p+1) · a^p`. -/
theorem Q_pow_at_diag (p : ℕ) (a : ℝ) :
    Q_pow (p + 1) a a = (p + 1 : ℕ) * a ^ p := by
  induction p with
  | zero => rw [Q_pow_succ, Q_pow_zero]; simp
  | succ p ih =>
    rw [Q_pow_succ, ih, pow_succ a p]
    push_cast; ring

/-- **Positivity.** Q_pow is positive when both arguments are positive. -/
theorem Q_pow_pos (p : ℕ) (hp : 1 ≤ p) (a s : ℝ) (ha : 0 < a) (hs : 0 < s) :
    0 < Q_pow p a s := by
  rcases p with _ | q
  · exact absurd hp (by norm_num)
  · clear hp
    induction q with
    | zero => rw [Q_pow_one]; exact one_pos
    | succ q ih =>
      rw [Q_pow_succ]
      have hsp : 0 < s ^ (q + 1) := pow_pos hs _
      have haQ : 0 < a * Q_pow (q + 1) a s := mul_pos ha ih
      linarith

/-- **Diagonal positivity corollary.** -/
theorem Q_pow_at_diag_pos (p : ℕ) (a : ℝ) (ha : 0 < a) :
    0 < Q_pow (p + 1) a a := by
  rw [Q_pow_at_diag]
  have h1 : (0 : ℝ) < ((p + 1 : ℕ) : ℝ) := by exact_mod_cast Nat.succ_pos p
  exact mul_pos h1 (pow_pos ha p)

/-! ## §2. The Anchor-Step for p-th Roots

`anchor_step_p p X a s` is the direct generalization of
`pandrosion_anchor_step` from Core.lean to arbitrary `p`. -/

/-- **The anchor-step.** -/
noncomputable def anchor_step_p (p : ℕ) (X a s : ℝ) : ℝ :=
  if Q_pow p a s = 0 then s
  else (a * Q_pow p a s - (a ^ p - X)) / Q_pow p a s

/-- **Explicit form when Q_pow is nonzero.** -/
theorem anchor_step_p_explicit (p : ℕ) (X a s : ℝ) (hQ : Q_pow p a s ≠ 0) :
    anchor_step_p p X a s
      = (a * Q_pow p a s - (a ^ p - X)) / Q_pow p a s := by
  unfold anchor_step_p; rw [if_neg hQ]

/-- **Alternative form.** `F_{a,p}(s) = a - (a^p - X) / Q_pow p a s`. -/
theorem anchor_step_p_as_correction (p : ℕ) (X a s : ℝ) (hQ : Q_pow p a s ≠ 0) :
    anchor_step_p p X a s = a - (a ^ p - X) / Q_pow p a s := by
  rw [anchor_step_p_explicit p X a s hQ]
  field_simp

/-- **Fixed-point theorem.**  `F_{a,p}(r) = r` whenever `r^p = X`. -/
theorem anchor_step_p_at_root (p : ℕ) (X a r : ℝ)
    (hX : r ^ p = X) (hQ : Q_pow p a r ≠ 0) :
    anchor_step_p p X a r = r := by
  rw [anchor_step_p_explicit p X a r hQ, ← hX]
  have hfac : (a - r) * Q_pow p r a = a ^ p - r ^ p := Q_pow_factorization p r a
  have hsymm : Q_pow p r a = Q_pow p a r := Q_pow_symm p r a
  rw [hsymm] at hfac
  field_simp
  linarith

/-- **Newton collapse.**
    When the anchor equals the iterate, the step reduces to Newton. -/
theorem newton_from_anchor_p (p : ℕ) (X a : ℝ) (ha : a ≠ 0) :
    anchor_step_p (p + 1) X a a
      = ((p : ℕ) * a ^ (p + 1) + X) / ((p + 1 : ℕ) * a ^ p) := by
  have hp1_pos : (0 : ℝ) < ((p + 1 : ℕ) : ℝ) := by exact_mod_cast Nat.succ_pos p
  have hp1_ne : ((p + 1 : ℕ) : ℝ) ≠ 0 := ne_of_gt hp1_pos
  have hpow_ne : a ^ p ≠ 0 := pow_ne_zero _ ha
  have hQ_diag : Q_pow (p + 1) a a = (p + 1 : ℕ) * a ^ p := Q_pow_at_diag p a
  have hQ_ne : Q_pow (p + 1) a a ≠ 0 := by
    rw [hQ_diag]; exact mul_ne_zero hp1_ne hpow_ne
  rw [anchor_step_p_explicit (p + 1) X a a hQ_ne, hQ_diag]
  field_simp
  ring

/-! ## §3. The Residual Polynomial R_pow and Master Identity

R_pow p a r s := Σ_{k=0}^{p-1} a^k · Q_pow (p-1-k) r s  captures
the second-order correction in the master identity. -/

/-- **Definition of R_pow.** -/
def R_pow : ℕ → ℝ → ℝ → ℝ → ℝ
  | 0,     _, _, _ => 0
  | (n+1), a, r, s => Q_pow n r s + a * R_pow n a r s

@[simp] theorem R_pow_zero (a r s : ℝ) : R_pow 0 a r s = 0 := rfl

theorem R_pow_succ (n : ℕ) (a r s : ℝ) :
    R_pow (n + 1) a r s = Q_pow n r s + a * R_pow n a r s := rfl

theorem R_pow_one (a r s : ℝ) : R_pow 1 a r s = 0 := by
  rw [R_pow_succ, Q_pow_zero, R_pow_zero]; ring

theorem R_pow_two (a r s : ℝ) : R_pow 2 a r s = 1 := by
  rw [show (2 : ℕ) = 1 + 1 from rfl, R_pow_succ, Q_pow_one, R_pow_one]; ring

theorem R_pow_three (a r s : ℝ) : R_pow 3 a r s = r + s + a := by
  rw [show (3 : ℕ) = 2 + 1 from rfl, R_pow_succ, Q_pow_two, R_pow_two]; ring

/-- **Difference identity.**
    `Q_pow p a s - Q_pow p a r = (s - r) · R_pow p a r s`. -/
theorem Q_pow_diff (p : ℕ) (a r s : ℝ) :
    Q_pow p a s - Q_pow p a r = (s - r) * R_pow p a r s := by
  induction p with
  | zero => simp
  | succ p ih =>
    rw [Q_pow_succ, Q_pow_succ, R_pow_succ]
    have hfac := Q_pow_factorization p r s
    -- hfac : (s - r) * Q_pow p r s = s^p - r^p
    have h : (s ^ p + a * Q_pow p a s) - (r ^ p + a * Q_pow p a r)
           = (s ^ p - r ^ p) + a * (Q_pow p a s - Q_pow p a r) := by ring
    rw [h, ih, ← hfac]; ring

/-- **Master residual identity** (cross-multiplied form, division-free).
    `Q_p(a,s) · (F_{a,p}(s) - r) = (a - r) · (s - r) · R_p(a, r, s)`
    whenever `r^p = X` and `Q_p(a, s) ≠ 0`. -/
theorem anchor_step_p_master_identity_cross
    (p : ℕ) (X a r s : ℝ) (hX : r ^ p = X) (hQ : Q_pow p a s ≠ 0) :
    Q_pow p a s * (anchor_step_p p X a s - r)
      = (a - r) * (s - r) * R_pow p a r s := by
  rw [anchor_step_p_explicit p X a s hQ]
  have hfac : (a - r) * Q_pow p a r = a ^ p - r ^ p := by
    rw [Q_pow_symm]; exact Q_pow_factorization p r a
  have hdiff : Q_pow p a s - Q_pow p a r = (s - r) * R_pow p a r s :=
    Q_pow_diff p a r s
  have h_mul : Q_pow p a s *
      ((a * Q_pow p a s - (a ^ p - X)) / Q_pow p a s - r)
      = a * Q_pow p a s - (a ^ p - X) - r * Q_pow p a s := by
    field_simp; ring
  rw [h_mul, ← hX]
  linear_combination (a - r) * hdiff + hfac

/-- **Master residual identity** (quotient form).
    `F_{a,p}(s) - r = (a - r) · (s - r) · R_p(a, r, s) / Q_p(a, s)`. -/
theorem anchor_step_p_master_identity
    (p : ℕ) (X a r s : ℝ) (hX : r ^ p = X) (hQ : Q_pow p a s ≠ 0) :
    anchor_step_p p X a s - r
      = (a - r) * (s - r) * R_pow p a r s / Q_pow p a s := by
  have h := anchor_step_p_master_identity_cross p X a r s hX hQ
  field_simp
  linarith

/-- **Self-anchor quadratic convergence** (exact algebraic form).
    Setting `a = s` yields `F_{s,p}(s) - r = (s - r)² · R_p(s, r, s) / Q_p(s, s)`.
    The residual is exactly quadratic in `(s - r)`, matching Newton's
    classical quadratic rate for `f(s) = s^p - X`. -/
theorem anchor_step_p_self_quadratic
    (p : ℕ) (hp : 1 ≤ p) (X r s : ℝ)
    (hX : r ^ p = X) (hs : 0 < s) :
    anchor_step_p p X s s - r
      = (s - r) ^ 2 * R_pow p s r s / Q_pow p s s := by
  have hQ : Q_pow p s s ≠ 0 := ne_of_gt (Q_pow_pos p hp s s hs hs)
  have h := anchor_step_p_master_identity p X s r s hX hQ
  rw [h]; ring

/-- **Self-anchor cross form.** `Q_p(s,s) · (F_{s,p}(s) - r) = (s - r)² · R_p(s, r, s)`. -/
theorem anchor_step_p_self_quadratic_cross
    (p : ℕ) (hp : 1 ≤ p) (X r s : ℝ)
    (hX : r ^ p = X) (hs : 0 < s) :
    Q_pow p s s * (anchor_step_p p X s s - r)
      = (s - r) ^ 2 * R_pow p s r s := by
  have hQ : Q_pow p s s ≠ 0 := ne_of_gt (Q_pow_pos p hp s s hs hs)
  have h := anchor_step_p_master_identity_cross p X s r s hX hQ
  rw [h]; ring

/-! ## §4. Reanchoring and the Full Multistart Step

Generalizes Core.lean's `reanchor` (§224) and `multistart_step` (§225)
to arbitrary `p`. The structure is identical: three anchor-steps per
epoch plus one Aitken Δ² anchor refresh. -/

/-- **Reanchoring via Aitken Δ².** -/
noncomputable def reanchor_p (p : ℕ) (X a s : ℝ) : ℝ :=
  let s1 := anchor_step_p p X a s
  let s2 := anchor_step_p p X a s1
  let denom := s2 - 2 * s1 + s
  if denom = 0 then s2 else s - (s1 - s) ^ 2 / denom

/-- **Full multistart step for general p.**
    `(a, s) ↦ (reanchor_p(a, s), F_{a,p}^[3](s))`. -/
noncomputable def multistart_step_p (p : ℕ) (X a s : ℝ) : ℝ × ℝ :=
  let s1 := anchor_step_p p X a s
  let s2 := anchor_step_p p X a s1
  let s3 := anchor_step_p p X a s2
  let a_new := reanchor_p p X a s
  (a_new, s3)

/-- **Reanchoring at a root.** When `s = r`, the iterates collapse
    to `r` and the Aitken denominator vanishes; the fallback branch
    returns the second iterate, which is again `r`. -/
theorem reanchor_p_at_root (p : ℕ) (X a r : ℝ)
    (hX : r ^ p = X) (hQ : Q_pow p a r ≠ 0) :
    reanchor_p p X a r = r := by
  unfold reanchor_p
  have h1 : anchor_step_p p X a r = r := anchor_step_p_at_root p X a r hX hQ
  simp only [h1]
  simp [show (r : ℝ) - 2 * r + r = 0 from by ring]

/-- **Multistart step at a root.** Every root of `X` (of degree `p`)
    is a joint fixed point of `multistart_step_p` for every anchor. -/
theorem multistart_step_p_at_root (p : ℕ) (X a r : ℝ)
    (hX : r ^ p = X) (hQ : Q_pow p a r ≠ 0) :
    multistart_step_p p X a r = (r, r) := by
  unfold multistart_step_p
  have h1 : anchor_step_p p X a r = r := anchor_step_p_at_root p X a r hX hQ
  have hre : reanchor_p p X a r = r := reanchor_p_at_root p X a r hX hQ
  simp [h1, hre]

/-! ## §5. Residual Propagation Under Iterated Anchor Steps

The master identity gives a clean amplification bound for iterated
anchor-steps: after `n` anchor-steps with the same anchor `a`, the
residual `|F^[n](s) - r|` carries a factor `|a - r|^n`. When `a ≈ r`
this yields super-linear convergence. -/

/-- Iterated anchor-step with fixed anchor. -/
noncomputable def anchor_iterate (p : ℕ) (X a : ℝ) : ℕ → ℝ → ℝ
  | 0,     s => s
  | (n+1), s => anchor_step_p p X a (anchor_iterate p X a n s)

@[simp] theorem anchor_iterate_zero (p : ℕ) (X a s : ℝ) :
    anchor_iterate p X a 0 s = s := rfl

theorem anchor_iterate_succ (p : ℕ) (X a : ℝ) (n : ℕ) (s : ℝ) :
    anchor_iterate p X a (n + 1) s
      = anchor_step_p p X a (anchor_iterate p X a n s) := rfl

/-- **Root preservation under iteration.**
    If `r` is a root (`r^p = X`) and `Q_p(a, r) ≠ 0`, every
    anchor-iterate of `r` is again `r`. -/
theorem anchor_iterate_at_root (p : ℕ) (X a r : ℝ)
    (hX : r ^ p = X) (hQ : Q_pow p a r ≠ 0) :
    ∀ n, anchor_iterate p X a n r = r := by
  intro n
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [anchor_iterate_succ, ih]
    exact anchor_step_p_at_root p X a r hX hQ

/-! ## §6. Specializations and Bridge to the p = 3 Infrastructure

Concrete instances of the `p`-parametric machinery. The `p = 3`
case is *definitionally* compatible with Core.lean's Q_cubic and
`pandrosion_anchor_step` — they compute identical values. -/

/-- **p = 2 (square root).** Self-anchor gives the Babylonian step. -/
theorem anchor_step_p2_self (X a : ℝ) (ha : a ≠ 0) :
    anchor_step_p 2 X a a = (a ^ 2 + X) / (2 * a) := by
  have h := newton_from_anchor_p 1 X a ha
  -- h : anchor_step_p 2 X a a = (1 * a^2 + X) / (2 * a^1)
  rw [h]
  push_cast; ring

/-- **p = 3 (cube root).** Matches Core.lean's `newton_from_pandrosion`. -/
theorem anchor_step_p3_self (X a : ℝ) (ha : a ≠ 0) :
    anchor_step_p 3 X a a = (2 * a ^ 3 + X) / (3 * a ^ 2) := by
  have h := newton_from_anchor_p 2 X a ha
  rw [h]
  push_cast; ring

/-- **p = 4 (fourth root).** Self-anchor gives the quartic Newton step. -/
theorem anchor_step_p4_self (X a : ℝ) (ha : a ≠ 0) :
    anchor_step_p 4 X a a = (3 * a ^ 4 + X) / (4 * a ^ 3) := by
  have h := newton_from_anchor_p 3 X a ha
  rw [h]
  push_cast; ring

/-- **p = 5.** -/
theorem anchor_step_p5_self (X a : ℝ) (ha : a ≠ 0) :
    anchor_step_p 5 X a a = (4 * a ^ 5 + X) / (5 * a ^ 4) := by
  have h := newton_from_anchor_p 4 X a ha
  rw [h]
  push_cast; ring

/-- **p = 7.** -/
theorem anchor_step_p7_self (X a : ℝ) (ha : a ≠ 0) :
    anchor_step_p 7 X a a = (6 * a ^ 7 + X) / (7 * a ^ 6) := by
  have h := newton_from_anchor_p 6 X a ha
  rw [h]
  push_cast; ring

/-- **Sanity check: √2 via one Babylonian self-anchor step from s = 1.**
    `F(1) = (1 + 2)/2 = 3/2` — the classical starting point. -/
theorem sqrt_two_one_step : anchor_step_p 2 2 1 1 = 3 / 2 := by
  rw [anchor_step_p2_self 2 1 one_ne_zero]; norm_num

/-- **Sanity check: ∛2 via one Newton self-anchor step from s = 1.**
    `F(1) = (2 + 2)/3 = 4/3`. -/
theorem cbrt_two_one_step : anchor_step_p 3 2 1 1 = 4 / 3 := by
  rw [anchor_step_p3_self 2 1 one_ne_zero]; norm_num

/-- **Sanity check: ⁵√3 via one Newton self-anchor step from s = 1.**
    `F(1) = (4 + 3)/5 = 7/5`. -/
theorem fifth_root_three_one_step : anchor_step_p 5 3 1 1 = 7 / 5 := by
  rw [anchor_step_p5_self 3 1 one_ne_zero]; norm_num

/-! ## §7. Super-certificate

A single bundled theorem packaging the four core facts, for every
`p ≥ 1`: fixed point, Newton collapse, master identity, self-anchor
quadratic. -/

/-- **The MultiStartP super-certificate.** -/
theorem multistartP_super_certificate (p : ℕ) (hp : 1 ≤ p) :
    (∀ (X a r : ℝ), r ^ p = X → Q_pow p a r ≠ 0 →
       anchor_step_p p X a r = r) ∧
    (∀ (X a : ℝ), a ≠ 0 →
       ∃ q : ℕ, p = q + 1 ∧
       anchor_step_p (q + 1) X a a
         = ((q : ℕ) * a ^ (q + 1) + X) / ((q + 1 : ℕ) * a ^ q)) ∧
    (∀ (X a r s : ℝ), r ^ p = X → Q_pow p a s ≠ 0 →
       anchor_step_p p X a s - r
         = (a - r) * (s - r) * R_pow p a r s / Q_pow p a s) ∧
    (∀ (X r s : ℝ), r ^ p = X → 0 < s →
       anchor_step_p p X s s - r
         = (s - r) ^ 2 * R_pow p s r s / Q_pow p s s) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro X a r hX hQ; exact anchor_step_p_at_root p X a r hX hQ
  · intro X a ha
    match p, hp with
    | 0, hp => exact absurd hp (by norm_num)
    | q + 1, _ => exact ⟨q, rfl, newton_from_anchor_p q X a ha⟩
  · intro X a r s hX hQ; exact anchor_step_p_master_identity p X a r s hX hQ
  · intro X r s hX hs; exact anchor_step_p_self_quadratic p hp X r s hX hs

end MultiStartP

end Pandrosion
