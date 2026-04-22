/-
  Universitas Pandrosion — Legacy / Anchor step.

  Algebra of the derivative-free *anchored* Pandrosion step for the
  cube-root problem `z³ = X`:

      F_a(s) = a − (a³ − X) / Q(a, s),     Q(a, s) = s² + a·s + a².

  Content:

    §L27  Cubic divided difference `Q_cubic`.
    §L28  Anchor step `pandrosion_anchor_step`, fixed-point property.
    §L29  Newton's step as the `a = s` degenerate case.
    §L30  Ratio form, residual numerator.
    §L31  Aitken reanchoring, idempotence at roots.
    §L32  Aitken exactness on a pure geometric sequence (`lam, lam²`).

  Ported from `AnchorStep.lean` on the `main_previous` branch.
  The speculative closing section ("§227 Structural Non-Holomorphicity"
  and "McMullen Impossibility Circumvention") was dropped as it
  offered no mathematical content beyond sloganeering.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §L27. Cubic divided difference `Q_cubic` -/

/-- `s³ − a³ = (s − a)(s² + a·s + a²)`. -/
theorem cubic_factorization (a s : ℝ) :
    s ^ 3 - a ^ 3 = (s - a) * (s ^ 2 + a * s + a ^ 2) := by
  ring

/-- Divided difference `Q(a, s) = s² + a·s + a²` — a polynomial
    (so no division needed), making `F_a(s)` derivative-free. -/
def Q_cubic (a s : ℝ) : ℝ := s ^ 2 + a * s + a ^ 2

/-- `Q(a, a) = 3a²` — coalescence limit is the Newton derivative. -/
theorem Q_selfevaluation (a : ℝ) :
    Q_cubic a a = 3 * a ^ 2 := by
  unfold Q_cubic; ring

/-- `Q(a, s) > 0` whenever `a > 0` and `s > 0`. -/
theorem Q_cubic_pos (a s : ℝ) (ha : a > 0) (_hs : s > 0) :
    Q_cubic a s > 0 := by
  unfold Q_cubic
  nlinarith [sq_nonneg (s + a / 2)]

/-! ## §L28. Anchor step and fixed-point property -/

/-- Anchored Pandrosion step: `F_a(s) = a − (a³ − X)/Q(a, s)`. -/
noncomputable def pandrosion_anchor_step (X a s : ℝ) : ℝ :=
  a - (a ^ 3 - X) / Q_cubic a s

/-- Cross-multiplied fixed-point identity:
    `a · Q(a,r) − (a³ − X) = r · Q(a,r)` when `r³ = X`. -/
theorem anchor_fixed_point_cross (X a r : ℝ) (hX : r ^ 3 = X) :
    a * Q_cubic a r - (a ^ 3 - X) = r * Q_cubic a r := by
  unfold Q_cubic; subst hX; ring

/-- **Every cube root of `X` is a fixed point of `F_a`** (for any anchor
    `a` with `Q(a, r) ≠ 0`). -/
theorem anchor_fixed_point (X a r : ℝ) (hX : r ^ 3 = X)
    (hQ : Q_cubic a r ≠ 0) :
    pandrosion_anchor_step X a r = r := by
  unfold pandrosion_anchor_step
  have h := anchor_fixed_point_cross X a r hX
  field_simp
  linarith

/-! ## §L29. Newton's step as the `a = s` degenerate case -/

/-- Cross-multiplied Newton identity:
    `a · 3a² − (a³ − X) = 2a³ + X`. -/
theorem newton_is_pandrosion_cross (X a : ℝ) :
    a * (3 * a ^ 2) - (a ^ 3 - X) = 2 * a ^ 3 + X := by
  ring

/-- **Self-anchoring collapses to Newton.** `F_a(a) = (2a³ + X)/(3a²)`. -/
theorem newton_from_pandrosion (X a : ℝ) (ha : a ≠ 0) :
    pandrosion_anchor_step X a a = (2 * a ^ 3 + X) / (3 * a ^ 2) := by
  unfold pandrosion_anchor_step Q_cubic
  have h3a2 : a ^ 2 + a * a + a ^ 2 ≠ 0 := by
    have : a ^ 2 + a * a + a ^ 2 = 3 * a ^ 2 := by ring
    rw [this]; positivity
  field_simp
  ring

/-! ## §L30. Ratio form of the anchor step -/

/-- Numerator identity: `a · Q(a,s) − (a³ − X) = a·s² + a²·s + X`. -/
theorem anchor_step_numerator (X a s : ℝ) :
    a * Q_cubic a s - (a ^ 3 - X) = a * s ^ 2 + a ^ 2 * s + X := by
  unfold Q_cubic; ring

/-- Explicit form: `F_a(s) = (a·s² + a²·s + X) / Q(a, s)`. -/
theorem anchor_step_explicit (X a s : ℝ) (hQ : Q_cubic a s ≠ 0) :
    pandrosion_anchor_step X a s =
    (a * s ^ 2 + a ^ 2 * s + X) / Q_cubic a s := by
  unfold pandrosion_anchor_step
  have h := anchor_step_numerator X a s
  field_simp
  linarith

/-- Pure re-statement of `newton_from_pandrosion` for pedagogical clarity. -/
theorem standard_vs_anchor (X a : ℝ) (ha : a ≠ 0) :
    pandrosion_anchor_step X a a = (2 * a ^ 3 + X) / (3 * a ^ 2) :=
  newton_from_pandrosion X a ha

/-! ## §L31. Aitken Δ² reanchoring and idempotence at roots -/

/-- Reanchoring: given `s`, compute `F_a(s)` and `F_a²(s)`, then
    apply Aitken Δ²; fall back to `s₂` when the denominator vanishes
    (covering the fixed-point case). -/
noncomputable def reanchor (X a s : ℝ) : ℝ :=
  let s1 := pandrosion_anchor_step X a s
  let s2 := pandrosion_anchor_step X a s1
  let denom := s2 - 2 * s1 + s
  if denom = 0 then s2 else s - (s1 - s) ^ 2 / denom

/-- One multi-start step: three anchor iterations plus a reanchor. -/
noncomputable def multistart_step (X a s : ℝ) : ℝ × ℝ :=
  let s1 := pandrosion_anchor_step X a s
  let s2 := pandrosion_anchor_step X a s1
  let s3 := pandrosion_anchor_step X a s2
  let a_new := reanchor X a s
  (a_new, s3)

/-- `F_a` is idempotent at every cube root. -/
theorem anchor_step_at_root (X a r : ℝ) (hX : r ^ 3 = X)
    (hQ : Q_cubic a r ≠ 0) :
    pandrosion_anchor_step X a r = r :=
  anchor_fixed_point X a r hX hQ

/-- Three consecutive anchor steps starting from a root return the root. -/
theorem three_steps_at_root (X a r : ℝ) (hX : r ^ 3 = X)
    (hQ : Q_cubic a r ≠ 0) :
    let s1 := pandrosion_anchor_step X a r
    let s2 := pandrosion_anchor_step X a s1
    s1 = r ∧ s2 = r := by
  refine ⟨anchor_fixed_point X a r hX hQ, ?_⟩
  have h1 := anchor_fixed_point X a r hX hQ
  rw [h1]
  exact anchor_fixed_point X a r hX hQ

/-- The Aitken denominator vanishes at a fixed point: `r − 2r + r = 0`. -/
theorem aitken_denom_at_root (X a r : ℝ) (hX : r ^ 3 = X)
    (hQ : Q_cubic a r ≠ 0) :
    let s1 := pandrosion_anchor_step X a r
    let s2 := pandrosion_anchor_step X a s1
    s2 - 2 * s1 + r = 0 := by
  have ⟨h1, h2⟩ := three_steps_at_root X a r hX hQ
  simp only
  rw [h2, h1]; ring

/-- Reanchoring at a root returns the root (fallback branch). -/
theorem reanchor_at_root (X a r : ℝ) (hX : r ^ 3 = X)
    (hQ : Q_cubic a r ≠ 0) :
    reanchor X a r = r := by
  unfold reanchor
  have h1 : pandrosion_anchor_step X a r = r := anchor_fixed_point X a r hX hQ
  simp only [h1]
  simp [show r - 2 * r + r = 0 from by ring]

/-- **Full multi-start step is a fixed point at every cube root.**
    `multistart_step(X, a, r) = (r, r)` for all anchors `a` with
    `Q(a, r) ≠ 0`. -/
theorem multistart_step_at_root (X a r : ℝ) (hX : r ^ 3 = X)
    (hQ : Q_cubic a r ≠ 0) :
    multistart_step X a r = (r, r) := by
  unfold multistart_step
  have h1 : pandrosion_anchor_step X a r = r := anchor_fixed_point X a r hX hQ
  have h2 : reanchor X a r = r := reanchor_at_root X a r hX hQ
  simp only [h1, h2]

/-! ## §L32. Aitken exactness on a pure geometric sequence -/

/-- **Aitken Δ² is exact on `s, r + lam·e, r + lam²·e`.**
    Variant of `Legacy.Basic.aitken_perfect_extrapolation`, phrased
    with the explicit `lam²·e` form (rather than a free `s₂`). -/
theorem aitken_exact_geometric (r e lam : ℝ) (hlam : lam ≠ 1) (he : e ≠ 0) :
    let s := r + e
    let s1 := r + lam * e
    let s2 := r + lam ^ 2 * e
    let denom := s2 - 2 * s1 + s
    denom ≠ 0 ∧ s - (s1 - s) ^ 2 / denom = r := by
  simp only
  refine ⟨?_, ?_⟩
  · intro h
    have hkey : (lam - 1) ^ 2 * e = 0 := by linarith
    rcases mul_eq_zero.mp hkey with h1 | h2
    · exact hlam (by nlinarith [sq_eq_zero_iff.mp h1])
    · exact he h2
  · have hdenom :
        r + lam ^ 2 * e - 2 * (r + lam * e) + (r + e) = (lam - 1) ^ 2 * e := by
      ring
    rw [hdenom]
    have hdenom_ne : (lam - 1) ^ 2 * e ≠ 0 :=
      mul_ne_zero (pow_ne_zero 2 (sub_ne_zero.mpr hlam)) he
    have hnum : (r + lam * e - (r + e)) ^ 2 = (lam - 1) ^ 2 * e ^ 2 := by ring
    rw [hnum]
    field_simp
    ring

/-- Aitken Δ² denominator at the Pandrosion rate `lam = −1/5`:
    `(lam − 1)² = (−6/5)² = 36/25`. Shows the accelerator has a
    well-controlled denominator bounded away from zero. -/
theorem aitken_pandrosion_denominator :
    (-(1 : ℝ) / 5 - 1) ^ 2 = 36 / 25 := by norm_num

end Pandrosion
