/-
  Universitas Pandrosion — Lean 4 Formalization
  Q_CUBIC IS POSITIVE DEFINITE ON ℝ² — UNCONDITIONAL

  The anchor-step divided difference
      Q(a, s) = s² + a·s + a²
  was shown in AnchorStep.Q_cubic_pos to be positive only under
  the restrictive hypothesis a > 0 ∧ s > 0. That hypothesis is
  unnecessarily strong: Q is positive definite on ℝ², vanishing
  only at the origin.

  Proof certificate (completion of the square):
      Q(a, s) = (s + a/2)² + (3/4)·a²

  As a consequence, whenever the anchor a is nonzero, the
  Pandrosion anchor step
      F_a(s) = a − (a³ − X) / Q(a, s)
  is well-defined for EVERY real s — no positivity hypothesis
  on s required. This cleans up the hypothesis tower of several
  downstream theorems in AnchorStep / MultiStart.

  No scaffold, no sorry — a direct nonlinear arithmetic certificate.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic
import Pandrosion.AnchorStep

namespace Pandrosion

/-! ## §304. Completion of the Square

The identity
    s² + a·s + a² = (s + a/2)² + (3/4)·a²
is proved by `ring`. Both summands on the right are nonnegative,
and their sum is strictly positive unless both vanish — which
forces (a, s) = (0, 0).
-/

/-- **Completion of the square for Q_cubic.** -/
theorem Q_cubic_complete_square (a s : ℝ) :
    Q_cubic a s = (s + a / 2) ^ 2 + 3 * a ^ 2 / 4 := by
  unfold Q_cubic; ring

/-! ## §305. Unconditional Nonnegativity

Q_cubic a s ≥ 0 for every real (a, s).
-/

/-- **Q_cubic is nonnegative on all of ℝ².** -/
theorem Q_cubic_nonneg (a s : ℝ) : 0 ≤ Q_cubic a s := by
  rw [Q_cubic_complete_square]
  have h1 : (0 : ℝ) ≤ (s + a / 2) ^ 2 := sq_nonneg _
  have h2 : (0 : ℝ) ≤ 3 * a ^ 2 / 4 := by positivity
  linarith

/-! ## §306. Zero Locus of Q_cubic

Q_cubic a s = 0 ⟺ a = 0 ∧ s = 0. So the zero set of the
divided difference is the single point (0, 0) ∈ ℝ².
-/

/-- **The zero locus of Q_cubic is exactly the origin.** -/
theorem Q_cubic_eq_zero_iff (a s : ℝ) :
    Q_cubic a s = 0 ↔ a = 0 ∧ s = 0 := by
  rw [Q_cubic_complete_square]
  constructor
  · intro h
    have h1 : (0 : ℝ) ≤ (s + a / 2) ^ 2 := sq_nonneg _
    have h2 : (0 : ℝ) ≤ 3 * a ^ 2 / 4 := by positivity
    have h_sq1 : (s + a / 2) ^ 2 = 0 := by linarith
    have h_a2 : a ^ 2 = 0 := by linarith
    have ha : a = 0 := sq_eq_zero_iff.mp h_a2
    have hsa : s + a / 2 = 0 := sq_eq_zero_iff.mp h_sq1
    rw [ha] at hsa
    refine ⟨ha, ?_⟩
    linarith
  · rintro ⟨ha, hs⟩
    rw [ha, hs]; ring

/-! ## §307. Strict Positivity off the Origin

Off the single degenerate point, Q_cubic is strictly positive.
-/

/-- **Q_cubic is strictly positive away from the origin.** -/
theorem Q_cubic_pos_of_ne_origin (a s : ℝ) (h : ¬(a = 0 ∧ s = 0)) :
    0 < Q_cubic a s := by
  rcases lt_or_eq_of_le (Q_cubic_nonneg a s) with hp | hp
  · exact hp
  · exfalso; exact h ((Q_cubic_eq_zero_iff a s).mp hp.symm)

/-- **Nonzero anchor ⇒ Q_cubic strictly positive for every s.**
    This is the practical form: as long as the anchor a is nonzero,
    the divided-difference denominator cannot vanish, regardless
    of the current iterate s. -/
theorem Q_cubic_pos_of_anchor_ne_zero (a s : ℝ) (ha : a ≠ 0) :
    0 < Q_cubic a s := by
  apply Q_cubic_pos_of_ne_origin
  rintro ⟨ha', _⟩
  exact ha ha'

/-- **Nonzero iterate ⇒ Q_cubic strictly positive for every anchor.** -/
theorem Q_cubic_pos_of_iterate_ne_zero (a s : ℝ) (hs : s ≠ 0) :
    0 < Q_cubic a s := by
  apply Q_cubic_pos_of_ne_origin
  rintro ⟨_, hs'⟩
  exact hs hs'

/-! ## §308. Anchor-Step Well-Definedness

The denominator Q(a, s) of the Pandrosion anchor step never
vanishes when a ≠ 0. This removes the Q-nonvanishing hypothesis
from the downstream fixed-point and multistart theorems.
-/

/-- **Anchor-step denominator never vanishes for nonzero anchor.** -/
theorem anchor_step_denominator_ne_zero (a s : ℝ) (ha : a ≠ 0) :
    Q_cubic a s ≠ 0 :=
  ne_of_gt (Q_cubic_pos_of_anchor_ne_zero a s ha)

/-- **Cleaned anchor fixed-point theorem.**
    For any nonzero anchor a ≠ 0 and any real cube root r of X
    (i.e. r³ = X), the Pandrosion anchor step fixes r.
    The Q-nonvanishing hypothesis of `anchor_fixed_point` is
    automatic here, as a consequence of Q's unconditional positivity. -/
theorem anchor_fixed_point_of_anchor_ne_zero
    (X a r : ℝ) (ha : a ≠ 0) (hX : r ^ 3 = X) :
    pandrosion_anchor_step X a r = r :=
  anchor_fixed_point X a r hX (anchor_step_denominator_ne_zero a r ha)

/-- **Cleaned multistart step at root.**
    For any nonzero anchor a ≠ 0 and any real cube root r of X,
    one epoch of the multistart step maps (a, r) ↦ (r, r). -/
theorem multistart_step_at_root_of_anchor_ne_zero
    (X a r : ℝ) (ha : a ≠ 0) (hX : r ^ 3 = X) :
    multistart_step X a r = (r, r) :=
  multistart_step_at_root X a r hX (anchor_step_denominator_ne_zero a r ha)

/-! ## §309. Lower Bound on Q_cubic

For any anchor a ≠ 0 and any real s, Q_cubic ≥ (3/4)·a². This
explicit quantitative lower bound is useful for Lipschitz and
contraction estimates on the anchor step.
-/

/-- **Explicit lower bound: Q_cubic a s ≥ (3/4)·a².**
    Follows directly from the completion of the square
    Q(a, s) = (s + a/2)² + (3/4)·a² and nonnegativity of squares. -/
theorem Q_cubic_ge_threequarter_anchor_sq (a s : ℝ) :
    3 * a ^ 2 / 4 ≤ Q_cubic a s := by
  rw [Q_cubic_complete_square]
  have : (0 : ℝ) ≤ (s + a / 2) ^ 2 := sq_nonneg _
  linarith

/-- **At the self-anchoring diagonal a = s, Q_cubic matches the
    derivative Q(a,a) = 3a² — consistent with Newton's scaling.** -/
theorem Q_cubic_diagonal (a : ℝ) : Q_cubic a a = 3 * a ^ 2 := by
  unfold Q_cubic; ring

end Pandrosion
