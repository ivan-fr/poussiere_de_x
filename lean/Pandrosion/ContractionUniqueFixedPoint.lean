/-
  Universitas Pandrosion — Lean 4 Formalization
  CONTRACTION ⇒ UNIQUE REAL FIXED POINT

  Real theorem connecting NoCycles and AnchorStep.

  If a map F : ℝ → ℝ satisfies the pointwise contraction
      |F s - r| ≤ c · |s - r|    (0 ≤ c < 1)
  toward an anchor point r, then EVERY real fixed point of F
  coincides with r. In particular, F has at most one real
  fixed point — the contraction target.

  Applied to the Pandrosion anchor step F_a (AnchorStep.lean):
  whenever F_a contracts toward a cube root r of X with rate c < 1,
  r is the unique real fixed point of F_a. Together with
  `anchor_fixed_point` (which furnishes existence), this characterises
  the anchor step's fixed-point set as the singleton {r}.

  No scaffold, no sorry — proved directly from the contraction
  inequality and the strict monotonicity c·t < t for t > 0, c < 1.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic
import Pandrosion.AnchorStep

namespace Pandrosion

/-! ## §301. One-Step Fixed Point Uniqueness

A single application of the contraction inequality at a
fixed point s (F s = s) already forces |s - r| = 0: the
chain `|s - r| = |F s - r| ≤ c·|s - r|` combined with
`c < 1` collapses the distance.
-/

/-- **Fixed-point uniqueness under pointwise contraction.**
    If `|F s - r| ≤ c · |s - r|` for every s and some c ∈ [0, 1),
    then any fixed point of F equals r.

    Proof: at a fixed point, the inequality becomes
    `|s - r| ≤ c · |s - r|`. If `|s - r| > 0` this forces `c ≥ 1`,
    contradicting `c < 1`. So `|s - r| = 0`, i.e. `s = r`. -/
theorem fixed_point_unique_under_contraction
    (F : ℝ → ℝ) (r : ℝ) (c : ℝ)
    (_hc_nn : 0 ≤ c) (hc_lt : c < 1)
    (h_contract : ∀ s, |F s - r| ≤ c * |s - r|)
    {s : ℝ} (h_fix : F s = s) : s = r := by
  by_contra h_ne
  have h_pos : (0 : ℝ) < |s - r| := abs_pos.mpr (sub_ne_zero.mpr h_ne)
  have h1 : |F s - r| ≤ c * |s - r| := h_contract s
  rw [h_fix] at h1
  have h_strict : c * |s - r| < |s - r| := by
    calc c * |s - r|
      _ < 1 * |s - r| := mul_lt_mul_of_pos_right hc_lt h_pos
      _ = |s - r| := one_mul _
  linarith

/-! ## §302. Anchor-Step Fixed-Point Characterisation

Specialisation of the general result to the Pandrosion anchor
step F_a (AnchorStep.pandrosion_anchor_step). Combined with
`AnchorStep.anchor_fixed_point` (which proves r IS a fixed point
whenever r³ = X and Q(a, r) ≠ 0), we obtain a full `↔`
characterisation of the fixed-point set of F_a.
-/

/-- **Anchor-step: at most one real fixed point under contraction.**
    If the Pandrosion anchor step F_a contracts toward r with rate
    c ∈ [0, 1), then every real s with F_a(s) = s satisfies s = r. -/
theorem anchor_step_fixed_point_unique
    (X a r : ℝ) (c : ℝ)
    (hc_nn : 0 ≤ c) (hc_lt : c < 1)
    (h_contract : ∀ s, |pandrosion_anchor_step X a s - r| ≤ c * |s - r|)
    {s : ℝ} (h_fix : pandrosion_anchor_step X a s = s) : s = r :=
  fixed_point_unique_under_contraction
    (pandrosion_anchor_step X a) r c hc_nn hc_lt h_contract h_fix

/-- **Anchor-step fixed-point set = {r}, under contraction.**
    If r³ = X, Q(a, r) ≠ 0, and F_a contracts toward r with rate c < 1,
    then for every s ∈ ℝ: `F_a(s) = s ↔ s = r`.

    This closes the characterisation: existence (from `anchor_fixed_point`)
    and uniqueness (from contraction) combine to pin the fixed-point set. -/
theorem anchor_step_fixed_point_iff
    (X a r : ℝ) (c : ℝ)
    (hX : r ^ 3 = X) (hQ : Q_cubic a r ≠ 0)
    (hc_nn : 0 ≤ c) (hc_lt : c < 1)
    (h_contract : ∀ s, |pandrosion_anchor_step X a s - r| ≤ c * |s - r|)
    (s : ℝ) :
    pandrosion_anchor_step X a s = s ↔ s = r := by
  constructor
  · intro h_fix
    exact anchor_step_fixed_point_unique X a r c hc_nn hc_lt h_contract h_fix
  · intro h_eq
    rw [h_eq]
    exact anchor_fixed_point X a r hX hQ

/-! ## §303. Multi-Start Coherence

In the Pandrosion multi-start architecture, d different orbits
run from d different anchors a_0, …, a_{d-1}. Each orbit, if it
contracts toward SOME cube root r_i of X, has r_i as its unique
real fixed point.

The next theorem states the cross-orbit coherence: two orbits
cannot silently converge to the same iterate unless they share
the same target root. This is the algebraic skeleton behind the
Voronoï-coverage guarantee in MultiStart.§215.
-/

/-- **Distinct contraction targets force distinct fixed points.**
    If F_{a₁} contracts toward r₁ with rate c₁ < 1, F_{a₂} contracts
    toward r₂ with rate c₂ < 1, and some real s is simultaneously a
    fixed point of BOTH anchor steps, then r₁ = r₂.

    Proof: uniqueness under contraction gives s = r₁ and s = r₂. -/
theorem shared_fixed_point_forces_equal_targets
    (X a₁ a₂ r₁ r₂ : ℝ) (c₁ c₂ : ℝ)
    (hc₁_nn : 0 ≤ c₁) (hc₁_lt : c₁ < 1)
    (hc₂_nn : 0 ≤ c₂) (hc₂_lt : c₂ < 1)
    (h_contract₁ : ∀ s, |pandrosion_anchor_step X a₁ s - r₁| ≤ c₁ * |s - r₁|)
    (h_contract₂ : ∀ s, |pandrosion_anchor_step X a₂ s - r₂| ≤ c₂ * |s - r₂|)
    (s : ℝ)
    (h_fix₁ : pandrosion_anchor_step X a₁ s = s)
    (h_fix₂ : pandrosion_anchor_step X a₂ s = s) :
    r₁ = r₂ := by
  have h1 : s = r₁ :=
    anchor_step_fixed_point_unique X a₁ r₁ c₁ hc₁_nn hc₁_lt h_contract₁ h_fix₁
  have h2 : s = r₂ :=
    anchor_step_fixed_point_unique X a₂ r₂ c₂ hc₂_nn hc₂_lt h_contract₂ h_fix₂
  rw [← h1, ← h2]

end Pandrosion
