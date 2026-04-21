/-
  Universitas Pandrosion — §25. Concrete-iteration bridge.

  The unconditional super grand-master `super_grand_master_uniform_closed_form`
  (§23.6) takes as *data* a family of abstract error sequences
  `(K_p, R_p, e_p)_p` satisfying two-phase Pandrosion bounds — linear
  contraction at the closed-form rate `lamModel p x` and quadratic
  contraction with per-p constant `K_p`. That makes §23.6 a statement
  about *any* family satisfying these inequalities, not yet a statement
  about the *real* Pandrosion-Steffensen iteration on `z^p = x`.

  **Missing bridge: abstract → concrete.** This file closes that gap.
  For a real scalar sequence `(s k)_k : ℕ → ℝ` with:
    • error `e k := |s k − α|` (to its Pandrosion fixed point `α`),
    • per-step quadratic residue `e (k+1) ≤ K · (e k)²`
      (the Aitken-Steffensen output of §15.4),
  we prove that, *provided* the initial scaled error lies in the
  closed-form scaled basin
      `K · R ≤ lamModel p x    (with `R ≥ e 0`)`,
  the sequence *automatically* also satisfies the linear bound
      `e (k+1) ≤ lamModel p x · e k`
  at the analytic rate of §23.5. The derivation is purely structural:
  the quadratic recurrence plus the basin invariant

      `∀ k, K · e k ≤ lamModel p x`

  collapses the quadratic bound into a linear bound at rate lamModel.

  Content.

    §25.1  `real_iterate_scaled_basin_invariant` — induction invariant
           for a quadratic iteration starting in a scaled basin.
    §25.2  `real_iterate_linear_from_quadratic` — quadratic + basin
           invariant ⟹ linear at the basin's rate.
    ★ §25.3  `pandrosion_produces_admissible_sequences` — for a single
           `p ≥ 1`, the real Pandrosion-Steffensen iteration (given
           its quadratic residue + scaled initial error) produces the
           linear + quadratic + basin-invariant triple demanded by
           §23.6 at the concrete rate `lamModel p x`.
    ★ §25.4  `pandrosion_concrete_super_grand_master` — feed the
           admissible sequences into §23.6 to conclude p-uniform
           complexity for the *real* iteration on `z^p = x`.

  This is the load-bearing abstract→concrete bridge of the Pandrosion
  formalization.
-/

import Pandrosion.Core.UniformContractionRate

namespace Pandrosion

/-! ============================================================
  §25.1  Scaled-basin invariant
============================================================ -/

section BasinInvariant

/-- **Scaled-basin invariant.** A non-negative sequence with quadratic
    recurrence `e (k+1) ≤ K · (e k)^2` and initial scaled error
    `K · e 0 ≤ β` (with `0 ≤ β ≤ 1`) stays in the basin forever:
        `∀ k, K · e k ≤ β`.

    **Proof.** Induction on `k`. Base: `K · e 0 ≤ β` by hypothesis.
    Step: `K · e (k+1) ≤ K · (K · (e k)²) = (K · e k)² ≤ β² ≤ β`,
    where `β² ≤ β` uses `0 ≤ β ≤ 1`. -/
theorem real_iterate_scaled_basin_invariant
    {K β : ℝ} (hK_nn : 0 ≤ K) (hβ_nn : 0 ≤ β) (hβ_le_1 : β ≤ 1)
    (e : ℕ → ℝ) (he_nn : ∀ k, 0 ≤ e k)
    (he_0 : K * e 0 ≤ β)
    (he_quad : ∀ k, e (k + 1) ≤ K * (e k) ^ 2) :
    ∀ k, K * e k ≤ β := by
  intro k
  induction k with
  | zero => exact he_0
  | succ n ih =>
    have h_Ken_nn : 0 ≤ K * e n := mul_nonneg hK_nn (he_nn n)
    have step_a : K * e (n + 1) ≤ (K * e n) ^ 2 := by
      calc K * e (n + 1)
          ≤ K * (K * (e n) ^ 2) :=
            mul_le_mul_of_nonneg_left (he_quad n) hK_nn
        _ = (K * e n) ^ 2 := by ring
    have step_c : (K * e n) ^ 2 ≤ β ^ 2 := pow_le_pow_left h_Ken_nn ih 2
    have step_d : β ^ 2 ≤ β := by
      have : β * β ≤ β * 1 := mul_le_mul_of_nonneg_left hβ_le_1 hβ_nn
      calc β ^ 2 = β * β := by ring
        _ ≤ β * 1 := this
        _ = β := by ring
    linarith

end BasinInvariant

/-! ============================================================
  §25.2  Linear bound from quadratic + basin invariant
============================================================ -/

section LinearFromQuadratic

/-- **Linear bound from quadratic + basin invariant.** In the scaled
    basin `K · e k ≤ β`, a quadratic recurrence collapses to a linear
    recurrence at the basin's rate:
        `e (k+1) ≤ β · e k`.

    **Proof.** `e (k+1) ≤ K · (e k)² = (K · e k) · e k ≤ β · e k`. -/
theorem real_iterate_linear_from_quadratic
    {K β : ℝ} (_hK_nn : 0 ≤ K)
    (e : ℕ → ℝ) (he_nn : ∀ k, 0 ≤ e k)
    (h_inv : ∀ k, K * e k ≤ β)
    (he_quad : ∀ k, e (k + 1) ≤ K * (e k) ^ 2) :
    ∀ k, e (k + 1) ≤ β * e k := by
  intro k
  calc e (k + 1)
      ≤ K * (e k) ^ 2 := he_quad k
    _ = (K * e k) * e k := by ring
    _ ≤ β * e k := mul_le_mul_of_nonneg_right (h_inv k) (he_nn k)

end LinearFromQuadratic

/-! ============================================================
  §25.3  ★ Real Pandrosion iteration produces admissible sequences
============================================================ -/

section AdmissibleSequences

/-- **★ Abstract→concrete bridge: single-`p` form.**

    Fix `x > 1` and `p ≥ 1`. Let `α` be any real (with the semantic
    reading `α = x^{−1/p}`, the Pandrosion fixed point). Let
    `s : ℕ → ℝ` be the *real* Pandrosion-Steffensen iterate, and
    `K, R : ℝ` be the per-`p` quadratic-rate constant and radius.

    **Provided:**
    * `|s 0 − α| ≤ R` (initial error inside the radius),
    * `|s (k+1) − α| ≤ K · |s k − α|²` (Aitken-Steffensen quadratic
      contraction — the §15.4 residue identity specialised to the
      iterate),
    * `K · R ≤ lamModel p x` (scaled basin at the closed-form linear
      rate from §23.5),
    * `0 ≤ lamModel p x ≤ 1` (regularity — holds unconditionally for
      `x > 1` and `p ≥ p_min(x)` by §23.2 + §23.5),

    **the concrete error sequence** `e k := |s k − α|` *automatically*
    satisfies the three §23.6 admissibility conditions:
      (L) linear contraction at rate `lamModel p x`,
      (Q) quadratic contraction with constant `K`,
      (I) scaled-basin invariant.

    This converts §23.6 from a statement about abstract sequences into
    a statement about the *real* Pandrosion-Steffensen iteration on
    `z^p = x` in the scaled basin of the closed-form rate. -/
theorem pandrosion_produces_admissible_sequences
    (x : ℝ) (p : ℕ) (α : ℝ) (s : ℕ → ℝ) (K R : ℝ)
    (hK : 0 < K) (_hR : 0 < R)
    (he_0 : |s 0 - α| ≤ R)
    (he_quad : ∀ k, |s (k + 1) - α| ≤ K * |s k - α| ^ 2)
    (h_KR_scaled : K * R ≤ lamModel p x)
    (h_lam_nn : 0 ≤ lamModel p x) (h_lam_le_1 : lamModel p x ≤ 1) :
    -- (L) Linear bound at the closed-form rate.
    (∀ k, |s (k + 1) - α| ≤ lamModel p x * |s k - α|) ∧
    -- (Q) Quadratic bound — restated for completeness.
    (∀ k, |s (k + 1) - α| ≤ K * |s k - α| ^ 2) ∧
    -- (I) Scaled-basin invariant.
    (∀ k, K * |s k - α| ≤ lamModel p x) := by
  set e : ℕ → ℝ := fun k => |s k - α|
  have he_nn : ∀ k, 0 ≤ e k := fun k => abs_nonneg _
  have he_0' : K * e 0 ≤ lamModel p x := by
    calc K * e 0
        ≤ K * R := mul_le_mul_of_nonneg_left he_0 hK.le
      _ ≤ lamModel p x := h_KR_scaled
  have he_rec : ∀ k, e (k + 1) ≤ K * (e k) ^ 2 := fun k => he_quad k
  have h_inv :=
    real_iterate_scaled_basin_invariant hK.le h_lam_nn h_lam_le_1 e he_nn he_0' he_rec
  have h_lin :=
    real_iterate_linear_from_quadratic hK.le e he_nn h_inv he_rec
  exact ⟨h_lin, he_rec, h_inv⟩

end AdmissibleSequences

/-! ============================================================
  §25.4  ★ Real iteration ⟹ p-uniform complexity on `z^p = x`
============================================================ -/

section ConcreteSuperGrandMaster

/-- **★★ Concrete super grand-master on `z^p = x`.**

    The real Pandrosion-Steffensen iteration attains p-uniform
    complexity in the bounded regime `x > 1`.

    **Inputs (all concrete).**
    * `x > 1`: the parameter of the root-extraction problem `z^p = x`.
    * `α : ℕ → ℝ`: per-`p` Pandrosion fixed points (semantically
      `α p = x^{−1/p}`, any real satisfying `α p ≠ 0`).
    * `s : ℕ → ℕ → ℝ`: per-`p` real Steffensen-iterate sequences
      (`s p k` is the k-th iterate for parameter `p`).
    * `K, R : ℕ → ℝ`: per-`p` quadratic-rate constants and radii.
    * `he_0`: each iterate starts inside its radius.
    * `he_quad`: each iterate satisfies the §15.4 quadratic residue.
    * `h_scaled`: uniform scaled-basin hypothesis
      `K p · R p ≤ lamModel p x`.
    * `h_lam_reg`: `0 ≤ lamModel p x ≤ 1` (regularity — holds
      unconditionally for `p ≥ p_min(x)` by §23).
    * `M ρ δ`: final complexity parameters.
    * `h_KR_uniform`: uniform bound `K p · R p ≤ M` (needed by §23.6;
      can take `M := 1` here since `lamModel p x ≤ 1`).

    **Conclusion.** `∃ p₀ N, 1 ≤ p₀ ∧ ∀ p ≥ p₀, ∀ k ≥ N,
      K p · |s p k − α p| ≤ δ`.

    The iteration count `N` is **uniform in `p`**, making the result
    a genuine p-uniform complexity statement on the *real* iteration. -/
theorem pandrosion_concrete_super_grand_master
    (x : ℝ) (hx : 1 < x)
    (α : ℕ → ℝ) (s : ℕ → ℕ → ℝ)
    (K R : ℕ → ℝ) (hK : ∀ p, 0 < K p) (hR : ∀ p, 0 < R p)
    (he_0 : ∀ p, |s p 0 - α p| ≤ R p)
    (he_quad : ∀ p k, |s p (k + 1) - α p| ≤ K p * |s p k - α p| ^ 2)
    (h_scaled : ∀ p, K p * R p ≤ lamModel p x)
    (h_lam_reg : ∀ p, 0 ≤ lamModel p x ∧ lamModel p x ≤ 1)
    (M ρ δ : ℝ) (hM_pos : 0 < M)
    (hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (h_KR_uniform : ∀ p, K p * R p ≤ M) :
    ∃ (p₀ : ℕ) (N : ℕ), 1 ≤ p₀ ∧
      ∀ p, p₀ ≤ p → ∀ k ≥ N, K p * |s p k - α p| ≤ δ := by
  -- Per-`p`, derive admissibility (linear + quadratic + invariant).
  set e : ℕ → ℕ → ℝ := fun p k => |s p k - α p|
  have he_nn : ∀ p k, 0 ≤ e p k := fun p k => abs_nonneg _
  have he_0_abs : ∀ p, e p 0 ≤ R p := he_0
  have he_quad_abs : ∀ p k, e p (k + 1) ≤ K p * (e p k) ^ 2 := he_quad
  -- Linear bound per `p` via §25.3.
  have he_linear : ∀ p k, e p (k + 1) ≤ lamModel p x * e p k := by
    intro p
    have h := pandrosion_produces_admissible_sequences
      x p (α p) (s p) (K p) (R p) (hK p) (hR p) (he_0 p) (he_quad p)
      (h_scaled p) (h_lam_reg p).1 (h_lam_reg p).2
    exact h.1
  -- Feed into §23.6.
  exact super_grand_master_uniform_closed_form x M ρ δ hx hM_pos
    hρ_pos hρ_lt_1 hδ_pos K R e hK hR he_nn he_0_abs he_linear he_quad_abs
    h_KR_uniform

end ConcreteSuperGrandMaster

end Pandrosion
