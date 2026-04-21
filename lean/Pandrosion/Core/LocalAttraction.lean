/-
  Universitas Pandrosion — §27. **Local dynamics near a Pandrosion fixed point.**

  Dynamical upgrade of §15 quadratic-rate: under a local quadratic
  contraction bound on a disk `B(α, r)`, the iterates of a map
  converge to its fixed point `α` whenever the start `z₀` lies in
  the disk. This is the "local attraction" content of a
  super-attracting fixed point in the Julia-dynamics sense.

  Three layers:

    §27.1  `local_contraction_iterate_converges` — abstract iteration:
           a genuine linear contraction `‖f z − α‖ ≤ q·‖z − α‖` with
           `q < 1` on `B(α, r)` drives iterates to `α`.

    §27.2  `local_quadratic_iterate_converges` — reducer: a quadratic
           bound `‖f z − α‖ ≤ K·‖z − α‖²` on `B(α, r)` with `K·r ≤ 1/2`
           specializes to a linear contraction at rate `1/2`.

    §27.3  `steffensen_local_attraction` — specialization to
           `steffensen_step_C` at a Pandrosion fixed point, conditional
           on the quadratic bound (the super-attracting content of
           Steffensen's method, deferred here as a hypothesis).
-/

import Pandrosion.Core.Foundations
import Pandrosion.Core.ComplexMultiplier
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Topology.MetricSpace.Basic

namespace Pandrosion

open Filter Topology

/-! ============================================================
  §27.1  Abstract local linear contraction ⇒ convergence
============================================================ -/

section LocalContraction

/-- **Abstract local attraction.** If `f` fixes `α` and satisfies
    a linear contraction `‖f z − α‖ ≤ q·‖z − α‖` with `q ∈ [0,1)`
    on an open disk `B(α, r)`, then the iterates of `f` starting
    from any `z₀ ∈ B(α, r)` converge to `α`. -/
theorem local_contraction_iterate_converges
    {f : ℂ → ℂ} {α : ℂ} {r q : ℝ}
    (_hr : 0 < r) (hq_nn : 0 ≤ q) (hq : q < 1)
    (h_fp : f α = α)
    (h_contract : ∀ z : ℂ, ‖z - α‖ < r → ‖f z - α‖ ≤ q * ‖z - α‖) :
    ∀ z₀ : ℂ, ‖z₀ - α‖ < r →
      Tendsto (fun k => f^[k] z₀) atTop (𝓝 α) := by
  intro z₀ hz₀
  -- Step 1: inductively, iterates stay in the disk and shrink geometrically.
  have h_stays : ∀ k : ℕ,
      ‖f^[k] z₀ - α‖ ≤ q ^ k * ‖z₀ - α‖ ∧ ‖f^[k] z₀ - α‖ < r := by
    intro k
    induction k with
    | zero =>
        refine ⟨?_, hz₀⟩
        simp
    | succ n ih =>
        obtain ⟨h_bd, h_in⟩ := ih
        have h_step : ‖f (f^[n] z₀) - α‖ ≤ q * ‖f^[n] z₀ - α‖ :=
          h_contract (f^[n] z₀) h_in
        have h_next_bd : ‖f^[n + 1] z₀ - α‖ ≤ q ^ (n + 1) * ‖z₀ - α‖ := by
          rw [Function.iterate_succ_apply']
          calc ‖f (f^[n] z₀) - α‖
              ≤ q * ‖f^[n] z₀ - α‖ := h_step
            _ ≤ q * (q ^ n * ‖z₀ - α‖) :=
                mul_le_mul_of_nonneg_left h_bd hq_nn
            _ = q ^ (n + 1) * ‖z₀ - α‖ := by ring
        have h_next_in : ‖f^[n + 1] z₀ - α‖ < r := by
          calc ‖f^[n + 1] z₀ - α‖
              ≤ q ^ (n + 1) * ‖z₀ - α‖ := h_next_bd
            _ ≤ 1 * ‖z₀ - α‖ := by
                apply mul_le_mul_of_nonneg_right _ (norm_nonneg _)
                exact pow_le_one _ hq_nn (le_of_lt hq)
            _ = ‖z₀ - α‖ := one_mul _
            _ < r := hz₀
        exact ⟨h_next_bd, h_next_in⟩
  -- Step 2: convergence via `q^k → 0` and the geometric sandwich.
  rw [Metric.tendsto_atTop]
  intro ε hε
  -- Pick `N` with `q^N * ‖z₀ - α‖ < ε`.
  by_cases hz0_zero : ‖z₀ - α‖ = 0
  · -- `z₀ = α`: iterates are all `α` (since `α` is a fixed point).
    refine ⟨0, fun n _ => ?_⟩
    have hz0_eq : z₀ = α := by
      have := (norm_sub_eq_zero_iff (a := z₀) (b := α)).mp hz0_zero
      exact this
    have : f^[n] z₀ = α := by
      rw [hz0_eq]
      exact Function.iterate_fixed h_fp n
    rw [this, dist_self]
    exact hε
  · -- Generic case: use `tendsto_pow_atTop_nhds_zero_of_lt_one` scaled by `‖z₀ - α‖`.
    have h_pow_tendsto : Tendsto (fun k : ℕ => q ^ k * ‖z₀ - α‖) atTop (𝓝 0) := by
      have h_pow : Tendsto (fun k : ℕ => q ^ k) atTop (𝓝 0) :=
        tendsto_pow_atTop_nhds_zero_of_lt_one hq_nn hq
      have := h_pow.mul_const (‖z₀ - α‖)
      simpa using this
    rw [Metric.tendsto_atTop] at h_pow_tendsto
    obtain ⟨N, hN⟩ := h_pow_tendsto ε hε
    refine ⟨N, fun n hn => ?_⟩
    have h_bd := (h_stays n).1
    have h_pow_nn : 0 ≤ q ^ n * ‖z₀ - α‖ :=
      mul_nonneg (pow_nonneg hq_nn _) (norm_nonneg _)
    have h_dist_eq : dist (f^[n] z₀) α = ‖f^[n] z₀ - α‖ := dist_eq_norm _ _
    have h_small : q ^ n * ‖z₀ - α‖ < ε := by
      have := hN n hn
      rw [Real.dist_eq, sub_zero, abs_of_nonneg h_pow_nn] at this
      exact this
    calc dist (f^[n] z₀) α
        = ‖f^[n] z₀ - α‖ := h_dist_eq
      _ ≤ q ^ n * ‖z₀ - α‖ := h_bd
      _ < ε := h_small

end LocalContraction

/-! ============================================================
  §27.2  Local quadratic bound ⇒ linear contraction ⇒ convergence
============================================================ -/

section LocalQuadratic

/-- **Quadratic ⇒ linear contraction.** Given a quadratic error bound
    `‖f z − α‖ ≤ K·‖z − α‖²` on `B(α, r)` with `K·r ≤ 1/2` and `K > 0`,
    iterates converge: the factor `K·‖z − α‖ ≤ K·r ≤ 1/2` turns the
    quadratic bound into a linear contraction at rate `1/2`. -/
theorem local_quadratic_iterate_converges
    {f : ℂ → ℂ} {α : ℂ} {r K : ℝ}
    (hr : 0 < r) (hK_pos : 0 < K) (h_small : K * r ≤ 1 / 2)
    (h_fp : f α = α)
    (h_quadratic : ∀ z : ℂ, ‖z - α‖ < r → ‖f z - α‖ ≤ K * ‖z - α‖ ^ 2) :
    ∀ z₀ : ℂ, ‖z₀ - α‖ < r →
      Tendsto (fun k => f^[k] z₀) atTop (𝓝 α) := by
  apply local_contraction_iterate_converges hr
      (by norm_num : (0 : ℝ) ≤ 1 / 2) (by norm_num : (1 : ℝ) / 2 < 1) h_fp
  intro z hz
  calc ‖f z - α‖
      ≤ K * ‖z - α‖ ^ 2 := h_quadratic z hz
    _ = (K * ‖z - α‖) * ‖z - α‖ := by ring
    _ ≤ (K * r) * ‖z - α‖ := by
        apply mul_le_mul_of_nonneg_right _ (norm_nonneg _)
        exact mul_le_mul_of_nonneg_left (le_of_lt hz) (le_of_lt hK_pos)
    _ ≤ (1 / 2) * ‖z - α‖ :=
        mul_le_mul_of_nonneg_right h_small (norm_nonneg _)

end LocalQuadratic

/-! ============================================================
  §27.3  ★ steffensen_local_attraction — specialization
============================================================ -/

section SteffensenLocalAttraction

/-- **★★ Steffensen local attraction.** The Pandrosion-Steffensen
    iterator `steffensen_step_C x p` has a locally attracting basin
    around every fixed point `α` satisfying `α^p = 1/x`, provided a
    local quadratic contraction bound is available on a disk `B(α, r)`
    with `K · r ≤ 1/2`.

    The quadratic bound is the super-attracting content of Steffensen's
    method (`h'(α) = 0` ⇒ `|h(z) − α| = O(|z − α|²)` near `α`). It is
    taken as a hypothesis here, deferring the analytic proof (which
    requires handling the `if` branch / removable singularity of
    `steffensen_step_C` at `α`).

    Under this hypothesis, any start `z₀` in the disk yields iterates
    converging to `α`. This is the dynamical upgrade of §15
    quadratic-rate and the load-bearing local theorem behind the
    global multi-start basin result. -/
theorem steffensen_local_attraction
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (α : ℂ)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    {r K : ℝ} (hr : 0 < r) (hK_pos : 0 < K) (h_small : K * r ≤ 1 / 2)
    (h_quadratic : ∀ z : ℂ, ‖z - α‖ < r →
      ‖steffensen_step_C x p z - α‖ ≤ K * ‖z - α‖ ^ 2) :
    ∀ z₀ : ℂ, ‖z₀ - α‖ < r →
      Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop (𝓝 α) := by
  have h_fp : steffensen_step_C x p α = α :=
    steffensen_step_C_fixed_point x hx p α hSp hfp
  exact local_quadratic_iterate_converges hr hK_pos h_small h_fp h_quadratic

end SteffensenLocalAttraction

end Pandrosion
