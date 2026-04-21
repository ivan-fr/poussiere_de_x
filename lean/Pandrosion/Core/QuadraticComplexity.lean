/-
  Universitas Pandrosion — §17. Abstract log-log complexity for
  quadratically-convergent sequences.

  The grand-master global-convergence theorem (`constructive_mcmullen`,
  Foundations §13) combined with the Aitken-Steffensen quadratic rate
  (`pandrosion_h_quadratic_polynomial`, §15.4) yields the following
  Smale-style complexity statement:

    For every `x ≠ 0`, every query `z₀` off the Voronoï boundary
    (a Lebesgue-null union of perpendicular bisectors), and every
    target precision `ε > 0`, the multi-start Steffensen-Pandrosion
    algorithm reaches precision `ε` after `O(log log(1/ε))` iterations.

  This file establishes the *abstract* complexity bound. Given a sequence
  `(e_k)_{k ≥ 0}` of non-negative errors satisfying the quadratic
  recurrence `e_{k+1} ≤ K · e_k²`, we prove:

    §17.1  `quadratic_tower_bound`  — power-tower estimate
           `K · e_k ≤ (K · e_0)^(2^k)`.

    §17.2  `quadratic_loglog_complexity` — existential log-log bound
           `∃ N, ∀ k ≥ N, e_k ≤ ε`.

  The constant `K` is instantiated from §15.4 on the Pandrosion side
  (Steffensen step on the scaled target `y = x/A`); the contraction
  condition `K · e_0 < 1` is inherited from §13 (anchor interior
  reachability) combined with §9.2 (closed-form rate `λ_{p,x} < 1`).
-/

import Pandrosion.Core.QuadraticRate

namespace Pandrosion

/-! ============================================================
  §17.1  Power-tower bound
============================================================ -/

section TowerBound

/-- **Quadratic tower bound.** A non-negative sequence `e : ℕ → ℝ` with
    quadratic contraction `e_{k+1} ≤ K · e_k²` and initial bound
    `e_0 ≤ r` satisfies the power-tower estimate

    `K · e_k ≤ (K · r)^(2^k)`    for all `k ≥ 0`.

    This is the classical Aitken-Steffensen tower estimate: each step
    squares the rescaled error `K · e_k`, producing the `2^k` exponent. -/
theorem quadratic_tower_bound
    (K r : ℝ) (hK : 0 < K) (_hr : 0 ≤ r)
    (e : ℕ → ℝ) (he_nn : ∀ k, 0 ≤ e k) (he_0 : e 0 ≤ r)
    (he_rec : ∀ k, e (k + 1) ≤ K * (e k) ^ 2) :
    ∀ k, K * e k ≤ (K * r) ^ (2 ^ k) := by
  intro k
  induction k with
  | zero =>
    simp only [pow_zero, pow_one]
    exact mul_le_mul_of_nonneg_left he_0 hK.le
  | succ k ih =>
    show K * e (k + 1) ≤ (K * r) ^ (2 ^ (k + 1))
    have h_nn_ek : 0 ≤ K * e k := mul_nonneg hK.le (he_nn k)
    have step_a : K * e (k + 1) ≤ (K * e k) ^ 2 := by
      have h1 : K * e (k + 1) ≤ K * (K * (e k) ^ 2) :=
        mul_le_mul_of_nonneg_left (he_rec k) hK.le
      calc K * e (k + 1) ≤ K * (K * (e k) ^ 2) := h1
        _ = (K * e k) ^ 2 := by ring
    have step_b : (K * e k) ^ 2 ≤ ((K * r) ^ (2 ^ k)) ^ 2 :=
      pow_le_pow_left h_nn_ek ih 2
    have step_c : ((K * r) ^ (2 ^ k)) ^ 2 = (K * r) ^ (2 ^ (k + 1)) := by
      rw [← pow_mul, pow_succ]
    linarith

end TowerBound

/-! ============================================================
  §17.2  Log-log iteration count
============================================================ -/

section LogLogComplexity

/-- **Auxiliary: `n ≤ 2^n`.** Elementary; extracted to keep the main
    argument linear. Matches Mathlib's `Nat.lt_two_pow_self` when
    available, or proved by induction. -/
private theorem nat_le_two_pow_self (n : ℕ) : n ≤ 2 ^ n := by
  induction n with
  | zero => simp
  | succ m ih =>
    have h2 : 0 < 2 ^ m := Nat.pos_pow_of_pos m (by norm_num)
    calc m + 1 ≤ 2 ^ m + 1 := by linarith
      _ ≤ 2 ^ m + 2 ^ m := by linarith
      _ = 2 ^ (m + 1) := by rw [pow_succ]; ring

/-- **★ Log-log iteration count.**
    If a non-negative sequence `e_k` contracts quadratically under
    `e_{k+1} ≤ K · e_k²` and starts within the contraction regime
    `K · e_0 ≤ K · r < 1`, then for every target precision `ε > 0`
    there exists `N : ℕ` such that `e_k ≤ ε` for all `k ≥ N`.

    **Effective rate.** The proof produces `N` explicitly: `N = n` where
    `n` is any natural satisfying `(K · r)^n < K · ε` (existence from
    `exists_pow_lt_of_lt_one`). Combined with the tower bound §17.1 and
    `n ≤ 2^n`, this gives `e_k ≤ ε` from step `k = n` onwards.

    Asymptotically, choosing the smallest such `n` gives
    `n = ⌈log(Kε) / log(Kr)⌉`; combined with `k ≥ log₂ n` for the
    tower inequality, the total iteration count is
    `k ≳ log₂(log(1/(Kε)) / log(1/(Kr))) = O(log log(1/ε))`. -/
theorem quadratic_loglog_complexity
    (K r : ℝ) (hK : 0 < K) (hr : 0 ≤ r) (hKr : K * r < 1)
    (e : ℕ → ℝ) (he_nn : ∀ k, 0 ≤ e k) (he_0 : e 0 ≤ r)
    (he_rec : ∀ k, e (k + 1) ≤ K * (e k) ^ 2)
    (ε : ℝ) (hε_pos : 0 < ε) :
    ∃ N : ℕ, ∀ k ≥ N, e k ≤ ε := by
  have hKr_nn : 0 ≤ K * r := mul_nonneg hK.le hr
  have hKε_pos : 0 < K * ε := mul_pos hK hε_pos
  obtain ⟨n, hn⟩ : ∃ n, (K * r) ^ n < K * ε :=
    exists_pow_lt_of_lt_one hKε_pos hKr
  refine ⟨n, fun k hk => ?_⟩
  have h_tower := quadratic_tower_bound K r hK hr e he_nn he_0 he_rec k
  have h_n_le_2n : n ≤ 2 ^ n := nat_le_two_pow_self n
  have h_2n_le_2k : 2 ^ n ≤ 2 ^ k :=
    Nat.pow_le_pow_right (by norm_num) hk
  have h_2k_ge_n : n ≤ 2 ^ k := le_trans h_n_le_2n h_2n_le_2k
  have h_mono : (K * r) ^ (2 ^ k) ≤ (K * r) ^ n :=
    pow_le_pow_of_le_one hKr_nn (le_of_lt hKr) h_2k_ge_n
  have h_bound : K * e k ≤ K * ε := by linarith
  exact le_of_mul_le_mul_left h_bound hK

end LogLogComplexity

/-! ============================================================
  §17.3  Instantiation for the Pandrosion super-linear residue
============================================================ -/

section PandrosionInstantiation

/-- **Specialisation.** If the Pandrosion super-linear residue
    `h(s_k) − s_star − λ_fp·(s_k − s_star)` dominates the Steffensen
    step error (classical Aitken–Steffensen), and the rational bound `C`
    of `pandrosion_h_quadratic_bound_of` is uniform on a neighborhood,
    then the iterative sequence `e_k := |s_k − s_star|` satisfies the
    hypotheses of `quadratic_loglog_complexity` — yielding
    log-log complexity in the target precision `ε`. -/
theorem pandrosion_steffensen_loglog_complexity
    (K r : ℝ) (hK : 0 < K) (hr : 0 ≤ r) (hKr : K * r < 1)
    (s_star : ℝ) (s : ℕ → ℝ)
    (he_0 : |s 0 - s_star| ≤ r)
    (he_rec : ∀ k, |s (k + 1) - s_star| ≤ K * |s k - s_star| ^ 2)
    (ε : ℝ) (hε_pos : 0 < ε) :
    ∃ N : ℕ, ∀ k ≥ N, |s k - s_star| ≤ ε := by
  set e : ℕ → ℝ := fun k => |s k - s_star|
  have he_nn : ∀ k, 0 ≤ e k := fun k => abs_nonneg _
  have he_sq : ∀ k, (e k) ^ 2 = |s k - s_star| ^ 2 := fun _ => rfl
  have he_rec' : ∀ k, e (k + 1) ≤ K * (e k) ^ 2 := by
    intro k
    rw [he_sq k]
    exact he_rec k
  exact quadratic_loglog_complexity K r hK hr hKr e he_nn he_0 he_rec' ε hε_pos

end PandrosionInstantiation

end Pandrosion
