/-
  Universitas Pandrosion — §22. Super multi-start grand-master
  (conditional two-phase uniform complexity).

  The uniform log-log complexity theorem `pandrosion_uniform_loglog`
  (§21.3) requires the initial error to already lie in the *scaled
  basin* `K_p · |s_p(0) − α_p| ≤ ρ < 1`. Outside this basin the
  quadratic recurrence `e_{k+1} ≤ K · e_k²` does NOT contract — it
  can grow.

  This file establishes a conditional "super" grand-master: given a
  joint linear+quadratic contraction structure
      `|s_p(k+1) − α_p| ≤ λ · |s_p(k) − α_p|`    (uniform linear),
      `|s_p(k+1) − α_p| ≤ K_p · |s_p(k) − α_p|²` (per-p quadratic),
  the two-phase analysis delivers a p-UNIFORM complexity bound for
  any family with UNIFORMLY bounded *scaled* initial error
  `K_p · R_p ≤ M` (arbitrary M, NOT required < 1).

  The theorems:

    §22.1  `linear_decay_bound` — `e_k ≤ λ^k · e_0` under linear
           contraction (standard).

    §22.2  `scaled_basin_entry_time` — entry into the scaled basin
           `K · e_k ≤ ρ` within `N_1 ≈ log(K·R/ρ)/log(1/λ)` steps.

    §22.3  `two_phase_complexity` — single-sequence total iteration
           count `N = N_1 + N_2`, linear transient + quadratic tail.

    ★ §22.4  `super_grand_master_uniform` — the p-UNIFORM two-phase
           bound. Under
                `λ ≤ λ_0 < 1`      (uniform linear rate),
                `K_p · R_p ≤ M`    (uniform scaled initial error,
                                    any bound M, not required < 1),
           the iteration count to reach scaled precision
           `K_p · e_p(k) ≤ δ` is bounded by a p-INDEPENDENT
           `N = N_1(λ_0, M, ρ) + N_2(ρ, δ)`.

  This EXTENDS `pandrosion_uniform_loglog` from "starting inside the
  scaled basin (K·R < 1)" to "starting anywhere with uniformly bounded
  scaled initial error (K·R ≤ M, any M)", absorbing the extra
  `O(log M / log(1/λ_0))` iterations into the p-uniform constant.

  **Caveat.** The uniform linear-contraction hypothesis `λ_0 < 1` is
  NOT automatic for Pandrosion-Steffensen on `z^p − x` over all
  `x ∈ ℂ \ {0}` — the closed-form rate `λ_{p,x}` approaches 1 as
  `p → ∞` when `x` is near 1 or as `|x| → ∞`. It holds uniformly in
  p in bounded regimes `x ∈ [x_min, x_max]` bounded away from 1 and
  ∞, which is a structural hypothesis about the parameter range, not
  a consequence of §13. This file provides the abstract two-phase
  machinery; applying it to a concrete x-regime is separate work.
-/

import Pandrosion.Core.UniformComplexity

namespace Pandrosion

/-! ============================================================
  §22.1  Linear decay bound
============================================================ -/

section LinearDecay

/-- **Geometric decay.** Under a linear-contraction recurrence
    `e_{k+1} ≤ λ · e_k` with `λ ≥ 0`, the iterate satisfies
    `e_k ≤ λ^k · e_0`. -/
theorem linear_decay_bound
    (lam : ℝ) (hlam_nn : 0 ≤ lam)
    (e : ℕ → ℝ) (_he_nn : ∀ k, 0 ≤ e k)
    (he_rec : ∀ k, e (k + 1) ≤ lam * e k) :
    ∀ k, e k ≤ lam ^ k * e 0 := by
  intro k
  induction k with
  | zero => simp
  | succ k ih =>
    calc e (k + 1) ≤ lam * e k := he_rec k
      _ ≤ lam * (lam ^ k * e 0) := mul_le_mul_of_nonneg_left ih hlam_nn
      _ = lam ^ (k + 1) * e 0 := by rw [pow_succ]; ring

end LinearDecay

/-! ============================================================
  §22.2  Scaled-basin entry time
============================================================ -/

section BasinEntry

/-- **Entry into the scaled basin.** Under linear contraction
    `e_{k+1} ≤ λ · e_k` with `λ < 1` and initial bound `e_0 ≤ R`,
    there exists `N_1 : ℕ` such that `K · e_{N_1} ≤ ρ`, i.e. the
    iterate enters the scaled basin `K · e ≤ ρ` by step `N_1`.

    The smallest such `N_1` is approximately
    `⌈log(K · R / ρ) / log(1/λ)⌉`. -/
theorem scaled_basin_entry_time
    (lam K R ρ : ℝ) (hlam_pos : 0 < lam) (hlam_lt_1 : lam < 1)
    (hK_pos : 0 < K) (hR_pos : 0 < R) (hρ_pos : 0 < ρ)
    (e : ℕ → ℝ) (he_nn : ∀ k, 0 ≤ e k) (he_0 : e 0 ≤ R)
    (he_rec : ∀ k, e (k + 1) ≤ lam * e k) :
    ∃ N1 : ℕ, K * e N1 ≤ ρ := by
  have hKR_pos : 0 < K * R := mul_pos hK_pos hR_pos
  have hρ_over_KR_pos : 0 < ρ / (K * R) := div_pos hρ_pos hKR_pos
  obtain ⟨N1, hN1⟩ : ∃ n, lam ^ n < ρ / (K * R) :=
    exists_pow_lt_of_lt_one hρ_over_KR_pos hlam_lt_1
  refine ⟨N1, ?_⟩
  have h_decay := linear_decay_bound lam hlam_pos.le e he_nn he_rec N1
  have hlam_N_nn : 0 ≤ lam ^ N1 := pow_nonneg hlam_pos.le N1
  have h_ek_R : e N1 ≤ lam ^ N1 * R :=
    le_trans h_decay (mul_le_mul_of_nonneg_left he_0 hlam_N_nn)
  calc K * e N1
      ≤ K * (lam ^ N1 * R) := mul_le_mul_of_nonneg_left h_ek_R hK_pos.le
    _ = lam ^ N1 * (K * R) := by ring
    _ ≤ (ρ / (K * R)) * (K * R) :=
        mul_le_mul_of_nonneg_right hN1.le hKR_pos.le
    _ = ρ := div_mul_cancel₀ _ hKR_pos.ne'

end BasinEntry

/-! ============================================================
  §22.3  Two-phase complexity — single sequence
============================================================ -/

section TwoPhaseSingle

/-- **Two-phase complexity (single sequence).**
    Under a joint linear+quadratic contraction structure
      `e_{k+1} ≤ λ · e_k`   (linear phase, transient)
      `e_{k+1} ≤ K · e_k²`  (quadratic phase, always valid),
    and any initial bound `e_0 ≤ R` (without requiring
    `K · R < 1`), the total iteration count to reach scaled
    precision `K · e_k ≤ δ` is finite, decomposing as
    `N = N_1 + N_2` where
      • `N_1` is the linear transient from `e_0 ≤ R` to
        `K · e_{N_1} ≤ ρ`,
      • `N_2` is the quadratic log-log from `K · e ≤ ρ` to
        `K · e ≤ δ` (§17.2). -/
theorem two_phase_complexity
    (lam K R ρ δ : ℝ) (hlam_pos : 0 < lam) (hlam_lt_1 : lam < 1)
    (hK_pos : 0 < K) (hR_pos : 0 < R)
    (hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (e : ℕ → ℝ) (he_nn : ∀ k, 0 ≤ e k) (he_0 : e 0 ≤ R)
    (he_linear : ∀ k, e (k + 1) ≤ lam * e k)
    (he_quad : ∀ k, e (k + 1) ≤ K * (e k)^2) :
    ∃ N : ℕ, ∀ k ≥ N, K * e k ≤ δ := by
  obtain ⟨N1, hN1⟩ := scaled_basin_entry_time lam K R ρ hlam_pos hlam_lt_1
    hK_pos hR_pos hρ_pos e he_nn he_0 he_linear
  set e' : ℕ → ℝ := fun k => e (N1 + k)
  have he'_nn : ∀ k, 0 ≤ e' k := fun _ => he_nn _
  have he'_0 : e' 0 ≤ ρ / K := by
    change e N1 ≤ ρ / K
    rw [le_div_iff hK_pos]
    linarith
  have he'_rec : ∀ k, e' (k + 1) ≤ K * (e' k) ^ 2 :=
    fun k => he_quad (N1 + k)
  have hρK_nn : 0 ≤ ρ / K := le_of_lt (div_pos hρ_pos hK_pos)
  have hK_ρK : K * (ρ / K) < 1 := by
    rw [mul_div_assoc', mul_div_cancel_left₀ _ hK_pos.ne']
    exact hρ_lt_1
  obtain ⟨N2, hN2⟩ := quadratic_loglog_complexity K (ρ / K) hK_pos hρK_nn hK_ρK
    e' he'_nn he'_0 he'_rec (δ / K) (div_pos hδ_pos hK_pos)
  refine ⟨N1 + N2, fun k hk => ?_⟩
  have hN1_le_k : N1 ≤ k := le_trans (Nat.le_add_right N1 N2) hk
  have h_diff : N2 ≤ k - N1 := by omega
  have h_eval : e' (k - N1) = e k := by
    show e (N1 + (k - N1)) = e k
    congr 1; omega
  have h_small : e' (k - N1) ≤ δ / K := hN2 (k - N1) h_diff
  rw [h_eval] at h_small
  calc K * e k ≤ K * (δ / K) := mul_le_mul_of_nonneg_left h_small hK_pos.le
    _ = δ := by rw [mul_div_assoc', mul_div_cancel_left₀ _ hK_pos.ne']

end TwoPhaseSingle

/-! ============================================================
  §22.4  ★ Super multi-start grand-master — uniform in p
============================================================ -/

section SuperGrandMaster

/-- **★ Super multi-start grand-master: uniform two-phase complexity.**

    Let `(K_p, R_p, e_p)_p` be a family of sequences with:
      • UNIFORM linear rate `λ_0 < 1`:
           `e_p(k+1) ≤ λ_0 · e_p(k)`    (∀p, ∀k),
      • per-p quadratic contraction:
           `e_p(k+1) ≤ K_p · e_p(k)²`   (∀p, ∀k),
      • UNIFORM scaled initial error:
           `K_p · R_p ≤ M`              (∀p)
        for some M > 0 (NOT required < 1 — may be arbitrary).

    Then there exists a p-UNIFORM iteration count `N` such that
    `K_p · e_p(k) ≤ δ` for every p and every `k ≥ N`.

    **Strict generalisation of `pandrosion_uniform_loglog` (§21.3).**
    §21 required `K_p · R_p ≤ ρ < 1` (inside the scaled basin). This
    theorem allows ANY uniform bound `K_p · R_p ≤ M`, using the
    linear phase to absorb the extra contraction needed to enter
    the scaled basin. The cost is `O(log M / log(1/λ_0))` additional
    iterations, uniform in p.

    **Structural content.** This is the closest one can get to a
    truly unconditional p-uniform grand-master with purely abstract
    methods. The remaining non-unconditionality is confined to the
    hypothesis `λ_0 < 1` uniform in p — a genuine constraint on the
    parameter regime of `z^p − x`, not an artefact of the proof. -/
theorem super_grand_master_uniform
    {ι : Type*}
    (lam_0 M ρ δ : ℝ)
    (hlam_pos : 0 < lam_0) (hlam_lt_1 : lam_0 < 1)
    (hM_pos : 0 < M)
    (hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (K : ι → ℝ) (R : ι → ℝ) (e : ι → ℕ → ℝ)
    (hK : ∀ p, 0 < K p) (_hR : ∀ p, 0 < R p)
    (he_nn : ∀ p k, 0 ≤ e p k)
    (he_0 : ∀ p, e p 0 ≤ R p)
    (he_linear : ∀ p k, e p (k + 1) ≤ lam_0 * e p k)
    (he_quad : ∀ p k, e p (k + 1) ≤ K p * (e p k)^2)
    (h_KR_uniform : ∀ p, K p * R p ≤ M) :
    ∃ N : ℕ, ∀ p, ∀ k ≥ N, K p * e p k ≤ δ := by
  -- Phase 1: p-uniform entry into the scaled basin.
  have hρ_over_M_pos : 0 < ρ / M := div_pos hρ_pos hM_pos
  obtain ⟨N1, hN1_lt⟩ : ∃ n, lam_0 ^ n < ρ / M :=
    exists_pow_lt_of_lt_one hρ_over_M_pos hlam_lt_1
  -- For every p, after N1 linear steps, K_p · e_p(N1) ≤ ρ (uniform!).
  have h_basin_uniform : ∀ p, K p * e p N1 ≤ ρ := by
    intro p
    have hlamN_nn : 0 ≤ lam_0 ^ N1 := pow_nonneg hlam_pos.le N1
    have h_decay := linear_decay_bound lam_0 hlam_pos.le (e p)
                    (he_nn p) (he_linear p) N1
    have h_ek_R : e p N1 ≤ lam_0 ^ N1 * R p :=
      le_trans h_decay (mul_le_mul_of_nonneg_left (he_0 p) hlamN_nn)
    calc K p * e p N1
        ≤ K p * (lam_0 ^ N1 * R p) :=
          mul_le_mul_of_nonneg_left h_ek_R (hK p).le
      _ = lam_0 ^ N1 * (K p * R p) := by ring
      _ ≤ lam_0 ^ N1 * M :=
          mul_le_mul_of_nonneg_left (h_KR_uniform p) hlamN_nn
      _ ≤ (ρ / M) * M :=
          mul_le_mul_of_nonneg_right hN1_lt.le hM_pos.le
      _ = ρ := div_mul_cancel₀ _ hM_pos.ne'
  -- Phase 2: apply §21.2 (uniform_quadratic_loglog_scaled) to the
  -- tail family e'(p,k) = e(p, k+N1).
  set e' : ι → ℕ → ℝ := fun p k => e p (N1 + k)
  set r' : ι → ℝ := fun p => e p N1
  have he'_nn : ∀ p k, 0 ≤ e' p k := fun p _ => he_nn p _
  have hr'_nn : ∀ p, 0 ≤ r' p := fun p => he_nn p N1
  have he'_0 : ∀ p, e' p 0 ≤ r' p := fun _ => le_refl _
  have he'_rec : ∀ p k, e' p (k + 1) ≤ K p * (e' p k) ^ 2 :=
    fun p k => he_quad p (N1 + k)
  have h_Kr'_uniform : ∀ p, K p * r' p ≤ ρ := h_basin_uniform
  obtain ⟨N2, hN2⟩ := uniform_quadratic_loglog_scaled ρ δ hρ_pos hρ_lt_1 hδ_pos
    K r' e' hK hr'_nn he'_nn he'_0 he'_rec h_Kr'_uniform
  refine ⟨N1 + N2, fun p k hk => ?_⟩
  have hN1_le_k : N1 ≤ k := le_trans (Nat.le_add_right N1 N2) hk
  have h_diff : N2 ≤ k - N1 := by omega
  have h_eval : e' p (k - N1) = e p k := by
    show e p (N1 + (k - N1)) = e p k
    congr 1; omega
  have h_small : K p * e' p (k - N1) ≤ δ := hN2 p (k - N1) h_diff
  rw [h_eval] at h_small
  exact h_small

end SuperGrandMaster

end Pandrosion
