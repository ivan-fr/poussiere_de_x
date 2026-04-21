/-
  Universitas Pandrosion — §21. Uniform log-log complexity in `p`.

  The log-log complexity theorem `quadratic_loglog_complexity` (§17.2)
  produces, for a SINGLE quadratic-contraction sequence with constants
  `(K, r, ε)`, an iteration count after which the error drops below ε.

  For the Pandrosion-Steffensen iteration on `z^p − x`, the quadratic
  constant `K_p = (p−1)/(2α) · (1 − λ_fp)` grows LINEARLY in p, and the
  pre-quadratic contraction rate `λ_fp` approaches 1 as p → ∞. Hence
  absolute precision `|s_k − α| < ε` is NOT achievable in a p-uniform
  number of iterations: the existential `N` of §17.2 depends on p.

  However, in the SCALED coordinate `ê_k := K_p · e_k`, the recurrence
  `|e_{k+1}| ≤ K_p · |e_k|²` becomes
      `ê_{k+1} ≤ ê_k²`,
  which is **p-INDEPENDENT**. The number of iterations to drive the
  scaled error below a fixed scaled-target `δ > 0` from a fixed
  scaled-initial `ρ < 1` is therefore bounded by a universal constant.

  This file formalises that insight:

    §21.1  `uniform_quadratic_loglog` — for any family of
           quadratic-contraction sequences `(K_p, r_p, ε_p, e_p)_p`
           with uniform scaled bounds
              `K_p · r_p ≤ ρ < 1`  and  `K_p · ε_p ≥ δ > 0`,
           there exists a p-UNIFORM iteration count `N` such that
           `e_p(k) ≤ ε_p` for every p and every `k ≥ N`.

    §21.2  `uniform_quadratic_loglog_scaled` — specialisation to the
           "scaled-precision" target `ε_p = δ/K_p`: drive the scaled
           error `K_p · e_p(k)` below `δ` in p-uniform steps.

    ★ §21.3  `pandrosion_uniform_loglog` — direct corollary for
           the Pandrosion-Steffensen family: under uniform scaled
           initial contraction `K_p · |s_p(0) − α_p| ≤ ρ < 1`, the
           scaled error `K_p · |s_p(k) − α_p|` reaches any fixed
           scaled-target `δ > 0` in a p-uniform number of iterations.

  **Interpretation.** The Newton constant `(p−1)/2` for the polynomial
  `u^p − 1` is the intrinsic curvature of a `p`-fold root; no iteration
  scheme can beat it. But the Pandrosion-Steffensen algorithm achieves
  that intrinsic rate in a **universally bounded** number of iterations
  per unit of scaled precision — this is the correct uniform statement.
-/

import Pandrosion.Core.QuadraticComplexity

namespace Pandrosion

/-! ============================================================
  §21.1  Uniform log-log complexity — core theorem
============================================================ -/

section UniformCore

/-- **Auxiliary: `n ≤ 2^n`.** Independent copy to avoid reaching into
    the private helper of §17.2. -/
private theorem nat_le_two_pow_self_u (n : ℕ) : n ≤ 2 ^ n := by
  induction n with
  | zero => simp
  | succ m ih =>
    have h2 : 0 < 2 ^ m := Nat.pos_pow_of_pos m (by norm_num)
    calc m + 1 ≤ 2 ^ m + 1 := by linarith
      _ ≤ 2 ^ m + 2 ^ m := by linarith
      _ = 2 ^ (m + 1) := by rw [pow_succ]; ring

/-- **★ Uniform log-log complexity across a family.**

    Let `(K_p, r_p, ε_p, e_p)_p` be a family of quadratic-contraction
    sequences:
      • `e_p(k+1) ≤ K_p · e_p(k)²`   (quadratic recurrence per p),
      • `e_p(0) ≤ r_p`                (initial bound per p),
    with uniform scaled bounds independent of p:
      • `K_p · r_p ≤ ρ < 1`           (uniform scaled contraction),
      • `K_p · ε_p ≥ δ > 0`           (uniform scaled target).

    Then there exists an iteration count `N : ℕ`, depending ONLY on
    `ρ` and `δ` (not on the parameter `p`), such that
      `e_p(k) ≤ ε_p`    for every `p` and every `k ≥ N`.

    **Proof idea.** The scaled error `ê_p(k) := K_p · e_p(k)` satisfies
    the p-independent recurrence `ê_p(k+1) ≤ ê_p(k)²`, with initial
    value `ê_p(0) ≤ ρ`. Hence `ê_p(k) ≤ ρ^{2^k}` (power-tower bound
    §17.1), and the smallest `N` with `ρ^N < δ` already ensures
    `ê_p(k) < δ ≤ K_p · ε_p` for `k ≥ N` (using `N ≤ 2^N`). -/
theorem uniform_quadratic_loglog
    {ι : Type*}
    (ρ δ : ℝ) (_hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (K r : ι → ℝ) (ε : ι → ℝ) (e : ι → ℕ → ℝ)
    (hK : ∀ p, 0 < K p) (hr : ∀ p, 0 ≤ r p) (_hε : ∀ p, 0 < ε p)
    (he_nn : ∀ p k, 0 ≤ e p k)
    (he_0 : ∀ p, e p 0 ≤ r p)
    (he_rec : ∀ p k, e p (k + 1) ≤ K p * (e p k) ^ 2)
    (h_Kr_uniform : ∀ p, K p * r p ≤ ρ)
    (h_Kε_uniform : ∀ p, δ ≤ K p * ε p) :
    ∃ N : ℕ, ∀ p, ∀ k ≥ N, e p k ≤ ε p := by
  -- Pick the smallest N with ρ^N < δ (exists since ρ < 1, δ > 0).
  obtain ⟨N, hN⟩ : ∃ n, ρ ^ n < δ :=
    exists_pow_lt_of_lt_one hδ_pos hρ_lt_1
  refine ⟨N, fun p k hk => ?_⟩
  have hK_pos := hK p
  have hKr_nn : 0 ≤ K p * r p := mul_nonneg hK_pos.le (hr p)
  have hKr_lt_1 : K p * r p < 1 := lt_of_le_of_lt (h_Kr_uniform p) hρ_lt_1
  -- Power-tower bound from §17.1.
  have h_tower := quadratic_tower_bound (K p) (r p) hK_pos (hr p) (e p)
                  (he_nn p) (he_0 p) (he_rec p) k
  -- N ≤ 2^N ≤ 2^k.
  have h_N_le_2N : N ≤ 2 ^ N := nat_le_two_pow_self_u N
  have h_2N_le_2k : 2 ^ N ≤ 2 ^ k :=
    Nat.pow_le_pow_right (by norm_num) hk
  have h_2k_ge_N : N ≤ 2 ^ k := le_trans h_N_le_2N h_2N_le_2k
  -- (K_p · r_p)^{2^k} ≤ (K_p · r_p)^N ≤ ρ^N < δ ≤ K_p · ε_p.
  have h_mono : (K p * r p) ^ (2 ^ k) ≤ (K p * r p) ^ N :=
    pow_le_pow_of_le_one hKr_nn (le_of_lt hKr_lt_1) h_2k_ge_N
  have h_rho : (K p * r p) ^ N ≤ ρ ^ N :=
    pow_le_pow_left hKr_nn (h_Kr_uniform p) N
  have h_bound : K p * e p k ≤ K p * ε p :=
    calc K p * e p k ≤ (K p * r p) ^ (2 ^ k) := h_tower
      _ ≤ (K p * r p) ^ N := h_mono
      _ ≤ ρ ^ N := h_rho
      _ ≤ δ := le_of_lt hN
      _ ≤ K p * ε p := h_Kε_uniform p
  exact le_of_mul_le_mul_left h_bound hK_pos

end UniformCore

/-! ============================================================
  §21.2  Scaled-precision corollary
============================================================ -/

section ScaledPrecision

/-- **★ Uniform scaled-precision log-log complexity.**

    Specialisation of `uniform_quadratic_loglog` to the scaled-target
    regime `ε_p := δ / K_p`: drive the *scaled* error `K_p · e_p(k)`
    below a p-independent threshold `δ > 0` in a p-UNIFORM number of
    iterations.

    This is the sharpest uniform statement available: the Pandrosion-
    Steffensen quadratic constant `K_p` grows with p, so absolute
    precision `e_p(k) < ε` requires `N ~ log log(1/(K_p · ε))` steps
    — a mild p-dependence inside a double logarithm. But in scaled
    coordinates, the complexity is completely p-free. -/
theorem uniform_quadratic_loglog_scaled
    {ι : Type*}
    (ρ δ : ℝ) (hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (K : ι → ℝ) (r : ι → ℝ) (e : ι → ℕ → ℝ)
    (hK : ∀ p, 0 < K p) (hr : ∀ p, 0 ≤ r p)
    (he_nn : ∀ p k, 0 ≤ e p k)
    (he_0 : ∀ p, e p 0 ≤ r p)
    (he_rec : ∀ p k, e p (k + 1) ≤ K p * (e p k) ^ 2)
    (h_Kr_uniform : ∀ p, K p * r p ≤ ρ) :
    ∃ N : ℕ, ∀ p, ∀ k ≥ N, K p * e p k ≤ δ := by
  -- Take ε_p := δ / K_p; then K_p · ε_p = δ ≥ δ.
  set ε : ι → ℝ := fun p => δ / K p
  have hε_pos : ∀ p, 0 < ε p := fun p => div_pos hδ_pos (hK p)
  have h_Kε_uniform : ∀ p, δ ≤ K p * ε p := fun p => by
    show δ ≤ K p * (δ / K p)
    rw [mul_div_assoc', mul_div_cancel_left₀ _ (hK p).ne']
  obtain ⟨N, hN⟩ := uniform_quadratic_loglog ρ δ hρ_pos hρ_lt_1 hδ_pos
    K r ε e hK hr hε_pos he_nn he_0 he_rec h_Kr_uniform h_Kε_uniform
  refine ⟨N, fun p k hk => ?_⟩
  have h_ek : e p k ≤ ε p := hN p k hk
  have h_mul : K p * e p k ≤ K p * ε p :=
    mul_le_mul_of_nonneg_left h_ek (hK p).le
  have h_eq : K p * ε p = δ := by
    show K p * (δ / K p) = δ
    rw [mul_div_assoc', mul_div_cancel_left₀ _ (hK p).ne']
  linarith

end ScaledPrecision

/-! ============================================================
  §21.3  Pandrosion-Steffensen specialisation
============================================================ -/

section PandrosionUniform

/-- **★ Uniform log-log complexity for the Pandrosion-Steffensen
    family.**

    Let `(s_p, α_p)_p` be a family of Pandrosion-Steffensen iterations
    (one per degree `p`) with p-dependent quadratic constants `K_p`:

      `|s_p(k+1) − α_p| ≤ K_p · |s_p(k) − α_p|²`    (per p).

    Under the UNIFORM scaled-contraction hypothesis on the initial
    guess,

      `K_p · |s_p(0) − α_p| ≤ ρ < 1`                (uniform in p),

    the scaled error `K_p · |s_p(k) − α_p|` reaches any target
    `δ > 0` in a p-UNIFORM number of iterations `N = N(ρ, δ)`.

    **Consequence.** Absolute precision `ε_p := δ / K_p` (i.e. absolute
    precision that tracks the intrinsic scale of a degree-p polynomial
    root) is achieved in a bounded number of iterations independent of
    p. This is the correct uniform form of the multi-start grand-master
    log-log complexity. -/
theorem pandrosion_uniform_loglog
    (ρ δ : ℝ) (hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (K : ℕ → ℝ) (s_star : ℕ → ℝ) (s : ℕ → ℕ → ℝ)
    (hK : ∀ p, 0 < K p)
    (he_rec : ∀ p k, |s p (k + 1) - s_star p| ≤ K p * |s p k - s_star p| ^ 2)
    (he_init : ∀ p, K p * |s p 0 - s_star p| ≤ ρ) :
    ∃ N : ℕ, ∀ p, ∀ k ≥ N, K p * |s p k - s_star p| ≤ δ := by
  set e : ℕ → ℕ → ℝ := fun p k => |s p k - s_star p|
  set r : ℕ → ℝ := fun p => |s p 0 - s_star p|
  have he_nn : ∀ p k, 0 ≤ e p k := fun _ _ => abs_nonneg _
  have hr_nn : ∀ p, 0 ≤ r p := fun _ => abs_nonneg _
  have he_0 : ∀ p, e p 0 ≤ r p := fun _ => le_refl _
  have he_rec' : ∀ p k, e p (k + 1) ≤ K p * (e p k) ^ 2 := he_rec
  have h_Kr_uniform : ∀ p, K p * r p ≤ ρ := he_init
  exact uniform_quadratic_loglog_scaled ρ δ hρ_pos hρ_lt_1 hδ_pos
    K r e hK hr_nn he_nn he_0 he_rec' h_Kr_uniform

end PandrosionUniform

end Pandrosion
