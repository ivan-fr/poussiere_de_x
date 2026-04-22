/-
  Universitas Pandrosion — §53. **Universal unconditional orbit
  convergence on `[0, 2α]` for all `p ≥ 2`, `x > 1`.**

  Upgrades the Legacy `contraction_general` (single-step contraction,
  ∀ p ≥ 2, ∀ x > 1, ∀ s ≥ 0: `|h_{p,x}(s) − α| ≤ (x−1)/x · |s − α|`)
  to a **full orbit convergence theorem**:

      For every `p ≥ 2`, every `x > 1` real, every `s₀ ∈ [0, 2α]`
      real, the Pandrosion iteration `h_{p,x}^[k](s₀) → α`
      (where `α = x^{-1/p}` is the unique positive real root),
      with explicit linear rate `λ = (x−1)/x < 1`:
          `|h^[k](s₀) − α| ≤ λ^k · |s₀ − α|`.

  The orbit restriction `s₀ ∈ [0, 2α]` ensures each iterate stays
  in this "contraction ball" (by the linear bound itself), avoiding
  edge-case concerns about Pandrosion-h rational-form non-negativity.
  For `x > 1`, `α = x^{-1/p} ∈ (0, 1)`, so the ball `[0, 2α] ⊂ [0, 2)`
  is a concrete compact subset of ℝ₊.

  **Inconditionnel** : aucune hypothèse McMullen / Julia / fractal.
  **Ouvre `p ≥ 3` sans discharge complexe-dynamique.**

  Contents.

    §53.1  `pandrosion_h_iterate_in_ball` — orbit stays in `[0, 2α]`.

    §53.2  `pandrosion_h_orbit_linear_rate` — linear rate bound
           `|h^[k](s₀) − α| ≤ λ^k · |s₀ − α|`.

    §53.3  `pandrosion_h_real_positive_converges` — orbit
           `Tendsto` to `α` as `k → ∞`.
-/

import Pandrosion.Core.Foundations
import Pandrosion.Legacy

namespace Pandrosion

open Filter Topology

/-! ============================================================
  §53.1  Iterate stays in the contraction ball `[0, 2α]`
============================================================ -/

section IterateInBallC

/-- **Contraction ball invariance.**

    For `x > 1`, `α ∈ (0, 1)` with `α^p = 1/x`, and any `s₀ ∈ [0, 2α]`,
    the entire forward orbit `h^[k](s₀)` stays in `[0, 2α]`.

    *Proof.* Induction on `k`. The base case is `hs₀ : s₀ ∈ [0, 2α]`.
    Inductive step: if `s_k ∈ [0, 2α]`, then `|s_k − α| ≤ α`, so by
    `contraction_general`, `|h(s_k) − α| ≤ (x−1)/x · α < α`, hence
    `h(s_k) ∈ (0, 2α) ⊂ [0, 2α]`. -/
theorem pandrosion_h_iterate_in_ball
    (p : ℕ) (hp : p ≥ 2) (x : ℝ) (hx : x > 1)
    (α : ℝ) (hα_pos : α > 0) (hα_lt : α < 1)
    (hα_eq : α ^ p = 1 / x)
    (s₀ : ℝ) (hs₀_nn : s₀ ≥ 0) (hs₀_le : s₀ ≤ 2 * α)
    (k : ℕ) :
    (pandrosion_h x p)^[k] s₀ ≥ 0 ∧ (pandrosion_h x p)^[k] s₀ ≤ 2 * α := by
  induction k with
  | zero => simpa using ⟨hs₀_nn, hs₀_le⟩
  | succ j ih =>
    obtain ⟨ih_nn, ih_le⟩ := ih
    rw [Function.iterate_succ_apply']
    set s_j := (pandrosion_h x p)^[j] s₀
    -- |s_j - α| ≤ α (from ball).
    have h_dist_sj : |s_j - α| ≤ α := by
      rw [abs_sub_le_iff]
      refine ⟨?_, ?_⟩ <;> linarith
    -- Apply contraction.
    by_cases h_s_eq : s_j = α
    · -- s_j = α: h(α) = α, already in ball.
      have h_fix : pandrosion_h x p α = α := by
        have hx_pos : x > 0 := by linarith
        exact (pandrosion_fixed_point_iff x hx_pos p (by omega) α
          (le_of_lt hα_pos)).mpr hα_eq
      rw [h_s_eq, h_fix]
      refine ⟨le_of_lt hα_pos, ?_⟩
      linarith
    · -- s_j ≠ α: apply contraction_general.
      have h_contract := contraction_general p hp x α hx hα_pos hα_lt
        hα_eq s_j ih_nn h_s_eq
      -- |h(s_j) - α| ≤ (x-1)/x · |s_j - α| ≤ (x-1)/x · α < α.
      have hxm1_nn : 0 ≤ (x - 1) / x := by
        apply div_nonneg
        · linarith
        · linarith
      have h_rate_lt_one : (x - 1) / x < 1 := by
        rw [div_lt_one (by linarith : (0:ℝ) < x)]; linarith
      have h_contract_bd : |pandrosion_h x p s_j - α| ≤ (x - 1) / x * α := by
        calc |pandrosion_h x p s_j - α|
            ≤ (x - 1) / x * |s_j - α| := h_contract
          _ ≤ (x - 1) / x * α :=
              mul_le_mul_of_nonneg_left h_dist_sj hxm1_nn
      have h_next_lt : |pandrosion_h x p s_j - α| < α := by
        calc |pandrosion_h x p s_j - α| ≤ (x - 1) / x * α := h_contract_bd
          _ < 1 * α := by
              apply mul_lt_mul_of_pos_right h_rate_lt_one hα_pos
          _ = α := by ring
      rw [abs_sub_lt_iff] at h_next_lt
      refine ⟨?_, ?_⟩ <;> linarith

end IterateInBallC

/-! ============================================================
  §53.2  Linear rate orbit bound
============================================================ -/

section LinearRateOrbitC

/-- **★★ Universal orbit linear-rate bound.**

    For every `p ≥ 2`, `x > 1`, positive fixed point `α` with
    `α^p = 1/x`, and starting point `s₀ ∈ [0, 2α]`:
        `|h^[k](s₀) − α| ≤ ((x−1)/x)^k · |s₀ − α|`  for every `k`.

    *Proof.* Induction on `k` applying `contraction_general` to each
    iterate (which stays in `[0, 2α]` by §53.1). -/
theorem pandrosion_h_orbit_linear_rate
    (p : ℕ) (hp : p ≥ 2) (x : ℝ) (hx : x > 1)
    (α : ℝ) (hα_pos : α > 0) (hα_lt : α < 1)
    (hα_eq : α ^ p = 1 / x)
    (s₀ : ℝ) (hs₀_nn : s₀ ≥ 0) (hs₀_le : s₀ ≤ 2 * α)
    (k : ℕ) :
    |(pandrosion_h x p)^[k] s₀ - α| ≤ ((x - 1) / x) ^ k * |s₀ - α| := by
  induction k with
  | zero =>
    simp
  | succ j ih =>
    rw [Function.iterate_succ_apply']
    set s_j := (pandrosion_h x p)^[j] s₀
    have ⟨hs_j_nn, _hs_j_le⟩ :=
      pandrosion_h_iterate_in_ball p hp x hx α hα_pos hα_lt hα_eq s₀
        hs₀_nn hs₀_le j
    have hx_nn : (0 : ℝ) ≤ x := by linarith
    have hxm1_nn_raw : (0 : ℝ) ≤ x - 1 := by linarith
    have hxm1_nn : (0 : ℝ) ≤ (x - 1) / x := div_nonneg hxm1_nn_raw hx_nn
    -- Apply contraction_general to the inner iterate.
    by_cases h_s_eq : s_j = α
    · -- s_j = α: h(α) = α; both sides give 0 ≤ something.
      have hx_pos : x > 0 := by linarith
      have h_fix : pandrosion_h x p α = α :=
        (pandrosion_fixed_point_iff x hx_pos p (by omega) α
          (le_of_lt hα_pos)).mpr hα_eq
      rw [h_s_eq, h_fix]
      simp only [sub_self, abs_zero]
      apply mul_nonneg (pow_nonneg hxm1_nn _) (abs_nonneg _)
    · have h_contract := contraction_general p hp x α hx hα_pos hα_lt
        hα_eq s_j hs_j_nn h_s_eq
      calc |pandrosion_h x p s_j - α|
          ≤ (x - 1) / x * |s_j - α| := h_contract
        _ ≤ (x - 1) / x * (((x - 1) / x) ^ j * |s₀ - α|) :=
            mul_le_mul_of_nonneg_left ih hxm1_nn
        _ = ((x - 1) / x) ^ (j + 1) * |s₀ - α| := by
            rw [pow_succ]; ring

end LinearRateOrbitC

/-! ============================================================
  §53.3  Orbit convergence
============================================================ -/

section OrbitConvergesC

/-- **★★★★ Universal unconditional Pandrosion convergence on ℝ₊
    for every `p ≥ 2`.**

    For every `p ≥ 2`, every `x > 1` real, every `α ∈ (0, 1)` with
    `α^p = 1/x` (the unique positive real p-th root of `1/x`), and
    every starting point `s₀ ∈ [0, 2α]`:
        `(pandrosion_h x p)^[k] s₀ → α`  as  `k → ∞`.

    **This is the p ≥ 3 unconditional real-line convergence theorem**
    obtained via pure linear contraction — no McMullen / Julia /
    complex-dynamics machinery. Works for cube roots (`p = 3`) and
    all higher roots. Opens the door to unconditional p ≥ 3 analysis
    on the real positive axis. -/
theorem pandrosion_h_real_positive_converges
    (p : ℕ) (hp : p ≥ 2) (x : ℝ) (hx : x > 1)
    (α : ℝ) (hα_pos : α > 0) (hα_lt : α < 1)
    (hα_eq : α ^ p = 1 / x)
    (s₀ : ℝ) (hs₀_nn : s₀ ≥ 0) (hs₀_le : s₀ ≤ 2 * α) :
    Tendsto (fun k => (pandrosion_h x p)^[k] s₀) atTop (𝓝 α) := by
  -- Squeeze via |h^[k] s₀ − α| ≤ λ^k · |s₀ − α| → 0.
  rw [Metric.tendsto_atTop]
  intro ε hε_pos
  have h_rate_pos : 0 < (x - 1) / x := by
    apply div_pos
    · linarith
    · linarith
  have h_rate_lt_one : (x - 1) / x < 1 := by
    rw [div_lt_one (by linarith : (0:ℝ) < x)]; linarith
  -- Find N such that λ^N · |s₀ − α| < ε.
  have h_tend_pow : Tendsto (fun k : ℕ => ((x - 1) / x) ^ k) atTop (𝓝 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (le_of_lt h_rate_pos) h_rate_lt_one
  by_cases hs₀_eq : s₀ = α
  · -- Trivial: starting at fixed point.
    refine ⟨0, fun k _hk => ?_⟩
    have hx_pos : x > 0 := by linarith
    have h_fix : pandrosion_h x p α = α :=
      (pandrosion_fixed_point_iff x hx_pos p (by omega) α
        (le_of_lt hα_pos)).mpr hα_eq
    rw [hs₀_eq, Function.iterate_fixed h_fix]
    simpa using hε_pos
  · have h_dist_pos : 0 < |s₀ - α| := abs_pos.mpr (sub_ne_zero.mpr hs₀_eq)
    -- Extract N from Tendsto: λ^N < ε/|s₀ - α|.
    set target : ℝ := ε / |s₀ - α|
    have h_target_pos : 0 < target := div_pos hε_pos h_dist_pos
    have h_event : ∀ᶠ k in atTop, ((x - 1) / x) ^ k < target := by
      have h_mem : Set.Iio target ∈ 𝓝 (0 : ℝ) := Iio_mem_nhds h_target_pos
      have h_ev := h_tend_pow h_mem
      filter_upwards [h_ev] with k hk using hk
    obtain ⟨N, hN⟩ := h_event.exists_forall_of_atTop
    refine ⟨N, fun k hk => ?_⟩
    have h_orbit := pandrosion_h_orbit_linear_rate p hp x hx α hα_pos hα_lt
      hα_eq s₀ hs₀_nn hs₀_le k
    have h_ratek_lt : ((x - 1) / x) ^ k < target := hN k hk
    rw [dist_eq_norm]
    show ‖(pandrosion_h x p)^[k] s₀ - α‖ < ε
    have h_norm_eq : ‖(pandrosion_h x p)^[k] s₀ - α‖
                   = |(pandrosion_h x p)^[k] s₀ - α| := Real.norm_eq_abs _
    rw [h_norm_eq]
    calc |(pandrosion_h x p)^[k] s₀ - α|
        ≤ ((x - 1) / x) ^ k * |s₀ - α| := h_orbit
      _ < target * |s₀ - α| := by
          apply mul_lt_mul_of_pos_right h_ratek_lt h_dist_pos
      _ = ε / |s₀ - α| * |s₀ - α| := rfl
      _ = ε := by field_simp

end OrbitConvergesC

end Pandrosion
