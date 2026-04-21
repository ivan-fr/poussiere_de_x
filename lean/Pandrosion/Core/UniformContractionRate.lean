/-
  Universitas Pandrosion — §23. Uniform-in-p contraction rate for
  Pandrosion-Steffensen on `z^p − x`.

  **The missing ingredient of §22.** The super multi-start grand-master
  `super_grand_master_uniform` (§22.4) is *conditional* on the hypothesis
  that the linear-contraction rate `λ_p := |h'(α_{p,x})|` of the
  Pandrosion map `h(s) = 1 − (x−1)/(x · S_p(s))` at its fixed point
  `α_{p,x} = x^{−1/p}` is uniformly bounded strictly below 1:
      `λ_p ≤ λ_0 < 1`    (uniform in p).

  This file discharges that hypothesis in bounded x-regimes `x > 1` by
  identifying the asymptotic limit of `λ_p` as `p → ∞`:

      `λ_∞(x) := 1 − log x / (x − 1)`.

  Numerical investigation (cf. `/tmp/pandrosion_uniform2.py`, EXP C+)
  confirms the convergence
      `|λ_{p,x} − λ_∞(x)| = O(1/p)`    (uniform in bounded x),
  with representative values `λ_∞(1.5) ≈ 0.189`, `λ_∞(2) ≈ 0.307`,
  `λ_∞(5) ≈ 0.598`, `λ_∞(10) ≈ 0.744`. For every `x > 1` the limit
  `λ_∞(x)` lies strictly in `(0, 1)`.

  **Content.**

    §23.1  `lamInfty` — definition `λ_∞(x) = 1 − log x / (x − 1)`.

    §23.2  `lamInfty_pos_of_one_lt`, `lamInfty_lt_one_of_one_lt` —
           the asymptotic rate lies in `(0, 1)` for every `x > 1`.
           Proved via Mathlib's `Real.log_pos` and
           `Real.log_lt_sub_one_of_pos`.

    §23.3  `uniform_contraction_from_rate_convergence` — abstract
           machinery: given a family `λ : ℕ → ℝ` with
           `|λ p − λ_∞| ≤ C / p`, there exists `p₀ : ℕ` such that
           for all `p ≥ p₀`, `λ p ≤ (1 + λ_∞) / 2 < 1`. This
           discharges the uniform-linear-rate hypothesis of §22.

    ★ §23.4  `super_grand_master_uniform_bounded_regime` — the
           p-unconditional super grand-master in a bounded x-regime.
           Combines §23.3 with `super_grand_master_uniform` (§22.4)
           to conclude that, for any `x > 1` with the Taylor input
           `|λ_{p,x} − λ_∞(x)| ≤ C / p` and uniform scaled initial
           error `K_p · R_p ≤ M`, the iteration count to reach
           scaled precision is bounded UNIFORMLY in `p`.

  **Scope of the Taylor input.** The bound `|λ_{p,x} − λ_∞(x)| ≤ C/p`
  is a first-order Taylor expansion of `x^{−1/p} = exp(−log x / p)`
  around `1/p = 0`, combined with the algebraic identity
      `h'(α) = [p (α − 1) + α (x − 1)] / [α (x − 1)]`
  (obtained by direct differentiation of `h(s) = 1 − (x−1)/(x · S_p(s))`
  at `s = α`). Both steps are classical real-analysis calculations;
  here they enter as a named hypothesis `rate_taylor` so the machinery
  is self-contained in Lean while the numerical-analysis input is
  cleanly isolated. This mirrors the style of §21/§22 where the
  per-p quadratic constants enter as data, not derived objects.
-/

import Pandrosion.Core.SuperGrandMaster
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Complex.Exponential

namespace Pandrosion

open Real

/-! ============================================================
  §23.1  The asymptotic contraction rate `λ_∞(x)`
============================================================ -/

section LamInfty

/-- **Asymptotic contraction rate** `λ_∞(x)` of Pandrosion-Steffensen
    on `z^p − x` as `p → ∞`:
        `λ_∞(x) = 1 − log x / (x − 1)`.

    Numerical values (cf. `/tmp/pandrosion_uniform2.py`, EXP C+):
    `λ_∞(1.5) ≈ 0.189`, `λ_∞(2) ≈ 0.307`, `λ_∞(5) ≈ 0.598`,
    `λ_∞(10) ≈ 0.744`. Strictly in `(0, 1)` for every `x > 1`. -/
noncomputable def lamInfty (x : ℝ) : ℝ :=
  1 - Real.log x / (x - 1)

end LamInfty

/-! ============================================================
  §23.2  `λ_∞(x) ∈ (0, 1)` for `x > 1`
============================================================ -/

section LamInftyBounds

/-- **Upper bound.** For every `x > 1`, the asymptotic rate satisfies
    `λ_∞(x) < 1`.

    **Proof.** `log x > 0` (Mathlib `Real.log_pos`) and `x − 1 > 0`,
    so `log x / (x − 1) > 0`, hence `1 − log x / (x − 1) < 1`. -/
theorem lamInfty_lt_one_of_one_lt (x : ℝ) (hx : 1 < x) :
    lamInfty x < 1 := by
  unfold lamInfty
  have h_log_pos : 0 < Real.log x := Real.log_pos hx
  have h_den_pos : 0 < x - 1 := by linarith
  have h_frac_pos : 0 < Real.log x / (x - 1) :=
    div_pos h_log_pos h_den_pos
  linarith

/-- **Lower bound.** For every `x > 1`, the asymptotic rate satisfies
    `0 < λ_∞(x)`.

    **Proof.** The sharp inequality `log x < x − 1` (Mathlib's
    `Real.log_lt_sub_one_of_pos`) gives `log x / (x − 1) < 1`, hence
    `1 − log x / (x − 1) > 0`. -/
theorem lamInfty_pos_of_one_lt (x : ℝ) (hx : 1 < x) :
    0 < lamInfty x := by
  unfold lamInfty
  have h_pos : 0 < x := by linarith
  have h_ne_one : x ≠ 1 := ne_of_gt hx
  have h_log_lt : Real.log x < x - 1 :=
    Real.log_lt_sub_one_of_pos h_pos h_ne_one
  have h_den_pos : 0 < x - 1 := by linarith
  have h_frac_lt_one : Real.log x / (x - 1) < 1 := by
    rw [div_lt_one h_den_pos]
    exact h_log_lt
  linarith

/-- **Combined bound.** For `x > 1`, `λ_∞(x) ∈ (0, 1)`. -/
theorem lamInfty_mem_Ioo (x : ℝ) (hx : 1 < x) :
    0 < lamInfty x ∧ lamInfty x < 1 :=
  ⟨lamInfty_pos_of_one_lt x hx, lamInfty_lt_one_of_one_lt x hx⟩

end LamInftyBounds

/-! ============================================================
  §23.3  Uniform contraction from rate-convergence
============================================================ -/

section UniformFromConvergence

/-- **Uniform linear-rate bound from rate-convergence.**

    Let `λ : ℕ → ℝ` be a family of linear-contraction rates (one per
    degree `p`) with the convergence bound
        `|λ p − λ_∞| ≤ C / p`  (for `p ≥ 1`),
    and `λ_∞ < 1`. Then for every `p` sufficiently large,
        `λ p ≤ (1 + λ_∞) / 2 < 1`.

    **Interpretation.** This discharges the uniform-linear-rate
    hypothesis `λ_p ≤ λ_0 < 1` of `super_grand_master_uniform` (§22.4).
    The new upper bound `λ_0 := (1 + λ_∞) / 2` is the midpoint of
    `λ_∞` and `1`, which is strictly less than 1 for any `λ_∞ < 1`.

    **Proof.** The gap `1 − λ_∞ > 0` is fixed. Taking `p₀` large enough
    that `C / p₀ ≤ (1 − λ_∞) / 2`, every `p ≥ p₀` satisfies
        `λ p ≤ λ_∞ + C/p ≤ λ_∞ + (1 − λ_∞)/2 = (1 + λ_∞)/2`. -/
theorem uniform_contraction_from_rate_convergence
    (lam_infty C : ℝ)
    (hlam_infty_lt_1 : lam_infty < 1)
    (hC_nn : 0 ≤ C)
    (p_min : ℕ)
    (lam : ℕ → ℝ)
    (h_taylor : ∀ p : ℕ, p_min ≤ p → 1 ≤ p → |lam p - lam_infty| ≤ C / (p : ℝ)) :
    ∃ p₀ : ℕ, 1 ≤ p₀ ∧ p_min ≤ p₀ ∧ ∀ p : ℕ, p₀ ≤ p →
      lam p ≤ (1 + lam_infty) / 2 := by
  -- The contraction "budget" (gap from 1) we need to stay inside.
  have h_gap_pos : 0 < 1 - lam_infty := by linarith
  have _h_gap_half_pos : 0 < (1 - lam_infty) / 2 := by linarith
  -- Need p₀ so that C / p₀ ≤ (1 − λ_∞) / 2, i.e. p₀ ≥ 2C / (1 − λ_∞),
  -- AND p₀ ≥ p_min, AND p₀ ≥ 1.
  set B : ℝ := 2 * C / (1 - lam_infty)
  obtain ⟨m, hm⟩ := exists_nat_gt B
  refine ⟨max (max 1 p_min) m, ?_, ?_, fun p hp => ?_⟩
  · exact le_trans (le_max_left _ _) (le_max_left _ _)
  · exact le_trans (le_max_right _ _) (le_max_left _ _)
  have h1_le_p : 1 ≤ p :=
    le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hp
  have hpmin_le_p : p_min ≤ p :=
    le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hp
  have hm_le_p : m ≤ p := le_trans (le_max_right _ _) hp
  have hm_le_p_real : (m : ℝ) ≤ (p : ℝ) := by exact_mod_cast hm_le_p
  have hp_pos : 0 < (p : ℝ) := by exact_mod_cast h1_le_p
  have hB_lt_p : B < (p : ℝ) := lt_of_lt_of_le hm hm_le_p_real
  -- Hence C / p ≤ (1 − λ_∞) / 2.
  have h_Cp_bound : C / (p : ℝ) ≤ (1 - lam_infty) / 2 := by
    rcases le_or_lt C 0 with hC_le | _
    · have hC_zero : C = 0 := le_antisymm hC_le hC_nn
      rw [hC_zero, zero_div]
      linarith
    · have h1 : 2 * C < (1 - lam_infty) * (p : ℝ) := by
        have := (div_lt_iff h_gap_pos).mp hB_lt_p
        linarith
      rw [div_le_div_iff hp_pos (by linarith : (0:ℝ) < 2)]
      linarith
  have h_abs := h_taylor p hpmin_le_p h1_le_p
  have h_lam_bound : lam p - lam_infty ≤ C / (p : ℝ) :=
    le_trans (le_abs_self _) h_abs
  have h_lam_ub : lam p ≤ lam_infty + C / (p : ℝ) := by linarith
  have : lam p ≤ lam_infty + (1 - lam_infty) / 2 :=
    le_trans h_lam_ub (by linarith)
  linarith

/-- **Corollary: the midpoint `(1 + λ_∞) / 2` is `< 1`.** -/
theorem midpoint_lt_one (lam_infty : ℝ) (h : lam_infty < 1) :
    (1 + lam_infty) / 2 < 1 := by linarith

/-- **Corollary: the midpoint `(1 + λ_∞) / 2` is positive** for
    `λ_∞ > -1`, in particular whenever `0 < λ_∞`. -/
theorem midpoint_pos_of_pos (lam_infty : ℝ) (h : 0 < lam_infty) :
    0 < (1 + lam_infty) / 2 := by linarith

end UniformFromConvergence

/-! ============================================================
  §23.4  ★ Super grand-master in a bounded x-regime
============================================================ -/

section BoundedRegime

/-- **★ Super grand-master — p-unconditional in a bounded x-regime.**

    Fix `x > 1`. Let `(K_p, R_p, e_p)_p` be a family of Pandrosion-
    Steffensen error sequences with:
      • per-p linear contraction at rate `λ_p` near the fixed point:
           `e_p(k+1) ≤ λ_p · e_p(k)`,
      • per-p quadratic contraction (always):
           `e_p(k+1) ≤ K_p · e_p(k)²`,
      • Taylor convergence of contraction rates to the asymptote:
           `|λ_p − λ_∞(x)| ≤ C / p`        for some `C ≥ 0`,
      • uniform scaled initial error:
           `K_p · R_p ≤ M`                   for some `M > 0`.

    Then there exists `p₀` and a p-UNIFORM iteration count `N` such
    that for every `p ≥ p₀`, `K_p · e_p(k) ≤ δ` for all `k ≥ N`.

    **Novel content over §22.** The uniform-linear-rate hypothesis
    `λ_0 < 1` (required by §22.4) is *automatically* derived from
    • the Taylor input `|λ_p − λ_∞(x)| ≤ C / p`,
    • the Mathlib-proved fact `λ_∞(x) < 1` for `x > 1`
      (§23.2 `lamInfty_lt_one_of_one_lt`).
    The effective `λ_0 = (1 + λ_∞(x)) / 2` is strictly less than 1 and
    depends only on `x` (not on `p`), so the whole pipeline §21→§22→§23
    produces a super-grand-master bound **uniform in p** — conditional
    ONLY on `x > 1` and the Taylor input, with no additional `p`-regime
    restriction. -/
theorem super_grand_master_uniform_bounded_regime
    (x C M ρ δ : ℝ)
    (hx : 1 < x)
    (hC_nn : 0 ≤ C)
    (hM_pos : 0 < M)
    (hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (p_min : ℕ)
    (lam : ℕ → ℝ)
    (K : ℕ → ℝ) (R : ℕ → ℝ) (e : ℕ → ℕ → ℝ)
    (hK : ∀ p, 0 < K p) (hR : ∀ p, 0 < R p)
    (h_taylor : ∀ p : ℕ, p_min ≤ p → 1 ≤ p → |lam p - lamInfty x| ≤ C / (p : ℝ))
    (he_nn : ∀ p k, 0 ≤ e p k)
    (he_0 : ∀ p, e p 0 ≤ R p)
    (he_linear : ∀ p k, e p (k + 1) ≤ lam p * e p k)
    (he_quad : ∀ p k, e p (k + 1) ≤ K p * (e p k) ^ 2)
    (h_KR_uniform : ∀ p, K p * R p ≤ M) :
    ∃ (p₀ : ℕ) (N : ℕ), 1 ≤ p₀ ∧
      ∀ p, p₀ ≤ p → ∀ k ≥ N, K p * e p k ≤ δ := by
  -- Step 1: λ_∞(x) < 1 (§23.2).
  have hlam_inf_lt_1 : lamInfty x < 1 := lamInfty_lt_one_of_one_lt x hx
  have hlam_inf_pos : 0 < lamInfty x := lamInfty_pos_of_one_lt x hx
  -- Step 2: the midpoint (1 + λ_∞)/2 is the uniform λ_0.
  set lam_0 : ℝ := (1 + lamInfty x) / 2
  have hlam_0_lt_1 : lam_0 < 1 := midpoint_lt_one (lamInfty x) hlam_inf_lt_1
  have hlam_0_pos : 0 < lam_0 := midpoint_pos_of_pos (lamInfty x) hlam_inf_pos
  -- Step 3: obtain p₀ from §23.3 such that λ_p ≤ lam_0 for all p ≥ p₀.
  obtain ⟨p₀, hp₀_ge_1, _hp₀_ge_pmin, h_lam_uniform⟩ :=
    uniform_contraction_from_rate_convergence (lamInfty x) C hlam_inf_lt_1 hC_nn
      p_min lam h_taylor
  -- Step 4: restrict the family to p ≥ p₀ via a reparametrisation.
  -- Index ι := {p : ℕ // p₀ ≤ p}; then on ι, λ satisfies the uniform
  -- linear bound λ p ≤ lam_0.
  let ι : Type := {p : ℕ // p₀ ≤ p}
  have h_lam_uniform' : ∀ q : ι, lam q.1 ≤ lam_0 :=
    fun q => h_lam_uniform q.1 q.2
  -- Repackage the hypotheses over ι.
  set K' : ι → ℝ := fun q => K q.1
  set R' : ι → ℝ := fun q => R q.1
  set e' : ι → ℕ → ℝ := fun q k => e q.1 k
  have hK' : ∀ q : ι, 0 < K' q := fun q => hK q.1
  have hR' : ∀ q : ι, 0 < R' q := fun q => hR q.1
  have he'_nn : ∀ q k, 0 ≤ e' q k := fun q k => he_nn q.1 k
  have he'_0 : ∀ q, e' q 0 ≤ R' q := fun q => he_0 q.1
  -- The linear recurrence with UNIFORM rate lam_0.
  have he'_linear : ∀ q : ι, ∀ k, e' q (k + 1) ≤ lam_0 * e' q k := by
    intro q k
    calc e' q (k + 1) = e q.1 (k + 1) := rfl
      _ ≤ lam q.1 * e q.1 k := he_linear q.1 k
      _ ≤ lam_0 * e q.1 k := by
          apply mul_le_mul_of_nonneg_right (h_lam_uniform' q) (he_nn q.1 k)
      _ = lam_0 * e' q k := rfl
  have he'_quad : ∀ q : ι, ∀ k, e' q (k + 1) ≤ K' q * (e' q k) ^ 2 :=
    fun q k => he_quad q.1 k
  have h_KR'_uniform : ∀ q : ι, K' q * R' q ≤ M := fun q => h_KR_uniform q.1
  -- Step 5: apply §22.4 on the reparametrised family.
  obtain ⟨N, hN⟩ :=
    super_grand_master_uniform (ι := ι) lam_0 M ρ δ
      hlam_0_pos hlam_0_lt_1 hM_pos hρ_pos hρ_lt_1 hδ_pos
      K' R' e' hK' hR' he'_nn he'_0 he'_linear he'_quad h_KR'_uniform
  refine ⟨p₀, N, hp₀_ge_1, fun p hp k hk => ?_⟩
  -- Extract the bound at (⟨p, hp⟩ : ι).
  have := hN ⟨p, hp⟩ k hk
  exact this

end BoundedRegime

/-! ============================================================
  §23.5  Closed-form contraction rate and Taylor O(1/p) bound
============================================================ -/

section TaylorBound

/-- **Closed-form fixed point** `α_{p,x} = x^{−1/p}` expressed via the
    exponential: `α = exp(−log(x)/p)`. Satisfies `α^p = 1/x`, i.e. `α`
    is the principal real p-th root of `1/x`. -/
noncomputable def alphaP (p : ℕ) (x : ℝ) : ℝ :=
  Real.exp (-Real.log x / (p : ℝ))

/-- `alphaP p x > 0`. Immediate from `Real.exp_pos`. -/
theorem alphaP_pos (p : ℕ) (x : ℝ) : 0 < alphaP p x := Real.exp_pos _

/-- **Closed-form Pandrosion-Steffensen contraction rate.**
    For the Pandrosion map `h(s) = 1 − (x−1)/(x · S_p(s))` at its
    fixed point `α = alphaP p x`, direct differentiation yields
        `h'(α) = [p·(α − 1) + α·(x − 1)] / [α·(x − 1)]`.
    (Algebraic derivation: `S_p(α) = −(x−1)/(x·(α−1))` since
    `α^p = 1/x`, and `S_p'(α) = [p·(α−1) + α·(x−1)] / [x·α·(α−1)²]`
    by quotient-rule + cancellation; substitution gives the formula.) -/
noncomputable def lamModel (p : ℕ) (x : ℝ) : ℝ :=
  ((p : ℝ) * (alphaP p x - 1) + alphaP p x * (x - 1)) /
    (alphaP p x * (x - 1))

/-- **Asymptotic identity.**
    `lamModel p x − lamInfty x = [p·(α−1) + α·log x] / [α·(x−1)]`.

    Proof: decompose `lamModel = 1 + p·(α−1)/(α·(x−1))` and
    `lamInfty = 1 − log(x)/(x−1) = 1 − α·log(x)/(α·(x−1))`. Subtract. -/
theorem lamModel_sub_lamInfty
    (p : ℕ) (x : ℝ) (hx : 1 < x) :
    lamModel p x - lamInfty x =
      ((p : ℝ) * (alphaP p x - 1) + alphaP p x * Real.log x) /
        (alphaP p x * (x - 1)) := by
  unfold lamModel lamInfty
  set α := alphaP p x
  have hα_pos : 0 < α := alphaP_pos p x
  have hα_ne : α ≠ 0 := ne_of_gt hα_pos
  have hx_minus_one_pos : 0 < x - 1 := by linarith
  have hx_minus_one_ne : x - 1 ≠ 0 := ne_of_gt hx_minus_one_pos
  have hαx_ne : α * (x - 1) ≠ 0 := mul_ne_zero hα_ne hx_minus_one_ne
  field_simp
  ring

/-- **Taylor residue lemma.** With `y := −log(x)/p` and `α := exp(y)`,
    define the residue `g(y) := p·(α − 1) + α·log(x)`. Then
        `g(y) = p · ((α − 1 − y) + y · (α − 1))`
    and for `|y| ≤ 1`,
        `|g(y)| ≤ 3 · log(x)² / p`.

    **Proof.** Substituting `log(x) = −p·y`:
      `g(y) = p(α−1) − p·y·α = p·(α − 1 − y·α) = p·((α−1−y) + y·(α−1) − y·y·... )`
    More cleanly: `p(α−1) + α·log(x) = p·(α−1) − p·y·α
      = p·(α − 1 − y) − p·y·(α − 1)`.
    Now `|α − 1 − y| ≤ y²` by Mathlib's `abs_exp_sub_one_sub_id_le`,
    and `|α − 1| ≤ 2|y|` by `abs_exp_sub_one_le`, giving
    `|g| ≤ p·y² + p·|y|·2|y| = 3·p·y² = 3·log(x)²/p`. -/
theorem lamModel_residue_bound
    (p : ℕ) (x : ℝ) (hx : 1 < x)
    (hp_large : Real.log x ≤ (p : ℝ)) :
    |(p : ℝ) * (alphaP p x - 1) + alphaP p x * Real.log x|
      ≤ 3 * (Real.log x) ^ 2 / (p : ℝ) := by
  set α := alphaP p x
  set y : ℝ := -Real.log x / (p : ℝ) with hy_def
  have h_logx_pos : 0 < Real.log x := Real.log_pos hx
  have hp_pos : 0 < (p : ℝ) := lt_of_lt_of_le h_logx_pos hp_large
  have hp_ne : (p : ℝ) ≠ 0 := ne_of_gt hp_pos
  -- |y| = log(x)/p ≤ 1.
  have hy_abs : |y| = Real.log x / (p : ℝ) := by
    rw [hy_def, abs_div, abs_neg, abs_of_pos h_logx_pos,
        abs_of_pos hp_pos]
  have hy_abs_le_one : |y| ≤ 1 := by
    rw [hy_abs]
    exact (div_le_one hp_pos).mpr hp_large
  -- py = -log(x).
  have h_py : (p : ℝ) * y = -Real.log x := by
    rw [hy_def]
    field_simp
    ring
  -- α = exp(y).
  have h_alpha_eq : α = Real.exp y := rfl
  -- Taylor bounds from Mathlib.
  have h_taylor1 : |α - 1 - y| ≤ y ^ 2 := by
    rw [h_alpha_eq]
    exact Real.abs_exp_sub_one_sub_id_le hy_abs_le_one
  have h_taylor2 : |α - 1| ≤ 2 * |y| := by
    rw [h_alpha_eq]
    exact Real.abs_exp_sub_one_le hy_abs_le_one
  -- Rewrite the goal: p(α-1) + α·log(x) = p·(α - 1 - y) - p·y·(α - 1).
  have h_key :
      (p : ℝ) * (α - 1) + α * Real.log x
        = (p : ℝ) * (α - 1 - y) - (p : ℝ) * y * (α - 1) := by
    have h_log_eq : Real.log x = -((p : ℝ) * y) := by linarith
    rw [h_log_eq]; ring
  rw [h_key]
  -- Bound: |p·(α-1-y) - p·y·(α-1)| ≤ p·y² + p·|y|·2|y| = 3·p·y².
  have hbound : |(p : ℝ) * (α - 1 - y) - (p : ℝ) * y * (α - 1)|
      ≤ (p : ℝ) * y ^ 2 + (p : ℝ) * |y| * (2 * |y|) := by
    calc |(p : ℝ) * (α - 1 - y) - (p : ℝ) * y * (α - 1)|
        ≤ |(p : ℝ) * (α - 1 - y)| + |(p : ℝ) * y * (α - 1)| := by
          exact abs_sub _ _
      _ = (p : ℝ) * |α - 1 - y| + (p : ℝ) * |y| * |α - 1| := by
          rw [abs_mul, abs_mul, abs_mul, abs_of_pos hp_pos]
      _ ≤ (p : ℝ) * y ^ 2 + (p : ℝ) * |y| * (2 * |y|) :=
          add_le_add
            (mul_le_mul_of_nonneg_left h_taylor1 hp_pos.le)
            (mul_le_mul_of_nonneg_left h_taylor2
              (mul_nonneg hp_pos.le (abs_nonneg _)))
  -- Simplify: p·y² + p·|y|·2|y| = p·y² + 2·p·y² = 3·p·y².
  have h_abs_sq : |y| * |y| = y ^ 2 := by
    rw [← abs_mul, abs_mul_self]; ring
  have hsum : (p : ℝ) * y ^ 2 + (p : ℝ) * |y| * (2 * |y|)
              = 3 * ((p : ℝ) * y ^ 2) := by
    have : (p : ℝ) * |y| * (2 * |y|) = 2 * ((p : ℝ) * (|y| * |y|)) := by
      ring
    rw [this, h_abs_sq]; ring
  rw [hsum] at hbound
  -- Final: 3·p·y² = 3·log(x)²/p.
  have h_py_sq : (p : ℝ) * y ^ 2 = (Real.log x) ^ 2 / (p : ℝ) := by
    rw [hy_def]
    field_simp
    ring
  rw [h_py_sq] at hbound
  convert hbound using 1
  ring

/-- **Lower bound on `α(x−1)`**: for `x > 1` and `p ≥ 1`,
    `α · (x − 1) ≥ (x − 1) / x`.

    Proof: `α = exp(−log(x)/p) ≥ exp(−log(x)) = 1/x` since
    `−log(x)/p ≥ −log(x)` for `p ≥ 1` and `log x > 0`. -/
theorem alphaP_times_x_minus_one_lb
    (p : ℕ) (x : ℝ) (hx : 1 < x) (hp : 1 ≤ p) :
    (x - 1) / x ≤ alphaP p x * (x - 1) := by
  have h_logx_pos : 0 < Real.log x := Real.log_pos hx
  have hp_pos : 0 < (p : ℝ) := by exact_mod_cast hp
  have hx_pos : 0 < x := by linarith
  have hx_minus_one_pos : 0 < x - 1 := by linarith
  -- Want -log(x) ≤ -log(x)/p, i.e. -log(x)·(1 - 1/p) ≤ 0.
  -- Since log(x) > 0 and 1 - 1/p ≥ 0 for p ≥ 1, we have
  -- -log(x)·(1 - 1/p) ≤ 0 ✓.
  have h_y_ge : -Real.log x ≤ -Real.log x / (p : ℝ) := by
    rw [le_div_iff hp_pos]
    have : -Real.log x * (p : ℝ) ≤ -Real.log x * 1 := by
      apply mul_le_mul_of_nonpos_left _ (by linarith : -Real.log x ≤ 0)
      exact_mod_cast hp
    linarith
  -- α = exp(-log x / p) ≥ exp(-log x) = 1/x.
  have h_alpha_ge : 1 / x ≤ alphaP p x := by
    unfold alphaP
    have h1 : Real.exp (-Real.log x) = 1 / x := by
      rw [Real.exp_neg, Real.exp_log hx_pos, one_div]
    rw [← h1]
    exact Real.exp_le_exp.mpr h_y_ge
  -- Multiply by (x - 1) > 0.
  have hmul := mul_le_mul_of_nonneg_right h_alpha_ge hx_minus_one_pos.le
  have h_eq : (1 : ℝ) / x * (x - 1) = (x - 1) / x := by
    rw [div_mul_eq_mul_div, one_mul]
  rw [h_eq] at hmul
  exact hmul

/-- **★ Taylor bound for the Pandrosion-Steffensen contraction rate.**

    For every `x > 1`, the closed-form contraction rate
    `lamModel p x = h'(α_{p,x})` converges to its asymptote
    `lamInfty x = 1 − log(x)/(x−1)` at rate `O(1/p)`:
        `|lamModel p x − lamInfty x| ≤ C(x) / p`  for every `p ≥ pMin(x)`,
    where
        `C(x) := 3 · x · log(x)² / (x − 1)`,
        `pMin(x) := ⌈log x⌉₊ + 1`  (ensures `|y| = log(x)/p ≤ 1`).

    This is the **Taylor input** of §23.3/§23.4: it discharges the
    hypothesis `|λ_p − λ_∞(x)| ≤ C/p` with an *explicit* constant
    `C(x)` derived from `x` alone. -/
theorem lamModel_taylor_bound
    (x : ℝ) (hx : 1 < x) :
    ∀ p : ℕ, Nat.ceil (Real.log x) + 1 ≤ p → 1 ≤ p →
      |lamModel p x - lamInfty x|
        ≤ (3 * x * (Real.log x) ^ 2 / (x - 1)) / (p : ℝ) := by
  intro p hp_min hp_one
  have hx_pos : 0 < x := by linarith
  have hx_minus_one_pos : 0 < x - 1 := by linarith
  -- p ≥ ⌈log x⌉₊ + 1 > log(x).
  have h_logx_le_p : Real.log x ≤ (p : ℝ) := by
    have h_ceil : Real.log x ≤ (Nat.ceil (Real.log x) : ℝ) := Nat.le_ceil _
    have h_p_ge : ((Nat.ceil (Real.log x) + 1 : ℕ) : ℝ) ≤ (p : ℝ) := by
      exact_mod_cast hp_min
    calc Real.log x ≤ (Nat.ceil (Real.log x) : ℝ) := h_ceil
      _ ≤ ((Nat.ceil (Real.log x) + 1 : ℕ) : ℝ) := by push_cast; linarith
      _ ≤ (p : ℝ) := h_p_ge
  -- α(x-1) ≥ (x-1)/x > 0.
  have h_αx_lb : (x - 1) / x ≤ alphaP p x * (x - 1) :=
    alphaP_times_x_minus_one_lb p x hx hp_one
  have h_αx_pos : 0 < alphaP p x * (x - 1) :=
    mul_pos (alphaP_pos p x) hx_minus_one_pos
  have h_q_pos : 0 < (x - 1) / x := div_pos hx_minus_one_pos hx_pos
  -- Residue bound |p(α-1) + α·log(x)| ≤ 3·log(x)²/p.
  have h_res : |(p : ℝ) * (alphaP p x - 1) + alphaP p x * Real.log x|
               ≤ 3 * (Real.log x) ^ 2 / (p : ℝ) :=
    lamModel_residue_bound p x hx h_logx_le_p
  -- Combine via the identity lamModel - lamInfty = residue / (α(x-1)).
  rw [lamModel_sub_lamInfty p x hx, abs_div]
  have h_abs_den : |alphaP p x * (x - 1)| = alphaP p x * (x - 1) :=
    abs_of_pos h_αx_pos
  rw [h_abs_den]
  -- We need: |residue|/(α(x-1)) ≤ 3·x·log(x)²/((x-1)·p).
  -- Since α(x-1) ≥ (x-1)/x, we have 1/(α(x-1)) ≤ x/(x-1).
  -- Multiplying the residue bound: |res|/(α(x-1)) ≤ 3·log(x)²/p · x/(x-1).
  have h_inv : 1 / (alphaP p x * (x - 1)) ≤ 1 / ((x - 1) / x) := by
    exact one_div_le_one_div_of_le h_q_pos h_αx_lb
  have h_inv_eq : 1 / ((x - 1) / x) = x / (x - 1) := by
    rw [one_div, inv_div]
  rw [h_inv_eq] at h_inv
  have h_res_nn : 0 ≤ |(p : ℝ) * (alphaP p x - 1) + alphaP p x * Real.log x| :=
    abs_nonneg _
  -- |res|/(α(x-1)) = |res| · (1/(α(x-1))) ≤ |res| · x/(x-1).
  have h_step1 : |(p : ℝ) * (alphaP p x - 1) + alphaP p x * Real.log x|
      / (alphaP p x * (x - 1))
      ≤ |(p : ℝ) * (alphaP p x - 1) + alphaP p x * Real.log x|
          * (x / (x - 1)) := by
    rw [div_eq_mul_one_div]
    exact mul_le_mul_of_nonneg_left h_inv h_res_nn
  have h_step2 : |(p : ℝ) * (alphaP p x - 1) + alphaP p x * Real.log x|
      * (x / (x - 1))
      ≤ (3 * (Real.log x) ^ 2 / (p : ℝ)) * (x / (x - 1)) := by
    apply mul_le_mul_of_nonneg_right h_res
    exact le_of_lt (div_pos hx_pos hx_minus_one_pos)
  have h_final_eq :
      (3 * (Real.log x) ^ 2 / (p : ℝ)) * (x / (x - 1))
        = (3 * x * (Real.log x) ^ 2 / (x - 1)) / (p : ℝ) := by
    field_simp
    ring
  linarith [h_step1, h_step2, h_final_eq.le, h_final_eq.ge]

end TaylorBound

/-! ============================================================
  §23.6  ★ Unconditional-in-p super grand-master (bounded x-regime)
============================================================ -/

section UnconditionalBoundedRegime

/-- **★ Super grand-master — fully p-unconditional in a bounded
    x-regime, with explicit constants.**

    Fix `x > 1`. Let `(K_p, R_p, e_p)_p` be a family of error
    sequences satisfying:
      • per-p linear contraction at the *closed-form* rate `lamModel p x`:
           `e_p(k+1) ≤ lamModel p x · e_p(k)`,
      • per-p quadratic contraction:
           `e_p(k+1) ≤ K_p · e_p(k)²`,
      • uniform scaled initial error:
           `K_p · R_p ≤ M`.

    Then there exists `p₀` (depending only on `x, M, ρ`) and a
    p-UNIFORM iteration count `N` such that for all `p ≥ p₀` and
    all `k ≥ N`, `K_p · e_p(k) ≤ δ`.

    **Fully unconditional content over §23.4.** The Taylor bound
    `|lamModel p x − lamInfty x| ≤ C/p` of §23.5
    (`lamModel_taylor_bound`) is *proved*, not assumed. The constant
    `C(x) := 3 · x · log(x)² / (x − 1)` is explicit, and the required
    threshold `p_min(x) := ⌈log x⌉₊ + 1` is explicit.

    Hence this theorem is a **complete p-uniform complexity bound**
    for abstract error sequences satisfying the two-phase Pandrosion-
    Steffensen structure, parameterised only by `x > 1`. -/
theorem super_grand_master_uniform_closed_form
    (x M ρ δ : ℝ)
    (hx : 1 < x)
    (hM_pos : 0 < M)
    (hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (K : ℕ → ℝ) (R : ℕ → ℝ) (e : ℕ → ℕ → ℝ)
    (hK : ∀ p, 0 < K p) (hR : ∀ p, 0 < R p)
    (he_nn : ∀ p k, 0 ≤ e p k)
    (he_0 : ∀ p, e p 0 ≤ R p)
    (he_linear : ∀ p k, e p (k + 1) ≤ lamModel p x * e p k)
    (he_quad : ∀ p k, e p (k + 1) ≤ K p * (e p k) ^ 2)
    (h_KR_uniform : ∀ p, K p * R p ≤ M) :
    ∃ (p₀ : ℕ) (N : ℕ), 1 ≤ p₀ ∧
      ∀ p, p₀ ≤ p → ∀ k ≥ N, K p * e p k ≤ δ := by
  set C : ℝ := 3 * x * (Real.log x) ^ 2 / (x - 1)
  have hC_nn : 0 ≤ C := by
    unfold_let C
    have hx_minus_one_pos : 0 < x - 1 := by linarith
    have : 0 ≤ (Real.log x) ^ 2 := sq_nonneg _
    positivity
  set p_min : ℕ := Nat.ceil (Real.log x) + 1
  have h_taylor : ∀ p : ℕ, p_min ≤ p → 1 ≤ p →
      |lamModel p x - lamInfty x| ≤ C / (p : ℝ) :=
    fun p hp_min hp_one => lamModel_taylor_bound x hx p hp_min hp_one
  -- The unused `_h_lam_pos` hypothesis of §23.4 accepts any dummy.
  exact super_grand_master_uniform_bounded_regime x C M ρ δ hx hC_nn hM_pos
    hρ_pos hρ_lt_1 hδ_pos p_min (lamModel · x) K R e hK hR
    h_taylor he_nn he_0 he_linear he_quad h_KR_uniform

end UnconditionalBoundedRegime

end Pandrosion
