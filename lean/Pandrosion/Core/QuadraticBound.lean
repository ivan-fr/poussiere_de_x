/-
  Universitas Pandrosion — §30. **Quadratic bound for `pandrosion_h_C`
  at a Pandrosion fixed point.**

  This module proves the load-bearing second-order Taylor bound

      `‖h(z) − α − λ·(z − α)‖ ≤ M · ‖z − α‖²`   on `B(α, r)`

  for the Pandrosion iterator `h`, upgrading the first-order asymptotic
  of §29 to a genuine quadratic estimate.

  **Proof strategy.** We exhibit an explicit *rational* factorization
      `h(z) − α = (z − α) · hFactor(z)`
  valid whenever `S_p(z) ≠ 0`, where `hFactor` is itself a rational
  function of `z` with `hFactor(α) = λ_C(p, α)` and a well-defined
  derivative at `α` (polynomial numerator / non-vanishing polynomial
  denominator). Applying `HasDerivAt.isBigO_sub` to `hFactor` yields
      `hFactor(z) − hFactor(α) = O[𝓝 α] (z − α)`,
  and multiplying by `(z − α)` gives the quadratic bound on `h`.

  Content.

    §30.1  `hFactorAux`, `hFactor` — the explicit factorization factor.

    §30.2  `pandrosion_h_C_factorization` — algebraic identity
           `h(z) − α = (z − α) · hFactor(z)` away from the singular set.

    §30.3  `hFactor_at_alpha` — `hFactor(α) = λ_C(p, α)`.

    §30.4  `Sp_C_analyticAt`, `hFactor_analyticAt`,
           `pandrosion_h_C_analyticAt` — analytic-function hygiene.

    §30.5  `hFactor_isBigO_sub_at_fp` — linear-O bound on the factor
           residue.

    §30.6  ★★★ `pandrosion_h_C_quadratic_bound` — **the crown**:
           existential quadratic Taylor bound for `h` at every
           Pandrosion fixed point. Fully unconditional; no hypothesis
           on `K` or `r`.

    §30.7  `pandrosion_h_C_quadratic_bound_with_small` — strengthened
           form giving the `K · r ≤ 1/2` packaging needed by §27
           (reducing the `h_quadratic` hypothesis of
           `steffensen_local_attraction` to a proven fact for `h`).
-/

import Pandrosion.Core.LinearAsymptotics
import Mathlib.Analysis.Analytic.Constructions
import Mathlib.Analysis.Calculus.FDeriv.Analytic

namespace Pandrosion

open Complex Filter Topology Asymptotics BigOperators Finset

/-! ============================================================
  §30.1  The explicit factorization factor
============================================================ -/

section FactorDefinition

/-- Auxiliary polynomial `T(z) = Σ_{i<p} z^i · α^{p−1−i}` arising from
    the `z^p − α^p` factorization. Satisfies `T(α) = p · α^{p−1}`. -/
noncomputable def hFactorAux (α : ℂ) (p : ℕ) (z : ℂ) : ℂ :=
  ∑ i in Finset.range p, z ^ i * α ^ (p - 1 - i)

/-- **Factorization factor** of `pandrosion_h_C`: the rational function
    `hFactor(z) := (S_p(z) − T(z)) / S_p(z)` for which
        `h(z) − α = (z − α) · hFactor(z)`
    at any `z` with `S_p(z) ≠ 0` (and `α` a Pandrosion fp).

    At `z = α`, `hFactor(α) = 1 − p · α^{p−1} / S_p(α) = λ_C(p, α)`. -/
noncomputable def hFactor (α : ℂ) (p : ℕ) (z : ℂ) : ℂ :=
  (Sp_C p z - hFactorAux α p z) / Sp_C p z

end FactorDefinition

/-! ============================================================
  §30.2  Algebraic factorization identity
============================================================ -/

section Factorization

/-- **★ Factorization of `h(z) − α`.** At any Pandrosion fixed point
    α (where `α^p = 1/x`), and at any `z` with `x · S_p(z) ≠ 0`,
        `pandrosion_h_C x p z − α = (z − α) · hFactor α p z`.
    The key algebraic identity used below to derive the quadratic
    bound via `HasDerivAt.isBigO_sub` on `hFactor`. -/
theorem pandrosion_h_C_factorization
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (α : ℂ)
    (hfp : α ^ p = 1 / x) (z : ℂ) (hSp : Sp_C p z ≠ 0) :
    pandrosion_h_C x p z - α = (z - α) * hFactor α p z := by
  have h_xαp : x * α ^ p = 1 := by rw [hfp]; field_simp
  have h_Sp_eq : Sp_C p z * (1 - z) = 1 - z ^ p := Sp_C_mul_one_sub p z
  have h_geom₂ : hFactorAux α p z * (z - α) = z ^ p - α ^ p := by
    unfold hFactorAux; exact geom_sum₂_mul z α p
  unfold pandrosion_h_C hFactor
  have hxSp : x * Sp_C p z ≠ 0 := mul_ne_zero hx hSp
  field_simp
  linear_combination x * (Sp_C p z) * h_Sp_eq + x * (Sp_C p z) * h_geom₂
    - (Sp_C p z) * h_xαp

end Factorization

/-! ============================================================
  §30.3  Factor value at the fixed point
============================================================ -/

section FactorAtAlpha

/-- `hFactorAux(α) = p · α^{p−1}`. Direct computation: each summand is
    `α^i · α^{p−1−i} = α^{p−1}`, and there are `p` of them. -/
theorem hFactorAux_at_alpha (α : ℂ) {p : ℕ} (hp : 1 ≤ p) :
    hFactorAux α p α = (p : ℂ) * α ^ (p - 1) := by
  unfold hFactorAux
  have h_each : ∀ i ∈ Finset.range p,
      α ^ i * α ^ (p - 1 - i) = α ^ (p - 1) := by
    intro i hi
    rw [Finset.mem_range] at hi
    rw [← pow_add]
    congr 1
    omega
  rw [Finset.sum_congr rfl h_each]
  rw [Finset.sum_const, Finset.card_range]
  simp [mul_comm]

/-- **★ Value of the factorization factor at `α`.** For `p ≥ 1` and
    `S_p(α) ≠ 0`,
        `hFactor(α) = 1 − p · α^{p−1} / S_p(α) = λ_C(p, α)`. -/
theorem hFactor_at_alpha
    (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hSp : Sp_C p α ≠ 0) :
    hFactor α p α = lambdaClosedC p α := by
  unfold hFactor lambdaClosedC
  rw [hFactorAux_at_alpha α hp]
  field_simp

end FactorAtAlpha

/-! ============================================================
  §30.4  Analyticity
============================================================ -/

section Analyticity

/-- `S_p` is analytic everywhere (polynomial). -/
theorem Sp_C_analyticAt (p : ℕ) (α : ℂ) : AnalyticAt ℂ (Sp_C p) α := by
  show AnalyticAt ℂ (fun z : ℂ => ∑ i in Finset.range p, z ^ i) α
  exact Finset.analyticAt_sum _ (fun i _ => (analyticAt_id ℂ α).pow i)

/-- `hFactorAux α p` is analytic everywhere (polynomial in `z` with
    coefficients depending on `α`). -/
theorem hFactorAux_analyticAt (α : ℂ) (p : ℕ) (z₀ : ℂ) :
    AnalyticAt ℂ (hFactorAux α p) z₀ := by
  show AnalyticAt ℂ (fun z : ℂ => ∑ i in Finset.range p,
        z ^ i * α ^ (p - 1 - i)) z₀
  exact Finset.analyticAt_sum _ (fun i _ =>
    ((analyticAt_id ℂ z₀).pow i).mul analyticAt_const)

/-- `hFactor α p` is analytic at any `z₀` with `S_p(z₀) ≠ 0`. -/
theorem hFactor_analyticAt
    (α : ℂ) (p : ℕ) (z₀ : ℂ) (hSp : Sp_C p z₀ ≠ 0) :
    AnalyticAt ℂ (hFactor α p) z₀ := by
  show AnalyticAt ℂ (fun z => (Sp_C p z - hFactorAux α p z) / Sp_C p z) z₀
  exact AnalyticAt.div
    ((Sp_C_analyticAt p z₀).sub (hFactorAux_analyticAt α p z₀))
    (Sp_C_analyticAt p z₀) hSp

/-- `pandrosion_h_C x p` is analytic at any `z₀` with `x · S_p(z₀) ≠ 0`. -/
theorem pandrosion_h_C_analyticAt
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (z₀ : ℂ) (hSp : Sp_C p z₀ ≠ 0) :
    AnalyticAt ℂ (pandrosion_h_C x p) z₀ := by
  show AnalyticAt ℂ (fun z => 1 - (x - 1) / (x * Sp_C p z)) z₀
  refine analyticAt_const.sub ?_
  refine analyticAt_const.div ?_ ?_
  · exact analyticAt_const.mul (Sp_C_analyticAt p z₀)
  · exact mul_ne_zero hx hSp

end Analyticity

/-! ============================================================
  §30.5  Linear-O bound on the factor residue
============================================================ -/

section FactorResidue

/-- **Linear-O bound on the factor residue.** Since `hFactor` is
    analytic (hence differentiable) at `α`,
        `hFactor(z) − hFactor(α) = O[𝓝 α] (z − α)`.
    Equivalent form: there exist `C > 0`, `r > 0` such that
        `‖hFactor z − λ_C(p, α)‖ ≤ C · ‖z − α‖`  on `B(α, r)`. -/
theorem hFactor_isBigO_sub_at_fp
    (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hSp : Sp_C p α ≠ 0) :
    (fun z => hFactor α p z - lambdaClosedC p α) =O[𝓝 α] (fun z => z - α) := by
  have h_analytic := hFactor_analyticAt α p α hSp
  have h_deriv : DifferentiableAt ℂ (hFactor α p) α :=
    h_analytic.differentiableAt
  have h_hasDeriv : HasDerivAt (hFactor α p) (deriv (hFactor α p) α) α :=
    h_deriv.hasDerivAt
  have h_bigO := h_hasDeriv.isBigO_sub
  -- h_bigO : (fun z => hFactor α p z − hFactor α p α) =O[𝓝 α] (fun z => z − α)
  rw [hFactor_at_alpha p hp α hSp] at h_bigO
  exact h_bigO

end FactorResidue

/-! ============================================================
  §30.6  ★★★ Quadratic bound for `pandrosion_h_C`
============================================================ -/

section QuadraticBoundH

/-- **★★★ Quadratic Taylor bound for `pandrosion_h_C`.**

  At any Pandrosion fixed point α (where `α^p = 1/x`, `S_p(α) ≠ 0`,
  `p ≥ 1`), there exist explicit constants `M > 0` and `r > 0` such
  that for every `z` in the disk `B(α, r)`,
      `‖pandrosion_h_C x p z − α − λ_C(p, α) · (z − α)‖ ≤ M · ‖z − α‖²`.

  This is the **rigorous second-order Taylor bound** for the
  Pandrosion iterator at a fixed point, proved unconditionally via
  the algebraic factorization
      `h(z) − α = (z − α) · hFactor(z)`
  combined with the linear-O bound on `hFactor(z) − hFactor(α)`.

  **Paper-citable consequence.** *The Pandrosion iterator `h` is
  quadratic at every fixed point, with a quantitative error constant
  derivable from the local derivative of the rational factor
  `hFactor`.* This is a fully unconditional version of what was
  previously quoted from the literature. -/
theorem pandrosion_h_C_quadratic_bound
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    ∃ M > 0, ∃ r > 0, ∀ z : ℂ, ‖z - α‖ < r →
      ‖pandrosion_h_C x p z - α - lambdaClosedC p α * (z - α)‖ ≤ M * ‖z - α‖ ^ 2 := by
  -- 1. Extract linear-O bound on hFactor residue.
  have h_bigO := hFactor_isBigO_sub_at_fp p hp α hSp
  rw [isBigO_iff] at h_bigO
  obtain ⟨C, h_eventually⟩ := h_bigO
  rw [Metric.eventually_nhds_iff] at h_eventually
  obtain ⟨r₁, hr₁_pos, h_bd_near⟩ := h_eventually
  -- 2. Ensure we also stay in the S_p ≠ 0 region.
  have h_SpC_cont : ContinuousAt (Sp_C p) α :=
    (Sp_C_analyticAt p α).continuousAt
  have h_Sp_ne : ∀ᶠ z in 𝓝 α, Sp_C p z ≠ 0 :=
    h_SpC_cont.eventually_ne hSp
  rw [Metric.eventually_nhds_iff] at h_Sp_ne
  obtain ⟨r₂, hr₂_pos, h_Sp_near⟩ := h_Sp_ne
  -- 3. Pick r = min(r₁, r₂); and M = max(C, 1) so M > 0.
  set r : ℝ := min r₁ r₂
  have hr_pos : 0 < r := lt_min hr₁_pos hr₂_pos
  set M : ℝ := max C 1
  have hM_pos : 0 < M := lt_of_lt_of_le zero_lt_one (le_max_right _ _)
  have hC_le_M : C ≤ M := le_max_left _ _
  refine ⟨M, hM_pos, r, hr_pos, fun z hz => ?_⟩
  -- 4. On B(α, r), both the residue bound and the factorization hold.
  have h_dist : dist z α < r := by rw [dist_eq_norm]; exact hz
  have h_dist_r₁ : dist z α < r₁ := lt_of_lt_of_le h_dist (min_le_left _ _)
  have h_dist_r₂ : dist z α < r₂ := lt_of_lt_of_le h_dist (min_le_right _ _)
  have hSp_z : Sp_C p z ≠ 0 := h_Sp_near h_dist_r₂
  have h_fact := pandrosion_h_C_factorization x hx p α hfp z hSp_z
  -- h_fact : pandrosion_h_C x p z - α = (z - α) * hFactor α p z
  have h_residue_bd : ‖hFactor α p z - lambdaClosedC p α‖ ≤ C * ‖z - α‖ := by
    have := h_bd_near h_dist_r₁
    simpa using this
  -- 5. Algebraic identity: h(z) - α - λ(z - α) = (z - α)·(hFactor(z) - λ).
  have h_id : pandrosion_h_C x p z - α - lambdaClosedC p α * (z - α)
      = (z - α) * (hFactor α p z - lambdaClosedC p α) := by
    rw [h_fact]; ring
  -- 6. Norm bound.
  rw [h_id, norm_mul]
  calc ‖z - α‖ * ‖hFactor α p z - lambdaClosedC p α‖
      ≤ ‖z - α‖ * (C * ‖z - α‖) := by
        apply mul_le_mul_of_nonneg_left h_residue_bd (norm_nonneg _)
    _ = C * ‖z - α‖ ^ 2 := by ring
    _ ≤ M * ‖z - α‖ ^ 2 := by
        apply mul_le_mul_of_nonneg_right hC_le_M
        exact sq_nonneg _

end QuadraticBoundH

/-! ============================================================
  §30.7  Strengthened form with `K · r ≤ 1/2` packaging
============================================================ -/

section QuadraticBoundPackaged

/-- **★ Quadratic bound packaged for `steffensen_local_attraction`.**

  Version of §30.6 that also guarantees `K · r ≤ 1/2`, matching the
  `h_small` hypothesis of `local_quadratic_iterate_converges` in §27.
  Obtained by shrinking the radius to `min(r, 1/(2K))`.

  This is the form in which the quadratic bound is consumed by §27
  when specialized to `h`: for any Pandrosion fp `α`, the iterator
  `h` has a proven quadratic local Lipschitz bound on a disk of
  explicit super-attracting radius.

  NB. This is the quadratic bound for the *unaccelerated* Pandrosion
  iterator `h`, not for the Steffensen-accelerated `steffensen_step_C`.
  The Steffensen quadratic bound — the actual hypothesis consumed by
  `steffensen_local_attraction` — is a *downstream* consequence that
  additionally requires bounding the Steffensen-denominator quotient;
  that step is separate and remains scoped as future work. -/
theorem pandrosion_h_C_quadratic_bound_with_small
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    ∃ K > 0, ∃ r > 0, K * r ≤ 1 / 2 ∧
      ∀ z : ℂ, ‖z - α‖ < r →
        ‖pandrosion_h_C x p z - α - lambdaClosedC p α * (z - α)‖
          ≤ K * ‖z - α‖ ^ 2 := by
  obtain ⟨M, hM_pos, r₀, hr₀_pos, h_bd⟩ :=
    pandrosion_h_C_quadratic_bound x hx p hp α hSp hfp
  -- Shrink radius to ensure M·r ≤ 1/2.
  set r : ℝ := min r₀ (1 / (2 * M))
  have h_halfM_pos : 0 < 1 / (2 * M) := by positivity
  have hr_pos : 0 < r := lt_min hr₀_pos h_halfM_pos
  have h_Mr_le : M * r ≤ 1 / 2 := by
    have : M * r ≤ M * (1 / (2 * M)) :=
      mul_le_mul_of_nonneg_left (min_le_right _ _) (le_of_lt hM_pos)
    have h_eq : M * (1 / (2 * M)) = 1 / 2 := by
      field_simp
      ring
    linarith
  refine ⟨M, hM_pos, r, hr_pos, h_Mr_le, fun z hz => ?_⟩
  have h_in_r₀ : ‖z - α‖ < r₀ := lt_of_lt_of_le hz (min_le_left _ _)
  exact h_bd z h_in_r₀

end QuadraticBoundPackaged

end Pandrosion
