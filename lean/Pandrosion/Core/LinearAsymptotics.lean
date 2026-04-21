/-
  Universitas Pandrosion — §29. **Linear asymptotics near a Pandrosion
  fixed point.**

  Structural bridge between §26 `ComplexMultiplier` (the `HasDerivAt`
  statement for the Pandrosion iterator `h`) and the quadratic bound
  needed by §27-§28 (`steffensen_local_attraction`,
  `steffensen_dynamical_convergence_ae`).

  Content.

    §29.1  `lambdaClosedC_ne_one_at_fp` — at every Pandrosion fixed
           point `α ≠ 0` with `p ≥ 1`, the closed-form multiplier
           `λ_C(p, α) = 1 − p·α^{p−1}/S_p(α)` is *not* 1. This is a
           purely algebraic fact and is a **load-bearing ingredient**
           for the Steffensen denominator non-vanishing argument:
           the leading-order behaviour
               `steffensen_denom_C(z) ~ (λ − 1)² · (z − α)`
           requires `λ ≠ 1`.

    §29.2  `cycAnchor_ne_zero` — cyclotomic anchors are non-zero when
           `α ≠ 0`.

    §29.3  `cycAnchor_lambdaClosedC_ne_one` — specialization of §29.1
           to the canonical cyclotomic anchor family.

    §29.4  `pandrosion_h_C_locally_lipschitz_at_fp` — **quantitative
           linear contraction for `h`**: for any `c > |λ_C(p, α)|` at
           a Pandrosion fixed point, there exists a radius `r > 0`
           such that
               `‖pandrosion_h_C x p z − α‖ ≤ c · ‖z − α‖`
           on `B(α, r)`. Direct consequence of §26's `HasDerivAt` via
           the `hasDerivAt_iff_isLittleO` Mathlib characterisation.

    §29.5  `pandrosion_h_C_strict_contraction_when_abs_lt_one` —
           immediate corollary: if `|λ_C(p, α)| < 1`, then `h` is a
           strict contraction in a neighbourhood of `α`. This is the
           linear (pre-Steffensen) attraction content — a genuine
           dynamical statement that *does not* require the full
           quadratic bound.
-/

import Pandrosion.Core.ComplexMultiplier
import Pandrosion.Core.CyclotomicMcMullen
import Mathlib.Analysis.Asymptotics.Asymptotics

namespace Pandrosion

open Complex Filter Topology Asymptotics

/-! ============================================================
  §29.1  Closed-form multiplier is never 1 at a Pandrosion fp
============================================================ -/

section LambdaNeOne

/-- **★ Closed-form multiplier ≠ 1 at a Pandrosion fixed point.**
    For `p ≥ 1`, `α ≠ 0`, and `S_p(α) ≠ 0`, the paper multiplier
        `λ_C(p, α) = 1 − p · α^{p−1} / S_p(α)`
    is not equal to 1.

    Proof: `λ_C = 1` iff `p · α^{p−1} / S_p(α) = 0` iff
    `p · α^{p−1} = 0` (using `S_p(α) ≠ 0`) iff `p = 0` or
    `α^{p−1} = 0`. Both are ruled out by `p ≥ 1` and `α ≠ 0`. -/
theorem lambdaClosedC_ne_one_at_fp
    {p : ℕ} (hp : 1 ≤ p) {α : ℂ} (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) :
    lambdaClosedC p α ≠ 1 := by
  intro heq
  unfold lambdaClosedC at heq
  -- heq : 1 - (p : ℂ) * α ^ (p - 1) / Sp_C p α = 1
  have hratio : (p : ℂ) * α ^ (p - 1) / Sp_C p α = 0 := by
    linear_combination -heq
  have h_num : (p : ℂ) * α ^ (p - 1) = 0 :=
    (div_eq_zero_iff.mp hratio).resolve_right hSp
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hα_pow : α ^ (p - 1) ≠ 0 := pow_ne_zero _ hα
  exact (mul_ne_zero hp_ne hα_pow) h_num

end LambdaNeOne

/-! ============================================================
  §29.2  Cyclotomic anchors are non-zero
============================================================ -/

section CycAnchorNonZero

/-- **Cyclotomic anchors are non-zero** whenever `α ≠ 0`.
    Follows from `cycAnchor α p k = α · exp(…)` and the fact that
    `Complex.exp` is never zero. -/
theorem cycAnchor_ne_zero {α : ℂ} (hα : α ≠ 0) (p : ℕ) (k : Fin p) :
    cycAnchor α p k ≠ 0 := by
  unfold cycAnchor
  exact mul_ne_zero hα (Complex.exp_ne_zero _)

end CycAnchorNonZero

/-! ============================================================
  §29.3  λ_C ≠ 1 at cyclotomic anchors
============================================================ -/

section CycLambdaNeOne

/-- **★ Closed-form multiplier is never 1 at a cyclotomic anchor.**
    Specialization of §29.1 to the canonical cyclotomic anchor
    family `γ_k := cycAnchor α p k` for `α ≠ 0`, `p ≥ 1`, assuming
    `S_p(γ_k) ≠ 0`. -/
theorem cycAnchor_lambdaClosedC_ne_one
    {p : ℕ} (hp : 1 ≤ p) {α : ℂ} (hα : α ≠ 0) (k : Fin p)
    (hSp : Sp_C p (cycAnchor α p k) ≠ 0) :
    lambdaClosedC p (cycAnchor α p k) ≠ 1 :=
  lambdaClosedC_ne_one_at_fp hp (cycAnchor_ne_zero hα p k) hSp

end CycLambdaNeOne

/-! ============================================================
  §29.4  Local Lipschitz bound for `pandrosion_h_C`
============================================================ -/

section LocalLipschitz

/-- **★★ Quantitative linear contraction for `h` at a fixed point.**
    For any `c > |λ_C(p, α)|` at a Pandrosion fixed point α, there
    exists a radius `r > 0` such that for every `z` in the disk
    `B(α, r)`,
        `‖pandrosion_h_C x p z − α‖ ≤ c · ‖z − α‖`.

    This is the rigorous linear-rate content of the iterator's
    derivative, extracted from §26 via
    `hasDerivAt_iff_isLittleO`. It is the **linear-attraction
    content** (pre-Steffensen) and a stepping stone to the
    full quadratic bound needed by §27. -/
theorem pandrosion_h_C_locally_lipschitz_at_fp
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (α : ℂ)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    {c : ℝ} (hc : ‖lambdaClosedC p α‖ < c) :
    ∃ r > 0, ∀ z : ℂ, ‖z - α‖ < r →
      ‖pandrosion_h_C x p z - α‖ ≤ c * ‖z - α‖ := by
  -- §26.5: h has derivative λ_C at the fixed point.
  have h_deriv := pandrosion_h_C_hasDerivAt_at_fp x hx p α hSp hfp
  -- h(α) = α (reverse direction of the fixed-point iff).
  have h_fp : pandrosion_h_C x p α = α :=
    pandrosion_h_C_fixed_point_reverse x hx p α hSp hfp
  -- Translate `HasDerivAt` into `isLittleO` form.
  rw [hasDerivAt_iff_isLittleO] at h_deriv
  rw [h_fp] at h_deriv
  -- h_deriv : (fun z => h z − α − (z − α) • λ) =o[𝓝 α] (fun z => z − α)
  -- Pick ε = c − |λ|.
  set ε : ℝ := c - ‖lambdaClosedC p α‖ with hε_def
  have hε_pos : 0 < ε := by rw [hε_def]; linarith
  -- From isLittleO get an eventually bound.
  have h_eventually := h_deriv.definition hε_pos
  -- h_eventually : ∀ᶠ z in 𝓝 α, ‖h z − α − (z − α) • λ‖ ≤ ε * ‖z − α‖
  rw [Metric.eventually_nhds_iff] at h_eventually
  obtain ⟨r, hr_pos, hr⟩ := h_eventually
  refine ⟨r, hr_pos, fun z hz => ?_⟩
  have h_dist : dist z α < r := by rw [dist_eq_norm]; exact hz
  have h_bd := hr h_dist
  -- h_bd : ‖h z − α − (z − α) • λ‖ ≤ ε · ‖z − α‖
  -- Convert smul to mul (on ℂ, `•` on ℂ-module = `*`).
  simp only [smul_eq_mul] at h_bd
  -- h_bd : ‖h z − α − (z − α) * λ‖ ≤ ε * ‖z − α‖
  have h_mul_norm : ‖(z - α) * lambdaClosedC p α‖
      = ‖lambdaClosedC p α‖ * ‖z - α‖ := by
    rw [norm_mul, mul_comm]
  calc ‖pandrosion_h_C x p z - α‖
      = ‖(pandrosion_h_C x p z - α - (z - α) * lambdaClosedC p α)
          + (z - α) * lambdaClosedC p α‖ := by congr 1; ring
    _ ≤ ‖pandrosion_h_C x p z - α - (z - α) * lambdaClosedC p α‖
        + ‖(z - α) * lambdaClosedC p α‖ := norm_add_le _ _
    _ ≤ ε * ‖z - α‖ + ‖lambdaClosedC p α‖ * ‖z - α‖ := by
        apply add_le_add h_bd
        rw [h_mul_norm]
    _ = c * ‖z - α‖ := by rw [hε_def]; ring

end LocalLipschitz

/-! ============================================================
  §29.5  Strict linear contraction when |λ_C| < 1
============================================================ -/

section StrictContraction

/-- **★ Strict linear contraction when `|λ_C| < 1`.**
    At a Pandrosion fixed point α with `|λ_C(p, α)| < 1`, there
    exists a contraction rate `q < 1` and a disk `B(α, r)` on which
    `h` is a genuine `q`-Lipschitz contraction around α.

    This is the honest pre-Steffensen linear attraction result:
    whenever the closed-form multiplier is strictly less than 1
    in magnitude, the iteration is locally attracting at rate `q`.
    Combined with §27.1 `local_contraction_iterate_converges`, it
    yields local convergence of the **unaccelerated** Pandrosion
    iterator `h` (not Steffensen-accelerated) to its fixed point —
    a strictly weaker but fully unconditional statement. -/
theorem pandrosion_h_C_strict_contraction_when_abs_lt_one
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (α : ℂ)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    (h_abs : ‖lambdaClosedC p α‖ < 1) :
    ∃ q : ℝ, 0 ≤ q ∧ q < 1 ∧ ∃ r > 0, ∀ z : ℂ, ‖z - α‖ < r →
      ‖pandrosion_h_C x p z - α‖ ≤ q * ‖z - α‖ := by
  -- Pick q = (|λ| + 1)/2, strictly between |λ| and 1.
  set q : ℝ := (‖lambdaClosedC p α‖ + 1) / 2 with hq_def
  have h_nn : 0 ≤ ‖lambdaClosedC p α‖ := norm_nonneg _
  have hq_nn : 0 ≤ q := by rw [hq_def]; linarith
  have hq_lt_one : q < 1 := by rw [hq_def]; linarith
  have hq_gt : ‖lambdaClosedC p α‖ < q := by rw [hq_def]; linarith
  obtain ⟨r, hr_pos, hr⟩ :=
    pandrosion_h_C_locally_lipschitz_at_fp x hx p α hSp hfp hq_gt
  exact ⟨q, hq_nn, hq_lt_one, r, hr_pos, hr⟩

end StrictContraction

end Pandrosion
