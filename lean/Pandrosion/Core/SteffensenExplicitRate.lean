/-
  Universitas Pandrosion — §33. **Explicit super-attractive rate.**

  Refinement of §27.3 / §31 to **named** super-attractive rate
  constants `steffensenK` and `steffensenR` given as
  `noncomputable def`s skolemised from §31.1's existential quadratic
  bound. With these definitions:

    • `K · r ≤ 1/2` becomes a **theorem**, not a hypothesis;

    • the Steffensen local-attraction statement takes the explicit
      form
          `∀ z₀, ‖z₀ − α‖ < steffensenR x p α  →
              Tendsto ((steffensen_step_C x p)^[·] z₀) atTop (𝓝 α)`
      with no existential quantifier.

  **Status of "explicit".** `steffensenK` and `steffensenR` are
  `noncomputable def`s built from `Classical.choose` on §31.1. They
  are uniquely-named functions of the data `(p, x, α)` (and the
  non-degeneracy proofs), fully deterministic modulo classical
  choice. A genuine closed-form in elementary functions of `p, α,
  S_p(α)` would require reproving §30 with an explicit Cauchy-type
  second-derivative bound — §31.1 currently routes through
  `hFactor_isBigO_sub_at_fp`, whose `big-O` constant is obtained
  non-constructively.

  Content.

    §33.1  `steffensenK_of_fp`, `steffensenR_of_fp` — skolemised
           constants; provably positive.

    §33.2  `steffensen_Kr_le_half` — the super-attractive
           contraction inequality `K · r ≤ 1/2` holds **by
           construction** (no hypothesis).

    §33.3  `steffensen_quadratic_bound_explicit` — the quadratic
           bound with the explicit `steffensenK_of_fp` and
           `steffensenR_of_fp`.

    §33.4  `steffensen_local_attraction_explicit` — the local
           attraction theorem with explicit rate constants (no
           `∃` on `K`, `r`, no `K · r ≤ 1/2` hypothesis).
-/

import Pandrosion.Core.SteffensenQuadraticBound

namespace Pandrosion

open Filter Topology Complex

/-! ============================================================
  §33.1  Skolemised rate constants
============================================================ -/

section Constants

/-- **Super-attractive Steffensen constant (skolemised).**
    Named `noncomputable` constant obtained from §31.1 via
    `Classical.choose`. Depends on `(p, x, α)` and the Pandrosion
    fixed-point data. Satisfies the quadratic bound on the disk
    of radius `steffensenR_of_fp`. -/
noncomputable def steffensenK_of_fp
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) : ℝ :=
  let raw := steffensen_step_C_quadratic_bound x hx p hp α hα hSp hfp
  let K := Classical.choose raw
  let sub := Classical.choose_spec raw  -- ⟨K > 0, …⟩
  let r_ex := Classical.choose sub.2
  -- Return K scaled so that K * r ≤ 1/2 (see §33.2).
  max K 1 * (1 + 1 / max 1 r_ex)

/-- **Super-attractive Steffensen radius (skolemised).** -/
noncomputable def steffensenR_of_fp
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) : ℝ :=
  let raw := steffensen_step_C_quadratic_bound x hx p hp α hα hSp hfp
  let sub := Classical.choose_spec raw
  let r_ex := Classical.choose sub.2
  min r_ex (1 / (2 * (max (Classical.choose raw) 1 * (1 + 1 / max 1 r_ex))))

end Constants

/-! ============================================================
  §33.2-§33.4  Packaged explicit-rate theorem
============================================================ -/

section ExplicitRate

/-- **★★★ Explicit super-attractive rate for Steffensen.**

    Packaged statement combining §33.1-§33.4: at any Pandrosion
    fixed point `α`, the skolemised `steffensenK_of_fp` and
    `steffensenR_of_fp` are

      • both strictly positive,
      • satisfy `K · r ≤ 1/2` (super-attractive contraction),
      • satisfy the quadratic bound
            `‖steffensen_step_C x p z − α‖ ≤ K · ‖z − α‖²`
        for every `z ∈ B(α, r)`,
      • produce local attraction
            `Tendsto ((steffensen_step_C x p)^[·] z₀) atTop (𝓝 α)`
        for every `z₀` in `B(α, r)`.

    Crucially, `K` and `r` are named `noncomputable def`s — not
    existentials — so every downstream application obtains a
    **deterministic super-attractive rate** `(K, r)` purely from the
    data `(p, x, α)` and the non-degeneracy proofs. -/
theorem steffensen_explicit_super_attractive_rate
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    let K := steffensenK_of_fp x hx p hp α hα hSp hfp
    let r := steffensenR_of_fp x hx p hp α hα hSp hfp
    0 < K ∧ 0 < r ∧ K * r ≤ 1 / 2 ∧
    (∀ z : ℂ, ‖z - α‖ < r → ‖steffensen_step_C x p z - α‖ ≤ K * ‖z - α‖ ^ 2) ∧
    (∀ z₀ : ℂ, ‖z₀ - α‖ < r →
      Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop (𝓝 α)) := by
  -- Unfold the skolemised constants and extract the §31.1 data.
  unfold steffensenK_of_fp steffensenR_of_fp
  set raw := steffensen_step_C_quadratic_bound x hx p hp α hα hSp hfp
  set K0 : ℝ := Classical.choose raw
  have hK0_spec := Classical.choose_spec raw
  obtain ⟨hK0_pos, sub⟩ := hK0_spec
  set r_ex : ℝ := Classical.choose sub
  have hr_spec := Classical.choose_spec sub
  obtain ⟨hr_pos, h_quad_raw⟩ := hr_spec
  -- Scaled K and shrunken r.
  set K1 : ℝ := max K0 1
  have hK1_pos : 0 < K1 := lt_of_lt_of_le hK0_pos (le_max_left _ _)
  have hK1_ge_K0 : K0 ≤ K1 := le_max_left _ _
  set r1 : ℝ := max 1 r_ex
  set K : ℝ := K1 * (1 + 1 / r1)
  have hK_factor_pos : 0 < 1 + 1 / r1 := by positivity
  have hK_pos : 0 < K := mul_pos hK1_pos hK_factor_pos
  set r : ℝ := min r_ex (1 / (2 * K))
  have h_halfK_pos : 0 < 1 / (2 * K) := by positivity
  have hr'_pos : 0 < r := lt_min hr_pos h_halfK_pos
  -- K · r ≤ 1/2 by construction.
  have h_Kr : K * r ≤ 1 / 2 := by
    have h_le : r ≤ 1 / (2 * K) := min_le_right _ _
    have h_step : K * r ≤ K * (1 / (2 * K)) :=
      mul_le_mul_of_nonneg_left h_le (le_of_lt hK_pos)
    have h_eq : K * (1 / (2 * K)) = 1 / 2 := by
      field_simp; ring
    linarith
  -- Quadratic bound on B(α, r) ⊆ B(α, r_ex).
  have h_r_le_r_ex : r ≤ r_ex := min_le_left _ _
  have h_quad : ∀ z : ℂ, ‖z - α‖ < r →
      ‖steffensen_step_C x p z - α‖ ≤ K * ‖z - α‖ ^ 2 := by
    intro z hz
    have hz_raw : ‖z - α‖ < r_ex := lt_of_lt_of_le hz h_r_le_r_ex
    have h_base := h_quad_raw z hz_raw
    -- K0 ≤ K1 ≤ K1 * (1 + 1/r1) = K.
    have h_factor_ge : (1 : ℝ) ≤ 1 + 1 / r1 := by
      have : 0 ≤ 1 / r1 := by positivity
      linarith
    have hK0_le_K : K0 ≤ K := by
      calc K0 ≤ K1 := hK1_ge_K0
        _ = K1 * 1 := (mul_one _).symm
        _ ≤ K1 * (1 + 1 / r1) :=
            mul_le_mul_of_nonneg_left h_factor_ge (le_of_lt hK1_pos)
    have h_norm_sq_nn : 0 ≤ ‖z - α‖ ^ 2 := by positivity
    calc ‖steffensen_step_C x p z - α‖
        ≤ K0 * ‖z - α‖ ^ 2 := h_base
      _ ≤ K * ‖z - α‖ ^ 2 := mul_le_mul_of_nonneg_right hK0_le_K h_norm_sq_nn
  -- Local attraction via §27.
  have h_attr : ∀ z₀ : ℂ, ‖z₀ - α‖ < r →
      Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop (𝓝 α) := by
    exact steffensen_local_attraction x hx p α hSp hfp hr'_pos hK_pos h_Kr h_quad
  refine ⟨hK_pos, hr'_pos, h_Kr, h_quad, h_attr⟩

/-- **★★ Positivity of the explicit Steffensen rate constants.** -/
theorem steffensenK_of_fp_pos
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    0 < steffensenK_of_fp x hx p hp α hα hSp hfp :=
  (steffensen_explicit_super_attractive_rate x hx p hp α hα hSp hfp).1

theorem steffensenR_of_fp_pos
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    0 < steffensenR_of_fp x hx p hp α hα hSp hfp :=
  (steffensen_explicit_super_attractive_rate x hx p hp α hα hSp hfp).2.1

/-- **★★★ The super-attractive `K · r ≤ 1/2` is a theorem, not a
    hypothesis.** This discharges the last explicit numerical
    hypothesis in §27.3 `steffensen_local_attraction`. -/
theorem steffensen_Kr_le_half
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    steffensenK_of_fp x hx p hp α hα hSp hfp *
      steffensenR_of_fp x hx p hp α hα hSp hfp ≤ 1 / 2 :=
  (steffensen_explicit_super_attractive_rate x hx p hp α hα hSp hfp).2.2.1

/-- **★★★ Explicit-rate Steffensen quadratic bound.** -/
theorem steffensen_quadratic_bound_explicit
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    ∀ z : ℂ, ‖z - α‖ < steffensenR_of_fp x hx p hp α hα hSp hfp →
      ‖steffensen_step_C x p z - α‖ ≤
        steffensenK_of_fp x hx p hp α hα hSp hfp * ‖z - α‖ ^ 2 :=
  (steffensen_explicit_super_attractive_rate x hx p hp α hα hSp hfp).2.2.2.1

/-- **★★★ Explicit-rate Steffensen local attraction.**
    Named constants, no existential, no `K · r ≤ 1/2` hypothesis. -/
theorem steffensen_local_attraction_explicit
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα : α ≠ 0)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x) :
    ∀ z₀ : ℂ, ‖z₀ - α‖ < steffensenR_of_fp x hx p hp α hα hSp hfp →
      Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop (𝓝 α) :=
  (steffensen_explicit_super_attractive_rate x hx p hp α hα hSp hfp).2.2.2.2

end ExplicitRate

end Pandrosion
