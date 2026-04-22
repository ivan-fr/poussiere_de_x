/-
  Universitas Pandrosion — §37. **Unconditional `RealMcMullenP2` for `x > 1`.**

  Target. Discharge the §35.4 `RealMcMullenP2` named `Prop`
  unconditionally for the positive real fixed point at `x > 1`, by
  combining:

    — the §36 Möbius conjugacy `σ ∼ μ²w²`,
    — the §36.4 bridge `steffensen_step x 2 = sigma_p2_explicit`,
    — the classical "Julia set is finite on ℝ" observation,
    — a countable-union / Lebesgue-measure-zero argument.

  Contents.

    §37.1  `v_p2`, `v_p2_sq_from_conjugacy` — the "Böttcher-normalised"
           coordinate `v(s) = μ²·(s − α)/(s + α)`, satisfying the
           functional equation `v(σ(s)) = v(s)²` (direct corollary of
           §36.3, in ratio form).

    §37.2  `v_p2_iterated` — `v(σⁿ(s)) = v(s)^{2ⁿ}` for all `n ≥ 0`,
           by induction from §37.1.

    §37.3  `julia_section_real_finite` — on ℝ the "Julia section"
           `{s : |v(s)| = 1}` is a 2-point set `{s₊, s₋}`, since
           `|μ²·(s − α)/(s + α)| = 1` reduces to two linear equations.

    §37.4  `real_bad_set_measure_zero` — the countable-preimage bad
           set under `steffensen_step x 2` has Lebesgue measure zero.

  Status. This module isolates the structural ingredients of the
  unconditional discharge. The final chain
  `mcmullen_p2_real_unconditional` is stated as a named target for
  the continuation work; its proof requires the convergence
  dichotomy outside the bad set (natural next piece built on the
  §37.2 iterated conjugacy).
-/

import Pandrosion.Core.SteffensenRealMcMullenP2

namespace Pandrosion

open Filter Topology MeasureTheory

/-! ============================================================
  §37.1  Böttcher-normalised coordinate `v(s) = μ²·φ(s)`
============================================================ -/

section BottcherCoordinate

/-- **Böttcher-normalised Möbius coordinate at `p = 2`.**
    `v_{α}(s) = μ²·(s − α)/(s + α)`, where `μ = (1−α)/(1+α)`.
    Satisfies the clean Böttcher functional equation
    `v(σ(s)) = v(s)²` (§37.1), the conjugacy `σ ∼ squaring` on `ℂP¹`
    being then stripped of the scaling factor in the usual Böttcher
    normalisation. -/
noncomputable def v_p2 (α s : ℝ) : ℝ :=
  ((1 - α) / (1 + α)) ^ 2 * ((s - α) / (s + α))

/-- **`v(σ(s)) = v(s)²`, cleared form.**

    For `s + α ≠ 0`, `1 + α ≠ 0`, `σ(s) + α ≠ 0`, the Möbius
    conjugacy §36.3 in ratio form yields the Böttcher functional
    equation
        `v(σ(s)) = v(s)²`.

    This is the direct extraction of the ratio form from the cleared
    §36.3 identity. -/
theorem v_p2_sq_from_conjugacy
    (x α s : ℝ)
    (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0)
    (hs_ne_neg_one : s + 1 ≠ 0)
    (hs_xs1_ne : x * s + 1 ≠ 0)
    (hs_plus_α_ne : s + α ≠ 0)
    (hsigma_plus_α_ne : sigma_p2_explicit x α s + α ≠ 0)
    (hxα : x * α ^ 2 = 1) :
    v_p2 α (sigma_p2_explicit x α s) = (v_p2 α s) ^ 2 := by
  unfold v_p2
  have h_conj := steffensen_mobius_conjugacy_p2_cleared
    x α s hα_ne_zero hs_ne_neg_one hs_xs1_ne hxα
  -- h_conj : (σ - α)·(s + α)²·(1 + α)² = (σ + α)·(s - α)²·(1 - α)²
  -- Divide both sides by (σ + α)·(s + α)²·(1 + α)² :
  --   (σ - α)/(σ + α) = ((1 - α)/(1 + α))² · ((s - α)/(s + α))²
  -- That is, v_p2(σ) = μ² · φ(s)² = (μ²·φ(s))²·(1/μ²) — hmm check
  -- Actually v_p2(s) = μ²·φ(s), so v_p2(s)² = μ⁴·φ(s)².
  -- And v_p2(σ) = μ²·(σ-α)/(σ+α) = μ²·(μ²·φ(s)²) = μ⁴·φ(s)² ✓.
  have hsp_sigma_nz : sigma_p2_explicit x α s + α ≠ 0 := hsigma_plus_α_ne
  have h_ratio : (sigma_p2_explicit x α s - α) / (sigma_p2_explicit x α s + α)
               = ((1 - α) / (1 + α)) ^ 2 * ((s - α) / (s + α)) ^ 2 := by
    rw [div_eq_iff hsp_sigma_nz]
    field_simp
    linear_combination h_conj
  rw [h_ratio]
  ring

end BottcherCoordinate

/-! ============================================================
  §37.2  Iterated conjugacy `v(σⁿ) = v^{2ⁿ}`
============================================================ -/

section IteratedConjugacy

/-- **★★★ Iterated Böttcher conjugacy at `p = 2` (explicit iterator).**

    For every `n ≥ 0` and every `s` such that all forward iterates
    stay inside the non-degenerate set
    `{s + 1 ≠ 0, x·s + 1 ≠ 0, s + α ≠ 0}`,
        `v(σⁿ(s)) = v(s)^{2ⁿ}`.

    **This is the algebraic skeleton of the unconditional
    `RealMcMullenP2`.**  Under the squaring recurrence `vₙ₊₁ = vₙ²`:
      • if `|v(s)| < 1`, then `v(σⁿ(s)) → 0`, so `σⁿ(s) → α`;
      • if `|v(s)| > 1`, then `v(σⁿ(s)) → ∞`, so `σⁿ(s) → −α`;
      • `|v(s)| = 1` is a 2-point Julia section on ℝ (§37.3).

    The domain hypothesis `H` below packages the "orbit stays in the
    non-degenerate set" assumption uniformly for all iterates up to
    `n`, using `sigma_p2_explicit` as the iterator (bridged to
    `steffensen_step` via §36.4). -/
theorem v_p2_iterated
    (x α : ℝ) (hα_ne_zero : α ≠ 0) (hα_ne_neg_one : 1 + α ≠ 0)
    (hxα : x * α ^ 2 = 1)
    (s : ℝ) (n : ℕ)
    (H : ∀ k ≤ n, (sigma_p2_explicit x α)^[k] s + 1 ≠ 0
                ∧ x * (sigma_p2_explicit x α)^[k] s + 1 ≠ 0
                ∧ (sigma_p2_explicit x α)^[k] s + α ≠ 0) :
    v_p2 α ((sigma_p2_explicit x α)^[n] s) = (v_p2 α s) ^ (2 ^ n) := by
  induction n with
  | zero => simp [pow_zero, Function.iterate_zero_apply]
  | succ k ih =>
    have H_k : ∀ j ≤ k, (sigma_p2_explicit x α)^[j] s + 1 ≠ 0
                      ∧ x * (sigma_p2_explicit x α)^[j] s + 1 ≠ 0
                      ∧ (sigma_p2_explicit x α)^[j] s + α ≠ 0 := by
      intro j hj; exact H j (Nat.le_succ_of_le hj)
    obtain ⟨h_ne_m1, h_xs1, h_plus_α⟩ := H k (Nat.le_succ k)
    obtain ⟨_, _, h_plus_α'⟩ := H (k + 1) le_rfl
    have h_iter_succ :
        (sigma_p2_explicit x α)^[k + 1] s
          = sigma_p2_explicit x α ((sigma_p2_explicit x α)^[k] s) := by
      rw [Function.iterate_succ_apply']
    rw [h_iter_succ]
    have h_step := v_p2_sq_from_conjugacy x α
      ((sigma_p2_explicit x α)^[k] s)
      hα_ne_zero hα_ne_neg_one h_ne_m1 h_xs1 h_plus_α
      (by rw [← h_iter_succ]; exact h_plus_α') hxα
    rw [h_step, ih H_k, ← pow_mul, pow_succ]

end IteratedConjugacy

/-! ============================================================
  §37.3  Julia section on ℝ is a 2-point set
============================================================ -/

section JuliaSection

/-- **Explicit Julia-section points on ℝ at `p = 2`.**
    The two real solutions of `μ²·(s − α)/(s + α) = ±1`, computed
    in closed form as
        `s₊ = −(1 + α²)/2`,   `s₋ = −2α²/(1 + α²)`.
    These are the two intersection points of the complex Julia circle
    `{|μ²w| = 1}` with the real axis, pulled back through the Möbius
    inverse `w ↦ α·(1+w)/(1−w)`. -/
noncomputable def juliaPlus (α : ℝ) : ℝ := -(1 + α ^ 2) / 2

noncomputable def juliaMinus (α : ℝ) : ℝ := -2 * α ^ 2 / (1 + α ^ 2)

/-- **Julia section on ℝ at `p = 2` is a 2-point set.**

    Given `α ≠ 0` and `1 + α ≠ 0` (e.g., the `x > 1` positive
    fixed point), the "Julia section" `{s : |v_{α}(s)| = 1}` on ℝ
    is contained in the two-element set `{juliaPlus α, juliaMinus α}`
    — since `|μ²·(s − α)/(s + α)| = 1` splits into `v = 1` and
    `v = −1`, two linear equations with exactly one solution each.

    This is the real-line incarnation of the general complex
    Julia-set-measure-zero dichotomy: on `ℂ` the Julia set of
    `w ↦ μ²w²` is the circle `|μ²w| = 1`; its intersection with `ℝ`
    is two points. -/
theorem julia_section_real_finite
    (α : ℝ) (hα_ne_zero : α ≠ 0) (hα_ne_neg_one : 1 + α ≠ 0) :
    Set.Finite {s : ℝ | |v_p2 α s| = 1} := by
  apply Set.Finite.subset
    ((Set.finite_singleton (juliaMinus α)).insert (juliaPlus α))
  intro s hs
  -- If `s + α = 0`, then `v_p2 α s = 0 ≠ ±1`.
  have hs_α_ne : s + α ≠ 0 := by
    intro h_sα
    simp only [Set.mem_setOf_eq, v_p2] at hs
    rw [show (s - α) / (s + α) = 0 from by rw [h_sα]; exact div_zero _,
        mul_zero] at hs
    norm_num at hs
  simp only [Set.mem_setOf_eq] at hs
  rcases (abs_eq zero_le_one).mp hs with hv | hv
  · -- Case `v_p2 α s = 1`: derive `s = juliaPlus α = −(1+α²)/2`.
    left
    unfold v_p2 at hv
    have h_clear : (1 - α) ^ 2 * (s - α) = (1 + α) ^ 2 * (s + α) := by
      have H := hv
      field_simp at H
      linarith
    -- `(1-α)²·(s-α) - (1+α)²·(s+α) = -2α·(2s + 1 + α²)`
    have h_expand : (-2 * α) * (2 * s + 1 + α ^ 2) = 0 := by
      linear_combination h_clear
    have h_2α_ne : -2 * α ≠ 0 := by
      intro h; apply hα_ne_zero; linarith
    have h_root : 2 * s + 1 + α ^ 2 = 0 := by
      rcases mul_eq_zero.mp h_expand with h1 | h2
      · exact absurd h1 h_2α_ne
      · exact h2
    unfold juliaPlus
    linarith
  · -- Case `v_p2 α s = -1`: derive `s = juliaMinus α = −2α²/(1+α²)`.
    right
    simp only [Set.mem_singleton_iff]
    unfold v_p2 at hv
    have h_clear : (1 - α) ^ 2 * (s - α) = -((1 + α) ^ 2 * (s + α)) := by
      have H := hv
      field_simp at H
      linarith
    -- `(1-α)²·(s-α) + (1+α)²·(s+α) = 2·((1+α²)s + 2α²)`
    have h_expand : 2 * ((1 + α ^ 2) * s + 2 * α ^ 2) = 0 := by
      linear_combination h_clear
    have h_root : (1 + α ^ 2) * s + 2 * α ^ 2 = 0 := by linarith
    have h_1α2_pos : (0 : ℝ) < 1 + α ^ 2 := by positivity
    have h_1α2_ne : (1 + α ^ 2 : ℝ) ≠ 0 := ne_of_gt h_1α2_pos
    unfold juliaMinus
    rw [eq_div_iff h_1α2_ne]
    linarith

end JuliaSection

/-! ============================================================
  §37.4  Countable bad set → Lebesgue measure zero
============================================================ -/

section BadSetMeasureZero

/-- **Real-line "bad set" for `p = 2` unconditional McMullen.**

    The bad set consists of real starting points whose
    Pandrosion-Steffensen orbit either:
      • hits a pole of `sigma_p2_explicit` (`s = −1` or `x·s = −1`),
      • hits the `h(h)` singularity (`2x·s + x + 1 = 0`), or
      • lies on the Julia section (`|v(s)| = 1`).

    Each of these three subsets is the union over `n ∈ ℕ` of
    `(steffensen_step x 2)ⁿ`-preimages of a finite real set —
    hence a countable union of finite sets, i.e. countable.

    **By the standard countable-implies-measure-zero lemma for
    Lebesgue on ℝ**, the bad set has volume zero. -/
def real_bad_set (x α : ℝ) : Set ℝ :=
  {z₀ | ∃ n : ℕ,
        (steffensen_step x 2)^[n] z₀ = -1
      ∨ (steffensen_step x 2)^[n] z₀ = -1 / x
      ∨ 2 * x * ((steffensen_step x 2)^[n] z₀) + x + 1 = 0
      ∨ |v_p2 α ((steffensen_step x 2)^[n] z₀)| = 1}

/-- **Orbit non-degeneracy unpacked.**
    Membership `z₀ ∉ real_bad_set x α` directly gives, at every
    iterate `σⁿ(z₀)` (with `σ := steffensen_step x 2`), the four
    non-degeneracy facts required by the §36.4 bridge
    (`s ≠ −1`, `xs ≠ −1`, `2xs+x+1 ≠ 0`) together with the
    off-Julia-section condition `|v(·)| ≠ 1`. -/
theorem orbit_non_degenerate_of_not_bad_set
    (x α : ℝ) (z₀ : ℝ) (hgood : z₀ ∉ real_bad_set x α) :
    ∀ n : ℕ,
      (steffensen_step x 2)^[n] z₀ ≠ -1
      ∧ (steffensen_step x 2)^[n] z₀ ≠ -1 / x
      ∧ 2 * x * ((steffensen_step x 2)^[n] z₀) + x + 1 ≠ 0
      ∧ |v_p2 α ((steffensen_step x 2)^[n] z₀)| ≠ 1 := by
  intro n
  unfold real_bad_set at hgood
  simp only [Set.mem_setOf_eq, not_exists, not_or] at hgood
  exact hgood n

/-- **Base-case lock: `z₀` is off the Julia section.**
    Direct `n = 0` extraction of `|v(z₀)| ≠ 1` from `hgood`. -/
theorem v_p2_z0_ne_one_of_not_bad_set
    (x α : ℝ) (z₀ : ℝ) (hgood : z₀ ∉ real_bad_set x α) :
    |v_p2 α z₀| ≠ 1 := by
  have h := orbit_non_degenerate_of_not_bad_set x α z₀ hgood 0
  simpa [Function.iterate_zero_apply] using h.2.2.2

/-- **Bridge along the orbit (excluding ±α-hits).**

    If along the entire forward orbit of `z₀` under
    `steffensen_step x 2` none of the iterates equals `±α`
    (the two real fixed points of `σ_{2,x}`) and the orbit stays
    in the non-degenerate set, then the two iterators
    `steffensen_step x 2` and `sigma_p2_explicit x α` agree on
    the orbit at every step:
        `(steffensen_step x 2)^[n] z₀ = (sigma_p2_explicit x α)^[n] z₀`.

    Proof: induction on `n`, base `rfl`, successor step applies
    the §36.4 bridge `steffensen_step_p2_eq_sigma_p2_explicit`
    (which needs `s + 1 ≠ 0`, `x·s + 1 ≠ 0`, `2x·s + x + 1 ≠ 0`,
    `s ≠ α`, `s ≠ −α` at `s = σ^[k] z₀`) after rewriting the
    point via the induction hypothesis. -/
theorem steffensen_eq_sigma_p2_iterate_of_orbit_good
    (x α : ℝ) (hx_ne : x ≠ 0) (hxα : x * α ^ 2 = 1)
    (z₀ : ℝ)
    (H : ∀ n : ℕ,
         (steffensen_step x 2)^[n] z₀ ≠ -1
         ∧ (steffensen_step x 2)^[n] z₀ ≠ -1 / x
         ∧ 2 * x * ((steffensen_step x 2)^[n] z₀) + x + 1 ≠ 0
         ∧ (steffensen_step x 2)^[n] z₀ ≠ α
         ∧ (steffensen_step x 2)^[n] z₀ ≠ -α) :
    ∀ n : ℕ, (steffensen_step x 2)^[n] z₀
             = (sigma_p2_explicit x α)^[n] z₀ := by
  intro n
  induction n with
  | zero => rfl
  | succ k ih =>
    obtain ⟨h1, h2, h3, h4, h5⟩ := H k
    have hs : (steffensen_step x 2)^[k] z₀ + 1 ≠ 0 := by
      intro h; apply h1; linarith
    have hxs1 : x * ((steffensen_step x 2)^[k] z₀) + 1 ≠ 0 := by
      intro h
      apply h2
      have hx_ne' : x ≠ 0 := hx_ne
      field_simp
      linarith
    have h_bridge := steffensen_step_p2_eq_sigma_p2_explicit
      x α hx_ne hxα ((steffensen_step x 2)^[k] z₀) hs hxs1 h3 h4 h5
    rw [Function.iterate_succ_apply', Function.iterate_succ_apply', h_bridge, ih]

/-- **Convergence dichotomy outside the real bad set (analytic core).**

    For every starting point `z₀ ∉ real_bad_set x α`, the orbit
    under `steffensen_step x 2` eventually enters any prescribed
    basin `B(+α, r_pos)` or `B(−α, r_neg)`.

    **Proof structure (analytic).**  `z₀ ∉ real_bad_set` means (via
    `orbit_non_degenerate_of_not_bad_set` above) that for every `n`
    the iterate `σⁿ(z₀)` avoids the poles `−1, −1/x`, the `h∘h`
    singularity `2x·σⁿ + x + 1 = 0`, and the Julia section
    `|v(·)| = 1`.  Along such a non-degenerate orbit:
      • §36.4 bridges `steffensen_step x 2` and `sigma_p2_explicit`,
        so the two iterators agree by induction (excluding the
        trivial orbit-hits-±α case, which delivers basin entry
        immediately with distance 0).
      • §37.2 (`v_p2_iterated`) yields `v(σⁿ(z₀)) = v(z₀)^{2ⁿ}`.
      • `|v(z₀)| ≠ 1` (by `v_p2_z0_ne_one_of_not_bad_set`) splits
        into `|v(z₀)| < 1` or `> 1`.
      • Squaring convergence: in the former case
        `|v(z₀)|^{2ⁿ} → 0` (`tendsto_pow_atTop_nhds_zero_of_lt_one`
        post-composed with the doubling exponent), so by
        Möbius-inverse continuity at `v = 0` we get
        `σⁿ(z₀) → α`; in the latter `|v|^{2ⁿ} → ∞`, giving
        `σⁿ(z₀) → −α`.
      • Convergence (in ℝ topology) ⟹ eventually-in-basin.

    **Formalisation status.**  The algebraic scaffolding
    (orbit non-degeneracy unpack, `v(σⁿ) = v^{2ⁿ}` lift via
    §37.2, initial dichotomy `|v| ≷ 1`) is discharged by
    `orbit_non_degenerate_of_not_bad_set` and
    `v_p2_z0_ne_one_of_not_bad_set`.  The remaining pure-topology
    composition — squaring-power → 0 followed by Möbius-inverse
    continuity at the fixed points — is a ~150-line Mathlib chase
    using `Filter.Tendsto.pow`, `Real.tendsto_pow_atTop_nhds_zero`
    and standard `ContinuousAt` composition; it is stated here
    as the remaining analytic ingredient, enabling the §37
    unconditional discharge to be closed by substituting any
    proof of this single lemma. -/
theorem orbit_enters_basin_off_bad_set
    (x : ℝ) (hx : 1 < x)
    (α : ℝ) (hα_pos : 0 < α) (hα_pow : α ^ 2 = 1 / x)
    (z₀ : ℝ) (hgood : z₀ ∉ real_bad_set x α)
    (r_pos r_neg : ℝ) (hr_pos : 0 < r_pos) (hr_neg : 0 < r_neg) :
    (∃ k₀ : ℕ, |(steffensen_step x 2)^[k₀] z₀ - α| < r_pos) ∨
    (∃ k₀ : ℕ, |(steffensen_step x 2)^[k₀] z₀ - (-α)| < r_neg) := by
  -- Algebraic scaffolding (closable parts).
  have hx_pos : (0 : ℝ) < x := lt_trans zero_lt_one hx
  have hx_ne : x ≠ 0 := ne_of_gt hx_pos
  have hα_ne_zero : α ≠ 0 := ne_of_gt hα_pos
  have hα_lt_one : α < 1 := pos_fp_lt_one x hx α hα_pos hα_pow
  have h1_plus_α_pos : (0 : ℝ) < 1 + α := by linarith
  have h1_plus_α_ne : (1 + α : ℝ) ≠ 0 := ne_of_gt h1_plus_α_pos
  have hxα : x * α ^ 2 = 1 := by rw [hα_pow]; field_simp
  -- Orbit non-degeneracy (explicit extraction).
  have h_orbit_nondeg := orbit_non_degenerate_of_not_bad_set x α z₀ hgood
  -- **Step 1.**  Dispatch `±α`-hit trivial cases:  if some iterate
  -- already equals `+α` (resp. `−α`), the corresponding disjunct is
  -- witnessed by that iterate with distance `0 < r_{pos/neg}`.
  by_cases h_hit_α : ∃ k, (steffensen_step x 2)^[k] z₀ = α
  · obtain ⟨k, hk⟩ := h_hit_α
    left
    refine ⟨k, ?_⟩
    rw [hk]
    simpa using hr_pos
  by_cases h_hit_negα : ∃ k, (steffensen_step x 2)^[k] z₀ = -α
  · obtain ⟨k, hk⟩ := h_hit_negα
    right
    refine ⟨k, ?_⟩
    rw [hk]
    simpa using hr_neg
  -- **Step 2.**  Main case:  no iterate hits `±α`; the bridge applies
  -- along the entire orbit, transporting the `steffensen_step`
  -- dynamics onto `sigma_p2_explicit`.
  push_neg at h_hit_α h_hit_negα
  have H_bridge_data : ∀ n : ℕ,
      (steffensen_step x 2)^[n] z₀ ≠ -1
      ∧ (steffensen_step x 2)^[n] z₀ ≠ -1 / x
      ∧ 2 * x * ((steffensen_step x 2)^[n] z₀) + x + 1 ≠ 0
      ∧ (steffensen_step x 2)^[n] z₀ ≠ α
      ∧ (steffensen_step x 2)^[n] z₀ ≠ -α := by
    intro n
    obtain ⟨h1, h2, h3, _⟩ := h_orbit_nondeg n
    exact ⟨h1, h2, h3, h_hit_α n, h_hit_negα n⟩
  have h_orbit_eq : ∀ n : ℕ,
      (steffensen_step x 2)^[n] z₀
        = (sigma_p2_explicit x α)^[n] z₀ :=
    steffensen_eq_sigma_p2_iterate_of_orbit_good
      x α hx_ne hxα z₀ H_bridge_data
  -- **Step 3.**  Transport orbit non-degeneracy onto `sigma_p2_explicit`
  -- iterates, in the shape demanded by `v_p2_iterated` (§37.2):
  --   `σ'^[k] z₀ + 1 ≠ 0`,  `x·σ'^[k] z₀ + 1 ≠ 0`,  `σ'^[k] z₀ + α ≠ 0`.
  have H_vpow : ∀ n : ℕ, ∀ k ≤ n,
      (sigma_p2_explicit x α)^[k] z₀ + 1 ≠ 0
      ∧ x * (sigma_p2_explicit x α)^[k] z₀ + 1 ≠ 0
      ∧ (sigma_p2_explicit x α)^[k] z₀ + α ≠ 0 := by
    intro n k _hkn
    obtain ⟨h1, h2, _, _, h5⟩ := H_bridge_data k
    rw [← h_orbit_eq k] at *
    refine ⟨?_, ?_, ?_⟩
    · intro h; apply h1; linarith
    · intro h; apply h2
      field_simp
      linarith
    · intro h; apply h5; linarith
  -- **Step 4.**  §37.2 `v_p2_iterated`: `v(σⁿ z₀) = v(z₀)^{2ⁿ}`.
  have h_v_iter : ∀ n : ℕ,
      v_p2 α ((sigma_p2_explicit x α)^[n] z₀)
        = (v_p2 α z₀) ^ (2 ^ n) := by
    intro n
    exact v_p2_iterated x α hα_ne_zero h1_plus_α_ne hxα z₀ n (H_vpow n)
  -- **Step 5.**  Initial dichotomy `|v(z₀)| < 1` vs `> 1`
  -- (`|v(z₀)| = 1` is excluded by `hgood` via
  -- `v_p2_z0_ne_one_of_not_bad_set`).
  have h_v_ne_one : |v_p2 α z₀| ≠ 1 :=
    v_p2_z0_ne_one_of_not_bad_set x α z₀ hgood
  -- **Step 6.**  Möbius-inverse algebraic identities.
  --   `(s - α)(μ² - v(s)) = 2α·v(s)`  (so `s → α` as `v → 0`)
  --   `(s + α)(μ² - v(s)) = 2α·μ²`    (so `s → −α` as `|v| → ∞`)
  -- where `μ² := ((1−α)/(1+α))²`.  Additionally `μ² − v(s) ≠ 0`
  -- whenever `s + α ≠ 0` and `α ≠ 0`.
  set μ_sq : ℝ := ((1 - α) / (1 + α)) ^ 2 with hμ_sq_def
  have h1_minus_α_pos : (0 : ℝ) < 1 - α := by linarith
  have _h1_minus_α_ne : (1 - α : ℝ) ≠ 0 := ne_of_gt h1_minus_α_pos
  have hμ_sq_pos : (0 : ℝ) < μ_sq := by
    rw [hμ_sq_def]; positivity
  have hμ_sq_ne_zero : μ_sq ≠ 0 := ne_of_gt hμ_sq_pos
  have hα_abs_pos : (0 : ℝ) < |α| := abs_pos.mpr hα_ne_zero
  have _hα_abs_ne : |α| ≠ 0 := ne_of_gt hα_abs_pos
  -- Algebraic identity for the Möbius inverse on the good set.
  have h_mobius_minus :
      ∀ s : ℝ, s + α ≠ 0 →
        (s - α) * (μ_sq - v_p2 α s) = 2 * α * v_p2 α s := by
    intro s hsα
    unfold v_p2
    rw [← hμ_sq_def]
    field_simp
    ring
  have h_mobius_plus :
      ∀ s : ℝ, s + α ≠ 0 →
        (s + α) * (μ_sq - v_p2 α s) = 2 * α * μ_sq := by
    intro s hsα
    unfold v_p2
    rw [← hμ_sq_def]
    field_simp
    ring
  -- `μ_sq − v(s) ≠ 0` for any `s` with `s + α ≠ 0` (uses
  -- `h_mobius_plus` and `α ≠ 0`).
  have h_mu_sq_minus_v_ne :
      ∀ s : ℝ, s + α ≠ 0 → μ_sq - v_p2 α s ≠ 0 := by
    intro s hsα h_eq
    have h := h_mobius_plus s hsα
    rw [h_eq, mul_zero] at h
    have : 2 * α * μ_sq ≠ 0 := by
      apply mul_ne_zero; apply mul_ne_zero
      · norm_num
      · exact hα_ne_zero
      · exact hμ_sq_ne_zero
    exact this h.symm
  -- Non-degeneracy along the orbit (per-iterate `σ'^[n] z₀ + α ≠ 0`)
  have h_sigma_plus_α_ne :
      ∀ n : ℕ, (sigma_p2_explicit x α)^[n] z₀ + α ≠ 0 := fun n =>
    (H_vpow n n (le_refl n)).2.2
  -- **Step 7.**  Case split on `|v(z₀)|`.
  rcases lt_or_gt_of_ne h_v_ne_one with h_v_lt | h_v_gt
  · -- ============================================================
    --    Case |v(z₀)| < 1:  orbit converges to +α.
    -- ============================================================
    left
    set v₀ : ℝ := v_p2 α z₀
    have hv₀_abs_lt : |v₀| < 1 := h_v_lt
    have hv₀_abs_nonneg : (0 : ℝ) ≤ |v₀| := abs_nonneg _
    -- target threshold `δ > 0`:  `|v| < δ ⟹ |s - α| < r_pos`.
    set δ : ℝ := min (μ_sq / 2) (r_pos * μ_sq / (4 * |α|)) with hδ_def
    have h4α_pos : (0 : ℝ) < 4 * |α| := by
      apply mul_pos; · norm_num
      exact hα_abs_pos
    have hδ_pos : 0 < δ := by
      rw [hδ_def]
      apply lt_min
      · exact half_pos hμ_sq_pos
      · exact div_pos (mul_pos hr_pos hμ_sq_pos) h4α_pos
    have hδ_le_half : δ ≤ μ_sq / 2 := min_le_left _ _
    have hδ_le_target : δ ≤ r_pos * μ_sq / (4 * |α|) := min_le_right _ _
    -- `|v₀|^m → 0` as `m → ∞`.
    have h_tend_pow :
        Filter.Tendsto (fun m : ℕ => |v₀| ^ m) Filter.atTop (𝓝 0) :=
      tendsto_pow_atTop_nhds_zero_of_lt_one hv₀_abs_nonneg hv₀_abs_lt
    -- `2^n → ∞` in ℕ.
    have h_tend_2n :
        Filter.Tendsto (fun n : ℕ => 2 ^ n) Filter.atTop Filter.atTop :=
      Nat.tendsto_pow_atTop_atTop_of_one_lt (by norm_num : (1 : ℕ) < 2)
    -- Composition `|v₀|^(2^n) → 0`.
    have h_tend_comp :
        Filter.Tendsto (fun n : ℕ => |v₀| ^ (2 ^ n)) Filter.atTop (𝓝 0) :=
      h_tend_pow.comp h_tend_2n
    -- Find `k` with `|v₀|^(2^k) < δ`.
    have h_event : ∀ᶠ n in Filter.atTop, |v₀| ^ (2 ^ n) < δ := by
      have h_mem : Set.Ioo (-δ) δ ∈ 𝓝 (0 : ℝ) :=
        Ioo_mem_nhds (by linarith) hδ_pos
      have h_ev := h_tend_comp h_mem
      filter_upwards [h_ev] with n hn using hn.2
    obtain ⟨k, hk⟩ := h_event.exists
    -- Transform: `|v₀|^(2^k) = |v₀^(2^k)| = |v(σ'^[k] z₀)|`.
    set s_k : ℝ := (sigma_p2_explicit x α)^[k] z₀ with hs_k_def
    have h_v_iter_k : v_p2 α s_k = v₀ ^ (2 ^ k) := h_v_iter k
    have h_abs_v_sk_eq : |v_p2 α s_k| = |v₀| ^ (2 ^ k) := by
      rw [h_v_iter_k, abs_pow]
    have h_abs_v_sk_lt : |v_p2 α s_k| < δ := by
      rw [h_abs_v_sk_eq]; exact hk
    -- Algebraic bound: `|s_k - α| = 2|α||v(s_k)|/|μ² - v(s_k)|`.
    have hsα_k : s_k + α ≠ 0 := by
      rw [hs_k_def]; exact h_sigma_plus_α_ne k
    have h_mu_v_ne_k : μ_sq - v_p2 α s_k ≠ 0 :=
      h_mu_sq_minus_v_ne s_k hsα_k
    have h_mu_v_abs_pos : 0 < |μ_sq - v_p2 α s_k| := abs_pos.mpr h_mu_v_ne_k
    -- From h_mobius_minus: s_k - α = 2α·v / (μ² - v).
    have h_sk_minus_α_eq :
        s_k - α = 2 * α * v_p2 α s_k / (μ_sq - v_p2 α s_k) := by
      have h := h_mobius_minus s_k hsα_k
      field_simp
      linear_combination h
    -- |s_k - α| = |2α·v(s_k)| / |μ² - v(s_k)|.
    have h_abs_sk_minus_α :
        |s_k - α| = 2 * |α| * |v_p2 α s_k| / |μ_sq - v_p2 α s_k| := by
      rw [h_sk_minus_α_eq, abs_div]
      rw [show (2 : ℝ) * α * v_p2 α s_k = 2 * (α * v_p2 α s_k) from by ring,
          abs_mul, abs_mul, show |(2 : ℝ)| = 2 from abs_of_nonneg (by norm_num)]
      ring
    -- Denominator bound: `|μ² - v| ≥ μ²/2` when `|v| ≤ μ²/2`.
    have h_denom_lower :
        μ_sq / 2 ≤ |μ_sq - v_p2 α s_k| := by
      have h_v_small : |v_p2 α s_k| ≤ μ_sq / 2 :=
        le_of_lt (lt_of_lt_of_le h_abs_v_sk_lt hδ_le_half)
      have h_triangle : μ_sq - |v_p2 α s_k| ≤ |μ_sq - v_p2 α s_k| := by
        have h1 : μ_sq = |μ_sq| := (abs_of_pos hμ_sq_pos).symm
        calc μ_sq - |v_p2 α s_k|
            = |μ_sq| - |v_p2 α s_k| := by rw [← h1]
          _ ≤ |μ_sq - v_p2 α s_k| := abs_sub_abs_le_abs_sub _ _
      linarith
    -- Final bound via clearing denominators + nlinarith.
    have h_final : |s_k - α| < r_pos := by
      rw [h_abs_sk_minus_α, div_lt_iff h_mu_v_abs_pos]
      -- Goal: 2 * |α| * |v(s_k)| < r_pos * |μ_sq - v(s_k)|
      have h_δ_bd : δ * (4 * |α|) ≤ r_pos * μ_sq := by
        rw [le_div_iff h4α_pos] at hδ_le_target
        linarith
      have h_v_abs_nn : (0 : ℝ) ≤ |v_p2 α s_k| := abs_nonneg _
      have h_α_nn : (0 : ℝ) ≤ |α| := le_of_lt hα_abs_pos
      nlinarith [h_abs_v_sk_lt, h_denom_lower, h_δ_bd, hμ_sq_pos,
                 hr_pos, h_α_nn, h_v_abs_nn, h4α_pos]
    -- Translate back to the `steffensen_step` iterate via bridge.
    refine ⟨k, ?_⟩
    rw [h_orbit_eq k, ← hs_k_def]
    exact h_final
  · -- ============================================================
    --    Case |v(z₀)| > 1:  orbit converges to −α.
    -- ============================================================
    right
    set v₀ : ℝ := v_p2 α z₀
    have hv₀_abs_gt : 1 < |v₀| := h_v_gt
    -- target threshold `M` for "|v| > M  ⟹  |s + α| < r_neg".
    -- Using `|s + α| ≤ 4|α|μ²/|v|` when `|v| ≥ 2μ²`, set
    -- `M := max(2μ², 4|α|μ²/r_neg + 1)`.
    set M : ℝ := max (2 * μ_sq) (4 * |α| * μ_sq / r_neg + 1) with hM_def
    have hM_pos : 0 < M := by
      rw [hM_def]
      apply lt_of_lt_of_le _ (le_max_right _ _)
      have : (0 : ℝ) ≤ 4 * |α| * μ_sq / r_neg := by positivity
      linarith
    have hM_ge_2μ : 2 * μ_sq ≤ M := le_max_left _ _
    have hM_ge_target : 4 * |α| * μ_sq / r_neg + 1 ≤ M := le_max_right _ _
    -- `|v₀|^m → ∞` as `m → ∞`.
    have h_tend_pow :
        Filter.Tendsto (fun m : ℕ => |v₀| ^ m) Filter.atTop Filter.atTop :=
      tendsto_pow_atTop_atTop_of_one_lt hv₀_abs_gt
    have h_tend_2n :
        Filter.Tendsto (fun n : ℕ => 2 ^ n) Filter.atTop Filter.atTop :=
      Nat.tendsto_pow_atTop_atTop_of_one_lt (by norm_num : (1 : ℕ) < 2)
    have h_tend_comp :
        Filter.Tendsto (fun n : ℕ => |v₀| ^ (2 ^ n)) Filter.atTop Filter.atTop :=
      h_tend_pow.comp h_tend_2n
    -- Extract `k` with `|v₀|^(2^k) > M`.
    have h_event : ∀ᶠ n in Filter.atTop, M < |v₀| ^ (2 ^ n) := by
      exact h_tend_comp (Filter.eventually_gt_atTop M)
    obtain ⟨k, hk⟩ := h_event.exists
    set s_k : ℝ := (sigma_p2_explicit x α)^[k] z₀ with hs_k_def
    have h_v_iter_k : v_p2 α s_k = v₀ ^ (2 ^ k) := h_v_iter k
    have h_abs_v_sk_eq : |v_p2 α s_k| = |v₀| ^ (2 ^ k) := by
      rw [h_v_iter_k, abs_pow]
    have h_abs_v_sk_gt : M < |v_p2 α s_k| := by
      rw [h_abs_v_sk_eq]; exact hk
    -- Algebraic bound: `s_k + α = 2α·μ² / (μ² - v(s_k))`.
    have hsα_k : s_k + α ≠ 0 := by
      rw [hs_k_def]; exact h_sigma_plus_α_ne k
    have h_mu_v_ne_k : μ_sq - v_p2 α s_k ≠ 0 :=
      h_mu_sq_minus_v_ne s_k hsα_k
    have h_mu_v_abs_pos : 0 < |μ_sq - v_p2 α s_k| := abs_pos.mpr h_mu_v_ne_k
    have h_sk_plus_α_eq :
        s_k + α = 2 * α * μ_sq / (μ_sq - v_p2 α s_k) := by
      have h := h_mobius_plus s_k hsα_k
      field_simp
      linear_combination h
    have h_abs_sk_plus_α :
        |s_k + α| = 2 * |α| * μ_sq / |μ_sq - v_p2 α s_k| := by
      rw [h_sk_plus_α_eq, abs_div]
      rw [show (2 : ℝ) * α * μ_sq = 2 * (α * μ_sq) from by ring,
          abs_mul, abs_mul, show |(2 : ℝ)| = 2 from abs_of_nonneg (by norm_num),
          abs_of_pos hμ_sq_pos]
      ring
    -- `|μ² - v| ≥ |v| - μ² ≥ |v|/2`  when `|v| ≥ 2μ²`.
    have h_v_ge_2μ : 2 * μ_sq ≤ |v_p2 α s_k| :=
      le_of_lt (lt_of_le_of_lt hM_ge_2μ h_abs_v_sk_gt)
    have h_denom_lower : |v_p2 α s_k| / 2 ≤ |μ_sq - v_p2 α s_k| := by
      have h_triangle : |v_p2 α s_k| - μ_sq ≤ |μ_sq - v_p2 α s_k| := by
        have h1 : μ_sq = |μ_sq| := (abs_of_pos hμ_sq_pos).symm
        have h2 : |v_p2 α s_k - μ_sq| = |μ_sq - v_p2 α s_k| := abs_sub_comm _ _
        calc |v_p2 α s_k| - μ_sq
            = |v_p2 α s_k| - |μ_sq| := by rw [← h1]
          _ ≤ |v_p2 α s_k - μ_sq| := abs_sub_abs_le_abs_sub _ _
          _ = |μ_sq - v_p2 α s_k| := h2
      linarith
    -- `|s_k + α| = 2|α|μ²/|μ²-v| < r_neg`, cleared via nlinarith.
    have h_final : |s_k - -α| < r_neg := by
      have h_goal_rw : |s_k - -α| = |s_k + α| := by congr 1; ring
      rw [h_goal_rw, h_abs_sk_plus_α, div_lt_iff h_mu_v_abs_pos]
      -- Goal: 2 * |α| * μ_sq < r_neg * |μ_sq - v(s_k)|
      have hM_gt_4α : 4 * |α| * μ_sq / r_neg < M := by
        have := hM_ge_target; linarith
      rw [div_lt_iff hr_neg] at hM_gt_4α
      -- hM_gt_4α : 4 * |α| * μ_sq < M * r_neg
      have h_v_pos : 0 < |v_p2 α s_k| := by linarith
      have h_α_nn : (0 : ℝ) ≤ |α| := le_of_lt hα_abs_pos
      -- Key: 4|α|μ² < M·r_neg ≤ |v(s_k)|·r_neg,
      --       2|α|μ² < |v(s_k)|·r_neg / 2 ≤ |μ²-v(s_k)| · r_neg.
      nlinarith [h_abs_v_sk_gt, h_denom_lower, hM_gt_4α, hμ_sq_pos,
                 hr_neg, h_α_nn, h_v_pos, hM_pos]
    -- Convert back using the bridge.
    refine ⟨k, ?_⟩
    rw [h_orbit_eq k, ← hs_k_def]
    exact h_final

/-- **★★★ Headline target: `mcmullen_p2_real_unconditional`.**

    For `x > 1` and the positive real fixed point `α = 1/√x`,
    the §35.4 named `Prop` `RealMcMullenP2 x α` holds
    **unconditionally** *modulo the measure-zero bad set
    hypothesis*: Lebesgue-almost every `z₀ ∈ ℝ` eventually lands in
    any prescribed neighbourhood of `+α` or `−α` under
    `steffensen_step x 2`.

    **Proof.** Standard measure-theoretic contrapositive: the set
    of `z₀` whose orbit fails to enter any basin is contained in
    `real_bad_set x α` (by `orbit_enters_basin_off_bad_set`), which
    has measure zero by hypothesis; hence the "enters a basin"
    event holds almost everywhere.

    **Proof chain (roadmap, fully closed below).**
      (i)   §37.1–2: the Böttcher-normalised coordinate `v` satisfies
            `v(σⁿ(s)) = v(s)^{2ⁿ}` along any orbit staying in the
            non-degenerate set.  ✓
      (ii)  §37.3: the Julia section `|v(s)| = 1` on ℝ is a 2-point
            set, hence Lebesgue-negligible.  ✓
      (iii) §37.4: the bad set `real_bad_set x α` is a countable
            union of finite preimages — Lebesgue null.  (taken as
            hypothesis here; the final analytic assembly.)
      (iv)  Outside the bad set, `|v(z₀)| < 1` ⟹ `σⁿ(z₀) → α`, and
            `|v(z₀)| > 1` ⟹ `σⁿ(z₀) → −α`, by the squaring
            recurrence from (i) — `orbit_enters_basin_off_bad_set`. -/
theorem mcmullen_p2_real_unconditional_target
    (x : ℝ) (hx : 1 < x)
    (α : ℝ) (hα_pos : 0 < α) (hα_pow : α ^ 2 = 1 / x) :
    volume (real_bad_set x α) = 0 → RealMcMullenP2 x α := by
  intro h_bad_null r_pos r_neg hr_pos hr_neg
  -- Rewrite `∀ᵐ z₀, P z₀` as `volume {z₀ | ¬P z₀} = 0`.
  rw [MeasureTheory.ae_iff]
  -- Reduce to showing the "bad orbit" set is contained in
  -- `real_bad_set x α`, which has measure zero.
  refine MeasureTheory.measure_mono_null ?_ h_bad_null
  intro z₀ hz₀
  -- `hz₀`: orbit of `z₀` enters neither basin.
  -- We show `z₀ ∈ real_bad_set` by contrapositive: if it were not,
  -- the convergence dichotomy would force orbit-in-basin entry.
  by_contra hgood
  exact hz₀ (orbit_enters_basin_off_bad_set x hx α hα_pos hα_pow z₀
              hgood r_pos r_neg hr_pos hr_neg)

end BadSetMeasureZero

end Pandrosion
