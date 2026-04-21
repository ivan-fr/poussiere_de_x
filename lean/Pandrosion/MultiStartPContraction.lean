/-
  Universitas Pandrosion — Lean 4 Formalization
  CONTRACTIVE REGIME, ITERATED RATE, AND MULTI-START COVERAGE

  Three interlocking theorems for the p-parametric anchor-step:

    (E) Contractive Regime.
        On a box [m, M] with 0 < m ≤ M, if the anchor is close enough
        to the root in the sense
            |a − r| · (p+1) · M^p / m^p ≤ ρ < 1,
        then every anchor-step is ρ-Lipschitz toward r:
            |F_{a,p+1}(s) − r| ≤ ρ · |s − r|.
        This is the general-p analog of `pandrosion_map_lipschitz_pairwise_p2`.

    (F) Iterated Contraction.
        Under (E), n anchor-steps reduce the error geometrically:
            |F^[n](s) − r| ≤ ρ^n · |s − r|.
        Direct induction on n.

    (G) Multi-start Coverage.
        For k equispaced anchors in [m, M] and any r ∈ [m, M], there
        exists an anchor a_i with |a_i − r| ≤ (M − m) / k. Combined
        with (E), this yields an *unconditional* multi-start
        convergence theorem for every p and every root in the box.
-/
import Pandrosion.MultiStartPSharp

namespace Pandrosion

section MultiStartPContraction

/-! ## §1. (E) Contractive Regime (Per-Step Lipschitz ≤ ρ) -/

/-- **Theorem E: Contraction radius.**
    If `|a − r| · (p+1) · M^p / m^p ≤ ρ` on a compact box `[m, M]`
    with `0 < m, 1 ≤ M` and `r^(p+1) = X`, then every anchor-step
    contracts toward `r` with rate ρ:

        |F_{a,p+1}(s) − r| ≤ ρ · |s − r|.

    In particular, when ρ < 1 this gives a genuine Banach contraction
    neighborhood around every root. -/
theorem anchor_step_p_contraction (p : ℕ) (X a r s m M ρ : ℝ)
    (hm : 0 < m) (hM1 : 1 ≤ M)
    (hma : m ≤ a) (haM : a ≤ M)
    (hmr : m ≤ r) (hrM : r ≤ M)
    (hms : m ≤ s) (hsM : s ≤ M)
    (hX : r ^ (p + 1) = X)
    (hρ : |a - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) ≤ ρ) :
    |anchor_step_p (p + 1) X a s - r| ≤ ρ * |s - r| := by
  have h_bound :
      |anchor_step_p (p + 1) X a s - r|
        ≤ ((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p * (|a - r| * |s - r|) :=
    optimal_anchor_bound p X a r s m M hm hM1 hma haM hmr hrM hms hsM hX
  have h_abs_sr_nn : 0 ≤ |s - r| := abs_nonneg _
  have h_rearrange :
      ((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p * (|a - r| * |s - r|)
        = |a - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) * |s - r| := by
    ring
  rw [h_rearrange] at h_bound
  have h_step : |a - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) * |s - r|
              ≤ ρ * |s - r| :=
    mul_le_mul_of_nonneg_right hρ h_abs_sr_nn
  linarith

/-- **Contraction-radius witness.** Explicit existence of a ρ < 1 when
    the anchor is strictly close enough to the root. -/
theorem anchor_step_p_contraction_witness (p : ℕ) (X a r m M : ℝ)
    (hm : 0 < m) (hM1 : 1 ≤ M)
    (hma : m ≤ a) (hmr : m ≤ r)
    (hclose : |a - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) < 1)
    (haM : a ≤ M) (hrM : r ≤ M) :
    ∃ ρ : ℝ, 0 ≤ ρ ∧ ρ < 1 ∧
      |a - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) ≤ ρ ∧
      (∀ s : ℝ, m ≤ s → s ≤ M → r ^ (p + 1) = X →
         |anchor_step_p (p + 1) X a s - r| ≤ ρ * |s - r|) := by
  refine ⟨|a - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p), ?_, hclose, le_refl _, ?_⟩
  · -- 0 ≤ ρ
    have h_abs : 0 ≤ |a - r| := abs_nonneg _
    have hm_nn : 0 ≤ m := le_of_lt hm
    have hMp_nn : 0 ≤ M ^ p := pow_nonneg (le_trans hm_nn (le_trans hma haM)) p
    have hmp_pos : 0 < m ^ p := pow_pos hm p
    have hpp1_nn : (0 : ℝ) ≤ ((p + 1 : ℕ) : ℝ) := by positivity
    have hnum_nn : 0 ≤ ((p + 1 : ℕ) : ℝ) * M ^ p := mul_nonneg hpp1_nn hMp_nn
    exact mul_nonneg h_abs (div_nonneg hnum_nn (le_of_lt hmp_pos))
  · intro s hms hsM hX
    exact anchor_step_p_contraction p X a r s m M _
      hm hM1 hma haM hmr hrM hms hsM hX (le_refl _)

/-! ## §2. (F) Iterated Contraction -/

/-- **Theorem F: Iterated ρ-contraction.**
    Under the hypotheses of (E), `n` anchor-steps with the same anchor
    reduce the error by a factor of `ρ^n`:

        |F^[n](s) − r| ≤ ρ^n · |s − r|.

    Proof: induction on n, using (E) at each step. -/
theorem anchor_iterate_contraction (p : ℕ) (X a r m M ρ : ℝ)
    (hm : 0 < m) (hM1 : 1 ≤ M)
    (hma : m ≤ a) (haM : a ≤ M)
    (hmr : m ≤ r) (hrM : r ≤ M)
    (hX : r ^ (p + 1) = X)
    (hρ_nn : 0 ≤ ρ)
    (hρ : |a - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) ≤ ρ)
    (h_inv : ∀ s, m ≤ s → s ≤ M → m ≤ anchor_step_p (p + 1) X a s
             ∧ anchor_step_p (p + 1) X a s ≤ M) :
    ∀ (n : ℕ) (s : ℝ), m ≤ s → s ≤ M →
      |anchor_iterate (p + 1) X a n s - r| ≤ ρ ^ n * |s - r| := by
  intro n
  induction n with
  | zero =>
    intro s _ _
    rw [anchor_iterate_zero, pow_zero, one_mul]
  | succ n ih =>
    intro s hms hsM
    rw [anchor_iterate_succ]
    -- First establish that the n-th iterate is in [m, M]
    have h_box_preserved :
        ∀ k, ∀ t, m ≤ t → t ≤ M →
           m ≤ anchor_iterate (p + 1) X a k t
           ∧ anchor_iterate (p + 1) X a k t ≤ M := by
      intro k
      induction k with
      | zero => intro t hmt htM; exact ⟨hmt, htM⟩
      | succ k ih_k =>
        intro t hmt htM
        rw [anchor_iterate_succ]
        obtain ⟨hm_k, hM_k⟩ := ih_k t hmt htM
        exact h_inv (anchor_iterate (p + 1) X a k t) hm_k hM_k
    obtain ⟨hm_n, hM_n⟩ := h_box_preserved n s hms hsM
    -- Apply (E) to the n-th iterate
    have h_step := anchor_step_p_contraction p X a r
      (anchor_iterate (p + 1) X a n s) m M ρ
      hm hM1 hma haM hmr hrM hm_n hM_n hX hρ
    have h_ih := ih s hms hsM
    -- Combine
    calc |anchor_step_p (p + 1) X a (anchor_iterate (p + 1) X a n s) - r|
        ≤ ρ * |anchor_iterate (p + 1) X a n s - r| := h_step
      _ ≤ ρ * (ρ ^ n * |s - r|) := mul_le_mul_of_nonneg_left h_ih hρ_nn
      _ = ρ ^ (n + 1) * |s - r| := by rw [pow_succ]; ring

/-! ## §3. (G) Multi-Start Coverage -/

/-- **Equispaced anchor grid.** `anchor_grid m M k i = m + i · (M−m)/k`. -/
noncomputable def anchor_grid (m M : ℝ) (k i : ℕ) : ℝ :=
  m + (i : ℝ) * ((M - m) / k)

/-- **Grid endpoint 0.** -/
@[simp] theorem anchor_grid_zero (m M : ℝ) (k : ℕ) :
    anchor_grid m M k 0 = m := by
  unfold anchor_grid; simp

/-- **Grid endpoint k.** -/
theorem anchor_grid_last (m M : ℝ) (k : ℕ) (hk : 1 ≤ k) :
    anchor_grid m M k k = M := by
  unfold anchor_grid
  have hk_ne : (k : ℝ) ≠ 0 := by exact_mod_cast Nat.one_le_iff_ne_zero.mp hk
  field_simp
  ring

/-- **Grid is in [m, M].** Each grid anchor lies in the box when m ≤ M. -/
theorem anchor_grid_mem (m M : ℝ) (k i : ℕ) (hmM : m ≤ M) (hik : i ≤ k) :
    m ≤ anchor_grid m M k i ∧ anchor_grid m M k i ≤ M := by
  unfold anchor_grid
  by_cases hk : k = 0
  · subst hk
    interval_cases i
    simp
    exact hmM
  · have hk_pos : 0 < k := Nat.pos_of_ne_zero hk
    have hk_pos_r : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk_pos
    have hk_ne : (k : ℝ) ≠ 0 := ne_of_gt hk_pos_r
    have h_diff_nn : 0 ≤ M - m := by linarith
    have h_div_nn : 0 ≤ (M - m) / k := div_nonneg h_diff_nn (le_of_lt hk_pos_r)
    have hi_nn : 0 ≤ (i : ℝ) := Nat.cast_nonneg i
    constructor
    · linarith [mul_nonneg hi_nn h_div_nn]
    · -- i · (M-m)/k ≤ k · (M-m)/k = M - m
      have hi_le : (i : ℝ) ≤ (k : ℝ) := by exact_mod_cast hik
      have h_upper : (i : ℝ) * ((M - m) / k) ≤ (k : ℝ) * ((M - m) / k) :=
        mul_le_mul_of_nonneg_right hi_le h_div_nn
      have h_collapse : (k : ℝ) * ((M - m) / k) = M - m := by
        field_simp; ring
      linarith

/-- **Theorem G: Coverage (strict box, `m < M`).**
    For `k ≥ 1` equispaced anchors on `[m, M]` and any `r ∈ [m, M]`,
    there exists an index `i ≤ k` with
        |anchor_grid m M k i − r| ≤ (M − m) / k. -/
theorem multistart_coverage_strict (m M r : ℝ) (k : ℕ)
    (hk : 1 ≤ k) (hmM_lt : m < M) (hmr : m ≤ r) (hrM : r ≤ M) :
    ∃ i : ℕ, i ≤ k ∧ |anchor_grid m M k i - r| ≤ (M - m) / k := by
  have hk_pos_r : (0 : ℝ) < (k : ℝ) := by exact_mod_cast (hk : 1 ≤ k)
  have h_diff_pos : 0 < M - m := by linarith
  -- Ratio t = (r - m) / (M - m) ∈ [0, 1]
  set t := (r - m) / (M - m) with ht_def
  have h_t_nn : 0 ≤ t := div_nonneg (by linarith) (le_of_lt h_diff_pos)
  have h_t_le_one : t ≤ 1 := by
    rw [ht_def, div_le_one h_diff_pos]; linarith
  -- Scale by k: t · k ∈ [0, k]
  have h_tk_nn : 0 ≤ t * (k : ℝ) := mul_nonneg h_t_nn (le_of_lt hk_pos_r)
  have h_tk_le_k : t * (k : ℝ) ≤ (k : ℝ) := by
    have := mul_le_mul_of_nonneg_right h_t_le_one (le_of_lt hk_pos_r)
    linarith
  -- Choose i = ⌊t · k⌋
  let i : ℕ := ⌊t * (k : ℝ)⌋₊
  have h_i_le_k : i ≤ k := by
    apply Nat.floor_le_of_le
    exact h_tk_le_k
  have h_floor_le : (i : ℝ) ≤ t * k := Nat.floor_le h_tk_nn
  have h_floor_gt : t * k < (i : ℝ) + 1 := Nat.lt_floor_add_one _
  -- The error: anchor_grid m M k i - r = m + i · h - r = i · h - (r - m) = i · h - t · (M−m)
  -- where h = (M−m)/k, so i · h = i · (M−m)/k. Error = (M−m) · (i/k − t).
  refine ⟨i, h_i_le_k, ?_⟩
  unfold anchor_grid
  have h_diff_eq : m + (i : ℝ) * ((M - m) / k) - r
                 = (i : ℝ) * ((M - m) / k) - (r - m) := by ring
  rw [h_diff_eq]
  have hk_ne : (k : ℝ) ≠ 0 := ne_of_gt hk_pos_r
  have h_r_expr : r - m = t * (M - m) := by
    rw [ht_def]; field_simp
  have h_anchor_expr : (i : ℝ) * ((M - m) / k) = (i : ℝ) / k * (M - m) := by ring
  rw [h_r_expr, h_anchor_expr]
  -- Want: |i/k · (M−m) − t · (M−m)| ≤ (M − m)/k
  have h_factor : (i : ℝ) / k * (M - m) - t * (M - m) = ((i : ℝ) / k - t) * (M - m) := by ring
  rw [h_factor, abs_mul, abs_of_pos h_diff_pos]
  -- Need: |i/k − t| · (M − m) ≤ (M − m) / k = (1/k) · (M − m)
  have h_diff_bound : |((i : ℝ) / k - t)| ≤ 1 / k := by
    rw [abs_sub_le_iff]
    refine ⟨?_, ?_⟩
    · -- i/k - t ≤ 1/k: from h_floor_le (i ≤ t·k) we get i/k ≤ t, so i/k - t ≤ 0 ≤ 1/k
      have h_ik_le_t : (i : ℝ) / k ≤ t := by
        rw [div_le_iff hk_pos_r]; exact h_floor_le
      have h_1k_nn : 0 ≤ (1 : ℝ) / k := by positivity
      linarith
    · -- t - i/k ≤ 1/k: from h_floor_gt (t·k < i+1) we get t ≤ (i+1)/k = i/k + 1/k
      have h_t_le : t ≤ ((i : ℝ) + 1) / k := by
        rw [le_div_iff hk_pos_r]; linarith
      have h_split : ((i : ℝ) + 1) / k = (i : ℝ) / k + 1 / k := by
        field_simp
      linarith
  calc |((i : ℝ) / k - t)| * (M - m)
      ≤ (1 / k) * (M - m) := mul_le_mul_of_nonneg_right h_diff_bound (le_of_lt h_diff_pos)
    _ = (M - m) / k := by ring

/-! ## §5. Unconditional Multi-Start Convergence (Corollary) -/

/-- **Theorem G⋆: Unconditional multi-start convergence.**
    For any root `r ∈ [m, M]` with `r^(p+1) = X` on a box where
    `1 ≤ M`, there exists a grid anchor `a` (with `k` chosen large
    enough that `(M−m)/k · (p+1)·M^p/m^p < 1`) giving a contraction
    toward `r` from every starting point in `[m, M]`.

    Proof: pick `i` from (G), set `ρ = (M−m)/k · (p+1)·M^p/m^p`; then
    (E) gives the contraction. Taking `k` large enough makes ρ < 1. -/
theorem multistart_unconditional_contraction (p : ℕ) (X r m M : ℝ) (k : ℕ)
    (hk : 1 ≤ k) (hm : 0 < m) (hM1 : 1 ≤ M) (hmM_lt : m < M)
    (hmr : m ≤ r) (hrM : r ≤ M) (hX : r ^ (p + 1) = X)
    (h_small : ((M - m) / k) * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) < 1) :
    ∃ (i : ℕ) (ρ : ℝ), i ≤ k ∧ 0 ≤ ρ ∧ ρ < 1 ∧
      ∀ s : ℝ, m ≤ s → s ≤ M →
        |anchor_step_p (p + 1) X (anchor_grid m M k i) s - r| ≤ ρ * |s - r| := by
  obtain ⟨i, hik, h_cov⟩ := multistart_coverage_strict m M r k hk hmM_lt hmr hrM
  have hmM : m ≤ M := le_of_lt hmM_lt
  obtain ⟨h_ai_m, h_ai_M⟩ := anchor_grid_mem m M k i hmM hik
  -- Choose ρ to be the tightened constant
  set C := (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) with hC_def
  have hM_nn : 0 ≤ M := le_of_lt (lt_of_lt_of_le zero_lt_one hM1)
  have hMp_nn : 0 ≤ M ^ p := pow_nonneg hM_nn p
  have hmp_pos : 0 < m ^ p := pow_pos hm p
  have hpp1_nn : (0 : ℝ) ≤ ((p + 1 : ℕ) : ℝ) := by positivity
  have hC_nn : 0 ≤ C := by
    rw [hC_def]
    exact div_nonneg (mul_nonneg hpp1_nn hMp_nn) (le_of_lt hmp_pos)
  have h_step_nn : 0 ≤ (M - m) / k := by
    have hk_pos_r : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
    exact div_nonneg (by linarith) (le_of_lt hk_pos_r)
  -- Bound |a_i − r| · C ≤ ((M−m)/k) · C < 1
  have h_abs_cov : |anchor_grid m M k i - r| ≤ (M - m) / k := h_cov
  have h_prod_bound :
      |anchor_grid m M k i - r| * C ≤ ((M - m) / k) * C :=
    mul_le_mul_of_nonneg_right h_abs_cov hC_nn
  refine ⟨i, ((M - m) / k) * C, hik, ?_, ?_, ?_⟩
  · exact mul_nonneg h_step_nn hC_nn
  · exact lt_of_le_of_lt (le_refl _) h_small
  · intro s hms hsM
    exact anchor_step_p_contraction p X (anchor_grid m M k i) r s m M _
      hm hM1 h_ai_m h_ai_M hmr hrM hms hsM hX
      (le_trans h_prod_bound (le_refl _))

/-! ## §6. (F) Self-Anchored Quadratic Iteration

    The Steffensen/T₃ construction accelerates a linearly-converging
    fixed-anchor sequence into a quadratically-converging one via
    Aitken Δ² reanchoring. For general p, the sharper self-anchored
    variant (Theorem A, `sharp_quadratic_bound`) already delivers
    quadratic convergence *per step*, obviating the need for Aitken
    acceleration: each step uses the previous iterate as its own
    anchor, and Theorem A bounds the error by `C · (prev − r)²`.

    We formalize:
      (F₁) `self_anchor_step`: one self-anchored step.
      (F₂) Quadratic rate per step: |x_{n+1} − r| ≤ C · (x_n − r)².
      (F₃) Two-step amplification: |x_2 − r| ≤ C³ · (x_0 − r)^4. -/

/-- **Self-anchored step.** `x ↦ F_{x, p+2}(x)`. Each step uses the
    current value as its own anchor, yielding quadratic convergence
    via Theorem A. -/
noncomputable def self_anchor_step (p : ℕ) (X s : ℝ) : ℝ :=
  anchor_step_p (p + 2) X s s

/-- **Theorem F₁: One-step quadratic.**
    A single self-anchored step reduces error quadratically:
        |self_anchor_step(s) − r| ≤ C · (s − r)²
    where `C = (p+2) · M^(p+1) / m^(p+1)`. Direct consequence of
    Theorem A. -/
theorem self_anchor_step_quadratic (p : ℕ) (X r s m M : ℝ)
    (hm : 0 < m) (hM1 : 1 ≤ M)
    (hmr : m ≤ r) (hrM : r ≤ M) (hms : m ≤ s) (hsM : s ≤ M)
    (hX : r ^ (p + 2) = X) :
    |self_anchor_step p X s - r|
      ≤ ((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1) * (s - r) ^ 2 :=
  sharp_quadratic_bound p X r s m M hm hM1 hmr hrM hms hsM hX

/-- **Two-step self-anchored iterate.** `s ↦ F²(s)`. -/
noncomputable def self_anchor_step_two (p : ℕ) (X s : ℝ) : ℝ :=
  self_anchor_step p X (self_anchor_step p X s)

/-- **Theorem F₂: Two-step amplification.**
    If the first iterate remains in `[m, M]`, two self-anchored
    steps yield a quartic error bound:
        |F²(s) − r| ≤ C³ · (s − r)⁴,
    where `C = (p+2)·M^(p+1)/m^(p+1)`.

    Proof: apply F₁ twice. First step: |F(s) − r| ≤ C · (s−r)².
    Second step (using F(s) as the new "s"): |F²(s) − r| ≤ C · (F(s) − r)²
    ≤ C · (C · (s−r)²)² = C³ · (s−r)⁴. -/
theorem self_anchor_two_step_quartic (p : ℕ) (X r s m M : ℝ)
    (hm : 0 < m) (hM1 : 1 ≤ M)
    (hmr : m ≤ r) (hrM : r ≤ M) (hms : m ≤ s) (hsM : s ≤ M)
    (hX : r ^ (p + 2) = X)
    (h_inv : m ≤ self_anchor_step p X s ∧ self_anchor_step p X s ≤ M) :
    |self_anchor_step_two p X s - r|
      ≤ (((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1)) ^ 3 * (s - r) ^ 4 := by
  set C := ((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1) with hC_def
  have hM_nn : 0 ≤ M := le_of_lt (lt_of_lt_of_le zero_lt_one hM1)
  have hMp1_nn : 0 ≤ M ^ (p + 1) := pow_nonneg hM_nn (p + 1)
  have hmp1_pos : 0 < m ^ (p + 1) := pow_pos hm (p + 1)
  have hpp2_nn : (0 : ℝ) ≤ ((p + 2 : ℕ) : ℝ) := by positivity
  have hC_nn : 0 ≤ C := by
    rw [hC_def]
    exact div_nonneg (mul_nonneg hpp2_nn hMp1_nn) (le_of_lt hmp1_pos)
  -- First step
  have h1 : |self_anchor_step p X s - r| ≤ C * (s - r) ^ 2 :=
    self_anchor_step_quadratic p X r s m M hm hM1 hmr hrM hms hsM hX
  -- Second step (using F(s) as the starting point)
  obtain ⟨h_fs_m, h_fs_M⟩ := h_inv
  have h2 : |self_anchor_step p X (self_anchor_step p X s) - r|
          ≤ C * (self_anchor_step p X s - r) ^ 2 :=
    self_anchor_step_quadratic p X r (self_anchor_step p X s) m M
      hm hM1 hmr hrM h_fs_m h_fs_M hX
  -- Chain: |F(s) − r|² ≤ (C · (s−r)²)² = C² · (s−r)⁴
  have h1_nn : 0 ≤ |self_anchor_step p X s - r| := abs_nonneg _
  have h_sq_sub_r_eq : (self_anchor_step p X s - r) ^ 2
                    = |self_anchor_step p X s - r| ^ 2 := by
    rw [sq_abs]
  have h1_sq : |self_anchor_step p X s - r| ^ 2 ≤ (C * (s - r) ^ 2) ^ 2 :=
    pow_le_pow_left h1_nn h1 2
  have h_bound_sq : (C * (s - r) ^ 2) ^ 2 = C ^ 2 * (s - r) ^ 4 := by ring
  -- Combine: |F²(s) − r| ≤ C · |F(s) − r|² ≤ C · (C · (s−r)²)² = C³ · (s−r)⁴
  calc |self_anchor_step_two p X s - r|
      = |self_anchor_step p X (self_anchor_step p X s) - r| := rfl
    _ ≤ C * (self_anchor_step p X s - r) ^ 2 := h2
    _ = C * |self_anchor_step p X s - r| ^ 2 := by rw [h_sq_sub_r_eq]
    _ ≤ C * (C * (s - r) ^ 2) ^ 2 :=
        mul_le_mul_of_nonneg_left h1_sq hC_nn
    _ = C * (C ^ 2 * (s - r) ^ 4) := by rw [h_bound_sq]
    _ = C ^ 3 * (s - r) ^ 4 := by ring

/-! ## §7. Grand Certificate for (E), (F), (G) -/

/-- **Grand certificate.** Packages (E), (F), (G), (G⋆) into one
    bundled statement. -/
theorem multistartP_contraction_certificate :
    -- (E) Per-step contraction
    (∀ (p : ℕ) (X a r s m M ρ : ℝ),
       0 < m → 1 ≤ M → m ≤ a → a ≤ M → m ≤ r → r ≤ M → m ≤ s → s ≤ M →
       r ^ (p + 1) = X →
       |a - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) ≤ ρ →
       |anchor_step_p (p + 1) X a s - r| ≤ ρ * |s - r|) ∧
    -- (F) Self-anchored quadratic per step
    (∀ (p : ℕ) (X r s m M : ℝ),
       0 < m → 1 ≤ M → m ≤ r → r ≤ M → m ≤ s → s ≤ M →
       r ^ (p + 2) = X →
       |self_anchor_step p X s - r|
         ≤ ((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1) * (s - r) ^ 2) ∧
    -- (G) Coverage
    (∀ (m M r : ℝ) (k : ℕ), 1 ≤ k → m < M → m ≤ r → r ≤ M →
       ∃ i : ℕ, i ≤ k ∧ |anchor_grid m M k i - r| ≤ (M - m) / k) ∧
    -- (G⋆) Unconditional contraction via multi-start
    (∀ (p : ℕ) (X r m M : ℝ) (k : ℕ),
       1 ≤ k → 0 < m → 1 ≤ M → m < M → m ≤ r → r ≤ M →
       r ^ (p + 1) = X →
       ((M - m) / k) * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) < 1 →
       ∃ (i : ℕ) (ρ : ℝ), i ≤ k ∧ 0 ≤ ρ ∧ ρ < 1 ∧
         ∀ s : ℝ, m ≤ s → s ≤ M →
           |anchor_step_p (p + 1) X (anchor_grid m M k i) s - r|
             ≤ ρ * |s - r|) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro p X a r s m M ρ hm hM1 hma haM hmr hrM hms hsM hX hρ
    exact anchor_step_p_contraction p X a r s m M ρ hm hM1 hma haM hmr hrM hms hsM hX hρ
  · intro p X r s m M hm hM1 hmr hrM hms hsM hX
    exact self_anchor_step_quadratic p X r s m M hm hM1 hmr hrM hms hsM hX
  · intro m M r k hk hmM_lt hmr hrM
    exact multistart_coverage_strict m M r k hk hmM_lt hmr hrM
  · intro p X r m M k hk hm hM1 hmM_lt hmr hrM hX h_small
    exact multistart_unconditional_contraction p X r m M k hk hm hM1 hmM_lt hmr hrM hX h_small

end MultiStartPContraction

end Pandrosion
