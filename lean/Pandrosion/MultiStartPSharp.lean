/-
  Universitas Pandrosion — Lean 4 Formalization
  SHARP QUADRATIC BOUNDS FOR THE p-TH ROOT MULTISTART

  Two explicit, constant-tracking convergence theorems for the
  `p`-parametric anchor-step of `MultiStartP.lean`:

    (A) Sharp Quadratic Bound (self-anchor).
        For r, s ∈ [m, M], 0 < m, 1 ≤ M, r^(p+2) = X:
            |F_{s,p+2}(s) - r| ≤ (p+2) · M^(p+1) / m^(p+1) · (s - r)².

    (D) Optimal Anchor Bound.
        For a, r, s ∈ [m, M], 0 < m, 1 ≤ M, r^(p+1) = X:
            |F_{a,p+1}(s) - r| ≤ (p+1) · M^p / m^p · |a - r| · |s - r|.

  Both constants are explicit, polynomial in (p, M/m), and optimal
  in their shape. Proofs combine the master identity
      F_{a,p}(s) - r = (a - r) · (s - r) · R_p(a,r,s) / Q_p(a,s)
  of MultiStartP with explicit polynomial bounds on Q_p and R_p
  over the compact box [m, M].
-/
import Pandrosion.MultiStartP

namespace Pandrosion

section MultiStartPSharp

/-! ## §1. Non-negativity of Q_pow and R_pow -/

theorem Q_pow_nonneg (p : ℕ) (a s : ℝ) (ha : 0 ≤ a) (hs : 0 ≤ s) :
    0 ≤ Q_pow p a s := by
  induction p with
  | zero => simp
  | succ p ih =>
    rw [Q_pow_succ]
    have h1 : 0 ≤ s ^ p := pow_nonneg hs p
    have h2 : 0 ≤ a * Q_pow p a s := mul_nonneg ha ih
    linarith

theorem R_pow_nonneg (p : ℕ) (a r s : ℝ)
    (ha : 0 ≤ a) (hr : 0 ≤ r) (hs : 0 ≤ s) :
    0 ≤ R_pow p a r s := by
  induction p with
  | zero => simp
  | succ p ih =>
    rw [R_pow_succ]
    have h1 : 0 ≤ Q_pow p r s := Q_pow_nonneg p r s hr hs
    have h2 : 0 ≤ a * R_pow p a r s := mul_nonneg ha ih
    linarith

/-! ## §2. Q_pow Bounds over a Compact Box -/

/-- **Upper bound.** `Q_pow (p+1) a s ≤ (p+1) · M^p` for `0 ≤ a, s ≤ M`. -/
theorem Q_pow_upper (p : ℕ) (a s M : ℝ)
    (ha : 0 ≤ a) (haM : a ≤ M) (hs : 0 ≤ s) (hsM : s ≤ M) :
    Q_pow (p + 1) a s ≤ ((p + 1 : ℕ) : ℝ) * M ^ p := by
  have hM_nn : 0 ≤ M := le_trans hs hsM
  induction p with
  | zero =>
    simp [Q_pow_one, pow_zero]
  | succ p ih =>
    rw [Q_pow_succ]
    have hQ_nn : 0 ≤ Q_pow (p + 1) a s := Q_pow_nonneg _ a s ha hs
    have hs_pow : s ^ (p + 1) ≤ M ^ (p + 1) :=
      pow_le_pow_left hs hsM (p + 1)
    have h_aQ_raw : a * Q_pow (p + 1) a s ≤ M * (((p + 1 : ℕ) : ℝ) * M ^ p) :=
      mul_le_mul haM ih hQ_nn hM_nn
    have h_aQ_simp : M * (((p + 1 : ℕ) : ℝ) * M ^ p)
                   = ((p + 1 : ℕ) : ℝ) * M ^ (p + 1) := by rw [pow_succ]; ring
    have h_aQ : a * Q_pow (p + 1) a s ≤ ((p + 1 : ℕ) : ℝ) * M ^ (p + 1) :=
      by rw [← h_aQ_simp]; exact h_aQ_raw
    calc s ^ (p + 1) + a * Q_pow (p + 1) a s
        ≤ M ^ (p + 1) + ((p + 1 : ℕ) : ℝ) * M ^ (p + 1) := add_le_add hs_pow h_aQ
      _ = ((p + 1 + 1 : ℕ) : ℝ) * M ^ (p + 1) := by push_cast; ring

/-- **Lower bound.** `(p+1) · m^p ≤ Q_pow (p+1) a s` when `0 ≤ m ≤ a, s`. -/
theorem Q_pow_lower (p : ℕ) (a s m : ℝ)
    (hm : 0 ≤ m) (hma : m ≤ a) (hms : m ≤ s) :
    ((p + 1 : ℕ) : ℝ) * m ^ p ≤ Q_pow (p + 1) a s := by
  have ha : 0 ≤ a := le_trans hm hma
  induction p with
  | zero =>
    simp [Q_pow_one, pow_zero]
  | succ p ih =>
    rw [Q_pow_succ]
    have hm_pow_nn : 0 ≤ m ^ p := pow_nonneg hm p
    have hs_pow_ge : m ^ (p + 1) ≤ s ^ (p + 1) :=
      pow_le_pow_left hm hms (p + 1)
    have hp1_nn : (0 : ℝ) ≤ ((p + 1 : ℕ) : ℝ) := by positivity
    have h_prod_nn : (0 : ℝ) ≤ ((p + 1 : ℕ) : ℝ) * m ^ p :=
      mul_nonneg hp1_nn hm_pow_nn
    have h_aQ_raw : m * (((p + 1 : ℕ) : ℝ) * m ^ p) ≤ a * Q_pow (p + 1) a s :=
      mul_le_mul hma ih h_prod_nn ha
    have h_aQ_simp : m * (((p + 1 : ℕ) : ℝ) * m ^ p)
                   = ((p + 1 : ℕ) : ℝ) * m ^ (p + 1) := by rw [pow_succ]; ring
    have h_aQ : ((p + 1 : ℕ) : ℝ) * m ^ (p + 1) ≤ a * Q_pow (p + 1) a s := by
      rw [← h_aQ_simp]; exact h_aQ_raw
    calc ((p + 1 + 1 : ℕ) : ℝ) * m ^ (p + 1)
        = m ^ (p + 1) + ((p + 1 : ℕ) : ℝ) * m ^ (p + 1) := by push_cast; ring
      _ ≤ s ^ (p + 1) + a * Q_pow (p + 1) a s := add_le_add hs_pow_ge h_aQ

/-! ## §3. R_pow Upper Bound (requires M ≥ 1) -/

/-- **Upper bound on R_pow.** For `0 ≤ a, r, s ≤ M` and `1 ≤ M`:
    `R_pow (p+1) a r s ≤ (p+1)² · M^p`. -/
theorem R_pow_upper (p : ℕ) (a r s M : ℝ)
    (ha : 0 ≤ a) (haM : a ≤ M) (hr : 0 ≤ r) (hrM : r ≤ M)
    (hs : 0 ≤ s) (hsM : s ≤ M) (hM1 : 1 ≤ M) :
    R_pow (p + 1) a r s ≤ ((p + 1 : ℕ) : ℝ) ^ 2 * M ^ p := by
  have hM_nn : 0 ≤ M := le_trans hs hsM
  induction p with
  | zero =>
    simp [R_pow_one, pow_zero]
  | succ p ih =>
    rw [R_pow_succ]
    have hR_nn : 0 ≤ R_pow (p + 1) a r s := R_pow_nonneg _ a r s ha hr hs
    have hMp_nn : 0 ≤ M ^ p := pow_nonneg hM_nn p
    have hMp1_nn : 0 ≤ M ^ (p + 1) := pow_nonneg hM_nn (p + 1)
    have hQ_ub : Q_pow (p + 1) r s ≤ ((p + 1 : ℕ) : ℝ) * M ^ p :=
      Q_pow_upper p r s M hr hrM hs hsM
    have h_aR_raw : a * R_pow (p + 1) a r s
                  ≤ M * (((p + 1 : ℕ) : ℝ) ^ 2 * M ^ p) :=
      mul_le_mul haM ih hR_nn hM_nn
    have h_aR_simp : M * (((p + 1 : ℕ) : ℝ) ^ 2 * M ^ p)
                   = ((p + 1 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) := by
      rw [pow_succ]; ring
    have h_aR : a * R_pow (p + 1) a r s
              ≤ ((p + 1 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) :=
      by rw [← h_aR_simp]; exact h_aR_raw
    have hp1_nn : (0 : ℝ) ≤ ((p + 1 : ℕ) : ℝ) := by positivity
    have hM_pow_mono : M ^ p ≤ M ^ (p + 1) := by
      rw [pow_succ]
      have : M ^ p = M ^ p * 1 := (mul_one _).symm
      conv_lhs => rw [this]
      exact mul_le_mul_of_nonneg_left hM1 hMp_nn
    have h_step1 : ((p + 1 : ℕ) : ℝ) * M ^ p
                 ≤ ((p + 1 : ℕ) : ℝ) * M ^ (p + 1) :=
      mul_le_mul_of_nonneg_left hM_pow_mono hp1_nn
    have h_coeff : ((p + 1 : ℕ) : ℝ) + ((p + 1 : ℕ) : ℝ) ^ 2
                 ≤ ((p + 1 + 1 : ℕ) : ℝ) ^ 2 := by
      have hp_nn : (0 : ℝ) ≤ (p : ℝ) := Nat.cast_nonneg p
      have key : ((p + 1 + 1 : ℕ) : ℝ) ^ 2
               = (((p + 1 : ℕ) : ℝ) + ((p + 1 : ℕ) : ℝ) ^ 2) + ((p : ℝ) + 2) := by
        push_cast; ring
      rw [key]
      have h_pos : (0 : ℝ) ≤ (p : ℝ) + 2 :=
        add_nonneg hp_nn (by norm_num)
      exact le_add_of_nonneg_right h_pos
    have h_mul : (((p + 1 : ℕ) : ℝ) + ((p + 1 : ℕ) : ℝ) ^ 2) * M ^ (p + 1)
               ≤ ((p + 1 + 1 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) :=
      mul_le_mul_of_nonneg_right h_coeff hMp1_nn
    calc Q_pow (p + 1) r s + a * R_pow (p + 1) a r s
        ≤ ((p + 1 : ℕ) : ℝ) * M ^ p + ((p + 1 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) :=
          add_le_add hQ_ub h_aR
      _ ≤ ((p + 1 : ℕ) : ℝ) * M ^ (p + 1) + ((p + 1 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) :=
          add_le_add_right h_step1 _
      _ = (((p + 1 : ℕ) : ℝ) + ((p + 1 : ℕ) : ℝ) ^ 2) * M ^ (p + 1) := by ring
      _ ≤ ((p + 1 + 1 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) := h_mul

/-! ## §4. (A) Sharp Quadratic Bound -/

/-- **Theorem A: Sharp Quadratic Bound.**
    For `p ≥ 0`, a compact box `[m, M]` with `0 < m, 1 ≤ M`, and
    `r, s ∈ [m, M]` with `r^(p+2) = X`, the self-anchor step obeys

        |F_{s,p+2}(s) - r| ≤ (p+2) · M^(p+1) / m^(p+1) · (s - r)².

    The `(s-r)²` factor exhibits the quadratic convergence of the
    self-anchor; the constant `(p+2) · (M/m)^(p+1)` is explicit and
    depends only on the box geometry, not on `r` or `s`. -/
theorem sharp_quadratic_bound (p : ℕ) (X r s m M : ℝ)
    (hm : 0 < m) (hM1 : 1 ≤ M)
    (hmr : m ≤ r) (hrM : r ≤ M) (hms : m ≤ s) (hsM : s ≤ M)
    (hX : r ^ (p + 2) = X) :
    |anchor_step_p (p + 2) X s s - r|
      ≤ ((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1) * (s - r) ^ 2 := by
  have hm_nn : 0 ≤ m := le_of_lt hm
  have hs_pos : 0 < s := lt_of_lt_of_le hm hms
  have hr_pos : 0 < r := lt_of_lt_of_le hm hmr
  have hs_nn : 0 ≤ s := le_of_lt hs_pos
  have hr_nn : 0 ≤ r := le_of_lt hr_pos
  -- Self-anchor quadratic identity
  have h_self : anchor_step_p (p + 2) X s s - r
              = (s - r) ^ 2 * R_pow (p + 2) s r s / Q_pow (p + 2) s s :=
    anchor_step_p_self_quadratic (p + 2) (by omega) X r s hX hs_pos
  -- Diagonal value of Q_pow at (s, s)
  have hQ_diag : Q_pow (p + 2) s s = ((p + 2 : ℕ) : ℝ) * s ^ (p + 1) :=
    Q_pow_at_diag (p + 1) s
  -- Polynomial bounds
  have hR_ub : R_pow (p + 2) s r s ≤ ((p + 2 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) :=
    R_pow_upper (p + 1) s r s M hs_nn hsM hr_nn hrM hs_nn hsM hM1
  have hR_nn : 0 ≤ R_pow (p + 2) s r s :=
    R_pow_nonneg _ s r s hs_nn hr_nn hs_nn
  have hmp1_pos : 0 < m ^ (p + 1) := pow_pos hm (p + 1)
  have hsp1_pos : 0 < s ^ (p + 1) := pow_pos hs_pos (p + 1)
  have hsp1_ge : m ^ (p + 1) ≤ s ^ (p + 1) :=
    pow_le_pow_left hm_nn hms (p + 1)
  have hpp2_pos : (0 : ℝ) < ((p + 2 : ℕ) : ℝ) := by positivity
  have hden_pos : 0 < ((p + 2 : ℕ) : ℝ) * s ^ (p + 1) :=
    mul_pos hpp2_pos hsp1_pos
  -- R / Q ratio bound
  have h_ratio : R_pow (p + 2) s r s / (((p + 2 : ℕ) : ℝ) * s ^ (p + 1))
               ≤ ((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1) := by
    rw [div_le_div_iff hden_pos hmp1_pos]
    have hpp2_sq_Mp1_nn : 0 ≤ ((p + 2 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) := by positivity
    calc R_pow (p + 2) s r s * m ^ (p + 1)
        ≤ ((p + 2 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) * m ^ (p + 1) :=
          mul_le_mul_of_nonneg_right hR_ub (le_of_lt hmp1_pos)
      _ ≤ ((p + 2 : ℕ) : ℝ) ^ 2 * M ^ (p + 1) * s ^ (p + 1) :=
          mul_le_mul_of_nonneg_left hsp1_ge hpp2_sq_Mp1_nn
      _ = ((p + 2 : ℕ) : ℝ) * M ^ (p + 1)
            * (((p + 2 : ℕ) : ℝ) * s ^ (p + 1)) := by ring
  -- Take absolute values and chain
  rw [h_self, hQ_diag]
  rw [abs_div, abs_of_pos hden_pos, abs_mul,
      abs_of_nonneg (sq_nonneg (s - r)), abs_of_nonneg hR_nn]
  calc (s - r) ^ 2 * R_pow (p + 2) s r s / (((p + 2 : ℕ) : ℝ) * s ^ (p + 1))
      = (s - r) ^ 2
          * (R_pow (p + 2) s r s / (((p + 2 : ℕ) : ℝ) * s ^ (p + 1))) :=
        mul_div_assoc _ _ _
    _ ≤ (s - r) ^ 2 * (((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1)) :=
        mul_le_mul_of_nonneg_left h_ratio (sq_nonneg _)
    _ = ((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1) * (s - r) ^ 2 := by ring

/-! ## §5. (D) Optimal Anchor Bound -/

/-- **Theorem D: Optimal Anchor Bound.**
    For `p ≥ 0`, compact box `[m, M]` with `0 < m, 1 ≤ M`, and
    `a, r, s ∈ [m, M]` with `r^(p+1) = X`:

        |F_{a,p+1}(s) - r| ≤ (p+1) · M^p / m^p · |a - r| · |s - r|.

    The bilinear factor `|a - r| · |s - r|` makes explicit the
    *optimality* of the anchor choice `a = r`: that choice drives the
    bound to zero, recovering the exactness of the anchor-step at
    the root. -/
theorem optimal_anchor_bound (p : ℕ) (X a r s m M : ℝ)
    (hm : 0 < m) (hM1 : 1 ≤ M)
    (hma : m ≤ a) (haM : a ≤ M)
    (hmr : m ≤ r) (hrM : r ≤ M)
    (hms : m ≤ s) (hsM : s ≤ M)
    (hX : r ^ (p + 1) = X) :
    |anchor_step_p (p + 1) X a s - r|
      ≤ ((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p * (|a - r| * |s - r|) := by
  have hm_nn : 0 ≤ m := le_of_lt hm
  have ha_pos : 0 < a := lt_of_lt_of_le hm hma
  have hr_pos : 0 < r := lt_of_lt_of_le hm hmr
  have hs_pos : 0 < s := lt_of_lt_of_le hm hms
  have ha_nn : 0 ≤ a := le_of_lt ha_pos
  have hr_nn : 0 ≤ r := le_of_lt hr_pos
  have hs_nn : 0 ≤ s := le_of_lt hs_pos
  have hM_pos : 0 < M := lt_of_lt_of_le zero_lt_one hM1
  have hM_nn : 0 ≤ M := le_of_lt hM_pos
  -- Q_pow positive ⇒ master identity applies
  have hQ_pos : 0 < Q_pow (p + 1) a s :=
    Q_pow_pos (p + 1) (by omega) a s ha_pos hs_pos
  have hQ_ne : Q_pow (p + 1) a s ≠ 0 := ne_of_gt hQ_pos
  have h_master : anchor_step_p (p + 1) X a s - r
                = (a - r) * (s - r) * R_pow (p + 1) a r s
                    / Q_pow (p + 1) a s :=
    anchor_step_p_master_identity (p + 1) X a r s hX hQ_ne
  -- Polynomial bounds
  have hQ_lb : ((p + 1 : ℕ) : ℝ) * m ^ p ≤ Q_pow (p + 1) a s :=
    Q_pow_lower p a s m hm_nn hma hms
  have hR_ub : R_pow (p + 1) a r s ≤ ((p + 1 : ℕ) : ℝ) ^ 2 * M ^ p :=
    R_pow_upper p a r s M ha_nn haM hr_nn hrM hs_nn hsM hM1
  have hR_nn : 0 ≤ R_pow (p + 1) a r s := R_pow_nonneg _ a r s ha_nn hr_nn hs_nn
  have hmp_pos : 0 < m ^ p := pow_pos hm p
  have hpp1_pos : (0 : ℝ) < ((p + 1 : ℕ) : ℝ) := by positivity
  -- Ratio bound
  have h_ratio : R_pow (p + 1) a r s / Q_pow (p + 1) a s
               ≤ ((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p := by
    rw [div_le_div_iff hQ_pos hmp_pos]
    have hMp_nn : 0 ≤ M ^ p := pow_nonneg hM_nn p
    have hmul_nn : 0 ≤ ((p + 1 : ℕ) : ℝ) * M ^ p :=
      mul_nonneg (le_of_lt hpp1_pos) hMp_nn
    calc R_pow (p + 1) a r s * m ^ p
        ≤ ((p + 1 : ℕ) : ℝ) ^ 2 * M ^ p * m ^ p :=
          mul_le_mul_of_nonneg_right hR_ub (le_of_lt hmp_pos)
      _ = ((p + 1 : ℕ) : ℝ) * M ^ p * (((p + 1 : ℕ) : ℝ) * m ^ p) := by ring
      _ ≤ ((p + 1 : ℕ) : ℝ) * M ^ p * Q_pow (p + 1) a s :=
          mul_le_mul_of_nonneg_left hQ_lb hmul_nn
  -- Take absolute values and chain
  rw [h_master]
  rw [abs_div, abs_of_pos hQ_pos, abs_mul, abs_mul, abs_of_nonneg hR_nn]
  have h_abs_prod_nn : 0 ≤ |a - r| * |s - r| :=
    mul_nonneg (abs_nonneg _) (abs_nonneg _)
  calc |a - r| * |s - r| * R_pow (p + 1) a r s / Q_pow (p + 1) a s
      = |a - r| * |s - r| * (R_pow (p + 1) a r s / Q_pow (p + 1) a s) :=
        mul_div_assoc _ _ _
    _ ≤ |a - r| * |s - r| * (((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p) :=
        mul_le_mul_of_nonneg_left h_ratio h_abs_prod_nn
    _ = ((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p * (|a - r| * |s - r|) := by ring

/-! ## §6. Sharp Bounds Certificate -/

/-- **Bundled certificate** for Theorems A and D. A single statement
    exposing both explicit-constant convergence bounds. -/
theorem sharp_bounds_certificate :
    (∀ (p : ℕ) (X r s m M : ℝ),
      0 < m → 1 ≤ M → m ≤ r → r ≤ M → m ≤ s → s ≤ M →
      r ^ (p + 2) = X →
      |anchor_step_p (p + 2) X s s - r|
        ≤ ((p + 2 : ℕ) : ℝ) * M ^ (p + 1) / m ^ (p + 1) * (s - r) ^ 2) ∧
    (∀ (p : ℕ) (X a r s m M : ℝ),
      0 < m → 1 ≤ M → m ≤ a → a ≤ M → m ≤ r → r ≤ M → m ≤ s → s ≤ M →
      r ^ (p + 1) = X →
      |anchor_step_p (p + 1) X a s - r|
        ≤ ((p + 1 : ℕ) : ℝ) * M ^ p / m ^ p * (|a - r| * |s - r|)) := by
  refine ⟨?_, ?_⟩
  · intros p X r s m M hm hM1 hmr hrM hms hsM hX
    exact sharp_quadratic_bound p X r s m M hm hM1 hmr hrM hms hsM hX
  · intros p X a r s m M hm hM1 hma haM hmr hrM hms hsM hX
    exact optimal_anchor_bound p X a r s m M hm hM1 hma haM hmr hrM hms hsM hX

end MultiStartPSharp

end Pandrosion
