/-
  Universitas Pandrosion — Basin / BabylonianContraction
  Split from Advanced2.lean: §I, IV, VI, VII, XI, XIII
  (Banach / Lipschitz / monotone descent / quadratic convergence)
-/
import Mathlib.Topology.MetricSpace.Basic
import Mathlib.Topology.MetricSpace.Contracting
import Mathlib.Topology.Instances.Real
import Mathlib.Topology.Order.Basic
import Mathlib.Tactic
import Pandrosion.Core
import Pandrosion.Diophantine

namespace Pandrosion

section BabylonianContraction

/-! ## §I. Global Basin Invariance (p=2) -/

/-- **One-step basin preservation.** -/
theorem basin_preserved_p2 (x r s : ℝ) (hr : r > 0)
    (hs : s > 0) (h_root : r ^ 2 = x) :
    pandrosion_map 2 x s > 0 ∧
    |pandrosion_map 2 x s - r| ≤ pandrosion_map 2 x s := by
  have hid : pandrosion_map 2 x s - r = (s - r) ^ 2 / (2 * s) :=
    pandrosion_p2_identity x r s hs h_root
  have h2s_pos : (0 : ℝ) < 2 * s := by linarith
  have h_diff_nn : pandrosion_map 2 x s - r ≥ 0 := by
    rw [hid]
    exact div_nonneg (sq_nonneg _) (le_of_lt h2s_pos)
  have hF_ge_r : pandrosion_map 2 x s ≥ r := by linarith
  have hF_pos : pandrosion_map 2 x s > 0 := by linarith
  refine ⟨hF_pos, ?_⟩
  rw [abs_of_nonneg h_diff_nn]
  linarith

/-- **Full basin invariance.** -/
theorem basin_invariant_p2 (x r s₀ : ℝ) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) :
    ∀ n, iterate 2 x s₀ n > 0 ∧
         |iterate 2 x s₀ n - r| ≤ iterate 2 x s₀ n := by
  intro n
  induction n with
  | zero =>
    simp only [iterate]
    have hs₀_pos : s₀ > 0 := by linarith
    refine ⟨hs₀_pos, ?_⟩
    have h_diff_nn : s₀ - r ≥ 0 := by linarith
    rw [abs_of_nonneg h_diff_nn]
    linarith
  | succ n ih =>
    simp only [iterate]
    obtain ⟨hn_pos, _⟩ := ih
    exact basin_preserved_p2 x r (iterate 2 x s₀ n) hr hn_pos h_root

/-- **Unconditional global convergence (p=2).** -/
theorem global_convergence_p2 (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ) :
    error (iterate 2 x s₀ n) r ≤ (1 / 2) ^ n * error s₀ r :=
  error_after_n_steps_p2 x r s₀ hx hr h_root
    (basin_invariant_p2 x r s₀ hr h_s₀ h_root) n

/-- **Limiting form.** -/
theorem global_convergence_eventually_p2
    (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (ε : ℝ) (hε : ε > 0) :
    ∃ N : ℕ, ∀ n ≥ N, error (iterate 2 x s₀ n) r < ε := by
  have h_err₀_nn : error s₀ r ≥ 0 := error_nonneg _ _
  by_cases h_err₀ : error s₀ r = 0
  · refine ⟨0, ?_⟩
    intro n _
    have hbound := global_convergence_p2 x r s₀ hx hr h_s₀ h_root n
    rw [h_err₀] at hbound
    simp at hbound
    linarith
  · have h_err₀_pos : error s₀ r > 0 := lt_of_le_of_ne h_err₀_nn (Ne.symm h_err₀)
    have h_half : (1 / 2 : ℝ) < 1 := by norm_num
    have h_half_nn : (0 : ℝ) ≤ 1 / 2 := by norm_num
    have h_ratio_pos : ε / error s₀ r > 0 := div_pos hε h_err₀_pos
    obtain ⟨N, hN⟩ := exists_pow_lt_of_lt_one h_ratio_pos h_half
    refine ⟨N, ?_⟩
    intro n hn
    have hbound := global_convergence_p2 x r s₀ hx hr h_s₀ h_root n
    have h_pow_le : (1 / 2 : ℝ) ^ n ≤ (1 / 2) ^ N :=
      pow_le_pow_of_le_one h_half_nn (le_of_lt h_half) hn
    calc error (iterate 2 x s₀ n) r
        ≤ (1 / 2) ^ n * error s₀ r := hbound
      _ ≤ (1 / 2) ^ N * error s₀ r :=
          mul_le_mul_of_nonneg_right h_pow_le h_err₀_nn
      _ < (ε / error s₀ r) * error s₀ r :=
          mul_lt_mul_of_pos_right hN h_err₀_pos
      _ = ε := by field_simp

/-! ## §IV. Contraction-Mapping Structure (p=2) -/

/-- **Unconditional Lipschitz-½ contraction inside the basin.** -/
theorem pandrosion_map_lipschitz_half_p2
    (x r s : ℝ) (hx : x > 0) (hr : r > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r) :
    error (pandrosion_map 2 x s) r ≤ (1 / 2) * error s r := by
  have hs_pos : s > 0 := by linarith
  have h_basin_abs : |s - r| ≤ s := by
    rw [abs_of_nonneg (by linarith : s - r ≥ 0)]; linarith
  exact contraction_step_p2 x r s hx hr hs_pos h_root h_basin_abs

/-- **Iterated Lipschitz-½: the basin is a contracting set.** -/
theorem pandrosion_iterate_lipschitz_half_p2
    (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ) :
    error (iterate 2 x s₀ n) r ≤ (1 / 2) ^ n * error s₀ r :=
  global_convergence_p2 x r s₀ hx hr h_s₀ h_root n

/-- **Uniqueness of the fixed point in the basin.** -/
theorem pandrosion_map_fixed_point_unique_p2
    (x r s : ℝ) (hx : x > 0) (hr : r > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r)
    (h_fix : pandrosion_map 2 x s = s) : s = r := by
  have h_lip := pandrosion_map_lipschitz_half_p2 x r s hx hr h_root h_basin
  rw [h_fix] at h_lip
  have h_err_nn : error s r ≥ 0 := error_nonneg _ _
  have h_half : error s r ≤ 0 := by linarith
  have h_err_zero : error s r = 0 := le_antisymm h_half h_err_nn
  unfold error at h_err_zero
  have : s - r = 0 := abs_eq_zero.mp h_err_zero
  linarith

/-! ## §VI. Monotone Descent of the Orbit -/

/-- **One-step descent.** -/
theorem pandrosion_map_le_self_p2 (x r s : ℝ) (hr : r > 0) (hs : s > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r) :
    pandrosion_map 2 x s ≤ s := by
  have hid : pandrosion_map 2 x s - r = (s - r) ^ 2 / (2 * s) :=
    pandrosion_p2_identity x r s hs h_root
  have h2s_pos : (0 : ℝ) < 2 * s := by linarith
  have h_num_bnd : (s - r) ^ 2 ≤ (s - r) * (2 * s) := by nlinarith
  have h_div_bnd : (s - r) ^ 2 / (2 * s) ≤ s - r := by
    rw [div_le_iff h2s_pos]; linarith
  linarith

/-- **Every iterate is ≥ r** (stronger than basin invariance). -/
theorem iterate_ge_r_p2 (x r s₀ : ℝ) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ) :
    iterate 2 x s₀ n ≥ r := by
  induction n with
  | zero => simp only [iterate]; exact h_s₀
  | succ n ih =>
    simp only [iterate]
    have hn_pos : iterate 2 x s₀ n > 0 := by linarith
    have hid : pandrosion_map 2 x (iterate 2 x s₀ n) - r
             = (iterate 2 x s₀ n - r) ^ 2 / (2 * iterate 2 x s₀ n) :=
      pandrosion_p2_identity x r (iterate 2 x s₀ n) hn_pos h_root
    have h2_pos : 2 * iterate 2 x s₀ n > 0 := by linarith
    have h_nn : (iterate 2 x s₀ n - r) ^ 2 / (2 * iterate 2 x s₀ n) ≥ 0 :=
      div_nonneg (sq_nonneg _) (le_of_lt h2_pos)
    linarith

/-- **Full orbit monotone descent.** -/
theorem iterate_monotone_descent_p2 (x r s₀ : ℝ) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ) :
    iterate 2 x s₀ (n + 1) ≤ iterate 2 x s₀ n := by
  have h_ge_r := iterate_ge_r_p2 x r s₀ hr h_s₀ h_root n
  have hn_pos : iterate 2 x s₀ n > 0 := by linarith
  show pandrosion_map 2 x (iterate 2 x s₀ n) ≤ iterate 2 x s₀ n
  exact pandrosion_map_le_self_p2 x r (iterate 2 x s₀ n) hr hn_pos h_root h_ge_r

/-! ## §VII. Quadratic Convergence -/

/-- **Quadratic error bound.** -/
theorem pandrosion_quadratic_convergence_p2 (x r s : ℝ) (hr : r > 0) (hs : s > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r) :
    |pandrosion_map 2 x s - r| ≤ (s - r) ^ 2 / (2 * r) := by
  have hid : pandrosion_map 2 x s - r = (s - r) ^ 2 / (2 * s) :=
    pandrosion_p2_identity x r s hs h_root
  rw [hid]
  have h2s_pos : (0 : ℝ) < 2 * s := by linarith
  have h2r_pos : (0 : ℝ) < 2 * r := by linarith
  have h_sqr_nn : 0 ≤ (s - r) ^ 2 := sq_nonneg _
  have h_diff_nn : (s - r) ^ 2 / (2 * s) ≥ 0 :=
    div_nonneg h_sqr_nn (le_of_lt h2s_pos)
  rw [abs_of_nonneg h_diff_nn]
  rw [div_le_div_iff h2s_pos h2r_pos]
  nlinarith [sq_nonneg (s - r)]

/-- **Iterated quadratic bound (single-step form).** -/
theorem pandrosion_error_squared_p2 (x r s : ℝ) (hr : r > 0) (hs : s > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r) :
    error (pandrosion_map 2 x s) r ≤ (error s r) ^ 2 / (2 * r) := by
  unfold error
  have h_sr_nn : s - r ≥ 0 := by linarith
  have : |s - r| = s - r := abs_of_nonneg h_sr_nn
  rw [this]
  exact pandrosion_quadratic_convergence_p2 x r s hr hs h_root h_basin

/-! ## §XI. Pairwise Banach Contraction -/

/-- **The difference identity for `pandrosion_map 2`.** -/
theorem pandrosion_map_diff_p2 (x s t : ℝ) (hs : s > 0) (ht : t > 0) :
    pandrosion_map 2 x s - pandrosion_map 2 x t =
    (s - t) * (s * t - x) / (2 * (s * t)) := by
  unfold pandrosion_map
  simp only
  have _hs_ne : s ≠ 0 := ne_of_gt hs
  have _ht_ne : t ≠ 0 := ne_of_gt ht
  have hs2_ne : (2 : ℝ) * s ^ 2 ≠ 0 := by positivity
  have ht2_ne : (2 : ℝ) * t ^ 2 ≠ 0 := by positivity
  have hden_s_ne :
      (↑(2 : ℕ) : ℝ) * s ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x ≠ 0 := by
    rw [show (↑(2 : ℕ) : ℝ) * s ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x =
          2 * s ^ 2 from by push_cast; ring]
    exact hs2_ne
  have hden_t_ne :
      (↑(2 : ℕ) : ℝ) * t ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x ≠ 0 := by
    rw [show (↑(2 : ℕ) : ℝ) * t ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x =
          2 * t ^ 2 from by push_cast; ring]
    exact ht2_ne
  rw [if_neg hden_s_ne, if_neg hden_t_ne]
  have h_sden :
      (↑(2 : ℕ) : ℝ) * s ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x = 2 * s ^ 2 := by
    push_cast; ring
  have h_tden :
      (↑(2 : ℕ) : ℝ) * t ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x = 2 * t ^ 2 := by
    push_cast; ring
  rw [h_sden, h_tden]
  field_simp
  push_cast
  ring

/-- **Pairwise Lipschitz-½ on the basin — full Banach contraction.** -/
theorem pandrosion_map_lipschitz_pairwise_p2
    (x r s t : ℝ) (hr : r > 0) (h_root : r ^ 2 = x)
    (hs : s ≥ r) (ht : t ≥ r) :
    |pandrosion_map 2 x s - pandrosion_map 2 x t| ≤ (1 / 2) * |s - t| := by
  have hs_pos : s > 0 := by linarith
  have ht_pos : t > 0 := by linarith
  have hx_nn : 0 ≤ x := by rw [← h_root]; exact sq_nonneg _
  have h_diff :
      pandrosion_map 2 x s - pandrosion_map 2 x t
        = (s - t) * (s * t - x) / (2 * (s * t)) :=
    pandrosion_map_diff_p2 x s t hs_pos ht_pos
  rw [h_diff]
  have hst_pos : 0 < s * t := mul_pos hs_pos ht_pos
  have h2st_pos : (0 : ℝ) < 2 * (s * t) := by linarith
  have hst_ge_x : x ≤ s * t := by
    rw [← h_root]
    have : r * r ≤ s * t := mul_le_mul hs ht (le_of_lt hr) (by linarith)
    nlinarith
  have h_st_x_nn : 0 ≤ s * t - x := by linarith
  rw [abs_div, abs_mul,
      abs_of_pos h2st_pos, abs_of_nonneg h_st_x_nn]
  have h_abs_nn : 0 ≤ |s - t| := abs_nonneg _
  have h_ratio : (s * t - x) / (2 * (s * t)) ≤ 1 / 2 := by
    rw [div_le_div_iff h2st_pos (by norm_num : (0:ℝ) < 2)]
    linarith
  calc |s - t| * (s * t - x) / (2 * (s * t))
      = |s - t| * ((s * t - x) / (2 * (s * t))) := by ring
    _ ≤ |s - t| * (1 / 2) :=
        mul_le_mul_of_nonneg_left h_ratio h_abs_nn
    _ = (1 / 2) * |s - t| := by ring

/-! ## §XIII. Banach Fixed-Point Theorem via Mathlib -/

section BanachFixedPoint

open Topology Filter

/-- **The basin as a subtype.** -/
abbrev PandrosionBasin (r : ℝ) := {s : ℝ // r ≤ s}

/-- The basin contains `r` itself. -/
instance instNonemptyPandrosionBasin (r : ℝ) : Nonempty (PandrosionBasin r) :=
  ⟨⟨r, le_refl r⟩⟩

/-- The basin is a closed subset of ℝ, hence complete. -/
instance instCompletePandrosionBasin (r : ℝ) : CompleteSpace (PandrosionBasin r) := by
  have h : IsClosed (Set.Ici r) := isClosed_Ici
  exact h.completeSpace_coe

/-- The distinguished basin element `r`. -/
def rBasin (r : ℝ) : PandrosionBasin r := ⟨r, le_refl r⟩

/-- **The Babylonian map lifted to the basin.** -/
noncomputable def babylonianOnBasin
    (x r : ℝ) (hr : r > 0) (h_root : r ^ 2 = x) :
    PandrosionBasin r → PandrosionBasin r :=
  fun s => ⟨pandrosion_map 2 x s.val, by
    have hs_pos : 0 < s.val := lt_of_lt_of_le hr s.property
    have hid := pandrosion_p2_identity x r s.val hs_pos h_root
    have h_nn : 0 ≤ (s.val - r) ^ 2 / (2 * s.val) :=
      div_nonneg (sq_nonneg _) (by linarith)
    linarith⟩

/-- **`r` is a fixed point** of the lifted map. -/
theorem babylonianOnBasin_fixed_r
    (x r : ℝ) (hr : r > 0) (h_root : r ^ 2 = x) :
    babylonianOnBasin x r hr h_root (rBasin r) = rBasin r := by
  apply Subtype.ext
  show pandrosion_map 2 x r = r
  have hid := pandrosion_p2_identity x r r hr h_root
  have h_zero : (r - r : ℝ) ^ 2 / (2 * r) = 0 := by
    rw [sub_self, pow_two, zero_mul, zero_div]
  linarith

/-- **Lipschitz-½ on the basin.** -/
theorem babylonianOnBasin_lipschitz
    (x r : ℝ) (hr : r > 0) (h_root : r ^ 2 = x) :
    LipschitzWith (1/2 : NNReal) (babylonianOnBasin x r hr h_root) := by
  apply LipschitzWith.of_dist_le_mul
  intro s t
  show dist (babylonianOnBasin x r hr h_root s).val
            (babylonianOnBasin x r hr h_root t).val
       ≤ ↑(1/2 : NNReal) * dist s.val t.val
  rw [Real.dist_eq, Real.dist_eq]
  have h_lip := pandrosion_map_lipschitz_pairwise_p2 x r s.val t.val
                  hr h_root s.property t.property
  have h_coe : ((1/2 : NNReal) : ℝ) = 1/2 := by
    simp [NNReal.coe_div]
  rw [h_coe]
  exact h_lip

/-- **The Babylonian map is a Mathlib `ContractingWith (1/2)`.** -/
theorem babylonianOnBasin_contracting
    (x r : ℝ) (hr : r > 0) (h_root : r ^ 2 = x) :
    ContractingWith (1/2 : NNReal) (babylonianOnBasin x r hr h_root) := by
  refine ⟨?_, babylonianOnBasin_lipschitz x r hr h_root⟩
  rw [← NNReal.coe_lt_coe]
  push_cast
  norm_num

/-- **Uniqueness of the fixed point (Banach).** -/
theorem babylonianOnBasin_fixed_unique
    (x r : ℝ) (hr : r > 0) (h_root : r ^ 2 = x) (s : PandrosionBasin r)
    (hs : babylonianOnBasin x r hr h_root s = s) : s = rBasin r := by
  have h_r_fixed := babylonianOnBasin_fixed_r x r hr h_root
  have hF_s : pandrosion_map 2 x s.val = s.val := congrArg Subtype.val hs
  have hF_r : pandrosion_map 2 x r = r := congrArg Subtype.val h_r_fixed
  have h_pair := pandrosion_map_lipschitz_pairwise_p2 x r s.val r
                   hr h_root s.property (le_refl r)
  rw [hF_s, hF_r] at h_pair
  have h_abs_nn : 0 ≤ |s.val - r| := abs_nonneg _
  have h_half : |s.val - r| ≤ 0 := by linarith
  have h_zero : |s.val - r| = 0 := le_antisymm h_half h_abs_nn
  have h_val_eq : s.val = r := by
    have := abs_eq_zero.mp h_zero
    linarith
  exact Subtype.ext h_val_eq

/-- **Convergence of iterates to `rBasin` (Banach).** -/
theorem babylonianOnBasin_tendsto_rBasin
    (x r : ℝ) (hr : r > 0) (h_root : r ^ 2 = x) (s : PandrosionBasin r) :
    Tendsto
      (fun n : ℕ => (babylonianOnBasin x r hr h_root)^[n] s)
      atTop (𝓝 (rBasin r)) := by
  set F := babylonianOnBasin x r hr h_root
  have hc : ContractingWith (1/2 : NNReal) F := babylonianOnBasin_contracting x r hr h_root
  have h_edist_ne : edist s (F s) ≠ ⊤ := by
    rw [edist_dist]; exact ENNReal.ofReal_ne_top
  have h_tendsto :
      Tendsto (fun n : ℕ => F^[n] s) atTop
        (𝓝 (ContractingWith.efixedPoint F hc s h_edist_ne)) :=
    ContractingWith.tendsto_iterate_efixedPoint hc h_edist_ne
  have h_fp :=
    ContractingWith.efixedPoint_isFixedPt (f := F) hc h_edist_ne
  have h_eq : ContractingWith.efixedPoint F hc s h_edist_ne = rBasin r :=
    babylonianOnBasin_fixed_unique x r hr h_root _ h_fp
  rw [h_eq] at h_tendsto
  exact h_tendsto

/-- **Existence-and-uniqueness wrapper.** -/
theorem babylonian_banach_fixed_point
    (x r : ℝ) (hr : r > 0) (h_root : r ^ 2 = x) :
    ∃! s : PandrosionBasin r, babylonianOnBasin x r hr h_root s = s := by
  refine ⟨rBasin r, babylonianOnBasin_fixed_r x r hr h_root, ?_⟩
  intro s hs
  exact babylonianOnBasin_fixed_unique x r hr h_root s hs

/-- **Geometric `edist` rate bound on the basin.** -/
theorem babylonianOnBasin_edist_iterate_rate
    (x r : ℝ) (hr : r > 0) (h_root : r ^ 2 = x) (s : PandrosionBasin r) (n : ℕ) :
    edist ((babylonianOnBasin x r hr h_root)^[n] s) (rBasin r)
      ≤ edist s (babylonianOnBasin x r hr h_root s)
          * ((1/2 : NNReal) : ENNReal) ^ n / (1 - (1/2 : NNReal)) := by
  set F := babylonianOnBasin x r hr h_root
  have hc : ContractingWith (1/2 : NNReal) F := babylonianOnBasin_contracting x r hr h_root
  have h_edist_ne : edist s (F s) ≠ ⊤ := by
    rw [edist_dist]; exact ENNReal.ofReal_ne_top
  have h_bound := ContractingWith.apriori_edist_iterate_efixedPoint_le hc h_edist_ne n
  have h_fp :=
    ContractingWith.efixedPoint_isFixedPt (f := F) hc h_edist_ne
  have h_eq : ContractingWith.efixedPoint F hc s h_edist_ne = rBasin r :=
    babylonianOnBasin_fixed_unique x r hr h_root _ h_fp
  rw [h_eq] at h_bound
  exact h_bound

end BanachFixedPoint

end BabylonianContraction

end Pandrosion
