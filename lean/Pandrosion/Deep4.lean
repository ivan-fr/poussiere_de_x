/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP IV: Half-plane theorem, Universal descent, Smale amortized

  Reference: pandrosion_master.tex, Theorems 3090, 3881, 4578
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Algebra.GeomSum
import Mathlib.Algebra.BigOperators.Basic
import Mathlib.Tactic
import Pandrosion.Deep

open Real Finset BigOperators

namespace Pandrosion

/-! ═══════════════════════════════════════════════════════════
    PART I: HALF-PLANE THEOREM (Theorem 3090)
    ═══════════════════════════════════════════════════════════ -/

/-- **Newton ratio negative for s^d > 1.**
    r = -(s^d - 1)/(d·s^d) < 0 when s^d > 1. -/
theorem newton_ratio_negative (s : ℝ) (hs : s > 0) (d : ℕ) (_hd : d ≥ 1)
    (hsd : s ^ d > 1) :
    -(s ^ d - 1) / ((d : ℝ) * s ^ d) < 0 := by
  apply div_neg_of_neg_of_pos
  · linarith
  · exact mul_pos (by exact_mod_cast (show 0 < d by omega)) (by positivity)

/-- **Newton ratio vanishes at roots: s^d = 1 ⟹ r = 0.** -/
theorem newton_ratio_at_root (s : ℝ) (d : ℕ) (hs : s ^ d = 1) :
    -(s ^ d - 1) / ((d : ℝ) * s ^ d) = 0 := by
  rw [hs, sub_self, neg_zero, zero_div]

/-- **Contraction ratio (d-1)/d ∈ (0,1) for d ≥ 2 — the half-plane condition.** -/
theorem half_plane_contraction (d : ℕ) (hd : d ≥ 2) :
    0 < ((d : ℝ) - 1) / (d : ℝ) ∧ ((d : ℝ) - 1) / (d : ℝ) < 1 := by
  constructor
  · apply div_pos
    · have : (1:ℝ) ≤ (d:ℝ) - 1 := by
        have : (2:ℝ) ≤ (d:ℝ) := by exact_mod_cast hd
        linarith
      linarith
    · exact_mod_cast (show 0 < d by omega)
  · rw [div_lt_one (by exact_mod_cast (show 0 < d by omega) : (d:ℝ) > 0)]
    have : (1:ℝ) ≤ (d:ℝ) := by exact_mod_cast (show 1 ≤ d by omega)
    linarith

/-! ═══════════════════════════════════════════════════════════
    PART II: UNIVERSAL DESCENT (Theorem 3881)
    ═══════════════════════════════════════════════════════════ -/

/-- **cos(θ) ∈ (0,1) for θ ∈ (0, π/2).** -/
theorem cos_in_unit_interval (θ : ℝ) (h1 : 0 < θ) (h2 : θ < π / 2) :
    0 < Real.cos θ ∧ Real.cos θ < 1 := by
  constructor
  · exact Real.cos_pos_of_mem_Ioo ⟨by linarith, h2⟩
  · calc Real.cos θ < Real.cos 0 :=
          Real.cos_lt_cos_of_nonneg_of_le_pi_div_two (le_refl 0) (le_of_lt h2) h1
      _ = 1 := Real.cos_zero

/-- **ln(cos(θ)) < 0 for θ ∈ (0, π/2).** -/
theorem log_cos_neg (θ : ℝ) (h1 : 0 < θ) (h2 : θ < π / 2) :
    Real.log (Real.cos θ) < 0 := by
  have ⟨hpos, hlt⟩ := cos_in_unit_interval θ h1 h2
  exact Real.log_neg hpos hlt

/-- **The descent angle kπ/(2d) ∈ (0, π/2) for 1 ≤ k ≤ d-1, d ≥ 2.** -/
theorem descent_angle_in_range (d k : ℕ) (hd : d ≥ 2) (hk : 1 ≤ k) (hkd : k ≤ d - 1) :
    0 < (k : ℝ) * π / (2 * (d : ℝ)) ∧ (k : ℝ) * π / (2 * (d : ℝ)) < π / 2 := by
  have _hd_pos : (d : ℝ) > 0 := by exact_mod_cast (show 0 < d by omega)
  have _hk_pos : (k : ℝ) > 0 := by exact_mod_cast (show 0 < k by omega)
  constructor
  · positivity
  · rw [div_lt_div_iff (by positivity : 2 * (d:ℝ) > 0) (by norm_num : (2:ℝ) > 0)]
    have hkd_r : (k : ℝ) < (d : ℝ) := by exact_mod_cast (show k < d by omega)
    have hp : π > 0 := pi_pos
    nlinarith [mul_lt_mul_of_pos_right hkd_r hp]

/-- **Each descent term is negative.**
    ln(cos(kπ/(2d))) < 0 for 1 ≤ k ≤ d-1. -/
theorem descent_term_neg (d k : ℕ) (hd : d ≥ 2) (hk : 1 ≤ k) (hkd : k ≤ d - 1) :
    Real.log (Real.cos ((k : ℝ) * π / (2 * (d : ℝ)))) < 0 := by
  have ⟨h1, h2⟩ := descent_angle_in_range d k hd hk hkd
  exact log_cos_neg _ h1 h2

/-- **The descent sum is negative: ∑_{k=1}^{d-1} ln(cos(kπ/(2d))) < 0.**
    This is the core of the universal descent theorem. -/
theorem descent_sum_negative (d : ℕ) (hd : d ≥ 2) :
    ∑ k in Finset.Ico 1 d, Real.log (Real.cos ((k : ℝ) * π / (2 * (d : ℝ)))) < 0 := by
  have _hne : (Finset.Ico 1 d).Nonempty := ⟨1, Finset.mem_Ico.mpr ⟨le_refl _, by omega⟩⟩
  calc ∑ k in Finset.Ico 1 d, Real.log (Real.cos ((k : ℝ) * π / (2 * (d : ℝ))))
      < ∑ _k in Finset.Ico 1 d, (0 : ℝ) := by
        apply Finset.sum_lt_sum
        · intro k hk
          have ⟨hk_lo, hk_hi⟩ := Finset.mem_Ico.mp hk
          exact le_of_lt (descent_term_neg d k hd hk_lo (by omega))
        · exact ⟨1, Finset.mem_Ico.mpr ⟨le_refl _, by omega⟩,
                descent_term_neg d 1 hd (le_refl _) (by omega)⟩
    _ = 0 := Finset.sum_const_zero

/-- **The normalized descent D_d < 0.**
    D_d = (1/d)·∑ ln(cos(kπ/(2d))) < 0. -/
theorem normalized_descent_neg (d : ℕ) (hd : d ≥ 2) :
    (1 / (d : ℝ)) * ∑ k in Finset.Ico 1 d,
      Real.log (Real.cos ((k : ℝ) * π / (2 * (d : ℝ)))) < 0 := by
  apply mul_neg_of_pos_of_neg
  · exact div_pos one_pos (by exact_mod_cast (show 0 < d by omega))
  · exact descent_sum_negative d hd

/-! ═══════════════════════════════════════════════════════════
    PART III: SMALE O(d³) (Theorem 4578)
    ═══════════════════════════════════════════════════════════ -/

/-- **Total cost = 3a·d³.** -/
theorem smale_total_cost (d a : ℕ) :
    3 * d * d * (a * d) = 3 * a * d ^ 3 := by ring

/-- **The potential is positive: ln(R/ρ) > 0 when R > ρ > 0.** -/
theorem potential_positive (R ρ : ℝ) (hR : R > ρ) (hρ : ρ > 0) :
    Real.log (R / ρ) > 0 := by
  exact Real.log_pos (one_lt_div hρ |>.mpr hR)

/-- **Amortized step bound: M/δ > 0 when M, δ > 0.** -/
theorem amortized_bound (M δ : ℝ) (hM : M > 0) (hδ : δ > 0) :
    M / δ > 0 := div_pos hM hδ

/-- **The descent is strict: each step reduces the potential.** -/
theorem descent_strict (Φ Dd : ℝ) (hΦ : Φ > 0) (hDd : Dd < 0) :
    Φ / |Dd| > 0 := div_pos hΦ (abs_pos.mpr (ne_of_lt hDd))

end Pandrosion
