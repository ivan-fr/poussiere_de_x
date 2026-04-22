/-
  Universitas Pandrosion — Legacy / Descent.

  The "universal descent" theorem: for every `p ≥ 2`,
    D_p := (1/p) · ∑_{k<p} log(cos(kπ / (2p)))  <  0.

  Proved via trigonometric bounds on each term:
    • `cos(θ) ∈ (0, 1)` for `θ ∈ (0, π/2)`
    • `log(cos(θ)) < 0` on the same range
    • every angle `kπ/(2p)` with `1 ≤ k ≤ p − 1` lies in `(0, π/2)`
    • the term at `k = 0` contributes `log 1 = 0`, and the rest are
      *strictly* negative.

  Ported from `HalfPlaneUniversal.lean` on the `main_previous`
  branch (Parts I + II). Part III of that file was trivial `d³`
  arithmetic and is omitted.
-/

import Pandrosion.Legacy.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic

namespace Pandrosion

open Real Finset BigOperators

/-! ## §L20. Newton-ratio sign at the cyclotomic fixed point -/

/-- Newton ratio `r = −(s^d − 1)/(d · s^d)` is negative when `s^d > 1`. -/
theorem newton_ratio_negative (s : ℝ) (hs : s > 0) (d : ℕ) (_hd : d ≥ 1)
    (hsd : s ^ d > 1) :
    -(s ^ d - 1) / ((d : ℝ) * s ^ d) < 0 := by
  apply div_neg_of_neg_of_pos
  · linarith
  · exact mul_pos (by exact_mod_cast (show 0 < d by omega)) (by positivity)

/-- Newton ratio vanishes at the root `s^d = 1`. -/
theorem newton_ratio_at_root (s : ℝ) (d : ℕ) (hs : s ^ d = 1) :
    -(s ^ d - 1) / ((d : ℝ) * s ^ d) = 0 := by
  rw [hs, sub_self, neg_zero, zero_div]

/-- The contraction ratio `(d − 1)/d` sits in `(0, 1)` for `d ≥ 2`. -/
theorem half_plane_contraction (d : ℕ) (hd : d ≥ 2) :
    0 < ((d : ℝ) - 1) / (d : ℝ) ∧ ((d : ℝ) - 1) / (d : ℝ) < 1 := by
  refine ⟨?_, ?_⟩
  · apply div_pos
    · have : (2 : ℝ) ≤ (d : ℝ) := by exact_mod_cast hd
      linarith
    · exact_mod_cast (show 0 < d by omega)
  · rw [div_lt_one (by exact_mod_cast (show 0 < d by omega) : (d : ℝ) > 0)]
    have : (1 : ℝ) ≤ (d : ℝ) := by exact_mod_cast (show 1 ≤ d by omega)
    linarith

/-! ## §L21. Trigonometric bounds at the descent angles -/

/-- `0 < cos θ < 1` for `θ ∈ (0, π/2)`. -/
theorem cos_in_unit_interval (θ : ℝ) (h1 : 0 < θ) (h2 : θ < π / 2) :
    0 < Real.cos θ ∧ Real.cos θ < 1 := by
  refine ⟨Real.cos_pos_of_mem_Ioo ⟨by linarith, h2⟩, ?_⟩
  calc Real.cos θ
      < Real.cos 0 :=
        Real.cos_lt_cos_of_nonneg_of_le_pi_div_two (le_refl 0) (le_of_lt h2) h1
    _ = 1 := Real.cos_zero

/-- `log(cos θ) < 0` for `θ ∈ (0, π/2)`. -/
theorem log_cos_neg (θ : ℝ) (h1 : 0 < θ) (h2 : θ < π / 2) :
    Real.log (Real.cos θ) < 0 :=
  let ⟨hpos, hlt⟩ := cos_in_unit_interval θ h1 h2
  Real.log_neg hpos hlt

/-- Every descent angle `kπ/(2p)` with `1 ≤ k ≤ p − 1` lies in `(0, π/2)`. -/
theorem descent_angle_in_range (p k : ℕ) (hp : p ≥ 2) (hk : 1 ≤ k) (hkp : k ≤ p - 1) :
    0 < (k : ℝ) * π / (2 * (p : ℝ)) ∧
    (k : ℝ) * π / (2 * (p : ℝ)) < π / 2 := by
  have _hp_pos : (p : ℝ) > 0 := by exact_mod_cast (show 0 < p by omega)
  have _hk_pos : (k : ℝ) > 0 := by exact_mod_cast (show 0 < k by omega)
  refine ⟨by positivity, ?_⟩
  rw [div_lt_div_iff (by positivity : 2 * (p : ℝ) > 0) (by norm_num : (2 : ℝ) > 0)]
  have hkp_r : (k : ℝ) < (p : ℝ) := by exact_mod_cast (show k < p by omega)
  have hπ : π > 0 := pi_pos
  nlinarith [mul_lt_mul_of_pos_right hkp_r hπ]

/-- Each term `log(cos(kπ/(2p)))` (`1 ≤ k ≤ p − 1`) is strictly negative. -/
theorem descent_term_neg (p k : ℕ) (hp : p ≥ 2) (hk : 1 ≤ k) (hkp : k ≤ p - 1) :
    Real.log (Real.cos ((k : ℝ) * π / (2 * (p : ℝ)))) < 0 :=
  let ⟨h1, h2⟩ := descent_angle_in_range p k hp hk hkp
  log_cos_neg _ h1 h2

/-! ## §L22. The universal descent theorem -/

/-- **Universal descent (sum form).** For every `p ≥ 2`,
    `∑_{k=1}^{p-1} log(cos(kπ/(2p))) < 0`. -/
theorem descent_sum_negative (p : ℕ) (hp : p ≥ 2) :
    ∑ k in Finset.Ico 1 p, Real.log (Real.cos ((k : ℝ) * π / (2 * (p : ℝ)))) < 0 := by
  have h1_mem : (1 : ℕ) ∈ Finset.Ico 1 p :=
    Finset.mem_Ico.mpr ⟨le_refl _, by omega⟩
  calc ∑ k in Finset.Ico 1 p, Real.log (Real.cos ((k : ℝ) * π / (2 * (p : ℝ))))
      < ∑ _k in Finset.Ico 1 p, (0 : ℝ) := by
        apply Finset.sum_lt_sum
        · intro k hk
          have ⟨hk_lo, hk_hi⟩ := Finset.mem_Ico.mp hk
          exact le_of_lt (descent_term_neg p k hp hk_lo (by omega))
        · exact ⟨1, h1_mem, descent_term_neg p 1 hp (le_refl _) (by omega)⟩
    _ = 0 := Finset.sum_const_zero

/-- **Normalised universal descent.** `D_p := (1/p) · ∑_{k=1}^{p-1} log(cos(kπ/(2p))) < 0`. -/
theorem normalized_descent_neg (p : ℕ) (hp : p ≥ 2) :
    (1 / (p : ℝ)) * ∑ k in Finset.Ico 1 p,
      Real.log (Real.cos ((k : ℝ) * π / (2 * (p : ℝ)))) < 0 := by
  apply mul_neg_of_pos_of_neg
  · exact div_pos one_pos (by exact_mod_cast (show 0 < p by omega))
  · exact descent_sum_negative p hp

end Pandrosion
