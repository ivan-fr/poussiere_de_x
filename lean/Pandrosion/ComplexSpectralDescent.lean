/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP XII: COMPLEX SPECTRAL DESCENT — Re(r) < 0

  The p-th roots of unity ω = exp(2πi/p) satisfy:
    Σ_{k=0}^{p-1} ω^k = 0     (from geom_sum with ω^p=1, ω≠1)
  Taking real parts:  Σ cos(2πk/p) = 0
  Subtracting k=0:    Σ_{k=1}^{p-1} cos(2πk/p) = -1
  Average = -1/(p-1) < 0: LEFT HALF-PLANE.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Data.Complex.Basic
import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Algebra.GeomSum
import Mathlib.Tactic

open Finset BigOperators Real Complex Filter

namespace Pandrosion

/-! ## §90. Roots of Unity -/

noncomputable def omega (p : ℕ) : ℂ :=
  Complex.exp (2 * Real.pi * Complex.I / (p : ℂ))

theorem omega_pow_eq_one (p : ℕ) (hp : p ≥ 1) : omega p ^ p = 1 := by
  unfold omega
  rw [← Complex.exp_nat_mul]
  have : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  rw [show (p : ℂ) * (2 * ↑π * I / (p : ℂ)) = 2 * ↑π * I from by field_simp; ring]
  exact Complex.exp_two_pi_mul_I

/-- ω ≠ 1 for p ≥ 2. Uses exp_eq_one_iff: exp(z)=1 ↔ z=2πni. -/
theorem omega_ne_one (p : ℕ) (hp : p ≥ 2) : omega p ≠ 1 := by
  unfold omega
  intro heq
  rw [Complex.exp_eq_one_iff] at heq
  obtain ⟨n, hn⟩ := heq
  -- hn : 2πI/p = n·(2πI)
  -- This gives 1/p = n, impossible for p ≥ 2, n ∈ ℤ
  have hpi : (2 : ℂ) * ↑π * I ≠ 0 := by
    apply mul_ne_zero (mul_ne_zero (by norm_num) _) Complex.I_ne_zero
    exact_mod_cast Real.pi_ne_zero
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  -- Rewrite hn: 2πI * (1/p) = 2πI * n
  have : (2 * ↑π * I) * (1 / (p : ℂ)) = (2 * ↑π * I) * (n : ℂ) := by
    rw [mul_one_div]
    rw [show (↑n : ℂ) * (2 * ↑π * I) = (2 * ↑π * I) * ↑n from by ring] at hn
    exact hn
  have h_eq : (1 : ℂ) / (p : ℂ) = (n : ℂ) := mul_left_cancel₀ hpi this
  -- Take imaginary parts: 0 = 0 (trivial). Take real parts: 1/p = n.
  -- Cast to ℝ: since both sides are real
  have h_real : (1 : ℝ) / (p : ℝ) = (n : ℝ) := by
    have h1 : ((1 : ℂ) / (p : ℂ)).re = (1 : ℝ) / (p : ℝ) := by
      push_cast
      simp [Complex.ofReal_div, Complex.ofReal_re]
    have h2 : ((n : ℂ) : ℂ).re = (n : ℝ) := Complex.int_cast_re n
    rw [← h1, ← h2]
    exact congr_arg Complex.re h_eq
  -- 0 < 1/p ≤ 1/2 but n ∈ ℤ: contradiction
  have : (0 : ℝ) < 1 / (p : ℝ) := by positivity
  have : (1 : ℝ) / (p : ℝ) ≤ 1 / 2 := by
    apply div_le_div_of_nonneg_left (le_of_lt one_pos) (by norm_num)
    exact_mod_cast hp
  have : (n : ℝ) ≥ 1 := by
    have : n ≥ 1 := by
      by_contra h
      push_neg at h
      have : (n : ℝ) ≤ 0 := by
        have : n ≤ 0 := by omega
        exact_mod_cast this
      linarith
    exact_mod_cast this
  linarith

/-- Σ_{k=0}^{p-1} ω^k = 0. -/
theorem roots_of_unity_sum (p : ℕ) (hp : p ≥ 2) :
    ∑ k in range p, omega p ^ k = 0 := by
  have hne : omega p - 1 ≠ 0 := sub_ne_zero.mpr (omega_ne_one p hp)
  have hpow : omega p ^ p = 1 := omega_pow_eq_one p (by omega)
  have := geom_sum_mul (omega p) p
  rw [hpow, sub_self] at this
  exact (mul_eq_zero.mp this).resolve_right hne

/-! ## §91. The Cosine Sum Identity -/

/-- Re(Σ ω^k) = 0. -/
theorem re_roots_sum_zero (p : ℕ) (hp : p ≥ 2) :
    (∑ k in range p, omega p ^ k).re = 0 := by
  rw [roots_of_unity_sum p hp]; simp

/-- **Σ_{k=1}^{p-1} Re(ω^k) = -1.** -/
theorem re_nontrivial_sum (p : ℕ) (hp : p ≥ 2) :
    (∑ k in Finset.Ico 1 p, (omega p ^ k).re) = -1 := by
  have htotal : (∑ k in range p, (omega p ^ k).re) = 0 := by
    rw [show (∑ k in range p, (omega p ^ k).re) =
        (∑ k in range p, omega p ^ k).re from
      (map_sum Complex.reAddGroupHom _ _).symm]
    exact re_roots_sum_zero p hp
  -- range p = {0} ∪ Ico 1 p, sum over {0} = ω^0.re = 1
  have hsplit : ∑ k in range p, (omega p ^ k).re =
      (omega p ^ 0).re + ∑ k in Finset.Ico 1 p, (omega p ^ k).re := by
    have h0 : range p = insert 0 (Finset.Ico 1 p) := by
      ext k; simp [Finset.mem_range, Finset.mem_Ico, Finset.mem_insert]; omega
    rw [h0, Finset.sum_insert (by simp [Finset.mem_Ico])]
  rw [pow_zero, Complex.one_re] at hsplit
  linarith

/-! ## §92. LEFT HALF-PLANE: Re(spectral average) < 0 -/

/-- **The spectral average over non-trivial roots has Re < 0.** -/
theorem spectral_left_half_plane (p : ℕ) (hp : p ≥ 2) :
    (1 / ((p : ℝ) - 1)) * (∑ k in Finset.Ico 1 p, (omega p ^ k).re) < 0 := by
  rw [re_nontrivial_sum p hp]
  have : (0 : ℝ) < (p : ℝ) - 1 := by
    have : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp
    linarith
  have : 1 / ((p : ℝ) - 1) > 0 := div_pos one_pos this
  linarith

/-- **Complex descent rate = -1/(p-1) < 0.** -/
theorem complex_descent_neg (p : ℕ) (hp : p ≥ 2) :
    (-1 : ℝ) / ((p : ℝ) - 1) < 0 := by
  apply div_neg_of_neg_of_pos
  · norm_num
  · have : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp
    linarith

end Pandrosion
