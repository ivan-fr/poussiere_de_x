/-
  Universitas Pandrosion — Legacy / Roots of unity and DFT.

  Two layers of `p`-th-roots-of-unity mathematics:

    §L33  The primitive `p`-th root `ω_p := exp(2πi/p)`; proofs that
          `ω_p^p = 1`, `ω_p ≠ 1` for `p ≥ 2`, and the geometric-sum
          consequence `Σ_{k<p} ω_p^k = 0`.
    §L34  Real-part sum `Σ_{k=1}^{p-1} Re(ω_p^k) = −1` and the
          half-plane corollary `(1/(p−1)) · Σ_{k=1}^{p-1} Re(ω_p^k) < 0`.
    §L35  Orthogonality of characters:
                 `Σ_{j<p} ω_p^{j·m} = p`  when  `p ∣ m`,
                 `Σ_{j<p} ω_p^{j·m} = 0`  when  `p ∤ m`.
          DFT matrix entry, forward/inverse DFT, periodicity,
          eigenvalue product.

  Ported from `ComplexSpectralDescent.lean` + `DFTDecomposition.lean`
  on the `main_previous` branch.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Data.Complex.Basic
import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Algebra.GeomSum
import Mathlib.Tactic

namespace Pandrosion

open Finset BigOperators Real Complex Filter

/-! ## §L33. The primitive `p`-th root of unity -/

/-- The primitive `p`-th root `ω_p = exp(2πi/p)`. -/
noncomputable def omega (p : ℕ) : ℂ :=
  Complex.exp (2 * Real.pi * Complex.I / (p : ℂ))

/-- `ω_p^p = 1`. -/
theorem omega_pow_eq_one (p : ℕ) (hp : p ≥ 1) : omega p ^ p = 1 := by
  unfold omega
  rw [← Complex.exp_nat_mul]
  have : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  rw [show (p : ℂ) * (2 * ↑π * I / (p : ℂ)) = 2 * ↑π * I from by
    field_simp; ring]
  exact Complex.exp_two_pi_mul_I

/-- `ω_p ≠ 1` for `p ≥ 2`. -/
theorem omega_ne_one (p : ℕ) (hp : p ≥ 2) : omega p ≠ 1 := by
  unfold omega
  intro heq
  rw [Complex.exp_eq_one_iff] at heq
  obtain ⟨n, hn⟩ := heq
  have hpi : (2 : ℂ) * ↑π * I ≠ 0 := by
    apply mul_ne_zero (mul_ne_zero (by norm_num) _) Complex.I_ne_zero
    exact_mod_cast Real.pi_ne_zero
  have _hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hstep :
      (2 * ↑π * I) * (1 / (p : ℂ)) = (2 * ↑π * I) * (n : ℂ) := by
    rw [mul_one_div]
    rw [show (↑n : ℂ) * (2 * ↑π * I) = (2 * ↑π * I) * ↑n from by ring] at hn
    exact hn
  have h_eq : (1 : ℂ) / (p : ℂ) = (n : ℂ) := mul_left_cancel₀ hpi hstep
  have h_real : (1 : ℝ) / (p : ℝ) = (n : ℝ) := by
    have h1 : ((1 : ℂ) / (p : ℂ)).re = (1 : ℝ) / (p : ℝ) := by
      push_cast
      simp [Complex.ofReal_div, Complex.ofReal_re]
    have h2 : ((n : ℂ) : ℂ).re = (n : ℝ) := Complex.int_cast_re n
    rw [← h1, ← h2]
    exact congr_arg Complex.re h_eq
  have h_pos : (0 : ℝ) < 1 / (p : ℝ) := by positivity
  have h_half : (1 : ℝ) / (p : ℝ) ≤ 1 / 2 := by
    apply div_le_div_of_nonneg_left (le_of_lt one_pos) (by norm_num)
    exact_mod_cast hp
  have h_nge : (n : ℝ) ≥ 1 := by
    have hnZ : n ≥ 1 := by
      by_contra hlt
      push_neg at hlt
      have : (n : ℝ) ≤ 0 := by
        have : n ≤ 0 := by omega
        exact_mod_cast this
      linarith
    exact_mod_cast hnZ
  linarith

/-- `Σ_{k<p} ω_p^k = 0` for `p ≥ 2`. -/
theorem roots_of_unity_sum (p : ℕ) (hp : p ≥ 2) :
    ∑ k in range p, omega p ^ k = 0 := by
  have hne : omega p - 1 ≠ 0 := sub_ne_zero.mpr (omega_ne_one p hp)
  have hpow : omega p ^ p = 1 := omega_pow_eq_one p (by omega)
  have hgm := geom_sum_mul (omega p) p
  rw [hpow, sub_self] at hgm
  exact (mul_eq_zero.mp hgm).resolve_right hne

/-! ## §L34. Real-part sum and left-half-plane corollary -/

/-- `Re(Σ ω_p^k) = 0`. -/
theorem re_roots_sum_zero (p : ℕ) (hp : p ≥ 2) :
    (∑ k in range p, omega p ^ k).re = 0 := by
  rw [roots_of_unity_sum p hp]; simp

/-- `Σ_{k=1}^{p-1} Re(ω_p^k) = −1`. -/
theorem re_nontrivial_sum (p : ℕ) (hp : p ≥ 2) :
    (∑ k in Finset.Ico 1 p, (omega p ^ k).re) = -1 := by
  have htotal : (∑ k in range p, (omega p ^ k).re) = 0 := by
    rw [show (∑ k in range p, (omega p ^ k).re) =
        (∑ k in range p, omega p ^ k).re from
      (map_sum Complex.reAddGroupHom _ _).symm]
    exact re_roots_sum_zero p hp
  have hsplit : ∑ k in range p, (omega p ^ k).re =
      (omega p ^ 0).re + ∑ k in Finset.Ico 1 p, (omega p ^ k).re := by
    have h0 : range p = insert 0 (Finset.Ico 1 p) := by
      ext k; simp [Finset.mem_range, Finset.mem_Ico, Finset.mem_insert]; omega
    rw [h0, Finset.sum_insert (by simp [Finset.mem_Ico])]
  rw [pow_zero, Complex.one_re] at hsplit
  linarith

/-- The averaged non-trivial real part sits strictly in the left
    half-plane for `p ≥ 2`. -/
theorem spectral_left_half_plane (p : ℕ) (hp : p ≥ 2) :
    (1 / ((p : ℝ) - 1)) * (∑ k in Finset.Ico 1 p, (omega p ^ k).re) < 0 := by
  rw [re_nontrivial_sum p hp]
  have hp_gt_one : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have h_den : (0 : ℝ) < (p : ℝ) - 1 := by linarith
  have : 1 / ((p : ℝ) - 1) > 0 := div_pos one_pos h_den
  linarith

/-- The closed-form complex descent rate `−1/(p−1) < 0`. -/
theorem complex_descent_neg (p : ℕ) (hp : p ≥ 2) :
    (-1 : ℝ) / ((p : ℝ) - 1) < 0 := by
  apply div_neg_of_neg_of_pos
  · norm_num
  · have : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp
    linarith

/-! ## §L35. DFT and character orthogonality -/

/-- DFT matrix entry `W_{jk} = ω_p^{jk}`. -/
noncomputable def dft_entry (p : ℕ) (j k : ℕ) : ℂ :=
  omega p ^ (j * k)

/-- Forward DFT of a signal `x : ℕ → ℂ` at frequency `k`. -/
noncomputable def dft (p : ℕ) (x : ℕ → ℂ) (k : ℕ) : ℂ :=
  ∑ j in range p, x j * dft_entry p j k

/-- Inverse DFT. -/
noncomputable def idft (p : ℕ) (X : ℕ → ℂ) (j : ℕ) : ℂ :=
  (1 / (p : ℂ)) * ∑ k in range p, X k * (dft_entry p j k)⁻¹

/-- **Character sum, `p ∣ m` case.** `Σ_{j<p} ω_p^{j·m} = p`. -/
theorem character_sum_dvd (p : ℕ) (m : ℕ) (_hp : p ≥ 2) (hdvd : p ∣ m) :
    ∑ j in range p, omega p ^ (j * m) = (p : ℂ) := by
  have h1 : ∀ j, omega p ^ (j * m) = 1 := by
    intro j
    obtain ⟨k, rfl⟩ := hdvd
    rw [show j * (p * k) = p * (j * k) from by ring]
    rw [pow_mul, omega_pow_eq_one p (by omega : p ≥ 1)]
    simp
  simp only [h1, Finset.sum_const, Finset.card_range, nsmul_eq_mul, mul_one]

/-- Helper: `ω_p^m ≠ 1` when `p ∤ m`. -/
theorem omega_pow_ne_one_of_not_dvd (p m : ℕ) (hp : p ≥ 2) (h : ¬ p ∣ m) :
    omega p ^ m ≠ 1 := by
  unfold omega
  intro heq
  rw [← Complex.exp_nat_mul] at heq
  rw [Complex.exp_eq_one_iff] at heq
  obtain ⟨n, hn⟩ := heq
  have hpi : (2 : ℂ) * ↑Real.pi * Complex.I ≠ 0 := by
    apply mul_ne_zero (mul_ne_zero (by norm_num) _) Complex.I_ne_zero
    exact_mod_cast Real.pi_ne_zero
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have h_eq : (m : ℂ) / (p : ℂ) = (n : ℂ) := by
    have h2 :
        (2 * ↑Real.pi * Complex.I) * ((m : ℂ) / (p : ℂ)) =
        (2 * ↑Real.pi * Complex.I) * (n : ℂ) := by
      rw [mul_div]
      rw [show (↑n : ℂ) * (2 * ↑Real.pi * Complex.I) =
          (2 * ↑Real.pi * Complex.I) * ↑n from by ring] at hn
      rw [show (↑m : ℂ) * (2 * ↑Real.pi * Complex.I / ↑p) =
          (2 * ↑Real.pi * Complex.I) * ↑m / ↑p from by ring] at hn
      exact hn
    exact mul_left_cancel₀ hpi h2
  apply h
  have hmp : (m : ℂ) = (n : ℂ) * (p : ℂ) := by
    have hdc : (m : ℂ) / (p : ℂ) * (p : ℂ) = (n : ℂ) * (p : ℂ) := by
      rw [h_eq]
    rwa [div_mul_cancel₀ _ hp_ne] at hdc
  have hmp_int : (m : ℤ) = n * (p : ℤ) := by
    have := congr_arg Complex.re hmp
    push_cast at this
    exact_mod_cast this
  have : (p : ℤ) ∣ (m : ℤ) := ⟨n, hmp_int.symm ▸ mul_comm n ↑p⟩
  exact Int.ofNat_dvd.mp this

/-- **Character sum, `p ∤ m` case.** `Σ_{j<p} ω_p^{j·m} = 0`. -/
theorem character_sum_not_dvd (p : ℕ) (m : ℕ) (hp : p ≥ 2) (hndvd : ¬ p ∣ m) :
    ∑ j in range p, omega p ^ (j * m) = 0 := by
  have hpow : (omega p ^ m) ^ p = 1 := by
    rw [← pow_mul]
    rw [show m * p = p * m from by ring]
    rw [pow_mul]
    rw [omega_pow_eq_one p (by omega : p ≥ 1)]
    simp
  have hne : omega p ^ m ≠ 1 := omega_pow_ne_one_of_not_dvd p m hp hndvd
  have hne' : omega p ^ m - 1 ≠ 0 := sub_ne_zero.mpr hne
  have hgm := geom_sum_mul (omega p ^ m) p
  rw [hpow, sub_self] at hgm
  have hsum : ∑ j in range p, (omega p ^ m) ^ j = 0 :=
    (mul_eq_zero.mp hgm).resolve_right hne'
  convert hsum using 1
  congr 1; ext j; rw [← pow_mul, mul_comm]

/-- Self-cancellation along the DFT diagonal. -/
theorem orthogonality_same (p : ℕ) (_hp : p ≥ 2) (k : ℕ) :
    ∑ j in range p, omega p ^ (j * k) * (omega p ^ (j * k))⁻¹ = (p : ℂ) := by
  have hne : ∀ j, omega p ^ (j * k) ≠ 0 :=
    fun j => pow_ne_zero _ (by unfold omega; exact Complex.exp_ne_zero _)
  simp only [mul_inv_cancel (hne _)]
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul, mul_one]

/-! ## §L36. Eigenvalue arithmetic and periodicity -/

/-- `ω_p^{k+p} = ω_p^k`: period-`p` eigenvalues. -/
theorem dft_periodicity (p : ℕ) (hp : p ≥ 1) (k : ℕ) :
    omega p ^ (k + p) = omega p ^ k := by
  rw [pow_add, omega_pow_eq_one p hp, mul_one]

/-- Eigenvalue product: `ω_p^a · ω_p^b = ω_p^{a+b}`. -/
theorem eigenvalue_product (p : ℕ) (a b : ℕ) :
    omega p ^ a * omega p ^ b = omega p ^ (a + b) := by
  rw [← pow_add]

/-- `(ω_p^k)² = ω_p^{2k}`. -/
theorem spectral_eigenvalue_sq (p : ℕ) (k : ℕ) :
    omega p ^ k * omega p ^ k = omega p ^ (2 * k) := by
  rw [← pow_add]; congr 1; ring

end Pandrosion
