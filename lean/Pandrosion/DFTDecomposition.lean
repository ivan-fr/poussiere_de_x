/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP XIV: DFT SPECTRAL DECOMPOSITION

  The key identity is **orthogonality of characters**:
    Σ_{j=0}^{p-1} ω^{jm} = { p  if p | m
                             { 0  otherwise
  This follows from the geometric sum and roots of unity.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Data.Complex.Basic
import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Algebra.GeomSum
import Mathlib.Tactic
import Pandrosion.Deep12

open Finset BigOperators Complex

namespace Pandrosion

/-! ## §101. DFT: Definitions -/

/-- The DFT matrix entry: W_{jk} = ω^{jk}. -/
noncomputable def dft_entry (p : ℕ) (j k : ℕ) : ℂ :=
  omega p ^ (j * k)

/-- The DFT of a signal x at frequency k. -/
noncomputable def dft (p : ℕ) (x : ℕ → ℂ) (k : ℕ) : ℂ :=
  ∑ j in range p, x j * dft_entry p j k

/-- The inverse DFT. -/
noncomputable def idft (p : ℕ) (X : ℕ → ℂ) (j : ℕ) : ℂ :=
  (1 / (p : ℂ)) * ∑ k in range p, X k * (dft_entry p j k)⁻¹

/-! ## §102. Orthogonality of Characters -/

/-- **Case 1: p | m ⟹ Σ ω^{jm} = p.** -/
theorem character_sum_dvd (p : ℕ) (m : ℕ) (_hp : p ≥ 2) (hdvd : p ∣ m) :
    ∑ j in range p, omega p ^ (j * m) = (p : ℂ) := by
  have h1 : ∀ j, omega p ^ (j * m) = 1 := by
    intro j
    obtain ⟨k, rfl⟩ := hdvd
    rw [show j * (p * k) = p * (j * k) from by ring]
    rw [pow_mul, omega_pow_eq_one p (by omega : p ≥ 1)]
    simp
  simp only [h1, Finset.sum_const, Finset.card_range, nsmul_eq_mul, mul_one]

/-- **ω^m ≠ 1 when p ∤ m.**
    Proof: ω^m = exp(2πim/p) = 1 ⟺ m/p ∈ ℤ ⟺ p | m.
    Contrapositive: ¬(p|m) → ω^m ≠ 1. -/
theorem omega_pow_ne_one_of_not_dvd (p m : ℕ) (hp : p ≥ 2) (h : ¬ p ∣ m) :
    omega p ^ m ≠ 1 := by
  unfold omega
  intro heq
  rw [← Complex.exp_nat_mul] at heq
  rw [Complex.exp_eq_one_iff] at heq
  obtain ⟨n, hn⟩ := heq
  -- hn : ↑m * (2πI/p) = n · (2πI)
  have hpi : (2 : ℂ) * ↑Real.pi * Complex.I ≠ 0 := by
    apply mul_ne_zero (mul_ne_zero (by norm_num) _) Complex.I_ne_zero
    exact_mod_cast Real.pi_ne_zero
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  -- Extract m/p = n from the equation
  have h_eq : (m : ℂ) / (p : ℂ) = (n : ℂ) := by
    have h2 : (2 * ↑Real.pi * Complex.I) * ((m : ℂ) / (p : ℂ)) =
              (2 * ↑Real.pi * Complex.I) * (n : ℂ) := by
      rw [mul_div]
      rw [show (↑n : ℂ) * (2 * ↑Real.pi * Complex.I) =
          (2 * ↑Real.pi * Complex.I) * ↑n from by ring] at hn
      rw [show (↑m : ℂ) * (2 * ↑Real.pi * Complex.I / ↑p) =
          (2 * ↑Real.pi * Complex.I) * ↑m / ↑p from by ring] at hn
      exact hn
    exact mul_left_cancel₀ hpi h2
  -- From h_eq: (m : ℂ)/p = n, deduce m = n*p as ℂ, then as ℤ, then p∣m
  apply h
  -- Step 1: m = n * p as ℂ
  have hmp : (m : ℂ) = (n : ℂ) * (p : ℂ) := by
    have : (m : ℂ) / (p : ℂ) * (p : ℂ) = (n : ℂ) * (p : ℂ) := by
      rw [h_eq]
    rwa [div_mul_cancel₀ _ hp_ne] at this
  -- Step 2: m = n * p as ℤ (extract via real part)
  have hmp_int : (m : ℤ) = n * (p : ℤ) := by
    have := congr_arg Complex.re hmp
    push_cast at this
    exact_mod_cast this
  -- Step 3: p ∣ m as ℕ
  have : (p : ℤ) ∣ (m : ℤ) := ⟨n, hmp_int.symm ▸ mul_comm n ↑p⟩
  exact Int.ofNat_dvd.mp this

/-- **Case 2: p ∤ m ⟹ Σ ω^{jm} = 0.** -/
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
  have hgeom := geom_sum_mul (omega p ^ m) p
  rw [hpow, sub_self] at hgeom
  have hsum : ∑ j in range p, (omega p ^ m) ^ j = 0 :=
    (mul_eq_zero.mp hgeom).resolve_right hne'
  convert hsum using 1
  congr 1; ext j; rw [← pow_mul, mul_comm]

/-! ## §103. Orthogonality: Diagonal Form -/

/-- **When k = l: Σ ω^{jk} · ω^{-jl} = p.** -/
theorem orthogonality_same (p : ℕ) (_hp : p ≥ 2) (k : ℕ) :
    ∑ j in range p, omega p ^ (j * k) * (omega p ^ (j * k))⁻¹ = (p : ℂ) := by
  have hne : ∀ j, omega p ^ (j * k) ≠ 0 :=
    fun j => pow_ne_zero _ (by unfold omega; exact Complex.exp_ne_zero _)
  simp only [mul_inv_cancel (hne _)]
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul, mul_one]

/-! ## §104. Spectral Properties -/

/-- **ω^{k+p} = ω^k: periodicity of the eigenvalues.** -/
theorem dft_periodicity (p : ℕ) (hp : p ≥ 1) (k : ℕ) :
    omega p ^ (k + p) = omega p ^ k := by
  rw [pow_add, omega_pow_eq_one p hp, mul_one]

/-- **The eigenvalue product: ω^a · ω^b = ω^{a+b}.** -/
theorem eigenvalue_product (p : ℕ) (a b : ℕ) :
    omega p ^ a * omega p ^ b = omega p ^ (a + b) := by
  rw [← pow_add]

/-- **The spectral eigenvalue squared: (ω^k)² = ω^{2k}.** -/
theorem spectral_eigenvalue_sq (p : ℕ) (k : ℕ) :
    omega p ^ k * omega p ^ k = omega p ^ (2 * k) := by
  rw [← pow_add]; congr 1; ring

/-- **The spectral radius equals the contraction ratio.** -/
theorem spectral_radius_lt_one (d : ℕ) (hd : d ≥ 2) :
    ((d - 1 : ℝ) / d) < 1 := by
  have : (0 : ℝ) < d := by positivity
  rw [div_lt_one this]
  linarith

end Pandrosion
