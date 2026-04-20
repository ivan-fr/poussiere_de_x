/-
  Universitas Pandrosion — Lean 4 Formalization
  ZSYGMONDY ORBITAL MODULE

  Zsygmondy's theorem (1892) states that for coprime integers a > b > 0,
  the sequence a^n - b^n has a primitive prime divisor for all n ≥ 2
  (except for a few specific base cases).
  
  For the Pandrosion norm orbit d_n = d_0 * ∏ Φ_k, the geometric escape
  |d_n| ≥ 2^n ensures that new magnitude is continuously generated.
  While proving the strict emergence of a *prime* requires deep prime
  distribution laws not yet in Mathlib 4, we formalize the ALGEBRAIC
  core: the multiplier Φ_n must contribute new absolute size (≥ 2) at
  every step, so |d_{n+1}| strictly exceeds |d_n|.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Tactic
import Pandrosion.DynMordellLang

namespace Pandrosion

/-! ## §3300. Zsygmondy Magnitude Increment

The absolute value of the norm sequence is strictly increasing,
guaranteeing that new prime factors (or at least strictly higher
multiplicities) MUST appear. This is the orbital Zsygmondy step.
-/

/-- **Non-zero orbit states.**
    If d 0 ≠ 0 and multipliers are non-zero, then all d n ≠ 0. -/
theorem zsygmondy_orbit_nonzero
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ∀ n, d n ≠ 0 := by
  intro n
  induction n with
  | zero => exact hd0
  | succ k ih =>
    rw [hd k]
    have h_phi_nz : Φ k ≠ 0 := by
      intro h_zero
      have h1 : 2 ≤ |Φ k| := hΦ k
      rw [h_zero] at h1
      norm_num at h1
    exact mul_ne_zero ih h_phi_nz

/-- **Zsygmondy orbital increment.**
    At each step, the absolute norm is strictly amplified. -/
theorem zsygmondy_orbital_increment
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (n : ℕ) :
    |d n| < |d (n + 1)| := by
  have hd_succ : d (n + 1) = d n * Φ n := hd n
  have h_ne_zero : d n ≠ 0 := zsygmondy_orbit_nonzero d Φ hd0 hΦ hd n
  have h_abs_pos : 0 < |d n| := abs_pos.mpr h_ne_zero
  have h_abs : |d (n + 1)| = |d n| * |Φ n| := by rw [hd_succ, abs_mul]
  have hΦ_lower : 2 ≤ |Φ n| := hΦ n
  calc |d n|
    < |d n| * 2 := by linarith
    _ ≤ |d n| * |Φ n| := mul_le_mul_of_nonneg_left hΦ_lower (by linarith)
    _ = |d (n + 1)| := h_abs.symm

end Pandrosion
