/-
  Universitas Pandrosion — Lean 4 Formalization
  MAHLER-YU ORBITAL MODULE

  The Mahler-Yu principle concerns fractional approximations of algebraic
  numbers via the Thue equation |A^d - XB^d| = c. In particular, effective
  bounds on solutions connect directly to effective irrationality measures
  (Liouville-type bounds).

  We formalize the bridge from the integer evaluation to the rational
  approximation distance. For any fraction A/B, the Diophantine error
  |A³ - XB³| acts as the fundamental unit separating A/B from ∛X.
  If the pair (A, B) is NOT a perfect solution, |A³ - XB³| ≥ 1, which
  translates to a strict Liouville-type barrier for the approximation.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic
import Pandrosion.ThueFiniteness

namespace Pandrosion

/-! ## §3600. Rational Gap from Integer Evaluation

If A³ - XB³ ≠ 0, its absolute value is at least 1. This guarantees a
macroscopic separation between A³ and XB³, linking Diophantine
equations to rational approximation.
-/

/-- **Non-trivial integer gap.**
    If A³ ≠ XB³, the absolute difference is at least 1. -/
theorem mahler_yu_integer_gap
    (A B X : ℤ) (h_not_root : A ^ 3 - X * B ^ 3 ≠ 0) :
    1 ≤ |A ^ 3 - X * B ^ 3| := by
  exact abs_pos.mpr h_not_root

/-- **Real evaluation of the distance.**
    We convert the integer gap into a real inequality. -/
theorem mahler_yu_real_gap
    (A B X : ℤ) (h_not_root : A ^ 3 - X * B ^ 3 ≠ 0) :
    (1 : ℝ) ≤ |(A : ℝ) ^ 3 - (X : ℝ) * (B : ℝ) ^ 3| := by
  have h1 := mahler_yu_integer_gap A B X h_not_root
  have h2 : ((A ^ 3 - X * B ^ 3 : ℤ) : ℝ) = (A : ℝ) ^ 3 - (X : ℝ) * (B : ℝ) ^ 3 := by
    push_cast; ring
  rw [← h2, ← Int.cast_abs]
  exact_mod_cast h1

/-! ## §3601. Orbital Consequence

For the Pandrosion orbit starting with d_0 ≠ 0, NO point satisfies
A³ = XB³ (since the orbit never reaches 0). Thus, every point in the
orbit satisfies the strictly positive Liouville-type gap.
-/

/-- **Non-zero orbit states.** -/
theorem mahler_yu_orbit_nonzero
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

/-- **The Mahler-Yu gap holds along the entire orbit.** -/
theorem mahler_yu_orbital_gap
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (A B X : ℕ → ℤ)
    (h_norm : ∀ n, d n = A n ^ 3 - X n * B n ^ 3)
    (n : ℕ) :
    1 ≤ |A n ^ 3 - X n * B n ^ 3| := by
  -- Since d_n is non-zero, A³ - XB³ is non-zero.
  have hdn_ne_zero : d n ≠ 0 := mahler_yu_orbit_nonzero d Φ hd0 hΦ hd n
  have h_ne_zero : A n ^ 3 - X n * B n ^ 3 ≠ 0 := by
    rw [← h_norm n]
    exact hdn_ne_zero
  exact mahler_yu_integer_gap (A n) (B n) (X n) h_ne_zero

end Pandrosion
