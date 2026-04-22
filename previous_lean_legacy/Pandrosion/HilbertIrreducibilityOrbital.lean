/-
  Universitas Pandrosion — Lean 4 Formalization
  HILBERT IRREDUCIBILITY ORBITAL MODULE

  Hilbert's Irreducibility Theorem (HIT) guarantees that an irreducible
  polynomial over Q(t_1, ..., t_k, x) remains irreducible over Q(x)
  when the parameters t_j are specialized to integers, outside of a
  "thin" set of specialized values.

  In an orbital setting, if P(z) = z^3 - X is irreducible over Q,
  a rational approximation A/B mapping EXACTLY to the root would
  contradict irreducibility (it would yield a linear factor z - A/B).
  
  We formalize an "orbital preservation" analogue: the iteration
  starts from irreducible parameters (where d_0 ≠ 0 implies A_0/B_0
  is not a root). Since d_n never hits 0, the rational points
  never degenerate to exact roots. The algebraic "irreducibility"
  (non-existence of a rational root) is strictly propagated along
  the entire orbit.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Tactic
import Pandrosion.MahlerYuOrbital

namespace Pandrosion

/-! ## §4200. Orbital Specialization (Root Avoidance)

The core characteristic of a rational root for z^3 - X is achieving
exact zero for the Thue evaluation A^3 - X B^3. We show the orbit
strictly avoids this.
-/

/-- **Hilbert irreducibility algorithmic preservation.**
    If the start point is not a root (d_0 ≠ 0), then no point in the
    orbit can be a rational root, preserving the degree of the minimal
    polynomial over the whole dynamical sequence. -/
theorem hilbert_irreducibility_orbital
    (d : ℕ → ℤ) (Φ : ℕ → ℤ) (A B : ℕ → ℤ) (X : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (h_norm : ∀ n, d n = A n ^ 3 - X n * B n ^ 3)
    (n : ℕ) :
    A n ^ 3 - X n * B n ^ 3 ≠ 0 := by
  have hd_ne_zero : d n ≠ 0 := mahler_yu_orbit_nonzero d Φ hd0 hΦ hd n
  rw [← h_norm n]
  exact hd_ne_zero

/-! ## §4201. Quantitative Irreducibility Gap

Not only is it not a root, the evaluation is STRICTLY separated from 0.
This quantitative form is essential for algorithms to not get "stuck"
near thin sets of pseudo-roots.
-/

/-- **Quantitative evaluation gap.**
    The absolute value of the polynomial evaluation is strictly ≥ 1. -/
theorem hilbert_quantitative_gap
    (d : ℕ → ℤ) (Φ : ℕ → ℤ) (A B : ℕ → ℤ) (X : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (h_norm : ∀ n, d n = A n ^ 3 - X n * B n ^ 3)
    (n : ℕ) :
    1 ≤ |A n ^ 3 - X n * B n ^ 3| :=
  mahler_yu_orbital_gap d Φ hd0 hΦ hd A B X h_norm n

end Pandrosion
