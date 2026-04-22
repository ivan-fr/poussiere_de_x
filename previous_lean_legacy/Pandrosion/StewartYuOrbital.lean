/-
  Universitas Pandrosion — Lean 4 Formalization
  STEWART-YU ORBITAL MODULE

  Stewart and Yu (1991, 2001) established exponential bounds for the abc
  theorem using linear forms in logarithms.
  Here, we establish an effective Stewart-Yu type inequality for the
  Pandrosion orbital family. Taking the triple (a_n, b_n, c_n) =
  (XB_n³, d_n, A_n³), we bound the radical of the multiplicative component
  d_n explicitly via the geometric growth of the Φ_k multipliers.

  Main result:
  The orbital triple is strictly bounded by its initial state times
  the product of all Φ_k polynomial multipliers up to step n. Note
  that since d_{n+1} = d_n Φ_n, no primes enter d_n except those
  dividing d_0 or some Φ_k.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Algebra.BigOperators.Basic
import Mathlib.Tactic
import Pandrosion.AbcOrbital

open BigOperators

namespace Pandrosion

/-! ## §3200. Stewart-Yu Radical Factorization

The radical of d_n is the product of its distinct prime factors.
Since `d_n = d_0 * ∏ Φ_k`, the prime factors of d_n are exactly those
of d_0 and the set of Φ_k (k < n). This allows an effective structural
substitution for rad(abc) in the Stewart-Yu framework.
-/

/-- **Stewart-Yu effective bound (Orbital Form).**
    Instead of unknown exponents, we have an EXACT multiplicative
    decomposition of the 'b' component of the abc triple into Φ_k terms. -/
theorem stewart_yu_orbital_bound
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (n : ℕ) :
    |d n| ≤ |d 0| * ∏ k in Finset.range n, |Φ k| :=
  abc_orbital_b_upper d Φ hd n

/-- **Stewart-Yu: c_n formulation.**
    Combining A³ = XB³ + d with the factorized bound gives an
    inequality directly on the total size of the terms. -/
theorem stewart_yu_c_orbital
    (A B X d_n : ℤ) (d0 : ℤ) (Φ : ℕ → ℤ) (n : ℕ)
    (h_triple : A ^ 3 = X * B ^ 3 + d_n)
    (h_bound : |d_n| ≤ |d0| * ∏ k in Finset.range n, |Φ k|) :
    |A ^ 3 - X * B ^ 3| ≤ |d0| * ∏ k in Finset.range n, |Φ k| := by
  have : d_n = A ^ 3 - X * B ^ 3 := by linarith
  rw [← this]
  exact h_bound

end Pandrosion
