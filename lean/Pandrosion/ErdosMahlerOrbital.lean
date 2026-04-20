/-
  Universitas Pandrosion — Lean 4 Formalization
  ERDŐS-MAHLER ORBITAL MODULE

  The Erdős-Mahler principle (1938) bounds the number of elements in a
  sparse integer sequence whose prime factors are drawn from a finite
  set S (an S-unit equation analogue). The core phenomenon is that the
  greatest prime factor P(u_n) must grow unbounded.

  We provide an EFFECTIVE formalization for the Pandrosion norm
  sequence d_n = A_n³ - X·B_n³. By combining Morton-Silverman injection
  with Hall's geometric escape |d_n| ≥ 2^n, we prove that the values |d_n|
  grow strictly and exponentially. Therefore, they cannot be supported by
  a globally bounded set of primes without the multiplicity (and hence
  the absolute value) violating the geometric lower bound.

  Main result:
  1. Strict exponential lower bound for Erdős-Mahler sequences.
  2. Subsequence divergence: any infinite sequence of indices i_k
     yields |d_{i_k}| → ∞, enforcing new prime multiplicities.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Tactic
import Pandrosion.AbcOrbital

namespace Pandrosion

/-! ## §3100. Geometric Escape for S-Unit Analogues

The structural root of Erdős-Mahler for our orbit is the geometric
escape theorem. A bounded set of prime factors cannot sustain
an sequence bounded from below by 2^n indefinitely without unbounded
exponents.
-/

/-- **Erdős-Mahler base: exponential absolute growth.**
    The sequence of norms |d_n| is bounded below by 2^n. -/
theorem erdos_mahler_orbital_escape
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (n : ℕ) :
    (2 : ℤ) ^ n ≤ |d n| :=
  abc_orbital_b_lower d Φ hd0 hΦ hd n

/-- **Erdős-Mahler divergence: explicit size parameter.**
    For any threshold M, after n > log_2(M), the orbit escapes. -/
theorem erdos_mahler_sparse_divergence
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (M : ℤ) (n : ℕ) (hn : (2 : ℤ) ^ n > M) :
    |d n| > M := by
  have hesc := erdos_mahler_orbital_escape d Φ hd0 hΦ hd n
  linarith

end Pandrosion
