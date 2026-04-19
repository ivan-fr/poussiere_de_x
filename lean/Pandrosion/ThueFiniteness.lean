/-
  Universitas Pandrosion — Lean 4 Formalization
  THUE FINITENESS MODULE

  Proves that the Pandrosion orbit visits at most ONE solution
  of each cubic Thue equation A³ - XB³ = c, and provides an
  explicit upper bound on the number of orbital solutions with
  bounded norm.

  Key results:
  1. Norm injectivity: d_m = d_n implies m = n
  2. Linear lower bound: |d_n| ≥ |d_0| + n
  3. Bounded norm → bounded index: |d_n| ≤ M implies n ≤ M - |d_0|
  4. Consequence: for X=2 starting at (1,1), the Thue equation
     A³ - 2B³ = c has at most |c| solutions in the Pandrosion orbit.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Tactic
import Pandrosion.ThueBridge

namespace Pandrosion

/-! ## §500. Integer Gap Principle

For integers, strict inequality implies a gap of at least 1.
This upgrades the strict monotonicity of |d_n| to a linear lower bound.
-/

/-- **Integer gap: a < b implies a + 1 ≤ b.** -/
theorem int_strict_gap (a b : ℤ) (h : a < b) : a + 1 ≤ b :=
  Int.add_one_le_of_lt h

/-! ## §501. Linear Lower Bound on Norm Growth

Since |d_n| is a strictly increasing sequence of nonneg integers,
consecutive terms differ by at least 1, giving |d_n| ≥ |d_0| + n.
This is STRONGER than the geometric bound for small n, and gives
an explicit bound on the index n in terms of |d_n|.
-/

/-- **Linear lower bound on the norm sequence.**
    |d_n| ≥ |d_0| + n for all n. This follows from the fact that
    |d_n| is a strictly increasing sequence of nonneg integers. -/
theorem norm_linear_lower_bound (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ∀ n : ℕ, |d 0| + (↑n : ℤ) ≤ |d n| := by
  intro n
  induction n with
  | zero => simp
  | succ k ih =>
    have hstrict := norm_values_distinct d Φ hd0 hΦ hd k (k + 1) (by omega)
    have hgap : |d k| + 1 ≤ |d (k + 1)| := int_strict_gap _ _ hstrict
    push_cast
    push_cast at ih
    linarith

/-! ## §502. Thue Value Injectivity

Each integer value c appears at most once in the norm sequence.
This means each Thue equation A³ - XB³ = c has at most ONE
solution in any given Pandrosion orbit.
-/

/-- **Thue orbit injectivity.**
    If d_m = d_n in the Pandrosion norm sequence, then m = n.
    Equivalently: each Thue equation A³ - XB³ = c has AT MOST ONE
    solution (A_n, B_n) in the orbit. -/
theorem thue_value_injective (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (m n : ℕ) (heq : d m = d n) : m = n := by
  by_contra hmn
  rcases Nat.lt_or_gt_of_ne hmn with h | h
  · have := norm_values_distinct d Φ hd0 hΦ hd m n h
    rw [heq] at this
    exact lt_irrefl _ this
  · have := norm_values_distinct d Φ hd0 hΦ hd n m h
    rw [heq] at this
    exact lt_irrefl _ this

/-- **Stronger: even |d_m| = |d_n| implies m = n.**
    Different Thue equations A³ - XB³ = c and A³ - XB³ = -c
    also have disjoint orbital solutions. -/
theorem thue_abs_value_injective (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (m n : ℕ) (heq : |d m| = |d n|) : m = n := by
  by_contra hmn
  rcases Nat.lt_or_gt_of_ne hmn with h | h
  · exact absurd heq (ne_of_lt (norm_values_distinct d Φ hd0 hΦ hd m n h))
  · exact absurd heq.symm (ne_of_lt (norm_values_distinct d Φ hd0 hΦ hd n m h))

/-! ## §503. Bounded Norm Implies Bounded Index

The linear lower bound |d_n| ≥ |d_0| + n directly gives:
if |d_n| ≤ M, then n ≤ M - |d_0|.

This is an EFFECTIVE finiteness result: it gives an explicit
computable upper bound on which orbital indices can satisfy
a given Thue equation.
-/

/-- **Bounded norm ⟹ bounded index.**
    If |d_n| ≤ M, then n ≤ M - |d_0|.
    This is the effective finiteness certificate. -/
theorem thue_orbital_bound (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (M : ℤ) (n : ℕ) (hn : |d n| ≤ M) :
    (n : ℤ) ≤ M - |d 0| := by
  have hbound := norm_linear_lower_bound d Φ hd0 hΦ hd n
  linarith

/-- **Escape theorem: eventually all norms exceed any bound.**
    For n > M - |d_0|, we have |d_n| > M.
    This is the contrapositive of the bounded-index theorem. -/
theorem thue_orbital_escape (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (M : ℤ) (n : ℕ) (hn : M - |d 0| < (n : ℤ)) :
    M < |d n| := by
  have hbound := norm_linear_lower_bound d Φ hd0 hΦ hd n
  linarith

/-! ## §504. Concrete Certificate for X = 2

For the Pandrosion orbit targeting ∛2, starting from (A₀, B₀) = (1, 1):
  d₀ = 1³ - 2·1³ = -1, so |d₀| = 1.

The bound gives: |d_n| ≤ M implies n ≤ M - 1.

Therefore: the Thue equation A³ - 2B³ = c with |c| ≤ M has at most M
solutions in the Pandrosion orbit starting from (1, 1).

In particular, A³ - 2B³ = 1 has at most ONE orbital solution (n = 0),
and A³ - 2B³ = 43 has at most ONE orbital solution (n = 1).
-/

/-- **For X=2: the initial norm has absolute value 1.** -/
theorem abs_initial_norm_x2 : |(1 : ℤ) ^ 3 - 2 * 1 ^ 3| = 1 := by norm_num

/-- **For X=2: the equation A³ - 2B³ = c with |c| ≤ M
    has at most M solutions in the Pandrosion orbit from (1,1).**
    Specifically: if d_n = A_n³ - 2·B_n³ satisfies |d_n| ≤ M,
    then n ≤ M - 1, which means n ∈ {0, 1, ..., M-1}. -/
theorem thue_x2_orbit_bound (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0_ne : d 0 ≠ 0)
    (hd0_val : |d 0| = 1)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (M : ℤ) (n : ℕ) (hn : |d n| ≤ M) :
    (n : ℤ) ≤ M - 1 := by
  have hbound := thue_orbital_bound d Φ hd0_ne hΦ hd M n hn
  rw [hd0_val] at hbound
  exact hbound

end Pandrosion
