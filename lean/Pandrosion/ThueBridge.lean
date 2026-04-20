/-
  Universitas Pandrosion — Lean 4 Formalization
  THUE BRIDGE MODULE: Norm Amplification and Diophantine Consequences

  This module proves the exact multiplicative law governing the cubic norm
  d_n = A_n³ - X·B_n³ under the Pandrosion iteration.

  The central identity:
    d_{n+1} = d_n · Φ(A_n, B_n, X)
  where Φ = A⁹ - 14·X·A⁶·B³ - 20·X²·A³·B⁶ + 8·X³·B⁹

  is proven purely by `ring`. Combined with the effective irrationality
  bound (EffectiveIrrationality.lean) and coprimality isolation
  (CoprimalityIsolation.lean), this provides the algebraic scaffolding
  connecting the Pandrosion iteration to the theory of Thue equations.

  Key consequences proven here:
  1. The norm amplification identity (ring certificate)
  2. The cross-determinant of consecutive approximants (ring certificate)
  3. Geometric growth of |d_n| under Φ-amplification (induction)
  4. Injectivity of the norm sequence (distinct Thue values)
-/
import Mathlib.Data.Int.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §400. The Norm Amplification Identity

The Pandrosion iteration maps (A, B) ↦ (A', B') where
  A' = A(A³ + 4XB³),    B' = B(3A³ + 2XB³).

The cubic norm d = A³ - XB³ transforms as d' = d · Φ(A,B,X)
where Φ is a degree-9 polynomial. This is the engine of the
Diophantine connection.
-/

/-- **The Pandrosion amplification factor Φ.**
    Φ(A, B, X) = A⁹ - 14·X·A⁶·B³ - 20·X²·A³·B⁶ + 8·X³·B⁹.
    This is the exact multiplicative factor governing norm evolution. -/
def pandrosion_phi (A B X : ℤ) : ℤ :=
  A ^ 9 - 14 * X * A ^ 6 * B ^ 3 - 20 * X ^ 2 * A ^ 3 * B ^ 6 + 8 * X ^ 3 * B ^ 9

/-- **Norm Amplification Identity.**
    After one Pandrosion step, the cubic norm transforms multiplicatively:
      A'³ - X·B'³ = (A³ - X·B³) · Φ(A, B, X).

    This is a pure polynomial identity, verified by `ring`.
    It is the algebraic engine connecting Pandrosion to Thue equations. -/
theorem norm_amplification (A B X : ℤ) :
    let A' := A * (A ^ 3 + 4 * X * B ^ 3)
    let B' := B * (3 * A ^ 3 + 2 * X * B ^ 3)
    A' ^ 3 - X * B' ^ 3 =
      (A ^ 3 - X * B ^ 3) * pandrosion_phi A B X := by
  intros
  unfold pandrosion_phi
  ring

/-! ## §401. Cross-Determinant of Consecutive Approximants

The determinant A_n·B_{n+1} - A_{n+1}·B_n encodes the arithmetic
separation between consecutive Pandrosion approximants. Its exact
factored form connects to the classical Thue equation theory.
-/

/-- **Cross-determinant identity.**
    A·B' - B·A' = 2·A·B·(A³ - X·B³) = 2·A·B·d.

    This shows that consecutive Pandrosion approximants
    have a determinant controlled by the norm d. Since d ≠ 0
    (by norm propagation), consecutive approximants are always
    distinct fractions. -/
theorem cross_determinant (A B X : ℤ) :
    let A' := A * (A ^ 3 + 4 * X * B ^ 3)
    let B' := B * (3 * A ^ 3 + 2 * X * B ^ 3)
    A * B' - B * A' = 2 * A * B * (A ^ 3 - X * B ^ 3) := by
  intros
  ring

/-! ## §402. Φ Factor Evaluations

Computing Φ at specific algebraic points connects the abstract
amplification to concrete Thue equation dynamics.
-/

/-- **Φ at the initial point for ∛2.**
    Φ(1, 1, 2) = 1 - 28 - 80 + 64 = -43. -/
theorem phi_at_x2_initial : pandrosion_phi 1 1 2 = -43 := by
  unfold pandrosion_phi; norm_num

/-- **The absolute value of Φ(1,1,2) is large, certifying amplification.** -/
theorem phi_at_x2_initial_abs : |pandrosion_phi 1 1 2| ≥ 2 := by
  unfold pandrosion_phi; norm_num

/-- **Φ at the initial point for ∛3.**
    Φ(1, 1, 3) = 1 - 42 - 180 + 216 = -5. -/
theorem phi_at_x3_initial : pandrosion_phi 1 1 3 = -5 := by
  unfold pandrosion_phi; norm_num

/-- **|Φ(1,1,3)| ≥ 2.** -/
theorem phi_at_x3_initial_abs : |pandrosion_phi 1 1 3| ≥ 2 := by
  unfold pandrosion_phi; norm_num

/-! ## §402b. Concrete Pandrosion Iterations

Explicit integer computations that ground the abstract theory.
-/

/-- **The first Pandrosion iterate for X=2, starting from (1,1).**
    A₁ = 1·(1³ + 4·2·1³) = 9,   B₁ = 1·(3·1³ + 2·2·1³) = 7. -/
theorem pandrosion_step_x2 :
    1 * (1 ^ 3 + 4 * 2 * 1 ^ 3) = (9 : ℤ) ∧
    1 * (3 * 1 ^ 3 + 2 * 2 * 1 ^ 3) = (7 : ℤ) := by
  constructor <;> norm_num

/-- **Verification: d₁ = d₀ · Φ(1,1,2) = (-1)·(-43) = 43.** -/
theorem norm_amplification_x2_verify :
    (-1 : ℤ) * pandrosion_phi 1 1 2 = 43 := by
  unfold pandrosion_phi; norm_num

/-! ## §403. Geometric Growth of the Norm

If the amplification factor |Φ_n| ≥ 2 at every step, the norm
|d_n| grows at least geometrically: |d_n| ≥ |d₀| · 2ⁿ.

This is the key number-theoretic consequence: the Pandrosion
orbit visits EXPONENTIALLY MANY distinct Thue equation values.
Since each d_n is a nonzero integer, they are well-separated,
giving the iteration genuine Diophantine content.
-/

/-- **Geometric norm growth under amplification.**
    If a sequence d satisfies d(n+1) = d(n) · Φ(n)
    with |Φ(n)| ≥ 2 for all n, then |d(n)| ≥ |d(0)| · 2ⁿ. -/
theorem norm_geometric_growth (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ∀ n, |d 0| * (2 : ℤ) ^ n ≤ |d n| := by
  intro n
  induction n with
  | zero => simp
  | succ k ih =>
    rw [hd k, abs_mul]
    calc |d 0| * 2 ^ (k + 1)
        = |d 0| * 2 ^ k * 2 := by ring
      _ ≤ |d k| * 2 := by linarith
      _ ≤ |d k| * |Φ k| := by
          exact Int.mul_le_mul_of_nonneg_left (hΦ k) (abs_nonneg _)

/-- **Corollary: the norm never vanishes under amplification.**
    If d₀ ≠ 0 and |Φ_n| ≥ 2 for all n, then d_n ≠ 0 for all n. -/
theorem norm_never_zero (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ∀ n, d n ≠ 0 := by
  intro n
  have hgrowth := norm_geometric_growth d Φ hΦ hd n
  have hd0_pos : 0 < |d 0| := abs_pos.mpr hd0
  have hpow_pos : (0 : ℤ) < 2 ^ n := by positivity
  have h_lower : 0 < |d n| := by linarith [mul_pos hd0_pos hpow_pos]
  exact abs_pos.mp h_lower

/-! ## §404. Injectivity of the Norm Sequence

Under amplification, the sequence |d_n| is strictly increasing
(after step 0), hence the d_n are all distinct. This means the
Pandrosion orbit produces infinitely many distinct solutions to
distinct Thue equations A³ - XB³ = d_n.
-/

/-- **Strict norm increase.**
    If |Φ_n| ≥ 2 and d_n ≠ 0, then |d_{n+1}| > |d_n|. -/
theorem norm_strict_increase (d_n : ℤ) (Φ_n : ℤ)
    (hd : d_n ≠ 0) (hΦ : 2 ≤ |Φ_n|) :
    |d_n| < |d_n * Φ_n| := by
  rw [abs_mul]
  have hd_pos : 0 < |d_n| := abs_pos.mpr hd
  calc |d_n| = |d_n| * 1 := by ring
    _ < |d_n| * 2 := by linarith
    _ ≤ |d_n| * |Φ_n| := by
        exact Int.mul_le_mul_of_nonneg_left hΦ (le_of_lt hd_pos)

/-- **The norm values are pairwise distinct.**
    If |Φ_n| ≥ 2 for all n and d₀ ≠ 0, then m ≠ n implies |d_m| ≠ |d_n|. -/
theorem norm_values_distinct (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ∀ m n, m < n → |d m| < |d n| := by
  intro m n hmn
  have h_never_zero := norm_never_zero d Φ hd0 hΦ hd
  induction n with
  | zero => exact absurd hmn (Nat.not_lt_zero m)
  | succ k ih =>
    by_cases hmk : m < k
    · have hk := ih hmk
      have hstep : |d k| < |d (k + 1)| := by
        rw [hd k]
        exact norm_strict_increase (d k) (Φ k) (h_never_zero k) (hΦ k)
      exact lt_trans hk hstep
    · have hmk_eq : m = k := by omega
      rw [hmk_eq, hd k]
      exact norm_strict_increase (d k) (Φ k) (h_never_zero k) (hΦ k)

/-! ## §405. The Effective Irrationality Strengthening

The standard Liouville bound gives |p/q - ∛X| ≥ 1/(3r²q³), i.e., μ = 3.
The norm amplification allows a STRENGTHENED bound for the specific
Pandrosion approximants: since |d_n| grows geometrically while the
quadratic form grows polynomially in B_n, the effective bound improves.

For a GENERAL rational p/q (not from the Pandrosion sequence),
improving μ below 3 requires the Thue-Siegel auxiliary polynomial
method or Baker's theory of linear forms in logarithms, which are
beyond the current Lean formalization scope. We formalize the bridge
identities that would feed into such an argument.
-/

/-- **Strengthened bound for Pandrosion approximants.**
    If |d| ≥ K (bigger than 1, from amplification), we get a better
    effective bound: |A - rB| ≥ K / (A² + A·rB + (rB)²). -/
theorem effective_distance_bound_amplified (A B r : ℝ) (d : ℤ) (K : ℤ)
    (hA : 0 < A) (hrB : 0 < r * B)
    (h_eq : (d : ℝ) = A ^ 3 - (r * B) ^ 3)
    (hK : K ≤ |d|) :
    |A - r * B| ≥ (K : ℝ) / (A ^ 2 + A * (r * B) + (r * B) ^ 2) := by
  have h_form_pos : 0 < A ^ 2 + A * (r * B) + (r * B) ^ 2 := by
    nlinarith [sq_nonneg A, sq_nonneg (r * B), mul_pos hA hrB]
  have h_factor : (d : ℝ) = (A - r * B) * (A ^ 2 + A * (r * B) + (r * B) ^ 2) := by
    rw [h_eq]; ring
  have h_abs : |(d : ℝ)| = |A - r * B| * (A ^ 2 + A * (r * B) + (r * B) ^ 2) := by
    rw [h_factor, abs_mul, abs_of_pos h_form_pos]
  have h_K_le : (K : ℝ) ≤ |(d : ℝ)| := by
    rw [← Int.cast_abs]
    exact_mod_cast hK
  have h_product : (K : ℝ) ≤ |A - r * B| * (A ^ 2 + A * (r * B) + (r * B) ^ 2) := by
    linarith
  rw [ge_iff_le, div_le_iff h_form_pos]
  exact h_product

end Pandrosion
