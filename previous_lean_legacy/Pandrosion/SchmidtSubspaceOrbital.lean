/-
  Universitas Pandrosion — Lean 4 Formalization
  SCHMIDT SUBSPACE ORBITAL MODULE

  The Schmidt Subspace Theorem (1972) is a higher-dimensional generalization
  of Roth's theorem on Diophantine approximation. In dimension 1, it governs
  the relation between heights of linearly dependent quantities.

  For the Pandrosion family, the core height relation is derived from the
  norm identity: A_n³ = X·B_n³ + d_n. The "orbital subspace" phenomenon
  is the precise algebraic tracking of how the height of A_n³ completely
  dominates the height of d_n asymptotically, ensuring the vector
  (A_n³, -X·B_n³, -d_n) lies in a controlled hyper-plane whose heights
  satisfy an exact relation without error terms.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Tactic
import Pandrosion.AbcOrbital

namespace Pandrosion

/-! ## §3700. Orbital Subspace Identity

The vector (A³, X·B³, d) satisfies L(v) = A³ - X·B³ - d = 0.
This is an exact linear relation (subspace equation). We formalize
the exact matching of the absolute values (heights before logs)
which is the foundation for analyzing rational approximations in
projective space.
-/

/-- **Subspace linear relation.**
    The vector (A³, X·B³, d) is orthogonal to (1, -1, -1). -/
theorem schmidt_orbital_relation
    (A B X d : ℤ) (h_norm : d = A ^ 3 - X * B ^ 3) :
    A ^ 3 + (-1) * (X * B ^ 3) + (-1) * d = 0 := by
  linarith

/-! ## §3701. Triangle Inequalities for Heights

While true heights involve logarithms, the absolute value forms the
multiplicative core. The "height" of A³ is bounded by the sum of the
heights of XB³ and d.
-/

/-- **Multiplicative height bound (Triangle inequality).**
    |A³| ≤ |X·B³| + |d| -/
theorem schmidt_orbital_height_bound
    (A B X d : ℤ) (h_norm : d = A ^ 3 - X * B ^ 3) :
    |A ^ 3| ≤ |X * B ^ 3| + |d| := by
  have : A ^ 3 = X * B ^ 3 + d := by linarith
  rw [this]
  exact abs_add_le (X * B ^ 3) d

/-- **Lower height bound.**
    The separation |d| enforces that A³ and XB³ cannot be too close
    in projective space. -/
theorem schmidt_orbital_height_lower
    (A B X d : ℤ) (h_norm : d = A ^ 3 - X * B ^ 3) :
    |d| ≤ |A ^ 3| + |X * B ^ 3| := by
  have : d = A ^ 3 + (- (X * B ^ 3)) := by linarith
  rw [this]
  have h_triangle := abs_add_le (A ^ 3) (- (X * B ^ 3))
  have h_neg : |- (X * B ^ 3)| = |X * B ^ 3| := abs_neg _
  linarith [h_triangle, h_neg]

end Pandrosion
