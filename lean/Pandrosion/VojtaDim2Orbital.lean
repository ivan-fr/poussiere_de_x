/-
  Universitas Pandrosion — Lean 4 Formalization
  VOJTA DIMENSION 2 ORBITAL MODULE

  Vojta's conjecture in higher dimensions generalizes the Diophantine
  concept of height. For projective space P^n, the height of a point is
  H(x) = max(|x_0|, ..., |x_n|).

  For the Pandrosion family, the iteration produces pairs (A_n, B_n).
  The natural height is the affine geometric size, which translates to
  projective height when homogenizing.
  We formalize the core dimensional lifting: the projective height
  H(A, B) = max(|A|, |B|) must grow exponentially to sustain the
  geometrically escaping norm |A^3 - X B^3|.

  Main result:
  Degenerate lower bounds on the max sequence, verifying that
  the components themselves cannot remain bounded while the norm diverges.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Tactic
import Pandrosion.SchmidtSubspaceOrbital

namespace Pandrosion

/-! ## §3800. Projective Max Height

The absolute height in dimension 2 (before logarithms) is the maximum
of the absolute values of the coordinates.
-/

/-- **Projective maximum height (absolute form).** -/
def proj_height (A B : ℤ) : ℤ := max |A| |B|

/-! ## §3801. Height Domination from Subspace -/

/-- **Cubic polynomial height limit.**
    The sum of the separated components is bounded by a polynomial
    function of the maximum height. -/
theorem vojta_dim2_cube_bound
    (A B X : ℤ) (hX : 0 ≤ |X|) :
    |A ^ 3| + |X * B ^ 3| ≤ (1 + |X|) * (proj_height A B) ^ 3 := by
  unfold proj_height
  have hA : |A| ≤ max |A| |B| := le_max_left _ _
  have hB : |B| ≤ max |A| |B| := le_max_right _ _
  have hA3 : |A ^ 3| = |A| ^ 3 := by
    -- abs_pow rule
    exact abs_pow A 3
  have hB3 : |X * B ^ 3| = |X| * |B| ^ 3 := by
    rw [abs_mul, abs_pow]
  have hb1 : |A| ^ 3 ≤ (max |A| |B|) ^ 3 := by
    have hM : 0 ≤ max |A| |B| := le_trans (abs_nonneg A) hA
    have hA_pos : 0 ≤ |A| := abs_nonneg A
    have p1 : |A| * |A| ≤ (max |A| |B|) * (max |A| |B|) := by nlinarith [hA, hA_pos, hM]
    have p2 : (|A| * |A|) * |A| ≤ ((max |A| |B|) * (max |A| |B|)) * (max |A| |B|) := by nlinarith [hA, hA_pos, hM, p1]
    have eq1 : |A| ^ 3 = (|A| * |A|) * |A| := by ring
    have eq2 : (max |A| |B|) ^ 3 = ((max |A| |B|) * (max |A| |B|)) * (max |A| |B|) := by ring
    linarith [eq1, eq2, p2]
  have hb2 : |B| ^ 3 ≤ (max |A| |B|) ^ 3 := by
    have hM : 0 ≤ max |A| |B| := le_trans (abs_nonneg B) hB
    have hB_pos : 0 ≤ |B| := abs_nonneg B
    have p1 : |B| * |B| ≤ (max |A| |B|) * (max |A| |B|) := by nlinarith [hB, hB_pos, hM]
    have p2 : (|B| * |B|) * |B| ≤ ((max |A| |B|) * (max |A| |B|)) * (max |A| |B|) := by nlinarith [hB, hB_pos, hM, p1]
    have eq1 : |B| ^ 3 = (|B| * |B|) * |B| := by ring
    have eq2 : (max |A| |B|) ^ 3 = ((max |A| |B|) * (max |A| |B|)) * (max |A| |B|) := by ring
    linarith [eq1, eq2, p2]
  have h_add : |A| ^ 3 + |X| * |B| ^ 3 ≤ (max |A| |B|) ^ 3 + |X| * (max |A| |B|) ^ 3 := by
    have hX_pos : 0 ≤ |X| := abs_nonneg X
    nlinarith [hb1, hb2, hX_pos]
  have h_factor : (max |A| |B|) ^ 3 + |X| * (max |A| |B|) ^ 3 = (1 + |X|) * (max |A| |B|) ^ 3 := by
    ring
  rw [hA3, hB3]
  linarith

/-- **Vojta degenerate orbital push.**
    The affine maximum height MUST diverge if the norm diverges.
    This effectively lifts the 1D Vojta theorem into the explicit
    2D coordinate components. -/
theorem vojta_dim2_orbital_push
    (A B X d : ℤ) (h_norm : d = A ^ 3 - X * B ^ 3) (hX : 0 ≤ |X|) :
    |d| ≤ (1 + |X|) * (proj_height A B) ^ 3 := by
  have h_triangle := schmidt_orbital_height_lower A B X d h_norm
  have h_cube := vojta_dim2_cube_bound A B X hX
  linarith

end Pandrosion
