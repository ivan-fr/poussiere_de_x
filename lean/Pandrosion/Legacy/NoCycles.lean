/-
  Universitas Pandrosion — Legacy / No periodic orbits.

  Under a linear contraction `|f s − r| ≤ c · |s − r|` with `c < 1`,
  the fixed point `r` is the *only* recurrent point: no `n`-cycle
  exists for any `n ≥ 1`. Combined with `contraction_general`
  (`Legacy.UniversalContraction`), this rules out every periodic
  Pandrosion orbit outside the fixed point.

  Ported from `NoCycles.lean` on the `main_previous` branch.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace Pandrosion

/-- **No 2-cycle under contraction.** If `f` contracts toward `r` by
    factor `c < 1`, then `f(f(s)) = s` forces `s = r`. -/
theorem no_two_cycle (f : ℝ → ℝ) (r : ℝ) (c : ℝ)
    (hc_nn : 0 ≤ c) (hc_lt : c < 1)
    (h_contract : ∀ s, |f s - r| ≤ c * |s - r|)
    (s : ℝ) (h_period : f (f s) = s) :
    s = r := by
  by_contra h_ne
  have h_pos : (0 : ℝ) < |s - r| := abs_pos.mpr (sub_ne_zero.mpr h_ne)
  have h1 : |f s - r| ≤ c * |s - r| := h_contract s
  have h2 : |f (f s) - r| ≤ c * |f s - r| := h_contract (f s)
  rw [h_period] at h2
  have h3 : |s - r| ≤ c ^ 2 * |s - r| := by
    calc |s - r|
      _ ≤ c * |f s - r| := h2
      _ ≤ c * (c * |s - r|) := mul_le_mul_of_nonneg_left h1 hc_nn
      _ = c ^ 2 * |s - r| := by ring
  have hc2 : c ^ 2 < 1 := by
    nlinarith [sq_nonneg c, sq_nonneg (1 - c)]
  have h4 : c ^ 2 * |s - r| < |s - r| := by
    calc c ^ 2 * |s - r|
      _ < 1 * |s - r| := mul_lt_mul_of_pos_right hc2 h_pos
      _ = |s - r| := one_mul _
  linarith

/-- **Contraction iterates compose:** `|f^n(s) − r| ≤ cⁿ · |s − r|`. -/
theorem contraction_iterate (f : ℝ → ℝ) (r : ℝ) (c : ℝ)
    (hc_nn : 0 ≤ c)
    (h_contract : ∀ s, |f s - r| ≤ c * |s - r|)
    (s : ℝ) : ∀ n : ℕ, |f^[n] s - r| ≤ c ^ n * |s - r| := by
  intro n
  induction n with
  | zero => simp
  | succ k ih =>
    simp only [Function.iterate_succ']
    show |f (f^[k] s) - r| ≤ c ^ (k + 1) * |s - r|
    calc |f (f^[k] s) - r|
      _ ≤ c * |f^[k] s - r| := h_contract _
      _ ≤ c * (c ^ k * |s - r|) := mul_le_mul_of_nonneg_left ih hc_nn
      _ = c ^ (k + 1) * |s - r| := by ring

/-- **No periodic orbit under contraction.** If `f` contracts toward
    `r` by factor `c < 1`, and `f^[n] s = s` for some `n ≥ 1`, then
    `s = r`. The fixed point is the only recurrent point. -/
theorem no_periodic_orbit (f : ℝ → ℝ) (r : ℝ) (c : ℝ)
    (hc_nn : 0 ≤ c) (hc_lt : c < 1)
    (h_contract : ∀ s, |f s - r| ≤ c * |s - r|)
    (s : ℝ) (n : ℕ) (hn : n ≥ 1)
    (h_period : f^[n] s = s) :
    s = r := by
  by_contra h_ne
  have h_pos : (0 : ℝ) < |s - r| := abs_pos.mpr (sub_ne_zero.mpr h_ne)
  have h_iter := contraction_iterate f r c hc_nn h_contract s n
  rw [h_period] at h_iter
  have hcn : c ^ n < 1 := pow_lt_one hc_nn hc_lt (by omega)
  have h_strict : c ^ n * |s - r| < |s - r| := by
    calc c ^ n * |s - r|
      _ < 1 * |s - r| := mul_lt_mul_of_pos_right hcn h_pos
      _ = |s - r| := one_mul _
  linarith

end Pandrosion
