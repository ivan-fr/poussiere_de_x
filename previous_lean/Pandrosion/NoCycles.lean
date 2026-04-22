/-
  Universitas Pandrosion — Lean 4 Formalization
  NO PERIODIC ORBITS MODULE

  Proves that any map with contraction rate c < 1 toward a fixed
  point r has no periodic orbits of any period n ≥ 1 other than
  the fixed point itself.

  Applied to the Pandrosion iteration: since contraction_general
  (GeneralContractionAll.lean) establishes |h(s) - s*| ≤ c|s - s*|
  with c = (x-1)/x < 1 for all p ≥ 2, this excludes 2-cycles,
  3-cycles, and all higher periodic orbits.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace Pandrosion

/-! ## §221. Exclusion of 2-Cycles Under Contraction

If a map f satisfies |f(s) - r| ≤ c·|s - r| with c < 1,
then f has no 2-cycles: f(f(s)) = s implies s = r.

Proof: Apply contraction twice:
  |s - r| = |f(f(s)) - r| ≤ c·|f(s) - r| ≤ c²·|s - r|
Since c² < 1 and |s - r| ≥ 0, this forces |s - r| = 0.
-/

/-- **No 2-cycles under contraction.**
    If f contracts toward r by factor c < 1, then f(f(s)) = s implies s = r.
    This is a standard consequence of Banach's contraction principle applied
    to f∘f, but we prove it directly with an elementary argument. -/
theorem no_two_cycle (f : ℝ → ℝ) (r : ℝ) (c : ℝ)
    (hc_nn : 0 ≤ c) (hc_lt : c < 1)
    (h_contract : ∀ s, |f s - r| ≤ c * |s - r|)
    (s : ℝ) (h_period : f (f s) = s) :
    s = r := by
  -- By contradiction: assume s ≠ r, so |s - r| > 0
  by_contra h_ne
  have h_pos : (0 : ℝ) < |s - r| := abs_pos.mpr (sub_ne_zero.mpr h_ne)
  -- Apply contraction to s: |f(s) - r| ≤ c·|s - r|
  have h1 : |f s - r| ≤ c * |s - r| := h_contract s
  -- Apply contraction to f(s): |f(f(s)) - r| ≤ c·|f(s) - r|
  have h2 : |f (f s) - r| ≤ c * |f s - r| := h_contract (f s)
  -- Since f(f(s)) = s: |s - r| ≤ c·|f(s) - r| ≤ c²·|s - r|
  rw [h_period] at h2
  have h3 : |s - r| ≤ c ^ 2 * |s - r| := by
    calc |s - r|
      _ ≤ c * |f s - r| := h2
      _ ≤ c * (c * |s - r|) := by
          exact mul_le_mul_of_nonneg_left h1 hc_nn
      _ = c ^ 2 * |s - r| := by ring
  -- But c² < 1, so c²·|s - r| < |s - r|, contradiction
  have hc2 : c ^ 2 < 1 := by nlinarith [sq_nonneg c, sq_nonneg (1 - c)]
  have h4 : c ^ 2 * |s - r| < |s - r| := by
    calc c ^ 2 * |s - r|
      _ < 1 * |s - r| := by exact mul_lt_mul_of_pos_right hc2 h_pos
      _ = |s - r| := one_mul _
  linarith

/-! ## §222. Exclusion of All Periodic Orbits

Generalisation: no finite period n ≥ 1 is possible, since
iterating the contraction n times gives |f^n(s) - r| ≤ cⁿ|s - r|,
and cⁿ < 1 for c < 1.

We prove this for n-fold iteration via induction on the number
of contraction applications.
-/

/-- **Contraction composes:** n applications multiply the rate.
    |f^n(s) - r| ≤ cⁿ · |s - r| for any n. -/
theorem contraction_iterate (f : ℝ → ℝ) (r : ℝ) (c : ℝ)
    (hc_nn : 0 ≤ c)
    (h_contract : ∀ s, |f s - r| ≤ c * |s - r|)
    (s : ℝ) : ∀ n : ℕ, |f^[n] s - r| ≤ c ^ n * |s - r| := by
  intro n
  induction n with
  | zero => simp
  | succ k ih =>
    simp only [Function.iterate_succ']
    -- Goal type: |(f ∘ f^[k]) s - r| ≤ c ^ (k+1) * |s - r|
    -- which is |f (f^[k] s) - r| ≤ ...
    show |f (f^[k] s) - r| ≤ c ^ (k + 1) * |s - r|
    calc |f (f^[k] s) - r|
      _ ≤ c * |f^[k] s - r| := h_contract _
      _ ≤ c * (c ^ k * |s - r|) := by
          exact mul_le_mul_of_nonneg_left ih hc_nn
      _ = c ^ (k + 1) * |s - r| := by ring

/-- **No periodic orbits under contraction.**
    If f contracts toward r by factor c < 1, and f^n(s) = s
    for some n ≥ 1, then s = r.

    This completely characterises the dynamics: the fixed point
    is the only recurrent point. -/
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
  -- |s - r| ≤ cⁿ · |s - r| with cⁿ < 1
  have hcn : c ^ n < 1 := by
    exact pow_lt_one hc_nn hc_lt (by omega)
  have h_strict : c ^ n * |s - r| < |s - r| := by
    calc c ^ n * |s - r|
      _ < 1 * |s - r| := mul_lt_mul_of_pos_right hcn h_pos
      _ = |s - r| := one_mul _
  linarith

end Pandrosion
