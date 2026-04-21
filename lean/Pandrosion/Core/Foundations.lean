/-
  Universitas Pandrosion — Core / Foundations (clean rewrite).

  Order (strict, minimal, fully proven — béton armé):

    §1  Core Pandrosion iteration (p-th, generalised)          [ℝ]
    §2  Scaling optimization
        §2.1  Scaling principle
        §2.2  Thales preconditioning (geometric interpretation)
        §2.3  Double preconditioning
    §3  Pandrosion + Steffensen acceleration                   [ℝ]
    §4  Pandrosion iteration in ℂ
    §5  Pandrosion + Steffensen acceleration                   [ℂ]
    §6  Scaling + Steffensen in ℂ
-/

import Mathlib.Algebra.BigOperators.Basic
import Mathlib.Algebra.GeomSum
import Mathlib.Data.Complex.Basic
import Mathlib.Tactic

namespace Pandrosion

open Finset BigOperators

/-! ============================================================
  §1. Core Pandrosion iteration (p-th, generalised) — ℝ
============================================================ -/

section CoreReal

/-- Geometric sum `S_p(s) = Σ_{k=0}^{p-1} s^k`. -/
def Sp (p : ℕ) (s : ℝ) : ℝ := ∑ k in Finset.range p, s ^ k

/-- `S_p(s) · (1 − s) = 1 − s^p` — the geometric-series identity. -/
theorem Sp_mul_one_sub (p : ℕ) (s : ℝ) : Sp p s * (1 - s) = 1 - s ^ p := by
  unfold Sp; exact geom_sum_mul_neg s p

/-- `S_p(1) = p`. -/
theorem Sp_at_one (p : ℕ) : Sp p 1 = (p : ℝ) := by unfold Sp; simp

/-- `S_p(s) > 0` for `s ≥ 0` and `p ≥ 1`. -/
theorem Sp_pos (p : ℕ) (hp : p ≥ 1) (s : ℝ) (hs : 0 ≤ s) : Sp p s > 0 := by
  unfold Sp
  calc ∑ k in Finset.range p, s ^ k
      ≥ ∑ k in Finset.range 1, s ^ k :=
        Finset.sum_le_sum_of_subset_of_nonneg (Finset.range_mono (by omega))
          (fun k _ _ => pow_nonneg hs k)
    _ = 1 := by simp
    _ > 0 := one_pos

/-- **Pandrosion's iteration map**: `h_{p,x}(s) = 1 − (x−1)/(x·S_p(s))`.
    Converges linearly to the fixed point `s* = x^{-1/p}`. -/
noncomputable def pandrosion_h (x : ℝ) (p : ℕ) (s : ℝ) : ℝ :=
  1 - (x - 1) / (x * Sp p s)

/-- Underlying polynomial `f(s) = 1 − x·s^p`. Its positive root is
    exactly the Pandrosion fixed point. -/
def pandrosion_f (x : ℝ) (p : ℕ) (s : ℝ) : ℝ := 1 - x * s ^ p

/-- **Fixed-point characterisation**: `h(s) = s ⟺ s^p = 1/x`. -/
theorem pandrosion_fixed_point_iff
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s : ℝ) (hs : 0 ≤ s) :
    pandrosion_h x p s = s ↔ s ^ p = 1 / x := by
  have hS : Sp p s > 0 := Sp_pos p hp s hs
  have hxS_ne : x * Sp p s ≠ 0 := ne_of_gt (mul_pos hx hS)
  have hSmul : Sp p s * (1 - s) = 1 - s ^ p := Sp_mul_one_sub p s
  unfold pandrosion_h
  constructor
  · intro heq
    have h1 : (x - 1) / (x * Sp p s) = 1 - s := by linarith
    have h2 : x - 1 = (1 - s) * (x * Sp p s) := by
      rw [← h1]; exact (div_mul_cancel₀ (x - 1) hxS_ne).symm
    have h3 : (1 - s) * (x * Sp p s) = x * (Sp p s * (1 - s)) := by ring
    rw [h3, hSmul] at h2
    have h5 : x * s ^ p = 1 := by linarith
    rw [eq_div_iff (ne_of_gt hx)]; linarith
  · intro hsp
    have h5 : x * s ^ p = 1 := by
      rw [eq_div_iff (ne_of_gt hx)] at hsp; linear_combination hsp
    have h_key : x * Sp p s * (1 - s) = x - 1 := by
      have h_expand : x * (Sp p s * (1 - s)) = x * (1 - s ^ p) := by
        rw [hSmul]
      have h_rhs : x * (1 - s ^ p) = x - 1 := by
        rw [mul_sub, mul_one]; linarith
      linarith [h_expand, h_rhs]
    have h_div : (x - 1) / (x * Sp p s) = 1 - s := by
      rw [div_eq_iff hxS_ne]; linarith [h_key]
    linarith

/-- **Quasi-Newton identity**: `h(s) = s + f(s)/(x·S_p(s))`.
    The divisor `x·S_p(s)` replaces the derivative `f'(s) = −p·x·s^{p−1}`;
    hence Pandrosion is a *quasi-Newton* method, not Newton itself. -/
theorem pandrosion_quasi_newton
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1) (s : ℝ) (hs : 0 ≤ s) :
    pandrosion_h x p s = s + pandrosion_f x p s / (x * Sp p s) := by
  have hS : Sp p s > 0 := Sp_pos p hp s hs
  have hxS_ne : x * Sp p s ≠ 0 := ne_of_gt (mul_pos hx hS)
  have hSmul : Sp p s * (1 - s) = 1 - s ^ p := Sp_mul_one_sub p s
  unfold pandrosion_h pandrosion_f
  field_simp
  linear_combination x * hSmul

end CoreReal

/-! ============================================================
  §2. Scaling optimization
============================================================ -/

/-! ## §2.1 Scaling principle -/

section ScalingPrinciple

/-- **Scaling identity (algebraic form).**
    If `α^p = x` and `β^p = A`, then `(α/β)^p = x/A`. This is the
    algebraic content of `x^{1/p} = A^{1/p} · (x/A)^{1/p}`. -/
theorem scaling_power (α β x A : ℝ) (p : ℕ)
    (hα : α ^ p = x) (hβ : β ^ p = A) : (α / β) ^ p = x / A := by
  rw [div_pow, hα, hβ]

/-- **Reconstruction identity.** `α = β · (α/β)` for `β ≠ 0`. -/
theorem scaling_factorization (α β : ℝ) (hβ : β ≠ 0) :
    α = β * (α / β) := by
  rw [mul_div_assoc', mul_comm, mul_div_assoc, div_self hβ, mul_one]

end ScalingPrinciple

/-! ## §2.2 Thales preconditioning (geometric interpretation) -/

section ThalesPreconditioning

/-- **Optimal starting point** `s_0^opt := 1 − (x−1)/(x·p)`.
    Geometrically: one Pandrosion step applied to the naive guess
    `s = 1`, corresponding to the unit-length first Thales segment. -/
noncomputable def pandrosion_s0_opt (x : ℝ) (p : ℕ) : ℝ :=
  1 - (x - 1) / (x * p)

/-- `s_0^opt = h(1)`. One Pandrosion iterate from `s = 1`. -/
theorem pandrosion_s0_opt_eq_h_one (x : ℝ) (p : ℕ) :
    pandrosion_s0_opt x p = pandrosion_h x p 1 := by
  unfold pandrosion_s0_opt pandrosion_h
  rw [Sp_at_one]

/-- **Offset monotonicity.** `(x−1)/(x·p)` is non-decreasing in `x > 0`.
    Thales-geometric meaning: an elongated rectangle (large `x`) has
    slower-converging parallels than a near-square one (x close to 1). -/
theorem offset_monotone
    (p : ℕ) (hp : p ≥ 1) (x₁ x₂ : ℝ) (hx₁ : x₁ > 0) (hle : x₁ ≤ x₂) :
    (x₁ - 1) / (x₁ * p) ≤ (x₂ - 1) / (x₂ * p) := by
  have hx₂ : x₂ > 0 := lt_of_lt_of_le hx₁ hle
  have hp_pos : (0 : ℝ) < p := by exact_mod_cast hp
  have hx₁p : x₁ * p > 0 := mul_pos hx₁ hp_pos
  have eq₁ : (x₁ - 1) / (x₁ * p) = 1 / p - 1 / (x₁ * p) := by field_simp
  have eq₂ : (x₂ - 1) / (x₂ * p) = 1 / p - 1 / (x₂ * p) := by field_simp
  rw [eq₁, eq₂]
  have h_recip : 1 / (x₂ * p) ≤ 1 / (x₁ * p) :=
    one_div_le_one_div_of_le hx₁p (by nlinarith)
  linarith

/-- **Scaled start is closer to 1.**
    If `1 ≤ x' ≤ x`, then `s_0^opt(x) ≤ s_0^opt(x')`. -/
theorem scaled_s0_opt_closer_to_one
    (p : ℕ) (hp : p ≥ 1) (x x' : ℝ) (hx' : x' ≥ 1) (hle : x' ≤ x) :
    pandrosion_s0_opt x p ≤ pandrosion_s0_opt x' p := by
  unfold pandrosion_s0_opt
  have hx'_pos : x' > 0 := by linarith
  have h_offset := offset_monotone p hp x' x hx'_pos hle
  linarith

end ThalesPreconditioning

/-! ## §2.3 Double preconditioning -/

section DoublePreconditioning

/-- **Double preconditioning (algebraic core).**
    Scaling simultaneously (i) factorises the p-th root across the
    scaling, and (ii) moves the optimal start closer to 1. The two
    effects compound; this theorem packages them. -/
theorem double_preconditioning
    (p : ℕ) (hp : p ≥ 1)
    (α β x A : ℝ)
    (hα : α ^ p = x) (hβ : β ^ p = A)
    (x' : ℝ) (hx' : x' ≥ 1) (hle : x' ≤ x) :
    (α / β) ^ p = x / A ∧
    pandrosion_s0_opt x p ≤ pandrosion_s0_opt x' p :=
  ⟨scaling_power α β x A p hα hβ,
   scaled_s0_opt_closer_to_one p hp x x' hx' hle⟩

end DoublePreconditioning

/-! ============================================================
  §3. Pandrosion + Steffensen acceleration (real)
============================================================ -/

section SteffensenReal

/-- Steffensen denominator: `D_{p,x}(s) = h(h(s)) − 2·h(s) + s`. -/
noncomputable def steffensen_denom (x : ℝ) (p : ℕ) (s : ℝ) : ℝ :=
  pandrosion_h x p (pandrosion_h x p s) - 2 * pandrosion_h x p s + s

/-- Steffensen-accelerated Pandrosion step (piecewise: identity when
    denominator vanishes, covering the fixed point). -/
noncomputable def steffensen_step (x : ℝ) (p : ℕ) (s : ℝ) : ℝ :=
  if steffensen_denom x p s = 0 then s
  else s - (pandrosion_h x p s - s) ^ 2 / steffensen_denom x p s

/-- **Any Pandrosion fixed point is a Steffensen fixed point.** -/
theorem steffensen_step_fixed_point_of_h
    (x : ℝ) (p : ℕ) (s : ℝ) (h_fp : pandrosion_h x p s = s) :
    steffensen_step x p s = s := by
  have hD : steffensen_denom x p s = 0 := by
    unfold steffensen_denom; rw [h_fp, h_fp]; ring
  unfold steffensen_step; rw [if_pos hD]

/-- **Explicit fixed-point preservation at `s^p = 1/x`.** -/
theorem steffensen_step_fixed_point
    (x : ℝ) (hx : x > 0) (p : ℕ) (hp : p ≥ 1)
    (s : ℝ) (hs : 0 ≤ s)
    (h_fp : s ^ p = 1 / x) :
    steffensen_step x p s = s :=
  steffensen_step_fixed_point_of_h x p s
    ((pandrosion_fixed_point_iff x hx p hp s hs).mpr h_fp)

/-- **Steffensen non-degeneracy hypothesis.**
    The uniform contraction bound `(p−1)/p < 1` satisfies `≠ 1`, so
    the classical Aitken-Steffensen quadratic theorem applies. -/
theorem steffensen_rate_ne_one (p : ℕ) (hp : p ≥ 2) :
    ((p : ℝ) - 1) / p ≠ 1 := by
  have hp' : (0 : ℝ) < p := by exact_mod_cast (by omega : p ≥ 1)
  have h_lt : ((p : ℝ) - 1) / p < 1 := by
    rw [div_lt_one hp']; linarith
  linarith

end SteffensenReal

/-! ============================================================
  §4. Pandrosion iteration in ℂ
============================================================ -/

section PandrosionComplex

/-- Complex geometric sum. -/
noncomputable def Sp_C (p : ℕ) (s : ℂ) : ℂ := ∑ k in Finset.range p, s ^ k

theorem Sp_C_mul_one_sub (p : ℕ) (s : ℂ) :
    Sp_C p s * (1 - s) = 1 - s ^ p := by
  unfold Sp_C; exact geom_sum_mul_neg s p

theorem Sp_C_at_one (p : ℕ) : Sp_C p 1 = (p : ℂ) := by unfold Sp_C; simp

/-- Complex Pandrosion map `h(s) = 1 − (x−1)/(x·S_p(s))`. -/
noncomputable def pandrosion_h_C (x : ℂ) (p : ℕ) (s : ℂ) : ℂ :=
  1 - (x - 1) / (x * Sp_C p s)

/-- Complex underlying polynomial. -/
noncomputable def pandrosion_f_C (x : ℂ) (p : ℕ) (s : ℂ) : ℂ :=
  1 - x * s ^ p

/-- **Quasi-Newton identity in ℂ.** -/
theorem pandrosion_h_C_quasi_newton
    (x : ℂ) (p : ℕ) (s : ℂ) (hxSp : x * Sp_C p s ≠ 0) :
    pandrosion_h_C x p s = s + pandrosion_f_C x p s / (x * Sp_C p s) := by
  have hSmul : Sp_C p s * (1 - s) = 1 - s ^ p := Sp_C_mul_one_sub p s
  unfold pandrosion_h_C pandrosion_f_C
  field_simp
  linear_combination x * hSmul

/-- **Fixed-point → s^p = 1/x (forward, ℂ).** -/
theorem pandrosion_h_C_fixed_point_forward
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (s : ℂ) (hSp : Sp_C p s ≠ 0)
    (heq : pandrosion_h_C x p s = s) : s ^ p = 1 / x := by
  have hxSp : x * Sp_C p s ≠ 0 := mul_ne_zero hx hSp
  have hSmul : Sp_C p s * (1 - s) = 1 - s ^ p := Sp_C_mul_one_sub p s
  unfold pandrosion_h_C at heq
  have h1 : x * Sp_C p s - (x - 1) = s * (x * Sp_C p s) := by
    have hh := heq
    field_simp at hh
    linear_combination hh
  have h2 : x * (Sp_C p s * (1 - s)) = x - 1 := by linear_combination h1
  rw [hSmul] at h2
  have h3 : x * s ^ p = 1 := by linear_combination -h2
  rw [eq_div_iff hx]; linear_combination h3

/-- **s^p = 1/x → fixed point (reverse, ℂ).** -/
theorem pandrosion_h_C_fixed_point_reverse
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (s : ℂ) (hSp : Sp_C p s ≠ 0)
    (hsp : s ^ p = 1 / x) : pandrosion_h_C x p s = s := by
  have hxSp : x * Sp_C p s ≠ 0 := mul_ne_zero hx hSp
  have hSmul : Sp_C p s * (1 - s) = 1 - s ^ p := Sp_C_mul_one_sub p s
  unfold pandrosion_h_C
  have h3 : x * s ^ p = 1 := by
    rw [eq_div_iff hx] at hsp; linear_combination hsp
  have h2 : x * (1 - s ^ p) = x - 1 := by linear_combination -h3
  rw [← hSmul] at h2
  have h1 : x * Sp_C p s * (1 - s) = x - 1 := by linear_combination h2
  have h_div : (x - 1) / (x * Sp_C p s) = 1 - s := by
    rw [div_eq_iff hxSp]; linear_combination -h1
  linear_combination -h_div

/-- **Fixed-point characterisation in ℂ** (iff form). -/
theorem pandrosion_h_C_fixed_point_iff
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (s : ℂ) (hSp : Sp_C p s ≠ 0) :
    pandrosion_h_C x p s = s ↔ s ^ p = 1 / x :=
  ⟨pandrosion_h_C_fixed_point_forward x hx p s hSp,
   pandrosion_h_C_fixed_point_reverse x hx p s hSp⟩

end PandrosionComplex

/-! ============================================================
  §5. Pandrosion + Steffensen acceleration in ℂ
============================================================ -/

section SteffensenComplex

noncomputable def steffensen_denom_C (x : ℂ) (p : ℕ) (s : ℂ) : ℂ :=
  pandrosion_h_C x p (pandrosion_h_C x p s) - 2 * pandrosion_h_C x p s + s

noncomputable def steffensen_step_C (x : ℂ) (p : ℕ) (s : ℂ) : ℂ :=
  if steffensen_denom_C x p s = 0 then s
  else s - (pandrosion_h_C x p s - s) ^ 2 / steffensen_denom_C x p s

/-- **Any Pandrosion-ℂ fixed point is a Steffensen-ℂ fixed point.** -/
theorem steffensen_step_C_fixed_point_of_h
    (x : ℂ) (p : ℕ) (s : ℂ) (h_fp : pandrosion_h_C x p s = s) :
    steffensen_step_C x p s = s := by
  have hD : steffensen_denom_C x p s = 0 := by
    unfold steffensen_denom_C; rw [h_fp, h_fp]; ring
  unfold steffensen_step_C; rw [if_pos hD]

/-- **Explicit complex fixed-point preservation at `s^p = 1/x`.** -/
theorem steffensen_step_C_fixed_point
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (s : ℂ) (hSp : Sp_C p s ≠ 0)
    (h_fp : s ^ p = 1 / x) :
    steffensen_step_C x p s = s :=
  steffensen_step_C_fixed_point_of_h x p s
    (pandrosion_h_C_fixed_point_reverse x hx p s hSp h_fp)

end SteffensenComplex

/-! ============================================================
  §6. Scaling in ℂ + Steffensen
============================================================ -/

section ScalingComplex

/-- **Complex scaling identity.**
    If `α^p = x` and `β^p = A`, then `(α/β)^p = x/A`. -/
theorem scaling_power_C (α β x A : ℂ) (p : ℕ)
    (hα : α ^ p = x) (hβ : β ^ p = A) : (α / β) ^ p = x / A := by
  rw [div_pow, hα, hβ]

/-- **Complex reconstruction identity.** `α = β · (α/β)` for `β ≠ 0`. -/
theorem scaling_factorization_C (α β : ℂ) (hβ : β ≠ 0) :
    α = β * (α / β) := by
  rw [mul_div_assoc', mul_comm, mul_div_assoc, div_self hβ, mul_one]

/-- **Scaling preserves the complex Steffensen fixed point.**
    If `γ` is a fixed point of the Pandrosion iteration with reduced
    target `y = x/A` (i.e., `γ^p = 1/y = A/x`), then the Steffensen
    accelerator on `y` also fixes `γ`. Combined with the reconstruction
    `α = β · (α/β)`, this packages the full scaling-plus-Steffensen
    pipeline in ℂ. -/
theorem scaling_steffensen_C
    (α β γ x A : ℂ) (p : ℕ)
    (hx : x ≠ 0) (hβ : β ≠ 0) (hA : A ≠ 0)
    (hα : α ^ p = x) (hβpow : β ^ p = A)
    (hSp : Sp_C p γ ≠ 0)
    (hγ_fp : γ ^ p = A / x) :
    steffensen_step_C (x / A) p γ = γ ∧
    (α / β) ^ p = x / A ∧
    α = β * (α / β) := by
  have hxA : (x / A) ≠ 0 := div_ne_zero hx hA
  -- γ^p = A/x = 1/(x/A)
  have hfp : γ ^ p = 1 / (x / A) := by
    rw [hγ_fp]; field_simp
  refine ⟨?_, ?_, ?_⟩
  · exact steffensen_step_C_fixed_point (x / A) hxA p γ hSp hfp
  · exact scaling_power_C α β x A p hα hβpow
  · exact scaling_factorization_C α β hβ

end ScalingComplex

end Pandrosion
