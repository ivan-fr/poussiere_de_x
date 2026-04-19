/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP XIV: COMPLEX DYNAMICS AND JULIA SET CONJECTURE

  Explores the Complex Dynamics of the Pandrosion iteration.
  Unlike generic Newton-type maps of degree ≥ 4 which have chaotic
  fractal boundaries (McMullen 1987), numerical evidence suggests
  P_X has smooth basins (a non-fractal Julia set).
-/
import Mathlib.Data.Complex.Basic
import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.Calculus.FDeriv.Basic
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.MeasureTheory.Measure.Lebesgue.Complex
import Mathlib.Topology.Instances.Complex
import Mathlib.Tactic

open Complex MeasureTheory

namespace Pandrosion

/-! ## §300. The Rational Map -/

/-- The Pandrosion iteration over ℂ for p=3. -/
noncomputable def P3_complex (X : ℂ) (s : ℂ) : ℂ :=
  (s ^ 4 + 4 * X * s) / (3 * s ^ 3 + 2 * X)

/-- The induced one-variable map on the normalized cubic coordinate `y = s^3 / X`. -/
noncomputable def H3_complex (y : ℂ) : ℂ :=
  y * (y + 4) ^ 3 / (3 * y + 2) ^ 3

/-- The real restriction of the normalized cubic-coordinate map. -/
noncomputable def H3_real (y : ℝ) : ℝ :=
  y * (y + 4) ^ 3 / (3 * y + 2) ^ 3

/-- The real restriction commutes with the natural embedding into `ℂ`. -/
theorem H3_real_coe (y : ℝ) :
    H3_complex (y : ℂ) = (H3_real y : ℂ) := by
  norm_num [H3_complex, H3_real]

/-- Exact real error identity around the attracting fixed point `1`. -/
theorem H3_real_error_identity (y : ℝ) (hden : 3 * y + 2 ≠ 0) :
    H3_real y - 1 =
      (y - 1) * (y ^ 3 - 14 * y ^ 2 - 20 * y + 8) / (3 * y + 2) ^ 3 := by
  unfold H3_real
  field_simp [hden]
  ring

/-- A real free critical point is not the pole. -/
theorem H3_real_free_critical_den_ne_zero (y : ℝ)
    (hcrit : 3 * y ^ 2 - 16 * y + 8 = 0) :
    3 * y + 2 ≠ 0 := by
  intro hden
  have hcert : 3 * y ^ 2 - 16 * y + 8 = (3 * y + 2) * (y - 6) + 20 := by
    ring
  rw [hcert, hden] at hcrit
  norm_num at hcrit

/-- Real free critical points are positive. -/
theorem H3_real_free_critical_pos (y : ℝ)
    (hcrit : 3 * y ^ 2 - 16 * y + 8 = 0) :
    0 < y := by
  by_contra h
  have hy : y ≤ 0 := le_of_not_gt h
  linarith [sq_nonneg y]

/-- Real free critical points are bounded above by `16/3`. -/
theorem H3_real_free_critical_le (y : ℝ)
    (hcrit : 3 * y ^ 2 - 16 * y + 8 = 0) :
    y ≤ 16 / 3 := by
  by_contra h
  have hy : 16 / 3 < y := lt_of_not_ge h
  have hpos : 0 < y * (3 * y - 16) + 8 := by
    have hypos : 0 < y := by linarith
    have hfactor : 0 < 3 * y - 16 := by linarith
    linarith [mul_pos hypos hfactor]
  have hrewrite : 3 * y ^ 2 - 16 * y + 8 = y * (3 * y - 16) + 8 := by ring
  rw [hrewrite] at hcrit
  linarith

/-- Exact image of a real free critical point under `H3_real`. -/
theorem H3_real_free_critical_image (y : ℝ)
    (hcrit : 3 * y ^ 2 - 16 * y + 8 = 0) :
    H3_real y = (808 - 56 * y) / 729 := by
  have hden : 3 * y + 2 ≠ 0 := H3_real_free_critical_den_ne_zero y hcrit
  unfold H3_real
  field_simp [hden]
  nlinarith

/-- After one step, each real free critical point lands in `[1/2, 3/2]`. -/
theorem H3_real_free_critical_image_mem_interval (y : ℝ)
    (hcrit : 3 * y ^ 2 - 16 * y + 8 = 0) :
    (1 : ℝ) / 2 ≤ H3_real y ∧ H3_real y ≤ (3 : ℝ) / 2 := by
  rw [H3_real_free_critical_image y hcrit]
  have hyle : y ≤ 16 / 3 := H3_real_free_critical_le y hcrit
  constructor <;> linarith

/-- Complex free critical points of `H3_complex` are real. -/
theorem H3_complex_free_critical_im_zero (z : ℂ)
    (hcrit : 3 * z ^ 2 - 16 * z + 8 = 0) :
    z.im = 0 := by
  let a : ℝ := z.re
  let b : ℝ := z.im
  have hpow_re : (z ^ 2).re = a ^ 2 - b ^ 2 := by
    dsimp [a, b]
    rw [pow_two, Complex.mul_re]
    ring
  have hpow_im : (z ^ 2).im = 2 * a * b := by
    dsimp [a, b]
    rw [pow_two, Complex.mul_im]
    ring
  have hre : 3 * (a ^ 2 - b ^ 2) - 16 * a + 8 = 0 := by
    have h := congr_arg Complex.re hcrit
    norm_num at h
    rw [hpow_re] at h
    dsimp [a] at h
    ring_nf at h ⊢
    exact h
  have him : (6 * a - 16) * b = 0 := by
    have h := congr_arg Complex.im hcrit
    norm_num at h
    rw [hpow_im] at h
    dsimp [b] at h
    ring_nf at h ⊢
    exact h
  by_contra hb
  have hb_ne : b ≠ 0 := by
    exact hb
  have ha : a = 8 / 3 := by
    have hlin : 6 * a - 16 = 0 := (mul_eq_zero.mp him).resolve_right hb_ne
    linarith
  rw [ha] at hre
  ring_nf at hre
  linarith [sq_nonneg b]

/-- The real part of a complex free critical point satisfies the real critical equation. -/
theorem H3_complex_free_critical_re_eq (z : ℂ)
    (hcrit : 3 * z ^ 2 - 16 * z + 8 = 0) :
    3 * z.re ^ 2 - 16 * z.re + 8 = 0 := by
  let a : ℝ := z.re
  let b : ℝ := z.im
  have him : b = 0 := by
    dsimp [b]
    exact H3_complex_free_critical_im_zero z hcrit
  have hpow_re : (z ^ 2).re = a ^ 2 - b ^ 2 := by
    dsimp [a, b]
    rw [pow_two, Complex.mul_re]
    ring
  have hre : 3 * (a ^ 2 - b ^ 2) - 16 * a + 8 = 0 := by
    have h := congr_arg Complex.re hcrit
    norm_num at h
    rw [hpow_re] at h
    dsimp [a] at h
    ring_nf at h ⊢
    exact h
  dsimp [a, b] at hre him
  rw [him] at hre
  ring_nf at hre
  exact hre

/-- A complex free critical point is exactly its real part embedded in `ℂ`. -/
theorem H3_complex_free_critical_eq_re (z : ℂ)
    (hcrit : 3 * z ^ 2 - 16 * z + 8 = 0) :
    z = (z.re : ℂ) := by
  apply Complex.ext
  · simp
  · simpa using H3_complex_free_critical_im_zero z hcrit

/-- The image of every complex free critical point is the embedded real image
    of its real part. -/
theorem H3_complex_free_critical_image_real (z : ℂ)
    (hcrit : 3 * z ^ 2 - 16 * z + 8 = 0) :
    H3_complex z = (H3_real z.re : ℂ) := by
  rw [H3_complex_free_critical_eq_re z hcrit]
  exact H3_real_coe z.re

/-- Consequently, every complex free critical point maps into the real
    interval `[1/2, 3/2]` in the normalized coordinate. -/
theorem H3_complex_free_critical_image_mem_real_interval (z : ℂ)
    (hcrit : 3 * z ^ 2 - 16 * z + 8 = 0) :
    (1 : ℝ) / 2 ≤ (H3_complex z).re ∧ (H3_complex z).re ≤ (3 : ℝ) / 2 := by
  rw [H3_complex_free_critical_image_real z hcrit]
  simp
  simpa [one_div] using H3_real_free_critical_image_mem_interval z.re
    (H3_complex_free_critical_re_eq z hcrit)

/-- Zero is a fixed point of the cubic-coordinate map. -/
theorem H3_zero_fixed : H3_complex 0 = 0 := by
  simp [H3_complex]

/-- The target cubic coordinate `1` is fixed. -/
theorem H3_one_fixed : H3_complex 1 = 1 := by
  norm_num [H3_complex]

/-- Algebraic factorization of the fixed-point equation for `H3_complex`. -/
theorem H3_fixed_sub_factor (y : ℂ) :
    y * (y + 4) ^ 3 - y * (3 * y + 2) ^ 3 =
      2 * y * (1 - y) * (13 * y ^ 2 + 34 * y + 28) := by
  ring

/-- Away from the pole, the normalized map has exactly the listed fixed points. -/
theorem H3_fixed_iff (y : ℂ) (hden : (3 * y + 2) ^ 3 ≠ 0) :
    H3_complex y = y ↔ y = 0 ∨ y = 1 ∨ 13 * y ^ 2 + 34 * y + 28 = 0 := by
  unfold H3_complex
  rw [div_eq_iff hden]
  constructor
  · intro h
    have hsub : y * (y + 4) ^ 3 - y * (3 * y + 2) ^ 3 = 0 := by
      rw [h]
      ring
    have hfac : 2 * y * (1 - y) * (13 * y ^ 2 + 34 * y + 28) = 0 := by
      rw [← H3_fixed_sub_factor y]
      exact hsub
    have hcases :
        y = 0 ∨ 1 - y = 0 ∨ 13 * y ^ 2 + 34 * y + 28 = 0 := by
      have hmul₁ :
          2 * y * (1 - y) = 0 ∨ 13 * y ^ 2 + 34 * y + 28 = 0 :=
        mul_eq_zero.mp hfac
      cases hmul₁ with
      | inl hleft =>
          have hmul₂ : 2 * y = 0 ∨ 1 - y = 0 := mul_eq_zero.mp hleft
          cases hmul₂ with
          | inl hy =>
              left
              cases mul_eq_zero.mp hy with
              | inl htwo =>
                  norm_num at htwo
              | inr hy0 =>
                  exact hy0
          | inr hy =>
              right
              left
              exact hy
      | inr hq =>
          right
          right
          exact hq
    cases hcases with
    | inl hy =>
        exact Or.inl hy
    | inr hrest =>
        cases hrest with
        | inl hy =>
            right
            left
            rw [sub_eq_zero] at hy
            exact hy.symm
        | inr hq =>
            right
            right
            exact hq
  · intro h
    rw [← sub_eq_zero]
    rw [H3_fixed_sub_factor y]
    cases h with
    | inl hy =>
        rw [hy]
        ring
    | inr hrest =>
        cases hrest with
        | inl hy =>
            rw [hy]
            ring
        | inr hq =>
            rw [hq]
            ring

/-- The two extra fixed points are not poles of `H3_complex`. -/
theorem H3_extra_fixed_den_ne_zero (y : ℂ)
    (hq : 13 * y ^ 2 + 34 * y + 28 = 0) :
    3 * y + 2 ≠ 0 := by
  intro hden
  have hcert :
      9 * (13 * y ^ 2 + 34 * y + 28) = (39 * y + 76) * (3 * y + 2) + 100 := by
    ring
  rw [hq, hden] at hcert
  norm_num at hcert

/-- The extra fixed points are not the attracting target `1`. -/
theorem H3_extra_fixed_ne_one (y : ℂ)
    (hq : 13 * y ^ 2 + 34 * y + 28 = 0) :
    y ≠ 1 := by
  intro hy
  rw [hy] at hq
  norm_num at hq

/-- The extra fixed points are not the fixed exceptional point `0`. -/
theorem H3_extra_fixed_ne_zero (y : ℂ)
    (hq : 13 * y ^ 2 + 34 * y + 28 = 0) :
    y ≠ 0 := by
  intro hy
  rw [hy] at hq
  norm_num at hq

/-- Any root of the extra quadratic is indeed a fixed point of `H3_complex`. -/
theorem H3_extra_fixed_is_fixed (y : ℂ)
    (hq : 13 * y ^ 2 + 34 * y + 28 = 0) :
    H3_complex y = y := by
  have hden : (3 * y + 2) ^ 3 ≠ 0 := pow_ne_zero _ (H3_extra_fixed_den_ne_zero y hq)
  exact (H3_fixed_iff y hden).mpr (Or.inr (Or.inr hq))

/-- Convergence in the normalized cubic coordinate. -/
def H3_converges_to_one (y : ℂ) : Prop :=
  Filter.Tendsto (fun n => (H3_complex)^[n] y) Filter.atTop (nhds (1 : ℂ))

/-- A fixed point different from `1` cannot converge to `1`. -/
theorem H3_fixed_not_converges_to_one (y : ℂ)
    (hyfix : H3_complex y = y) (hyne : y ≠ 1) :
    ¬ H3_converges_to_one y := by
  intro ht
  have hconst_orbit :
      (fun n : ℕ => (H3_complex)^[n] y) = fun _ : ℕ => y := by
    funext n
    exact Function.iterate_fixed hyfix n
  have ht1 : Filter.Tendsto (fun _ : ℕ => y) Filter.atTop (nhds (1 : ℂ)) := by
    simpa [H3_converges_to_one, hconst_orbit] using ht
  have hty : Filter.Tendsto (fun _ : ℕ => y) Filter.atTop (nhds y) := tendsto_const_nhds
  have hlim : (1 : ℂ) = y := tendsto_nhds_unique ht1 hty
  exact hyne hlim.symm

/-- The normalized fixed point `0` is exceptional for convergence to `1`. -/
theorem H3_zero_not_converges_to_one : ¬ H3_converges_to_one 0 :=
  H3_fixed_not_converges_to_one 0 H3_zero_fixed (by norm_num)

/-- The extra quadratic fixed points are exceptional for convergence to `1`. -/
theorem H3_extra_fixed_not_converges_to_one (y : ℂ)
    (hq : 13 * y ^ 2 + 34 * y + 28 = 0) :
    ¬ H3_converges_to_one y :=
  H3_fixed_not_converges_to_one y
    (H3_extra_fixed_is_fixed y hq) (H3_extra_fixed_ne_one y hq)

/-- Algebraic numerator of the derivative of `H3_complex`. -/
theorem H3_derivative_numerator (y : ℂ) :
    let u := y * (y + 4) ^ 3
    let v := (3 * y + 2) ^ 3
    let u' := (y + 4) ^ 3 + y * (3 * (y + 4) ^ 2)
    let v' := 9 * (3 * y + 2) ^ 2
    u' * v - u * v' =
      (y + 4) ^ 2 * (3 * y + 2) ^ 2 * (3 * y ^ 2 - 16 * y + 8) := by
  intros
  ring

/-- The formal multiplier obtained from the quotient-rule numerator. -/
noncomputable def H3_formal_multiplier (y : ℂ) : ℂ :=
  ((y + 4) ^ 2 * (3 * y ^ 2 - 16 * y + 8)) / (3 * y + 2) ^ 4

/-- The normalized target fixed point is attracting with multiplier `-1/5`. -/
theorem H3_formal_multiplier_at_one : H3_formal_multiplier 1 = -(1 : ℂ) / 5 := by
  norm_num [H3_formal_multiplier]

/-- At the two non-target fixed points, the formal multiplier satisfies
    `25 m² + 235 m + 559 = 0`.

Over `ℂ`, these two multipliers are `-47/10 ± i * sqrt(27) / 10`, hence
repelling. The polynomial relation is the algebraic part needed before the
metric `Complex.normSq` statement.
-/
theorem H3_formal_multiplier_extra_fixed_poly (y : ℂ)
    (hden : 3 * y + 2 ≠ 0)
    (hq : 13 * y ^ 2 + 34 * y + 28 = 0) :
    25 * (H3_formal_multiplier y) ^ 2 + 235 * H3_formal_multiplier y + 559 = 0 := by
  unfold H3_formal_multiplier
  field_simp [hden]
  have hcert :
      25 * ((y + 4) ^ 2 * (3 * y ^ 2 - 16 * y + 8)) ^ 2 * (3 * y + 2) ^ 4 +
            235 * ((y + 4) ^ 2 * (3 * y ^ 2 - 16 * y + 8)) * ((3 * y + 2) ^ 4) ^ 2 +
          559 * (((3 * y + 2) ^ 4) ^ 2 * (3 * y + 2) ^ 4) =
      (13 * y ^ 2 + 34 * y + 28) *
        (36928 + 49952 * y + 328320 * y ^ 2 + 397640 * y ^ 3 +
          793720 * y ^ 4 + 778782 * y ^ 5 + 286533 * y ^ 6) * (3 * y + 2) ^ 4 := by
    ring
  rw [hcert, hq]
  ring

/-- The quadratic satisfied by the extra fixed-point multipliers forces their
    squared modulus to be exactly `559 / 25`, hence greater than `1`. -/
theorem extra_fixed_multiplier_normSq (m : ℂ)
    (hpoly : 25 * m ^ 2 + 235 * m + 559 = 0) :
    Complex.normSq m = 559 / 25 := by
  let a : ℝ := m.re
  let b : ℝ := m.im
  have hpow_re : (m ^ 2).re = a ^ 2 - b ^ 2 := by
    dsimp [a, b]
    rw [pow_two, Complex.mul_re]
    ring
  have hpow_im : (m ^ 2).im = 2 * a * b := by
    dsimp [a, b]
    rw [pow_two, Complex.mul_im]
    ring
  have hre : 25 * (a ^ 2 - b ^ 2) + 235 * a + 559 = 0 := by
    have h := congr_arg Complex.re hpoly
    norm_num at h
    rw [hpow_re] at h
    dsimp [a] at h
    nlinarith
  have him : (50 * a + 235) * b = 0 := by
    have h := congr_arg Complex.im hpoly
    norm_num at h
    rw [hpow_im] at h
    dsimp [b] at h
    nlinarith
  have hb_ne : b ≠ 0 := by
    intro hb
    have hreal : 25 * a ^ 2 + 235 * a + 559 = 0 := by
      nlinarith
    have hsquare : (10 * a + 47) ^ 2 + 27 = 0 := by
      nlinarith
    nlinarith [sq_nonneg (10 * a + 47)]
  have ha : a = -47 / 10 := by
    have hlin : 50 * a + 235 = 0 := by
      exact mul_eq_zero.mp him |>.resolve_right hb_ne
    nlinarith
  have hb_sq : b ^ 2 = 27 / 100 := by
    rw [ha] at hre
    norm_num at hre ⊢
    nlinarith
  rw [Complex.normSq_apply]
  change a * a + b * b = 559 / 25
  rw [ha]
  rw [show b * b = b ^ 2 by ring, hb_sq]
  norm_num

/-- The two non-target fixed-point multipliers are repelling in squared norm. -/
theorem extra_fixed_multiplier_normSq_gt_one (m : ℂ)
    (hpoly : 25 * m ^ 2 + 235 * m + 559 = 0) :
    1 < Complex.normSq m := by
  rw [extra_fixed_multiplier_normSq m hpoly]
  norm_num

/-- Same repulsion statement in the usual complex norm. -/
theorem extra_fixed_multiplier_norm_gt_one (m : ℂ)
    (hpoly : 25 * m ^ 2 + 235 * m + 559 = 0) :
    1 < ‖m‖ := by
  have hsq : 1 < ‖m‖ ^ 2 := by
    rw [← Complex.normSq_eq_norm_sq]
    exact extra_fixed_multiplier_normSq_gt_one m hpoly
  by_contra hnot
  have hle : ‖m‖ ≤ 1 := le_of_not_gt hnot
  have hnonneg : 0 ≤ ‖m‖ := norm_nonneg _
  nlinarith

/-- The two non-target fixed points of `H3_complex` are repelling in the
    formal multiplier sense. -/
theorem H3_extra_fixed_multiplier_normSq_gt_one (y : ℂ)
    (hden : 3 * y + 2 ≠ 0)
    (hq : 13 * y ^ 2 + 34 * y + 28 = 0) :
    1 < Complex.normSq (H3_formal_multiplier y) := by
  exact extra_fixed_multiplier_normSq_gt_one _
    (H3_formal_multiplier_extra_fixed_poly y hden hq)

/-- The two non-target fixed points of `H3_complex` are repelling in the
    usual complex norm. -/
theorem H3_extra_fixed_multiplier_norm_gt_one (y : ℂ)
    (hden : 3 * y + 2 ≠ 0)
    (hq : 13 * y ^ 2 + 34 * y + 28 = 0) :
    1 < ‖H3_formal_multiplier y‖ := by
  exact extra_fixed_multiplier_norm_gt_one _
    (H3_formal_multiplier_extra_fixed_poly y hden hq)

/-- The origin is a fixed point of the totalized complex Pandrosion map. -/
theorem P3_complex_zero_fixed (X : ℂ) : P3_complex X 0 = 0 := by
  simp [P3_complex]

/-- Every cubic root of `X` is a fixed point of `P3_complex`. -/
theorem P3_complex_root_fixed (X r : ℂ) (hX : X ≠ 0) (hr : r ^ 3 = X) :
    P3_complex X r = r := by
  have hr_ne : r ≠ 0 := by
    intro hr0
    apply hX
    rw [← hr, hr0]
    norm_num
  unfold P3_complex
  rw [← hr]
  have hden : 3 * r ^ 3 + 2 * r ^ 3 = 5 * r ^ 3 := by ring
  have hnum : r ^ 4 + 4 * r ^ 3 * r = 5 * r ^ 4 := by ring
  rw [hden, hnum]
  field_simp [hr_ne]
  ring

/-- Away from the pole, the only finite fixed points are `0` and the cubic roots. -/
theorem P3_complex_fixed_iff (X z : ℂ) (hden : 3 * z ^ 3 + 2 * X ≠ 0) :
    P3_complex X z = z ↔ z = 0 ∨ z ^ 3 = X := by
  unfold P3_complex
  rw [div_eq_iff hden]
  constructor
  · intro h
    have hsub : z ^ 4 + 4 * X * z - z * (3 * z ^ 3 + 2 * X) = 0 := by
      rw [h]
      ring
    have hfac : 2 * z * (X - z ^ 3) = 0 := by
      rw [← hsub]
      ring
    have hcases : z = 0 ∨ X - z ^ 3 = 0 := by
      have hmul : 2 * z = 0 ∨ X - z ^ 3 = 0 := mul_eq_zero.mp hfac
      cases hmul with
      | inl hz =>
          left
          cases mul_eq_zero.mp hz with
          | inl htwo =>
              norm_num at htwo
          | inr hz0 =>
              exact hz0
      | inr hx =>
          right
          exact hx
    cases hcases with
    | inl hz => exact Or.inl hz
    | inr hx =>
        right
        rw [sub_eq_zero] at hx
        exact hx.symm
  · intro h
    cases h with
    | inl hz =>
        rw [hz]
        ring
    | inr hroot =>
        rw [← hroot]
        ring

/-- Away from the pole, the normalized cubic coordinate is semi-conjugate to `H3_complex`.

If `y = s^3 / X`, then the next normalized coordinate is `H3_complex y`.
This is the algebraic reduction that turns the complex-plane conjecture into
the study of a single rational map on the Riemann sphere.
-/
theorem P3_complex_cubic_coordinate (X s : ℂ) (hX : X ≠ 0)
    (hden : 3 * s ^ 3 + 2 * X ≠ 0) :
    (P3_complex X s) ^ 3 / X = H3_complex (s ^ 3 / X) := by
  have hyden : 3 * (s ^ 3 / X) + 2 ≠ 0 := by
    intro h
    apply hden
    have hmul : (3 * (s ^ 3 / X) + 2) * X = 0 := by rw [h, zero_mul]
    field_simp [hX] at hmul
    simpa [mul_comm, mul_left_comm, mul_assoc] using hmul
  unfold P3_complex H3_complex
  field_simp [hX, hden, hyden]
  ring

/-! ## §301. Critical Points Identity -/

/-- **Derivative Algebraic Form.**
    The critical points of P_X are governed by the numerator of its derivative.
    Since P_X is a rational function, P'_X(s) = 0 iff 3s^6 - 16Xs^3 + 8X^2 = 0.
    Unlike super-attracting Newton maps, P'_X(∛X) = -1/5 ≠ 0.
-/
theorem P3_derivative_numerator (X s : ℂ) :
    let u := s ^ 4 + 4 * X * s
    let v := 3 * s ^ 3 + 2 * X
    let u' := 4 * s ^ 3 + 4 * X
    let v' := 9 * s ^ 2
    u' * v - u * v' = 3 * s ^ 6 - 16 * X * s ^ 3 + 8 * X ^ 2 := by
  intros
  calc (4 * s ^ 3 + 4 * X) * (3 * s ^ 3 + 2 * X) - (s ^ 4 + 4 * X * s) * (9 * s ^ 2)
      = 12 * s ^ 6 + 8 * X * s ^ 3 + 12 * X * s ^ 3 + 8 * X ^ 2 - 9 * s ^ 6 - 36 * X * s ^ 3 := by ring
    _ = 3 * s ^ 6 - 16 * X * s ^ 3 + 8 * X ^ 2 := by ring

/-! ## §302. The Formal Conjecture -/

/-- Definition of convergence to a target root of X. -/
def converges_to_root (X : ℂ) (s : ℂ) : Prop :=
  ∃ r : ℂ, r ^ 3 = X ∧ Filter.Tendsto (fun n => (P3_complex X)^[n] s) Filter.atTop (nhds r)

/-- For `X ≠ 0`, the exceptional set is genuinely nonempty: `0` is fixed,
    so it cannot converge to a nonzero cubic root. -/
theorem not_converges_to_root_zero (X : ℂ) (hX : X ≠ 0) :
    ¬ converges_to_root X 0 := by
  rintro ⟨r, hr, ht⟩
  have hconst_orbit :
      (fun n : ℕ => (P3_complex X)^[n] (0 : ℂ)) = fun _ : ℕ => (0 : ℂ) := by
    funext n
    exact Function.iterate_fixed (P3_complex_zero_fixed X) n
  have ht0 : Filter.Tendsto (fun _ : ℕ => (0 : ℂ)) Filter.atTop (nhds r) := by
    simpa [hconst_orbit] using ht
  have h0 : Filter.Tendsto (fun _ : ℕ => (0 : ℂ)) Filter.atTop (nhds (0 : ℂ)) :=
    tendsto_const_nhds
  have hr0 : r = 0 := tendsto_nhds_unique ht0 h0
  apply hX
  rw [← hr, hr0]
  norm_num

/-- The fixed exceptional point `0` has zero Lebesgue measure. -/
theorem volume_singleton_zero_complex : MeasureTheory.volume ({0} : Set ℂ) = 0 := by
  simp

/-- **Open Problem: Absence of Chaos (McMullen Exemption)**
    It is conjectured that the Fatou set of P_X has full measure in ℂ.
    In other words, almost all starting points converge to a finite cycle
    (and empirically, only to the roots of X).
    Because Mathlib currently lacks Montel's theorem and Riemann surface dynamics,
    we encode this as a formal axiom serving as an open research question.
-/
axiom pandrosion_fatou_full_measure (X : ℂ) (hX : X ≠ 0) :
  MeasureTheory.volume {s : ℂ | ¬ converges_to_root X s} = 0

end Pandrosion
