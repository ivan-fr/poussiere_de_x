/-
  Universitas Pandrosion — §61. **Unconditional building block for
  `ComplexJuliaP3MeasureZero 2 α`.**

  Empiriquement (scans Python 360 000 pts + zoom ε = 10⁻⁹) au paramètre
  canonique `x = 2, p = 3` : l'ensemble divergent se réduit exactement à
  l'orbite singulière (itérés arrière des zéros de `Sp_C 3`).

  **Contribution inconditionnelle** de ce module :

    (i)  `Sp_C 3 z = 0 ↔ z³ = 1 ∧ z ≠ 1` (identité algébrique) ;
    (ii) l'ensemble des zéros est `{ω, ω²}` où `ω = e^{2πi/3}` ;
    (iii) il est fini, de cardinalité ≤ 2 ;
    (iv) il est de mesure de Lebesgue nulle.

  Ces briques élémentaires sont les **zéros directs** du dénominateur
  de `h = pandrosion_h_C 2 3`. La montée vers `ComplexJuliaP3MeasureZero`
  (= `McMullenAEEntry 3 2 α`) nécessiterait la dénombrabilité des
  **préimages itérées** + un argument de dynamique (Fatou-exhaustion),
  que nous laissons en conjecture résiduelle dans un module aval.

  Contents.

    §61.1  `Sp_C_p3_zero_iff_cube_root` — factorisation clé.

    §61.2  `Sp_C_p3_omega3_zero`, `Sp_C_p3_omega3_sq_zero`.

    §61.3  `Sp_C_p3_zeros_subset` — ⊆ `{ω, ω²}`.

    §61.4  `Sp_C_p3_zeros_finite` — fini.

    §61.5  `Sp_C_p3_zeros_measure_zero` — mesure nulle.
-/

import Pandrosion.Core.UniversalConjectureC

namespace Pandrosion

open Complex MeasureTheory

/-! ============================================================
  §61.1  Factorisation algébrique clé
============================================================ -/

/-- **`Sp_C 3 z = 0 ↔ z³ = 1 ∧ z ≠ 1`.** -/
theorem Sp_C_p3_zero_iff_cube_root (z : ℂ) :
    Sp_C 3 z = 0 ↔ z ^ 3 = 1 ∧ z ≠ 1 := by
  rw [Sp_C_p3]
  constructor
  · intro h
    have h_factor : z ^ 3 - 1 = (z - 1) * (1 + z + z ^ 2) := by ring
    have h_z3_sub : z ^ 3 - 1 = 0 := by rw [h_factor, h]; ring
    refine ⟨?_, ?_⟩
    · linear_combination h_z3_sub
    · intro h_eq
      rw [h_eq] at h
      norm_num at h
  · rintro ⟨h_cube, h_ne⟩
    have h_factor : z ^ 3 - 1 = (z - 1) * (1 + z + z ^ 2) := by ring
    have h_prod : (z - 1) * (1 + z + z ^ 2) = 0 := by
      rw [← h_factor]
      linear_combination h_cube
    have h_sub_ne : z - 1 ≠ 0 := sub_ne_zero.mpr h_ne
    exact (mul_eq_zero.mp h_prod).resolve_left h_sub_ne

/-! ============================================================
  §61.2  `ω³ = 1`, `ω ≠ 1`, et zeros spécifiques
============================================================ -/

/-- **`ω³ = 1`.** -/
theorem omega3_cube : omega3 ^ 3 = 1 := by
  unfold omega3
  rw [show ((Complex.exp (2 * (Real.pi : ℂ) * Complex.I / 3)) ^ 3 : ℂ)
      = Complex.exp (3 * (2 * (Real.pi : ℂ) * Complex.I / 3)) from by
        rw [← Complex.exp_nat_mul]; push_cast; ring_nf]
  have h_simp : (3 : ℂ) * (2 * (Real.pi : ℂ) * Complex.I / 3)
              = 2 * (Real.pi : ℂ) * Complex.I := by ring
  rw [h_simp, Complex.exp_two_pi_mul_I]

/-- **`ω ≠ 1`** (car `Im ω = sin(2π/3) > 0`). -/
theorem omega3_ne_one : omega3 ≠ 1 := by
  intro h
  have h_one_im : (1 : ℂ).im = 0 := by simp
  have h_eq_im : omega3.im = 0 := by rw [h]; exact h_one_im
  have h_re0 : (2 * ((Real.pi : ℝ) : ℂ) * Complex.I / 3).re = 0 := by
    simp [Complex.div_re, Complex.mul_re, Complex.mul_im,
          Complex.I_re, Complex.I_im, Complex.normSq_apply,
          Complex.ofReal_re, Complex.ofReal_im]
  have h_im_val : (2 * ((Real.pi : ℝ) : ℂ) * Complex.I / 3).im
                = 2 * Real.pi / 3 := by
    simp [Complex.div_im, Complex.mul_re, Complex.mul_im,
          Complex.I_re, Complex.I_im, Complex.normSq_apply,
          Complex.ofReal_re, Complex.ofReal_im]
  have h_omega_im : omega3.im = Real.sin (2 * Real.pi / 3) := by
    unfold omega3
    rw [Complex.exp_im, h_re0, h_im_val, Real.exp_zero, one_mul]
  rw [h_omega_im] at h_eq_im
  have h_sin_pos : 0 < Real.sin (2 * Real.pi / 3) := by
    apply Real.sin_pos_of_pos_of_lt_pi
    · linarith [Real.pi_pos]
    · linarith [Real.pi_pos]
  linarith

/-- **`Sp_C 3 ω = 0`**. -/
theorem Sp_C_p3_omega3_zero : Sp_C 3 omega3 = 0 :=
  (Sp_C_p3_zero_iff_cube_root _).mpr ⟨omega3_cube, omega3_ne_one⟩

/-- **`Sp_C 3 (ω²) = 0`**. -/
theorem Sp_C_p3_omega3_sq_zero : Sp_C 3 (omega3 ^ 2) = 0 := by
  apply (Sp_C_p3_zero_iff_cube_root _).mpr
  refine ⟨?_, ?_⟩
  · rw [show (omega3 ^ 2) ^ 3 = (omega3 ^ 3) ^ 2 from by ring, omega3_cube]
    ring
  · intro h_eq_one
    -- ω² = 1, ω³ = 1 ⟹ ω = ω³ / ω² = 1 (contradiction).
    have h_omega_eq_one : omega3 = 1 := by
      have h_cube := omega3_cube
      have : omega3 * omega3 ^ 2 = omega3 ^ 3 := by ring
      rw [h_eq_one, mul_one] at this
      rw [this, h_cube]
    exact omega3_ne_one h_omega_eq_one

/-! ============================================================
  §61.3  Zéros de `Sp_C 3` ⊆ `{ω, ω²}`
============================================================ -/

/-- **`{z : Sp_C 3 z = 0} ⊆ {ω, ω²}`.** -/
theorem Sp_C_p3_zeros_subset :
    {z : ℂ | Sp_C 3 z = 0} ⊆ ({omega3, omega3 ^ 2} : Set ℂ) := by
  intro z hz
  rw [Set.mem_setOf_eq, Sp_C_p3_zero_iff_cube_root] at hz
  obtain ⟨h_cube, h_ne⟩ := hz
  -- Factorisation complète : z³ − 1 = (z − 1)(z − ω)(z − ω²).
  have h_sum_omega_zero : (1 : ℂ) + omega3 + omega3 ^ 2 = 0 := by
    have h_sp := Sp_C_p3_omega3_zero
    rw [Sp_C_p3] at h_sp
    linear_combination h_sp
  have h_prod_omega : omega3 * omega3 ^ 2 = 1 := by
    rw [show omega3 * omega3 ^ 2 = omega3 ^ 3 from by ring, omega3_cube]
  have h_sum_one : omega3 + omega3 ^ 2 + 1 = 0 := by
    linear_combination h_sum_omega_zero
  have h_full_factor :
      z ^ 3 - 1 = (z - 1) * (z - omega3) * (z - omega3 ^ 2) := by
    have h_expand :
        (z - 1) * (z - omega3) * (z - omega3 ^ 2)
      = z ^ 3 - ((1 : ℂ) + omega3 + omega3 ^ 2) * z ^ 2
          + (omega3 + omega3 ^ 2 + omega3 * omega3 ^ 2) * z
          - omega3 * omega3 ^ 2 := by ring
    rw [h_expand, h_sum_omega_zero, h_prod_omega]
    have h_mid : omega3 + omega3 ^ 2 + 1 = 0 := h_sum_one
    linear_combination -z * h_mid
  have h_prod_zero :
      (z - 1) * (z - omega3) * (z - omega3 ^ 2) = 0 := by
    rw [← h_full_factor]; linear_combination h_cube
  have h_sub_one_ne : z - 1 ≠ 0 := sub_ne_zero.mpr h_ne
  -- On factorise : (A·B·C) = 0 et A ≠ 0 ⟹ B·C = 0.
  have h_BC_zero : (z - omega3) * (z - omega3 ^ 2) = 0 := by
    have := h_prod_zero
    rcases mul_eq_zero.mp this with hAB | hC
    · rcases mul_eq_zero.mp hAB with hA | hB
      · exact absurd hA h_sub_one_ne
      · rw [hB]; ring
    · rw [hC]; ring
  rcases mul_eq_zero.mp h_BC_zero with h1 | h2
  · left; exact sub_eq_zero.mp h1
  · right; exact sub_eq_zero.mp h2

/-! ============================================================
  §61.4  Finitude
============================================================ -/

/-- **★ L'ensemble des zéros de `Sp_C 3` est fini.** -/
theorem Sp_C_p3_zeros_finite : ({z : ℂ | Sp_C 3 z = 0}).Finite :=
  (Set.toFinite ({omega3, omega3 ^ 2} : Set ℂ)).subset Sp_C_p3_zeros_subset

/-! ============================================================
  §61.5  Mesure nulle
============================================================ -/

/-- **★★ Les zéros de `Sp_C 3` sont de mesure de Lebesgue nulle.**

    Conséquence directe de la finitude + `Set.Countable.measure_zero`.
    Brique inconditionnelle vers `ComplexJuliaP3MeasureZero 2 α` :
    la partie mesure-théorique du "bad set" direct est acquise. -/
theorem Sp_C_p3_zeros_measure_zero :
    volume ({z : ℂ | Sp_C 3 z = 0}) = 0 :=
  Sp_C_p3_zeros_finite.countable.measure_zero _

end Pandrosion
