/-
  Universitas Pandrosion — §80. **Borne polynomiale `‖Sp_C 4 z‖²` à p=4.**

  Généralisation partielle de §54 (p=3) à p=4. Cible originale du
  ROADMAP : "même structure algébrique, degrés bumped". **Découverte**
  empirique en cours d'écriture : la borne `‖Sp_C 4 z‖² ≥
  (1+Re z+Re z²+Re z³)²` ne tient PAS sur tout `Re z ≥ 1/2` à p=4.

  **Contre-exemple p=4** : à `z = 1/2 + i`, `‖Sp_C 4 z‖² ≈ 3.453` mais
  `(1+1/2+1/4+1/8)² = (15/8)² ≈ 3.516`. Échec ~ 0.06.

  **Threshold corrigé** : la borne tient pour `Re z ≥ 1` (où les
  coefficients de l'identité polynomiale en `b² = (Im z)²` sont tous
  ≥ 0). Pour `1/2 ≤ Re z < 1`, le coefficient de b² peut être négatif,
  compensé par b⁴, b⁶ pour |b| grand mais pas pour |b| petit.

  Contents.

    §80.1  `Sp_C_p4` — closed form `Sp_C 4 z = 1 + z + z² + z³`.
    §80.2  `Sp_C_p4_norm_sq_identity` — identité polynomiale exacte
           pour `‖Sp_C 4 z‖²`.
    §80.3  `Sp_C_p4_norm_sq_lower_bound_at_one` — borne
           `‖Sp_C 4 z‖² ≥ (1 + a + a² + a³)²` pour `Re z ≥ 1`.
-/

import Pandrosion.Core.HalfPlaneContractionP3

namespace Pandrosion

open Complex

/-! §80.1  Closed form `Sp_C 4` -/

theorem Sp_C_p4 (z : ℂ) : Sp_C 4 z = 1 + z + z ^ 2 + z ^ 3 := by
  unfold Sp_C
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_succ,
      Finset.sum_range_succ, Finset.sum_range_zero]
  ring

/-! §80.2  Polynomial identity for `‖Sp_C 4 z‖²` -/

/-- **Identité algébrique polynomiale pour `‖Sp_C 4 z‖²`**.

    Avec `a = Re z, b = Im z, A = 1+a+a²+a³, P = 1+3a, Q = 1+2a+3a²` :

      `‖Sp_C 4 z‖² − A² = b²·(Q² − 2·A·P) + b⁴·(P² − 2·Q) + b⁶`. -/
theorem Sp_C_p4_norm_sq_identity (z : ℂ) :
    ‖Sp_C 4 z‖ ^ 2 - (1 + z.re + z.re ^ 2 + z.re ^ 3) ^ 2
      = z.im ^ 2 * ((1 + 2*z.re + 3*z.re^2)^2
                      - 2 * (1+z.re+z.re^2+z.re^3) * (1+3*z.re))
        + z.im ^ 4 * ((1+3*z.re)^2 - 2*(1+2*z.re+3*z.re^2))
        + z.im ^ 6 := by
  have h_Sp : Sp_C 4 z = 1 + z + z ^ 2 + z ^ 3 := Sp_C_p4 z
  have h_norm_sq : ‖Sp_C 4 z‖ ^ 2 = (Sp_C 4 z).re ^ 2 + (Sp_C 4 z).im ^ 2 := by
    rw [Complex.norm_eq_abs, Complex.sq_abs, Complex.normSq_apply]; ring
  have h_z2_re : (z ^ 2).re = z.re ^ 2 - z.im ^ 2 := by
    rw [sq]; simp [Complex.mul_re]; ring
  have h_z2_im : (z ^ 2).im = 2 * z.re * z.im := by
    rw [sq]; simp [Complex.mul_im]; ring
  have h_z3_re : (z ^ 3).re = z.re ^ 3 - 3 * z.re * z.im ^ 2 := by
    have h_eq : z ^ 3 = z * z ^ 2 := by ring
    rw [h_eq]
    simp [Complex.mul_re, h_z2_re, h_z2_im]
    ring
  have h_z3_im : (z ^ 3).im = 3 * z.re ^ 2 * z.im - z.im ^ 3 := by
    have h_eq : z ^ 3 = z * z ^ 2 := by ring
    rw [h_eq]
    simp [Complex.mul_im, h_z2_re, h_z2_im]
    ring
  have h_Sp_re : (Sp_C 4 z).re = 1 + z.re + (z.re^2 - z.im^2)
                                   + (z.re^3 - 3*z.re*z.im^2) := by
    rw [h_Sp]
    simp only [Complex.add_re, Complex.one_re, h_z2_re, h_z3_re]
  have h_Sp_im : (Sp_C 4 z).im = z.im + 2*z.re*z.im + (3*z.re^2*z.im - z.im^3) := by
    rw [h_Sp]
    simp only [Complex.add_im, Complex.one_im, h_z2_im, h_z3_im]
    ring
  rw [h_norm_sq, h_Sp_re, h_Sp_im]
  ring

/-! §80.3  Lower bound for Re z ≥ 1 -/

/-- **Borne `‖Sp_C 4 z‖² ≥ (1+a+a²+a³)²` pour `Re z ≥ 1`** (PAS pour
    `Re z ≥ 1/2` — voir contre-exemple à `z = 1/2 + i` dans le header).

    *Preuve.* Sur `a ≥ 1`, les coefficients `(Q² − 2·A·P)` et
    `(P² − 2·Q)` de l'identité §80.2 sont tous deux ≥ 0 (chacun
    ≥ 4), donc la somme `b²·... + b⁴·... + b⁶ ≥ 0`. -/
theorem Sp_C_p4_norm_sq_lower_bound_at_one (z : ℂ) (hre : z.re ≥ 1) :
    ‖Sp_C 4 z‖ ^ 2 ≥ (1 + z.re + z.re ^ 2 + z.re ^ 3) ^ 2 := by
  have h_id := Sp_C_p4_norm_sq_identity z
  set a := z.re
  set b := z.im
  -- Coeff de b² : `(1+2a+3a²)² − 2·(1+a+a²+a³)·(1+3a) = 3a⁴+4a³+2a²−4a−1`.
  -- Identité : 3a⁴+4a³+2a²−4a−1 = 4 + (a−1)·(3a³+7a²+9a+5).
  have h_coeff_b2 :
      (1 + 2*a + 3*a^2)^2 - 2 * (1+a+a^2+a^3) * (1+3*a)
      = 3*a^4 + 4*a^3 + 2*a^2 - 4*a - 1 := by ring
  have h_coeff_b2_nn :
      0 ≤ (1 + 2*a + 3*a^2)^2 - 2 * (1+a+a^2+a^3) * (1+3*a) := by
    rw [h_coeff_b2]
    have h_a_m1_nn : 0 ≤ a - 1 := by linarith
    have h_poly_nn : 0 ≤ 3*a^3 + 7*a^2 + 9*a + 5 := by
      have h_a_nn : 0 ≤ a := by linarith
      linarith [pow_nonneg h_a_nn 3, sq_nonneg a]
    have h_prod_nn : 0 ≤ (a - 1) * (3*a^3 + 7*a^2 + 9*a + 5) :=
      mul_nonneg h_a_m1_nn h_poly_nn
    linarith [show 3*a^4 + 4*a^3 + 2*a^2 - 4*a - 1
                  = 4 + (a - 1) * (3*a^3 + 7*a^2 + 9*a + 5) from by ring]
  -- Coeff de b⁴ : `(1+3a)² − 2·(1+2a+3a²) = 3a²+2a−1`.
  -- Identité : 3a²+2a−1 = 4 + (a−1)·(3a+5).
  have h_coeff_b4 :
      (1+3*a)^2 - 2*(1+2*a+3*a^2) = 3*a^2 + 2*a - 1 := by ring
  have h_coeff_b4_nn :
      0 ≤ (1+3*a)^2 - 2*(1+2*a+3*a^2) := by
    rw [h_coeff_b4]
    have h_a_m1_nn : 0 ≤ a - 1 := by linarith
    have h_poly_nn : 0 ≤ 3*a + 5 := by linarith
    linarith [mul_nonneg h_a_m1_nn h_poly_nn,
              show 3*a^2 + 2*a - 1 = 4 + (a - 1) * (3*a + 5) from by ring]
  -- Combine.
  have h_b2_term_nn : 0 ≤ b^2 * ((1 + 2*a + 3*a^2)^2
                                    - 2 * (1+a+a^2+a^3) * (1+3*a)) :=
    mul_nonneg (sq_nonneg _) h_coeff_b2_nn
  have h_b4_term_nn : 0 ≤ b^4 * ((1+3*a)^2 - 2*(1+2*a+3*a^2)) :=
    mul_nonneg (by positivity) h_coeff_b4_nn
  have h_b6_nn : 0 ≤ b^6 := by positivity
  linarith

end Pandrosion
