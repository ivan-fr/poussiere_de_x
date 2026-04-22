/-
  Universitas Pandrosion — §72. **Forme fermée rationnelle de σ à
  `x = 2, p = 3`** : définitions et identités algébriques clefs.

  Par calcul symbolique (sympy validé), l'itération Steffensen σ à
  `x = 2, p = 3` s'exprime comme :

      `σ(z) = (20z⁶ + 56z⁵ + 92z⁴ + 90z³ + 62z² + 26z + 7) /
              (4·(6z⁶ + 17z⁵ + 29z⁴ + 29z³ + 20z² + 8z + 2))`.

  Le dénominateur factorise :
      `6z⁶+17z⁵+29z⁴+29z³+20z²+8z+2 = (1+z+z²)·(6z⁴+11z³+12z²+6z+2)`.

  Et le dénominateur de Steffensen (`h(h(z)) − 2h(z) + z`) s'exprime :
      `denom_C = (2z³−1)·(6z⁴+11z³+12z²+6z+2) /
                 ((1+z+z²)·(12z⁴+24z³+30z²+18z+7))`.

  La chaîne algébrique complète (équivalence `steffensen_step_C 2 3 z
  = sigma_x2_closed z`) repose sur une `field_simp` de degré 18+, qui
  dépasse la tolérance build (>50 s). On garde ici les **identités
  polynomiales élémentaires** prouvées via `ring`, qui capturent
  l'algèbre sans passer par field_simp.

  Contents.

    §72.1  `sigma_x2_closed` — forme rationnelle explicite.
    §72.2  `sigma_x2_closed_at_alpha` — `σ_closed(α₀) = α₀`.
    §72.3  `sigma_x2_closed_denom_factor` — factorisation clef.
    §72.4  `sigma_x2_closed_pivot_identity` — identité polynomiale
           pivot `z·4·BigPoly − (2z³−1)·P₄ = N`.
    §72.5  `pandrosion_h_C_p3_x2_sub_id` — `h(z)−z = (1−2z³)/(2·Sp)`.
    §72.6  `Sp_of_h_p3_x2` — forme fermée de `1+h+h²`.
-/

import Pandrosion.Core.PrincipalDominanceP3X2

namespace Pandrosion

open Complex

/-! §72.1  Forme rationnelle -/

/-- **Forme fermée de σ à `x = 2, p = 3`** (validée par sympy). -/
noncomputable def sigma_x2_closed (z : ℂ) : ℂ :=
  (20*z^6 + 56*z^5 + 92*z^4 + 90*z^3 + 62*z^2 + 26*z + 7) /
    (4 * (6*z^6 + 17*z^5 + 29*z^4 + 29*z^3 + 20*z^2 + 8*z + 2))

/-! §72.2  `sigma_x2_closed α₀ = α₀` -/

theorem sigma_x2_closed_at_alpha :
    sigma_x2_closed ((alphaX2 : ℝ) : ℂ) = ((alphaX2 : ℝ) : ℂ) := by
  unfold sigma_x2_closed
  have h_cube : ((alphaX2 : ℝ) : ℂ) ^ 3 = 1/2 := alphaX2_cube_C
  rw [div_eq_iff]
  · ring_nf
    linear_combination (-24 * ((alphaX2 : ℝ) : ℂ)^4
                        - 48 * ((alphaX2 : ℝ) : ℂ)^3
                        - 60 * ((alphaX2 : ℝ) : ℂ)^2
                        - 36 * ((alphaX2 : ℝ) : ℂ)
                        - 14) * h_cube
  · have h_denom_val :
        4 * (6 * ((alphaX2 : ℝ) : ℂ)^6 + 17 * ((alphaX2 : ℝ) : ℂ)^5
             + 29 * ((alphaX2 : ℝ) : ℂ)^4 + 29 * ((alphaX2 : ℝ) : ℂ)^3
             + 20 * ((alphaX2 : ℝ) : ℂ)^2 + 8 * ((alphaX2 : ℝ) : ℂ) + 2)
        = 72 + 114 * ((alphaX2 : ℝ) : ℂ)^2 + 90 * ((alphaX2 : ℝ) : ℂ) := by
      have h6 : ((alphaX2 : ℝ) : ℂ)^6 = (((alphaX2 : ℝ) : ℂ)^3)^2 := by ring
      have h5 : ((alphaX2 : ℝ) : ℂ)^5 = ((alphaX2 : ℝ) : ℂ)^3 * ((alphaX2 : ℝ) : ℂ)^2 := by ring
      have h4 : ((alphaX2 : ℝ) : ℂ)^4 = ((alphaX2 : ℝ) : ℂ)^3 * ((alphaX2 : ℝ) : ℂ) := by ring
      rw [h6, h5, h4, h_cube]
      ring
    rw [h_denom_val]
    have h_real : (72 : ℂ) + 114 * ((alphaX2 : ℝ) : ℂ)^2 + 90 * ((alphaX2 : ℝ) : ℂ)
                = (((72 + 114 * alphaX2^2 + 90 * alphaX2) : ℝ) : ℂ) := by
      push_cast; ring
    rw [h_real, Ne, Complex.ofReal_eq_zero]
    have h_pos : 0 < 72 + 114 * alphaX2^2 + 90 * alphaX2 := by
      have := sq_nonneg alphaX2
      linarith [alphaX2_pos]
    linarith

/-! §72.3  Factorisation du dénominateur -/

theorem sigma_x2_closed_denom_factor (z : ℂ) :
    (6 : ℂ)*z^6 + 17*z^5 + 29*z^4 + 29*z^3 + 20*z^2 + 8*z + 2
      = (1 + z + z^2) * (6*z^4 + 11*z^3 + 12*z^2 + 6*z + 2) := by ring

/-! §72.4  Identité-pivot polynomiale -/

/-- **Identité clef** : `z · 4·BigPoly − (2z³−1)·P₄ = N`.

    C'est la forme polynomiale pure de la conclusion algébrique
    `sigma_x2_closed z = z − (2z³−1)·P₄ / (4·BigPoly)`, transposable
    à l'analyse σ de demi-plan sans passer par `field_simp` coûteux. -/
theorem sigma_x2_closed_pivot_identity (z : ℂ) :
    z * (4 * (6*z^6 + 17*z^5 + 29*z^4 + 29*z^3 + 20*z^2 + 8*z + 2))
      - (2*z^3 - 1) * (12*z^4 + 24*z^3 + 30*z^2 + 18*z + 7)
      = 20*z^6 + 56*z^5 + 92*z^4 + 90*z^3 + 62*z^2 + 26*z + 7 := by ring

/-! §72.5  Forme fermée de h(z) − z -/

/-- `h(z) − z = (1 − 2z³) / (2·Sp(z))` à x=2, p=3, si `Sp(z) ≠ 0`. -/
theorem pandrosion_h_C_p3_x2_sub_id (z : ℂ) (h_Sp : 1 + z + z^2 ≠ 0) :
    pandrosion_h_C (2 : ℂ) 3 z - z = (1 - 2*z^3) / (2 * (1 + z + z^2)) := by
  unfold pandrosion_h_C
  rw [Sp_C_p3]
  have h_two_Sp_ne : (2 : ℂ) * (1 + z + z^2) ≠ 0 :=
    mul_ne_zero (by norm_num) h_Sp
  field_simp
  ring

/-! §72.6  Forme fermée de 1 + h(z) + h(z)² -/

/-- `1 + h(z) + h(z)² = (12z⁴+24z³+30z²+18z+7) / (4·(1+z+z²)²)` à x=2, p=3. -/
theorem Sp_of_h_p3_x2 (z : ℂ) (h_Sp : 1 + z + z^2 ≠ 0) :
    1 + pandrosion_h_C (2 : ℂ) 3 z + (pandrosion_h_C (2 : ℂ) 3 z)^2
      = (12*z^4 + 24*z^3 + 30*z^2 + 18*z + 7) / (4 * (1 + z + z^2)^2) := by
  unfold pandrosion_h_C
  rw [Sp_C_p3]
  have h_2S_ne : (2 : ℂ) * (1 + z + z^2) ≠ 0 :=
    mul_ne_zero (by norm_num) h_Sp
  field_simp
  ring

end Pandrosion
