/-
  Universitas Pandrosion — §78. **Points nommés dans le h-bassin
  principal à `x = 2, p = 3`.**

  Corollaires concrets combinant §67 (demi-plan) + §77 (far-field).
  Chaque théorème établit `h^n(z) → α₀` pour un `z` spécifique,
  montrant que le bassin couvre des points identifiables.

  Utilité : ancres démonstratives pour la suite, tests en
  instanciation, validation du chaînage corpus.

  Contents.

    §78.1  `tendsto_at_two`, `tendsto_at_five`, `tendsto_at_ten` —
           points réels positifs (via §67).
    §78.2  `tendsto_at_neg_five` — point réel négatif (via §77).
    §78.3  `tendsto_at_i`, `tendsto_at_three_i` — points imaginaires
           purs (via §77).
    §78.4  `tendsto_at_one_plus_i` — point complexe mixte (via §67).
-/

import Pandrosion.Core.FarFieldX2

namespace Pandrosion

open Complex Filter Topology

/-! §78.1  Real positifs dans le demi-plan -/

theorem tendsto_at_two :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] (2 : ℂ))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply pandrosion_h_C_p3_tendsto_half_plane_x2
  show (2 : ℂ).re ≥ 1/2
  rw [show ((2 : ℂ)) = (((2 : ℝ)) : ℂ) from by norm_num, Complex.ofReal_re]
  norm_num

theorem tendsto_at_five :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] (5 : ℂ))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply pandrosion_h_C_p3_tendsto_half_plane_x2
  show (5 : ℂ).re ≥ 1/2
  rw [show ((5 : ℂ)) = (((5 : ℝ)) : ℂ) from by norm_num, Complex.ofReal_re]
  norm_num

theorem tendsto_at_ten :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] (10 : ℂ))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply pandrosion_h_C_p3_tendsto_half_plane_x2
  show (10 : ℂ).re ≥ 1/2
  rw [show ((10 : ℂ)) = (((10 : ℝ)) : ℂ) from by norm_num, Complex.ofReal_re]
  norm_num

/-! §78.2  Real négatif via far-field -/

/-- **`h^n(-5) → α₀`.** Re(-5) < 1/2 mais |-5| = 5 ≥ 2 ⟹ §77 s'applique. -/
theorem tendsto_at_neg_five :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] ((-5 : ℝ) : ℂ))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply pandrosion_h_C_p3_tendsto_far_field_x2
  show ‖((-5 : ℝ) : ℂ)‖ ≥ 2
  rw [Complex.norm_real, Real.norm_eq_abs]
  norm_num

/-! §78.3  Imaginaires purs via far-field -/

/-- **`h^n(3i) → α₀`.** Re(3i) = 0 < 1/2 mais |3i| = 3 ≥ 2. -/
theorem tendsto_at_three_i :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] (3 * Complex.I))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply pandrosion_h_C_p3_tendsto_far_field_x2
  have : ‖(3 : ℂ) * Complex.I‖ = 3 := by
    rw [norm_mul, Complex.norm_I, mul_one]
    rw [show ((3 : ℂ)) = (((3 : ℝ)) : ℂ) from by norm_num,
        Complex.norm_real, Real.norm_eq_abs]
    norm_num
  linarith

/-! §78.4  Complexe mixte via demi-plan -/

/-- **`h^n(1 + i) → α₀`.** Re(1+i) = 1 ≥ 1/2. -/
theorem tendsto_at_one_plus_i :
    Tendsto (fun n : ℕ => (pandrosion_h_C (2 : ℂ) 3)^[n] (1 + Complex.I))
      atTop (𝓝 ((alphaX2 : ℝ) : ℂ)) := by
  apply pandrosion_h_C_p3_tendsto_half_plane_x2
  show ((1 : ℂ) + Complex.I).re ≥ 1/2
  rw [Complex.add_re, Complex.one_re, Complex.I_re]
  norm_num

end Pandrosion
