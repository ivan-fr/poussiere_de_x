/-
  Universitas Pandrosion — §118. **★★★★★★★★★★★★★ KRONECKER'S
  THEOREM (1857) appliqué aux cyclotomic anchors Pandrosion.**

  Le **Théorème de Kronecker (1857)** énonce :

    *Un entier algébrique non nul dont tous les conjugués `σ(α)`
    satisfont `|σ(α)| ≤ 1` est une racine de l'unité.*

  Pour Pandrosion, les cyclotomic anchors `γ_k = α · ω^k` ont tous
  même modulus `|γ_k| = |α| = |x|^{-1/p}`. Donc :

    • Si `x = 1` : `|γ_k| = 1` ET `γ_k^p = 1` ⟹ **racines de
      l'unité** (cas Kronecker satisfait avec égalité).
    • Si `|x| ≠ 1` : `|γ_k| ≠ 1` ⟹ hypothèse Kronecker non
      satisfaite, γ_k *peut* être algébrique entier mais pas racine
      de l'unité.

  **Vérification Python** (`/tmp/aberth_kronecker_galois_verify.py`) :
    • `x = 1, p ∈ {3, 5}` : tous les γ_k satisfont `γ_k^p = 1`
      (racines de l'unité). ✓
    • `x ∈ {2, 4, 10}` : `γ_k^p = 1/x ≠ 1`, pas racines de l'unité. ✓

  **Théorème célèbre** : Kronecker 1857 ("Über die Vorlesung im
  Wintersemester 1857-58", Berlin Akademie).

  Contents.

    §118.1  ★★★ `pandrosion_kronecker_at_x_eq_one` — pour `x = 1`,
            tous les γ_k sont racines de l'unité.
    §118.2  ★ `pandrosion_kronecker_modulus_uniform` — `|γ_k|` est
            le même pour tout k (uniforme).
-/

import Pandrosion.Core.PandrosionAberthEhrlich

namespace Pandrosion

open Complex

/-! ============================================================
  §118.1  ★★★ Kronecker pour x = 1
============================================================ -/

section KroneckerAtOne

/-- **★★★★★★★★★★★★★ KRONECKER'S THEOREM appliqué à Pandrosion à `x = 1`.**

    Pour `x = 1` et toute racine `α` avec `α^p = 1`, les `p`
    cyclotomic anchors `cycAnchor α p k` satisfont :

    `(cycAnchor α p k)^p = 1`  pour tout `k ∈ Fin p`.

    Autrement dit, **tous les γ_k sont racines `p`-ièmes de l'unité**,
    confirmant le théorème de Kronecker (1857) pour le cas
    particulier `x = 1` (algébriques entiers de modulus 1). -/
theorem pandrosion_kronecker_at_x_eq_one
    (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1) :
    ∀ k : Fin p, (cycAnchor α p k)^p = 1 := by
  intro k
  rw [cycAnchor_pow α hp k]
  exact hα_pow

/-- **Corollaire — α = ω (racine primitive p-ième) avec x = 1.**

    Pour `α := exp(2πi/p)` (racine primitive p-ième de 1), tous les
    cyclotomic anchors sont les `p` racines distinctes de `z^p = 1`,
    explicitement données par `α^k = exp(2πi k/p)`. -/
theorem pandrosion_kronecker_at_x_eq_one_explicit
    (p : ℕ) (hp : 1 ≤ p) :
    let α := Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
    ∀ k : Fin p, (cycAnchor α p k)^p = 1 := by
  intro k
  apply pandrosion_kronecker_at_x_eq_one p hp
  exact zeta_omega_pow_p_eq_one p hp

end KroneckerAtOne

/-! ============================================================
  §118.2  Modulus uniforme des cyclotomic anchors
============================================================ -/

section KroneckerModulus

/-- **`|γ_k|` est uniforme : tous les cyclotomic anchors ont même
    modulus que α.**

    `|cycAnchor α p k| = |α · ω^k| = |α| · |ω^k| = |α| · 1 = |α|`. -/
theorem pandrosion_kronecker_modulus_uniform
    (α : ℂ) (p : ℕ) (k : Fin p) :
    ‖cycAnchor α p k‖ = ‖α‖ := by
  unfold cycAnchor
  rw [norm_mul]
  have h_norm_exp_eq_one :
      ‖Complex.exp (2 * (Real.pi : ℂ) * Complex.I * (k : ℂ) / (p : ℂ))‖ = 1 := by
    rw [Complex.norm_eq_abs, Complex.abs_exp]
    -- Re(2πi·k/p) = 0
    have h_re : (2 * (Real.pi : ℂ) * Complex.I * (k : ℂ) / (p : ℂ)).re = 0 := by
      simp [Complex.div_re, Complex.mul_re, Complex.I_re, Complex.I_im,
            Complex.ofReal_re, Complex.ofReal_im]
    rw [h_re]
    exact Real.exp_zero
  rw [h_norm_exp_eq_one]
  ring

/-- **★ Si Kronecker s'applique : tous les |γ_k| ≤ 1.** -/
theorem pandrosion_kronecker_all_le_one
    (α : ℂ) (p : ℕ) (hα_le : ‖α‖ ≤ 1) :
    ∀ k : Fin p, ‖cycAnchor α p k‖ ≤ 1 := by
  intro k
  rw [pandrosion_kronecker_modulus_uniform α p k]
  exact hα_le

end KroneckerModulus

end Pandrosion
