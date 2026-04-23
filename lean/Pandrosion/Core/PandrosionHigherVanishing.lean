/-
  Universitas Pandrosion — §108. **★★★★★★★★★★★★ HIGHER VANISHING
  RELATIONS pour la zeta Pandrosion (PDF Théorème 2.2).**

  Pour `P(z) = z^p − 1/x` avec racines cyclotomiques `α_k`,
  l'article "The Pandrosion Zeta Function" (Théorème 2.2) énonce
  les `d − 1` *relations de vanishing supérieures* :

    `Σ_{k ∈ Fin p} α_k^m / P'(α_k) = {0 si 0 ≤ m ≤ p−2; 1 si m = p−1}`

  Ces relations encodent le **système de biorthogonalité** :

    `(α_1^m, …, α_p^m) ⊥ (1/P'(α_1), …, 1/P'(α_p))` pour `m ≤ p−2`,

  avec normalisation à `m = p−1` qui fixe l'échelle.

  **Vérification Python** (`/tmp/pandrosion_higher_vanishing_verify.py`) :
  tous les cas `p ∈ {2,…,7}`, `x ∈ {1,2,5,10}`, `m ∈ {0,…,p−1}`
  testés et confirmés à précision flottante.

  **Stratégie de preuve** :
    1. `1/P'(α_k) = x · α_k / p` ⟹ `α_k^m / P'(α_k) = x · α_k^{m+1} / p`.
    2. `Σ_k α_k^{m+1} = α^{m+1} · Σ_k ω^{k(m+1)}`.
    3. Pour `1 ≤ m+1 ≤ p−1` : `Σ_k ω^{k(m+1)} = 0` (geom_sum + ω^r ≠ 1).
    4. Pour `m+1 = p` : `α_k^{p-1} / P'(α_k) = 1/p`, somme `= 1`.

  Contents.

    §108.1  `omega_pow_r_ne_one` — `ω^r ≠ 1` pour `0 < r < p`.
    §108.2  `sum_omega_pow_kr_eq_zero` — `Σ ω^{kr} = 0` (`0 < r < p`).
    §108.3  `sum_cycAnchor_pow_eq_zero` — `Σ α_k^{m+1} = 0`
            pour `0 ≤ m ≤ p-2`.
    §108.4  ★ `pandrosion_higher_vanishing` — vanishing `m ≤ p-2`.
    §108.5  ★ `pandrosion_normalization` — `Σ α_k^{p-1}/P'(α_k) = 1`.
-/

import Pandrosion.Core.PandrosionSpectralDeterminant

namespace Pandrosion

open Complex Finset

/-! ============================================================
  §108.1  `ω^r ≠ 1` pour `0 < r < p`
============================================================ -/

section OmegaPowNonOne

/-- **`ω^r ≠ 1` pour `0 < r < p`** (via `IsPrimitiveRoot`). -/
theorem omega_pow_r_ne_one (p r : ℕ) (hp : 2 ≤ p) (hr_pos : 1 ≤ r) (hr_lt : r < p) :
    (Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ)))^r ≠ 1 := by
  intro h
  have h_primitive : IsPrimitiveRoot
      (Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))) p :=
    Complex.isPrimitiveRoot_exp p (by omega)
  have h_dvd : p ∣ r := (h_primitive.pow_eq_one_iff_dvd r).mp h
  have h_le : p ≤ r := Nat.le_of_dvd hr_pos h_dvd
  omega

end OmegaPowNonOne

/-! ============================================================
  §108.2  Σ_{k ∈ range p} ω^{k·r} = 0 pour `0 < r < p`
============================================================ -/

section SumOmegaPowZero

/-- **`Σ_{k ∈ range p} (ω^r)^k = 0` pour `0 < r < p`** (geom_sum). -/
theorem sum_omega_pow_kr_eq_zero
    (p r : ℕ) (hp : 2 ≤ p) (hr_pos : 1 ≤ r) (hr_lt : r < p) :
    (Finset.range p).sum
      (fun k : ℕ => (Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ)))^(k * r))
      = 0 := by
  set ω := Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
  have hζ_ne_one : ω^r ≠ 1 := omega_pow_r_ne_one p r hp hr_pos hr_lt
  have hω_p_eq_one : ω^p = 1 := zeta_omega_pow_p_eq_one p (by omega)
  have hζ_p_eq_one : (ω^r)^p = 1 := by
    rw [← pow_mul, mul_comm, pow_mul, hω_p_eq_one, one_pow]
  -- Rewrite ω^(k·r) = (ω^r)^k.
  have h_rewrite : (Finset.range p).sum (fun k : ℕ => ω^(k * r))
                 = (Finset.range p).sum (fun k : ℕ => (ω^r)^k) := by
    apply Finset.sum_congr rfl
    intro k _
    rw [show k * r = r * k from Nat.mul_comm k r, pow_mul]
  rw [h_rewrite, geom_sum_eq hζ_ne_one, hζ_p_eq_one]
  ring

end SumOmegaPowZero

/-! ============================================================
  §108.3  Σ_{k ∈ Fin p} α_k^{m+1} = 0 pour `0 ≤ m ≤ p-2`
============================================================ -/

section SumCycAnchorPowZero

/-- **`Σ_{k ∈ Fin p} (cycAnchor α p k)^{m+1} = 0` pour `0 ≤ m ≤ p-2`.** -/
theorem sum_cycAnchor_pow_eq_zero
    (α : ℂ) (p : ℕ) (hp : 2 ≤ p) (m : ℕ) (hm : m ≤ p - 2) :
    (Finset.univ : Finset (Fin p)).sum
      (fun k : Fin p => (cycAnchor α p k)^(m + 1)) = 0 := by
  set ω := Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
  have h_rewrite : ∀ k : Fin p, (cycAnchor α p k)^(m + 1)
                              = α^(m + 1) * ω^((k : ℕ) * (m + 1)) := by
    intro k
    rw [cycAnchor_eq_alpha_omega_pow α p k]
    rw [mul_pow, ← pow_mul]
  rw [Finset.sum_congr rfl (fun k _ => h_rewrite k)]
  rw [← Finset.mul_sum]
  have hr_pos : 1 ≤ m + 1 := by omega
  have hr_lt : m + 1 < p := by omega
  -- Reindex Fin p → range p.
  have h_reindex : (Finset.univ : Finset (Fin p)).sum
                     (fun k : Fin p => ω^((k : ℕ) * (m + 1)))
                 = (Finset.range p).sum (fun k : ℕ => ω^(k * (m + 1))) := by
    rw [← Fin.sum_univ_eq_sum_range (fun k : ℕ => ω^(k * (m + 1)))]
  rw [h_reindex]
  rw [sum_omega_pow_kr_eq_zero p (m + 1) hp hr_pos hr_lt]
  ring

end SumCycAnchorPowZero

/-! ============================================================
  §108.4  ★ Higher Vanishing : Σ α_k^m / P'(α_k) = 0 pour `m ≤ p-2`
============================================================ -/

section HigherVanishing

/-- **★★★★★★★★★★★★ HIGHER VANISHING (PDF Théorème 2.2).**

    Pour `0 ≤ m ≤ p − 2` :

    `Σ_{k ∈ Fin p} (cycAnchor α p k)^m / P'(cycAnchor α p k) = 0`

    Encode l'orthogonalité : `(α_k^m)_k ⊥ (1/P'(α_k))_k`. -/
theorem pandrosion_higher_vanishing
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 2 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x) (m : ℕ) (hm : m ≤ p - 2) :
    (Finset.univ : Finset (Fin p)).sum
      (fun k : Fin p => (cycAnchor α p k)^m / pandrosionPDeriv p (cycAnchor α p k))
      = 0 := by
  have hp_ge_one : 1 ≤ p := by omega
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  -- Local non-zero helper.
  have hαk_ne : ∀ k : Fin p, cycAnchor α p k ≠ 0 := by
    intro k h_zero
    have h_pow : (cycAnchor α p k) ^ p = 0 := by
      rw [h_zero]; exact zero_pow (by omega)
    have h_eq : (cycAnchor α p k) ^ p = 1 / x :=
      (cycAnchor_pow α hp_ge_one k).trans hα_pow
    rw [h_eq] at h_pow
    exact one_div_ne_zero hx h_pow
  -- α_k^m / P'(α_k) = x · α_k^(m+1) / p
  have h_rewrite : ∀ k : Fin p,
      (cycAnchor α p k)^m / pandrosionPDeriv p (cycAnchor α p k)
        = x * (cycAnchor α p k)^(m + 1) / (p : ℂ) := by
    intro k
    have h_pow_eq : (cycAnchor α p k)^p = 1 / x :=
      (cycAnchor_pow α hp_ge_one k).trans hα_pow
    have hαk : cycAnchor α p k ≠ 0 := hαk_ne k
    -- α_k^(p-1) = (1/x) / α_k via α_k^p = α_k^(p-1) · α_k
    have h_pow_pminus1 : (cycAnchor α p k)^(p - 1) = (1 / x) / cycAnchor α p k := by
      have h_succ : (cycAnchor α p k)^(p - 1 + 1)
                  = (cycAnchor α p k)^(p - 1) * cycAnchor α p k := pow_succ _ _
      rw [show p - 1 + 1 = p from Nat.sub_add_cancel hp_ge_one] at h_succ
      rw [h_pow_eq] at h_succ
      rw [h_succ]
      field_simp
    unfold pandrosionPDeriv
    rw [h_pow_pminus1]
    field_simp
    ring
  rw [Finset.sum_congr rfl (fun k _ => h_rewrite k)]
  -- Σ x · α_k^(m+1) / p = (x/p) · Σ α_k^(m+1) = (x/p) · 0 = 0
  rw [show (fun k : Fin p => x * (cycAnchor α p k)^(m + 1) / (p : ℂ))
       = (fun k : Fin p => (x / (p : ℂ)) * (cycAnchor α p k)^(m + 1)) from by
    funext k; ring]
  rw [← Finset.mul_sum]
  rw [sum_cycAnchor_pow_eq_zero α p hp m hm]
  ring

end HigherVanishing

/-! ============================================================
  §108.5  ★ Normalization : Σ α_k^{p-1} / P'(α_k) = 1
============================================================ -/

section Normalization

/-- **★ Normalization à `m = p − 1` (PDF Thm 2.2 case `m = d − 1`).**

    `Σ_{k ∈ Fin p} (cycAnchor α p k)^{p-1} / P'(cycAnchor α p k) = 1`.

    Preuve directe : `α_k^{p-1} / (p · α_k^{p-1}) = 1/p`, somme = `1`. -/
theorem pandrosion_normalization
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x) :
    (Finset.univ : Finset (Fin p)).sum
      (fun k : Fin p => (cycAnchor α p k)^(p - 1) / pandrosionPDeriv p (cycAnchor α p k))
      = 1 := by
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hαk_ne : ∀ k : Fin p, cycAnchor α p k ≠ 0 := by
    intro k h_zero
    have h_pow : (cycAnchor α p k) ^ p = 0 := by
      rw [h_zero]; exact zero_pow (by omega)
    have h_eq : (cycAnchor α p k) ^ p = 1 / x :=
      (cycAnchor_pow α hp k).trans hα_pow
    rw [h_eq] at h_pow
    exact one_div_ne_zero hx h_pow
  have h_rewrite : ∀ k : Fin p,
      (cycAnchor α p k)^(p - 1) / pandrosionPDeriv p (cycAnchor α p k)
        = (1 : ℂ) / (p : ℂ) := by
    intro k
    unfold pandrosionPDeriv
    have hαk_pow_ne : (cycAnchor α p k)^(p - 1) ≠ 0 := pow_ne_zero _ (hαk_ne k)
    rw [show (p : ℂ) * (cycAnchor α p k)^(p - 1)
         = (cycAnchor α p k)^(p - 1) * (p : ℂ) from by ring]
    field_simp
  rw [Finset.sum_congr rfl (fun k _ => h_rewrite k)]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
  field_simp

end Normalization

end Pandrosion
