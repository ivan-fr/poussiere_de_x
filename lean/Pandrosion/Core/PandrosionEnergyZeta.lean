/-
  Universitas Pandrosion — §111. **★★★★★★★★★★★ PANDROSION ENERGY
  ZETA FUNCTION (PDF §4).**

  Définit la **fonction zeta d'énergie** (PDF §4) :

    `Z_P(s) := Σ_{k} E_P(α_k)^{-s}`

  où `E_P(α_k) := P'(α_k)^2` (self-energy au k-ième zéro).

  Pour Pandrosion `P(z) = z^p − 1/x` :
    `E_P(α_k) = (p · α_k^{p-1})^2`,
    `1/E_P(α_k) = (x · α_k / p)^2 = x²·α_k²/p²`.

  **Vérification Python** (`/tmp/pandrosion_energy_zeta_verify.py`) :
    • `Z_P(0) = p` pour tout `p, x` ✓
    • `Z_P(1) = 0` pour `p ≥ 3` ✓
    • `Z_P(-1) = 0` pour `p ≥ 3` ✓
    • `Z_P(s) = ζ_P(2s)` à `s` entier ✓

  Contents.

    §111.1  `pandrosionEnergy` — `E_P(α_k) := P'(α_k)^2`.
    §111.2  `sum_cycAnchor_pow_eq_zero_general` — généralisation
            de §108 pour tout exposant `e` avec `p ∤ e`.
    §111.3  Valeurs spéciales `Z_P(0)`, `Z_P(1)`, `Z_P(-1)`.
    §111.4  ★★★ `pandrosion_energy_zeta_master` — capstone à 3
            valeurs (PDF §4 Table).
-/

import Pandrosion.Core.PandrosionHigherVanishing

namespace Pandrosion

open Complex Finset

/-! ============================================================
  §111.1  Définition `pandrosionEnergy`
============================================================ -/

section EnergyDefinition

/-- **Self-energy au k-ième zéro** : `E_P(α_k) := P'(α_k)^2`. -/
noncomputable def pandrosionEnergy (p : ℕ) (z : ℂ) : ℂ :=
  (pandrosionPDeriv p z)^2

end EnergyDefinition

/-! ============================================================
  §111.2  Généralisation : Σ α_k^e = 0 pour tout `e` avec `p ∤ e`
============================================================ -/

section SumCycAnchorPowGeneral

/-- **`Σ_{k ∈ Fin p} α_k^e = 0` pour tout exposant `e ≥ 1` avec `p ∤ e`.**

    Généralisation de §108.3 `sum_cycAnchor_pow_eq_zero` (qui couvrait
    le cas `e = m + 1` avec `m + 1 ≤ p - 1`). -/
theorem sum_cycAnchor_pow_eq_zero_general
    (α : ℂ) (p : ℕ) (hp : 2 ≤ p) (e : ℕ) (_he_pos : 1 ≤ e) (he_not_dvd : ¬ p ∣ e) :
    (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => (cycAnchor α p k)^e) = 0 := by
  set ω := Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
  have h_rewrite : ∀ k : Fin p, (cycAnchor α p k)^e = α^e * ω^((k : ℕ) * e) := by
    intro k
    rw [cycAnchor_eq_alpha_omega_pow α p k]
    rw [mul_pow, ← pow_mul]
  rw [Finset.sum_congr rfl (fun k _ => h_rewrite k)]
  rw [← Finset.mul_sum]
  have h_reindex : (Finset.univ : Finset (Fin p)).sum
                     (fun k : Fin p => ω^((k : ℕ) * e))
                 = (Finset.range p).sum (fun k : ℕ => ω^(k * e)) := by
    rw [← Fin.sum_univ_eq_sum_range (fun k : ℕ => ω^(k * e))]
  rw [h_reindex]
  -- ω^e ≠ 1 since p ∤ e.
  have hζ_ne_one : ω^e ≠ 1 := by
    intro h
    have h_primitive : IsPrimitiveRoot ω p := Complex.isPrimitiveRoot_exp p (by omega)
    have h_dvd : p ∣ e := (h_primitive.pow_eq_one_iff_dvd e).mp h
    exact he_not_dvd h_dvd
  have hω_p_eq_one : ω^p = 1 := zeta_omega_pow_p_eq_one p (by omega)
  have hζ_p_eq_one : (ω^e)^p = 1 := by
    rw [← pow_mul, mul_comm, pow_mul, hω_p_eq_one, one_pow]
  -- Rewrite as (ω^e)^k.
  have h_rewrite_inner : (Finset.range p).sum (fun k : ℕ => ω^(k * e))
                       = (Finset.range p).sum (fun k : ℕ => (ω^e)^k) := by
    apply Finset.sum_congr rfl
    intro k _
    rw [show k * e = e * k from Nat.mul_comm k e, pow_mul]
  rw [h_rewrite_inner]
  rw [geom_sum_eq hζ_ne_one, hζ_p_eq_one]
  ring

end SumCycAnchorPowGeneral

/-! ============================================================
  §111.3  Valeurs spéciales `Z_P(0)`, `Z_P(1)`, `Z_P(-1)`
============================================================ -/

section EnergyZetaSpecialValues

/-- **Helper : `cycAnchor α p k ≠ 0`.** -/
private theorem cycAnchor_ne_zero_ez
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x) (k : Fin p) :
    cycAnchor α p k ≠ 0 := by
  intro h_zero
  have h_pow : (cycAnchor α p k) ^ p = 0 := by
    rw [h_zero]; exact zero_pow (by omega)
  have h_eq : (cycAnchor α p k) ^ p = 1 / x :=
    (cycAnchor_pow α hp k).trans hα_pow
  rw [h_eq] at h_pow
  exact one_div_ne_zero hx h_pow

/-- **`Z_P(0) = p`** : valeur triviale (compte des racines). -/
theorem pandrosion_energy_zeta_at_zero (p : ℕ) :
    (Finset.univ : Finset (Fin p)).sum (fun _ : Fin p => (1 : ℂ)) = (p : ℂ) := by
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
  simp

/-- **`Z_P(1) = Σ 1/E_P(α_k) = 0` pour `p ≥ 3`.**

    Preuve : `1/E_P(α_k) = x²·α_k²/p²` ; `Σ α_k² = 0` (cas `e = 2`,
    `p ∤ 2` pour `p ≥ 3`). -/
theorem pandrosion_energy_zeta_at_one_vanishing
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 3 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x) :
    (Finset.univ : Finset (Fin p)).sum
      (fun k : Fin p => 1 / pandrosionEnergy p (cycAnchor α p k)) = 0 := by
  have hp_ge_one : 1 ≤ p := by omega
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hαk_ne : ∀ k : Fin p, cycAnchor α p k ≠ 0 := fun k =>
    cycAnchor_ne_zero_ez x hx p hp_ge_one α hα_pow k
  -- 1/E_P(α_k) = x² · α_k² / p²
  have h_rewrite : ∀ k : Fin p,
      1 / pandrosionEnergy p (cycAnchor α p k) = x^2 * (cycAnchor α p k)^2 / (p : ℂ)^2 := by
    intro k
    unfold pandrosionEnergy pandrosionPDeriv
    have hαk : cycAnchor α p k ≠ 0 := hαk_ne k
    -- α_k^(p-1) = (1/x) / α_k
    have h_pow_pminus1 : (cycAnchor α p k)^(p - 1) = (1 / x) / cycAnchor α p k := by
      have h_succ : (cycAnchor α p k)^(p - 1 + 1)
                  = (cycAnchor α p k)^(p - 1) * cycAnchor α p k := pow_succ _ _
      rw [show p - 1 + 1 = p from Nat.sub_add_cancel hp_ge_one] at h_succ
      have h_pow_eq : (cycAnchor α p k)^p = 1 / x :=
        (cycAnchor_pow α hp_ge_one k).trans hα_pow
      rw [h_pow_eq] at h_succ
      rw [h_succ]
      field_simp
    rw [h_pow_pminus1]
    field_simp
    ring
  rw [Finset.sum_congr rfl (fun k _ => h_rewrite k)]
  rw [show (fun k : Fin p => x^2 * (cycAnchor α p k)^2 / (p : ℂ)^2)
       = (fun k : Fin p => (x^2 / (p : ℂ)^2) * (cycAnchor α p k)^2) from by
    funext k; ring]
  rw [← Finset.mul_sum]
  -- Σ α_k^2 = 0 (e = 2, p ∤ 2 for p ≥ 3)
  have h_not_dvd_2 : ¬ p ∣ 2 := by
    intro h
    have : p ≤ 2 := Nat.le_of_dvd (by norm_num) h
    omega
  rw [sum_cycAnchor_pow_eq_zero_general α p (by omega) 2 (by norm_num) h_not_dvd_2]
  ring

/-- **`Z_P(-1) = Σ E_P(α_k) = 0` pour `p ≥ 3`.**

    Preuve : `E_P(α_k) = p² · α_k^{2(p-1)}` ; `2(p-1) mod p = p - 2 ≠ 0`
    pour `p ≥ 3`, donc `Σ α_k^{2(p-1)} = 0`. -/
theorem pandrosion_energy_zeta_at_neg_one_vanishing
    (α : ℂ) (p : ℕ) (hp : 3 ≤ p) :
    (Finset.univ : Finset (Fin p)).sum
      (fun k : Fin p => pandrosionEnergy p (cycAnchor α p k)) = 0 := by
  -- E_P(α_k) = (p · α_k^(p-1))^2 = p² · α_k^(2(p-1))
  have h_rewrite : ∀ k : Fin p,
      pandrosionEnergy p (cycAnchor α p k) = (p : ℂ)^2 * (cycAnchor α p k)^(2 * (p - 1)) := by
    intro k
    unfold pandrosionEnergy pandrosionPDeriv
    rw [mul_pow, ← pow_mul,
        show (p - 1) * 2 = 2 * (p - 1) from Nat.mul_comm _ _]
  rw [Finset.sum_congr rfl (fun k _ => h_rewrite k)]
  rw [← Finset.mul_sum]
  -- Σ α_k^(2(p-1)) = 0 since p ∤ 2(p-1) for p ≥ 3
  have h_not_dvd : ¬ p ∣ 2 * (p - 1) := by
    intro h
    -- Lift to ℤ to handle subtraction cleanly.
    have hp_ge_one : 1 ≤ p := by omega
    have h_int : (p : ℤ) ∣ (2 * ((p : ℤ) - 1)) := by
      have h_cast : (2 * ((p : ℤ) - 1)) = ((2 * (p - 1) : ℕ) : ℤ) := by
        rw [Nat.cast_mul, Nat.cast_sub hp_ge_one]
        push_cast; ring
      rw [h_cast]
      exact_mod_cast h
    have h_2p : (p : ℤ) ∣ 2 * (p : ℤ) := dvd_mul_left _ _
    have h_diff : (p : ℤ) ∣ (2 * (p : ℤ) - 2 * ((p : ℤ) - 1)) := dvd_sub h_2p h_int
    have h_simp : 2 * (p : ℤ) - 2 * ((p : ℤ) - 1) = 2 := by ring
    rw [h_simp] at h_diff
    have h_p_le : (p : ℤ) ≤ 2 := Int.le_of_dvd (by norm_num) h_diff
    have _hp_ge_int : (3 : ℤ) ≤ (p : ℤ) := by exact_mod_cast hp
    omega
  have h_pos : 1 ≤ 2 * (p - 1) := by omega
  rw [sum_cycAnchor_pow_eq_zero_general α p (by omega) (2 * (p - 1)) h_pos h_not_dvd]
  ring

end EnergyZetaSpecialValues

/-! ============================================================
  §111.4  ★★★★★★★★★★★ Capstone : Energy Zeta Master
============================================================ -/

section EnergyZetaMaster

/-- **★★★★★★★★★★★ PANDROSION ENERGY ZETA MASTER (PDF §4).**

    Pour `p ≥ 3`, **trois valeurs spéciales** de `Z_P(s)` :

    **(A) Trivial value** :  `Z_P(0) = p` (degree).

    **(B) Vanishing à s = 1** : `Z_P(1) = Σ 1/E_P(α_k) = 0`
                                 ("capacitance" nulle).

    **(C) Vanishing à s = -1** : `Z_P(-1) = Σ E_P(α_k) = 0`
                                  (total self-energy nulle).

    **Interprétation physique (PDF §4 Table)** : pour `p ≥ 3`, les
    bassins cyclotomiques de Pandrosion forment un **système
    électrostatique parfaitement équilibré** — la self-energy
    totale (Z_P(-1)) ET la capacitance (Z_P(1)) s'annulent
    simultanément par cancellation cyclotomique. Seule la
    cardinalité (degree, Z_P(0)) persiste.

    Pour `p = 2`, ces vanishings échouent (cas dégénéré). -/
theorem pandrosion_energy_zeta_master
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 3 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x) :
    -- (A) Z_P(0) = p
    (Finset.univ : Finset (Fin p)).sum (fun _ : Fin p => (1 : ℂ)) = (p : ℂ) ∧
    -- (B) Z_P(1) = Σ 1/E_P(α_k) = 0
    (Finset.univ : Finset (Fin p)).sum
      (fun k : Fin p => 1 / pandrosionEnergy p (cycAnchor α p k)) = 0 ∧
    -- (C) Z_P(-1) = Σ E_P(α_k) = 0
    (Finset.univ : Finset (Fin p)).sum
      (fun k : Fin p => pandrosionEnergy p (cycAnchor α p k)) = 0 := by
  refine ⟨?_, ?_, ?_⟩
  · exact pandrosion_energy_zeta_at_zero p
  · exact pandrosion_energy_zeta_at_one_vanishing x hx p hp α hα_pow
  · exact pandrosion_energy_zeta_at_neg_one_vanishing α p hp

end EnergyZetaMaster

end Pandrosion
