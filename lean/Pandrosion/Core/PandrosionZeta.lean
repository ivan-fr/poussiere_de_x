/-
  Universitas Pandrosion — §106. **★★★★★★★★★★★ PANDROSION ZETA
  FUNCTION + MULTI-START SPECTRAL BALANCE THEOREM.**

  Combinaison du corpus multi-start Pandrosion avec la **Pandrosion
  Zeta Function** introduite dans l'article "The Pandrosion Zeta
  Function and the Discriminant as Spectral Determinant" (avril 2026).

  La zeta Pandrosion est :

    `ζ_P(s) := Σ_{k=1}^d P'(α_k)^{-s}`

  Pour Pandrosion : `P(z) = z^p − 1/x`, `P'(α_k) = p/(x·α_k)`, donc
  `1/P'(α_k) = x·α_k/p` — le **poids spectral** `ν_k`.

  Résultats principaux :
    • `ζ_P(1) = Σ_k ν_k = 0` pour `p ≥ 2` (vanishing identity).
    • Combiné avec §95-§105 : **balance spectrale multi-start**.

  Contents.

    §106.1  `pandrosionSpectralWeight` — `ν_k := x · α_k / p`.
    §106.2  `sum_cycAnchor_eq_zero` — Σ γ_k = 0 pour p ≥ 2.
    §106.3  `pandrosion_zeta_vanishing_identity` — ζ_P(1) = 0.
    §106.4  ★★★ `pandrosion_multi_start_spectral_balance` — théorème
            principal à 5 conjuncts.
-/

import Pandrosion.Core.FatouContinuity

namespace Pandrosion

open Complex Finset MeasureTheory

/-! ============================================================
  §106.1  Poids spectral Pandrosion
============================================================ -/

section SpectralWeight

/-- **Poids spectral Pandrosion** au `k`-ième cyclotomic anchor :
    `ν_k := x · α_k / p = 1/P'(α_k)` où `P(z) = z^p − 1/x`. -/
noncomputable def pandrosionSpectralWeight
    (x : ℂ) (p : ℕ) (α : ℂ) (k : Fin p) : ℂ :=
  x * cycAnchor α p k / (p : ℂ)

end SpectralWeight

/-! ============================================================
  §106.2  Somme des anchors cyclotomiques est zéro pour `p ≥ 2`
============================================================ -/

section SumCycAnchor

/-- **Identité clé : `ω^p = 1` avec `ω := exp(2πi/p)`.** -/
theorem zeta_omega_pow_p_eq_one (p : ℕ) (hp : 1 ≤ p) :
    (Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ)))^p = 1 := by
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  rw [← Complex.exp_nat_mul]
  rw [show ((p : ℕ) : ℂ) * (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
       = 2 * (Real.pi : ℂ) * Complex.I from by field_simp; ring]
  exact Complex.exp_two_pi_mul_I

/-- **Identité clé : `ω ≠ 1` pour `p ≥ 2`** (via `IsPrimitiveRoot`). -/
theorem zeta_omega_ne_one (p : ℕ) (hp : 2 ≤ p) :
    Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ)) ≠ 1 := by
  have h_primitive : IsPrimitiveRoot
      (Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))) p :=
    Complex.isPrimitiveRoot_exp p (by omega)
  intro h_eq
  have h_one_pow : (Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ)))^1 = 1 := by
    rw [pow_one, h_eq]
  have h_dvd : p ∣ 1 := (h_primitive.pow_eq_one_iff_dvd 1).mp h_one_pow
  have h_le : p ≤ 1 := Nat.le_of_dvd (by omega) h_dvd
  omega

/-- **`cycAnchor α p k = α · ω^k` où `ω := exp(2πi/p)`.** -/
theorem cycAnchor_eq_alpha_omega_pow (α : ℂ) (p : ℕ) (k : Fin p) :
    cycAnchor α p k
      = α * (Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ)))^(k : ℕ) := by
  unfold cycAnchor
  congr 1
  rw [← Complex.exp_nat_mul]
  congr 1
  push_cast
  ring

/-- **★ Somme des anchors cyclotomiques est zéro pour `p ≥ 2`.** -/
theorem sum_cycAnchor_eq_zero
    (α : ℂ) (p : ℕ) (hp : 2 ≤ p) :
    (Finset.univ : Finset (Fin p)).sum (fun k => cycAnchor α p k) = 0 := by
  set ω := Complex.exp (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
  have hω_ne_one : ω ≠ 1 := zeta_omega_ne_one p hp
  have hω_p_eq_one : ω^p = 1 := zeta_omega_pow_p_eq_one p (by omega)
  have h_rewrite : ∀ k : Fin p, cycAnchor α p k = α * ω^(k : ℕ) := by
    intro k
    exact cycAnchor_eq_alpha_omega_pow α p k
  rw [Finset.sum_congr rfl (fun k _ => h_rewrite k)]
  rw [← Finset.mul_sum]
  suffices h_sum : (Finset.univ : Finset (Fin p)).sum (fun k => ω^(k : ℕ)) = 0 by
    rw [h_sum]; ring
  rw [show (Finset.univ : Finset (Fin p)).sum (fun k => ω^(k : ℕ))
       = (Finset.range p).sum (fun k => ω^k) from by
    rw [← Fin.sum_univ_eq_sum_range]]
  rw [geom_sum_eq hω_ne_one]
  rw [hω_p_eq_one]
  ring

end SumCycAnchor

/-! ============================================================
  §106.3  ζ_P(1) = 0 — vanishing identity
============================================================ -/

section ZetaVanishing

/-- **★★★ Pandrosion zeta vanishing identity : `ζ_P(1) = 0` pour
    `p ≥ 2`.**

    La somme des poids spectraux `ν_k := x · α_k / p` sur toutes les
    racines cyclotomiques est zéro. Equivalent du Théorème 2.1
    "Pandrosion vanishing" de l'article zeta, spécialisé à
    `P(z) = z^p − 1/x`. -/
theorem pandrosion_zeta_vanishing_identity
    (x : ℂ) (p : ℕ) (hp : 2 ≤ p) (α : ℂ) :
    (Finset.univ : Finset (Fin p)).sum
      (fun k => pandrosionSpectralWeight x p α k) = 0 := by
  unfold pandrosionSpectralWeight
  rw [show (fun k : Fin p => x * cycAnchor α p k / (p : ℂ))
       = (fun k => (x / (p : ℂ)) * cycAnchor α p k) from by
    funext k; ring]
  rw [← Finset.mul_sum]
  rw [sum_cycAnchor_eq_zero α p hp]
  ring

end ZetaVanishing

/-! ============================================================
  §106.4  ★★★★★★★★★★★ Multi-Start Spectral Balance Theorem
============================================================ -/

section MultiStartSpectralBalance

/-- **★★★★★★★★★★★ PANDROSION MULTI-START SPECTRAL BALANCE.**

  Théorème combinant le corpus multi-start Pandrosion (§95-§105)
  avec la zeta Pandrosion (PDF "The Pandrosion Zeta Function"), via
  les racines cyclotomiques communes `γ_k = cycAnchor α p k`.

  **Cinq propriétés simultanées** pour tout `(x, p, α)` avec `p ≥ 2`,
  les hypothèses standard Pandrosion :

  **(A) Vanishing identity.** Somme des poids spectraux sur les `p`
  anchors est zéro : `Σ_k ν_k = 0`.

  **(B) Multi-start image.** Pour tout `z ∈ ℂ`, la cible du multi-
  start est un des `p` anchors cyclotomiques — donc a un poids
  spectral bien défini.

  **(C) Unique target a.e.** Pour a.e. `z`, la cible est uniquement
  déterminée (§99).

  **(D) Julia-null.** L'ensemble des `z` où la cible est ambiguë
  a mesure zéro (§103).

  **(E) Continuité a.e.** La fonction cible (version Borel) est
  continue a.e. (§105).

  **Interprétation physique.** Les bassins multi-start forment un
  *système spectral équilibré* : la somme (signée) des poids
  spectraux sur toutes les cibles est zéro. C'est une contrainte
  globale structurelle qui généralise le principe d'équilibre de la
  zeta Pandrosion au cadre algorithmique.

  **Nouveau contenu.** Le corpus Pandrosion (algorithme + dynamique
  + mesure) s'étend maintenant à la théorie spectrale : les cibles
  multi-start sont les "valeurs propres" d'un opérateur Pandrosion
  avec une structure zeta analogue à Minakshisundaram-Pleijel. -/
theorem pandrosion_multi_start_spectral_balance
    (x : ℂ) (p : ℕ) (hp : 2 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0) :
    let hp_ge_one : 1 ≤ p := by omega
    -- (A) Vanishing identity : Σ poids spectraux = 0.
    (Finset.univ : Finset (Fin p)).sum
      (fun k => pandrosionSpectralWeight x p α k) = 0 ∧
    -- (B) Multi-start image dans les cyclotomic anchors.
    (∀ z : ℂ, ∃ k : Fin p,
      niveau5_multi_start_target α p hp_ge_one z = cycAnchor α p k) ∧
    -- (C) Unique target a.e.
    (∀ᵐ z : ℂ ∂volume,
      ∃! s : Fin p, ∀ t : Fin p,
        ‖cycAnchor α p s - z‖ ≤ ‖cycAnchor α p t - z‖) ∧
    -- (D) Julia-null (via §103).
    volume (PandrosionJuliaSet α p) = 0 ∧
    -- (E) Continuité a.e. (via §105).
    (∀ᵐ z : ℂ ∂volume,
      ContinuousAt (niveau5_multi_start_target_borel α p hp_ge_one) z) := by
  have hp_ge_one : 1 ≤ p := by omega
  have hx_ne_zero : x ≠ 0 := by
    intro h_zero
    rw [h_zero] at hα_pow
    have h1 : (1 : ℂ) / (0 : ℂ) = 0 := div_zero 1
    rw [h1] at hα_pow
    have h_alpha_zero : α = 0 :=
      (pow_eq_zero_iff (show p ≠ 0 by omega)).mp hα_pow
    exact hα h_alpha_zero
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · exact pandrosion_zeta_vanishing_identity x p hp α
  · intro z
    exact ⟨nearest_anchor_selector α p hp_ge_one z, rfl⟩
  · exact niveau5_multi_start_target_is_nearest_ae x hx_ne_zero p hp_ge_one α hα hα_pow hSp
  · exact PandrosionJuliaSet_volume_zero α p hα hp_ge_one
  · exact niveau5_multi_start_target_borel_ae_continuous α p hp_ge_one hα

end MultiStartSpectralBalance

end Pandrosion
