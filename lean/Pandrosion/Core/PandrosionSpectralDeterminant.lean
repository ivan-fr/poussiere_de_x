/-
  Universitas Pandrosion — §107. **★★★★★★★★★★★★ PANDROSION
  SPECTRAL DETERMINANT.**

  Pour `P(z) = z^p − 1/x` avec racines cyclotomiques
  `α_k := cycAnchor α p k`, on prouve la **formule produit closed-form** :

    `Π_{k ∈ Fin p} P'(α_k) = (-1)^{p-1} · p^p / x^{p-1}`

  C'est l'analogue de l'égalité `exp(−ζ'_P(0)) = (-1)^{d(d-1)/2} D(P)^2`
  du Théorème 3.1 de l'article "The Pandrosion Zeta Function and the
  Discriminant as Spectral Determinant", spécialisée au polynôme
  `z^p − 1/x`.

  **Vérification Python** : tests p=2..7 avec divers x, ✓.

  Contents.

    §107.1  `fin_sum_coe_complex_two_mul` — `2·∑ k = p(p-1)` en ℂ.
    §107.2  `prod_exp_cycAnchor_eq` — `Π exp(2πi k/p) = (-1)^{p-1}`.
    §107.3  `prod_cycAnchor_eq` — `Π cycAnchor α p k = (-1)^{p-1}·α^p`.
    §107.4  `pandrosionPDeriv` — définition de `P'(z) = p·z^{p-1}`.
    §107.5  ★★★ `pandrosion_spectral_determinant` — théorème principal.
-/

import Pandrosion.Core.PandrosionZeta

namespace Pandrosion

open Complex Finset

/-! ============================================================
  §107.1  Somme Gauss des indices
============================================================ -/

section GaussSum

/-- **`2 · Σ_{k ∈ Fin p} (k : ℂ) = p · (p-1)` en ℂ.** -/
theorem fin_sum_coe_complex_two_mul (p : ℕ) :
    2 * (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => ((k : ℕ) : ℂ))
      = (p : ℂ) * ((p : ℂ) - 1) := by
  rcases Nat.eq_zero_or_pos p with hp_zero | hp_pos
  · subst hp_zero; simp
  · -- Reindex Fin p → range p.
    have h_reindex : (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => ((k : ℕ) : ℂ))
                   = (Finset.range p).sum (fun k : ℕ => (k : ℂ)) := by
      rw [← Fin.sum_univ_eq_sum_range (fun k : ℕ => (k : ℂ))]
    rw [h_reindex]
    -- Σ k en ℕ.
    have h_nat : (Finset.range p).sum (fun k : ℕ => k) * 2 = p * (p - 1) :=
      Finset.sum_range_id_mul_two p
    -- Cast en ℂ.
    have h_sum_cast : ((Finset.range p).sum (fun k : ℕ => (k : ℂ)))
                    = (((Finset.range p).sum (fun k : ℕ => k) : ℕ) : ℂ) := by
      push_cast; rfl
    rw [h_sum_cast]
    have h_complex : 2 * (((Finset.range p).sum (fun k : ℕ => k) : ℕ) : ℂ)
                  = ((p * (p - 1) : ℕ) : ℂ) := by
      have : (((Finset.range p).sum (fun k : ℕ => k) * 2 : ℕ) : ℂ)
           = ((p * (p - 1) : ℕ) : ℂ) := by exact_mod_cast h_nat
      push_cast at this ⊢
      linear_combination this
    rw [h_complex]
    rw [Nat.cast_mul, Nat.cast_sub hp_pos]
    push_cast; ring

end GaussSum

/-! ============================================================
  §107.2  Produit des exponentielles cyclotomiques
============================================================ -/

section ProdExpCyclotomic

/-- **`exp(πi · n) = (-1)^n` pour tout `n : ℕ`.** -/
theorem exp_pi_mul_I_nat (n : ℕ) :
    Complex.exp ((n : ℂ) * ((Real.pi : ℂ) * Complex.I)) = (-1 : ℂ)^n := by
  rw [Complex.exp_nat_mul, Complex.exp_pi_mul_I]

/-- **`Π_{k ∈ Fin p} exp(2πi k/p) = (-1)^{p-1}` pour `p ≥ 1`.** -/
theorem prod_exp_cycAnchor_eq
    (p : ℕ) (hp : 1 ≤ p) :
    (Finset.univ : Finset (Fin p)).prod
      (fun k : Fin p => Complex.exp (2 * (Real.pi : ℂ) * Complex.I * (k : ℂ) / (p : ℂ)))
      = (-1 : ℂ)^(p - 1) := by
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  -- Π exp(z_k) = exp(Σ z_k)
  rw [← Complex.exp_sum]
  -- Σ (2πi k/p) = (p-1) · πi
  have h_sum_eq : (Finset.univ : Finset (Fin p)).sum
        (fun k : Fin p => 2 * (Real.pi : ℂ) * Complex.I * (k : ℂ) / (p : ℂ))
      = ((p - 1 : ℕ) : ℂ) * ((Real.pi : ℂ) * Complex.I) := by
    have h_factor : (Finset.univ : Finset (Fin p)).sum
          (fun k : Fin p => 2 * (Real.pi : ℂ) * Complex.I * (k : ℂ) / (p : ℂ))
        = (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
            * (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => ((k : ℕ) : ℂ)) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro k _
      push_cast
      ring
    rw [h_factor]
    have h_2S := fin_sum_coe_complex_two_mul p
    have h_two_ne : (2 : ℂ) ≠ 0 := two_ne_zero
    -- Multiplier par 2 et utiliser 2·S = p(p-1).
    have h_sub : ((p - 1 : ℕ) : ℂ) = (p : ℂ) - 1 := by
      have h1 : 1 ≤ p := hp
      rw [Nat.cast_sub h1]; push_cast; ring
    rw [h_sub]
    -- Goal: (2πi/p) · S = (p - 1) · πi
    -- 2 · LHS = (2πi/p) · 2S = (2πi/p) · p(p-1) = 2(p-1)πi
    have h_dbl : 2 * ((2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
                  * (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => ((k : ℕ) : ℂ)))
              = 2 * (((p : ℂ) - 1) * ((Real.pi : ℂ) * Complex.I)) := by
      have step1 : 2 * ((2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
                    * (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => ((k : ℕ) : ℂ)))
                = (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
                    * (2 * (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => ((k : ℕ) : ℂ))) := by
        ring
      rw [step1, h_2S]
      field_simp
      ring
    -- Diviser par 2.
    have h_factor_2 :
        (2 * (Real.pi : ℂ) * Complex.I / (p : ℂ))
            * (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => ((k : ℕ) : ℂ))
          = ((p : ℂ) - 1) * ((Real.pi : ℂ) * Complex.I) := by
      have := h_dbl
      have h := mul_left_cancel₀ h_two_ne this
      exact h
    exact h_factor_2
  rw [h_sum_eq]
  exact exp_pi_mul_I_nat (p - 1)

end ProdExpCyclotomic

/-! ============================================================
  §107.3  Produit des cyclotomic anchors
============================================================ -/

section ProdCycAnchor

/-- **★ Produit des anchors cyclotomiques.** -/
theorem prod_cycAnchor_eq
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) :
    (Finset.univ : Finset (Fin p)).prod (fun k => cycAnchor α p k)
      = (-1 : ℂ)^(p - 1) * α^p := by
  unfold cycAnchor
  rw [Finset.prod_mul_distrib]
  rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  rw [prod_exp_cycAnchor_eq p hp]
  ring

/-- **Corollaire : `Π cycAnchor = (-1)^{p-1}/x` sous `α^p = 1/x`.** -/
theorem prod_cycAnchor_eq_of_pow
    (x : ℂ) (p : ℕ) (hp : 1 ≤ p) (α : ℂ) (hα_pow : α ^ p = 1 / x) :
    (Finset.univ : Finset (Fin p)).prod (fun k => cycAnchor α p k)
      = (-1 : ℂ)^(p - 1) / x := by
  rw [prod_cycAnchor_eq α p hp, hα_pow]
  ring

end ProdCycAnchor

/-! ============================================================
  §107.4  Définition `P'`
============================================================ -/

section PandrosionDeriv

/-- **Dérivée `P'(z) = p·z^{p-1}` pour `P(z) = z^p − 1/x`.** -/
noncomputable def pandrosionPDeriv (p : ℕ) (z : ℂ) : ℂ :=
  (p : ℂ) * z ^ (p - 1)

end PandrosionDeriv

/-! ============================================================
  §107.5  ★★★★★★★★★★★★ Pandrosion Spectral Determinant
============================================================ -/

section SpectralDeterminant

/-- **`(-1)^(n*n) = (-1)^n` pour tout `n : ℕ`** (parité préservée). -/
theorem neg_one_pow_mul_self (n : ℕ) :
    (-1 : ℂ)^(n * n) = (-1 : ℂ)^n := by
  rcases Nat.even_or_odd n with h | h
  · have h_sq : Even (n * n) := h.mul_right n
    rw [Even.neg_pow h, Even.neg_pow h_sq, one_pow, one_pow]
  · have h_sq : Odd (n * n) := h.mul h
    rw [Odd.neg_pow h, Odd.neg_pow h_sq, one_pow, one_pow]

/-- **★★★★★★★★★★★★ PANDROSION SPECTRAL DETERMINANT.**

    Pour `P(z) = z^p − 1/x` avec racines cyclotomiques
    `α_k := cycAnchor α p k`, on a :

    `Π_{k ∈ Fin p} P'(α_k) = (-1)^{p-1} · p^p / x^{p-1}`.

    **Interprétation (PDF Théorème 3.1)** :

      `exp(−ζ'_P(0)) = Π_k P'(α_k) = (-1)^{d(d-1)/2} · D(P)^2`

    où `D(P)` est le discriminant. Ceci donne au discriminant de
    `z^p − 1/x` une interprétation **spectrale** via la zeta
    Pandrosion. Combinée avec §106 (vanishing identity `ζ_P(1) = 0`),
    on obtient l'image analytique complète du root spectrum
    Pandrosion.

    **Chaîne de preuve** :
      1. `Π_k P'(α_k) = Π_k (p · α_k^{p-1}) = p^p · (Π_k α_k)^{p-1}`.
      2. `Π_k α_k = (-1)^{p-1} / x` (§107.3).
      3. `((-1)^{p-1}/x)^{p-1} = (-1)^{(p-1)·(p-1)} / x^{p-1}
                              = (-1)^{p-1} / x^{p-1}` (parité).
      4. Combiner : `p^p · (-1)^{p-1} / x^{p-1}`. -/
theorem pandrosion_spectral_determinant
    (x : ℂ) (_hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x) :
    (Finset.univ : Finset (Fin p)).prod
      (fun k => pandrosionPDeriv p (cycAnchor α p k))
      = (-1 : ℂ)^(p - 1) * (p : ℂ)^p / x^(p - 1) := by
  unfold pandrosionPDeriv
  -- Π (p · α_k^(p-1)) = p^p · Π α_k^(p-1) = p^p · (Π α_k)^(p-1)
  rw [Finset.prod_mul_distrib]
  rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  rw [Finset.prod_pow]
  rw [prod_cycAnchor_eq_of_pow x p hp α hα_pow]
  -- p^p · ((-1)^(p-1)/x)^(p-1) = p^p · (-1)^((p-1)*(p-1)) / x^(p-1)
  --                            = p^p · (-1)^(p-1) / x^(p-1)
  rw [div_pow, ← pow_mul, neg_one_pow_mul_self]
  ring

end SpectralDeterminant

end Pandrosion
