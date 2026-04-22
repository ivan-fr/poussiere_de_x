/-
  Universitas Pandrosion — §62. **Countability of the iterated
  `h`-preimages of the Sp-zeros at `x = 2, p = 3`.**

  Ce module consolide la contribution mesure-théorique inconditionnelle
  de §61 en passant des zéros directs aux **préimages itérées** de
  l'itération Pandrosion `h`. La structure clé :

    • `h` est une fraction rationnelle de degré ≤ 2. Pour tout `w ∈ ℂ`,
      `h⁻¹({w})` est l'ensemble des racines d'un polynôme de degré 2,
      donc **fini** (cardinalité ≤ 2).
    • Par induction, `h⁻ⁿ(F)` est fini pour tout `F` fini.
    • Une union dénombrable d'ensembles finis est **dénombrable**.
    • Un ensemble dénombrable dans ℂ est de mesure de Lebesgue nulle.

  Ce résultat **réduit** la conjecture résiduelle
  `ComplexJuliaP3MeasureZero 2 α` à :

      *"Tout point hors de l'orbite singulière h-itérée des Sp-zéros
      converge vers l'une des trois racines cyclotomiques."*

  C'est une **conjecture purement dynamique** (Fatou-exhaustion),
  sans composante mesure-théorique résiduelle.

  Contents.

    §62.1  `h2_singleton_preimage_finite` — `h⁻¹({w}).Finite`.

    §62.2  `h2_preimage_finite_of_finite` — extension aux F finis.

    §62.3  `h2_iterate_preimage_finite` — finitude en n itérations.

    §62.4  `h2_singular_orbit_countable` — dénombrabilité totale.

    §62.5  `h2_singular_orbit_measure_zero` — mesure nulle.
-/

import Pandrosion.Core.JuliaNullX2P3
import Mathlib.Data.Polynomial.RingDivision

namespace Pandrosion

open Complex MeasureTheory

/-- Shorthand for the Pandrosion `h` map at `x = 2, p = 3`. -/
noncomputable def h2 : ℂ → ℂ := pandrosion_h_C (2 : ℂ) 3

/-! ============================================================
  §62.1  Single-preimage is finite
============================================================ -/

/-- **★ Preimage of a single point is finite (cardinalité ≤ 2).**

    Pour tout `w ∈ ℂ`, `{z : h₂(z) = w}` est fini.

    *Idée.* `h(z) = w` équivaut (algébriquement) à
    `2(1-w)·(1 + z + z²) = 1` — polynôme de degré 2 (si `w ≠ 1`),
    ou à `Sp(z) = 0` (si `w = 1`, grâce à la convention Lean `1/0 = 0`).
    Dans les deux cas on a une **équation polynomiale de degré ≤ 2**,
    donc ≤ 2 solutions. -/
theorem h2_singleton_preimage_finite (w : ℂ) :
    ({z : ℂ | h2 z = w}).Finite := by
  by_cases hw : w = 1
  · -- Case w = 1: preimage = Sp-zeros ∪ (points with Sp ≠ 0 and h(z) = 1).
    -- Second set: (x-1)/(x·Sp) = 0 with Sp ≠ 0 ⟹ x-1 = 0, but x = 2 so false.
    -- So preimage = {z : Sp_C 3 z = 0}, finite by §61.
    subst hw
    have h_set_eq :
        {z : ℂ | h2 z = 1} = {z : ℂ | Sp_C 3 z = 0} := by
      ext z
      simp only [Set.mem_setOf_eq, h2, pandrosion_h_C]
      constructor
      · intro hz
        -- 1 - (2-1)/(2·Sp) = 1 ⟹ 1/(2·Sp) = 0
        have h_div_zero : ((2 : ℂ) - 1) / ((2 : ℂ) * Sp_C 3 z) = 0 := by
          linear_combination -hz
        -- `a/b = 0` ⟹ `a = 0 ∨ b = 0`. Here `2 - 1 = 1 ≠ 0`, so `2·Sp = 0`.
        rw [div_eq_zero_iff] at h_div_zero
        rcases h_div_zero with h_num | h_den
        · exfalso; norm_num at h_num
        · have h_two_ne : (2 : ℂ) ≠ 0 := by norm_num
          rcases mul_eq_zero.mp h_den with h1 | h2
          · exact absurd h1 h_two_ne
          · exact h2
      · intro hSp
        rw [hSp]
        simp
    rw [h_set_eq]
    exact Sp_C_p3_zeros_finite
  · -- Case w ≠ 1: preimage ⊆ roots of the quadratic
    --   Q(z) := 2(1-w)·(1 + z + z²) - 1 = 0.
    have h1w_ne : (1 - w) ≠ 0 := sub_ne_zero.mpr (Ne.symm hw)
    have h_coeff_ne : (2 : ℂ) * (1 - w) ≠ 0 := mul_ne_zero (by norm_num) h1w_ne
    -- Build polynomial Q = 2(1-w)·X² + 2(1-w)·X + (2(1-w) - 1).
    set a : ℂ := (2 : ℂ) * (1 - w)
    have ha_ne : a ≠ 0 := h_coeff_ne
    let Q : Polynomial ℂ :=
      Polynomial.C a * Polynomial.X^2 + Polynomial.C a * Polynomial.X
        + Polynomial.C (a - 1)
    have hQ_ne : Q ≠ 0 := by
      intro h
      have h_deg : Q.natDegree = 2 := by
        show (Polynomial.C a * Polynomial.X^2 + Polynomial.C a * Polynomial.X
              + Polynomial.C (a - 1)).natDegree = 2
        compute_degree!
      rw [h] at h_deg
      simp at h_deg
    have h_subset : {z : ℂ | h2 z = w} ⊆ {z : ℂ | Q.IsRoot z} := by
      intro z hz
      simp only [Set.mem_setOf_eq, h2, pandrosion_h_C] at hz
      -- Case Sp_C 3 z = 0.
      by_cases hSp : Sp_C 3 z = 0
      · -- h(z) with Sp=0 evaluates to 1 (Lean convention 1/0=0).
        exfalso
        apply hw
        rw [← hz, hSp]; simp
      · -- Sp ≠ 0. Clear denominators.
        have h2Sp_ne : (2 : ℂ) * Sp_C 3 z ≠ 0 :=
          mul_ne_zero (by norm_num) hSp
        -- 1 - (x-1)/(x·Sp) = w ⟹ (x-1) = (1-w)·x·Sp ⟹
        -- 2·(1-w)·Sp = 1, i.e., a·Sp = 1.
        have h_aSp : a * Sp_C 3 z = 1 := by
          have h_frac : ((2 : ℂ) - 1) / ((2 : ℂ) * Sp_C 3 z) = 1 - w := by
            linear_combination -hz
          field_simp [h2Sp_ne] at h_frac
          unfold_let a
          linear_combination -h_frac
        -- Sp(z) = 1 + z + z², so a·(1 + z + z²) = 1, i.e. Q(z) = 0.
        rw [Sp_C_p3] at h_aSp
        show Q.IsRoot z
        unfold_let Q
        simp [Polynomial.IsRoot, Polynomial.eval_add, Polynomial.eval_mul,
              Polynomial.eval_X, Polynomial.eval_pow, Polynomial.eval_C]
        linear_combination h_aSp
    exact (Polynomial.finite_setOf_isRoot hQ_ne).subset h_subset

/-! ============================================================
  §62.2  Preimage of a finite set
============================================================ -/

/-- **★ Preimage of a finite set is finite.** -/
theorem h2_preimage_finite_of_finite (F : Set ℂ) (hF : F.Finite) :
    ({z : ℂ | h2 z ∈ F}).Finite := by
  have h_union : {z : ℂ | h2 z ∈ F} = ⋃ w ∈ F, {z : ℂ | h2 z = w} := by
    ext z; simp [Set.mem_iUnion]
  rw [h_union]
  exact hF.biUnion (fun w _ => h2_singleton_preimage_finite w)

/-! ============================================================
  §62.3  n-fold preimage is finite
============================================================ -/

/-- n-fold preimage sequence. -/
noncomputable def h2_iterPreimage (n : ℕ) (F : Set ℂ) : Set ℂ :=
  {z : ℂ | (h2)^[n] z ∈ F}

/-- **★ n-fold preimage of a finite set is finite.** -/
theorem h2_iterPreimage_finite (F : Set ℂ) (hF : F.Finite) (n : ℕ) :
    (h2_iterPreimage n F).Finite := by
  induction n with
  | zero =>
    simp [h2_iterPreimage, Function.iterate_zero, id]
    exact hF
  | succ k ih =>
    have h_rewrite : h2_iterPreimage (k + 1) F = {z : ℂ | h2 z ∈ h2_iterPreimage k F} := by
      ext z
      simp only [h2_iterPreimage, Set.mem_setOf_eq]
      rw [Function.iterate_succ, Function.comp_apply]
    rw [h_rewrite]
    exact h2_preimage_finite_of_finite _ ih

/-! ============================================================
  §62.4  Union over n is countable
============================================================ -/

/-- **★★★ L'orbite h-singulière est dénombrable.**

    `⋃_n h⁻ⁿ({ω, ω²})` est une union dénombrable d'ensembles finis,
    donc dénombrable. -/
theorem h2_singular_orbit_countable :
    (⋃ n : ℕ, h2_iterPreimage n ({omega3, omega3 ^ 2} : Set ℂ)).Countable := by
  apply Set.countable_iUnion
  intro n
  exact (h2_iterPreimage_finite _ (Set.toFinite _) n).countable

/-! ============================================================
  §62.5  Measure zero
============================================================ -/

/-- **★★★★ L'orbite h-singulière est de mesure de Lebesgue nulle.**

    Brique inconditionnelle : la **partie mesure-théorique** de
    `ComplexJuliaP3MeasureZero 2 α` est prouvée. Il ne reste qu'une
    conjecture dynamique (Fatou-exhaustion en dehors de l'orbite). -/
theorem h2_singular_orbit_measure_zero :
    volume (⋃ n : ℕ, h2_iterPreimage n ({omega3, omega3 ^ 2} : Set ℂ)) = 0 :=
  h2_singular_orbit_countable.measure_zero _

end Pandrosion
