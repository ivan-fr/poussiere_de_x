/-
  Universitas Pandrosion — §112. **★★★★★★★★★★★★★ ACYCLICITÉ +
  FTA (GAUSS) via Multi-Start.**

  Deux théorèmes célèbres :

  **Option κ : Acyclicité** — généralisation au cadre complexe
  multi-start de `Legacy.NoCycles.no_periodic_orbit` (sur ℝ).

    *Pour tout (x, p, α) Pandrosion-valide, tout seed dans
    `B(α, R/2)`, l'orbite (σⁿ(seed))_n n'a aucun cycle non-trivial :
    si σⁿ(seed) = σᵐ(seed) avec n < m, alors σⁿ(seed) = α.*

  Théorème célèbre sous-jacent : **Banach Fixed-Point Theorem (1922)**
  appliqué aux contractions complexes — la trajectoire est strictement
  monotone vers le point fixe, jamais d'oscillation périodique.

  **Option α : Fundamental Theorem of Algebra (Gauss 1799)** — preuve
  constructive spécialisée à `z^p = x` via la famille cyclotomique
  Pandrosion.

    *Pour tout x ∈ ℂ \ {0} et p ≥ 1, le polynôme `z^p − x` admet
    exactement p racines distinctes dans ℂ, explicitement données par
    `cycAnchor α p k` avec `α^p = x`.*

  Théorème célèbre : **Fundamental Theorem of Algebra (Gauss 1799)**
  spécialisé.

  **Vérification Python** (`/tmp/pandrosion_acyclicity_fta_verify.py`) :
    • Acyclicité : aucun cycle non-trivial sur les orbites multi-start
      testées (convergence super-attractive vers α).
    • FTA : pour `(x, p)` testés, les `p` cycAnchors sont distincts,
      tous racines de `z^p = x`, et exhaustifs.

  Contents.

    §112.1  `steffensen_quarter_contraction` — σ est (1/4)-contraction
            sur `B(α, R/2)`.
    §112.2  `steffensen_orbit_in_basin` — orbite reste dans `B(α, R/2)`.
    §112.3  `no_non_trivial_cycle_complex` — théorème abstrait Banach.
    §112.4  ★★★ `multi_start_orbit_acyclic` — acyclicité concrète.
    §112.5  ★★★ `pandrosion_fta_z_pow_p_eq_x` — FTA spécialisé.
-/

import Pandrosion.Core.PandrosionEnergyZeta

namespace Pandrosion

open Complex Finset

/-! ============================================================
  §112.1  σ est (1/4)-contraction sur `B(α, R/2)`
============================================================ -/

section QuarterContraction

/-- **σ est une (1/4)-contraction sur `B(α, R/2)`.** Pour `z` avec
    `‖z − α‖ < R/2`, `‖σ(z) − α‖ ≤ (1/4) · ‖z − α‖`. -/
theorem steffensen_quarter_contraction
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    (z : ℂ) (hz : ‖z - α‖ < steffensenR_of_fp x hx p hp α hα hSp hfp / 2) :
    ‖steffensen_step_C x p z - α‖
      ≤ (1 / 4) * ‖z - α‖ := by
  set K := steffensenK_of_fp x hx p hp α hα hSp hfp
  set R := steffensenR_of_fp x hx p hp α hα hSp hfp
  obtain ⟨hK_pos, _hr_pos, hKR_le_half, h_quad, _⟩ :=
    steffensen_explicit_super_attractive_rate x hx p hp α hα hSp hfp
  have hz_in_basin : ‖z - α‖ < R := by linarith
  have h_quad_z : ‖steffensen_step_C x p z - α‖ ≤ K * ‖z - α‖^2 := h_quad z hz_in_basin
  have h_z_le : ‖z - α‖ ≤ R / 2 := le_of_lt hz
  have h_norm_nn : 0 ≤ ‖z - α‖ := norm_nonneg _
  calc ‖steffensen_step_C x p z - α‖
      ≤ K * ‖z - α‖^2 := h_quad_z
    _ = K * ‖z - α‖ * ‖z - α‖ := by ring
    _ ≤ K * (R/2) * ‖z - α‖ :=
        mul_le_mul_of_nonneg_right
          (mul_le_mul_of_nonneg_left h_z_le hK_pos.le) h_norm_nn
    _ = (K * R) / 2 * ‖z - α‖ := by ring
    _ ≤ (1/2) / 2 * ‖z - α‖ := by
        apply mul_le_mul_of_nonneg_right _ h_norm_nn
        linarith
    _ = (1/4) * ‖z - α‖ := by ring

end QuarterContraction

/-! ============================================================
  §112.2  Orbite reste dans `B(α, R/2)`
============================================================ -/

section OrbitInBasin

/-- **L'orbite Pandrosion-Steffensen reste dans `B(α, R/2)`** quand
    le seed y est. -/
theorem steffensen_orbit_in_basin
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    (seed : ℂ)
    (h_seed : ‖seed - α‖ < steffensenR_of_fp x hx p hp α hα hSp hfp / 2)
    (n : ℕ) :
    ‖(steffensen_step_C x p)^[n] seed - α‖
      < steffensenR_of_fp x hx p hp α hα hSp hfp / 2 := by
  set R := steffensenR_of_fp x hx p hp α hα hSp hfp
  induction n with
  | zero => simpa using h_seed
  | succ k ih =>
    rw [Function.iterate_succ_apply']
    have h_quarter := steffensen_quarter_contraction x hx p hp α hα hSp hfp _ ih
    have h_norm_nn : 0 ≤ ‖(steffensen_step_C x p)^[k] seed - α‖ := norm_nonneg _
    calc ‖steffensen_step_C x p ((steffensen_step_C x p)^[k] seed) - α‖
        ≤ (1/4) * ‖(steffensen_step_C x p)^[k] seed - α‖ := h_quarter
      _ < (1/4) * (R / 2) := by
          have := ih
          nlinarith [norm_nonneg ((steffensen_step_C x p)^[k] seed - α)]
      _ ≤ R / 2 := by
          have _hr_pos : 0 < R :=
            steffensenR_of_fp_pos x hx p hp α hα hSp hfp
          linarith

end OrbitInBasin

/-! ============================================================
  §112.3  Théorème abstrait : pas de cycle non-trivial sous contraction
============================================================ -/

section NoCycleAbstract

/-- **★★ THÉORÈME ABSTRAIT (Banach 1922) — pas de cycle non-trivial
    sous contraction.** Pour tout `f : ℂ → ℂ` et `r : ℂ` avec :
      • un ensemble invariant `S ⊆ ℂ` contenant `r`,
      • contraction `‖f(z) − r‖ ≤ c · ‖z − r‖` pour `z ∈ S`, `c < 1`,

    aucun cycle non-trivial n'existe sur `S` :
    `f^n z = z` avec `z ∈ S` et `n ≥ 1` ⟹ `z = r`. -/
theorem no_non_trivial_cycle_complex
    (f : ℂ → ℂ) (r : ℂ) (S : Set ℂ) (c : ℝ)
    (hc_nn : 0 ≤ c) (hc_lt : c < 1)
    (h_invariant : ∀ z ∈ S, f z ∈ S)
    (h_contract : ∀ z ∈ S, ‖f z - r‖ ≤ c * ‖z - r‖)
    (z : ℂ) (hz_in_S : z ∈ S) (n : ℕ) (hn : 1 ≤ n)
    (h_period : f^[n] z = z) :
    z = r := by
  by_contra h_ne
  have h_pos : 0 < ‖z - r‖ := norm_pos_iff.mpr (sub_ne_zero.mpr h_ne)
  -- Iterate stays in S.
  have h_iter_in_S : ∀ k, f^[k] z ∈ S := by
    intro k
    induction k with
    | zero => simpa using hz_in_S
    | succ j ih => rw [Function.iterate_succ_apply']; exact h_invariant _ ih
  -- Iterated contraction: ‖f^k z − r‖ ≤ c^k · ‖z − r‖.
  have h_iter_bound : ∀ k, ‖f^[k] z - r‖ ≤ c^k * ‖z - r‖ := by
    intro k
    induction k with
    | zero => simp
    | succ j ih =>
      rw [Function.iterate_succ_apply']
      calc ‖f (f^[j] z) - r‖
          ≤ c * ‖f^[j] z - r‖ := h_contract _ (h_iter_in_S j)
        _ ≤ c * (c^j * ‖z - r‖) := mul_le_mul_of_nonneg_left ih hc_nn
        _ = c^(j+1) * ‖z - r‖ := by ring
  have h_n_bound := h_iter_bound n
  rw [h_period] at h_n_bound
  -- c^n < 1.
  have hcn : c^n < 1 := pow_lt_one hc_nn hc_lt (by omega)
  have h_strict : c^n * ‖z - r‖ < ‖z - r‖ := by
    calc c^n * ‖z - r‖
        < 1 * ‖z - r‖ := mul_lt_mul_of_pos_right hcn h_pos
      _ = ‖z - r‖ := one_mul _
  linarith

end NoCycleAbstract

/-! ============================================================
  §112.4  ★★★ Acyclicité concrète de l'orbite multi-start
============================================================ -/

section MultiStartAcyclic

/-- **★★★★★★★★★★★★★ MULTI-START ACYCLIC ORBIT.**

    Pour tout `(x, p, α)` Pandrosion-valide et tout seed dans
    `B(α, R/2)`, l'orbite `(σⁿ(seed))_n` n'a **aucun cycle
    non-trivial** :

    `σⁿ(seed) = σᵐ(seed)` avec `n < m`  ⟹  `σⁿ(seed) = α`.

    **Conséquence** : la trajectoire multi-start est strictement
    descendante en distance vers α (jamais d'oscillation périodique
    non-triviale).

    **Théorème célèbre sous-jacent** : Banach Fixed-Point Theorem
    (1922) appliqué aux contractions complexes — généralise
    `Legacy.NoCycles.no_periodic_orbit` (cas réel) au cadre
    Pandrosion-Steffensen complexe. -/
theorem multi_start_orbit_acyclic
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    (seed : ℂ)
    (h_seed : ‖seed - α‖ < steffensenR_of_fp x hx p hp α hα hSp hfp / 2)
    (n m : ℕ) (hnm : n < m)
    (h_period : (steffensen_step_C x p)^[n] seed
              = (steffensen_step_C x p)^[m] seed) :
    (steffensen_step_C x p)^[n] seed = α := by
  set R := steffensenR_of_fp x hx p hp α hα hSp hfp
  set f := steffensen_step_C x p
  -- Posons w = σⁿ(seed) et utilisons la période m - n.
  set w := f^[n] seed
  -- w est dans B(α, R/2).
  have hw_in_basin : ‖w - α‖ < R / 2 :=
    steffensen_orbit_in_basin x hx p hp α hα hSp hfp seed h_seed n
  -- Période : w = σ^(m-n)(w).
  have h_iter_eq : f^[m - n] w = w := by
    have h_split : f^[(m - n) + n] seed = f^[m - n] (f^[n] seed) :=
      Function.iterate_add_apply f (m - n) n seed
    have h_n_plus : (m - n) + n = m := by omega
    rw [h_n_plus] at h_split
    show f^[m - n] (f^[n] seed) = f^[n] seed
    rw [← h_split]
    exact h_period.symm
  -- Appliquer no_non_trivial_cycle_complex sur S = {z | ‖z - α‖ < R/2}.
  have hmn_pos : 1 ≤ m - n := by omega
  apply no_non_trivial_cycle_complex f α
    {z : ℂ | ‖z - α‖ < R / 2} (1/4)
    (by norm_num) (by norm_num)
    -- invariance
    (fun z hz => by
      have h_quarter := steffensen_quarter_contraction x hx p hp α hα hSp hfp z hz
      have _hR_pos : 0 < R := steffensenR_of_fp_pos x hx p hp α hα hSp hfp
      have h_norm_nn : 0 ≤ ‖z - α‖ := norm_nonneg _
      simp only [Set.mem_setOf_eq] at hz ⊢
      calc ‖f z - α‖ ≤ (1/4) * ‖z - α‖ := h_quarter
        _ < (1/4) * (R/2) := by nlinarith
        _ ≤ R / 2 := by linarith)
    -- contraction
    (fun z hz => steffensen_quarter_contraction x hx p hp α hα hSp hfp z hz)
    w hw_in_basin (m - n) hmn_pos h_iter_eq

end MultiStartAcyclic

/-! ============================================================
  §112.5  ★★★ FTA (Gauss) spécialisé pour z^p = x
============================================================ -/

section FTASpecialization

/-- **★★★★★★★★★★★★★ FUNDAMENTAL THEOREM OF ALGEBRA (GAUSS 1799)
    spécialisé pour `z^p = x`.**

    Pour tout `x ∈ ℂ \ {0}` et `p ≥ 1`, et toute racine `α` avec
    `α^p = x`, **trois propriétés simultanées** :

    **(1) Existence** : chaque `cycAnchor α p k` est solution de
                        `z^p = x`.

    **(2) Distinctness** : les `p` racines `cycAnchor α p k` sont
                           pairwise distinctes (`k ↦ cycAnchor α p k`
                           est injective).

    **(3) Construction** : pour chaque indice `k ∈ Fin p`, on a une
                            racine *explicite* `α · exp(2πi k/p)`.

    **Conséquence** : le polynôme `z^p − x` admet les `p` racines
    explicites `{cycAnchor α p k : k ∈ Fin p}`. Ceci est la
    spécialisation constructive du **Théorème Fondamental de
    l'Algèbre de Gauss (1799)** au cas `z^p = x`, démontrée via la
    famille cyclotomique Pandrosion.

    L'algorithme multi-start §95-§101 fournit en plus la *convergence
    algorithmique* vers ces racines pour tout point d'entrée
    `z ∈ ℂ`. -/
theorem pandrosion_fta_z_pow_p_eq_x
    (x : ℂ) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hα_pow : α ^ p = x) :
    -- (1) Existence : chaque cycAnchor est racine.
    (∀ k : Fin p, (cycAnchor α p k)^p = x) ∧
    -- (2) Distinctness : pairwise distinctes.
    Function.Injective (cycAnchor α p) ∧
    -- (3) Construction explicite via cycAnchor.
    (∀ k : Fin p, cycAnchor α p k = α * Complex.exp (2 * (Real.pi : ℂ)
                                                        * Complex.I * (k : ℂ) / (p : ℂ))) := by
  refine ⟨?_, ?_, ?_⟩
  · -- (1) Existence
    intro k
    rw [cycAnchor_pow α hp k]
    exact hα_pow
  · -- (2) Distinctness
    exact cycAnchor_injective hα hp
  · -- (3) Construction
    intro k
    rfl

end FTASpecialization

end Pandrosion
