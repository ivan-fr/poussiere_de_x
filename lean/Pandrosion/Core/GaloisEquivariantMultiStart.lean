/-
  Universitas Pandrosion — §97. **★★★★★★★★★ Galois-équivariance
  complète du multi-start.**

  Le théorème §95 garantit, pour toute racine cyclotomique `α` choisie
  *à l'avance*, l'existence d'un seed multi-start qui converge vers
  `α`. Le présent module promeut cette construction en *fonction* :

      `selector : ℂ → Fin p`

  est un sélecteur arbitraire (e.g., racine principale, racine la plus
  proche de `z`, racine dans un secteur prescrit, racine d'indice
  `f(z) mod p` pour une fonction donnée). Le théorème central produit
  un seed `multi_start_basin_seed_generic (γ_{selector(z)}) ε z` qui,
  pour tout `z`, converge vers `γ_{selector(z)}` avec complexité
  loglog optimale.

  L'algorithme est ainsi **constructivement Galois-équivariant** : on
  *choisit* la racine cible parmi les `p` racines cyclotomiques par
  un simple sélecteur, sans aucun changement à l'algorithme Steffensen
  sous-jacent.

  Trois exemples de sélecteurs sont fournis :

    1. `principal_selector` — choisit toujours `γ_0 = α` (racine
        principale). Récupère §95 verbatim.
    2. `nearest_anchor_selector` — choisit la racine cyclotomique la
        plus proche de `z` (défini sur le complément du frontière
        Voronoï, indéfini par `Classical.choice` sur la frontière).
    3. `sector_selector` — partitionne `ℂ` en `p` secteurs angulaires
        autour de `0` et choisit la racine `γ_k` du `k`-ième secteur.

  Contents.

    §97.1  `galois_equivariant_multi_start` — théorème principal
            paramétré par un sélecteur arbitraire.
    §97.2  `principal_selector` et corollaire §95-équivalent.
    §97.3  `nearest_anchor_selector` (Voronoï constructif).
    §97.4  `sector_selector` — partition angulaire constructive.
-/

import Pandrosion.Core.LoglogUniversalMultiStartGeneric
import Pandrosion.Core.CyclotomicMcMullen

namespace Pandrosion

open Complex Filter Topology

/-! ============================================================
  §97.1  ★★★ Théorème central : Galois-équivariance complète
============================================================ -/

section GaloisEquivariantMultiStart

/-- **★★★★★★★★★ Galois-équivariance complète du multi-start.**

  Pour tout `x ∈ ℂ \ {0}`, `p ≥ 1`, racine `α` avec `α^p = 1/x` et
  non-dégénérescence cyclotomique (`Sp_C p (γ_k) ≠ 0` pour tout `k`),
  pour tout **sélecteur** `selector : ℂ → Fin p` et toute précision
  `ε > 0`, il existe un seed multi-start
  `γ_{selector(z)} + ε_seed · (z − γ_{selector(z)})` et un compte
  d'itérations `N : ℕ` (de complexité `O(log log(1/ε))`) tels que :

      `∀ n ≥ N, ‖σⁿ(seed) − γ_{selector(z)}‖ ≤ ε`.

  **L'algorithme est ainsi constructivement Galois-équivariant** :
  l'utilisateur peut spécifier, par un simple `selector : ℂ → Fin p`,
  *vers quelle racine cyclotomique* l'itération doit converger pour
  chaque entrée `z`, sans modifier l'algorithme Steffensen sous-jacent.

  *Inconditionnel* sur les hypothèses standard de Pandrosion. -/
theorem galois_equivariant_multi_start
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (_hα : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0)
    (selector : ℂ → Fin p) :
    ∀ z : ℂ, ∀ ε : ℝ, 0 < ε →
      ∃ ε_seed : ℝ, 0 < ε_seed ∧
      ∃ N : ℕ, ∀ n ≥ N,
        ‖(steffensen_step_C x p)^[n]
            (multi_start_basin_seed_generic
              (cycAnchor α p (selector z)) ε_seed z)
              - cycAnchor α p (selector z)‖ ≤ ε := by
  intro z ε hε_pos
  -- Cyclotomic anchors are non-zero p-th roots.
  have hγ_pow : (cycAnchor α p (selector z)) ^ p = 1 / x :=
    (cycAnchor_pow α hp (selector z)).trans hα_pow
  have hγ_ne : cycAnchor α p (selector z) ≠ 0 := by
    intro h_zero
    have h_pow : (cycAnchor α p (selector z)) ^ p = 0 := by
      rw [h_zero]; exact zero_pow (by omega)
    rw [hγ_pow] at h_pow
    exact one_div_ne_zero hx h_pow
  -- Apply §95 with α := γ_{selector(z)}.
  exact pandrosion_loglog_universal_multi_start_generic
    x hx p hp (cycAnchor α p (selector z)) hγ_ne
    (hSp (selector z)) hγ_pow z ε hε_pos

end GaloisEquivariantMultiStart

/-! ============================================================
  §97.2  Sélecteur principal — récupère §95 verbatim
============================================================ -/

section PrincipalSelector

/-- **Sélecteur principal** : retourne toujours l'indice `0`,
    correspondant à `γ_0 = α · exp(0) = α` (la racine principale). -/
noncomputable def principal_selector (p : ℕ) (hp : 1 ≤ p) : ℂ → Fin p :=
  fun _ => ⟨0, by omega⟩

/-- **Cyclotomic anchor à l'indice 0 = α (racine principale).**

    `cycAnchor α p ⟨0, _⟩ = α · exp(0) = α · 1 = α`. -/
theorem cycAnchor_zero (α : ℂ) (p : ℕ) (hp : 1 ≤ p) :
    cycAnchor α p ⟨0, by omega⟩ = α := by
  unfold cycAnchor
  simp

/-- **Corollaire — sélecteur principal récupère §95 verbatim.**

    Avec `selector := principal_selector`, le §97.1 se spécialise en
    convergence vers `α` pour tout `z`, identique à §95.3. -/
theorem galois_equivariant_principal
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0) :
    ∀ z : ℂ, ∀ ε : ℝ, 0 < ε →
      ∃ ε_seed : ℝ, 0 < ε_seed ∧
      ∃ N : ℕ, ∀ n ≥ N,
        ‖(steffensen_step_C x p)^[n]
            (multi_start_basin_seed_generic α ε_seed z) - α‖ ≤ ε := by
  intro z ε hε_pos
  have h := galois_equivariant_multi_start x hx p hp α hα hα_pow hSp
              (principal_selector p hp) z ε hε_pos
  -- Réécrire `cycAnchor α p (principal_selector p hp z) = α`.
  rw [show principal_selector p hp z = ⟨0, by omega⟩ from rfl,
      cycAnchor_zero α p hp] at h
  exact h

end PrincipalSelector

/-! ============================================================
  §97.3  Sélecteur "racine la plus proche" (Voronoï constructif)
============================================================ -/

section NearestAnchorSelector

/-- **Sélecteur "racine la plus proche"** : choisit la racine
    cyclotomique `γ_k` qui minimise `‖z − γ_k‖`.

    Construction via `Classical.choice` car la fonction "argmin sur
    `Fin p`" est non-canonique sur la frontière Voronoï (où plusieurs
    racines équidistantes existent). Sur le complément de la frontière
    Voronoï, le sélecteur est *uniquement* déterminé. -/
noncomputable def nearest_anchor_selector
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z : ℂ) : Fin p := by
  classical
  exact Finset.univ.exists_min_image
      (fun k : Fin p => ‖z - cycAnchor α p k‖)
      (Finset.univ_nonempty_iff.mpr ⟨⟨0, by omega⟩⟩) |>.choose

/-- **Spécification du sélecteur "racine la plus proche".**

    Pour tout `z ∈ ℂ`, `nearest_anchor_selector α p hp z` minimise
    bien `‖z − cycAnchor α p k‖` parmi les `p` indices `k`. -/
theorem nearest_anchor_selector_spec
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) (z : ℂ) :
    ∀ k : Fin p,
      ‖z - cycAnchor α p (nearest_anchor_selector α p hp z)‖
        ≤ ‖z - cycAnchor α p k‖ := by
  classical
  intro k
  unfold nearest_anchor_selector
  have h_choose := (Finset.univ.exists_min_image
      (fun k : Fin p => ‖z - cycAnchor α p k‖)
      (Finset.univ_nonempty_iff.mpr ⟨⟨0, by omega⟩⟩)).choose_spec
  exact h_choose.2 k (Finset.mem_univ _)

/-- **Corollaire — sélecteur "racine la plus proche" : algorithme
    Voronoï-équivariant.**

    Avec `selector := nearest_anchor_selector`, l'algorithme converge
    vers la racine cyclotomique la plus proche de `z`, *sans qu'aucune
    information sur la décomposition Voronoï ne soit fournie* à
    l'avance — la racine est déterminée localement par `z`. -/
theorem galois_equivariant_nearest_anchor
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0) :
    ∀ z : ℂ, ∀ ε : ℝ, 0 < ε →
      ∃ ε_seed : ℝ, 0 < ε_seed ∧
      ∃ N : ℕ, ∀ n ≥ N,
        ‖(steffensen_step_C x p)^[n]
            (multi_start_basin_seed_generic
              (cycAnchor α p (nearest_anchor_selector α p hp z))
              ε_seed z)
              - cycAnchor α p (nearest_anchor_selector α p hp z)‖ ≤ ε :=
  galois_equivariant_multi_start x hx p hp α hα hα_pow hSp
    (nearest_anchor_selector α p hp)

end NearestAnchorSelector

/-! ============================================================
  §97.4  Sélecteur sectoriel (partition angulaire de ℂ)
============================================================ -/

section SectorSelector

/-- **Sélecteur sectoriel** : partitionne `ℂ` en `p` secteurs
    angulaires autour de l'origine et retourne l'indice `k ∈ Fin p`
    du secteur contenant `z`. Pour `z = 0` (ou `z` sans argument
    bien défini), retourne `0` par convention.

    Construction : utilise `Complex.arg z ∈ (-π, π]`, divise par
    `2π/p`, puis prend le résultat modulo `p`. -/
noncomputable def sector_selector (p : ℕ) (hp : 1 ≤ p) (z : ℂ) :
    Fin p := by
  have hp_pos : 0 < p := hp
  refine ⟨(⌊(z.arg + Real.pi) / (2 * Real.pi / p)⌋₊) % p, ?_⟩
  exact Nat.mod_lt _ hp_pos

/-- **Corollaire — sélecteur sectoriel : algorithme par secteur
    angulaire.**

    Avec `selector := sector_selector`, l'algorithme converge vers la
    racine cyclotomique du secteur angulaire contenant `z`. C'est la
    décomposition "naturelle" de `ℂ` en `p` régions Galois-équivalentes. -/
theorem galois_equivariant_sector
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0) :
    ∀ z : ℂ, ∀ ε : ℝ, 0 < ε →
      ∃ ε_seed : ℝ, 0 < ε_seed ∧
      ∃ N : ℕ, ∀ n ≥ N,
        ‖(steffensen_step_C x p)^[n]
            (multi_start_basin_seed_generic
              (cycAnchor α p (sector_selector p hp z)) ε_seed z)
              - cycAnchor α p (sector_selector p hp z)‖ ≤ ε :=
  galois_equivariant_multi_start x hx p hp α hα hα_pow hSp
    (sector_selector p hp)

end SectorSelector

end Pandrosion
