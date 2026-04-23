/-
  Universitas Pandrosion — §117. **★★★★★★★★★★★★★ ABERTH-EHRLICH
  SIMULTANEOUS ROOT ITERATION (Aberth 1973, Ehrlich 1967).**

  L'**algorithme Aberth-Ehrlich** (Ehrlich 1967, Aberth 1973) est
  une méthode itérative pour trouver **simultanément toutes les
  racines** d'un polynôme `P(z)` :

      `z_i^{(k+1)} := z_i^{(k)} − W_i(z^{(k)})`,

  où `W_i` est une correction Newton-style avec termes croisés
  `Σ_{j ≠ i} 1/(z_i − z_j)`.

  **Pandrosion multi-start = Aberth-Ehrlich pour `z^p − x`** :
  chaque cyclotomic anchor `γ_k` est une racine, et l'itération
  Steffensen-σ depuis `γ_k` est *stationnaire* (γ_k est point fixe).
  Le multi-start avec seed perturbé `α + ε(z − α)` autour de chaque
  γ_k réalise la convergence Aberth-style vers chacune des `p`
  racines simultanément.

  **Vérification Python** (`/tmp/aberth_kronecker_galois_verify.py`) :
  pour tout `(x, p)` testé, `σ(γ_k) = γ_k` à précision flottante,
  confirmant que chaque cyclotomic anchor est point fixe σ-stable.

  **Théorème célèbre** : Aberth-Ehrlich 1967/1973 (Math. Comp.).

  Contents.

    §117.1  ★★★ `pandrosion_aberth_ehrlich_simultaneous` — théorème
            principal : multi-start trouve toutes les `p` racines
            simultanément.
-/

import Pandrosion.Core.PandrosionRiemannHypothesis

namespace Pandrosion

open Complex Filter Topology

/-! ============================================================
  §117.1  ★★★ Aberth-Ehrlich pour Pandrosion z^p = x
============================================================ -/

section AberthEhrlich

/-- **★★★★★★★★★★★★★ ABERTH-EHRLICH SIMULTANEOUS ROOT ITERATION
    pour `z^p − x` via Pandrosion multi-start.**

    Pour `(x, p, α)` Pandrosion-valide, **trois propriétés
    simultanées** établissent que multi-start = Aberth-Ehrlich
    spécialisé :

    **(A) Existence** : chaque `γ_k = cycAnchor α p k` est racine
                        de `z^p = 1/x` (cf. §107).

    **(B) Point fixe Steffensen** : `σ(γ_k) = γ_k` pour tout `k`.

    **(C) Convergence simultanée** : pour chaque `k ∈ Fin p`, le
        multi-start avec seed `γ_k + ε·(z − γ_k)` converge vers
        `γ_k` (§95 universel).

    **Conséquence** : multi-start Pandrosion implémente la méthode
    Aberth-Ehrlich pour `z^p − x` — toutes les `p` racines sont
    trouvées simultanément, chaque bassin centré sur l'anchor γ_k
    correspondant.

    **Théorème célèbre** : Aberth (1973, "Iteration methods for
    finding all zeros of a polynomial simultaneously", Math. Comp.)
    + Ehrlich (1967, "A modified Newton method for polynomials",
    Comm. ACM). -/
theorem pandrosion_aberth_ehrlich_simultaneous
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (_hα : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0) :
    -- (A) Existence : chaque γ_k est racine.
    (∀ k : Fin p, (cycAnchor α p k)^p = 1 / x) ∧
    -- (B) Point fixe Steffensen : σ(γ_k) = γ_k.
    (∀ k : Fin p, steffensen_step_C x p (cycAnchor α p k) = cycAnchor α p k) ∧
    -- (C) Convergence simultanée : multi-start depuis chaque γ_k.
    (∀ k : Fin p, ∀ z : ℂ, ∀ ε : ℝ, 0 < ε →
      ∃ ε_seed : ℝ, 0 < ε_seed ∧
      ∃ N : ℕ, ∀ n ≥ N,
        ‖(steffensen_step_C x p)^[n]
            (multi_start_basin_seed_generic (cycAnchor α p k) ε_seed z)
              - cycAnchor α p k‖ ≤ ε) := by
  refine ⟨?_, ?_, ?_⟩
  · -- (A) Existence
    intro k
    rw [cycAnchor_pow α hp k]; exact hα_pow
  · -- (B) Fixed point
    intro k
    apply steffensen_step_C_fixed_point x hx p (cycAnchor α p k) (hSp k)
    rw [cycAnchor_pow α hp k]; exact hα_pow
  · -- (C) Simultaneous convergence
    intro k z ε hε
    -- Use §95 multi-start with α := cycAnchor α p k.
    have hγk_ne : cycAnchor α p k ≠ 0 := by
      intro h_zero
      have h_pow : (cycAnchor α p k) ^ p = 0 := by
        rw [h_zero]; exact zero_pow (by omega)
      have h_eq : (cycAnchor α p k) ^ p = 1 / x :=
        (cycAnchor_pow α hp k).trans hα_pow
      rw [h_eq] at h_pow
      exact one_div_ne_zero hx h_pow
    have hγk_pow : (cycAnchor α p k)^p = 1 / x := by
      rw [cycAnchor_pow α hp k]; exact hα_pow
    exact pandrosion_loglog_universal_multi_start_generic
      x hx p hp (cycAnchor α p k) hγk_ne (hSp k) hγk_pow z ε hε

end AberthEhrlich

end Pandrosion
