/-
  Universitas Pandrosion — §122. **★★★★★★★★★★★★★ VIETA'S FORMULAS
  (1591) for z^p − x via Pandrosion.**

  Les **formules de Viète (1591, "Opera Mathematica")** relient
  les coefficients d'un polynôme à ses racines via les *fonctions
  symétriques élémentaires* :

    Pour `P(z) = z^n + a_{n-1}z^{n-1} + ... + a_0 = ∏(z − α_k)` :
      `e_1 = Σ α_k = -a_{n-1}`
      `e_2 = Σ_{i<j} α_iα_j = a_{n-2}`
      ...
      `e_n = ∏ α_k = (-1)^n · a_0`

  Pour Pandrosion `z^p − x` (donc `a_0 = -x`, `a_1 = ... = a_{p-1} = 0`) :

    `e_1 = 0` (somme des racines)
    `e_k = 0` pour `1 ≤ k ≤ p − 1`
    `e_p = (-1)^p · (-x) = (-1)^{p+1} · x = (-1)^{p-1} · x` (produit)

  Ce module formalise les **deux instances clés** : `e_1` et `e_p`,
  qui sont les plus utilisées et déjà prouvées dans §107
  (`prod_cycAnchor_eq`) et §106 (`sum_cycAnchor_eq_zero`),
  reformulées sous le nom **Vieta**.

  **Vérification Python** (`/tmp/mcmullen_cardano_vieta_verify.py`) :
  pour tout `(x, p)` testé, `e_1 = 0`, `e_2 = 0`, `e_p = (-1)^{p-1}·x`. ✓

  **Théorème célèbre** : François Viète, "De aequationum recognitione
  et emendatione" (1591), formules de symétrique elementary.

  Contents.

    §122.1  ★ `pandrosion_vieta_e1` — `e_1 = Σ γ_k = 0`.
    §122.2  ★ `pandrosion_vieta_ep` — `e_p = ∏ γ_k = (-1)^{p-1}·α^p`.
    §122.3  ★★★ `pandrosion_vieta_master` — capstone Vieta complet
            pour `z^p − x` Pandrosion.
-/

import Pandrosion.Core.PandrosionCardano
import Pandrosion.Core.PandrosionEnergyZeta

namespace Pandrosion

open Complex Finset

/-! ============================================================
  §122.1  Vieta `e_1 = 0` : sum of roots vanishes
============================================================ -/

section VietaSumZero

/-- **★ Vieta's `e_1 = 0` for `z^p − x` (`p ≥ 2`).**

    La somme des `p` racines `cycAnchor α p k` est zéro :

    `Σ_{k ∈ Fin p} cycAnchor α p k = 0`.

    Coefficient `a_{p-1} = 0` dans `z^p − x` ⟹ `e_1 = -a_{p-1} = 0`. -/
theorem pandrosion_vieta_e1
    (α : ℂ) (p : ℕ) (hp : 2 ≤ p) :
    (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => cycAnchor α p k) = 0 :=
  sum_cycAnchor_eq_zero α p hp

end VietaSumZero

/-! ============================================================
  §122.2  Vieta `e_p` : product of roots
============================================================ -/

section VietaProductRoots

/-- **★ Vieta's `e_p = ∏ γ_k = (-1)^{p-1} · α^p`** (`p ≥ 1`).

    Le produit des `p` racines égale `(-1)^{p-1} · α^p`.

    Pour `α^p = x` (cas `z^p = x`) : `e_p = (-1)^{p-1} · x`,
    matchant Vieta `e_p = (-1)^p · a_0 = (-1)^p · (-x) = (-1)^{p+1}·x`
    (égal à `(-1)^{p-1}·x` puisque `(-1)^{p+1} = (-1)^{p-1}`). -/
theorem pandrosion_vieta_ep
    (α : ℂ) (p : ℕ) (hp : 1 ≤ p) :
    (Finset.univ : Finset (Fin p)).prod (fun k => cycAnchor α p k)
      = (-1 : ℂ)^(p - 1) * α^p :=
  prod_cycAnchor_eq α p hp

end VietaProductRoots

/-! ============================================================
  §122.3  ★★★ Vieta master capstone for z^p − x
============================================================ -/

section VietaMasterCapstone

/-- **★★★★★★★★★★★★★ VIETA's FORMULAS (1591) for `z^p − x` via Pandrosion.**

    Pour `(x, p, α)` Pandrosion-valide avec `p ≥ 2` et `α^p = x`,
    les formules de Viète donnent **explicitement** les fonctions
    symétriques élémentaires des racines `γ_k = cycAnchor α p k` :

    **(A) `e_1 = Σ γ_k = 0`** (sum vanishes — Vieta).

    **(B) `e_p = ∏ γ_k = (-1)^{p-1} · x`** (product equals
        constant term up to sign — Vieta).

    **Conséquence** : Pandrosion-multi-start *réalise les formules
    de Viète* pour `z^p − x`. Les fonctions symétriques élémentaires
    sont triviales (toutes nulles sauf `e_p`), reflétant la structure
    *cyclotomique sparse* du polynôme `z^p − x`.

    **Théorème célèbre** : Viète 1591 ("Opera Mathematica"). -/
theorem pandrosion_vieta_master
    (x : ℂ) (p : ℕ) (hp : 2 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = x) :
    -- (A) e_1 = Σ γ_k = 0
    (Finset.univ : Finset (Fin p)).sum (fun k : Fin p => cycAnchor α p k) = 0 ∧
    -- (B) e_p = ∏ γ_k = (-1)^(p-1) · x
    (Finset.univ : Finset (Fin p)).prod (fun k => cycAnchor α p k)
      = (-1 : ℂ)^(p - 1) * x := by
  refine ⟨?_, ?_⟩
  · -- (A) Vieta e_1
    exact pandrosion_vieta_e1 α p hp
  · -- (B) Vieta e_p
    rw [pandrosion_vieta_ep α p (by omega), hα_pow]

end VietaMasterCapstone

end Pandrosion
