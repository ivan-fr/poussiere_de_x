/-
  Universitas Pandrosion — §60. **Universal Pandrosion convergence
  conjecture on ℂ.**

  Le corpus établit :
    (a) p = 2 : convergence a.e. inconditionnelle vers `{α, −α}` (§40).
    (b) p = 3 : convergence a.e. **conditionnelle** modulo Julia-mesure-
        nulle (§41).
    (c) p ≥ 4 : non formalisé à ce jour.

  Ce module **cristallise** la conjecture universelle comme `Prop`
  explicite, chaîne les cas prouvés, et laisse p ≥ 4 comme hypothèse
  portée (à discharger dans un développement futur).

  Trois formes de plus en plus fortes :

    * **faible (Pandrosion-McMullen universel)** : `∀ p ≥ 2,
      a.e. convergence vers une racine cyclotomique` — exactement
      `McMullenAEEntry p x α` vue comme conjecture universelle.
    * **forte (bassin principal global)** : `a.e. convergence vers
      α` (racine principale) — motivée par l'observation empirique
      Python (bassin de γ₀ ≫ bassins de γ_s, s ≠ 0).
    * **ultra-forte (uniformité en p)** : taux loglog uniforme en p
      (lien avec §22.4).

  Contents.

    §60.1  `UniversalPandrosionConvergence` — conjecture faible,
           forme `Prop`.

    §60.2  `StrongPandrosionConvergence` — conjecture forte
           (bassin principal).

    §60.3  `universal_implies_some_root` — réciproque triviale.

    §60.4  `universal_convergence_p2` — cas p=2 prouvé inconditionnel.

    §60.5  `universal_convergence_p3_mod_julia` — cas p=3 prouvé
           modulo Julia-mesure-nulle.

    §60.6  `pandrosion_universal_conjecture_statement` — l'énoncé
           formel de la conjecture universelle à p ≥ 2.
-/

import Pandrosion.Core.ComplexMcMullenP2Unconditional
import Pandrosion.Core.SteffensenP3Loglog

namespace Pandrosion

open Filter Topology MeasureTheory Complex

/-! ============================================================
  §60.1  Conjecture faible (Pandrosion-McMullen universel)
============================================================ -/

/-- **Convergence universelle faible à p donnés.**

    Pour a.e. `z₀ ∈ ℂ`, l'itérée Pandrosion-Steffensen converge vers
    *quelque* racine cyclotomique `γ_s = ω^s · α`. -/
def UniversalPandrosionConvergence (p : ℕ) (x α : ℂ) : Prop :=
  ∀ᵐ z₀ : ℂ ∂volume,
    ∃ s : Fin p,
      Filter.Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) Filter.atTop
        (nhds (cycAnchor α p s))

/-! ============================================================
  §60.2  Conjecture forte (bassin principal)
============================================================ -/

/-- **Convergence forte vers la racine principale.**

    Pour a.e. `z₀ ∈ ℂ`, l'itérée converge vers `α` spécifiquement —
    pas vers une autre racine cyclotomique. Motivé par les
    expériences numériques : le bassin de `γ₀ = α` domine
    (plus de 90 % du plan pour `x = 2, p = 3`). -/
def StrongPandrosionConvergence (p : ℕ) (x α : ℂ) : Prop :=
  ∀ᵐ z₀ : ℂ ∂volume,
    Filter.Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) Filter.atTop
      (nhds α)

/-! ============================================================
  §60.3  Implication triviale forte ⇒ faible
============================================================ -/

/-- **Forte implique faible.** Si l'itération converge vers `α`
    a.e., elle converge en particulier vers *une* racine cyclotomique
    a.e. (prendre `s = 0`). -/
theorem strong_implies_universal
    (p : ℕ) (hp : 1 ≤ p) (x α : ℂ)
    (h_strong : StrongPandrosionConvergence p x α) :
    UniversalPandrosionConvergence p x α := by
  have h_cyc0 : cycAnchor α p ⟨0, hp⟩ = α := by
    simp [cycAnchor]
  unfold StrongPandrosionConvergence UniversalPandrosionConvergence at *
  filter_upwards [h_strong] with z₀ hz₀
  refine ⟨⟨0, hp⟩, ?_⟩
  rw [h_cyc0]; exact hz₀

/-! ============================================================
  §60.4  Cas p=2 prouvé inconditionnel
============================================================ -/

/-- **★★★ Cas p=2 : conjecture universelle faible prouvée.**

    Chaîne `steffensen_p2_solves_complex_unconditional` (§40.9). -/
theorem universal_convergence_p2
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0)
    (hα_ne_neg_one : 1 + α ≠ 0) (hα_ne_one : 1 - α ≠ 0)
    (hα_pow : α ^ 2 = 1 / x)
    (hSp : ∀ s : Fin 2, Sp_C 2 (cycAnchor α 2 s) ≠ 0) :
    UniversalPandrosionConvergence 2 x α :=
  steffensen_p2_solves_complex_unconditional x hx_ne α hα_ne_zero
    hα_ne_neg_one hα_ne_one hα_pow hSp

/-! ============================================================
  §60.5  Cas p=3 prouvé conditionnel
============================================================ -/

/-- **★★★★ Cas p=3 : conjecture universelle faible modulo
    `ComplexJuliaP3MeasureZero`.**

    Chaîne `steffensen_p3_solves_complex_mod_julia_null` (§41.5). -/
theorem universal_convergence_p3_mod_julia
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0) (hα_pow : α ^ 3 = 1 / x)
    (hSp : ∀ s : Fin 3, Sp_C 3 (cycAnchor α 3 s) ≠ 0)
    (hJulia : ComplexJuliaP3MeasureZero x α) :
    UniversalPandrosionConvergence 3 x α :=
  steffensen_p3_solves_complex_mod_julia_null x hx_ne α hα_ne_zero hα_pow hSp
    hJulia

/-! ============================================================
  §60.6  Énoncé formel de la conjecture universelle
============================================================ -/

/-- **★★★★★ Conjecture universelle de Pandrosion.**

    Pour tout `p ≥ 2`, tout `x ∈ ℂ \ {0}`, tout `α` de racine `α^p = 1/x`
    avec les non-dégénérescences cyclotomiques nécessaires,
    `UniversalPandrosionConvergence p x α` est vraie.

    **Statut :**
    - `p = 2` : **prouvé** inconditionnel (§60.4).
    - `p = 3` : **conditionnel** modulo `ComplexJuliaP3MeasureZero`
      (§60.5).
    - `p ≥ 4` : **ouvert**.

    L'énoncé formel général est porté comme hypothèse `McMullenAEEntry p x α`
    (une instance de laquelle suffit pour §32.2 `steffensen_solves_ae_mod_mcmullen`
    à fournir la convergence a.e.).

    **Chaînage final (forme universelle) :**
    Si `McMullenAEEntry p x α` est disponible pour un `p`, alors
    `UniversalPandrosionConvergence p x α` l'est aussi. -/
theorem pandrosion_universal_conjecture_statement
    (p : ℕ) (hp : 1 ≤ p)
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ s : Fin p, Sp_C p (cycAnchor α p s) ≠ 0)
    (hMcM : McMullenAEEntry p x α) :
    UniversalPandrosionConvergence p x α :=
  steffensen_solves_ae_mod_mcmullen p hp x hx_ne α hα_ne_zero hα_pow hSp hMcM

end Pandrosion
