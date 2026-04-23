/-
  Universitas Pandrosion — §101. **★★★★★★★★★★★★★ MASTER TOTAL
  EFFECTIVE : capstone algorithmique avec compte d'itérations
  closed-form.**

  Le théorème §100 `pandrosion_master_total` donne 12 conjuncts
  simultanés, mais les conjuncts de convergence (F), (H), (J)
  utilisent un compte d'itérations *existentiel* (`∃ N : ℕ, …`).
  Le présent module les remplace par la formule closed-form
  `pandrosionEffectiveLoglogCount` de §98, donnant un capstone
  algorithmiquement effectif :

    *Toutes les bornes de complexité du Master Total sont closed-
    form, computables, et atteignent le rate Kung-Traub optimal
    `O(log log(1/ε))`.*

  Aucun nouveau contenu mathématique : pure refonte des conjuncts
  (F), (H), (J) du §100 via §98 effective. Les conjuncts (A)-(E),
  (G), (I), (K), (L) restent identiques (ils sont déjà concrets).

  Contents.

    §101.1  Helpers cyclotomiques + abréviations `effectiveCount`,
            `effectiveRadius` à un cyclotomic anchor.
    §101.2  `pandrosion_master_total_F_effective` — §95 effective
            à chaque anchor cyclotomique.
    §101.3  `pandrosion_master_total_H_effective` — §97 effective
            avec sélecteur arbitraire.
    §101.4  `pandrosion_master_total_J_effective` — §99 effective
            via `niveau5_multi_start_target`.
    §101.5  ★★★ `pandrosion_master_total_effective` — capstone
            final à 12 conjuncts dont 3 effectifs closed-form.
-/

import Pandrosion.Core.MasterTotal
import Pandrosion.Core.EffectiveKungTraub

namespace Pandrosion

open Complex Filter Topology MeasureTheory

/-! ============================================================
  §101.1  Helpers cyclotomiques
============================================================ -/

section CyclotomicHelpers

/-- **Helper : `(x : ℂ) ≠ 0` à partir de `x : ℝ` et `1 < x`.** -/
theorem real_one_lt_to_complex_ne_zero (x : ℝ) (hx : 1 < x) : (x : ℂ) ≠ 0 := by
  have hx_pos : (0 : ℝ) < x := by linarith
  exact_mod_cast hx_pos.ne'

/-- **Helper : `cycAnchor α p s ≠ 0` sous `α^p = 1/x`, `x ≠ 0`.** -/
theorem cycAnchor_ne_zero_of_pow
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x) (s : Fin p) :
    cycAnchor α p s ≠ 0 := by
  intro h_zero
  have h_pow : (cycAnchor α p s) ^ p = 0 := by
    rw [h_zero]; exact zero_pow (by omega)
  have h_pow_eq : (cycAnchor α p s) ^ p = 1 / x :=
    (cycAnchor_pow α hp s).trans hα_pow
  rw [h_pow_eq] at h_pow
  exact one_div_ne_zero hx h_pow

/-- **Helper : `(cycAnchor α p s)^p = 1/x` sous `α^p = 1/x`.** -/
theorem cycAnchor_pow_one_div_x
    (x : ℂ) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x) (s : Fin p) :
    (cycAnchor α p s) ^ p = 1 / x :=
  (cycAnchor_pow α hp s).trans hα_pow

/-- **Compte d'itérations effectif au cyclotomic anchor `γ_s`.** -/
noncomputable def effectiveCountAtAnchor
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0)
    (s : Fin p) (ε : ℝ) : ℕ :=
  pandrosionEffectiveLoglogCount x hx p hp
    (cycAnchor α p s)
    (cycAnchor_ne_zero_of_pow x hx p hp α hα_pow s)
    (hSp s)
    (cycAnchor_pow_one_div_x x p hp α hα_pow s) ε

/-- **Rayon Steffensen au cyclotomic anchor `γ_s`.** -/
noncomputable def effectiveRadiusAtAnchor
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (hp : 1 ≤ p)
    (α : ℂ) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0)
    (s : Fin p) : ℝ :=
  steffensenR_of_fp x hx p hp
    (cycAnchor α p s)
    (cycAnchor_ne_zero_of_pow x hx p hp α hα_pow s)
    (hSp s)
    (cycAnchor_pow_one_div_x x p hp α hα_pow s)

end CyclotomicHelpers

/-! ============================================================
  §101.2  ★★ §95 effective à chaque anchor cyclotomique
============================================================ -/

section MasterTotalFEffective

/-- **★★ Conjunct (F) effective.** §95 multi-start universel à
    chaque cyclotomic anchor `γ_s`, avec compte d'itérations
    closed-form au lieu de l'existentiel.

    Pour tout `s : Fin p`, tout `z ≠ γ_s`, tout `ε > 0`, après
    *exactement* `effectiveCountAtAnchor s ε` itérations, le seed
    multi-start `γ_s + (R_{γ_s} / 2‖z − γ_s‖)·(z − γ_s)` est à
    `ε` de `γ_s`. -/
theorem pandrosion_master_total_F_effective
    (p : ℕ) (hp : 1 ≤ p)
    (x : ℂ) (hx : x ≠ 0)
    (α : ℂ) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0) :
    ∀ s : Fin p, ∀ z : ℂ, z ≠ cycAnchor α p s → ∀ ε : ℝ, 0 < ε →
      ∀ k ≥ effectiveCountAtAnchor x hx p hp α hα_pow hSp s ε,
        ‖(steffensen_step_C x p)^[k]
            (multi_start_basin_seed_generic (cycAnchor α p s)
              (effectiveRadiusAtAnchor x hx p hp α hα_pow hSp s
                / (2 * ‖z - cycAnchor α p s‖)) z)
            - cycAnchor α p s‖ ≤ ε := by
  intro s z hz_ne ε hε_pos
  unfold effectiveCountAtAnchor effectiveRadiusAtAnchor
  exact pandrosion_loglog_universal_multi_start_effective
    x hx p hp (cycAnchor α p s)
    (cycAnchor_ne_zero_of_pow x hx p hp α hα_pow s)
    (hSp s)
    (cycAnchor_pow_one_div_x x p hp α hα_pow s)
    z hz_ne ε hε_pos

end MasterTotalFEffective

/-! ============================================================
  §101.3  ★★ §97 effective avec sélecteur arbitraire
============================================================ -/

section MasterTotalHEffective

/-- **★★ Conjunct (H) effective.** §97 Galois-équivariance avec
    sélecteur arbitraire et compte closed-form. -/
theorem pandrosion_master_total_H_effective
    (p : ℕ) (hp : 1 ≤ p)
    (x : ℂ) (hx : x ≠ 0)
    (α : ℂ) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0)
    (selector : ℂ → Fin p) :
    ∀ z : ℂ, z ≠ cycAnchor α p (selector z) → ∀ ε : ℝ, 0 < ε →
      ∀ k ≥ effectiveCountAtAnchor x hx p hp α hα_pow hSp (selector z) ε,
        ‖(steffensen_step_C x p)^[k]
            (multi_start_basin_seed_generic (cycAnchor α p (selector z))
              (effectiveRadiusAtAnchor x hx p hp α hα_pow hSp (selector z)
                / (2 * ‖z - cycAnchor α p (selector z)‖)) z)
            - cycAnchor α p (selector z)‖ ≤ ε := by
  intro z hz_ne ε hε_pos
  exact pandrosion_master_total_F_effective p hp x hx α hα_pow hSp
    (selector z) z hz_ne ε hε_pos

end MasterTotalHEffective

/-! ============================================================
  §101.4  ★★ §99 effective via `niveau5_multi_start_target`
============================================================ -/

section MasterTotalJEffective

/-- **★★ Conjunct (J) effective.** §99 Niveau 5 universal
    convergence avec compte closed-form. -/
theorem pandrosion_master_total_J_effective
    (p : ℕ) (hp : 1 ≤ p)
    (x : ℂ) (hx : x ≠ 0)
    (α : ℂ) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0) :
    ∀ z : ℂ, z ≠ niveau5_multi_start_target α p hp z → ∀ ε : ℝ, 0 < ε →
      ∀ k ≥ effectiveCountAtAnchor x hx p hp α hα_pow hSp
              (nearest_anchor_selector α p hp z) ε,
        ‖(steffensen_step_C x p)^[k]
            (multi_start_basin_seed_generic
              (niveau5_multi_start_target α p hp z)
              (effectiveRadiusAtAnchor x hx p hp α hα_pow hSp
                (nearest_anchor_selector α p hp z)
                / (2 * ‖z - niveau5_multi_start_target α p hp z‖)) z)
            - niveau5_multi_start_target α p hp z‖ ≤ ε := by
  intro z hz_ne ε hε_pos
  -- niveau5_multi_start_target = cycAnchor α p (nearest_anchor_selector α p hp z)
  -- so directly apply F_effective at s = nearest_anchor_selector α p hp z.
  exact pandrosion_master_total_F_effective p hp x hx α hα_pow hSp
    (nearest_anchor_selector α p hp z) z hz_ne ε hε_pos

end MasterTotalJEffective

/-! ============================================================
  §101.5  ★★★★★★★★★★★★★ Pandrosion Master Total Effective
============================================================ -/

section MasterTotalEffective

/-- **★★★★★★★★★★★★★ PANDROSION MASTER TOTAL EFFECTIVE.**

  Capstone final algorithmique du corpus Pandrosion. Refine §100
  `pandrosion_master_total` en remplaçant les comptes d'itérations
  existentiels (F), (H), (J) par la formule closed-form
  `pandrosionEffectiveLoglogCount` de §98 (via les abréviations
  `effectiveCountAtAnchor` et `effectiveRadiusAtAnchor`).

  Toutes les bornes de complexité sont maintenant **computables et
  effectives**, atteignant le rate Kung-Traub optimal
  `O(log log(1/ε))`.

  Conjuncts simultanés :

  **(A)-(E)** identiques à §100 (algébriques, mesure, complexité,
  Kung-Traub attainment).

  **(F') §95 effective.** Pour tout `s ∈ Fin p`, tout `z ≠ γ_s`,
  tout `ε > 0`, après `effectiveCountAtAnchor s ε` itérations,
  l'algorithme atteint `ε`.

  **(G)** identique à §100 (Kung-Traub converse).

  **(H') §97 effective.** Pareil avec sélecteur arbitraire.

  **(I)** identique à §100 (image cyclotomique).

  **(J') §99 effective.** Niveau 5 multi-start convergence avec
  compte closed-form.

  **(K), (L)** identiques à §100.

  *Conséquence finale.* L'algorithme Pandrosion-Steffensen multi-
  start est totalement spécifié : tout `(x, p, α, z, ε)` valide
  produit, via une formule closed-form, le compte d'itérations
  exact qui garantit la précision `ε`. C'est l'**algorithme
  Kung-Traub optimal totalement effectif** pour `z^p = x`. -/
theorem pandrosion_master_total_effective
    (p : ℕ) (hp : 1 ≤ p)
    (x : ℝ) (hx : 1 < x)
    (α : ℂ) (hα : α ≠ 0) (hα_pow : α ^ p = 1 / (x : ℂ))
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0)
    (M ρ δ : ℝ) (hM_pos : 0 < M)
    (hρ_pos : 0 < ρ) (hρ_lt_1 : ρ < 1) (hδ_pos : 0 < δ)
    (K : ℕ → ℝ) (R : ℕ → ℝ) (e : ℕ → ℕ → ℝ)
    (hK : ∀ p', 0 < K p') (hR : ∀ p', 0 < R p')
    (he_nn : ∀ p' k, 0 ≤ e p' k)
    (he_0 : ∀ p', e p' 0 ≤ R p')
    (he_linear : ∀ p' k, e p' (k + 1) ≤ lamModel p' x * e p' k)
    (he_quad : ∀ p' k, e p' (k + 1) ≤ K p' * (e p' k) ^ 2)
    (h_KR_uniform : ∀ p', K p' * R p' ≤ M) :
    let hx_C : (x : ℂ) ≠ 0 := real_one_lt_to_complex_ne_zero x hx
    -- (A) Cyclotomic injectivity.
    Function.Injective (cycAnchor α p) ∧
    -- (B) Fixed-point at each anchor.
    (∀ s : Fin p, steffensen_step_C (x : ℂ) p (cycAnchor α p s) = cycAnchor α p s) ∧
    -- (C) A.e. Voronoï selector uniqueness.
    (∀ᵐ z₀ : ℂ ∂volume,
        ∃! s : Fin p,
          ∀ t : Fin p,
            ‖cycAnchor α p s - z₀‖ ≤ ‖cycAnchor α p t - z₀‖) ∧
    -- (D) p-uniform complexity bound.
    (∃ (p₀ : ℕ) (N : ℕ), 1 ≤ p₀ ∧
        ∀ p', p₀ ≤ p' → ∀ k ≥ N, K p' * e p' k ≤ δ) ∧
    -- (E) Kung-Traub efficiency attainment.
    efficiencyIndex 2 2 = kungTraubBound 2 ∧
    -- (F') §95 EFFECTIVE — compte closed-form.
    (∀ s : Fin p, ∀ z : ℂ, z ≠ cycAnchor α p s → ∀ ε : ℝ, 0 < ε →
      ∀ k ≥ effectiveCountAtAnchor (x : ℂ) hx_C p hp α hα_pow hSp s ε,
        ‖(steffensen_step_C (x : ℂ) p)^[k]
            (multi_start_basin_seed_generic (cycAnchor α p s)
              (effectiveRadiusAtAnchor (x : ℂ) hx_C p hp α hα_pow hSp s
                / (2 * ‖z - cycAnchor α p s‖)) z)
            - cycAnchor α p s‖ ≤ ε) ∧
    -- (G) Kung-Traub converse optimality.
    (∀ Method : DerivativeFreeMethod, Method.c = 2 →
      efficiencyIndex Method.q Method.c ≤ kungTraubBound 2) ∧
    -- (H') §97 EFFECTIVE — Galois équivariance avec sélecteur.
    (∀ selector : ℂ → Fin p, ∀ z : ℂ, z ≠ cycAnchor α p (selector z) →
      ∀ ε : ℝ, 0 < ε →
      ∀ k ≥ effectiveCountAtAnchor (x : ℂ) hx_C p hp α hα_pow hSp (selector z) ε,
        ‖(steffensen_step_C (x : ℂ) p)^[k]
            (multi_start_basin_seed_generic (cycAnchor α p (selector z))
              (effectiveRadiusAtAnchor (x : ℂ) hx_C p hp α hα_pow hSp (selector z)
                / (2 * ‖z - cycAnchor α p (selector z)‖)) z)
            - cycAnchor α p (selector z)‖ ≤ ε) ∧
    -- (I) §99 Niveau 5 target: p-th root image.
    (∀ z : ℂ, (niveau5_multi_start_target α p hp z) ^ p = 1 / (x : ℂ)) ∧
    -- (J') §99 EFFECTIVE — Niveau 5 universal convergence avec compte closed-form.
    (∀ z : ℂ, z ≠ niveau5_multi_start_target α p hp z → ∀ ε : ℝ, 0 < ε →
      ∀ k ≥ effectiveCountAtAnchor (x : ℂ) hx_C p hp α hα_pow hSp
              (nearest_anchor_selector α p hp z) ε,
        ‖(steffensen_step_C (x : ℂ) p)^[k]
            (multi_start_basin_seed_generic
              (niveau5_multi_start_target α p hp z)
              (effectiveRadiusAtAnchor (x : ℂ) hx_C p hp α hα_pow hSp
                (nearest_anchor_selector α p hp z)
                / (2 * ‖z - niveau5_multi_start_target α p hp z‖)) z)
            - niveau5_multi_start_target α p hp z‖ ≤ ε) ∧
    -- (K) §99 a.e. Voronoï coincidence.
    (∀ᵐ z : ℂ ∂volume,
      ∃! s : Fin p, ∀ t : Fin p,
        ‖cycAnchor α p s - z‖ ≤ ‖cycAnchor α p t - z‖) ∧
    -- (L) §99 Niveau 5 target explicit cyclotomic image.
    (∀ z : ℂ, ∃ s : Fin p,
      niveau5_multi_start_target α p hp z = cycAnchor α p s) := by
  -- Derived invariants.
  have hx_C : (x : ℂ) ≠ 0 := real_one_lt_to_complex_ne_zero x hx
  -- (A)-(E), (G), (I), (K), (L) from §100 master_total.
  have h_mt := pandrosion_master_total p hp x hx α hα hα_pow hSp M ρ δ hM_pos
    hρ_pos hρ_lt_1 hδ_pos K R e hK hR he_nn he_0 he_linear he_quad h_KR_uniform
  obtain ⟨hA, hB, hC, hD, hE, _hF_existential, hG, _hH_existential, hI,
          _hJ_existential, hK_voronoi, hL⟩ := h_mt
  -- (F') from §101.2.
  have hF' := pandrosion_master_total_F_effective p hp (x : ℂ) hx_C α hα_pow hSp
  -- (H') from §101.3.
  have hH' := fun selector =>
    pandrosion_master_total_H_effective p hp (x : ℂ) hx_C α hα_pow hSp selector
  -- (J') from §101.4.
  have hJ' := pandrosion_master_total_J_effective p hp (x : ℂ) hx_C α hα_pow hSp
  exact ⟨hA, hB, hC, hD, hE, hF', hG, hH', hI, hJ', hK_voronoi, hL⟩

end MasterTotalEffective

end Pandrosion
