/-
  Universitas Pandrosion — §87. **Stub : série formelle de Böttcher
  pour σ à `x = 2, p = 3`.**

  Cible originale : tenter une version formelle de `ψ` (§68
  `BottcherCoordinateP3X2`) via série de puissances tronquée,
  preuve de convergence locale.

  **Constat** : la formalisation Mathlib de séries holomorphes
  convergentes (`HasFPowerSeriesAt` etc.) est non triviale et hors
  scope d'une session autonome.

  **Cette session ajoute des cadres Prop** pour les coefficients de
  Böttcher tronqués (Taylor série d'ordre `n`) :

    `ψ_n(z) = ∑_{k=1}^{n} c_k · (z − α₀)^k`

  satisfait `ψ_n(σ(z)) = ψ_n(z)² + O((z − α₀)^{n+1})`.

  Coefficient leading `c_1 = σ''(α₀)/2` (calcul classique pour ancre
  super-attractive de degré 2).

  Contents.

    §87.1  `BottcherTaylorCoefficient k` — Prop pour coeff k-ième.
    §87.2  `BottcherTruncatedExists n` — existence de série tronquée
           d'ordre n.
-/

import Pandrosion.Core.BottcherP3Framework

namespace Pandrosion

open Complex

/-! §87.1  Conjecture coefficient leading `c_1` -/

/-- **Conjecture coefficient leading `c_1` de la série de Böttcher**
    à `x = 2, p = 3`. Relié à la dérivée seconde de σ au point fixe
    par `c_1 = σ''(α₀)/2`. -/
def BottcherLeadingCoefficient (c : ℂ) : Prop :=
  ∃ (ψ : ℂ → ℂ) (r : ℝ), 0 < r ∧ BottcherEquationP3X2 ψ r ∧
    ∀ z : ℂ, ‖z - ((alphaX2 : ℝ) : ℂ)‖ < r →
      ‖ψ z - c * (z - ((alphaX2 : ℝ) : ℂ))‖
        ≤ ‖z - ((alphaX2 : ℝ) : ℂ)‖ ^ 2

/-! §87.2  Existence de série tronquée d'ordre `n` -/

/-- **Conjecture** : pour tout ordre `n`, il existe une série de
    Böttcher tronquée à l'ordre `n` satisfaisant l'équation
    fonctionnelle modulo `O(‖z - α₀‖^{n+1})`. -/
def BottcherTruncatedExists (n : ℕ) : Prop :=
  ∃ (ψ : ℂ → ℂ) (r : ℝ), 0 < r ∧
    ψ ((alphaX2 : ℝ) : ℂ) = 0 ∧
    ∀ z : ℂ, ‖z - ((alphaX2 : ℝ) : ℂ)‖ < r →
      ‖ψ (steffensen_step_C (2 : ℂ) 3 z) - (ψ z)^2‖
        ≤ ‖z - ((alphaX2 : ℝ) : ℂ)‖ ^ (n + 1)

/-! §87.3  Implication `BottcherCoordinateP3X2` ⟹ truncated -/

/-- Si la coordonnée Böttcher complète existe, alors la série
    tronquée d'ordre `n` existe trivialement (avec error = 0). -/
theorem bottcher_truncated_of_full
    (hB : BottcherCoordinateP3X2) (n : ℕ) :
    BottcherTruncatedExists n := by
  obtain ⟨ψ, r, hr, hEq⟩ := hB
  refine ⟨ψ, r, hr, hEq.1, fun z hz => ?_⟩
  rw [hEq.2 z hz]
  rw [sub_self, norm_zero]
  positivity

end Pandrosion
