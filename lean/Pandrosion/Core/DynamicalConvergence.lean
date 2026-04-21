/-
  Universitas Pandrosion — §28. **Dynamical convergence almost-everywhere.**

  The dynamical-systems refinement of §19 McMullen and §18 Voronoï:
  for Lebesgue-almost every initial point `z₀ ∈ ℂ`, the Pandrosion-
  Steffensen orbit `(z_k)_k := (steffensen_step_C x p)^[k] z₀`
  *converges in ℂ* to one of the `p` cyclotomic anchors
  `γ_0, …, γ_{p-1}`.

  This upgrades the static "a.e. unique nearest anchor" statement
  (§18.3 `voronoiBoundary_volume_zero`) to an honest dynamical
  statement about the trajectory.

  Pieces stitched:

    • §27 `steffensen_local_attraction` — local basin around each
      anchor, conditional on a local quadratic bound per anchor
      (the super-attracting content of Steffensen's method).

    • A.e. trajectory-entry hypothesis — for a.e. `z₀`, the orbit
      eventually enters one of the local disks `B(γ_s, r_s)`. This
      is the global dynamical content (a Fatou-type McMullen
      statement at the measure level), taken here as a hypothesis.

  Content.

    §28.1  `tendsto_iterate_of_tail` — tail-shift invariance of
           `Tendsto` for iteration: a limit on the shifted orbit
           lifts to a limit on the full orbit.

    §28.2  `steffensen_dynamical_convergence_pointwise` — pointwise
           statement at a single `z₀`: orbit enters `B(α, r)` at some
           time `k₀` ⇒ orbit converges to `α`.

    §28.3  ★★★ `steffensen_dynamical_convergence_ae` — the crown
           jewel: almost-everywhere dynamical convergence of the
           Steffensen orbit to a cyclotomic anchor, conditional on
           per-anchor local quadratic bounds and a.e. trajectory
           entry.
-/

import Pandrosion.Core.LocalAttraction
import Pandrosion.Core.CyclotomicMcMullen

namespace Pandrosion

open Filter Topology MeasureTheory Complex

/-! ============================================================
  §28.1  Tail-shift invariance of iteration limits
============================================================ -/

section TailShift

/-- Iteration limits are invariant under a finite initial shift:
    if the shifted orbit `(f^[k] (f^[k₀] z₀))_k` converges to `α`,
    then so does the full orbit `(f^[k] z₀)_k`. -/
theorem tendsto_iterate_of_tail
    {f : ℂ → ℂ} {α : ℂ} (z₀ : ℂ) (k₀ : ℕ)
    (h_tail : Tendsto (fun k => f^[k] (f^[k₀] z₀)) atTop (𝓝 α)) :
    Tendsto (fun k => f^[k] z₀) atTop (𝓝 α) := by
  have h_add : (fun k => f^[k] (f^[k₀] z₀)) = (fun k => f^[k + k₀] z₀) := by
    funext k; rw [← Function.iterate_add_apply]
  rw [h_add] at h_tail
  exact (tendsto_add_atTop_iff_nat k₀).mp h_tail

end TailShift

/-! ============================================================
  §28.2  Pointwise dynamical convergence from trajectory entry
============================================================ -/

section PointwiseConvergence

/-- **Pointwise dynamical convergence.** If at some time `k₀` the
    Steffensen orbit of `z₀` enters the local basin `B(α, r)` of a
    Pandrosion fixed point `α`, and the map satisfies a local
    quadratic contraction bound on that basin, then the full orbit
    converges to `α`. -/
theorem steffensen_dynamical_convergence_pointwise
    (x : ℂ) (hx : x ≠ 0) (p : ℕ) (α : ℂ)
    (hSp : Sp_C p α ≠ 0) (hfp : α ^ p = 1 / x)
    {r K : ℝ} (hr : 0 < r) (hK_pos : 0 < K) (h_small : K * r ≤ 1 / 2)
    (h_quadratic : ∀ z : ℂ, ‖z - α‖ < r →
      ‖steffensen_step_C x p z - α‖ ≤ K * ‖z - α‖ ^ 2)
    (z₀ : ℂ) (k₀ : ℕ)
    (h_enter : ‖(steffensen_step_C x p)^[k₀] z₀ - α‖ < r) :
    Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop (𝓝 α) := by
  have h_local := steffensen_local_attraction x hx p α hSp hfp
      hr hK_pos h_small h_quadratic
  exact tendsto_iterate_of_tail z₀ k₀
      (h_local ((steffensen_step_C x p)^[k₀] z₀) h_enter)

end PointwiseConvergence

/-! ============================================================
  §28.3  ★★★ steffensen_dynamical_convergence_ae — the crown jewel
============================================================ -/

section DynamicalConvergenceAE

/-- **★★★ Steffensen dynamical convergence a.e.**

  Fix `p ≥ 1`, a real `x > 0` (promoted to `ℂ`, `x ≠ 0`), and a
  complex `p`-th root `α` of `1/x`. For each cyclotomic anchor
  `γ_s := cycAnchor α p s`, assume:

    • `Sp_C p γ_s ≠ 0` (non-degeneracy).

    • A local quadratic contraction bound on `B(γ_s, r_s)`:
          `‖steffensen_step_C x p z − γ_s‖ ≤ K_s · ‖z − γ_s‖²`
      with `K_s · r_s ≤ 1/2` (the super-attracting content of
      Steffensen's method; proven analytically in the literature,
      taken as hypothesis here).

    • A.e. trajectory entry: for a.e. `z₀`, the Steffensen orbit
      eventually enters some local basin `B(γ_s, r_s)`.

  Conclusion. For Lebesgue-almost every `z₀ ∈ ℂ`, there exists an
  anchor index `s` such that
      `Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop (𝓝 γ_s)`,
  i.e., the Steffensen orbit converges *in ℂ* to one of the
  cyclotomic `p`-th roots of `1/x`.

  **Paper-citable consequence.** *Pandrosion-Steffensen almost-
  everywhere solves `z^p = x` in the dynamical sense: the iterates
  converge, not merely live in the correct basin.* -/
theorem steffensen_dynamical_convergence_ae
    (p : ℕ) (hp : 1 ≤ p)
    (x : ℂ) (hx : x ≠ 0)
    (α : ℂ) (hα_pow : α ^ p = 1 / x)
    (hSp : ∀ s : Fin p, Sp_C p (cycAnchor α p s) ≠ 0)
    (r : Fin p → ℝ) (K : Fin p → ℝ)
    (hr : ∀ s, 0 < r s) (hK_pos : ∀ s, 0 < K s)
    (h_small : ∀ s, K s * r s ≤ 1 / 2)
    (h_quadratic : ∀ s : Fin p, ∀ z : ℂ,
      ‖z - cycAnchor α p s‖ < r s →
      ‖steffensen_step_C x p z - cycAnchor α p s‖ ≤
        K s * ‖z - cycAnchor α p s‖ ^ 2)
    (h_entry : ∀ᵐ z₀ : ℂ ∂volume,
      ∃ s : Fin p, ∃ k₀ : ℕ,
        ‖(steffensen_step_C x p)^[k₀] z₀ - cycAnchor α p s‖ < r s) :
    ∀ᵐ z₀ : ℂ ∂volume,
      ∃ s : Fin p,
        Tendsto (fun k => (steffensen_step_C x p)^[k] z₀) atTop
          (𝓝 (cycAnchor α p s)) := by
  filter_upwards [h_entry] with z₀ hz₀
  obtain ⟨s, k₀, h_enter⟩ := hz₀
  refine ⟨s, ?_⟩
  have hfp_s : (cycAnchor α p s) ^ p = 1 / x := by
    rw [cycAnchor_pow α hp s]; exact hα_pow
  exact steffensen_dynamical_convergence_pointwise x hx p (cycAnchor α p s)
    (hSp s) hfp_s (hr s) (hK_pos s) (h_small s) (h_quadratic s)
    z₀ k₀ h_enter

end DynamicalConvergenceAE

end Pandrosion
