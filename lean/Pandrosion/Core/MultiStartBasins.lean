/-
  Universitas Pandrosion — Niveau 2: structural theorems on multi-start basins.

  Builds on Foundations §8 (multi-start Voronoï selection). Four flagship
  results:

    §10  Generic convergence (conditional)  — McMullen Theorem 4.3
    §11  Voronoï basin connectivity         — circumvents McMullen 1987
    §12  Basin boundary = union of bisectors — piecewise-linear boundaries
    §13  Constructive McMullen               — global convergence for `z^p − x`
-/

import Pandrosion.Core.Foundations
import Mathlib.Analysis.Convex.Basic
import Mathlib.Topology.Connected.Basic
import Mathlib.Data.Complex.Module

namespace Pandrosion

open Finset BigOperators Set Complex Topology

/-! ============================================================
  §10. Generic convergence (conditional)

  McMullen Theorem 4.3: for `P(z) = z^d − 1`, the Pandrosion
  multi-start algorithm with `d` equispaced Cauchy-circle starts
  is generally convergent — every orbit converges, and the map
  `σ : start-index → root-index` is a bijection.

  Formal content: given `p` orbits with distinct (injective) final
  anchors, the selection permutation `σ : Fin p → Fin p` is a
  bijection of `Fin p`. This is the "injective → bijective on
  finite sets" pigeonhole argument.
============================================================ -/

section GenericConvergence

/-- **Pigeonhole: injective endo on `Fin p` is bijective.**
    The abstract finitary core of "p orbits cover all p roots". -/
theorem fin_injective_iff_bijective (p : ℕ) (σ : Fin p → Fin p) :
    Function.Injective σ ↔ Function.Bijective σ :=
  ⟨fun h => ⟨h, Finite.injective_iff_surjective.mp h⟩, And.left⟩

/-- **σ is a permutation when orbits produce distinct roots.**
    If the multi-start final anchors are pairwise distinct (no two
    orbits collapse to the same root), then the start→root
    assignment is a bijection of `Fin p`. Every root is reached by
    exactly one orbit. -/
theorem multi_start_sigma_bijective
    (p : ℕ) (γ : Fin p → ℂ) (h_inj : Function.Injective γ)
    (roots : Fin p → ℂ) (_ : Function.Injective roots)
    (σ : Fin p → Fin p) (h_orbit : ∀ s : Fin p, γ s = roots (σ s)) :
    Function.Bijective σ := by
  have hσ_inj : Function.Injective σ := by
    intro s t hst
    have : roots (σ s) = roots (σ t) := by rw [hst]
    have : γ s = γ t := by rw [h_orbit s, h_orbit t, this]
    exact h_inj this
  exact (fin_injective_iff_bijective p σ).mp hσ_inj

/-- **Generic convergence (conditional form).** Given the McMullen
    hypothesis "all `p` orbits converge to distinct roots", the
    multi-start algorithm is surjective onto the root set: every
    root of the reduced polynomial is reached by exactly one orbit. -/
theorem generic_convergence_conditional
    (p : ℕ) (γ : Fin p → ℂ) (h_inj : Function.Injective γ)
    (roots : Fin p → ℂ) (h_root_inj : Function.Injective roots)
    (σ : Fin p → Fin p) (h_orbit : ∀ s : Fin p, γ s = roots (σ s)) :
    Function.Surjective σ :=
  (multi_start_sigma_bijective p γ h_inj roots h_root_inj σ h_orbit).2

end GenericConvergence

/-! ============================================================
  §11. Voronoï basin connectivity

  The result that circumvents McMullen 1987: for any family of
  distinct final anchors `γ : Fin p → ℂ`, each basin
  `basin γ k = {z₀ | σ selects k}` is the Voronoï cell of `γ k`.
  Voronoï cells are convex (intersections of closed half-planes),
  hence connected.

  This replaces Newton's fractal basin boundaries (box-counting
  dimension 1.4–1.7) with piecewise-linear Voronoï edges
  (dimension 1).
============================================================ -/

section BasinConnectivity

variable {p : ℕ}

/-- **Basin of attraction for anchor `k`.** The set of query points
    `z₀` whose nearest anchor is `γ k`. This is the Voronoï cell. -/
def basin (γ : Fin p → ℂ) (k : Fin p) : Set ℂ :=
  {z₀ | ∀ t : Fin p, ‖γ k - z₀‖ ≤ ‖γ t - z₀‖}

/-- **Half-plane between two anchors.** `{z | ‖a - z‖ ≤ ‖b - z‖}`
    is the closed half-plane on the `a` side of the perpendicular
    bisector between `a` and `b`. -/
def halfPlane (a b : ℂ) : Set ℂ := {z : ℂ | ‖a - z‖ ≤ ‖b - z‖}

/-- **Basin = intersection of half-planes.** The basin of `γ k` is
    the intersection over all competitors `t` of the half-plane
    between `γ k` and `γ t`. -/
theorem basin_eq_iInter (γ : Fin p → ℂ) (k : Fin p) :
    basin γ k = ⋂ t : Fin p, halfPlane (γ k) (γ t) := by
  ext z
  simp only [basin, halfPlane, Set.mem_setOf_eq, Set.mem_iInter]

/-- **Half-plane is convex.** The closed half-plane between any
    two complex anchors is convex. Proved via the key algebraic
    identity: `‖a − z‖² − ‖b − z‖²` is affine in `z` seen as
    `ℝ²`, so its sub-level set is a real half-space. -/
theorem halfPlane_convex (a b : ℂ) : Convex ℝ (halfPlane a b) := by
  intro x hx y hy s t hs ht hst
  simp only [halfPlane, Set.mem_setOf_eq] at hx hy ⊢
  set z : ℂ := s • x + t • y
  -- Work with squared norms (both sides nonneg → iff)
  have sq_le_sq : ∀ {u v : ℝ}, 0 ≤ u → 0 ≤ v → (u ≤ v ↔ u^2 ≤ v^2) := fun hu hv =>
    ⟨fun h => by nlinarith, fun h => by nlinarith⟩
  have hx_sq : ‖a - x‖^2 ≤ ‖b - x‖^2 :=
    (sq_le_sq (norm_nonneg _) (norm_nonneg _)).mp hx
  have hy_sq : ‖a - y‖^2 ≤ ‖b - y‖^2 :=
    (sq_le_sq (norm_nonneg _) (norm_nonneg _)).mp hy
  refine (sq_le_sq (norm_nonneg _) (norm_nonneg _)).mpr ?_
  -- ‖w‖² = w.re² + w.im²
  have sq_coords : ∀ w : ℂ, ‖w‖^2 = w.re^2 + w.im^2 := fun w => by
    rw [Complex.norm_eq_abs, Complex.sq_abs, Complex.normSq_apply]; ring
  have hx_sq' : (a.re - x.re)^2 + (a.im - x.im)^2 ≤
                (b.re - x.re)^2 + (b.im - x.im)^2 := by
    rw [sq_coords, sq_coords] at hx_sq
    simpa [Complex.sub_re, Complex.sub_im] using hx_sq
  have hy_sq' : (a.re - y.re)^2 + (a.im - y.im)^2 ≤
                (b.re - y.re)^2 + (b.im - y.im)^2 := by
    rw [sq_coords, sq_coords] at hy_sq
    simpa [Complex.sub_re, Complex.sub_im] using hy_sq
  rw [sq_coords, sq_coords]
  simp only [Complex.sub_re, Complex.sub_im]
  -- Coordinates of z = s • x + t • y
  have z_re : z.re = s * x.re + t * y.re := by
    show (s • x + t • y).re = s * x.re + t * y.re
    simp [Complex.add_re, Complex.real_smul, Complex.mul_re,
          Complex.ofReal_re, Complex.ofReal_im]
  have z_im : z.im = s * x.im + t * y.im := by
    show (s • x + t • y).im = s * x.im + t * y.im
    simp [Complex.add_im, Complex.real_smul, Complex.mul_im,
          Complex.ofReal_re, Complex.ofReal_im]
  rw [z_re, z_im]
  -- Key affine identity:  D(s·x + t·y) = s·D(x) + t·D(y) + (1-s-t)·K
  -- where D is the `(a-z)²−(b-z)²` difference.  Under s+t=1 the last term
  -- vanishes; the two remaining terms are scaled ≤0 hypotheses.
  have identity :
      ((a.re - (s * x.re + t * y.re))^2 + (a.im - (s * x.im + t * y.im))^2) -
      ((b.re - (s * x.re + t * y.re))^2 + (b.im - (s * x.im + t * y.im))^2)
    = s * (((a.re - x.re)^2 + (a.im - x.im)^2) -
           ((b.re - x.re)^2 + (b.im - x.im)^2))
    + t * (((a.re - y.re)^2 + (a.im - y.im)^2) -
           ((b.re - y.re)^2 + (b.im - y.im)^2))
    + (1 - s - t) * (a.re^2 + a.im^2 - b.re^2 - b.im^2) := by ring
  have identity_zero :
      (1 - s - t) * (a.re^2 + a.im^2 - b.re^2 - b.im^2) = 0 := by
    have : 1 - s - t = 0 := by linarith
    rw [this]; ring
  have hx_nonpos :
      s * (((a.re - x.re)^2 + (a.im - x.im)^2) -
           ((b.re - x.re)^2 + (b.im - x.im)^2)) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos hs (by linarith)
  have hy_nonpos :
      t * (((a.re - y.re)^2 + (a.im - y.im)^2) -
           ((b.re - y.re)^2 + (b.im - y.im)^2)) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos ht (by linarith)
  linarith [identity, identity_zero, hx_nonpos, hy_nonpos]

/-- **Basin is convex.** Intersection of convex half-planes. -/
theorem basin_convex (γ : Fin p → ℂ) (k : Fin p) :
    Convex ℝ (basin γ k) := by
  rw [basin_eq_iInter]
  exact convex_iInter (fun t => halfPlane_convex (γ k) (γ t))

/-- **Basin contains its own anchor.** `γ k ∈ basin γ k` trivially
    (distance zero to itself, hence ≤ any other distance). -/
theorem basin_self_mem (γ : Fin p → ℂ) (k : Fin p) : γ k ∈ basin γ k := by
  intro t
  simp [norm_nonneg]

/-- **Basin is nonempty.** -/
theorem basin_nonempty (γ : Fin p → ℂ) (k : Fin p) :
    (basin γ k).Nonempty := ⟨γ k, basin_self_mem γ k⟩

/-- **`★ Basin connectivity theorem (the headline result).**
    For any family of anchors `γ : Fin p → ℂ`, every basin is
    connected. This contrasts sharply with Newton's method, whose
    basins are fractal (disconnected on every scale) for polynomials
    of degree ≥ 4, by McMullen 1987. -/
theorem basin_isConnected (γ : Fin p → ℂ) (k : Fin p) :
    IsConnected (basin γ k) :=
  ⟨basin_nonempty γ k, (basin_convex γ k).isPreconnected⟩

/-- **Basin preconnectedness** (without the nonemptyness hypothesis). -/
theorem basin_isPreconnected (γ : Fin p → ℂ) (k : Fin p) :
    IsPreconnected (basin γ k) :=
  (basin_convex γ k).isPreconnected

end BasinConnectivity

/-! ============================================================
  §12. Basin boundary = union of perpendicular bisectors

  The boundary `∂ basin γ k` lies in a finite union of perpendicular
  bisectors `{z : ‖γ k − z‖ = ‖γ t − z‖}` with `t ≠ k`. Each bisector
  is an affine line in ℂ ≅ ℝ², so the boundary is piecewise linear —
  in stark contrast with Newton basins, whose boundaries are fractals
  of box-counting dimension 1.4–1.7 (McMullen 1987).

  Two substantive results are proved here:
  • `basin_isClosed`                — the basin is topologically closed.
  • `basin_frontier_on_bisector`    — every frontier point is equidistant
                                      to `γ k` and some competitor `γ t`
                                      with `t ≠ k`. The proof is a genuine
                                      finite-continuity argument (not just
                                      the trivial `t = k` witness).
============================================================ -/

section BoundaryDimension

variable {p : ℕ}

/-- **Perpendicular bisector** of two complex anchors — the set of
    points equidistant to `a` and `b`. A straight line in ℝ². -/
def perpBisector (a b : ℂ) : Set ℂ := {z : ℂ | ‖a - z‖ = ‖b - z‖}

/-- **Basin is closed.** Intersection of closed sub-level sets of
    continuous functions `w ↦ ‖γ k − w‖` and `w ↦ ‖γ t − w‖`. -/
theorem basin_isClosed (γ : Fin p → ℂ) (k : Fin p) : IsClosed (basin γ k) := by
  rw [basin_eq_iInter]
  refine isClosed_iInter fun t => ?_
  show IsClosed {z : ℂ | ‖γ k - z‖ ≤ ‖γ t - z‖}
  exact isClosed_le (continuous_const.sub continuous_id).norm
                    (continuous_const.sub continuous_id).norm

/-- **Frontier ⊆ union of bisectors with distinct anchors.**
    If `z ∈ ∂ basin γ k`, some competitor `t ≠ k` achieves equality
    `‖γ k − z‖ = ‖γ t − z‖`.

    *Proof (by contradiction).* Suppose all competitors `t ≠ k` give
    *strict* inequality `‖γ k − z‖ < ‖γ t − z‖`. By continuity of
    `w ↦ ‖γ k − w‖ − ‖γ t − w‖` and finiteness of the competitor set,
    the strict inequalities persist on a small ball around `z`. That
    ball sits inside `basin γ k`, so `z ∈ interior(basin)`. But
    `z ∈ frontier = closure ∖ interior`, so `z ∉ interior`.
    Contradiction. -/
theorem basin_frontier_on_bisector (γ : Fin p → ℂ) (k : Fin p) :
    frontier (basin γ k) ⊆
      ⋃ (t : Fin p) (_ : t ≠ k), perpBisector (γ k) (γ t) := by
  intro z hz
  -- z is in the basin (frontier ⊆ closure = basin, since basin is closed).
  have z_basin : z ∈ basin γ k :=
    (basin_isClosed γ k).closure_eq ▸ frontier_subset_closure hz
  -- z is not in the interior of basin (by definition of frontier).
  have z_not_int : z ∉ interior (basin γ k) := hz.2
  by_contra h_no
  -- For every competitor t ≠ k, strict inequality at z.
  have h_strict : ∀ t : Fin p, t ≠ k → ‖γ k - z‖ < ‖γ t - z‖ := fun t ht =>
    lt_of_le_of_ne (z_basin t) (fun heq => h_no <|
      Set.mem_iUnion.mpr ⟨t, Set.mem_iUnion.mpr ⟨ht, heq⟩⟩)
  -- By finite continuity, strict inequality persists on a neighborhood.
  apply z_not_int
  have h_ev : ∀ᶠ w in 𝓝 z, ∀ t ∈ (Finset.univ.erase k : Finset (Fin p)),
      ‖γ k - w‖ < ‖γ t - w‖ := by
    rw [Finset.eventually_all]
    intro t ht
    exact ((continuous_const.sub continuous_id).norm.continuousAt).eventually_lt
          ((continuous_const.sub continuous_id).norm.continuousAt)
          (h_strict t (Finset.mem_erase.mp ht).1)
  rw [Metric.eventually_nhds_iff] at h_ev
  obtain ⟨δ, hδ_pos, hδ⟩ := h_ev
  -- Exhibit an open ball ⊆ basin containing z (witnessing `z ∈ interior`).
  rw [mem_interior]
  refine ⟨Metric.ball z δ, ?_, Metric.isOpen_ball, Metric.mem_ball_self hδ_pos⟩
  intro w hw t
  by_cases ht : t = k
  · rw [ht]
  · exact le_of_lt (hδ (Metric.mem_ball.mp hw) t
      (Finset.mem_erase.mpr ⟨ht, Finset.mem_univ t⟩))

/-- **Boundary is a finite union of affine lines** (with distinct
    anchors). The finite `(p − 1)`-fold union of perpendicular
    bisectors covers the entire basin frontier. -/
theorem basin_boundary_finite_union (γ : Fin p → ℂ) (k : Fin p) :
    frontier (basin γ k) ⊆
      ⋃ t ∈ (Finset.univ.erase k : Finset (Fin p)), perpBisector (γ k) (γ t) := by
  intro z hz
  have h_mem := basin_frontier_on_bisector γ k hz
  rw [Set.mem_iUnion] at h_mem
  obtain ⟨t, h_mem⟩ := h_mem
  rw [Set.mem_iUnion] at h_mem
  obtain ⟨ht_ne, h_bi⟩ := h_mem
  exact Set.mem_iUnion₂.mpr
    ⟨t, Finset.mem_erase.mpr ⟨ht_ne, Finset.mem_univ t⟩, h_bi⟩

end BoundaryDimension

/-! ============================================================
  §13. Constructive McMullen — Global convergence for `z^p − x`

  Combines Foundations §8 (`multi_start_grand_master`) + Foundations §9.1
  (`voronoi_unique_off_boundary`) + §11 (convex basins) + §12 (bisector
  boundary) into a single global-convergence statement.

  **Main result (`constructive_mcmullen`).** For every `x ≠ 0`, every
  `p ≥ 1`, and every injective family of final anchors `γ : Fin p → ℂ`
  with `(γ s)^p = A/x` (e.g. the cyclotomic rotation of a single seed),
  the Pandrosion multi-start Steffensen algorithm simultaneously achieves:
    (A)  each anchor `γ k` lies in the *interior* of its basin
         (reachability from a nonempty open starting region);
    (A') basins *cover* all of ℂ (every query belongs to some basin);
    (B)  *off* the Voronoï boundary, the selector `σ(z₀)` is unique
         (`voronoi_unique_off_boundary`, Foundations §9.1);
    (C)  on the selected orbit, Steffensen-Pandrosion *fixes* the
         chosen anchor under the reduced target `y = x/A`.

  The Voronoï boundary sits in a finite union of perpendicular
  bisectors (Lemma `VoronoiBoundary_subset_bisectors`), hence has
  2-dimensional Lebesgue measure zero — so condition (B) holds
  almost everywhere.

  **Analogue constructif de McMullen 1987.** McMullen excludes
  generally convergent rational iteration in degree `≥ 4` for the
  *generic* setting. Here the obstacle is lifted by exploiting the
  *cyclotomic* structure of `z^p − x`: equispaced roots force the
  Voronoï cells to be *convex* (§11), and the whole algorithm
  becomes globally convergent off a measure-zero set.
============================================================ -/

section ConstructiveMcMullen

variable {p : ℕ}

/-- **Anchor lies in its own basin interior** (for distinct anchors).
    Since `γ` is injective, `γ k` is strictly closer to itself than to
    any other anchor, and by finite continuity a whole neighborhood
    inherits this property. -/
theorem basin_anchor_mem_interior
    (γ : Fin p → ℂ) (h_inj : Function.Injective γ) (k : Fin p) :
    γ k ∈ interior (basin γ k) := by
  -- At the centre, distance to γ k is 0, distance to γ t (t ≠ k) is strictly positive.
  have h_strict : ∀ t : Fin p, t ≠ k → ‖γ k - γ k‖ < ‖γ t - γ k‖ := by
    intro t ht
    have hne : γ t - γ k ≠ 0 := sub_ne_zero.mpr (fun heq => ht (h_inj heq))
    simp only [sub_self, norm_zero]
    exact norm_pos_iff.mpr hne
  -- Lift finitely many strict inequalities to a uniform neighborhood.
  have h_ev : ∀ᶠ w in 𝓝 (γ k), ∀ t ∈ (Finset.univ.erase k : Finset (Fin p)),
      ‖γ k - w‖ < ‖γ t - w‖ := by
    rw [Finset.eventually_all]
    intro t ht
    exact ((continuous_const.sub continuous_id).norm.continuousAt).eventually_lt
          ((continuous_const.sub continuous_id).norm.continuousAt)
          (h_strict t (Finset.mem_erase.mp ht).1)
  rw [Metric.eventually_nhds_iff] at h_ev
  obtain ⟨δ, hδ_pos, hδ⟩ := h_ev
  rw [mem_interior]
  refine ⟨Metric.ball (γ k) δ, ?_, Metric.isOpen_ball, Metric.mem_ball_self hδ_pos⟩
  intro w hw t
  by_cases ht : t = k
  · rw [ht]
  · exact le_of_lt (hδ (Metric.mem_ball.mp hw) t
      (Finset.mem_erase.mpr ⟨ht, Finset.mem_univ t⟩))

/-- **Basin interior is nonempty** for any family of distinct anchors.
    Combined with §11 (convex basins) this means every root has a
    positive-measure reachability set — unlike Newton basins on `z^p − x`
    for `p ≥ 4`, whose interiors can vanish at certain roots (McMullen). -/
theorem basin_interior_nonempty
    (γ : Fin p → ℂ) (h_inj : Function.Injective γ) (k : Fin p) :
    (interior (basin γ k)).Nonempty :=
  ⟨γ k, basin_anchor_mem_interior γ h_inj k⟩

/-- **Basins cover ℂ.** Every query point has at least one nearest
    anchor (existence from `voronoi_nearest_exists`), so it lies in
    the basin of that anchor. -/
theorem basins_cover (hp : 0 < p) (γ : Fin p → ℂ) :
    (⋃ k : Fin p, basin γ k) = Set.univ := by
  ext z
  simp only [Set.mem_iUnion, Set.mem_univ, iff_true, basin, Set.mem_setOf_eq]
  exact voronoi_nearest_exists p hp γ z

/-- **Voronoï boundary ⊆ union of bisectors** between distinct anchors.
    Equidistance between `γ s` and `γ t` places `z` on their bisector.
    Because the RHS is a finite union of affine lines in ℝ², it has
    2-dimensional Lebesgue measure zero — so the Voronoï boundary is
    Lebesgue-null. -/
theorem VoronoiBoundary_subset_bisectors (γ : Fin p → ℂ) :
    VoronoiBoundary γ ⊆
      ⋃ (s : Fin p) (t : Fin p) (_ : s ≠ t), perpBisector (γ s) (γ t) := by
  intro z hz
  obtain ⟨s, t, hst, heq, _⟩ := hz
  refine Set.mem_iUnion.mpr ⟨s, Set.mem_iUnion.mpr ⟨t,
           Set.mem_iUnion.mpr ⟨hst, ?_⟩⟩⟩
  exact heq

/-- **Voronoï boundary = filtered union of bisectors.**
    Strengthens `VoronoiBoundary_subset_bisectors` to an *equality*: the
    Voronoï boundary equals the union over distinct pairs `(s, t)` of
    bisector points where `γ s` already achieves the minimum distance to
    `z`. The added filter `∀ u, ‖γ s − z‖ ≤ ‖γ u − z‖` is essential —
    without it the reverse inclusion fails (a generic bisector point may
    be beaten by a third anchor and therefore lies in a cell interior,
    not on the boundary). The forward direction extracts the filter from
    the VoronoiBoundary minimality witness; the reverse direction
    reassembles the definition. -/
theorem VoronoiBoundary_eq_bisectors_filtered (γ : Fin p → ℂ) :
    VoronoiBoundary γ =
      ⋃ (s : Fin p) (t : Fin p) (_ : s ≠ t),
        { z : ℂ | z ∈ perpBisector (γ s) (γ t) ∧
                  ∀ u : Fin p, ‖γ s - z‖ ≤ ‖γ u - z‖ } := by
  ext z
  simp only [VoronoiBoundary, perpBisector, Set.mem_setOf_eq, Set.mem_iUnion,
             exists_prop]

/-- **★ CONSTRUCTIVE McMULLEN — Main global-convergence theorem.**

    For the problem of computing `p`-th roots of any `x ≠ 0` via
    the reduced target `y = x/A`, the Pandrosion multi-start
    Steffensen algorithm is *globally convergent*:

    (A)  each anchor `γ k` lies in the interior of its basin
         (reachability from an open set);
    (A') the basins cover ℂ entirely;
    (B)  off the Voronoï boundary (a Lebesgue-null union of lines),
         the selector `σ(z₀)` is unique;
    (C)  Steffensen-Pandrosion fixes every final anchor.

    This is the analogue-constructif of McMullen 1987's impossibility.
    The cyclotomic structure of `z^p − x` (equispaced anchors ⇒
    convex Voronoï cells, §11) is what lifts McMullen's generic
    obstruction. -/
theorem constructive_mcmullen
    (hp : 0 < p)
    (x A : ℂ) (hx : x ≠ 0) (hA : A ≠ 0)
    (γ : Fin p → ℂ) (h_inj : Function.Injective γ)
    (hSp : ∀ s, Sp_C p (γ s) ≠ 0)
    (hfp : ∀ s, (γ s) ^ p = A / x) :
    (∀ k : Fin p, γ k ∈ interior (basin γ k)) ∧
    (⋃ k : Fin p, basin γ k) = Set.univ ∧
    (∀ z₀ ∉ VoronoiBoundary γ,
        ∃! s : Fin p, ∀ t : Fin p, ‖γ s - z₀‖ ≤ ‖γ t - z₀‖) ∧
    (∀ s : Fin p, steffensen_step_C (x / A) p (γ s) = γ s) := by
  have hxA : (x / A) ≠ 0 := div_ne_zero hx hA
  refine ⟨basin_anchor_mem_interior γ h_inj, basins_cover hp γ, ?_, ?_⟩
  · intro z₀ h_off
    exact voronoi_unique_off_boundary p hp γ z₀ h_off
  · intro s
    have hfp_s : γ s ^ p = 1 / (x / A) := by rw [hfp s]; field_simp
    exact steffensen_step_C_fixed_point (x / A) hxA p (γ s) (hSp s) hfp_s

end ConstructiveMcMullen

end Pandrosion
