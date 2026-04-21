/-
  Universitas Pandrosion — Lean 4 Formalization
  ADVANCED² — Non-Trivial Extensions

  Four non-trivial results, each requiring genuine analytical,
  algebraic, or topological content beyond what Core and Advanced
  establish:

    (§I)   Global basin invariance (p=2) — removes the `h_inv`
           hypothesis from `error_after_n_steps_p2`, yielding
           *unconditional* global convergence of the Babylonian
           iteration from any `s₀ ≥ r`.

    (§II)  Pell-cone preservation — the Pandrosion-p2 iteration
           preserves the cone `{(A,B) : A² ≥ X·B²}` for X ≥ 0,
           via `norm_amplification_p2` and `Φ₂ ≥ 0`.

    (§III) Voronoï interior is open — the strict Voronoï cell
           `{z ∈ ℝ² : |z-r₁|² < |z-r₂|²}` is open in the
           product topology, with explicit δ-radius bound.

    (§IV)  Unconditional contraction on the basin — the Babylonian
           map is Lipschitz-½ on `[r, ∞)`, and the iteration has
           a unique fixed point at `r = √x` inside the basin.

  Each result is proved unconditionally with no sorry, no admit,
  and depends only on Core, Advanced, and Mathlib primitives.
-/
import Mathlib.Topology.MetricSpace.Basic
import Mathlib.Topology.Instances.Real
import Mathlib.Topology.Order.Basic
import Mathlib.Tactic
import Pandrosion.Core
import Pandrosion.Advanced

namespace Pandrosion

section Advanced2

/-! ## §I. Global Basin Invariance (p=2)

`error_after_n_steps_p2` in Core.lean assumes a hypothesis
  `h_inv : ∀ n, iterate 2 x s₀ n > 0 ∧ |iterate 2 x s₀ n - r| ≤ iterate 2 x s₀ n`
The original proof does not establish `h_inv`.

Here we prove `h_inv` holds automatically from the single
initial condition `s₀ ≥ r > 0`. The key is:

  **AM-GM**: for s > 0 and r² = x,  F(s) = (s²+x)/(2s) ≥ r.

Equivalently: `pandrosion_p2_identity` gives `F(s) − r = (s−r)²/(2s) ≥ 0`.
So every iterate satisfies `iterate n ≥ r`, preserving the basin
condition `|iterate n − r| ≤ iterate n` (since iterate ≥ r > 0).
-/

/-- **One-step basin preservation.**
    If `s > 0` and `r² = x` with `r > 0`, then after one Pandrosion-p2
    step, the new iterate is positive and remains in the basin
    `|F(s) − r| ≤ F(s)`. -/
theorem basin_preserved_p2 (x r s : ℝ) (hr : r > 0)
    (hs : s > 0) (h_root : r ^ 2 = x) :
    pandrosion_map 2 x s > 0 ∧
    |pandrosion_map 2 x s - r| ≤ pandrosion_map 2 x s := by
  have hid : pandrosion_map 2 x s - r = (s - r) ^ 2 / (2 * s) :=
    pandrosion_p2_identity x r s hs h_root
  have h2s_pos : (0 : ℝ) < 2 * s := by linarith
  have h_diff_nn : pandrosion_map 2 x s - r ≥ 0 := by
    rw [hid]
    exact div_nonneg (sq_nonneg _) (le_of_lt h2s_pos)
  have hF_ge_r : pandrosion_map 2 x s ≥ r := by linarith
  have hF_pos : pandrosion_map 2 x s > 0 := by linarith
  refine ⟨hF_pos, ?_⟩
  rw [abs_of_nonneg h_diff_nn]
  linarith

/-- **Full basin invariance.**
    Under `s₀ ≥ r > 0` and `r² = x`, every iterate satisfies the
    positivity + basin condition needed for `error_after_n_steps_p2`. -/
theorem basin_invariant_p2 (x r s₀ : ℝ) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) :
    ∀ n, iterate 2 x s₀ n > 0 ∧
         |iterate 2 x s₀ n - r| ≤ iterate 2 x s₀ n := by
  intro n
  induction n with
  | zero =>
    simp only [iterate]
    have hs₀_pos : s₀ > 0 := by linarith
    refine ⟨hs₀_pos, ?_⟩
    have h_diff_nn : s₀ - r ≥ 0 := by linarith
    rw [abs_of_nonneg h_diff_nn]
    linarith
  | succ n ih =>
    simp only [iterate]
    obtain ⟨hn_pos, _⟩ := ih
    exact basin_preserved_p2 x r (iterate 2 x s₀ n) hr hn_pos h_root

/-- **Unconditional global convergence (p=2).**
    For any `x > 0`, `r = √x > 0`, and any `s₀ ≥ r`, the Babylonian
    iteration converges at rate ≤ ½ per step:
        `|sₙ − r| ≤ (½)ⁿ · |s₀ − r|`. -/
theorem global_convergence_p2 (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ) :
    error (iterate 2 x s₀ n) r ≤ (1 / 2) ^ n * error s₀ r :=
  error_after_n_steps_p2 x r s₀ hx hr h_root
    (basin_invariant_p2 x r s₀ hr h_s₀ h_root) n

/-- **Limiting form.**
    The error converges to zero: for every `ε > 0`, some finite N
    suffices to drive `|sₙ − r|` below ε. -/
theorem global_convergence_eventually_p2
    (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (ε : ℝ) (hε : ε > 0) :
    ∃ N : ℕ, ∀ n ≥ N, error (iterate 2 x s₀ n) r < ε := by
  have h_err₀_nn : error s₀ r ≥ 0 := error_nonneg _ _
  by_cases h_err₀ : error s₀ r = 0
  · refine ⟨0, ?_⟩
    intro n _
    have hbound := global_convergence_p2 x r s₀ hx hr h_s₀ h_root n
    rw [h_err₀] at hbound
    simp at hbound
    linarith
  · have h_err₀_pos : error s₀ r > 0 := lt_of_le_of_ne h_err₀_nn (Ne.symm h_err₀)
    have h_half : (1 / 2 : ℝ) < 1 := by norm_num
    have h_half_nn : (0 : ℝ) ≤ 1 / 2 := by norm_num
    have h_ratio_pos : ε / error s₀ r > 0 := div_pos hε h_err₀_pos
    obtain ⟨N, hN⟩ := exists_pow_lt_of_lt_one h_ratio_pos h_half
    refine ⟨N, ?_⟩
    intro n hn
    have hbound := global_convergence_p2 x r s₀ hx hr h_s₀ h_root n
    have h_pow_le : (1 / 2 : ℝ) ^ n ≤ (1 / 2) ^ N :=
      pow_le_pow_of_le_one h_half_nn (le_of_lt h_half) hn
    calc error (iterate 2 x s₀ n) r
        ≤ (1 / 2) ^ n * error s₀ r := hbound
      _ ≤ (1 / 2) ^ N * error s₀ r :=
          mul_le_mul_of_nonneg_right h_pow_le h_err₀_nn
      _ < (ε / error s₀ r) * error s₀ r :=
          mul_lt_mul_of_pos_right hN h_err₀_pos
      _ = ε := by field_simp

/-! ## §II. Pell-Cone Preservation

For X ≥ 0, the set `C_X = {(A, B) ∈ ℤ² : A² ≥ X·B²}` is the
non-negative Pell cone. The Pandrosion-p2 iteration
    A' = A(A² + 2XB²),  B' = B(2A² + XB²)
preserves `C_X` because of two facts from Advanced.lean:

  * `norm_amplification_p2 : A'² − X·B'² = (A² − XB²) · Φ₂(A,B,X)`
  * `phi_p2_nonneg_squares : Φ₂(A,B,X) ≥ 0 whenever X ≥ 0`

Combined: the sign of `A² − XB²` propagates along the orbit.
If the initial state is in `C_X` (non-negative norm), the entire
orbit stays in `C_X`.

This is strictly weaker than the classical Pell group law
(composition of solutions), but follows *directly* from the
corpus' existing algebraic identities without new structural work.
-/

/-- **One-step Pell-cone preservation.** -/
theorem pell_cone_step (A B X : ℤ) (hX : 0 ≤ X) (h_cone : A ^ 2 ≥ X * B ^ 2) :
    (A * (A ^ 2 + 2 * X * B ^ 2)) ^ 2 ≥
    X * (B * (2 * A ^ 2 + X * B ^ 2)) ^ 2 := by
  have h_amp : (A * (A ^ 2 + 2 * X * B ^ 2)) ^ 2
               - X * (B * (2 * A ^ 2 + X * B ^ 2)) ^ 2
             = (A ^ 2 - X * B ^ 2) * pandrosion_phi_p2 A B X :=
    norm_amplification_p2 A B X
  have h_d_nn : A ^ 2 - X * B ^ 2 ≥ 0 := by linarith
  have h_phi_nn : 0 ≤ pandrosion_phi_p2 A B X := phi_p2_nonneg_squares A B X hX
  have h_nn : (A * (A ^ 2 + 2 * X * B ^ 2)) ^ 2
              - X * (B * (2 * A ^ 2 + X * B ^ 2)) ^ 2 ≥ 0 := by
    rw [h_amp]; exact mul_nonneg h_d_nn h_phi_nn
  linarith

/-- **The Pandrosion-p2 integer step.** -/
def pell_step (X : ℤ) : ℤ × ℤ → ℤ × ℤ :=
  fun AB => (AB.1 * (AB.1 ^ 2 + 2 * X * AB.2 ^ 2),
             AB.2 * (2 * AB.1 ^ 2 + X * AB.2 ^ 2))

/-- **n-fold iteration of the Pell step.** -/
def pell_iterate (X : ℤ) (AB : ℤ × ℤ) : ℕ → ℤ × ℤ
  | 0 => AB
  | n + 1 => pell_step X (pell_iterate X AB n)

/-- **Pell-cone invariance for the full orbit.**
    If `(A₀, B₀)` is in the non-negative cone `A₀² ≥ X·B₀²` for
    `X ≥ 0`, then every iterate `(Aₙ, Bₙ) = pell_iterate^n (A₀,B₀)`
    stays in the cone. -/
theorem pell_cone_invariant (X : ℤ) (hX : 0 ≤ X) (AB : ℤ × ℤ)
    (h_init : AB.1 ^ 2 ≥ X * AB.2 ^ 2) (n : ℕ) :
    (pell_iterate X AB n).1 ^ 2 ≥ X * (pell_iterate X AB n).2 ^ 2 := by
  induction n with
  | zero => exact h_init
  | succ n ih =>
    show (pell_step X (pell_iterate X AB n)).1 ^ 2 ≥
         X * (pell_step X (pell_iterate X AB n)).2 ^ 2
    unfold pell_step
    exact pell_cone_step (pell_iterate X AB n).1 (pell_iterate X AB n).2 X hX ih

/-- **Sign-propagation corollary: trivial Pell solutions stay trivial.**
    If the initial norm is zero (a trivial Pell solution), the orbit
    preserves that zero — `A² = X·B²` is preserved. -/
theorem pell_norm_zero_preserved (X : ℤ) (AB : ℤ × ℤ)
    (h_init : AB.1 ^ 2 = X * AB.2 ^ 2) (n : ℕ) :
    (pell_iterate X AB n).1 ^ 2 = X * (pell_iterate X AB n).2 ^ 2 := by
  induction n with
  | zero => exact h_init
  | succ n ih =>
    set A := (pell_iterate X AB n).1
    set B := (pell_iterate X AB n).2
    show (A * (A ^ 2 + 2 * X * B ^ 2)) ^ 2 =
         X * (B * (2 * A ^ 2 + X * B ^ 2)) ^ 2
    have h_amp : (A * (A ^ 2 + 2 * X * B ^ 2)) ^ 2
                 - X * (B * (2 * A ^ 2 + X * B ^ 2)) ^ 2
               = (A ^ 2 - X * B ^ 2) * pandrosion_phi_p2 A B X :=
      norm_amplification_p2 A B X
    have h_d_zero : A ^ 2 - X * B ^ 2 = 0 := by linarith
    rw [h_d_zero, zero_mul] at h_amp
    linarith

/-! ## §III. The Strict Voronoï Interior is Open

The Pandrosion multi-start relies on assigning each starting
point to its *nearest* anchor (Voronoï assignment). For
robustness, we need the interior of each Voronoï cell to be
*open* in the ambient topology — a small perturbation of z
should keep it in the same cell.

**Key observation:** the squared-distance-difference function
    D(z) = |z − r₁|² − |z − r₂|²
is a polynomial in `z.1, z.2, r₁.1, r₁.2, r₂.1, r₂.2`, hence
continuous in z. The strict Voronoï cell is `D⁻¹(−∞, 0)`, the
preimage of an open interval under a continuous map — therefore
open.

We also provide an explicit δ-radius witnessing openness.
-/

/-- **The squared-distance difference is continuous in z.** -/
theorem squared_distance_diff_continuous (r1x r1y r2x r2y : ℝ) :
    Continuous (fun p : ℝ × ℝ =>
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 -
      ((p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2)) := by
  have c1 : Continuous (fun p : ℝ × ℝ => p.1) := continuous_fst
  have c2 : Continuous (fun p : ℝ × ℝ => p.2) := continuous_snd
  exact ((c1.sub continuous_const).pow 2).add ((c2.sub continuous_const).pow 2)
        |>.sub (((c1.sub continuous_const).pow 2).add ((c2.sub continuous_const).pow 2))

/-- **The strict Voronoï half-plane is open in ℝ².** -/
theorem voronoi_strict_open (r1x r1y r2x r2y : ℝ) :
    IsOpen {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2} := by
  let D : ℝ × ℝ → ℝ := fun p =>
    (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 -
    ((p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2)
  have hcont : Continuous D := squared_distance_diff_continuous r1x r1y r2x r2y
  have h_eq : {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2} = D ⁻¹' Set.Iio 0 := by
    ext p
    simp only [D, Set.mem_setOf_eq, Set.mem_preimage, Set.mem_Iio]
    constructor
    · intro h; linarith
    · intro h; linarith
  rw [h_eq]
  exact isOpen_Iio.preimage hcont

/-- **Explicit δ-neighborhood witness.**
    If `(z₁, z₂)` is strictly closer to `r₁` than to `r₂` in ℝ²,
    there is a radius δ > 0 such that every point within sup-distance
    δ of `(z₁, z₂)` also lies strictly closer to `r₁`. -/
theorem voronoi_strict_delta_exists (r1x r1y r2x r2y z1 z2 : ℝ)
    (h : (z1 - r1x) ^ 2 + (z2 - r1y) ^ 2 <
         (z1 - r2x) ^ 2 + (z2 - r2y) ^ 2) :
    ∃ δ > 0, ∀ w1 w2 : ℝ, |w1 - z1| < δ → |w2 - z2| < δ →
      (w1 - r1x) ^ 2 + (w2 - r1y) ^ 2 <
      (w1 - r2x) ^ 2 + (w2 - r2y) ^ 2 := by
  have hopen := voronoi_strict_open r1x r1y r2x r2y
  have hmem : (z1, z2) ∈ {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2} := by
    show (z1 - r1x) ^ 2 + (z2 - r1y) ^ 2 < (z1 - r2x) ^ 2 + (z2 - r2y) ^ 2
    exact h
  obtain ⟨δ, hδ_pos, hδ⟩ := Metric.isOpen_iff.mp hopen (z1, z2) hmem
  refine ⟨δ, hδ_pos, ?_⟩
  intro w1 w2 hw1 hw2
  have hball : (w1, w2) ∈ Metric.ball (z1, z2) δ := by
    rw [Metric.mem_ball, Prod.dist_eq, max_lt_iff]
    refine ⟨?_, ?_⟩
    · simpa [Real.dist_eq] using hw1
    · simpa [Real.dist_eq] using hw2
  exact hδ hball

/-! ## §IV. Contraction-Mapping Structure (p=2)

Mathlib's `ContractingWith K f` formally requires `f : α → α`
with `LipschitzWith K f` and `K < 1`. The Babylonian map is not
globally a contraction on `ℝ` (it is only defined for `s > 0`
and only contracting toward `r` in the half-line `s ≥ r`).

The content of door #13 is the **unconditional** Lipschitz-½
bound inside the basin, packaged as a one-shot lemma that no
longer requires the basin hypothesis `h_inv` or the positivity
hypothesis `hs` separately — `s ≥ r > 0` suffices.
-/

/-- **Unconditional Lipschitz-½ contraction inside the basin.**
    For any `s ≥ r` with `r > 0` and `r² = x`,
        `|F(s) − r| ≤ (1/2) · |s − r|`.
    No separate basin / positivity / continuity hypotheses needed. -/
theorem pandrosion_map_lipschitz_half_p2
    (x r s : ℝ) (hx : x > 0) (hr : r > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r) :
    error (pandrosion_map 2 x s) r ≤ (1 / 2) * error s r := by
  have hs_pos : s > 0 := by linarith
  have h_basin_abs : |s - r| ≤ s := by
    rw [abs_of_nonneg (by linarith : s - r ≥ 0)]; linarith
  exact contraction_step_p2 x r s hx hr hs_pos h_root h_basin_abs

/-- **Iterated Lipschitz-½: the basin is a contracting set.**
    The Babylonian iteration, restricted to `s ≥ r`, is a strict
    contraction with rate ≤ 1/2 per step, for every orbit. -/
theorem pandrosion_iterate_lipschitz_half_p2
    (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ) :
    error (iterate 2 x s₀ n) r ≤ (1 / 2) ^ n * error s₀ r :=
  global_convergence_p2 x r s₀ hx hr h_s₀ h_root n

/-- **Uniqueness of the fixed point in the basin.**
    On the closed half-line `[r, ∞)`, the Babylonian map has a
    unique fixed point at `r = √x`. -/
theorem pandrosion_map_fixed_point_unique_p2
    (x r s : ℝ) (hx : x > 0) (hr : r > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r)
    (h_fix : pandrosion_map 2 x s = s) : s = r := by
  have h_lip := pandrosion_map_lipschitz_half_p2 x r s hx hr h_root h_basin
  rw [h_fix] at h_lip
  -- Now: error s r ≤ (1/2) * error s r, hence error s r = 0, hence s = r.
  have h_err_nn : error s r ≥ 0 := error_nonneg _ _
  have h_half : error s r ≤ 0 := by linarith
  have h_err_zero : error s r = 0 := le_antisymm h_half h_err_nn
  unfold error at h_err_zero
  have : s - r = 0 := abs_eq_zero.mp h_err_zero
  linarith

/-! ## §V. Unified Advanced² Certificate

A single statement packaging §I–§IV: the Pandrosion corpus
admits (i) unconditional global convergence for the Babylonian
iteration, (ii) Pell-cone preservation under the integer-norm
iteration, (iii) topological openness of Voronoï interiors, and
(iv) uniqueness of the Babylonian fixed point in the basin.
-/

/-- **Unified Advanced² certificate.**
    Four unconditional facts, each non-trivial, assembled from
    §I, §II, §III, §IV. -/
theorem advanced2_grand_certificate
    (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ)
    (A B X : ℤ) (hX : 0 ≤ X) (h_cone : A ^ 2 ≥ X * B ^ 2)
    (r1x r1y r2x r2y : ℝ) :
    -- (I) Unconditional convergence:
    (error (iterate 2 x s₀ n) r ≤ (1 / 2) ^ n * error s₀ r) ∧
    -- (II) Pell-cone preservation:
    ((pell_iterate X (A, B) n).1 ^ 2 ≥ X * (pell_iterate X (A, B) n).2 ^ 2) ∧
    -- (III) Voronoï openness:
    (IsOpen {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2}) ∧
    -- (IV) Lipschitz-½ contraction on the basin:
    (error (pandrosion_map 2 x s₀) r ≤ (1 / 2) * error s₀ r) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · exact global_convergence_p2 x r s₀ hx hr h_s₀ h_root n
  · exact pell_cone_invariant X hX (A, B) h_cone n
  · exact voronoi_strict_open r1x r1y r2x r2y
  · exact pandrosion_map_lipschitz_half_p2 x r s₀ hx hr h_root h_s₀

/-! ## §VI. Monotone Descent of the Orbit (#A)

For `s ≥ r > 0` with `r² = x`, each Babylonian step is a **descent**:
`F(s) ≤ s`. Combined with §I (basin invariance), this gives full
monotone descent of the orbit, so `(iterate 2 x s₀ n)` is a
non-increasing sequence bounded below by `r`.

**Proof**: `F(s) - r = (s-r)²/(2s)` and `F(s) - s = -(s-r)(s+r)/(2s)
≤ 0` for `s ≥ r > 0`.
-/

/-- **One-step descent.** -/
theorem pandrosion_map_le_self_p2 (x r s : ℝ) (hr : r > 0) (hs : s > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r) :
    pandrosion_map 2 x s ≤ s := by
  have hid : pandrosion_map 2 x s - r = (s - r) ^ 2 / (2 * s) :=
    pandrosion_p2_identity x r s hs h_root
  have h2s_pos : (0 : ℝ) < 2 * s := by linarith
  have h_num_bnd : (s - r) ^ 2 ≤ (s - r) * (2 * s) := by nlinarith
  have h_div_bnd : (s - r) ^ 2 / (2 * s) ≤ s - r := by
    rw [div_le_iff h2s_pos]; linarith
  linarith

/-- **Every iterate is ≥ r** (stronger than basin invariance). -/
theorem iterate_ge_r_p2 (x r s₀ : ℝ) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ) :
    iterate 2 x s₀ n ≥ r := by
  induction n with
  | zero => simp only [iterate]; exact h_s₀
  | succ n ih =>
    simp only [iterate]
    have hn_pos : iterate 2 x s₀ n > 0 := by linarith
    have hid : pandrosion_map 2 x (iterate 2 x s₀ n) - r
             = (iterate 2 x s₀ n - r) ^ 2 / (2 * iterate 2 x s₀ n) :=
      pandrosion_p2_identity x r (iterate 2 x s₀ n) hn_pos h_root
    have h2_pos : 2 * iterate 2 x s₀ n > 0 := by linarith
    have h_nn : (iterate 2 x s₀ n - r) ^ 2 / (2 * iterate 2 x s₀ n) ≥ 0 :=
      div_nonneg (sq_nonneg _) (le_of_lt h2_pos)
    linarith

/-- **Full orbit monotone descent.** -/
theorem iterate_monotone_descent_p2 (x r s₀ : ℝ) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ) :
    iterate 2 x s₀ (n + 1) ≤ iterate 2 x s₀ n := by
  have h_ge_r := iterate_ge_r_p2 x r s₀ hr h_s₀ h_root n
  have hn_pos : iterate 2 x s₀ n > 0 := by linarith
  show pandrosion_map 2 x (iterate 2 x s₀ n) ≤ iterate 2 x s₀ n
  exact pandrosion_map_le_self_p2 x r (iterate 2 x s₀ n) hr hn_pos h_root h_ge_r

/-! ## §VII. Quadratic Convergence (#B)

The contraction of §IV was `|F(s) - r| ≤ (1/2)|s - r|` (linear).
In fact, the Babylonian method converges **quadratically**:
    `|F(s) - r| = (s-r)² / (2s) ≤ (s-r)² / (2r)`
for `s ≥ r > 0`. The error is squared at each step (up to a
factor `1/(2r)`).
-/

/-- **Quadratic error bound.** -/
theorem pandrosion_quadratic_convergence_p2 (x r s : ℝ) (hr : r > 0) (hs : s > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r) :
    |pandrosion_map 2 x s - r| ≤ (s - r) ^ 2 / (2 * r) := by
  have hid : pandrosion_map 2 x s - r = (s - r) ^ 2 / (2 * s) :=
    pandrosion_p2_identity x r s hs h_root
  rw [hid]
  have h2s_pos : (0 : ℝ) < 2 * s := by linarith
  have h2r_pos : (0 : ℝ) < 2 * r := by linarith
  have h_sqr_nn : 0 ≤ (s - r) ^ 2 := sq_nonneg _
  have h_diff_nn : (s - r) ^ 2 / (2 * s) ≥ 0 :=
    div_nonneg h_sqr_nn (le_of_lt h2s_pos)
  rw [abs_of_nonneg h_diff_nn]
  rw [div_le_div_iff h2s_pos h2r_pos]
  nlinarith [sq_nonneg (s - r)]

/-- **Iterated quadratic bound (single-step form).**
    After one step, the error is squared — not merely halved. -/
theorem pandrosion_error_squared_p2 (x r s : ℝ) (hr : r > 0) (hs : s > 0)
    (h_root : r ^ 2 = x) (h_basin : s ≥ r) :
    error (pandrosion_map 2 x s) r ≤ (error s r) ^ 2 / (2 * r) := by
  unfold error
  have h_sr_nn : s - r ≥ 0 := by linarith
  have : |s - r| = s - r := abs_of_nonneg h_sr_nn
  rw [this]
  exact pandrosion_quadratic_convergence_p2 x r s hr hs h_root h_basin

/-! ## §VIII. Strict Pell-Cone Preservation (#C)

§II showed the weak cone `A² ≥ X·B²` is preserved. A stronger
fact: the **strict** cone `A² > X·B²` is also preserved (i.e.
no orbit can degenerate from strict to weak).

**Proof**: `A'² - X·B'² = (A² - X·B²) · Φ₂(A,B,X)`. If the left
factor is > 0, and Φ₂ > 0 (which follows from X ≥ 0 and `A ≠ 0`,
which follows from A² > XB² ≥ 0), then the product is > 0.
-/

/-- **One-step strict Pell-cone preservation.** -/
theorem pell_cone_step_strict (A B X : ℤ) (hX : 0 ≤ X)
    (h_cone : A ^ 2 > X * B ^ 2) :
    (A * (A ^ 2 + 2 * X * B ^ 2)) ^ 2 >
    X * (B * (2 * A ^ 2 + X * B ^ 2)) ^ 2 := by
  have h_amp : (A * (A ^ 2 + 2 * X * B ^ 2)) ^ 2
               - X * (B * (2 * A ^ 2 + X * B ^ 2)) ^ 2
             = (A ^ 2 - X * B ^ 2) * pandrosion_phi_p2 A B X :=
    norm_amplification_p2 A B X
  have hXB2_nn : 0 ≤ X * B ^ 2 := mul_nonneg hX (sq_nonneg _)
  have hA2_pos : 0 < A ^ 2 := by linarith
  have hA4_pos : 0 < A ^ 4 := by
    have : A ^ 4 = (A ^ 2) ^ 2 := by ring
    rw [this]; exact pow_pos hA2_pos 2
  have h_phi_pos : 0 < pandrosion_phi_p2 A B X := by
    unfold pandrosion_phi_p2
    have h2 : 0 ≤ X * A ^ 2 * B ^ 2 :=
      mul_nonneg (mul_nonneg hX (sq_nonneg _)) (sq_nonneg _)
    have h3 : 0 ≤ X ^ 2 * B ^ 4 := by positivity
    linarith
  have h_diff_pos : 0 < A ^ 2 - X * B ^ 2 := by linarith
  have h_prod_pos : 0 < (A ^ 2 - X * B ^ 2) * pandrosion_phi_p2 A B X :=
    mul_pos h_diff_pos h_phi_pos
  linarith

/-- **Strict Pell-cone invariance for the full orbit.** -/
theorem pell_cone_invariant_strict (X : ℤ) (hX : 0 ≤ X) (AB : ℤ × ℤ)
    (h_init : AB.1 ^ 2 > X * AB.2 ^ 2) (n : ℕ) :
    (pell_iterate X AB n).1 ^ 2 > X * (pell_iterate X AB n).2 ^ 2 := by
  induction n with
  | zero => exact h_init
  | succ n ih =>
    show (pell_step X (pell_iterate X AB n)).1 ^ 2 >
         X * (pell_step X (pell_iterate X AB n)).2 ^ 2
    unfold pell_step
    exact pell_cone_step_strict (pell_iterate X AB n).1
      (pell_iterate X AB n).2 X hX ih

/-! ## §IX. Finite Voronoï Intersection (#D)

For a finite set of anchor points `anchors`, the strict Voronoï
cell of a point `r₀` is the intersection of the open half-planes
`{z : |z - r₀|² < |z - r|²}` over all `r ∈ anchors`. A finite
intersection of open sets is open, so the full cell is open.
-/

/-- **Finite Voronoï cell is open in ℝ².** -/
theorem voronoi_cell_finite_open (r0 : ℝ × ℝ) (anchors : Finset (ℝ × ℝ)) :
    IsOpen {p : ℝ × ℝ | ∀ a ∈ anchors,
      (p.1 - r0.1) ^ 2 + (p.2 - r0.2) ^ 2 <
      (p.1 - a.1) ^ 2 + (p.2 - a.2) ^ 2} := by
  have h_eq : {p : ℝ × ℝ | ∀ a ∈ anchors,
      (p.1 - r0.1) ^ 2 + (p.2 - r0.2) ^ 2 <
      (p.1 - a.1) ^ 2 + (p.2 - a.2) ^ 2} =
      ⋂ a ∈ (anchors : Set (ℝ × ℝ)), {p : ℝ × ℝ |
        (p.1 - r0.1) ^ 2 + (p.2 - r0.2) ^ 2 <
        (p.1 - a.1) ^ 2 + (p.2 - a.2) ^ 2} := by
    ext p
    simp only [Set.mem_setOf_eq, Set.mem_iInter, Finset.mem_coe]
  rw [h_eq]
  refine Set.Finite.isOpen_biInter (Finset.finite_toSet _) ?_
  intro a _
  exact voronoi_strict_open r0.1 r0.2 a.1 a.2

/-! ## §X. Pell Orbit Strict Growth (#E)

For `X ≥ 1` and `A₀, B₀ ≥ 1`, every iterate stays in the cone
`A, B ≥ 1` AND strictly grows: `A_{n+1} > A_n` and `B_{n+1} > B_n`.

**Consequence**: the orbit never cycles — each iterate is strictly
larger than the last, hence all pairs `(A_n, B_n)` are distinct.
-/

/-- **One-step positivity preservation.** -/
theorem pell_step_pos_preserved (X : ℤ) (AB : ℤ × ℤ) (hX : 1 ≤ X)
    (hA : 1 ≤ AB.1) (hB : 1 ≤ AB.2) :
    1 ≤ (pell_step X AB).1 ∧ 1 ≤ (pell_step X AB).2 := by
  unfold pell_step
  have hA2 : 1 ≤ AB.1 ^ 2 := by nlinarith [sq_nonneg (AB.1 - 1)]
  have hB2 : 1 ≤ AB.2 ^ 2 := by nlinarith [sq_nonneg (AB.2 - 1)]
  have hXB2 : 1 ≤ X * AB.2 ^ 2 := by
    have h1 : (1 : ℤ) * 1 ≤ X * AB.2 ^ 2 :=
      mul_le_mul hX hB2 zero_le_one (by linarith)
    linarith
  have h_fac1 : 1 ≤ AB.1 ^ 2 + 2 * X * AB.2 ^ 2 := by linarith
  have h_fac2 : 1 ≤ 2 * AB.1 ^ 2 + X * AB.2 ^ 2 := by linarith
  have h_fac1_nn : 0 ≤ AB.1 ^ 2 + 2 * X * AB.2 ^ 2 := by linarith
  have h_fac2_nn : 0 ≤ 2 * AB.1 ^ 2 + X * AB.2 ^ 2 := by linarith
  refine ⟨?_, ?_⟩
  · show 1 ≤ AB.1 * (AB.1 ^ 2 + 2 * X * AB.2 ^ 2)
    have step : AB.1 ^ 2 + 2 * X * AB.2 ^ 2
                ≤ AB.1 * (AB.1 ^ 2 + 2 * X * AB.2 ^ 2) :=
      le_mul_of_one_le_left h_fac1_nn hA
    linarith
  · show 1 ≤ AB.2 * (2 * AB.1 ^ 2 + X * AB.2 ^ 2)
    have step : 2 * AB.1 ^ 2 + X * AB.2 ^ 2
                ≤ AB.2 * (2 * AB.1 ^ 2 + X * AB.2 ^ 2) :=
      le_mul_of_one_le_left h_fac2_nn hB
    linarith

/-- **One-step strict growth.** -/
theorem pell_step_strict_growth (X : ℤ) (AB : ℤ × ℤ) (hX : 1 ≤ X)
    (hA : 1 ≤ AB.1) (hB : 1 ≤ AB.2) :
    AB.1 < (pell_step X AB).1 ∧ AB.2 < (pell_step X AB).2 := by
  unfold pell_step
  have hA2 : 1 ≤ AB.1 ^ 2 := by nlinarith [sq_nonneg (AB.1 - 1)]
  have hB2 : 1 ≤ AB.2 ^ 2 := by nlinarith [sq_nonneg (AB.2 - 1)]
  have hXB2 : 1 ≤ X * AB.2 ^ 2 := by
    have h1 : (1 : ℤ) * 1 ≤ X * AB.2 ^ 2 :=
      mul_le_mul hX hB2 zero_le_one (by linarith)
    linarith
  have h_fac1 : 2 ≤ AB.1 ^ 2 + 2 * X * AB.2 ^ 2 := by linarith
  have h_fac2 : 2 ≤ 2 * AB.1 ^ 2 + X * AB.2 ^ 2 := by linarith
  refine ⟨?_, ?_⟩
  · show AB.1 < AB.1 * (AB.1 ^ 2 + 2 * X * AB.2 ^ 2)
    have hA_pos : 0 < AB.1 := by linarith
    have : AB.1 * 1 < AB.1 * (AB.1 ^ 2 + 2 * X * AB.2 ^ 2) :=
      mul_lt_mul_of_pos_left (by linarith) hA_pos
    linarith
  · show AB.2 < AB.2 * (2 * AB.1 ^ 2 + X * AB.2 ^ 2)
    have hB_pos : 0 < AB.2 := by linarith
    have : AB.2 * 1 < AB.2 * (2 * AB.1 ^ 2 + X * AB.2 ^ 2) :=
      mul_lt_mul_of_pos_left (by linarith) hB_pos
    linarith

/-- **Full orbit positivity invariance.** -/
theorem pell_iterate_pos_preserved (X : ℤ) (hX : 1 ≤ X) (AB : ℤ × ℤ)
    (hA : 1 ≤ AB.1) (hB : 1 ≤ AB.2) (n : ℕ) :
    1 ≤ (pell_iterate X AB n).1 ∧ 1 ≤ (pell_iterate X AB n).2 := by
  induction n with
  | zero => exact ⟨hA, hB⟩
  | succ n ih =>
    exact pell_step_pos_preserved X (pell_iterate X AB n) hX ih.1 ih.2

/-- **Full orbit strict growth.** -/
theorem pell_iterate_strict_growth (X : ℤ) (hX : 1 ≤ X) (AB : ℤ × ℤ)
    (hA : 1 ≤ AB.1) (hB : 1 ≤ AB.2) (n : ℕ) :
    (pell_iterate X AB n).1 < (pell_iterate X AB (n + 1)).1 ∧
    (pell_iterate X AB n).2 < (pell_iterate X AB (n + 1)).2 := by
  have ⟨hAn, hBn⟩ := pell_iterate_pos_preserved X hX AB hA hB n
  show (pell_iterate X AB n).1 <
         (pell_step X (pell_iterate X AB n)).1 ∧
       (pell_iterate X AB n).2 <
         (pell_step X (pell_iterate X AB n)).2
  exact pell_step_strict_growth X (pell_iterate X AB n) hX hAn hBn

/-! ## §XI. Pairwise Banach Contraction (#F)

Door #13 only proved contraction **to the fixed point**:
    `|F(s) - r| ≤ (1/2)|s - r|`.
The full **Banach contraction** condition is pairwise:
    `|F(s) - F(t)| ≤ (1/2)|s - t|` for all `s, t ≥ r`.

**Key identity**: for `s, t > 0`,
    `F(s) - F(t) = (s - t)(st - x) / (2st)`.
For `s, t ≥ r` with `r² = x`, we have `st ≥ r² = x ≥ 0`, so
`(st - x)/(2st) ≤ 1/2`, giving the Lipschitz-½ bound.

This is exactly the hypothesis Banach's fixed-point theorem
requires; combined with basin invariance (§I), the Babylonian
map `F : [r, ∞) → [r, ∞)` is a full contraction mapping with
unique fixed point `r`.
-/

/-- **The difference identity for `pandrosion_map 2`.** -/
theorem pandrosion_map_diff_p2 (x s t : ℝ) (hs : s > 0) (ht : t > 0) :
    pandrosion_map 2 x s - pandrosion_map 2 x t =
    (s - t) * (s * t - x) / (2 * (s * t)) := by
  unfold pandrosion_map
  simp only
  have _hs_ne : s ≠ 0 := ne_of_gt hs
  have _ht_ne : t ≠ 0 := ne_of_gt ht
  have hs2_ne : (2 : ℝ) * s ^ 2 ≠ 0 := by positivity
  have ht2_ne : (2 : ℝ) * t ^ 2 ≠ 0 := by positivity
  have hden_s_ne :
      (↑(2 : ℕ) : ℝ) * s ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x ≠ 0 := by
    rw [show (↑(2 : ℕ) : ℝ) * s ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x =
          2 * s ^ 2 from by push_cast; ring]
    exact hs2_ne
  have hden_t_ne :
      (↑(2 : ℕ) : ℝ) * t ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x ≠ 0 := by
    rw [show (↑(2 : ℕ) : ℝ) * t ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x =
          2 * t ^ 2 from by push_cast; ring]
    exact ht2_ne
  rw [if_neg hden_s_ne, if_neg hden_t_ne]
  have h_sden :
      (↑(2 : ℕ) : ℝ) * s ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x = 2 * s ^ 2 := by
    push_cast; ring
  have h_tden :
      (↑(2 : ℕ) : ℝ) * t ^ 2 + ((↑(2 : ℕ) : ℝ) - 1) * x - x = 2 * t ^ 2 := by
    push_cast; ring
  rw [h_sden, h_tden]
  field_simp
  push_cast
  ring

/-- **Pairwise Lipschitz-½ on the basin — full Banach contraction.** -/
theorem pandrosion_map_lipschitz_pairwise_p2
    (x r s t : ℝ) (hr : r > 0) (h_root : r ^ 2 = x)
    (hs : s ≥ r) (ht : t ≥ r) :
    |pandrosion_map 2 x s - pandrosion_map 2 x t| ≤ (1 / 2) * |s - t| := by
  have hs_pos : s > 0 := by linarith
  have ht_pos : t > 0 := by linarith
  have hx_nn : 0 ≤ x := by rw [← h_root]; exact sq_nonneg _
  have h_diff :
      pandrosion_map 2 x s - pandrosion_map 2 x t
        = (s - t) * (s * t - x) / (2 * (s * t)) :=
    pandrosion_map_diff_p2 x s t hs_pos ht_pos
  rw [h_diff]
  have hst_pos : 0 < s * t := mul_pos hs_pos ht_pos
  have h2st_pos : (0 : ℝ) < 2 * (s * t) := by linarith
  have hst_ge_x : x ≤ s * t := by
    rw [← h_root]
    have : r * r ≤ s * t := mul_le_mul hs ht (le_of_lt hr) (by linarith)
    nlinarith
  have h_st_x_nn : 0 ≤ s * t - x := by linarith
  rw [abs_div, abs_mul,
      abs_of_pos h2st_pos, abs_of_nonneg h_st_x_nn]
  have h_abs_nn : 0 ≤ |s - t| := abs_nonneg _
  have h_ratio : (s * t - x) / (2 * (s * t)) ≤ 1 / 2 := by
    rw [div_le_div_iff h2st_pos (by norm_num : (0:ℝ) < 2)]
    linarith
  calc |s - t| * (s * t - x) / (2 * (s * t))
      = |s - t| * ((s * t - x) / (2 * (s * t))) := by ring
    _ ≤ |s - t| * (1 / 2) :=
        mul_le_mul_of_nonneg_left h_ratio h_abs_nn
    _ = (1 / 2) * |s - t| := by ring

/-! ## §XII. Super-Certificate — Ten Unconditional Facts

All results of §I–§XI, bundled in a single statement. -/

/-- **Super-certificate.** Packages every non-trivial fact in
    Advanced2 into one compound theorem. -/
theorem advanced2_super_certificate
    (x r s₀ : ℝ) (hx : x > 0) (hr : r > 0)
    (h_s₀ : s₀ ≥ r) (h_root : r ^ 2 = x) (n : ℕ)
    (A B X : ℤ) (hX : 1 ≤ X) (hA : 1 ≤ A) (hB : 1 ≤ B)
    (h_cone_strict : A ^ 2 > X * B ^ 2)
    (r1x r1y r2x r2y : ℝ) (anchors : Finset (ℝ × ℝ))
    (r0 : ℝ × ℝ) (s t : ℝ) (hs : s ≥ r) (ht : t ≥ r) :
    -- (I) Unconditional convergence:
    (error (iterate 2 x s₀ n) r ≤ (1 / 2) ^ n * error s₀ r) ∧
    -- (II) Pell-cone preservation (weak):
    ((pell_iterate X (A, B) n).1 ^ 2 ≥ X * (pell_iterate X (A, B) n).2 ^ 2) ∧
    -- (III) Voronoï openness (pairwise):
    (IsOpen {p : ℝ × ℝ |
      (p.1 - r1x) ^ 2 + (p.2 - r1y) ^ 2 <
      (p.1 - r2x) ^ 2 + (p.2 - r2y) ^ 2}) ∧
    -- (IV) Lipschitz-½ to fixed point:
    (error (pandrosion_map 2 x s₀) r ≤ (1 / 2) * error s₀ r) ∧
    -- (VI) Monotone descent:
    (iterate 2 x s₀ (n + 1) ≤ iterate 2 x s₀ n) ∧
    -- (VII) Quadratic convergence:
    (|pandrosion_map 2 x s₀ - r| ≤ (s₀ - r) ^ 2 / (2 * r)) ∧
    -- (VIII) Strict Pell-cone preservation:
    ((pell_iterate X (A, B) n).1 ^ 2 > X * (pell_iterate X (A, B) n).2 ^ 2) ∧
    -- (IX) Finite Voronoï openness:
    (IsOpen {p : ℝ × ℝ | ∀ a ∈ anchors,
      (p.1 - r0.1) ^ 2 + (p.2 - r0.2) ^ 2 <
      (p.1 - a.1) ^ 2 + (p.2 - a.2) ^ 2}) ∧
    -- (X) Pell orbit strict growth:
    ((pell_iterate X (A, B) n).1 < (pell_iterate X (A, B) (n + 1)).1) ∧
    -- (XI) Pairwise Banach contraction:
    (|pandrosion_map 2 x s - pandrosion_map 2 x t| ≤ (1 / 2) * |s - t|) := by
  have hX0 : 0 ≤ X := by linarith
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact global_convergence_p2 x r s₀ hx hr h_s₀ h_root n
  · exact pell_cone_invariant X hX0 (A, B) (by
      show A ^ 2 ≥ X * B ^ 2; linarith) n
  · exact voronoi_strict_open r1x r1y r2x r2y
  · exact pandrosion_map_lipschitz_half_p2 x r s₀ hx hr h_root h_s₀
  · exact iterate_monotone_descent_p2 x r s₀ hr h_s₀ h_root n
  · exact pandrosion_quadratic_convergence_p2 x r s₀ hr (by linarith) h_root h_s₀
  · exact pell_cone_invariant_strict X hX0 (A, B) h_cone_strict n
  · exact voronoi_cell_finite_open r0 anchors
  · exact (pell_iterate_strict_growth X hX (A, B) hA hB n).1
  · exact pandrosion_map_lipschitz_pairwise_p2 x r s t hr h_root hs ht

end Advanced2

end Pandrosion
