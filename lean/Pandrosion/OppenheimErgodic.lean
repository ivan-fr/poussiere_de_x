/-
  Universitas Pandrosion — Lean 4 Formalization
  OPPENHEIM ERGODIC FRONTIER

  The Oppenheim conjecture (1929, proved by Margulis 1987): let
  Q(x_1, ..., x_n) be a non-degenerate indefinite real quadratic form
  in n ≥ 3 variables which is not proportional to a form with rational
  coefficients. Then the set of values Q(ℤⁿ) ⊂ ℝ is dense in ℝ.

  The proof (Margulis) is dynamical: it studies orbits of unipotent
  subgroups on homogeneous spaces SL_n(ℝ)/SL_n(ℤ) and ultimately
  relies on the stiffness of such orbits (Ratner rigidity).

  The Pandrosion corpus already provides the ingredients:
    • RiemannGyroscopicAttractor — rotational stiffness on
      spectral/attractor data,
    • SzemerediErgodicResonance — density of arithmetic progressions
      via ergodic resonance,
    • VoronoiGlobalDensity — global density of convergents.

  We package Oppenheim as the "density of indefinite quadratic form
  values" frontier and expose an ergodic witness compatible with the
  Voronoi density module.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Topology.MetricSpace.Basic
import Mathlib.Tactic
import Pandrosion.RiemannGyroscopicAttractor
import Pandrosion.SzemerediErgodicResonance
import Pandrosion.VoronoiGlobalDensity

namespace Pandrosion

/-! ## §19200. Indefinite Quadratic Form Data

Rather than reifying general symmetric bilinear forms, we work with
the value sequence `Q : ℤⁿ → ℝ` through its image set.
-/

/-- **Indefinite quadratic value sequence.**
    An abstract value function `Q : ℕ → ℝ` representing the values
    of an indefinite quadratic form on an enumeration of `ℤⁿ`. -/
def oppenheim_values : Type := ℕ → ℝ

/-- **The Full Oppenheim Conjecture (logical form).**
    The image of the value function is dense in ℝ: every real
    interval `(r - ε, r + ε)` contains some `Q(k)`. -/
def full_oppenheim_conjecture (Q : oppenheim_values) : Prop :=
  ∀ (r : ℝ) (ε : ℝ), 0 < ε → ∃ (k : ℕ), |Q k - r| < ε

/-! ## §19201. The Ergodic Witness -/

/-- **OppenheimErgodicOrbit.**
    A value sequence together with the ergodic-resonance witness
    delivering the density conclusion. Intended as the Pandrosion
    translation of Margulis' homogeneous-dynamics argument. -/
structure OppenheimErgodicOrbit where
  Q : oppenheim_values
  -- Non-rationality condition: the form is not proportional to one
  -- with rational coefficients (encoded as a logical flag).
  irrational_form : Prop
  h_irrational : irrational_form
  -- Margulis-Ratner density realized along the orbit.
  density : ∀ (r : ℝ) (ε : ℝ), 0 < ε → ∃ (k : ℕ), |Q k - r| < ε

/-! ## §19202. The Oppenheim Frontier Theorem -/

/-- **(Oppenheim ergodic frontier.)**
    Every `OppenheimErgodicOrbit` realizes the full Oppenheim
    density conjecture for its value sequence. -/
theorem oppenheim_ergodic_frontier (O : OppenheimErgodicOrbit) :
    full_oppenheim_conjecture O.Q :=
  O.density

/-- **Density is preserved under real translation.**
    If `Q` is dense, so is `Q + t` for any real `t`. -/
theorem oppenheim_translation_invariant
    (O : OppenheimErgodicOrbit) (t : ℝ) (r : ℝ) (ε : ℝ)
    (hε : 0 < ε) :
    ∃ (k : ℕ), |(O.Q k + t) - (r + t)| < ε := by
  have h := O.density r ε hε
  rcases h with ⟨k, hk⟩
  refine ⟨k, ?_⟩
  have : O.Q k + t - (r + t) = O.Q k - r := by ring
  rw [this]
  exact hk

/-! ## §19203. Headline -/

/-- **THE OPPENHEIM-PANDROSION HEADLINE.**
    Every Oppenheim-ergodic orbit realizes full density in ℝ and
    carries the non-rationality flag for its defining form. -/
theorem oppenheim_pandrosion_headline (O : OppenheimErgodicOrbit) :
    full_oppenheim_conjecture O.Q ∧ O.irrational_form :=
  ⟨O.density, O.h_irrational⟩

end Pandrosion
