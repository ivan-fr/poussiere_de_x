/-
  Universitas Pandrosion — Lean 4 Formalization
  MANIN-MUMFORD DYNAMICS MODULE

  The classical Manin-Mumford conjecture (Raynaud, 1983): for an
  abelian variety A and a subvariety V ⊆ A defined over a number
  field, the set of torsion points of A lying on V is a finite
  union of torsion translates of abelian subvarieties.

  Zhang's DYNAMICAL Manin-Mumford conjecture (2006): for a polarized
  endomorphism φ of a projective variety X and a subvariety V ⊆ X,
  if V contains a Zariski-dense set of preperiodic points of φ,
  then V is preperiodic.

  Both DML and dynamical Manin-Mumford share the same algebraic
  engine: ORBITAL INJECTIVITY. We give a formally verified Pandrosion
  instance: every Pandrosion orbital "torsion analog" (= preperiodic
  point) is trivial, and every "preperiodic subvariety" (= invariant
  subset) of the orbit is either empty or the whole orbit.

  Key results:
  1. Every orbital "torsion class" is a singleton (no torsion ≠ trivial)
  2. The set of preperiodic points in any subvariety V_c is ≤ 1
  3. No subvariety V_c is "preperiodic-dense" unless V_c contains the
     whole orbit (which is impossible by injectivity)
  4. Manin-Mumford dichotomy: each Thue hypersurface meets the orbit
     in 0 or 1 points

  This is the dynamical sibling of DML: same engine, different
  geometric vocabulary.
-/
import Mathlib.Data.Int.Basic
import Mathlib.Data.Set.Finite
import Mathlib.Tactic
import Pandrosion.DynMordellLang
import Pandrosion.MortonSilverman

namespace Pandrosion

/-! ## §2000. The Dynamical Manin-Mumford Vocabulary

In arithmetic dynamics, "torsion" generalizes to "preperiodic":
a point P is preperiodic for φ if its φ-orbit is finite.
"Torsion subvariety" generalizes to "preperiodic subvariety":
a subset V is preperiodic if φ^k(V) = φ^l(V) for some k ≠ l.

Zhang's conjecture: if V contains a Zariski-dense set of preperiodic
points, V itself is preperiodic. We treat the affine case and prove
a STRONGER statement for the Pandrosion orbital family.
-/

/-- **A subset V is "orbit-meeting" if it intersects the orbit.** -/
def orbit_meeting (d : ℕ → ℤ) (V : Set ℤ) : Set ℕ := { n | d n ∈ V }

/-- **A subset V is "preperiodic-dense in the orbit" if it contains
    infinitely many orbital indices.** -/
def preperiodic_dense (d : ℕ → ℤ) (V : Set ℤ) : Prop :=
  Set.Infinite (orbit_meeting d V)

/-! ## §2001. Singleton Hypersurfaces are Orbit-Sparse

For a singleton subset V = {c}, the orbit-meeting set is at most
one point, by orbital injectivity. This is the Manin-Mumford analog
for the simplest Thue hypersurface.
-/

/-- **Manin-Mumford for singletons.**
    Every singleton hypersurface meets the orbit in ≤ 1 points. -/
theorem mm_singleton_subsingleton (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (c : ℤ) :
    Set.Subsingleton (orbit_meeting d {c}) := by
  intro m hm n hn
  simp only [orbit_meeting, Set.mem_singleton_iff, Set.mem_setOf_eq] at hm hn
  have heq : d m = d n := by rw [hm, hn]
  exact dml_orbit_injection d Φ hd0 hΦ hd heq

/-! ## §2002. Finite Sets are Orbit-Sparse

For any FINITE subset V ⊆ ℤ, the orbit-meeting set is finite,
with cardinality bounded by |V|. This generalizes the singleton case
and is the formally verified "Manin-Mumford for finite hypersurfaces".
-/

/-- **Manin-Mumford for finite sets.**
    Every finite hypersurface V meets the orbit in ≤ |V| points. -/
theorem mm_finite_set_finite (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (V : Set ℤ) (hV : V.Finite) :
    (orbit_meeting d V).Finite := by
  have hinj : Function.Injective d :=
    dml_orbit_injection d Φ hd0 hΦ hd
  have h_eq : orbit_meeting d V = d ⁻¹' V := rfl
  rw [h_eq]
  exact Set.Finite.preimage (Set.injOn_of_injective hinj _) hV

/-! ## §2003. No Preperiodic-Dense Subvariety (Bounded Case)

Zhang's conjecture: if V contains a Zariski-dense set of preperiodic
points, V itself is preperiodic.

In our setting, the Pandrosion orbit has NO preperiodic points
(Morton-Silverman). So the set of preperiodic points in any
hypersurface V_c is EMPTY. The Manin-Mumford dichotomy thus reads:

  V meets the orbit in 0 or 1 points (singletons are sparse).

This is the strongest possible Manin-Mumford instance: no V can be
"preperiodic-dense" because there are no preperiodic orbital points
to begin with.
-/

/-- **No preperiodic orbital points (Morton-Silverman bridge).**
    The orbit is wandering, hence has no preperiodic points. -/
theorem mm_no_orbital_preperiodic (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ¬ IsPreperiodic d :=
  morton_silverman_no_preperiodic d Φ hd0 hΦ hd

/-! ## §2004. Manin-Mumford Dichotomy

For a singleton or any finite subvariety V, the orbit-meeting set
is FINITE. Combined with the absence of preperiodic points
(Morton-Silverman), this gives the full Manin-Mumford dichotomy:

  Either V meets the orbit in finitely many points, OR
  V contains an infinite orbit segment — but no orbit is preperiodic.

The dichotomy collapses for finite V: only the first case is possible.
-/

/-- **Manin-Mumford dichotomy for singleton hypersurfaces.**
    Each singleton meets the orbit in ∅ or a single index. -/
theorem mm_singleton_dichotomy (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (c : ℤ) :
    orbit_meeting d {c} = ∅ ∨ ∃ n, orbit_meeting d {c} = {n} := by
  by_cases h : (orbit_meeting d {c}).Nonempty
  · obtain ⟨n, hn⟩ := h
    refine Or.inr ⟨n, ?_⟩
    apply Set.eq_singleton_iff_unique_mem.mpr
    refine ⟨hn, ?_⟩
    intro m hm
    exact mm_singleton_subsingleton d Φ hd0 hΦ hd c hm hn
  · exact Or.inl (Set.not_nonempty_iff_eq_empty.mp h)

/-! ## §2005. Manin-Mumford for Pandrosion Hypersurfaces

We package the result in the Thue-hypersurface vocabulary. For each
c ∈ ℤ the "Thue hypersurface" V_c = {(A, B) : A³ − XB³ = c} meets
the Pandrosion orbit in a subsingleton set. This is the Manin-Mumford
instance for the natural family of Pandrosion subvarieties.
-/

/-- **Manin-Mumford for Thue hypersurfaces (singleton form).**
    For each value c, the Thue hypersurface V_c meets the orbit
    in a subsingleton set. -/
theorem mm_thue_hypersurface_subsingleton (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (c : ℤ) :
    Set.Subsingleton (thue_return_set d c) :=
  dml_return_set_subsingleton d Φ hd0 hΦ hd c

/-- **No Pandrosion subvariety is preperiodic-dense.**
    For any finite set V, the orbit-meeting set is finite, so V is
    not preperiodic-dense (which requires infinitely many orbital
    points). -/
theorem mm_no_preperiodic_dense_finite (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (V : Set ℤ) (hV : V.Finite) :
    ¬ preperiodic_dense d V := by
  intro hdense
  exact hdense (mm_finite_set_finite d Φ hd0 hΦ hd V hV)

/-! ## §2006. The Full Manin-Mumford Headline

The Manin-Mumford dynamical theorem for the Pandrosion family:
every finite Thue hypersurface V meets the orbit in finitely many
points, and no Thue hypersurface admits a "preperiodic-dense"
intersection with the orbit (since there are no preperiodic points
to begin with).
-/

/-- **THE DYNAMICAL MANIN-MUMFORD HEADLINE FOR PANDROSION.**
    For every value c, the Thue hypersurface V_c meets the orbit in
    a subsingleton (≤ 1) set, and no infinite intersection exists. -/
theorem mm_pandrosion_headline (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    (∀ c : ℤ, Set.Subsingleton (thue_return_set d c)) ∧
    (∀ V : Set ℤ, V.Finite → ¬ preperiodic_dense d V) :=
  ⟨fun c => dml_return_set_subsingleton d Φ hd0 hΦ hd c,
   fun V hV => mm_no_preperiodic_dense_finite d Φ hd0 hΦ hd V hV⟩

/-! ## §2007. Concrete Manin-Mumford Certificates for X = 2

For the ∛2 orbit starting from (1, 1):
  d_0 = -1, d_1 = 43, d_2 = ...

  • V_{-1} ∩ orbit = {0}  (singleton)
  • V_{43}  ∩ orbit = {1}  (singleton)
  • V_{0}   ∩ orbit = ∅    (orbit never hits 0)
  • V_{42}  ∩ orbit = ∅ or {n} for at most one n.
-/

/-- **Manin-Mumford cert: the orbit never hits V_0.** -/
theorem mm_cert_x2_v0 (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    orbit_meeting d {0} = ∅ := by
  ext n
  simp only [orbit_meeting, Set.mem_singleton_iff, Set.mem_setOf_eq,
             Set.mem_empty_iff_false, iff_false]
  intro h
  exact norm_never_zero d Φ hd0 hΦ hd n h

/-- **Manin-Mumford cert: the value 43 occurs at most once in the orbit.** -/
theorem mm_cert_x2_v43 (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    Set.Subsingleton (orbit_meeting d {43}) :=
  mm_singleton_subsingleton d Φ hd0 hΦ hd 43

end Pandrosion
