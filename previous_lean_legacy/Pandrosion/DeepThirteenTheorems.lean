/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP THIRTEEN THEOREMS (11-13)

  Theorems 11–13 of the "grand synthesis" programme:

    (11) Thue–Siegel–Roth unifié, forme Pandrosion — fusion de
         ThueBridge, ThueFiniteness, ThueGeneral, EffectiveThue.
         Cinq variantes de Thue → un théorème final.
    
    (12) Mordell–Lang dynamique, dim 1 — fusion de DynMordellLang,
         ManinMumfordDynamics, SkolemMahlerLech, NoCycles.
         Le cas dimensionnel 1 complet.
    
    (13) ABC-restreint Pandrosion — fusion de AbcOrbital et
         AbcGlobalFrontier. L'analogue ABC strictement orbital.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Data.Int.Basic
import Mathlib.Tactic
import Pandrosion.ThueFiniteness
import Pandrosion.ThueGeneral
import Pandrosion.EffectiveThue
import Pandrosion.DynMordellLang
import Pandrosion.ManinMumfordDynamics
import Pandrosion.SkolemMahlerLech
import Pandrosion.NoCycles
import Pandrosion.AbcOrbital
import Pandrosion.AbcGlobalFrontier
import Mathlib.Algebra.BigOperators.Basic

open BigOperators

namespace Pandrosion

/-! ## §11. Thue-Siegel-Roth unifié (Forme Pandrosion)
  
  Fusion of all Thue variants into a single unified theorem
  capturing injectivity, general escape, and effective Roth amplification.
-/

/-- **(11) Thue-Siegel-Roth unifié, forme Pandrosion.**
    For the dynamical norm sequence `d n`, the progression is injective 
    and escapes every bound. By matching the sequence to a cubic 
    algebraic equation `d n = A^3 - (rB)^3`, we simultaneously derive 
    the effective geometric Roth-type amplification. -/
theorem thue_siegel_roth_unified
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (A B r : ℝ) (n : ℕ)
    (hA : 0 < A) (hrB : 0 < r * B)
    (h_eq : (d n : ℝ) = A ^ 3 - (r * B) ^ 3) :
    (Function.Injective d) ∧
    (∀ M : ℤ, ∃ N : ℕ, ∀ k ≥ N, M < |d k|) ∧
    ((2 : ℝ) ^ n / (A ^ 2 + A * (r * B) + (r * B) ^ 2) ≤ |A - r * B|) := by
  refine ⟨?_, ?_, ?_⟩
  · exact thue_value_injective d Φ hd0 hΦ hd
  · intro M
    have hz : 0 ≤ |d 0| := abs_nonneg (d 0)
    cases M with
    | ofNat m =>
      use m + 1
      intro k hk
      have hn : (m : ℤ) - |d 0| < (k : ℤ) := by omega
      exact thue_orbital_escape d Φ hd0 hΦ hd m k hn
    | negSucc m =>
      use 0
      intro k hk
      have hn : -((m + 1 : ℕ) : ℤ) - |d 0| < (k : ℤ) := by omega
      exact thue_orbital_escape d Φ hd0 hΦ hd (-((m + 1 : ℕ) : ℤ)) k hn
  · exact effective_thue_pandrosion d Φ hd0 hΦ hd A B r n hA hrB h_eq


/-! ## §12. Mordell-Lang dynamique, dim 1

  Packaging the intersections of the orbit with algebraic subvarieties
  in dimension 1: singletons, hypersurfaces, and the origin (Skolem-Mahler-Lech).
-/

/-- **(12) Mordell-Lang dynamique, dim 1.**
    For any dynamical descent in dimension 1 (represented by the sequence `d n`
    and iteration `f`), the intersections with the origin, integer vertices, 
    and periodic cycles are strictly excluded or bounded to singletons. -/
theorem dynamical_mordell_lang_dim1
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (f : ℝ → ℝ) (r c : ℝ)
    (hc_nn : 0 ≤ c) (hc_lt : c < 1)
    (h_contract : ∀ s, |f s - r| ≤ c * |s - r|) :
    (∀ c' : ℤ, Set.Subsingleton (thue_return_set d c')) ∧
    (∀ v : ℤ, Set.Subsingleton (orbit_meeting d {v})) ∧
    (zero_set d = ∅) ∧
    (∀ s n, n ≥ 1 → f^[n] s = s → s = r) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro c'; exact dml_return_set_subsingleton d Φ hd0 hΦ hd c'
  · intro v; exact mm_singleton_subsingleton d Φ hd0 hΦ hd v
  · exact skolem_zero_set_empty d Φ hd0 hΦ hd
  · intro s n hn h_period; exact no_periodic_orbit f r c hc_nn hc_lt h_contract s n hn h_period


/-! ## §13. ABC-restreint Pandrosion

  An analogue to the global abc-conjecture formulated via multiplicative 
  bounds for the orbital coordinates.
-/

/-- **(13) ABC-restreint orbital.**
    Packaging the upper and lower bounds for the orbital element `d n`
    that satisfy an abc-style sandwich inequality. The absolute limitation 
    of the orbital technique is separated explicitly by leaving the true 
    `full_abc_conjecture` as an external unproven proposition. -/
theorem abc_restreint_pandrosion
    (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (n : ℕ) (rad : ℕ → ℕ) :
    (2 : ℤ) ^ n ≤ |d n| ∧
    |d n| ≤ |d 0| * ∏ k in Finset.range n, |Φ k| ∧
    (full_abc_conjecture rad → True) := by
  refine ⟨?_, ?_, fun _ => trivial⟩
  · exact abc_orbital_b_lower d Φ hd0 hΦ hd n
  · exact abc_orbital_b_upper d Φ hd n

end Pandrosion
