/-
  Universitas Pandrosion — Legacy / Basic.

  Ported from the `main_previous` branch (pre-rewrite). All definitions
  of `Sp`, `pandrosion_h`, the fixed-point iff, etc. are re-used from
  `Pandrosion.Core.Foundations`. This file only adds legacy lemmas that
  the clean rewrite did not port:

    • `Sp_one_eq` — `S_1(s) = 1` (helper).
    • `h_lt_one` — `h(s) < 1` for `x > 1`, `s ≥ 0`.
    • `h_pos`    — `h(s) > 0` for `x > 1`, `p ≥ 2`, `s ∈ [0,1]`.
    • `iter`, `orbit_in_interval` — `p ≥ 2` orbit stays in `(0,1)`.
    • `aitken_perfect_extrapolation` — Aitken Δ² is *exact* on a
      geometric error sequence (one-step convergence recovery).
-/

import Pandrosion.Core.Foundations
import Mathlib.Tactic

namespace Pandrosion

open Finset BigOperators

/-! ## §L1. `Sp` at `p = 1` -/

/-- `S_1(s) = 1`. -/
theorem Sp_one_eq (s : ℝ) : Sp 1 s = 1 := by
  simp [Sp]

/-! ## §L2. Range invariance of `h` on `[0,1]` for `x > 1` -/

/-- `h(s) < 1` for `s ≥ 0`, `x > 1`. Since `(x-1)/(x·S_p(s)) > 0`, we
    have `h(s) = 1 − positive < 1`. -/
theorem h_lt_one (x : ℝ) (hx : x > 1) (p : ℕ) (hp : p ≥ 1) (s : ℝ) (hs : 0 ≤ s) :
    pandrosion_h x p s < 1 := by
  unfold pandrosion_h
  simp only [sub_lt_self_iff]
  exact div_pos (by linarith) (mul_pos (by linarith) (Sp_pos p hp s hs))

/-- `h(s) > 0` for `s ∈ [0,1]`, `x > 1`, `p ≥ 2`. Requires
    `(x−1)/(x·S_p(s)) < 1`, which follows from `S_p(s) ≥ 1`. -/
theorem h_pos (x : ℝ) (hx : x > 1) (p : ℕ) (hp : p ≥ 2) (s : ℝ)
    (hs : 0 ≤ s) (_hs1 : s ≤ 1) :
    pandrosion_h x p s > 0 := by
  unfold pandrosion_h
  have hS := Sp_pos p (by omega) s hs
  rw [gt_iff_lt, sub_pos, div_lt_one (mul_pos (by linarith) hS)]
  have hS1 : Sp p s ≥ 1 := by
    have : Sp p s ≥ Sp 1 s := by
      unfold Sp
      exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.range_mono (by omega))
        (fun k _ _ => pow_nonneg hs k)
    linarith [Sp_one_eq s]
  nlinarith

/-! ## §L3. Orbit invariance -/

/-- The n-th iterate of the Pandrosion map. -/
noncomputable def iter (x : ℝ) (p : ℕ) : ℕ → ℝ → ℝ
  | 0, s => s
  | n + 1, s => pandrosion_h x p (iter x p n s)

/-- For `p ≥ 2`, `x > 1`, `s₀ ∈ [0,1]`: every iterate `n ≥ 1` is in `(0,1)`. -/
theorem orbit_in_interval (x : ℝ) (hx : x > 1) (p : ℕ) (hp : p ≥ 2)
    (s₀ : ℝ) (hs₀ : 0 ≤ s₀) (hs₀1 : s₀ ≤ 1)
    (n : ℕ) (hn : n ≥ 1) :
    0 < iter x p n s₀ ∧ iter x p n s₀ < 1 := by
  induction n with
  | zero => exfalso; exact Nat.not_succ_le_zero 0 hn
  | succ m ih =>
    unfold iter
    cases m with
    | zero =>
        simp [iter]
        exact ⟨h_pos x hx p hp s₀ hs₀ hs₀1,
               h_lt_one x hx p (by omega) s₀ hs₀⟩
    | succ k =>
        have ⟨hlo, hhi⟩ := ih (by omega)
        exact ⟨h_pos x hx p hp _ (le_of_lt hlo) (le_of_lt hhi),
               h_lt_one x hx p (by omega) _ (le_of_lt hlo)⟩

/-! ## §L4. Exact Aitken extrapolation -/

/-- **Aitken's Δ² is exact on a pure geometric error.**
    If `s₁ − r = λ(s − r)` and `s₂ − r = λ(s₁ − r)` with `λ ≠ 1` and
    `s ≠ r`, then the Aitken-Steffensen formula
    `s − (s₁ − s)² / (s₂ − 2s₁ + s)` recovers the limit `r` in a
    single step. The denominator is non-zero on the non-degenerate
    branch `(s − r)(λ − 1)² ≠ 0`. -/
theorem aitken_perfect_extrapolation (r s lam : ℝ) (h_lam : lam ≠ 1) (hs : s ≠ r) :
    let s1 := r + lam * (s - r)
    let s2 := r + lam * (s1 - r)
    let denom := s2 - 2 * s1 + s
    denom ≠ 0 ∧ s - (s1 - s) ^ 2 / denom = r := by
  intros s1 s2 denom
  have h_denom : denom = (s - r) * (lam - 1) ^ 2 := by
    dsimp [denom, s2, s1]; ring
  have h_sr : s - r ≠ 0 := sub_ne_zero.mpr hs
  have h_lam1 : lam - 1 ≠ 0 := sub_ne_zero.mpr h_lam
  have hn_lam1 : (lam - 1) ^ 2 ≠ 0 := pow_ne_zero 2 h_lam1
  have hz : denom ≠ 0 := by rw [h_denom]; exact mul_ne_zero h_sr hn_lam1
  refine ⟨hz, ?_⟩
  have h_num : (s1 - s) ^ 2 = (s - r) ^ 2 * (lam - 1) ^ 2 := by
    dsimp [s1]; ring
  rw [h_num, h_denom]
  have h_decomp : (s - r) ^ 2 * (lam - 1) ^ 2
      = (s - r) * ((s - r) * (lam - 1) ^ 2) := by ring
  rw [h_decomp, mul_div_cancel_right₀ (s - r)
        (mul_ne_zero h_sr hn_lam1)]
  ring

end Pandrosion
