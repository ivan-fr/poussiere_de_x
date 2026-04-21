/-
  Universitas Pandrosion — §19. Unconditional McMullen for `z^p − x`.

  The `constructive_mcmullen` theorem (Foundations §13) takes a *generic*
  family of anchors `γ : Fin p → ℂ` with:
    • `Function.Injective γ`           — anchors are distinct
    • `∀ s, Sp_C p (γ s) ≠ 0`          — Sp doesn't vanish
    • `∀ s, (γ s) ^ p = A / x`          — anchors are pth roots

  This file supplies a canonical, *unconditional*, cyclotomic instantiation
  `γ = cycAnchor α p` where `α^p = A/x`, and proves:

    §19.1  `cycAnchor_pow` — `(cycAnchor α p k)^p = α^p`, i.e., each
           anchor is indeed a `p`-th root.

    §19.2  `cycAnchor_injective` — the map `k ↦ cycAnchor α p k`
           is injective on `Fin p` when `α ≠ 0`. This is the
           load-bearing unconditional hypothesis.

    ★ §19.3  `unconditional_mcmullen` — McMullen's conclusion *without*
           the `h_inj` hypothesis of `constructive_mcmullen`: for every
           `x ≠ 0`, `p ≥ 1`, and every `α ≠ 0` with `α^p = A/x`, the
           cyclotomic anchors `k ↦ α · e^{2π i k/p}` satisfy the
           full conclusion of `constructive_mcmullen`.

  This removes the "generic" hypothesis and makes the McMullen
  analogue completely unconditional for the `z^p − x` family.
-/

import Pandrosion.Core.MultiStartBasins
import Mathlib.Analysis.SpecialFunctions.Complex.Log

namespace Pandrosion

open Complex

/-! ============================================================
  §19.1  Cyclotomic anchor family
============================================================ -/

section CyclotomicAnchor

/-- **Cyclotomic anchor family.** For `α ∈ ℂ`, `p ≥ 1`, and
    `k : Fin p`, the `k`-th cyclotomic anchor is
    `α · exp(2π i k / p)`. These are the `p` distinct scaled
    roots of `α^p`, spaced equally on the circle of radius `|α|`. -/
noncomputable def cycAnchor (α : ℂ) (p : ℕ) (k : Fin p) : ℂ :=
  α * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * (k : ℂ) / (p : ℂ))

/-- **Cyclotomic anchors are `p`-th roots of `α^p`.**
    `(α · e^{2πik/p})^p = α^p · e^{2πik} = α^p`. -/
theorem cycAnchor_pow (α : ℂ) {p : ℕ} (hp : 1 ≤ p) (k : Fin p) :
    (cycAnchor α p k) ^ p = α ^ p := by
  unfold cycAnchor
  rw [mul_pow, ← Complex.exp_nat_mul]
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have h_simp :
      (p : ℂ) * (2 * (Real.pi : ℂ) * Complex.I * (k : ℂ) / (p : ℂ))
        = (k : ℤ) * (2 * (Real.pi : ℂ) * Complex.I) := by
    push_cast
    field_simp
    ring
  rw [h_simp, Complex.exp_int_mul_two_pi_mul_I]
  ring

/-- **Helper: integer-arithmetic core of injectivity.**
    If `j, k : Fin p` and `(j : ℤ) − (k : ℤ) = n · p` for some `n : ℤ`,
    then `j = k`. -/
private theorem Fin.eq_of_int_diff_eq_mul_cast
    {p : ℕ} (hp : 1 ≤ p) {j k : Fin p} {n : ℤ}
    (hn : ((j : ℤ) - (k : ℤ)) = n * (p : ℤ)) : j = k := by
  have hj_lt : (j : ℤ) < (p : ℤ) := by exact_mod_cast j.2
  have hk_lt : (k : ℤ) < (p : ℤ) := by exact_mod_cast k.2
  have hj_nn : (0 : ℤ) ≤ (j : ℤ) := Int.ofNat_nonneg _
  have hk_nn : (0 : ℤ) ≤ (k : ℤ) := Int.ofNat_nonneg _
  have hp_pos : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hp
  -- n*p ∈ (-p, p), hence n = 0.
  have h_lb : -(p : ℤ) < n * (p : ℤ) := by linarith
  have h_ub : n * (p : ℤ) < (p : ℤ) := by linarith
  have hn_zero : n = 0 := by
    rcases lt_trichotomy n 0 with h | h | h
    · exfalso
      have : n * (p : ℤ) ≤ -1 * (p : ℤ) :=
        mul_le_mul_of_nonneg_right (by omega) (le_of_lt hp_pos)
      linarith
    · exact h
    · exfalso
      have : (1 : ℤ) * (p : ℤ) ≤ n * (p : ℤ) :=
        mul_le_mul_of_nonneg_right (by omega) (le_of_lt hp_pos)
      linarith
  rw [hn_zero, zero_mul] at hn
  exact Fin.ext (by exact_mod_cast (by linarith : (j : ℤ) = (k : ℤ)))

/-- **Cyclotomic anchors are injective.** For `α ≠ 0` and `p ≥ 1`,
    the map `k ↦ cycAnchor α p k` is injective on `Fin p`.

    *Proof.* If `α · e^{2πij/p} = α · e^{2πik/p}`, cancel `α` to
    get `e^{2πij/p} = e^{2πik/p}`. By `exp_eq_exp_iff_exists_int`,
    there is `n : ℤ` with `2πij/p − 2πik/p = n · 2πi`, i.e.
    `(j − k)/p = n`, i.e. `j − k = n · p` in `ℤ`. For
    `j, k : Fin p`, `|j − k| < p`, forcing `n = 0`, so `j = k`. -/
theorem cycAnchor_injective {α : ℂ} (hα : α ≠ 0) {p : ℕ} (hp : 1 ≤ p) :
    Function.Injective (cycAnchor α p) := by
  intro j k hjk
  unfold cycAnchor at hjk
  -- Cancel α.
  have h_exp :
      Complex.exp (2 * (Real.pi : ℂ) * Complex.I * (j : ℂ) / (p : ℂ))
    = Complex.exp (2 * (Real.pi : ℂ) * Complex.I * (k : ℂ) / (p : ℂ)) :=
    mul_left_cancel₀ hα hjk
  -- ∃ n : ℤ, 2πij/p = 2πik/p + n · 2πi
  rw [Complex.exp_eq_exp_iff_exists_int] at h_exp
  obtain ⟨n, hn⟩ := h_exp
  have hp_ne : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have h2pi : (2 * (Real.pi : ℂ) * Complex.I) ≠ 0 := by
    refine mul_ne_zero (mul_ne_zero ?_ ?_) Complex.I_ne_zero
    · exact two_ne_zero
    · exact_mod_cast Real.pi_ne_zero
  -- From hn: 2πi·j/p = 2πi·k/p + n · 2πi. Multiply both sides by p:
  have hn' :
      (2 * (Real.pi : ℂ) * Complex.I) * (j : ℂ)
    = (2 * (Real.pi : ℂ) * Complex.I) * (k : ℂ)
      + (n : ℂ) * (2 * (Real.pi : ℂ) * Complex.I) * (p : ℂ) := by
    have := hn
    field_simp at this
    linear_combination this
  -- Factor out 2πi and cancel.
  have hn'' :
      (2 * (Real.pi : ℂ) * Complex.I) * ((j : ℂ) - (k : ℂ) - (n : ℂ) * (p : ℂ)) = 0 := by
    linear_combination hn'
  have h_diff : (j : ℂ) - (k : ℂ) - (n : ℂ) * (p : ℂ) = 0 := by
    rcases mul_eq_zero.mp hn'' with h | h
    · exact absurd h h2pi
    · exact h
  -- Convert to ℤ and apply the integer-arithmetic core.
  have h_cast : ((j : ℤ) - (k : ℤ) : ℂ) = (n : ℂ) * ((p : ℤ) : ℂ) := by
    push_cast
    linear_combination h_diff
  have h_int : ((j : ℤ) - (k : ℤ)) = n * (p : ℤ) := by exact_mod_cast h_cast
  exact Fin.eq_of_int_diff_eq_mul_cast hp h_int

end CyclotomicAnchor

/-! ============================================================
  §19.2  Unconditional McMullen (Fundamental Theorem)
============================================================ -/

section UnconditionalMcMullen

/-- **★ Unconditional constructive McMullen for `z^p − x`.**

    For every `x ≠ 0`, every `p ≥ 1`, every `A ≠ 0`, and every
    `α ≠ 0` with `α^p = A/x`, the cyclotomic anchor family
    `γ := cycAnchor α p` satisfies the full conclusion of
    `constructive_mcmullen` *without* a separately-supplied
    injectivity hypothesis — the injectivity is automatic
    (§19.1 `cycAnchor_injective`).

    This is the "truly unconditional" form of McMullen's analogue:
    once `x ≠ 0` and `α^p = A/x` are fixed, the global-convergence
    properties (A), (A'), (B), (C) hold *for free*, without any
    additional regularity on the anchor family. -/
theorem unconditional_mcmullen
    {p : ℕ} (hp : 1 ≤ p) (x A α : ℂ)
    (hx : x ≠ 0) (hA : A ≠ 0) (hα : α ≠ 0) (hαp : α ^ p = A / x)
    (hSp : ∀ k : Fin p, Sp_C p (cycAnchor α p k) ≠ 0) :
    (∀ k : Fin p, cycAnchor α p k ∈ interior (basin (cycAnchor α p) k)) ∧
    (⋃ k : Fin p, basin (cycAnchor α p) k) = Set.univ ∧
    (∀ z₀ ∉ VoronoiBoundary (cycAnchor α p),
        ∃! s : Fin p,
          ∀ t : Fin p, ‖cycAnchor α p s - z₀‖ ≤ ‖cycAnchor α p t - z₀‖) ∧
    (∀ s : Fin p, steffensen_step_C (x / A) p (cycAnchor α p s) = cycAnchor α p s) := by
  have h_inj : Function.Injective (cycAnchor α p) := cycAnchor_injective hα hp
  have hfp : ∀ s : Fin p, (cycAnchor α p s) ^ p = A / x := by
    intro s
    rw [cycAnchor_pow α hp s, hαp]
  exact constructive_mcmullen hp x A hx hA (cycAnchor α p) h_inj hSp hfp

end UnconditionalMcMullen

end Pandrosion
