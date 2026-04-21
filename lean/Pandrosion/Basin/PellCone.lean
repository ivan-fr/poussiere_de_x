/-
  Universitas Pandrosion — Basin / PellCone
  Split from Advanced2.lean: §II, VIII, X
  (Pell cone preservation, strict cone, strict growth)
-/
import Mathlib.Tactic
import Pandrosion.Core
import Pandrosion.Advanced

namespace Pandrosion

section PellCone

/-! ## §II. Pell-Cone Preservation -/

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

/-- **Pell-cone invariance for the full orbit.** -/
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

/-- **Sign-propagation corollary: trivial Pell solutions stay trivial.** -/
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

/-! ## §VIII. Strict Pell-Cone Preservation -/

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

/-! ## §X. Pell Orbit Strict Growth -/

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

end PellCone

end Pandrosion
