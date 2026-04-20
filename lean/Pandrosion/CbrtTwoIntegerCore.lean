/-
  Universitas Pandrosion — Lean 4 Formalization
  CUBE ROOT OF 2 IS IRRATIONAL — INTEGER CORE (no scaffold, no sorry)

  Real theorem: for all natural numbers p, q with q > 0,
      p^3 ≠ 2 · q^3.
  This is the integer core of the classical result that ∛2 is irrational:
  if ∛2 = p/q in lowest terms with q > 0, then p^3 = 2·q^3, which this
  theorem refutes.

  The proof uses no structure-field trick: it appeals only to the
  2-adic divisibility chain 2 ∣ p^3 → 2 ∣ p via primality of 2,
  iterated through Fermat's infinite descent on q.

  This complements LogTwoThreeIrrational by attacking the cube-root
  case directly — the central object of the Pandrosion ∛X framework.
-/
import Mathlib.Data.Nat.Prime
import Mathlib.Tactic

namespace Pandrosion

/-- **Integer core of ∛2 irrationality.**
    For every natural number `q > 0` and every natural `p`, `p^3 ≠ 2 · q^3`.

    Proof sketch (infinite descent on q):
    * `2 ∣ 2·q^3 = p^3`, and `2` prime gives `2 ∣ p`.
    * Write `p = 2p'`. Then `8·p'^3 = 2·q^3`, i.e. `q^3 = 4·p'^3`.
    * Hence `2 ∣ q^3`, so `2 ∣ q`. Write `q = 2q'`.
    * Then `8·q'^3 = 4·p'^3`, i.e. `p'^3 = 2·q'^3`.
    * Strong induction on `q` applies to `q' < q`, a contradiction. -/
theorem cube_ne_two_times_cube :
    ∀ q : ℕ, 0 < q → ∀ p : ℕ, p ^ 3 ≠ 2 * q ^ 3 := by
  intro q
  induction q using Nat.strong_induction_on with
  | _ q ih =>
    intro hq p heq
    -- `2 ∣ p^3` follows from `heq : p^3 = 2 * q^3`.
    have h2_dvd_p3 : (2 : ℕ) ∣ p ^ 3 := ⟨q ^ 3, heq⟩
    -- Primality of 2 transfers divisibility from `p^3` to `p`.
    have h2_p : (2 : ℕ) ∣ p := Nat.prime_two.dvd_of_dvd_pow h2_dvd_p3
    obtain ⟨p', rfl⟩ := h2_p
    -- Expand `(2 p')^3 = 8 p'^3`.
    have expand_p : (2 * p') ^ 3 = 8 * p' ^ 3 := by ring
    rw [expand_p] at heq
    -- From `8 · p'^3 = 2 · q^3`, derive `q^3 = 4 · p'^3`.
    have hq3_eq : q ^ 3 = 4 * p' ^ 3 := by omega
    -- Hence `2 ∣ q^3`, and primality of 2 gives `2 ∣ q`.
    have h2_dvd_q3 : (2 : ℕ) ∣ q ^ 3 := ⟨2 * p' ^ 3, by omega⟩
    have h2_q : (2 : ℕ) ∣ q := Nat.prime_two.dvd_of_dvd_pow h2_dvd_q3
    obtain ⟨q', rfl⟩ := h2_q
    -- Expand `(2 q')^3 = 8 q'^3`.
    have expand_q : (2 * q') ^ 3 = 8 * q' ^ 3 := by ring
    rw [expand_q] at hq3_eq
    -- From `8 · q'^3 = 4 · p'^3`, derive `p'^3 = 2 · q'^3`.
    have hpc : p' ^ 3 = 2 * q' ^ 3 := by omega
    -- Strict descent: `q' < 2·q' = q`, and `q' > 0` since `q > 0`.
    have hq'_pos : 0 < q' := by omega
    have hq'_lt : q' < 2 * q' := by omega
    exact ih q' hq'_lt hq'_pos p' hpc

/-- **Convenient form.** For every natural `p` and every positive natural `q`,
    `p^3 ≠ 2 · q^3`. -/
theorem cube_ne_two_times_cube_of_pos (p q : ℕ) (hq : 0 < q) :
    p ^ 3 ≠ 2 * q ^ 3 :=
  cube_ne_two_times_cube q hq p

/-- **Symmetric form.** -/
theorem two_times_cube_ne_cube (p q : ℕ) (hq : 0 < q) :
    2 * q ^ 3 ≠ p ^ 3 :=
  fun h => cube_ne_two_times_cube_of_pos p q hq h.symm

/-- **Positive-positive form.** -/
theorem no_pos_cube_doubles_cube (p q : ℕ) (_ : 0 < p) (hq : 0 < q) :
    p ^ 3 ≠ 2 * q ^ 3 :=
  cube_ne_two_times_cube_of_pos p q hq

end Pandrosion
