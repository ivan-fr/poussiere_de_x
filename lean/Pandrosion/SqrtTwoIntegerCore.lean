/-
  Universitas Pandrosion — Lean 4 Formalization
  SQUARE ROOT OF 2 IS IRRATIONAL — INTEGER CORE (no scaffold, no sorry)

  Real theorem: for all natural numbers p, q with q > 0,
      p^2 ≠ 2 · q^2.
  This is the classical integer core of √2's irrationality:
  if √2 = p/q with q > 0, then p^2 = 2·q^2, which this theorem refutes.

  The proof uses no structure-field trick: it appeals only to the
  2-adic divisibility chain 2 ∣ p^2 → 2 ∣ p via primality of 2,
  iterated through Fermat's infinite descent on q.

  This mirrors CbrtTwoIntegerCore for the quadratic case, closing
  the √2/∛2 pair in the Pandrosion integer-core corpus.
-/
import Mathlib.Data.Nat.Prime
import Mathlib.Tactic

namespace Pandrosion

/-- **Integer core of √2 irrationality.**
    For every natural number `q > 0` and every natural `p`, `p^2 ≠ 2 · q^2`.

    Proof sketch (infinite descent on q):
    * `2 ∣ 2·q^2 = p^2`, and `2` prime gives `2 ∣ p`.
    * Write `p = 2p'`. Then `4·p'^2 = 2·q^2`, i.e. `q^2 = 2·p'^2`.
    * Hence `2 ∣ q^2`, so `2 ∣ q`. Write `q = 2q'`.
    * Then `4·q'^2 = 2·p'^2`, i.e. `p'^2 = 2·q'^2`.
    * Strong induction on `q` applies to `q' < q`, a contradiction. -/
theorem sq_ne_two_times_sq :
    ∀ q : ℕ, 0 < q → ∀ p : ℕ, p ^ 2 ≠ 2 * q ^ 2 := by
  intro q
  induction q using Nat.strong_induction_on with
  | _ q ih =>
    intro hq p heq
    -- `2 ∣ p^2` follows from `heq : p^2 = 2 * q^2`.
    have h2_dvd_p2 : (2 : ℕ) ∣ p ^ 2 := ⟨q ^ 2, heq⟩
    -- Primality of 2 transfers divisibility from `p^2` to `p`.
    have h2_p : (2 : ℕ) ∣ p := Nat.prime_two.dvd_of_dvd_pow h2_dvd_p2
    obtain ⟨p', rfl⟩ := h2_p
    -- Expand `(2 p')^2 = 4 p'^2`.
    have expand_p : (2 * p') ^ 2 = 4 * p' ^ 2 := by ring
    rw [expand_p] at heq
    -- From `4 · p'^2 = 2 · q^2`, derive `q^2 = 2 · p'^2`.
    have hq2_eq : q ^ 2 = 2 * p' ^ 2 := by omega
    -- Hence `2 ∣ q^2`, and primality of 2 gives `2 ∣ q`.
    have h2_dvd_q2 : (2 : ℕ) ∣ q ^ 2 := ⟨p' ^ 2, by omega⟩
    have h2_q : (2 : ℕ) ∣ q := Nat.prime_two.dvd_of_dvd_pow h2_dvd_q2
    obtain ⟨q', rfl⟩ := h2_q
    -- Expand `(2 q')^2 = 4 q'^2`.
    have expand_q : (2 * q') ^ 2 = 4 * q' ^ 2 := by ring
    rw [expand_q] at hq2_eq
    -- From `4 · q'^2 = 2 · p'^2`, derive `p'^2 = 2 · q'^2`.
    have hpc : p' ^ 2 = 2 * q' ^ 2 := by omega
    -- Strict descent: `q' < 2·q' = q`, and `q' > 0` since `q > 0`.
    have hq'_pos : 0 < q' := by omega
    have hq'_lt : q' < 2 * q' := by omega
    exact ih q' hq'_lt hq'_pos p' hpc

/-- **Convenient form.** For every natural `p` and every positive natural `q`,
    `p^2 ≠ 2 · q^2`. -/
theorem sq_ne_two_times_sq_of_pos (p q : ℕ) (hq : 0 < q) :
    p ^ 2 ≠ 2 * q ^ 2 :=
  sq_ne_two_times_sq q hq p

/-- **Symmetric form.** -/
theorem two_times_sq_ne_sq (p q : ℕ) (hq : 0 < q) :
    2 * q ^ 2 ≠ p ^ 2 :=
  fun h => sq_ne_two_times_sq_of_pos p q hq h.symm

/-- **Positive-positive form.** -/
theorem no_pos_square_doubles_square (p q : ℕ) (_ : 0 < p) (hq : 0 < q) :
    p ^ 2 ≠ 2 * q ^ 2 :=
  sq_ne_two_times_sq_of_pos p q hq

end Pandrosion
