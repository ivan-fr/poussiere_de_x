/-
  Universitas Pandrosion — Lean 4 Formalization
  PRIME-ROOT IRRATIONALITY — INTEGER CORE (GENERAL)

  A uniform generalization of:
    • SqrtTwoIntegerCore      (Advanced.lean, §SqrtTwoIntegerCore)
    • CbrtTwoIntegerCore      (Advanced.lean, §CbrtTwoIntegerCore)
    • LogTwoThreeIrrational   (Advanced.lean, §LogTwoThreeIrrational)

  We prove two universal integer-core facts that subsume all three
  cases above:

    (A) ⁿ√p  is irrational for every prime p and every n ≥ 2.
        Integer core:  ∀ a b, b > 0 → a^n ≠ p · b^n.
        Method: infinite descent via `Nat.Prime.dvd_of_dvd_pow`.

    (B) log_p(q) is irrational for distinct primes p ≠ q.
        Integer core:  ∀ a ≥ 1, ∀ b, p^a ≠ q^b.
        Method: primality transfer of divisibility.

  Together, these give irrationality certificates for
      √2, √3, √5, √7, …, ∛2, ∛3, …, ⁿ√p
  and for log₂(3), log₂(5), log₃(5), log₅(7), … — in one shot.
-/
import Mathlib.Data.Nat.Prime
import Mathlib.Tactic

namespace Pandrosion

section PrimeRootIrrationality

/-! ## §A. n-th Root of a Prime is Irrational (Integer Core)

For every prime `p` and every `n ≥ 2`, no pair of naturals `(a, b)`
with `b > 0` can witness `a^n = p · b^n`. The proof is a direct
infinite descent on `b`:

  * `p ∣ a^n  ⇒  p ∣ a`   (prime divides power)
  * substitute `a = p·a'`; then `p^(n-1) · a'^n = b^n`
  * since `n ≥ 2`, we have `p ∣ b^n`, hence `p ∣ b`
  * substitute `b = p·b'`; get `a'^n = p · b'^n` with `b' < b`
  * strong induction on `b` yields the contradiction.
-/

/-- **Core descent lemma.**
    For every prime `p` and every `n ≥ 2`, no positive natural `b`
    and natural `a` can satisfy `a^n = p · b^n`. -/
theorem pow_ne_prime_mul_pow
    (p : ℕ) (hp : p.Prime) (n : ℕ) (hn : 2 ≤ n) :
    ∀ b : ℕ, 0 < b → ∀ a : ℕ, a ^ n ≠ p * b ^ n := by
  intro b
  induction b using Nat.strong_induction_on with
  | _ b ih =>
    intro hb a heq
    -- (1) p prime and p ∣ a^n ⇒ p ∣ a.
    have hp_dvd_an : p ∣ a ^ n := ⟨b ^ n, heq⟩
    have hp_a : p ∣ a := hp.dvd_of_dvd_pow hp_dvd_an
    obtain ⟨a', rfl⟩ := hp_a
    -- (2) Rewrite heq : p^n · a'^n = p · b^n, then strip one p.
    rw [mul_pow] at heq
    have hp_pos : 0 < p := hp.pos
    have h_split : p ^ n = p * p ^ (n - 1) := by
      conv_lhs => rw [show n = (n - 1) + 1 from by omega]
      rw [pow_succ, mul_comm]
    rw [h_split, mul_assoc] at heq
    have heq2 : p ^ (n - 1) * a' ^ n = b ^ n :=
      Nat.eq_of_mul_eq_mul_left hp_pos heq
    -- (3) Since n ≥ 2, p ∣ p^(n-1) · a'^n = b^n, so p ∣ b.
    have hp_dvd_bn : p ∣ b ^ n := by
      rw [← heq2]
      have hsplit2 : p ^ (n - 1) = p * p ^ (n - 2) := by
        conv_lhs => rw [show n - 1 = (n - 2) + 1 from by omega]
        rw [pow_succ, mul_comm]
      rw [hsplit2, mul_assoc]
      exact ⟨p ^ (n - 2) * a' ^ n, rfl⟩
    have hp_b : p ∣ b := hp.dvd_of_dvd_pow hp_dvd_bn
    obtain ⟨b', rfl⟩ := hp_b
    -- (4) From p^(n-1) · a'^n = (p·b')^n = p^n · b'^n, get a'^n = p · b'^n.
    rw [mul_pow] at heq2
    have h_split_n : p ^ n = p ^ (n - 1) * p := by
      conv_lhs => rw [show n = (n - 1) + 1 from by omega, pow_succ]
    rw [h_split_n] at heq2
    have heq3 : a' ^ n = p * b' ^ n := by
      have hpow_pos : 0 < p ^ (n - 1) := pow_pos hp_pos _
      have step : p ^ (n - 1) * a' ^ n = p ^ (n - 1) * (p * b' ^ n) := by
        rw [heq2]; ring
      exact Nat.eq_of_mul_eq_mul_left hpow_pos step
    -- (5) b' > 0 and b' < p·b' = b: strict descent.
    have hp_ge : 2 ≤ p := hp.two_le
    have hb'_pos : 0 < b' := by
      rcases Nat.eq_zero_or_pos b' with h | h
      · simp [h] at hb
      · exact h
    have hb'_lt : b' < p * b' := by nlinarith
    exact ih b' hb'_lt hb'_pos a' heq3

/-- **Convenient symmetric form.** -/
theorem prime_mul_pow_ne_pow
    (p : ℕ) (hp : p.Prime) (n : ℕ) (hn : 2 ≤ n)
    (a b : ℕ) (hb : 0 < b) : p * b ^ n ≠ a ^ n :=
  fun h => pow_ne_prime_mul_pow p hp n hn b hb a h.symm

/-! ## §A.1 Concrete specializations.

Each of the following is a one-line corollary, providing the
integer-core irrationality witness for a specific (prime, degree)
pair. They recover (and strictly extend) the theorems in
`SqrtTwoIntegerCore`, `CbrtTwoIntegerCore`.
-/

/-- **√2 is irrational (integer core).** -/
theorem sqrt_two_integer_core (a b : ℕ) (hb : 0 < b) :
    a ^ 2 ≠ 2 * b ^ 2 :=
  pow_ne_prime_mul_pow 2 Nat.prime_two 2 (le_refl _) b hb a

/-- **∛2 is irrational (integer core).** -/
theorem cbrt_two_integer_core (a b : ℕ) (hb : 0 < b) :
    a ^ 3 ≠ 2 * b ^ 3 :=
  pow_ne_prime_mul_pow 2 Nat.prime_two 3 (by omega) b hb a

/-- **√3 is irrational (integer core).** -/
theorem sqrt_three_integer_core (a b : ℕ) (hb : 0 < b) :
    a ^ 2 ≠ 3 * b ^ 2 :=
  pow_ne_prime_mul_pow 3 Nat.prime_three 2 (le_refl _) b hb a

/-- **√5 is irrational (integer core).** -/
theorem sqrt_five_integer_core (a b : ℕ) (hb : 0 < b) :
    a ^ 2 ≠ 5 * b ^ 2 :=
  pow_ne_prime_mul_pow 5 (by decide) 2 (le_refl _) b hb a

/-- **⁴√2 is irrational (integer core).** -/
theorem fourth_root_two_integer_core (a b : ℕ) (hb : 0 < b) :
    a ^ 4 ≠ 2 * b ^ 4 :=
  pow_ne_prime_mul_pow 2 Nat.prime_two 4 (by omega) b hb a

/-- **⁵√2 is irrational (integer core).** -/
theorem fifth_root_two_integer_core (a b : ℕ) (hb : 0 < b) :
    a ^ 5 ≠ 2 * b ^ 5 :=
  pow_ne_prime_mul_pow 2 Nat.prime_two 5 (by omega) b hb a

/-! ## §B. log_p(q) is Irrational for Distinct Primes (Integer Core)

For distinct primes `p ≠ q` and any `a ≥ 1`, `b`, we have
`p^a ≠ q^b`. The proof does not use parity / descent, only the
divisibility transfer for primes:

  * `p ∣ p^a` (since `a ≥ 1`),
  * if `p^a = q^b`, then `p ∣ q^b`, so `p ∣ q`,
  * but `q` prime and `p ≠ 1` force `p = q`, contradicting `p ≠ q`.

This generalizes `two_pow_ne_three_pow` (the log₂(3) argument)
from the specific primes `(2, 3)` to all distinct prime pairs.
-/

/-- **Distinct-primes integer core.**
    For any two distinct primes `p, q` and any `a ≥ 1`, `b`,
    `p^a ≠ q^b`. -/
theorem prime_pow_ne_prime_pow
    (p q : ℕ) (hp : p.Prime) (hq : q.Prime) (hpq : p ≠ q)
    (a : ℕ) (ha : 1 ≤ a) (b : ℕ) : p ^ a ≠ q ^ b := by
  intro h
  have h_p_dvd_pa : p ∣ p ^ a := dvd_pow_self p (by omega)
  rw [h] at h_p_dvd_pa
  have h_p_dvd_q : p ∣ q := hp.dvd_of_dvd_pow h_p_dvd_pa
  rcases hq.eq_one_or_self_of_dvd p h_p_dvd_q with h1 | h1
  · exact hp.ne_one h1
  · exact hpq h1

/-! ## §B.1 Concrete specializations. -/

/-- **log₂(3) is irrational (integer core).** Recovers the result
    of `LogTwoThreeIrrational`. -/
theorem two_pow_ne_three_pow' (a : ℕ) (ha : 1 ≤ a) (b : ℕ) :
    (2 : ℕ) ^ a ≠ 3 ^ b :=
  prime_pow_ne_prime_pow 2 3 Nat.prime_two Nat.prime_three (by decide) a ha b

/-- **log₂(5) is irrational (integer core).** -/
theorem two_pow_ne_five_pow (a : ℕ) (ha : 1 ≤ a) (b : ℕ) :
    (2 : ℕ) ^ a ≠ 5 ^ b :=
  prime_pow_ne_prime_pow 2 5 Nat.prime_two (by decide) (by decide) a ha b

/-- **log₃(5) is irrational (integer core).** -/
theorem three_pow_ne_five_pow (a : ℕ) (ha : 1 ≤ a) (b : ℕ) :
    (3 : ℕ) ^ a ≠ 5 ^ b :=
  prime_pow_ne_prime_pow 3 5 Nat.prime_three (by decide) (by decide) a ha b

/-- **log₅(7) is irrational (integer core).** -/
theorem five_pow_ne_seven_pow (a : ℕ) (ha : 1 ≤ a) (b : ℕ) :
    (5 : ℕ) ^ a ≠ 7 ^ b :=
  prime_pow_ne_prime_pow 5 7 (by decide) (by decide) (by decide) a ha b

/-! ## §C. Unified Statement — The Integer-Core Irrationality Corpus

A single theorem that packages §A and §B: for every prime `p`
and every prime `q ≠ p`, and for every exponent `n ≥ 2`,

  * `ⁿ√p ∉ ℚ`   (integer core: `a^n ≠ p · b^n` for `b > 0`)
  * `log_p(q) ∉ ℚ`   (integer core: `p^a ≠ q^b` for `a ≥ 1`)

both hold simultaneously.
-/

/-- **Unified integer-core certificate.** -/
theorem irrationality_frontier_certificate
    (p q : ℕ) (hp : p.Prime) (hq : q.Prime) (hpq : p ≠ q)
    (n : ℕ) (hn : 2 ≤ n) :
    (∀ a b : ℕ, 0 < b → a ^ n ≠ p * b ^ n) ∧
    (∀ a : ℕ, 1 ≤ a → ∀ b : ℕ, p ^ a ≠ q ^ b) := by
  refine ⟨?_, ?_⟩
  · intro a b hb; exact pow_ne_prime_mul_pow p hp n hn b hb a
  · intro a ha b; exact prime_pow_ne_prime_pow p q hp hq hpq a ha b

end PrimeRootIrrationality

end Pandrosion
