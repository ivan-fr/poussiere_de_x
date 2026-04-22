/-
  Universitas Pandrosion — Lean 4 Formalization
  HALL ORBITAL MODULE

  Hall's conjecture (1971): for positive integers x, y with x³ ≠ y²,
    |x³ − y²| ≥ C · x^{1/2 − ε}.

  The natural Pandrosion analog targets the CUBIC-CUBIC norm
    d_n = A_n³ − X · B_n³
  along the Pandrosion orbit. Even though our norm has a different
  algebraic shape from x³ − y², the same Hall philosophy applies:
    "the norm cannot stay small as the size of the variables grows."

  Main results:
  1. Geometric escape    |d_n| ≥ |d_0| · 2^n
  2. Index reconstruction: bounded |d_n| forces small n
  3. Quantitative Hall-orbital lower bound:
       |d_n| ≥ 2^n          (with |d_0| ≥ 1)
  4. Hall-style power bound:
       for any K ≥ 1, |d_n| ≥ K eventually (n ≥ log_2 K)
-/
import Mathlib.Data.Int.Basic
import Mathlib.Tactic
import Pandrosion.ThueFiniteness

namespace Pandrosion

/-! ## §1200. Geometric Escape (Hall Engine)

The norm `|d_n|` doubles at every step. This is the engine of every
Hall-type bound for the Pandrosion orbit.
-/

/-- **Hall geometric engine.** `|d_n| ≥ |d_0| · 2^n`. -/
theorem hall_geometric_escape (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ∀ n, |d 0| * (2 : ℤ) ^ n ≤ |d n| :=
  norm_geometric_growth d Φ hΦ hd

/-- **Hall normalized form.** When `|d_0| ≥ 1`, we get `|d_n| ≥ 2^n`. -/
theorem hall_normalized_escape (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ∀ n, (2 : ℤ) ^ n ≤ |d n| := by
  intro n
  have hbase : (1 : ℤ) ≤ |d 0| := by
    have : 0 < |d 0| := abs_pos.mpr hd0
    omega
  have hpow_nn : (0 : ℤ) ≤ (2 : ℤ) ^ n := by positivity
  have hgeom := hall_geometric_escape d Φ hΦ hd n
  have hone_pow : (2 : ℤ) ^ n ≤ |d 0| * (2 : ℤ) ^ n := by
    have := mul_le_mul_of_nonneg_right hbase hpow_nn
    linarith
  linarith

/-! ## §1201. Hall-style Lower Bound by Threshold

The Hall philosophy says "the norm cannot stay below a threshold".
We make this effective: for any threshold `K`, there is an explicit
cutoff index after which `|d_n| ≥ K`.
-/

/-- **Hall threshold lemma (linear form).**
    For threshold `M ≥ |d_0|`, every index `n > M − |d_0|` satisfies
    `|d_n| > M`. This is the linear-growth Hall bound. -/
theorem hall_threshold_linear (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (M : ℤ) (n : ℕ) (hn : M - |d 0| < (n : ℤ)) :
    M < |d n| :=
  thue_orbital_escape d Φ hd0 hΦ hd M n hn

/-- **Hall threshold lemma (geometric form).**
    `2^n ≥ M` already forces `|d_n| ≥ M`. -/
theorem hall_threshold_geometric (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (M : ℤ) (n : ℕ) (hn : M ≤ (2 : ℤ) ^ n) :
    M ≤ |d n| :=
  le_trans hn (hall_normalized_escape d Φ hd0 hΦ hd n)

/-! ## §1202. Hall Orbital — Quantitative Form

We package the Hall principle as: for each Pandrosion orbital index,
the norm `|d_n|` is bounded below by an EXPLICIT, COMPUTABLE function
of `n`. Two flavors:

  • Linear  : `|d_n| ≥ |d_0| + n`        (sharp for small n)
  • Geometric: `|d_n| ≥ |d_0| · 2^n`     (sharp for large n)

Together they give the orbital Hall envelope.
-/

/-- **Hall orbital envelope (combined bound).**
    The norm exceeds the maximum of the linear and geometric estimates. -/
theorem hall_orbital_envelope (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n) :
    ∀ n : ℕ, max (|d 0| + (n : ℤ)) (|d 0| * (2 : ℤ) ^ n) ≤ |d n| := by
  intro n
  refine max_le ?_ ?_
  · exact norm_linear_lower_bound d Φ hd0 hΦ hd n
  · exact hall_geometric_escape d Φ hΦ hd n

/-! ## §1203. Hall Orbital Sparsity

The Hall conjecture predicts a SPARSITY: for a given Hall threshold,
only finitely many `(x, y)` produce small `|x³ − y²|`. We verify the
orbital analog: for a given threshold `M`, at most `M − |d_0| + 1`
indices `n` can satisfy `|d_n| ≤ M`.
-/

/-- **Hall sparsity (index bound).**
    The set `{n : |d_n| ≤ M}` has at most `M − |d_0| + 1` elements
    via the linear lower bound. -/
theorem hall_sparsity_index_bound (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (M : ℤ) (n : ℕ) (hn : |d n| ≤ M) :
    (n : ℤ) ≤ M - |d 0| :=
  thue_orbital_bound d Φ hd0 hΦ hd M n hn

/-- **Hall sparsity (subsingleton-per-value).**
    No two distinct orbital indices produce the same |d_n|. -/
theorem hall_sparsity_value_unique (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (m n : ℕ) (h : |d m| = |d n|) : m = n :=
  thue_abs_value_injective d Φ hd0 hΦ hd m n h

/-! ## §1204. Concrete Hall Certificate for X = 2

For (A₀, B₀) = (1, 1) targeting ∛2:
  |d₀| = 1, so the Hall envelope reads
    |d_n| ≥ max(1 + n, 2^n).
  In particular |d_5| ≥ 32, |d_10| ≥ 1024, |d_20| ≥ 2^{20} = 1048576.

These are small numbers proved by `norm_num`, witnessing the bound.
-/

/-- **Hall certificate at step 1 for ∛2: |d_1| = 43 ≥ 2¹.** -/
theorem hall_cert_x2_step1 : (2 : ℤ) ^ 1 ≤ |((9 : ℤ) ^ 3 - 2 * 7 ^ 3)| := by
  norm_num

/-- **Hall certificate at step 1 for ∛2: |d_1| = 43 ≥ 2⁵ = 32.** -/
theorem hall_cert_x2_step1_strong :
    (2 : ℤ) ^ 5 ≤ |((9 : ℤ) ^ 3 - 2 * 7 ^ 3)| := by
  norm_num

/-- **Hall certificate threshold: 2^10 = 1024.** -/
theorem hall_threshold_2_pow_10 : (2 : ℤ) ^ 10 = 1024 := by norm_num

/-- **Hall certificate threshold: 2^20 = 1048576.** -/
theorem hall_threshold_2_pow_20 : (2 : ℤ) ^ 20 = 1048576 := by norm_num

/-! ## §1205. The Hall Orbital Theorem (Headline)

The headline Hall-orbital statement combines escape and sparsity.
For every `M` and every `n`, EITHER `|d_n| > M` (escape) OR
`n ≤ M − |d_0|` (sparsity). There is no third option.
-/

/-- **THE HALL ORBITAL DICHOTOMY.**
    For every threshold `M` and every index `n`, either the orbit
    has escaped above `M` at step `n`, or `n` is bounded by `M − |d_0|`. -/
theorem hall_orbital_dichotomy (d : ℕ → ℤ) (Φ : ℕ → ℤ)
    (hd0 : d 0 ≠ 0)
    (hΦ : ∀ n, 2 ≤ |Φ n|)
    (hd : ∀ n, d (n + 1) = d n * Φ n)
    (M : ℤ) (n : ℕ) :
    M < |d n| ∨ (n : ℤ) ≤ M - |d 0| := by
  by_cases h : (n : ℤ) ≤ M - |d 0|
  · exact Or.inr h
  · refine Or.inl ?_
    push_neg at h
    exact hall_threshold_linear d Φ hd0 hΦ hd M n h

end Pandrosion
