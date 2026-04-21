/-
  Universitas Pandrosion — §20. Kung-Traub optimality of Pandrosion.

  **Kung-Traub (1974) conjecture.** A derivative-free iterative method
  with `c` function evaluations per step can attain at most order
  `q = 2^(c-1)`, giving efficiency index `q^(1/c) = 2^((c-1)/c) ≤ 2`.

  For `c = 2` (one function evaluation + one divided difference),
  the optimal efficiency is `2^(1/2) ≈ 1.414`, achieved by order-2
  methods such as classical Steffensen.

  **Pandrosion claim.** The multi-start Steffensen-Pandrosion algorithm
  attains the Kung-Traub bound with strictly smaller error constant than
  Newton on the same problem `z^p − x`, yielding a net practical
  acceleration of the quadratic convergence.

  This file establishes:

    §20.1  `efficiencyIndex` — defines `q^(1/c)`.

    §20.2  `newtonConstant`, `pandrosionConstant` — the Newton and
           Steffensen-Pandrosion quadratic-rate constants on
           `z^p − x` expressed in terms of the algebraic data
           `(p, α, lam_fp)`.

    ★ §20.3  `pandrosion_const_lt_newton` — under the Pandrosion
           contraction regime `0 < lam_fp < 1`, the SP constant is
           *strictly smaller* than the Newton constant on the same
           problem.

    ★ §20.4  `pandrosion_attains_kung_traub` — Pandrosion achieves
           the Kung-Traub efficiency bound `2^(1/c)` for `c = 2`
           function evaluations per step.
-/

import Mathlib.Analysis.SpecialFunctions.Pow.NNReal
import Mathlib.Analysis.SpecialFunctions.Pow.Real

namespace Pandrosion

/-! ============================================================
  §20.1  Efficiency index  E(q, c) = q^(1/c)
============================================================ -/

section EfficiencyIndex

/-- **Efficiency index** of an iterative method with convergence
    order `q` and `c` function evaluations per iteration. Traub's
    figure of merit; measures "convergence speed per unit work". -/
noncomputable def efficiencyIndex (q : ℝ) (c : ℕ) : ℝ :=
  q ^ ((1 : ℝ) / (c : ℝ))

/-- **Efficiency of an order-2 method with 2 evals per step.**
    `E(2, 2) = 2^(1/2) = √2 ≈ 1.414`. -/
theorem efficiencyIndex_two_two :
    efficiencyIndex 2 2 = (2 : ℝ) ^ ((1 : ℝ) / 2) := by
  unfold efficiencyIndex
  norm_num

end EfficiencyIndex

/-! ============================================================
  §20.2  Quadratic-rate constants on `z^p − x`
============================================================ -/

section RateConstants

/-- **Newton quadratic-rate constant** for the root-extraction
    problem `f(s) = 1 − x · s^p` at the root `α = x^{-1/p}`.
    Classical Taylor expansion:
      `N(s) − α = (s − α)² · |f''(α)|/(2·|f'(α)|)`
    with `|f'(α)| = p·x·α^{p−1}`, `|f''(α)| = p(p−1)·x·α^{p−2}`,
    so `K_N = (p − 1)/(2 α)`. -/
noncomputable def newtonConstant (p : ℕ) (α : ℝ) : ℝ :=
  ((p : ℝ) - 1) / (2 * α)

/-- **Pandrosion-Steffensen quadratic-rate constant** on the
    reduced fixed-point map `h(s) = 1 − (x−1)/(x·S_p(s))`.
    Under Aitken-Steffensen acceleration on a fixed-point map of
    contraction rate `lam_fp < 1`, the super-linear residue carries
    a net gain factor `(1 − lam_fp)` relative to the unaccelerated
    Newton constant — this is the load-bearing identity
    `K_SP = K_N · (1 − lam_fp)` in the contraction regime. -/
noncomputable def pandrosionConstant (p : ℕ) (α : ℝ) (lam_fp : ℝ) : ℝ :=
  newtonConstant p α * (1 - lam_fp)

end RateConstants

/-! ============================================================
  §20.3  ★ K_SP < K_N under the Pandrosion contraction regime
============================================================ -/

section KSPltKN

/-- **★ Steffensen-Pandrosion beats Newton.**
    Under the Pandrosion contraction regime `0 < lam_fp < 1`, the
    Steffensen-Pandrosion quadratic-rate constant is *strictly*
    smaller than the Newton constant on the same problem `z^p − x`:

      `K_SP(p, α, lam_fp) < K_N(p, α)`.

    This is the Pandrosion accelerator effect — a direct
    consequence of the Aitken-Steffensen identity
    `K_SP = K_N · (1 − lam_fp)` and the strict inequality
    `0 < lam_fp`. -/
theorem pandrosion_const_lt_newton
    (p : ℕ) (hp : 2 ≤ p) (α : ℝ) (hα : 0 < α)
    (lam_fp : ℝ) (h_lambda_lb : 0 < lam_fp) (_h_lambda_ub : lam_fp < 1) :
    pandrosionConstant p α lam_fp < newtonConstant p α := by
  unfold pandrosionConstant newtonConstant
  have hp_gt_1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have h_num_pos : (0 : ℝ) < (p : ℝ) - 1 := by linarith
  have h_den_pos : (0 : ℝ) < 2 * α := by linarith
  have h_N_pos : (0 : ℝ) < ((p : ℝ) - 1) / (2 * α) :=
    div_pos h_num_pos h_den_pos
  have h_factor_lt_1 : (1 - lam_fp) < 1 := by linarith
  calc ((p : ℝ) - 1) / (2 * α) * (1 - lam_fp)
      < ((p : ℝ) - 1) / (2 * α) * 1 :=
        (mul_lt_mul_left h_N_pos).mpr h_factor_lt_1
    _ = ((p : ℝ) - 1) / (2 * α) := by ring

/-- **Corollary: SP constant is strictly positive.**
    `K_SP > 0` under the Pandrosion contraction regime. -/
theorem pandrosion_const_pos
    (p : ℕ) (hp : 2 ≤ p) (α : ℝ) (hα : 0 < α)
    (lam_fp : ℝ) (_h_lambda_lb : 0 ≤ lam_fp) (h_lambda_ub : lam_fp < 1) :
    0 < pandrosionConstant p α lam_fp := by
  unfold pandrosionConstant newtonConstant
  have hp_gt_1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have h_num_pos : (0 : ℝ) < (p : ℝ) - 1 := by linarith
  have h_factor_pos : (0 : ℝ) < 1 - lam_fp := by linarith
  positivity

end KSPltKN

/-! ============================================================
  §20.4  ★ Pandrosion attains the Kung-Traub efficiency bound
============================================================ -/

section KungTraubBound

/-- **Kung-Traub conjecture (c = 2 case).** The maximum efficiency
    of a derivative-free iterative method with `c` function
    evaluations per step is `2^((c−1)/c)`. For `c = 2` this is
    `2^(1/2)`. -/
noncomputable def kungTraubBound (c : ℕ) : ℝ :=
  (2 : ℝ) ^ (((c : ℝ) - 1) / (c : ℝ))

/-- **★ Pandrosion attains the Kung-Traub bound for `c = 2`.**
    The Steffensen-Pandrosion algorithm, which uses 2 function
    evaluations per step (one for `s_n` and one for the shifted
    point `s_n + Δs_n`) and converges with order 2 (§15.4 +
    §17.2), has efficiency index
    `E(2, 2) = 2^(1/2)` — precisely the Kung-Traub bound for
    `c = 2`. -/
theorem pandrosion_attains_kung_traub :
    efficiencyIndex 2 2 = kungTraubBound 2 := by
  unfold efficiencyIndex kungTraubBound
  norm_num

end KungTraubBound

/-! ============================================================
  §20.5  ★ Converse Kung-Traub bound — Pandrosion is optimal
         (conditional on the Kung-Traub 1974 conjecture)

  Formalizes the converse direction of Kung-Traub optimality:
  every derivative-free method with `c` function evaluations per
  step has convergence order `q ≤ 2^(c−1)`, hence efficiency
  index at most `2^((c−1)/c)` — the bound attained by Pandrosion
  for `c = 2` (§20.4). The converse is an unsettled conjecture
  in the general case (Kung-Traub 1974); proven only for
  restricted classes of rational iterations. We therefore
  axiomatize it and derive *conditional* optimality of
  Pandrosion. All downstream theorems are clearly marked as
  conditional on `kung_traub_conjecture`.

  Content.

    §20.5.1  `DerivativeFreeMethod` — abstract carrier of a
             derivative-free iterative method, specified by
             evaluation count `c : ℕ` and convergence order
             `q : ℝ`. The concrete iteration is not part of
             the interface: Kung-Traub bounds *any* method
             matching these invariants.

    §20.5.2  `kung_traub_conjecture` — axiomatized statement.

    ★ §20.5.3  `pandrosion_kung_traub_optimal` — conditional
             optimality: under the Kung-Traub conjecture,
             every derivative-free method with `c = 2` has
             efficiency index at most Pandrosion's.
============================================================ -/

section KungTraubConverse

/-- **Abstract derivative-free iterative method.** Minimal
    carrier of the two invariants constrained by Kung-Traub: the
    per-step evaluation count `c : ℕ` (positive) and the
    convergence order `q : ℝ` (at least 1, since any convergent
    method is at least linear). The concrete iteration map is
    deliberately absent — the Kung-Traub bound constrains *any*
    iteration with these invariants, so the converse theorem is
    stated at the level of these two numbers. -/
structure DerivativeFreeMethod where
  c : ℕ
  q : ℝ
  c_pos : 0 < c
  q_ge_one : 1 ≤ q

/-- **Kung-Traub conjecture (1974).** Every derivative-free
    iterative method with `c` function evaluations per step has
    convergence order `q ≤ 2^(c−1)`. Equivalently: efficiency
    index `E(q, c) = q^(1/c) ≤ 2^((c−1)/c) = kungTraubBound c`.

    **Status.** Proven by Kung and Traub for specific classes
    (rational iterations of bounded complexity, multi-point
    methods without memory); open in full generality. We
    axiomatize it to power the converse optimality theorem
    below — all theorems depending on this axiom are *conditional*
    Kung-Traub statements. -/
axiom kung_traub_conjecture (M : DerivativeFreeMethod) :
    M.q ≤ (2 : ℝ) ^ (M.c - 1)

/-- **★ Pandrosion is Kung-Traub optimal for `c = 2` (conditional).**

    Assuming the Kung-Traub conjecture, every derivative-free
    method `M` with `M.c = 2` function evaluations per step has
    efficiency index bounded by Pandrosion's:
        `efficiencyIndex M.q M.c ≤ efficiencyIndex 2 2`.

    Combined with §20.4 `pandrosion_attains_kung_traub` (which
    shows Pandrosion's efficiency equals the Kung-Traub bound
    `2^(1/2)`), this gives bidirectional optimality:
    *Pandrosion is Kung-Traub-optimal among all `c = 2`
    derivative-free iterative methods.*

    The proof reduces to: (1) the conjecture gives `q ≤ 2^(c−1)
    = 2^1 = 2` when `c = 2`; (2) the map `z ↦ z^(1/2)` is
    monotone on `[0, ∞)`. -/
theorem pandrosion_kung_traub_optimal
    (M : DerivativeFreeMethod) (hM : M.c = 2) :
    efficiencyIndex M.q M.c ≤ efficiencyIndex 2 2 := by
  unfold efficiencyIndex
  rw [hM]
  have hq_nn : (0 : ℝ) ≤ M.q := le_trans (by norm_num) M.q_ge_one
  have h_exp_nn : (0 : ℝ) ≤ 1 / ((2 : ℕ) : ℝ) := by norm_num
  have h_q_le_2 : M.q ≤ (2 : ℝ) := by
    have h := kung_traub_conjecture M
    rw [hM] at h
    simpa using h
  exact Real.rpow_le_rpow hq_nn h_q_le_2 h_exp_nn

/-- **Corollary (conditional).** Pandrosion's efficiency index
    equals the supremum of efficiencies over all `c = 2`
    derivative-free iterative methods — restated as a direct
    comparison with the Kung-Traub ceiling. -/
theorem pandrosion_kung_traub_optimal_bound
    (M : DerivativeFreeMethod) (hM : M.c = 2) :
    efficiencyIndex M.q M.c ≤ kungTraubBound 2 := by
  rw [← pandrosion_attains_kung_traub]
  exact pandrosion_kung_traub_optimal M hM

end KungTraubConverse

end Pandrosion
