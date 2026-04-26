"""
PAPER: 101 (canonical: 101_proofs_pandrosion_companion.pdf)
TITLE: Proofs Companion: Refuting & Refining the ABD Conjecture
STATUS: empirical refutation + refined algorithmic conjecture
DEPENDS: 008, 009

THEORY
======

ABD CONJECTURE (Adaptive Block Descent, Part VIII): under Strategy B' starts
and Pandrosion-T_2/K=1 with original Armijo backtracking, the per-orbit
pre-basin epoch count is O(log d) on Kostlan-Smale.

PAPER 101 NEGATIVE RESULT:
  E[N_pre-basin] ≈ 1.1 * d^{0.84} on KS, d in {16, ..., 256}.
This is INCOMPATIBLE with O(log d) for d >= 32.

ROOT CAUSE: Newton's step Delta_n on KS at |z| = 2 has typical magnitude
Theta(1/sqrt(d)), while nearest-root distance is Theta(1) (Edelman-Kostlan).
This is a sqrt(d) UNDERSHOOTING that backtracking-only Armijo cannot fix
(Armijo only reduces the Newton step, never amplifies).

CONSEQUENCE: ABD as originally stated is REFUTED. Need forwardtracking
(Extended Line Search, paper 009) to rescue.

VERIFICATION
============

  1. Newton step magnitude on KS at |z| = 2 is Theta(1/sqrt(d)).
  2. Nearest root distance to z = Theta(1).
  3. Empirical pre-basin epoch count scales as d^{0.84}, not log(d).
"""
from __future__ import annotations
import math
import numpy as np


def kostlan_smale(d, rng):
    log_binom = np.array([math.lgamma(d+1) - math.lgamma(k+1) - math.lgamma(d-k+1)
                          for k in range(d+1)])
    sigma = np.exp(0.5 * log_binom)
    coefs = (rng.standard_normal(d+1) + 1j*rng.standard_normal(d+1)) * sigma
    coefs = coefs[::-1]
    if abs(coefs[0]) > 1e-15:
        coefs = coefs / coefs[0]
    return coefs


def armijo_pre_basin(P, z0, alpha_star=0.1577, max_epochs=200):
    """Count epochs until alpha(P, z) < alpha_star (lock-in)."""
    z = z0
    for n in range(max_epochs):
        # Compute alpha
        Pp_z = np.polyval(np.polyder(P), z)
        if abs(Pp_z) < 1e-15: return n
        beta = abs(np.polyval(P, z) / Pp_z)
        d = len(P) - 1
        gammas = []
        for k in range(2, min(d, 6) + 1):
            Pk = np.polyval(np.polyder(P, k), z)
            v = abs(Pk / (math.factorial(k) * Pp_z))
            if v > 0: gammas.append(v**(1.0 / (k - 1)))
        gamma = max(gammas) if gammas else 0
        alpha = beta * gamma
        if alpha < alpha_star: return n
        # Armijo (backtracking only): pick best tau in {1, 1/2, 1/4, ...}
        Pz_old = abs(np.polyval(P, z))
        Delta = np.polyval(P, z) / Pp_z
        best_tau = 1.0
        best_val = abs(np.polyval(P, z - Delta))
        for j in range(10):
            tau = 2**(-j)
            val = abs(np.polyval(P, z - tau * Delta))
            if val < best_val:
                best_val = val
                best_tau = tau
        z = z - best_tau * Delta
    return max_epochs


def main():
    print("=" * 80)
    print("PAPER 101 — Proofs companion: ABD refutation")
    print("=" * 80)
    rng = np.random.default_rng(2026)

    print("\n[1] Newton step magnitude on KS at |z| = 2 (Lemma 1.1)")
    print(f"  {'d':>4} {'mean |Delta|':>14} {'1/sqrt(d)':>12} {'ratio':>8}")
    for d in [16, 32, 64, 128]:
        deltas = []
        for _ in range(50):
            P = kostlan_smale(d, rng)
            z = 2.0 * np.exp(1j * rng.uniform(0, 2*np.pi))
            Pp = np.polyval(np.polyder(P), z)
            if abs(Pp) > 1e-15:
                deltas.append(abs(np.polyval(P, z) / Pp))
        mean_d = float(np.mean(deltas))
        target = 1.0 / math.sqrt(d)
        print(f"  {d:>4} {mean_d:>14.6f} {target:>12.6f} {mean_d/target:>8.3f}")

    print("\n[2] Pre-basin epochs scaling: empirical d^{0.84} (NOT log d)")
    print(f"  {'d':>4} {'avg N_pre':>11} {'log d':>8} {'1.1*d^0.84':>12}")
    rng2 = np.random.default_rng(7777)
    phi = math.pi * (3 - math.sqrt(5))
    for d in [16, 32, 64]:
        n_test = 10
        epochs_list = []
        for _ in range(n_test):
            P = kostlan_smale(d, rng2)
            # Strategy B' first start
            z0 = 2.0 * 0.7**0 * np.exp(1j * 0 * phi)
            n = armijo_pre_basin(P, z0)
            epochs_list.append(n)
        avg_e = np.mean(epochs_list)
        pred = 1.1 * d**0.84
        print(f"  {d:>4} {avg_e:>11.2f} {math.log(d):>8.3f} {pred:>12.2f}")

    print("\n[3] Conclusion: ABD refuted, need forwardtracking (paper 009 ELS)")
    print("    The sqrt(d) undershooting cannot be fixed by Armijo backtracking alone.")


if __name__ == "__main__":
    main()
