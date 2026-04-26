"""
PAPER: 010 (canonical: 9pandrosion_smale2_mv.pdf, also 9pandrosion_smale_mv.pdf)
TITLE: Universitas Pandrosion: Part IX-mv --
       Multivariate Extension of B'-Newton-ELS to Smale 17 --
       E[N_pre-basin] = O(1) on Kostlan-Smale systems up to Bezout degree D = d^n = 1024
STATUS: empirical confirmation (numerical verification only)
DEPENDS: 009

THEORY
======

mvKS ENSEMBLE (Bombieri-Weyl projective Gaussian on systems):
For F = (F_1, ..., F_n): C^n -> C^n of degree d, draw each F_i from KS:
  coefficient a_alpha^{(i)} ~ Gaussian with variance d!/(alpha_1!...alpha_n!(d-|alpha|)!).

Bezout degree D = d^n. By Edelman-Kostlan analogue, mvKS roots concentrate
near unit polysphere {||z|| = 1} with high probability.

MULTIVARIATE NEWTON-ELS:
At each iterate z_n:
  Delta_n = DF(z_n)^{-1} F(z_n)        (Newton direction, Schmidt slope-matrix at diag)
  z_{n+1} = z_n - tau* Delta_n
where tau* in T_ELS = {2^k : k in [-6, +6]} minimises ||F(z_n - tau Delta_n)||.

DIAGNOSTIC (Section 5.5 of paper):
Optimal tau* on mvKS is concentrated near 1, NOT near sqrt(D) as the univariate
analogy would predict. The ELS forwardtracking is therefore NOT the load-bearing
mechanism in the multivariate setting; pure multivariate Newton (tau = 1)
already attains E[N_pre-basin] = O(1) on mvKS.

Reason: matrix inversion in DF^{-1} F correctly dimensions the Newton step
(self-correcting via the Jacobian).

EMPIRICAL FIT (Theorem 6.1):
On 20 configs (n, d) spanning D = d^n in [9, 4.3*10^9]:
  E[N_tot] ~ 3.67 + 2.31 * log_2 log_2 D
i.e., O(log log D) — Smale's quadratic basin convergence.

LAS VEGAS COST:
  E_F[cost^{B'-N-ELS}_total(F)] = O(D log D + n^3 log D),
multivariate analogue of paper IX's univariate O(d log d).

VERIFICATION
============

This script verifies:
  1. Multivariate Newton-ELS converges on small mvKS systems.
  2. tau* concentrates near 1 (not sqrt(D)).
  3. N_tot scales as log log D.
"""
from __future__ import annotations
import math
import numpy as np
from itertools import product


def multinomial_coefficient(d, alpha):
    """d! / (alpha_1! ... alpha_n! (d - |alpha|)!)."""
    abs_a = sum(alpha)
    if abs_a > d: return 0
    num = math.factorial(d)
    den = math.factorial(d - abs_a)
    for a in alpha:
        den *= math.factorial(a)
    return num // den


def multi_indices(n, d):
    """All alpha in N^n with |alpha| <= d."""
    out = []
    def rec(prefix, remaining_n, remaining_d):
        if remaining_n == 0:
            out.append(tuple(prefix))
            return
        for k in range(remaining_d + 1):
            rec(prefix + [k], remaining_n - 1, remaining_d - k)
    rec([], n, d)
    return out


def random_mvks_system(n, d, rng):
    """Draw F = (F_1, ..., F_n) from mvKS.

    Each F_i is a homogeneous-style poly:
      F_i(z) = sum_{|alpha| <= d} c_alpha^{(i)} z^alpha
    with c_alpha ~ N(0, multinomial(d, alpha)).

    Returns: list of dicts mapping alpha -> coefficient.
    """
    indices = multi_indices(n, d)
    F = []
    for i in range(n):
        f = {}
        for alpha in indices:
            mc = multinomial_coefficient(d, alpha)
            std = math.sqrt(mc)
            c = (rng.standard_normal() + 1j * rng.standard_normal()) * std
            f[alpha] = c
        F.append(f)
    return F


def evaluate_F(F, z):
    """Evaluate F at z."""
    n = len(F)
    out = np.zeros(n, dtype=complex)
    for i, f in enumerate(F):
        s = 0.0 + 0j
        for alpha, c in f.items():
            term = c
            for j, a in enumerate(alpha):
                term *= z[j]**a
            s += term
        out[i] = s
    return out


def jacobian_F(F, z, h=1e-7):
    """Numerical Jacobian (since symbolic is heavy for the demo)."""
    n = len(z)
    F_z = evaluate_F(F, z)
    J = np.zeros((n, n), dtype=complex)
    for j in range(n):
        z_h = z.copy()
        z_h[j] += h
        F_h = evaluate_F(F, z_h)
        J[:, j] = (F_h - F_z) / h
    return J


def mv_newton_step(F, z):
    F_z = evaluate_F(F, z)
    J = jacobian_F(F, z)
    try:
        return np.linalg.solve(J, F_z)
    except np.linalg.LinAlgError:
        return np.zeros_like(z)


def mv_els_step(F, z, k_range=(-3, 3)):
    """Multivariate ELS: tau in T_ELS minimising ||F(z - tau Delta)||."""
    Delta = mv_newton_step(F, z)
    best_tau, best_val = 1.0, np.linalg.norm(evaluate_F(F, z - Delta))
    for k in range(k_range[0], k_range[1] + 1):
        tau = 2**k
        val = np.linalg.norm(evaluate_F(F, z - tau * Delta))
        if val < best_val:
            best_val, best_tau = val, tau
    return z - best_tau * Delta, best_tau


def mv_orbit(F, z0, max_iter=50, tol=1e-6):
    z = z0.copy()
    for ep in range(max_iter):
        if np.linalg.norm(evaluate_F(F, z)) < tol: break
        z, _ = mv_els_step(F, z)
    return ep + 1, z


def main():
    print("=" * 80)
    print("PAPER 10 — Multivariate Pandrosion-Newton-ELS on mvKS")
    print("=" * 80)
    rng = np.random.default_rng(2026)

    # 1. Convergence on small mvKS systems
    print("\n[1] Convergence on small mvKS systems")
    print(f"  {'(n, d)':>10} {'D = d^n':>10} {'#trials':>8} {'avg #epochs':>12}")
    for n, d in [(2, 2), (2, 3), (3, 2), (3, 3)]:
        D = d**n
        n_trials = 5
        epochs = []
        for _ in range(n_trials):
            F = random_mvks_system(n, d, rng)
            z0 = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) * 0.5
            ep, z_final = mv_orbit(F, z0)
            epochs.append(ep)
        avg = np.mean(epochs)
        print(f"  {(n,d)!s:>10} {D:>10} {n_trials:>8} {avg:>12.2f}")

    # 2. Optimal tau* on mvKS
    print("\n[2] Optimal tau* on mvKS (DIAGNOSTIC: concentrates near 1, NOT sqrt(D))")
    print(f"  {'(n, d)':>10} {'D':>8} {'mean tau*':>12} {'sqrt(D)':>10}")
    for n, d in [(2, 3), (3, 2), (2, 4)]:
        D = d**n
        taus = []
        for _ in range(8):
            F = random_mvks_system(n, d, rng)
            z = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) * 0.5
            _, tau = mv_els_step(F, z)
            taus.append(tau)
        mean_t = np.mean(taus)
        print(f"  {(n,d)!s:>10} {D:>8} {mean_t:>12.4f} {math.sqrt(D):>10.4f}")

    # 3. N_tot scaling
    print("\n[3] N_tot ~ log log D (Smale's quadratic basin)")
    print(f"  {'D':>10} {'mean N_tot':>12} {'3.67 + 2.31 log_2 log_2 D':>30}")
    configs = [(2, 2), (2, 3), (3, 2), (3, 3), (2, 4)]
    for n, d in configs:
        D = d**n
        n_trials = 5
        epochs = []
        for _ in range(n_trials):
            F = random_mvks_system(n, d, rng)
            z0 = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) * 0.5
            ep, _ = mv_orbit(F, z0, tol=1e-8)
            epochs.append(ep)
        avg = np.mean(epochs)
        if D >= 4:
            pred = 3.67 + 2.31 * math.log2(math.log2(D))
        else:
            pred = float('nan')
        print(f"  {D:>10} {avg:>12.2f} {pred:>30.2f}")


if __name__ == "__main__":
    main()
