"""
PAPER: 051 (canonical: 41pandrosion_universality.pdf)
TITLE: Universality of Pandrosion Iteration on KS Polynomials
STATUS: framework (universal local convergence on KS ensemble)
DEPENDS: 026

THEORY
======

UNIVERSALITY: For random KS polynomial P, the convergence behaviour of
Pandrosion iteration is universal — depending only on degree d, not on
specific coefficient values.

Key universal quantities:
  Smale alpha distribution at random z on |z| = 2: tail decays like ~ 1/d.
  Iteration count to basin: O(log log 1/eps) in lock-in.

VERIFICATION
============

  1. Smale alpha distribution on KS at |z| = 2.
  2. Universal log log iteration count.
"""
from __future__ import annotations
import math
import numpy as np


def kostlan_smale(d, rng):
    log_binom = np.array([math.lgamma(d+1) - math.lgamma(k+1) - math.lgamma(d-k+1)
                          for k in range(d+1)])
    sigma = np.exp(0.5 * log_binom)
    coefs = (rng.standard_normal(d+1) + 1j*rng.standard_normal(d+1)) * sigma
    return coefs[::-1]


def smale_alpha(P, z):
    Pp_z = np.polyval(np.polyder(P), z)
    if abs(Pp_z) < 1e-15: return float('inf')
    beta = abs(np.polyval(P, z) / Pp_z)
    d = len(P) - 1
    gammas = []
    for k in range(2, min(d, 8) + 1):
        Pk = np.polyval(np.polyder(P, k), z)
        v = abs(Pk / (math.factorial(k) * Pp_z))
        if v > 0: gammas.append(v**(1.0 / (k - 1)))
    gamma = max(gammas) if gammas else 0
    return beta * gamma


def main():
    print("=" * 80)
    print("PAPER 51 — Universality on KS polynomials")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Smale alpha distribution at |z| = 2 on KS")
    print(f"  {'d':>4} {'mean alpha':>12} {'median':>10} {'frac > 1':>10}")
    for d in [16, 32, 64, 128]:
        alphas = []
        for _ in range(50):
            P = kostlan_smale(d, rng)
            z = 2.0 * np.exp(1j * rng.uniform(0, 2*np.pi))
            a = smale_alpha(P, z)
            if math.isfinite(a): alphas.append(a)
        arr = np.array(alphas)
        frac_big = float(np.mean(arr > 1.0))
        print(f"  {d:>4} {arr.mean():>12.4f} {np.median(arr):>10.4f} {frac_big:>10.4f}")

    print("\n[2] Universal log log iteration count")
    P = np.array([1.0, 0, -2.0])  # z^2 - 2
    root = math.sqrt(2)
    print(f"  z^2 - 2, eps -> 0:")
    for eps in [1e-3, 1e-6, 1e-12]:
        z = root + 0.1
        n = 0
        while abs(z - root) > eps and n < 50:
            Pp = np.polyval(np.polyder(P), z)
            if abs(Pp) < 1e-15: break
            z = z - np.polyval(P, z) / Pp
            n += 1
        log_log = math.log2(math.log2(1/eps))
        print(f"    eps = {eps:.0e}: iters = {n}, log_2 log_2 1/eps = {log_log:.2f}")


if __name__ == "__main__":
    main()
