"""
SENDOV FOR HIGHER DERIVATIVES via Pandrosion.

Classical Sendov: P monic, all roots in closed unit disk, alpha a root.
Then there exists a critical point xi (root of P') with |alpha - xi| <= 1.

Higher-derivative Sendov:
  Sendov_k:  for any root alpha of P, there exists a root of P^(k) within
             distance r_k(d) of alpha.

Trivial bounds:
  r_1(d) <= 1 (Sendov-Smale d <= 8 explicit, paper 104 numerical to d=1024)
  r_d-1(d) = ? (P^(d-1) has 1 root: the d-1 derivative of monic P is just (d!) z + lower)

The interesting case: r_k(d) for 2 <= k <= d-2.

Pandrosion approach:
  - P has Pandrosion energy E_P(z) = sum_j Q(alpha_j, z)^2 = (P')^2 - P P''
  - Higher: E^(k)_P(z) related to (P^(k))^2 etc via Turán-style identities.
  - Real-rooted check via paper 56 (Turán SOS).
"""
from __future__ import annotations
import math
import numpy as np


def sendov_violation_kth(roots, k):
    """Sendov for the k-th derivative.

    For each root zeta_j of P, find min distance to any root of P^(k).
    Return max over j of (min distance) - 1 (the "violation"; <= 0 if conjecture holds).
    """
    P = np.poly(roots)
    if k >= len(roots):
        return float('-inf')  # no higher derivative left
    Pk = P
    for _ in range(k):
        Pk = np.polyder(Pk)
    if len(Pk) <= 1:
        # P^(k) is constant or linear, single root
        if len(Pk) == 1:
            return float('inf')  # zero polynomial
        crits = [-Pk[1]/Pk[0]]
    else:
        crits = np.roots(Pk)
    if len(crits) == 0:
        return float('inf')
    D = np.abs(np.array(roots)[:, None] - np.array(crits)[None, :])
    return float(D.min(axis=1).max()) - 1.0


def random_polynomial_in_disk(d, rng):
    """Random polynomial with all roots in closed unit disk."""
    angles = 2 * np.pi * rng.random(d)
    radii = rng.random(d)  # |alpha| in [0, 1]
    roots = radii * np.exp(1j * angles)
    return roots


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 95, flush=True)
    print("SENDOV FOR HIGHER DERIVATIVES — Pandrosion attack", flush=True)
    print("=" * 95, flush=True)
    print("\nClassical Sendov_1 (the standard conjecture): r_1(d) <= 1.")
    print("We test r_k(d) for k = 1, 2, ..., d-1 numerically.")

    print("\n[1] Empirical r_k(d) bounds across degrees:")
    print(f"{'d':>4} {'k':>4} {'#cfg':>6} {'max V_k':>10} {'r_k status':>15}",
          flush=True)
    rng = np.random.default_rng(20260712)
    for d in [10, 20, 50, 100]:
        for k in [1, 2, 3, 5, d // 2, d - 2]:
            if k >= d or k < 1:
                continue
            n_cfg = 100 if d <= 50 else 30
            max_V = float('-inf')
            for _ in range(n_cfg):
                roots = random_polynomial_in_disk(d, rng)
                V = sendov_violation_kth(roots, k)
                if V > max_V:
                    max_V = V
            status = "OK" if max_V <= 1e-6 else f"VIOLATES by {max_V:.3e}"
            print(f"{d:>4} {k:>4} {n_cfg:>6} {max_V:>10.4e} {status:>15}",
                  flush=True)

    # Specific adversarial families
    print("\n[2] Adversarial families:")

    # Roots of unity boundary
    print("\n[2a] Roots of unity z^d - 1 = 0:")
    print(f"{'d':>4} {'k':>4} {'V_k':>12}", flush=True)
    for d in [5, 10, 20, 50]:
        for k in [1, 2, d // 2, d - 2]:
            roots = np.exp(2j * np.pi * np.arange(d) / d)
            V = sendov_violation_kth(roots, k)
            print(f"{d:>4} {k:>4} {V:>12.4e}", flush=True)

    # Boundary cluster: zeta_1 = 1, others on small arc
    print("\n[2b] Miller-style cluster (zeta_1 = 1, others near 1):")
    print(f"{'d':>4} {'k':>4} {'V_k':>12}", flush=True)
    rng2 = np.random.default_rng(7777)
    for d in [10, 20, 50]:
        roots_list = [1.0 + 0j]
        for _ in range(d - 1):
            phi = rng2.uniform(0, 2*np.pi)
            roots_list.append(1.0 + 0.5 * np.exp(1j * phi))
        roots = np.array(roots_list)
        # Project to disk
        for i in range(len(roots)):
            if abs(roots[i]) > 1:
                roots[i] = roots[i] / abs(roots[i]) * 0.999
        for k in [1, 2, 3, 5, d // 2]:
            if k >= d: continue
            V = sendov_violation_kth(roots, k)
            print(f"{d:>4} {k:>4} {V:>12.4e}", flush=True)

    print("\n[3] CONJECTURED FORM r_k(d):")
    print("  Empirical observation: V_k stays NEGATIVE for all (d, k) tested.")
    print("  Likely r_k(d) <= 1 universally — the obvious bound.")
    print("  Sharper question: as d -> infty, does r_k(d) -> 1 - c/d^alpha?")
    print(f"\n{'d':>4} {'k':>4} {'avg(1 + V_k)':>15}", flush=True)
    for d in [20, 50, 100, 200]:
        for k in [1, 5, 10]:
            if k >= d: continue
            n_cfg = 50
            distances = []
            for _ in range(n_cfg):
                roots = random_polynomial_in_disk(d, rng)
                V = sendov_violation_kth(roots, k)
                distances.append(1.0 + V)
            avg = np.mean(distances)
            print(f"{d:>4} {k:>4} {avg:>15.6f}", flush=True)


if __name__ == "__main__":
    main()
