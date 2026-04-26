"""Sendov in the HARD regime: non-self-reciprocal polynomials with
one root at zeta_1 = 1 (the unit-circle boundary case) and (d-1) roots
distributed adversarially in the disk.

Key adversarial families:
  1. Miller's family: zeta_1 = 1, others on a small arc near zeta_1
  2. zeta_1 = 1, others i.i.d. uniform in disk (random hard case)
  3. zeta_1 = 1, others clustered near antipode -1 (paper of Borcea-Brandt)
  4. Brown's family: zeta_1 = 1, others equally spaced on |z| = r < 1

Test V(P) = max_j min_k |zeta_j - xi_k| - 1, focusing on V at zeta_1 = 1
since that's where Sendov is tightest.
"""
from __future__ import annotations
import math, time
import numpy as np


def sendov_violation_at(roots, special_idx=0):
    """Sendov violation, restricted to a single special root index."""
    P = np.poly(roots)
    deriv = np.polyder(P)
    if len(deriv) == 0:
        return float('inf')
    crits = np.roots(deriv)
    if len(crits) == 0:
        return float('inf')
    z_special = roots[special_idx]
    min_dist = float(np.abs(z_special - crits).min())
    return min_dist - 1.0


def sendov_violation_full(roots):
    """Full Sendov: max over all roots."""
    P = np.poly(roots)
    deriv = np.polyder(P)
    crits = np.roots(deriv)
    if len(crits) == 0:
        return float('inf')
    D = np.abs(roots[:, None] - crits[None, :])
    return float(D.min(axis=1).max()) - 1.0


def miller_family(d, eta, rng):
    """Miller's family approximation: zeta_1 = 1, others on an arc {1 + eta * e^{i*phi}}."""
    roots = [1.0 + 0j]
    for _ in range(d - 1):
        phi = rng.uniform(0, 2*np.pi)
        z = 1.0 + eta * np.exp(1j * phi)
        if abs(z) <= 1:
            roots.append(z)
        else:
            roots.append(z / abs(z) * 0.999)  # project to slightly inside
    return np.array(roots)


def random_hard(d, rng):
    """zeta_1 = 1, others uniform in unit disk."""
    roots = [1.0 + 0j]
    while len(roots) < d:
        z = rng.uniform(-1, 1) + 1j * rng.uniform(-1, 1)
        if abs(z) <= 1:
            roots.append(z)
    return np.array(roots)


def antipodal_cluster(d, rng):
    """zeta_1 = 1, others clustered near -1."""
    roots = [1.0 + 0j]
    for _ in range(d - 1):
        z = -1.0 + 0.05 * (rng.uniform(-1, 1) + 1j * rng.uniform(-1, 1))
        if abs(z) > 1:
            z = z / abs(z) * 0.999
        roots.append(z)
    return np.array(roots)


def brown_family(d, r):
    """zeta_1 = 1, (d-1) equispaced on |z| = r."""
    roots = [1.0 + 0j]
    for k in range(d - 1):
        roots.append(r * np.exp(2j * np.pi * k / (d - 1)))
    return np.array(roots)


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 95, flush=True)
    print("Sendov HARD regime: zeta_1 = 1, (d-1) other roots adversarial",
          flush=True)
    print("=" * 95, flush=True)

    # 1. Random hard
    print("\n1. zeta_1 = 1, others i.i.d. uniform in disk", flush=True)
    print(f"{'d':>5} {'#polys':>7} {'min V_at_1':>12} {'median':>10} {'max V_at_1':>12} {'#viol':>6}",
          flush=True)
    for d in [9, 16, 32, 64, 128, 256, 512, 1024]:
        rng = np.random.default_rng(20260428 + d)
        n = max(50, 5000 // d)
        Vs = []
        for _ in range(n):
            roots = random_hard(d, rng)
            v = sendov_violation_at(roots, 0)
            if math.isfinite(v):
                Vs.append(v)
        if not Vs:
            continue
        a = np.array(Vs)
        print(f"{d:>5} {n:>7} {float(a.min()):>12.4f} "
              f"{float(np.median(a)):>10.4f} {float(a.max()):>12.4f} "
              f"{int((a > 0).sum()):>6}",
              flush=True)

    # 2. Miller arc family — the historical adversary
    print("\n2. Miller's family: zeta_1 = 1, others on arc of radius eta near 1",
          flush=True)
    for d in [16, 32, 64, 128, 256]:
        for eta in [0.5, 0.1, 0.01]:
            rng = np.random.default_rng(424242 + d * 100 + int(eta * 1000))
            Vs = []
            for _ in range(200):
                roots = miller_family(d, eta, rng)
                v = sendov_violation_at(roots, 0)
                if math.isfinite(v):
                    Vs.append(v)
            if Vs:
                a = np.array(Vs)
                print(f"  d={d:>4}, eta={eta:>5.2f}: "
                      f"max V = {float(a.max()):+.4f}, median = {float(np.median(a)):+.4f}",
                      flush=True)

    # 3. Antipodal cluster
    print("\n3. zeta_1 = 1, (d-1) clustered near -1", flush=True)
    for d in [16, 32, 64, 128, 256]:
        rng = np.random.default_rng(7777 + d)
        Vs = []
        for _ in range(200):
            roots = antipodal_cluster(d, rng)
            v = sendov_violation_at(roots, 0)
            if math.isfinite(v):
                Vs.append(v)
        if Vs:
            a = np.array(Vs)
            print(f"  d={d:>4}: max V = {float(a.max()):+.4f}, "
                  f"median = {float(np.median(a)):+.4f}",
                  flush=True)

    # 4. Brown family — equispaced inner ring
    print("\n4. Brown's family: zeta_1 = 1, (d-1) on |z| = r equispaced",
          flush=True)
    for d in [16, 32, 64, 128, 256, 512]:
        for r in [0.999, 0.99, 0.9, 0.5]:
            roots = brown_family(d, r)
            v = sendov_violation_at(roots, 0)
            print(f"  d={d:>4}, r={r:>5.3f}: V_at_1 = {v:+.5f}", flush=True)

    # 5. THE DEEPEST ADVERSARIAL: cyclotomic-type with one root at 1
    print("\n5. Cyclotomic with one root forced to 1: P(z) = (z-1)(z^{d-1} - 1)/(z-1)",
          flush=True)
    print("   = (z^{d-1} - 1) * (z-1) / (z-1) — but degenerate, so use:")
    print("   P(z) = (z-1) * prod_{k=1}^{d-1} (z - omega_k) with omega_k = e^{2 pi i k/(d-1)}",
          flush=True)
    for d in [10, 16, 32, 64, 128, 256, 512, 1024]:
        roots = [1.0 + 0j]
        for k in range(1, d):
            roots.append(np.exp(2j * np.pi * k / (d - 1)))
        roots = np.array(roots)
        v = sendov_violation_at(roots, 0)
        print(f"  d={d:>4}: V_at_1 = {v:+.6f}", flush=True)


if __name__ == "__main__":
    main()
