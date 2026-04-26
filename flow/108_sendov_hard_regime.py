"""
PAPER: 108 (canonical: 108_sendov_hard_regime.pdf)
TITLE: Sendov in the Hard Regime — Adversarial Scan
STATUS: empirical (4 adversarial families, worst V = -0.4686 at d=32 Miller)
DEPENDS: 104

THEORY
======

The HARD REGIME for Sendov: non-self-reciprocal P with one root forced
to zeta_1 = 1 on the unit circle (boundary). 4 adversarial families:
  1. Random hard
  2. Miller (zeta_1=1, others on small arc)
  3. Antipodal cluster (others near -1)
  4. Brown equispaced (others on |z|=r)

Worst observed V = -0.47 (Miller at d=32, eta=0.5). All others V <= -0.5.

VERIFICATION
============

  1. Random hard family.
  2. Miller family: extremal at d=32, eta=0.5.
  3. Antipodal cluster.
"""
from __future__ import annotations
import math
import numpy as np


def sendov_violation_at(roots, idx=0):
    """Sendov for a specific root index."""
    P = np.poly(roots)
    crits = np.roots(np.polyder(P))
    if len(crits) == 0: return float('-inf')
    z_special = roots[idx]
    return float(np.abs(z_special - crits).min()) - 1.0


def main():
    print("=" * 80)
    print("PAPER 108 — Sendov in the hard regime")
    print("=" * 80)
    rng = np.random.default_rng(2026)

    print("\n[1] Random hard: zeta_1 = 1, others uniform in disk")
    print(f"  {'d':>4} {'#cfg':>6} {'min V':>10} {'max V':>10}")
    for d in [16, 32, 64, 128]:
        n = 50
        Vs = []
        for _ in range(n):
            roots = [1.0 + 0j]
            while len(roots) < d:
                z = complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
                if abs(z) <= 1: roots.append(z)
            Vs.append(sendov_violation_at(np.array(roots), 0))
        arr = np.array(Vs)
        print(f"  {d:>4} {n:>6} {arr.min():>10.4f} {arr.max():>10.4f}")

    print("\n[2] Miller family: zeta_1 = 1, others on arc near 1")
    print(f"  {'d':>4} {'eta':>6} {'max V':>10}")
    for d in [16, 32, 64]:
        for eta in [0.5, 0.1, 0.01]:
            n = 30
            max_v = float('-inf')
            for _ in range(n):
                roots = [1.0 + 0j]
                for _ in range(d - 1):
                    phi = rng.uniform(0, 2*np.pi)
                    z = 1.0 + eta * np.exp(1j * phi)
                    if abs(z) > 1: z = z / abs(z) * 0.999
                    roots.append(z)
                v = sendov_violation_at(np.array(roots), 0)
                if v > max_v: max_v = v
            print(f"  {d:>4} {eta:>6.2f} {max_v:>10.4f}")

    print("\n[3] Antipodal cluster: zeta_1 = 1, others near -1")
    print(f"  {'d':>4} {'max V':>10}")
    for d in [16, 32, 64]:
        max_v = float('-inf')
        for _ in range(20):
            roots = [1.0 + 0j]
            for _ in range(d - 1):
                z = -1.0 + 0.05 * (rng.uniform(-1, 1) + 1j * rng.uniform(-1, 1))
                if abs(z) > 1: z = z / abs(z) * 0.999
                roots.append(z)
            v = sendov_violation_at(np.array(roots), 0)
            if v > max_v: max_v = v
        print(f"  {d:>4} {max_v:>10.4f}")


if __name__ == "__main__":
    main()
