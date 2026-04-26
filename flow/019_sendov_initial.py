"""
PAPER: 019 (canonical: 09pandrosion_sendov.pdf)
TITLE: Sendov's Conjecture: Initial Pandrosion Treatment
STATUS: empirical (refined in papers 053, 104, 108, 113)
DEPENDS: 011

THEORY
======

SENDOV'S CONJECTURE (1962):
For monic P of degree d >= 2 with all roots in {|z| <= 1}, every root zeta
has a critical point xi (with P'(xi) = 0) within distance 1:
  min_xi |zeta - xi| <= 1.

PANDROSION-FORM:
The Pandrosion field F_P = P'/P has poles at roots, zeros at critical points.
Sendov is a statement about the geometric "containment" of zeros within
unit disks of poles.

VERIFICATION
============

  1. Sendov holds on random configurations in unit disk.
  2. Boundary case: tight at z^d - 1 (roots of unity).
  3. Miller's adversarial family near z = 1.
"""
from __future__ import annotations
import math
import numpy as np


def sendov_violation(roots):
    P = np.poly(roots)
    crits = np.roots(np.polyder(P))
    if len(crits) == 0: return float('-inf')
    D = np.abs(np.array(roots)[:, None] - np.array(crits)[None, :])
    return float(D.min(axis=1).max()) - 1.0


def main():
    print("=" * 80)
    print("PAPER 19 — Sendov: initial Pandrosion treatment")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Random unit-disk configurations")
    print(f"  {'d':>4} {'#cfg':>6} {'max V':>10} {'status':>10}")
    for d in [3, 5, 8, 16, 32]:
        n = 100
        max_v = float('-inf')
        for _ in range(n):
            r = rng.uniform(0, 1, size=d)
            t = rng.uniform(0, 2*math.pi, size=d)
            roots = r * np.exp(1j * t)
            v = sendov_violation(roots)
            if v > max_v: max_v = v
        ok = "OK" if max_v <= 1e-6 else "VIOLATES"
        print(f"  {d:>4} {n:>6} {max_v:>10.4e} {ok:>10}")

    print("\n[2] Roots of unity z^d - 1 (boundary case)")
    print(f"  {'d':>4} {'V':>14}")
    for d in [3, 5, 8, 16, 32]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        v = sendov_violation(roots)
        print(f"  {d:>4} {v:>14.4e}")

    print("\n[3] Miller adversarial family (zeta_1 = 1, others on small arc near 1)")
    print(f"  {'d':>4} {'eta':>6} {'V':>10}")
    for d in [10, 20, 50]:
        for eta in [0.5, 0.1, 0.01]:
            roots = [1.0 + 0j]
            for k in range(d - 1):
                phi = rng.uniform(0, 2*np.pi)
                z = 1.0 + eta * np.exp(1j * phi)
                if abs(z) > 1: z = z / abs(z) * 0.999
                roots.append(z)
            v = sendov_violation(np.array(roots))
            print(f"  {d:>4} {eta:>6.2f} {v:>10.4e}")


if __name__ == "__main__":
    main()
