"""Sendov for self-reciprocal real polynomials of even degree d.

A self-reciprocal polynomial P satisfies z^d P(1/z) = P(z), so its roots come
in pairs (zeta, 1/conj(zeta)). For real self-reciprocal P, roots are either:
  - on the unit circle |zeta| = 1
  - in pairs (r, 1/r) with r real
  - in conjugate pairs (z, conj(z), 1/z, 1/conj(z))

Sendov's conjecture for self-reciprocal: every root has a critical point within distance 1.

Test:
  1. Sample real self-reciprocal P with all roots in disk (so all roots on unit circle).
  2. Compute V(P) = max_j min_k |zeta_j - xi_k| - 1.
  3. Push to d up to 64, large samples.
"""
from __future__ import annotations
import math, time
import numpy as np


def make_self_reciprocal_real(d, rng):
    """Build a real self-reciprocal polynomial of even degree d with all roots on
    the unit circle. Roots come in conjugate pairs e^{i theta}, e^{-i theta} for
    d/2 random angles theta.
    """
    n_pairs = d // 2
    thetas = rng.uniform(0, np.pi, size=n_pairs)
    roots = []
    for t in thetas:
        roots.append(np.exp(1j * t))
        roots.append(np.exp(-1j * t))
    return np.array(roots, dtype=complex)


def sendov_violation(roots):
    P = np.poly(roots)
    deriv = np.polyder(P)
    if len(deriv) == 0:
        return float('inf')
    crits = np.roots(deriv)
    if len(crits) == 0:
        return float('inf')
    D = np.abs(roots[:, None] - crits[None, :])
    return float(D.min(axis=1).max()) - 1.0


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 75, flush=True)
    print("Sendov for real self-reciprocal polynomials with all roots on |z|=1",
          flush=True)
    print("=" * 75, flush=True)
    print(f"{'d':>4} {'#polys':>7} {'min V':>10} {'median V':>10} {'max V':>10} "
          f"{'#viol':>6}", flush=True)
    print("-" * 75, flush=True)
    for d in [10, 12, 14, 16, 18, 20, 24, 32, 48, 64]:
        rng = np.random.default_rng(20260428 + d)
        Vs = []
        n_polys = max(50, 5000 // d)
        for _ in range(n_polys):
            roots = make_self_reciprocal_real(d, rng)
            V = sendov_violation(roots)
            if math.isfinite(V):
                Vs.append(V)
        if not Vs:
            print(f"{d:>4} no valid", flush=True)
            continue
        arr = np.array(Vs)
        print(f"{d:>4} {n_polys:>7} {float(arr.min()):>10.4f} "
              f"{float(np.median(arr)):>10.4f} {float(arr.max()):>10.4f} "
              f"{int((arr > 0).sum()):>6}",
              flush=True)

    # Adversarial: roots clustered near a single point on |z|=1
    print("\n" + "=" * 75, flush=True)
    print("Adversarial: roots clustered near z=1 + tiny imaginary perturbations",
          flush=True)
    print("=" * 75, flush=True)
    for d in [10, 16, 24, 48]:
        for delta in [1e-1, 1e-2, 1e-3]:
            rng = np.random.default_rng(424242 + d)
            n_pairs = d // 2
            thetas = delta * rng.uniform(-1, 1, size=n_pairs)
            roots = []
            for t in thetas:
                roots.append(np.exp(1j * t))
                roots.append(np.exp(-1j * t))
            roots = np.array(roots, dtype=complex)
            V = sendov_violation(roots)
            print(f"  d={d:>3}, delta={delta:.0e}: V = {V:+.4f}", flush=True)

    # Roots equispaced on circle (cyclotomic-like)
    print("\n" + "=" * 75, flush=True)
    print("Adversarial: equispaced roots of z^d - 1 (cyclotomic)", flush=True)
    print("=" * 75, flush=True)
    for d in [10, 16, 24, 48, 64]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        V = sendov_violation(roots)
        print(f"  d={d:>3}: V = {V:+.4f}", flush=True)


if __name__ == "__main__":
    main()
