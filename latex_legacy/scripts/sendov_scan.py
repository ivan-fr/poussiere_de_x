"""Scan Sendov's conjecture: for polynomials with all roots in unit disk,
verify that every root has a critical point within distance 1.

Test: max_j min_xi |zeta_j - xi| where zeta_j are roots and xi are critical points.
Sendov holds iff this max is <= 1.
We record V(P) = max_j min_xi |zeta_j - xi| - 1; Sendov ⇔ V(P) <= 0.
"""
from __future__ import annotations
import math, time
import numpy as np


def sample_disk_uniform(d, rng):
    """Sample d random complex numbers in the unit disk uniformly."""
    # Rejection sampling
    out = []
    while len(out) < d:
        z = rng.uniform(-1, 1) + 1j * rng.uniform(-1, 1)
        if abs(z) <= 1.0:
            out.append(z)
    return np.array(out, dtype=complex)


def sample_disk_BW_like(d, rng):
    """Approximate Bombieri-Weyl-on-disk: roots concentrate near unit circle.
    Sample r ~ Beta(d, 1) (CDF: r^d), theta ~ Uniform.
    """
    r = rng.beta(d, 1.0, size=d)  # concentrated near r=1
    theta = rng.uniform(0, 2*np.pi, size=d)
    return r * np.exp(1j * theta)


def sendov_violation(roots):
    """Build P from roots, compute critical points, return V(P)."""
    d = len(roots)
    # Build polynomial coefficients (descending order for numpy.poly)
    P = np.poly(roots)  # P[0] z^d + P[1] z^{d-1} + ... + P[d]
    # Derivative coefficients (descending)
    deriv = np.polyder(P)
    if len(deriv) == 0:
        return float('inf')
    crits = np.roots(deriv)
    if len(crits) == 0:
        return float('inf')
    # For each root, find min distance to a critical point
    # Vectorized: |roots[:, None] - crits[None, :]|
    D = np.abs(roots[:, None] - crits[None, :])
    min_dists = D.min(axis=1)
    max_min = float(min_dists.max())
    return max_min - 1.0


def scan(d, n_polys, sampler='uniform', rng=None):
    if rng is None:
        rng = np.random.default_rng(20260427 + d)
    Vs = []
    for _ in range(n_polys):
        if sampler == 'uniform':
            roots = sample_disk_uniform(d, rng)
        elif sampler == 'BW':
            roots = sample_disk_BW_like(d, rng)
        else:
            raise ValueError(sampler)
        try:
            V = sendov_violation(roots)
        except (np.linalg.LinAlgError, ValueError):
            continue
        if math.isfinite(V):
            Vs.append(V)
    return np.array(Vs)


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    CASES = [
        (9,    20000),
        (16,   20000),
        (32,   10000),
        (64,    5000),
        (128,   2000),
        (256,    500),
        (512,    100),
        (1024,    30),
    ]
    t0 = time.perf_counter()
    print("=" * 110, flush=True)
    print("Sendov scan on Bombieri-Weyl-like uniform-on-disk polynomials", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>5} {'#polys':>8} {'min V':>10} {'median V':>10} {'max V':>10} "
          f"{'#violations':>12} {'time':>9}", flush=True)
    print("-" * 110, flush=True)
    summary = {}
    for d, npolys in CASES:
        t1 = time.perf_counter()
        Vs = scan(d, npolys, sampler='BW')
        elapsed = time.perf_counter() - t1
        if len(Vs) == 0:
            print(f"{d:>5} {npolys:>8}  no valid", flush=True)
            continue
        min_V = float(Vs.min())
        med_V = float(np.median(Vs))
        max_V = float(Vs.max())
        n_viol = int((Vs > 0).sum())
        flag = "  <-- COUNTEREXAMPLE" if n_viol > 0 else ""
        print(f"{d:>5} {npolys:>8} {min_V:>10.4f} {med_V:>10.4f} {max_V:>10.4f} "
              f"{n_viol:>12} {elapsed:>8.1f}s{flag}",
              flush=True)
        summary[d] = (min_V, med_V, max_V, n_viol, npolys, elapsed)

    # Adversarial: P_d(z) = z^{d-1}(z - 1) — known asymptotic equality case
    print("\n" + "=" * 110, flush=True)
    print("Adversarial: P_d(z) = z^{d-1}(z - 1) — Sendov boundary case", flush=True)
    print("=" * 110, flush=True)
    for d in [9, 16, 32, 64, 128, 256, 512]:
        roots = np.array([0.0+0j]*(d-1) + [1.0+0j])
        V = sendov_violation(roots)
        print(f"  d = {d:>4}: V = {V:.6f}", flush=True)

    # Adversarial: roots of unity scaled by 0.999
    print("\n" + "=" * 110, flush=True)
    print("Adversarial: scaled roots of unity 0.999 * exp(2 pi i k/d)", flush=True)
    print("=" * 110, flush=True)
    for d in [9, 16, 32, 64, 128, 256]:
        roots = 0.999 * np.exp(2j*np.pi*np.arange(d)/d)
        V = sendov_violation(roots)
        print(f"  d = {d:>4}: V = {V:.6f}", flush=True)

    # LaTeX-ready table
    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table for paper 104 §3.2:", flush=True)
    print("=" * 70, flush=True)
    for d, _ in CASES:
        if d in summary:
            min_V, med_V, max_V, n_viol, npolys, elapsed = summary[d]
            status = ("\\textbf{COUNTEREXAMPLE}" if n_viol > 0 else "$\\le 0$ everywhere")
            print(f"  ${d}$  & ${npolys:,}$ & ${max_V:+.4f}$ & ${med_V:+.4f}$ & "
                  f"${elapsed:.1f}$\\,s & {status} \\\\",
                  flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
