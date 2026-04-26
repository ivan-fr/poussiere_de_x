"""
LONELY RUNNER CONJECTURE via Pandrosion-trigonometric polynomials.

Conjecture (Wills 1967, Bienia-Goddyn-Gvozdjak-Sebo-Tarsi 1998):
For any n distinct nonzero integer speeds v_1, ..., v_n, there exists time t
such that for every i, the position v_i*t (mod 1) is at distance >= 1/(n+1)
from 0 (and hence from every other v_j*t mod 1 by symmetry).

Equivalently: || v_i * t || >= 1/(n+1) for all i, where ||x|| = dist(x, Z).

PROVED for n <= 6 runners (so n+1 <= 7 speeds).
OPEN for n >= 7.

Pandrosion approach:
  - Encode positions on circle via z_i(t) = exp(2*pi*i * v_i * t).
  - Define the "loneliness polynomial":
       L(t) = product over i of (1 - cos(2*pi*v_i*t)).
    L(t) is small iff some runner is near 0; L(t) maximal iff all are loneliest.
  - Pandrosion-Hadamard: relate L(t) to a Gram-determinant on the speed-vectors.
  - Numerical scan: maximize min_i ||v_i*t|| over t in [0, T] for many speed sets.
"""
from __future__ import annotations
import math
import numpy as np


def loneliness_score(speeds, t):
    """The loneliness gap at time t: gap = min_i ||v_i * t||
    where ||x|| = dist(x, Z)."""
    positions = np.array(speeds, dtype=float) * t
    fractional = positions - np.round(positions)
    distances = np.abs(fractional)
    return float(distances.min())


def max_loneliness_over_t(speeds, T_max=None, n_grid=10000):
    """Compute max_t min_i ||v_i * t|| over t in [0, T_max].
    Default T_max = 1 / gcd(speeds)."""
    if T_max is None:
        T_max = 1.0  # by symmetry, period of system divides 1/gcd
    ts = np.linspace(0, T_max, n_grid)
    best = 0.0
    best_t = 0.0
    for t in ts[1:-1]:
        s = loneliness_score(speeds, t)
        if s > best:
            best = s
            best_t = t
    return best, best_t


def conjecture_target(n):
    """Lonely runner conjecture: max_t min_i ||v_i*t|| >= 1/(n+1) for n speeds."""
    return 1.0 / (n + 1)


def lehmer_search(n_speeds, n_trials, max_speed, rng):
    """Search for adversarial speed sets with smallest max-loneliness.
    Conjecture: max-loneliness >= 1/(n+1) always."""
    worst_score = float('inf')
    worst_speeds = None
    for _ in range(n_trials):
        # Pick n distinct integer speeds
        speeds = sorted(set(rng.integers(1, max_speed + 1, size=2*n_speeds)))
        if len(speeds) < n_speeds:
            continue
        speeds = speeds[:n_speeds]
        score, _ = max_loneliness_over_t(speeds, T_max=1.0, n_grid=10000)
        if score < worst_score:
            worst_score = score
            worst_speeds = speeds
    return worst_score, worst_speeds


def pandrosion_loneliness_polynomial(speeds, t):
    """L(t) := prod_i (1 - cos(2*pi*v_i*t)). 0 iff some runner at 0.

    Logarithm: log L(t) = sum_i log(1 - cos(2*pi*v_i*t)) = sum_i log(2 sin^2(pi v_i t))
                       = sum_i [log 2 + 2 log |sin(pi v_i t)|].
    """
    s = 0.0
    for v in speeds:
        ct = math.cos(2 * math.pi * v * t)
        s += math.log(max(1 - ct, 1e-300))
    return s


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 95, flush=True)
    print("LONELY RUNNER CONJECTURE — Pandrosion attack", flush=True)
    print("=" * 95, flush=True)

    # n+1 runners total, n other speeds (excluding the lonely runner who has speed 0).
    # Conjecture: max_t min_i ||v_i*t|| >= 1/(n+1).

    print("\n[1] KNOWN CASES (n_other = number of OTHER speeds, conjecture proven n_other <= 6):")
    print(f"{'n_other':>9} {'1/(n+1)':>12} {'speeds tested':>18} {'best gap':>12} {'status':>15}",
          flush=True)
    rng = np.random.default_rng(20260711)
    for n_other in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12]:
        target = conjecture_target(n_other)
        # adversarial scan
        worst, worst_speeds = lehmer_search(n_other, n_trials=100, max_speed=20, rng=rng)
        status = "OK (>= 1/(n+1))" if worst >= target - 1e-6 else f"VIOLATES by {target - worst:.4e}"
        print(f"{n_other:>9} {target:>12.6f} 100 random sets   {worst:>12.6f} {status:>15}",
              flush=True)

    print("\n[2] CLASSIC HARD CASES (literature):")
    classic = [
        ("Bohr 1924", [1, 2, 3]),
        ("Wills 1967 hard", [1, 2, 3, 4]),
        ("BGGST 1998 sharp", [1, 2, 3, 4, 5, 6]),
        ("equal-step n=7", [1, 2, 3, 4, 5, 6, 7]),
        ("equal-step n=8", [1, 2, 3, 4, 5, 6, 7, 8]),
        ("Fibonacci n=7", [1, 2, 3, 5, 8, 13, 21]),
        ("Fibonacci n=8", [1, 2, 3, 5, 8, 13, 21, 34]),
    ]
    print(f"{'name':>20} {'speeds':>30} {'gap':>10} {'1/(n+1)':>10}", flush=True)
    for name, speeds in classic:
        gap, t_opt = max_loneliness_over_t(speeds, n_grid=20000)
        target = 1.0 / (len(speeds) + 1)
        print(f"{name:>20} {str(speeds):>30} {gap:>10.6f} {target:>10.6f}",
              flush=True)

    print("\n[3] PANDROSION REFORMULATION via L(t):")
    print("  log L(t) = sum_i log(1 - cos(2*pi*v_i*t)) measures total proximity")
    print("  Lonely instant t* = argmin |log L(t)| (all runners far from each other)")
    print(f"\n{'speeds':>20} {'gap (t*)':>12} {'log L(t*)':>14}", flush=True)
    for speeds in [[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 3, 5]]:
        gap, t_opt = max_loneliness_over_t(speeds, n_grid=10000)
        logL = pandrosion_loneliness_polynomial(speeds, t_opt)
        print(f"{str(speeds):>20} {gap:>12.6f} {logL:>14.4f}", flush=True)

    # KEY DISCOVERY ATTEMPT: Pandrosion-Hadamard for Lonely Runner
    # The "Pandrosion" of the speed-set is via the polynomial
    #   P(z) = prod_i (z^{v_i} - 1)
    # The conjecture relates to the discriminant of P modulo cyclotomic factors.
    print("\n[4] KEY DISCOVERY ATTEMPT: speed polynomial P(z) = prod (z^v_i - 1):")
    print("  Roots are v_i-th roots of unity for each i.")
    print("  Conjecture conjecturally relates to common-rotation cycles.")
    for speeds in [[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6, 7]]:
        # Build polynomial coefficients
        P = np.array([1])
        for v in speeds:
            factor = np.zeros(v + 1)
            factor[0] = 1
            factor[-1] = -1
            P = np.convolve(P, factor)
        # Find total degree, common roots etc
        roots = np.roots(P)
        n_unit = sum(1 for r in roots if abs(abs(r) - 1) < 1e-6)
        n_total = len(roots)
        gap, _ = max_loneliness_over_t(speeds, n_grid=20000)
        print(f"  speeds={speeds}: deg P = {n_total}, |unit roots| = {n_unit}, gap = {gap:.6f}",
              flush=True)


if __name__ == "__main__":
    main()
