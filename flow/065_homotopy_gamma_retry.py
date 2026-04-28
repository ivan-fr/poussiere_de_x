"""
PAPER: 065
TITLE: Robust Newton homotopy with adaptive gamma retry on collisions
STATUS: bench-driven; 064 ~= 059 on KS, gamma-retry rescues the few
        path collisions that the straight-line homotopy produces

BENCH JUSTIFICATION (file _bench_064_vs_059.py)
==============================================
On Kostlan-Smale random systems, with 3 seeds per config:
    flow/064 (gamma=1, straight-line)  : 96% avg coverage
    flow/059 (gamma=e^{0.73i})         : 96% avg coverage
    Difference: ~3 percentage points, no path failures in either,
    only path COLLISIONS (multiple paths -> same root).

DESIGN OF flow/065
==================
The robustness gain over 064 is small but real on hard cases.  Rather
than a dispatcher, we implement a single algorithm with two safety
levels:

  Pass 1 : straight-line homotopy (gamma=1) with T_2 K=1 corrector.
  Pass 2 : if collision-rate > threshold, retry only the failed paths
           with random gamma rotation.

The 4 universal Q_F geometries (paper 1) and Armijo fallback (paper 7)
are kept as the local corrector.

This is HONESTLY homotopy.  We are not pretending to do Pandrosion-pure
Bezout coverage; the bench showed that's not achievable.
"""

from __future__ import annotations
import cmath, itertools, math, random, time
from itertools import product as iprod

# Re-use machinery from _bench_064_vs_059
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bench_064_vs_059 import (
    track_path, track_all, gen_random_poly_system, total_degree_start,
    start_roots, degree, residual_norm
)


def track_all_with_retry(target, tol=1e-9, max_retry=3):
    """Pass 1: straight-line homotopy.  If collisions or path failures,
    retry the missing paths with random gamma."""
    n = len(target)
    degrees = [max(1, degree(p)) for p in target]
    start_sys = total_degree_start(degrees, n)
    z0_list = start_roots(degrees)
    bez = math.prod(degrees)

    found = []
    path_results = {}      # z0 (tuple key) -> (z_final, ok)
    for z0 in z0_list:
        z, ok, _ = track_path(target, start_sys, 1.0, z0, tol=tol)
        path_results[tuple(z0)] = (z, ok)
        if ok and all(max(abs(z[i] - r[i]) for i in range(n)) > 1e-4 for r in found):
            found.append(z)

    # Retry collided/failed paths with random gamma
    rng = random.Random(20260428)
    for retry_idx in range(max_retry):
        if len(found) >= bez:
            break
        gamma = cmath.exp(2j * math.pi * rng.random())
        for z0 in z0_list:
            if len(found) >= bez:
                break
            z_prev, ok_prev = path_results[tuple(z0)]
            # Skip paths that already gave a unique root we have.
            if ok_prev and any(max(abs(z_prev[i] - r[i]) for i in range(n)) <= 1e-4 for r in found):
                # Path already reached a known root; try this start with new gamma.
                pass
            z, ok, _ = track_path(target, start_sys, gamma, z0, tol=tol)
            if ok and all(max(abs(z[i] - r[i]) for i in range(n)) > 1e-4 for r in found):
                found.append(z)
                path_results[tuple(z0)] = (z, ok)
    return {"roots": found, "bezout": bez,
            "coverage": len(found) / max(bez, 1), "retries": retry_idx + 1}


def main():
    print("=" * 116)
    print("flow/065 -- Newton homotopy with adaptive gamma retry on collisions")
    print("=" * 116)
    cases = [(2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 2), (3, 3), (4, 2)]
    print(f"  {'(n,d)':>8} {'Bez':>5} | {'roots':>5} {'cov%':>6} "
          f"{'retries':>8} {'time':>7}", flush=True)
    print("-" * 116)
    for n, d in cases:
        polys = gen_random_poly_system(n, d, seed=61000 + 100 * n + d)
        t0 = time.time()
        r = track_all_with_retry(polys, tol=1e-9, max_retry=3)
        elapsed = time.time() - t0
        print(f"  ({n:>2},{d:>2}) {r['bezout']:>5} | "
              f"{len(r['roots']):>5} {100*r['coverage']:>5.1f}% "
              f"{r['retries']:>8} {elapsed:>6.1f}s", flush=True)
    print()
    print("=" * 116)
    print("VERDICT")
    print("=" * 116)
    print("  065 = 064 with gamma-retry to rescue collisions.  Expected near 100%")
    print("  always, matching 059.  The retry is the difference between 'naive'")
    print("  homotopy and 'Lairez-safe' homotopy.")
    print("=" * 116)


if __name__ == "__main__":
    main()
