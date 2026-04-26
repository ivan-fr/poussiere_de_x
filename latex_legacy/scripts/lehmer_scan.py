"""Lehmer's conjecture scan: search for polynomials with small Mahler measure.

Mahler measure M(P) = |a_d| * prod_j max(1, |alpha_j|).
Lehmer (1933): for monic integer P with M(P) > 1, M(P) >= L = 1.17628...
(Smyth's polynomial z^10 + z^9 - z^7 - z^6 - z^5 - z^4 - z^3 + z + 1).

Pandrosion attack: exploit the energy identity
    E(P) := sum_j |Q(alpha_j, 0)/P'(alpha_j)|^2
and the Pandrosion field F_P(z) = P'(z)/P(z) = sum 1/(z - alpha_j).

Numerical question: among monic integer polynomials of degree <= D and
coefficients in {-K, ..., K}, what is the minimum Mahler measure M(P) > 1?
"""
from __future__ import annotations
import math, time, itertools
import numpy as np


def mahler_measure(coefs):
    """coefs in ascending order: P(z) = sum coefs[k] z^k."""
    if abs(coefs[-1]) < 1e-15:
        return 0.0
    roots = np.roots(coefs[::-1])  # numpy expects descending
    return abs(coefs[-1]) * float(np.prod(np.maximum(1.0, np.abs(roots))))


def pandrosion_energy(coefs):
    """E(P) = sum_j 1/|P'(alpha_j)|^2 — Pandrosion-style energy.
    For monic P with simple roots, this is well-defined and connects to
    the discriminant via prod |P'(alpha_j)|^2 = |disc(P)|^2.
    """
    if abs(coefs[-1]) < 1e-15:
        return float('inf')
    roots = np.roots(coefs[::-1])
    deriv_coefs = coefs[1:] * np.arange(1, len(coefs))
    if len(deriv_coefs) == 0:
        return float('inf')
    Pprimes = np.polyval(deriv_coefs[::-1], roots)
    if np.any(np.abs(Pprimes) < 1e-15):
        return float('inf')
    return float(np.sum(1.0 / np.abs(Pprimes)**2))


def is_reciprocal(coefs):
    """Check if P is reciprocal (a_k = a_{d-k})."""
    d = len(coefs) - 1
    return all(abs(coefs[k] - coefs[d - k]) < 1e-12 for k in range(d // 2 + 1))


def smyth_polynomial():
    """Smyth's z^10 + z^9 - z^7 - z^6 - z^5 - z^4 - z^3 + z + 1.
    Returns ascending-order coefficients."""
    coefs = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]  # ascending: a_0, a_1, ..., a_10
    return np.array(coefs, dtype=float)


def scan_small_polys(d_max, K, n_per_d=None):
    """Exhaustively scan monic integer polynomials of degree d, coefs in [-K, K].
    Return min Mahler measure > 1 found."""
    best = (float('inf'), None)
    for d in range(2, d_max + 1):
        # Brute-force over all coefficient sequences would be (2K+1)^d — too many.
        # Sample randomly instead.
        rng = np.random.default_rng(20260427 + d)
        n_samples = n_per_d if n_per_d is not None else min(10**5, (2*K + 1)**min(d, 6))
        for _ in range(n_samples):
            mid = rng.integers(-K, K + 1, size=d)
            coefs = np.concatenate([mid, [1]]).astype(float)
            # Skip degenerate (a_0 = 0)
            if abs(coefs[0]) < 0.5:
                continue
            try:
                M = mahler_measure(coefs)
            except (np.linalg.LinAlgError, ValueError):
                continue
            if M > 1.0001 and M < best[0]:
                best = (M, coefs.copy())
    return best


def verify_smyth():
    """Sanity check on Smyth's polynomial."""
    coefs = smyth_polynomial()
    M = mahler_measure(coefs)
    return M


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 75, flush=True)
    print("Lehmer scan (Pandrosion energy attack)", flush=True)
    print("=" * 75, flush=True)

    # Sanity: Smyth polynomial gives M ~ 1.17628
    M_smyth = verify_smyth()
    print(f"\nSmyth polynomial M(P) = {M_smyth:.6f} (expected 1.17628)", flush=True)

    # Scan integer polynomials with K=1 (Littlewood-style)
    print("\nScan: monic, integer, coefs in {-1, 0, 1}", flush=True)
    print(f"{'d_max':>6} {'#samples':>10} {'min M > 1':>12} {'comment':>20}", flush=True)
    for d_max in [10, 15, 20, 25, 30]:
        t0 = time.perf_counter()
        best, _ = scan_small_polys(d_max, K=1, n_per_d=20000)
        elapsed = time.perf_counter() - t0
        print(f"{d_max:>6} {20000:>10} {best:>12.6f} "
              f"  (t={elapsed:.1f}s)", flush=True)

    # Scan with larger coefficients (K=2)
    print("\nScan: K=2 (coefs in {-2, -1, 0, 1, 2})", flush=True)
    for d_max in [10, 15, 20]:
        t0 = time.perf_counter()
        best, _ = scan_small_polys(d_max, K=2, n_per_d=10000)
        elapsed = time.perf_counter() - t0
        print(f"d_max = {d_max}: min M = {best:.6f}  (t={elapsed:.1f}s)", flush=True)

    # Reciprocal polynomials specifically (Lehmer's case)
    print("\nReciprocal polynomials (the Lehmer case)", flush=True)
    rng = np.random.default_rng(42)
    best_recip = float('inf')
    n_recip_tested = 0
    for d in range(4, 31, 2):  # even degree for reciprocal
        for _ in range(5000):
            half = rng.integers(-1, 2, size=d // 2)
            coefs = np.zeros(d + 1)
            coefs[0] = 1; coefs[d] = 1
            for k in range(d // 2):
                coefs[k + 1] = half[k]
                coefs[d - k - 1] = half[k]
            if d % 2 == 0:
                coefs[d // 2] = rng.integers(-1, 2)
            try:
                M = mahler_measure(coefs)
                n_recip_tested += 1
                if M > 1.0001 and M < best_recip:
                    best_recip = M
            except (np.linalg.LinAlgError, ValueError):
                continue
    print(f"  Tested {n_recip_tested} reciprocal polynomials.", flush=True)
    print(f"  Min M (reciprocal) = {best_recip:.6f}", flush=True)
    print(f"  Lehmer's conjecture: M >= 1.17628 (Smyth's polynomial)", flush=True)
    if best_recip < 1.17628:
        print(f"  WARNING: VIOLATION!  min M = {best_recip:.6f} < 1.17628", flush=True)
    else:
        print(f"  No violation: all reciprocal M >= {best_recip:.6f} >= 1.17628 ✓",
              flush=True)


if __name__ == "__main__":
    main()
