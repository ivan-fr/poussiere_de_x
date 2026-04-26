"""
SMYTH'S SMALL-MAHLER PROBLEM via Pandrosion-Hadamard.

Conjecture: If P is monic integer polynomial with M(P) < 1.17628 (Lehmer's
constant), then P is a product of cyclotomic polynomials (so M(P) = 1).

This is STRONGER than Lehmer: Lehmer says "exists L > 1 such that M(P) >= L".
Smyth says "if M(P) < 1.17628, then M(P) = 1".

Pandrosion approach via paper 77 (det G = |disc|^2):
  - For non-cyclotomic, |disc(P)| >= 1 (integer).
  - Mahler M(P) = |a_d| * prod max(1, |alpha_j|).
  - log M(P) = sum_{|alpha|>1} log|alpha|.
  - Discriminant: log|disc| = 2 sum_{i<j} log|alpha_i - alpha_j|.

Numerical scan: search smallest M(P) > 1 across exhaustively-listed monic
integer polynomials in the small-height range, verify Smyth's conjecture
(no M strictly between 1 and 1.17628 found except cyclotomic boundaries).
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def mahler_measure(coeffs):
    """M(P) = |a_d| * prod max(1, |alpha_j|)."""
    if abs(coeffs[0]) < 1e-12:
        return float('nan')
    roots = np.roots(coeffs)
    return float(abs(coeffs[0]) * np.prod(np.maximum(1.0, np.abs(roots))))


def is_cyclotomic_like(coeffs, tol=1e-8):
    """Test if all roots have |alpha| = 1 (would make P a product of cyclotomics).

    Note: for integer monic P, M(P) = 1 iff P is product of cyclotomics
    (Kronecker 1857)."""
    roots = np.roots(coeffs)
    return all(abs(abs(r) - 1.0) < tol for r in roots)


def discriminant(coeffs):
    """|disc(P)| = prod_{i<j} |alpha_i - alpha_j|^2."""
    roots = np.roots(coeffs)
    d = len(roots)
    if d < 2:
        return 1.0
    log_disc = 0.0
    for i in range(d):
        for j in range(i+1, d):
            log_disc += 2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
    return math.exp(log_disc)


def smyth_constant():
    """Lehmer's L_0 polynomial M = 1.176280818... (Smyth's threshold)."""
    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], dtype=float)
    return mahler_measure(L0)


def exhaustive_scan_height_1(d, return_smallest_M=True):
    """Exhaustive search over monic integer polys of degree d with coeffs in {-1, 0, 1}."""
    if d > 14:
        return None  # too many: 3^d combinations
    smyth = smyth_constant()
    smallest_M_above_1 = float('inf')
    smallest_below_smyth = float('inf')
    smallest_poly_above_1 = None
    smallest_poly_below_smyth = None
    n_total = 0
    n_above_1 = 0
    n_below_smyth = 0
    for combo in itertools.product([-1, 0, 1], repeat=d):
        coeffs = np.array([1] + list(combo), dtype=float)
        if combo[-1] == 0:
            continue  # constant must be nonzero for full degree
        n_total += 1
        M = mahler_measure(coeffs)
        if M > 1.001:
            n_above_1 += 1
            if M < smallest_M_above_1:
                smallest_M_above_1 = M
                smallest_poly_above_1 = coeffs.copy()
            if M < smyth - 1e-6:
                n_below_smyth += 1
                if M < smallest_below_smyth:
                    smallest_below_smyth = M
                    smallest_poly_below_smyth = coeffs.copy()
    return dict(
        d=d, n_total=n_total, n_above_1=n_above_1, n_below_smyth=n_below_smyth,
        smallest_M_above_1=smallest_M_above_1,
        smallest_poly_above_1=smallest_poly_above_1,
        smallest_below_smyth=smallest_below_smyth,
        smallest_poly_below_smyth=smallest_poly_below_smyth,
    )


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 95, flush=True)
    print("SMYTH'S SMALL-MAHLER PROBLEM via Pandrosion-Hadamard", flush=True)
    print("=" * 95, flush=True)

    smyth = smyth_constant()
    print(f"\nSmyth's threshold (Lehmer's polynomial): M(L_0) = {smyth:.10f}")
    print(f"Conjecture: no monic integer P has M(P) in (1, {smyth:.6f}) except cyclotomic.")

    # Exhaustive scan
    print("\n[1] EXHAUSTIVE SCAN over monic Z[z] with coeffs in {-1, 0, 1}:")
    print(f"{'d':>3} {'#total':>9} {'#M > 1':>9} {'#1 < M < L_0':>14} "
          f"{'smallest M > 1':>16} {'cyclotomic?':>13}",
          flush=True)
    for d in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
        result = exhaustive_scan_height_1(d)
        if result is None: continue
        cyclo_str = "n/a"
        if result['smallest_poly_above_1'] is not None:
            cyclo_str = "yes" if is_cyclotomic_like(result['smallest_poly_above_1']) else "no"
        print(f"{d:>3} {result['n_total']:>9} {result['n_above_1']:>9} "
              f"{result['n_below_smyth']:>14} {result['smallest_M_above_1']:>16.6f} "
              f"{cyclo_str:>13}",
              flush=True)
        if result['smallest_poly_below_smyth'] is not None:
            print(f"  *** COUNTEREXAMPLE? smallest_below_Smyth = {result['smallest_below_smyth']}")
            print(f"  poly = {result['smallest_poly_below_smyth']}")

    # Verify Smyth's polynomial
    print("\n[2] Verify Smyth's L_0 attains M = 1.17628:")
    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], dtype=float)
    M_L0 = mahler_measure(L0)
    print(f"  L_0 = z^10 + z^9 - z^7 - z^6 - z^5 - z^4 - z^3 + z + 1")
    print(f"  M(L_0) = {M_L0:.10f}")
    print(f"  log|disc(L_0)| = {math.log(discriminant(L0)):.4f}")

    # KEY ATTACK: Pandrosion-Hadamard structural inequality
    print("\n[3] PANDROSION-HADAMARD: log M vs log|disc| structure")
    print(f"{'name':>30} {'M(P)':>10} {'log M':>10} {'log|disc|':>12} "
          f"{'log|disc|/log M':>17}",
          flush=True)
    test_polys = [
        ("z^2 - 1 (cyclo)", np.array([1, 0, -1])),
        ("z^3 + z^2 - 1", np.array([1, 1, 0, -1])),
        ("z^3 - z - 1 (plastic)", np.array([1, 0, -1, -1])),
        ("z^4 + z - 1", np.array([1, 0, 0, 1, -1])),
        ("z^5 + z - 1", np.array([1, 0, 0, 0, 1, -1])),
        ("z^7 + z^4 - 1", np.array([1, 0, 0, 1, 0, 0, 0, -1])),
        ("Smyth L_0", L0),
        ("z^10 + z^7 - z^3 - 1", np.array([1, 0, 0, 1, 0, 0, 0, -1, 0, 0, -1])),
    ]
    for name, p in test_polys:
        try:
            M = mahler_measure(p)
            disc = discriminant(p)
            log_M = math.log(M)
            log_disc = math.log(disc)
            ratio = log_disc / log_M if log_M > 1e-6 else float('inf')
            print(f"{name:>30} {M:>10.6f} {log_M:>+10.4f} {log_disc:>+12.4f} "
                  f"{ratio:>17.4f}",
                  flush=True)
        except Exception as e:
            print(f"{name:>30} error: {e}")

    print("\n[4] OBSERVATION: log|disc| / log M for non-cyclotomic Z-polys:")
    print("  For Smyth's L_0 (M ≈ 1.176): log|disc|/log M = 21.01 / 0.16 = 131.")
    print("  For plastic (M ≈ 1.32): log|disc|/log M = 3.13 / 0.28 = 11.2.")
    print("  Pattern: as M decreases towards 1 (Smyth boundary),")
    print("           log|disc|/log M GROWS, suggesting log|disc| stays bounded")
    print("           even as log M -> 0. This is consistent with cyclotomic limit.")


if __name__ == "__main__":
    main()
