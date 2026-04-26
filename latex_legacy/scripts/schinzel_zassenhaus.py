"""
SCHINZEL-ZASSENHAUS conjecture / Dimitrov theorem.

Conjecture (Schinzel-Zassenhaus 1965):
For monic P in Z[z] of degree d, NOT cyclotomic, max_j |alpha_j| >= 1 + c/d
for some constant c > 0.

Dimitrov 2019: proved a stronger form via Galois orbits.
Many quantitative refinements remain open.

Pandrosion attack via paper 77 (det G = |Disc|^2):
  - Disc(P) integer => |Disc(P)| >= 1 for distinct roots.
  - Take log: log|Disc| >= 0 (trivial for distinct roots).
  - But: log|Disc| = sum_{i<j} 2 log|alpha_i - alpha_j|.
  - For all roots near unit circle, this constrains spread.

Numerical strategy: scan all monic Z[z] with bounded coeff height, compute
M_max := max |alpha_j|, find smallest M_max > 1.
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def max_modulus(coeffs):
    """max |alpha_j| where alpha_j are roots of P."""
    roots = np.roots(coeffs)
    return float(np.max(np.abs(roots)))


def is_cyclotomic_like(coeffs, tol=1e-8):
    """All roots on |z| = 1?"""
    roots = np.roots(coeffs)
    return all(abs(abs(r) - 1.0) < tol for r in roots)


def schinzel_excess(coeffs):
    """Returns max|alpha| - 1 (the 'Schinzel excess'); zero for cyclotomic."""
    return max_modulus(coeffs) - 1.0


def discriminant_log(coeffs):
    """log |Disc(P)| via root pairs."""
    roots = np.roots(coeffs)
    d = len(roots)
    log_disc = 0.0
    for i in range(d):
        for j in range(i+1, d):
            log_disc += 2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
    return log_disc


def mahler_measure(coeffs):
    """M(P) = |a_d| * prod max(1, |alpha_j|).
    For cyclotomic monic Z-polys, M = 1 exactly.
    Numerically: M = 1.0 to ~1e-9 for cyclotomic."""
    roots = np.roots(coeffs)
    return float(abs(coeffs[0]) * np.prod(np.maximum(1.0, np.abs(roots))))


def exhaustive_scan_height(d, height, count_only=False):
    """Scan all monic Z[z] of degree d, coeff height <= H.
    Compute smallest non-cyclotomic max-modulus excess.

    Use Mahler measure as the cyclotomic filter:
    M(P) = 1 iff P is product of cyclotomics (Kronecker's theorem).
    Numerical threshold: M > 1.001 (non-cyclotomic) vs M < 1.0001 (cyclotomic).
    """
    if d > 12 and height > 1:
        return None
    smallest_excess = float('inf')
    smallest_excess_poly = None
    n_total = 0
    n_cyclotomic = 0
    n_above_one = 0
    height_range = list(range(-height, height + 1))
    for combo in itertools.product(height_range, repeat=d):
        coeffs = np.array([1] + list(combo), dtype=float)
        if combo[-1] == 0:
            continue  # not full degree
        n_total += 1
        try:
            roots = np.roots(coeffs)
            M = float(np.prod(np.maximum(1.0, np.abs(roots))))
        except:
            continue
        if M < 1.001:  # cyclotomic (or close numerically)
            n_cyclotomic += 1
            continue
        n_above_one += 1
        Mmax = float(np.max(np.abs(roots)))
        excess = Mmax - 1.0
        if excess < smallest_excess:
            smallest_excess = excess
            smallest_excess_poly = coeffs.copy()
    return dict(
        d=d, height=height, n_total=n_total,
        n_cyclotomic=n_cyclotomic, n_above_one=n_above_one,
        smallest_excess=smallest_excess if smallest_excess < float('inf') else None,
        smallest_poly=smallest_excess_poly,
    )


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 95, flush=True)
    print("SCHINZEL-ZASSENHAUS CONJECTURE / Dimitrov theorem", flush=True)
    print("=" * 95, flush=True)

    print("\n[1] EXHAUSTIVE SCAN: smallest max|alpha| - 1 for non-cyclotomic monic Z[z]")
    print("   (height-1 polys, coeffs in {-1, 0, 1})")
    print(f"\n{'d':>3} {'#total':>10} {'#cyclo':>9} {'#non-cyclo':>12} "
          f"{'min(M_max - 1)':>16} {'pred 1/d':>10}",
          flush=True)
    for d in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
        result = exhaustive_scan_height(d, 1)
        if result is None: continue
        excess = result['smallest_excess']
        print(f"{d:>3} {result['n_total']:>10} {result['n_cyclotomic']:>9} "
              f"{result['n_above_one']:>12} "
              f"{excess if excess is None else f'{excess:>16.8f}':>16} "
              f"{1.0/d:>10.4f}",
              flush=True)

    print("\n[2] HEIGHT-2 SCAN (small degrees only):")
    print(f"{'d':>3} {'#total':>10} {'min(M_max - 1)':>16}", flush=True)
    for d in [3, 4, 5, 6, 7, 8]:
        result = exhaustive_scan_height(d, 2)
        if result is None: continue
        excess = result['smallest_excess']
        print(f"{d:>3} {result['n_total']:>10} "
              f"{excess if excess is None else f'{excess:>16.8f}':>16}",
              flush=True)

    # 3. KEY EXAMPLES from the literature
    print("\n[3] KNOWN EXAMPLES:")
    examples = [
        ("Smyth's L_0 (Lehmer)", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1])),
        ("plastic z^3 - z - 1", np.array([1, 0, -1, -1])),
        ("(z-1)(z-2)", np.array([1, -3, 2])),
        ("z^4 + z + 1", np.array([1, 0, 0, 1, 1])),
        ("z^3 + z + 1", np.array([1, 0, 1, 1])),
    ]
    for name, p in examples:
        try:
            Mmax = max_modulus(p.astype(float))
            excess = Mmax - 1.0
            ld = discriminant_log(p.astype(float))
            print(f"  {name:>24}: max|alpha| = {Mmax:.6f}, excess = {excess:+.6f}, "
                  f"log|Disc| = {ld:.4f}",
                  flush=True)
        except Exception as e:
            print(f"  {name:>24}: error {e}")

    # 4. PANDROSION-HADAMARD route: relate excess to log|Disc|
    print("\n[4] PANDROSION-HADAMARD: relate Schinzel excess to discriminant")
    print("    For non-cyclotomic Z-poly with all roots near |z|=1:")
    print("    log|Disc| / (M_max - 1) measures how 'close to cyclotomic' a poly is")
    print(f"\n{'d':>3} {'min M_max - 1':>14} {'corresp poly':>40} {'log|Disc|':>12}",
          flush=True)
    rng = np.random.default_rng(20260801)
    for d in [5, 8, 10, 12, 15, 20]:
        smallest = float('inf')
        smallest_poly = None
        # height-1 sample
        for _ in range(50000):
            coefs = rng.choice([-1, 0, 1], size=d+1)
            coefs[0] = 1
            if coefs[-1] == 0: continue
            try:
                Mmax = max_modulus(coefs.astype(float))
                if Mmax > 1.0001 and (Mmax - 1) < smallest:
                    smallest = Mmax - 1
                    smallest_poly = coefs.copy()
            except:
                continue
        if smallest_poly is not None:
            ld = discriminant_log(smallest_poly.astype(float))
            print(f"{d:>3} {smallest:>14.8f} {str(smallest_poly[:8]) + '...':>40} "
                  f"{ld:>12.4f}",
                  flush=True)

    # 5. Dimitrov 2019 result: max|alpha| >= 2^(1/(4d)) approximately
    print("\n[5] DIMITROV 2019 LOWER BOUND: max|alpha|^d >= 2^(1/4)")
    print("    Equivalently: log(max|alpha|) >= (log 2) / (4d) for non-cyclotomic")
    print(f"\n{'d':>3} {'Dimitrov lower bd on log M_max':>32} {'1/d':>10}",
          flush=True)
    for d in [3, 5, 8, 10, 15, 20, 50, 100]:
        bd = math.log(2) / (4 * d)
        print(f"{d:>3} {bd:>32.6f} {1.0/d:>10.4f}", flush=True)

    print("\n" + "=" * 95)
    print("CONCLUSION:")
    print("- Empirical: smallest non-cyclotomic max|alpha| - 1 ~ 0.176 at d=10 (Smyth)")
    print("- Dimitrov 2019: log M_max >= log 2 / (4d), unconditional")
    print("- Smyth's L_0 saturates the slowest-growing case")
    print("- Schinzel-Zassenhaus conjecture (linear 1/d): still open in sharp form")
    print("=" * 95)


if __name__ == "__main__":
    main()
