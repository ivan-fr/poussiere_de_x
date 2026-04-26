"""
ERDŐS-LITTLEWOOD problem.

Open since 1950s: do there exist polynomials P_n(z) of degree n with
coefficients in {+1, -1} such that
    c_1 sqrt(n) <= |P_n(e^{i theta})| <= c_2 sqrt(n)  forall theta?

Equivalently: are flat Littlewood polynomials of arbitrarily high degree
possible? Beck (1991) showed max/min ratio >= 1.196 (so c_2/c_1 >= 1.196).

Explicit constructions (Rudin-Shapiro) give max/min ratio ~ sqrt(2).

Pandrosion approach via paper 23 (Pandrosion energy):
  - L^2 norm of P_n on |z|=1: sqrt(n+1) (Parseval).
  - Sup-norm: |P_n|_∞.
  - Inf-norm: |P_n|_inf.
  - Ratio (sup/L^2) measures flatness.

Numerical strategy: search small-n {±1} polynomials for smallest sup/L^2 ratio.
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def evaluate_on_circle(coeffs, n_pts=1024):
    """|P(e^{i theta})| at n_pts uniform theta values."""
    thetas = 2 * np.pi * np.arange(n_pts) / n_pts
    z = np.exp(1j * thetas)
    # Polyval with high coeffs first
    vals = np.polyval(coeffs, z)
    return np.abs(vals)


def littlewood_stats(coeffs, n_pts=2048):
    """Return (min |P|, max |P|, L^2-norm, max/min ratio, max/L^2 ratio)."""
    abs_vals = evaluate_on_circle(coeffs, n_pts)
    pmin = float(abs_vals.min())
    pmax = float(abs_vals.max())
    L2 = math.sqrt(np.mean(abs_vals**2))  # Parseval: = sqrt(sum |a_k|^2)
    # for ±1 polys with d+1 coeffs: L^2 = sqrt(d+1)
    if pmin < 1e-12:
        ratio_max_min = float('inf')
    else:
        ratio_max_min = pmax / pmin
    ratio_max_L2 = pmax / L2 if L2 > 0 else float('inf')
    ratio_min_L2 = pmin / L2 if L2 > 0 else 0.0
    return dict(
        pmin=pmin, pmax=pmax, L2=L2,
        ratio_max_min=ratio_max_min,
        ratio_max_L2=ratio_max_L2,
        ratio_min_L2=ratio_min_L2,
    )


def rudin_shapiro_polys(n):
    """Generate Rudin-Shapiro polynomials P_n, Q_n recursively up to degree 2^n."""
    P = np.array([1, 1])
    Q = np.array([1, -1])
    for k in range(n):
        Pn = np.concatenate([P, Q])
        Qn = np.concatenate([P, -Q])
        P, Q = Pn, Qn
    return P, Q


def fekete_poly(n):
    """Fekete polynomial F_p(z) = sum_{k=1}^{p-1} (k/p) z^k for prime p.
    (Legendre-symbol coefficients, take p = n + 1 if prime.)"""
    if not is_prime(n + 1):
        return None
    p = n + 1
    coeffs = [0]  # x^0 coefficient = 0 in Fekete
    for k in range(1, p):
        coeffs.append(legendre_symbol(k, p))
    return np.array(coeffs[::-1])  # high to low order for numpy


def is_prime(n):
    if n < 2: return False
    for k in range(2, int(math.sqrt(n)) + 1):
        if n % k == 0: return False
    return True


def legendre_symbol(a, p):
    """Legendre symbol (a/p)."""
    res = pow(a, (p - 1) // 2, p)
    return res if res < p / 2 else res - p


def search_smallest_ratio(n, n_random=10000, rng=None):
    """Search for ±1 polynomial of degree n with smallest max/L^2 ratio."""
    if rng is None: rng = np.random.default_rng(20260802)
    best = float('inf')
    best_coeffs = None
    for _ in range(n_random):
        coeffs = rng.choice([1, -1], size=n+1).astype(float)
        coeffs[0] = 1  # leading coefficient = 1 by convention
        stats = littlewood_stats(coeffs)
        if stats['ratio_max_L2'] < best:
            best = stats['ratio_max_L2']
            best_coeffs = coeffs.copy()
    return best, best_coeffs


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 95, flush=True)
    print("ERDŐS-LITTLEWOOD PROBLEM: search for flat ±1 polynomials", flush=True)
    print("=" * 95, flush=True)

    print("\n[1] RUDIN-SHAPIRO polynomials: max/L^2 = sqrt(2) approx")
    print(f"{'k':>3} {'deg':>5} {'min/L2':>10} {'max/L2':>10} {'max/min':>10}",
          flush=True)
    for k in range(1, 9):
        P, Q = rudin_shapiro_polys(k)
        n = len(P) - 1
        # Reverse to high-coeff-first for np.polyval
        stats = littlewood_stats(P[::-1])
        print(f"{k:>3} {n:>5} {stats['ratio_min_L2']:>10.6f} "
              f"{stats['ratio_max_L2']:>10.6f} {stats['ratio_max_min']:>10.6f}",
              flush=True)

    print("\n[2] FEKETE polynomials (Legendre-symbol, prime p=n+1):")
    print(f"{'p (prime)':>10} {'deg n':>7} {'min/L2':>10} {'max/L2':>10}",
          flush=True)
    for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        if not is_prime(p): continue
        coeffs = fekete_poly(p - 1)
        if coeffs is None: continue
        # Coeffs in {0, +1, -1}; make ±1 only by removing leading 0
        coeffs = coeffs[coeffs != 0] if 0 in coeffs else coeffs
        # Reverse to high-coeff-first
        stats = littlewood_stats(coeffs)
        print(f"{p:>10} {len(coeffs)-1:>7} {stats['ratio_min_L2']:>10.6f} "
              f"{stats['ratio_max_L2']:>10.6f}",
              flush=True)

    print("\n[3] EXHAUSTIVE search smallest max/L^2 ratio (small n):")
    print(f"{'n':>3} {'#tested':>9} {'best max/L2':>13} {'corresponding ±1 polys':>30}",
          flush=True)
    for n in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
        best = float('inf')
        best_coefs = None
        # Iterate all ±1 polys of degree n
        for combo in itertools.product([-1, 1], repeat=n+1):
            coefs = np.array(combo, dtype=float)
            coefs[0] = 1
            stats = littlewood_stats(coefs)
            if stats['ratio_max_L2'] < best:
                best = stats['ratio_max_L2']
                best_coefs = coefs.copy()
        signs = ''.join('+' if c > 0 else '-' for c in best_coefs)
        print(f"{n:>3} {2**n:>9} {best:>13.6f} {signs:>30}", flush=True)

    print("\n[4] RANDOM SEARCH for higher n:")
    print(f"{'n':>4} {'#random':>8} {'best max/L2':>13} {'sqrt(2) ref':>13}",
          flush=True)
    for n in [10, 20, 50, 100, 200, 500]:
        rng = np.random.default_rng(20260803 + n)
        best, _ = search_smallest_ratio(n, n_random=10000, rng=rng)
        print(f"{n:>4} {10000:>8} {best:>13.6f} {math.sqrt(2):>13.6f}", flush=True)

    print("\n[5] PANDROSION L^2 trace identity: sum |a_k|^2 = (n+1) for ±1 polys")
    print("    By paper 77 Hadamard-Gram: tr(G) = sum ||Q_i||^2.")
    print("    For ±1 P, L^2-norm fixed at sqrt(n+1). Conjecture concerns sup/L^2.")
    print(f"\n  Beck 1991 lower bound on max/min: 1.196")
    print(f"  Rudin-Shapiro: sqrt(2) = 1.414")
    print(f"  Conjecture: smallest max/min approaches 1 as n -> infty? (open)")


if __name__ == "__main__":
    main()
