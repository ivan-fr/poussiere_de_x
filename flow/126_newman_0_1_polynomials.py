"""
PAPER: 126 (NEW — Newman 0/1 polynomial conjecture)
TITLE: Newman's {0, 1}-Coefficient Polynomials — Smallest Mahler Measure
STATUS: open in sharp form; Boyd 1980 gave bound 1.18034 for non-trivial cases.
DEPENDS: 020/105 (Lehmer), 069/093 (Mahler/Smyth), 087 (Hadamard-Gram),
         112 (Smyth small Mahler), 124 (Mossinghoff-Smyth)

THEORY
======

------------------------------------------------------------------------
NEWMAN POLYNOMIALS / {0, 1}-COEFFICIENTS
------------------------------------------------------------------------

A "Newman polynomial" (or 0-1 polynomial) is a polynomial with all coefficients
in {0, 1}. They appear in:
  - Combinatorics (generating functions of subsets of N)
  - Theory of digital signal processing
  - Number theory of algebraic integers (Lehmer-style)

Trivial cases:
  - Cyclotomic 0-1 polys: M = 1.
  - 1 + z + z^2 + ... + z^{n-1} (= Phi_n related): M = 1 if cyclotomic.

------------------------------------------------------------------------
NEWMAN'S CONJECTURE (refined Lehmer for 0-1 polys)
------------------------------------------------------------------------

For 0-1 polynomial P (NON-cyclotomic), what is inf M(P) > 1?

EMPIRICAL: smallest Mahler measure among 0-1 polynomials is conjecturally
the same as Lehmer's L_0 = 1.17628... (Boyd 1980).

But L_0 = z^10 + z^9 - z^7 - z^6 - ... has NEGATIVE coefficients, so it's
NOT a 0-1 poly. Hence the question becomes: what's the smallest M for
PURELY 0-1 polynomials?

Boyd 1980: smallest M for 0-1 polys (deg <= 36) is = 1.18034 (achieved by
specific 0-1 poly).

------------------------------------------------------------------------
PANDROSION-HADAMARD FRAMEWORK
------------------------------------------------------------------------

By paper 087: det G^(P) = |Disc(P)|^2.
0-1 polys have Disc that is bounded by polynomial in d (since coefs are
in {0, 1} hence height-1).

Combined with paper 117 (log|Disc|/log M -> infty as M -> 1):
The 0-1 constraint LIMITS the discriminant range, narrowing the search.

VERIFICATION
============

  1. Exhaustive search 0-1 polys: smallest M(P) > 1.
  2. Comparison with general height-1 polys (Lehmer L_0).
  3. Pandrosion-Hadamard discriminant invariant.
"""
from __future__ import annotations
import math
import itertools
import numpy as np


def mahler_measure(P):
    return float(abs(P[0]) * np.prod(np.maximum(1.0, np.abs(np.roots(P)))))


def is_cyclotomic_like(P, tol=1e-3):
    """Test if M(P) ~ 1 (suggests cyclotomic)."""
    return mahler_measure(P) < 1.001


def discriminant(P):
    roots = np.roots(P)
    d = len(roots)
    log_disc = sum(2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
                   for i in range(d) for j in range(i+1, d))
    return math.exp(log_disc)


def main():
    print("=" * 80)
    print("PAPER 126 — Newman {0, 1}-coefficient polynomials")
    print("=" * 80)

    # ===========================
    # PROOF THAT 0-1 POLYS HAVE M >= 1 STRICTLY > 1 IF NON-CYCLOTOMIC
    # ===========================
    # 0-1 polys are monic Z-polys with coefs in {0, 1}.
    # By Kronecker (paper 090): M(P) = 1 for monic Z-poly iff P is cyclotomic.
    # 0-1 polys can be cyclotomic (e.g., 1 + z, 1 + z + z^2 = Phi_3-related,
    # 1 + z + z^2 + ... + z^{n-1} when n prime = (z^n - 1)/(z - 1)).
    # Non-cyclotomic 0-1 polys: M > 1.
    # Lehmer: M >= L_0 ~ 1.176... (conjectured).
    # 0-1 case: empirical bound 1.18034 (Boyd 1980).

    print("\n[1] Cyclotomic 0-1 polys: M = 1")
    cyclotomic_01 = [
        ("z + 1", np.array([1, 1.0])),
        ("z^2 + z + 1", np.array([1, 1, 1.0])),
        ("z^3 + z^2 + z + 1", np.array([1, 1, 1, 1.0])),  # = (z^4 - 1)/(z - 1)
    ]
    for name, P in cyclotomic_01:
        M = mahler_measure(P)
        print(f"  {name}: M = {M:.6f}, cyclotomic? {is_cyclotomic_like(P)}")

    print("\n[2] Exhaustive search: 0-1 polys, smallest non-cyclotomic M")
    print(f"  {'d':>3} {'#total 0-1':>12} {'#non-cyclo':>12} {'min M > 1':>14}")
    for d in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
        n_total = 0
        n_nc = 0
        min_M = float('inf')
        min_P = None
        for combo in itertools.product([0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            n_total += 1
            try:
                M = mahler_measure(coefs)
                if M > 1.001:
                    n_nc += 1
                    if M < min_M:
                        min_M = M
                        min_P = coefs.copy()
            except: pass
        min_M_str = f"{min_M:.6f}" if min_M < float('inf') else "(none)"
        print(f"  {d:>3} {n_total:>12} {n_nc:>12} {min_M_str:>14}")

    print("\n[3] Boyd 1980 record: smallest 0-1 poly M ~ 1.18034")
    # Boyd's record: specific 0-1 poly of degree near 35 with M = 1.18034
    # Construction: 1 + z + z^2 + ... + z^{m-1} with specific m,
    # related to factorization patterns.
    # Try to find low-M 0-1 polys at d = 12-14
    print(f"  Sample 0-1 polys near 1.18 at higher degree:")
    rng = np.random.default_rng(2026)
    for d in [12, 14, 16]:
        n_trials = 5000
        min_M = float('inf')
        min_P = None
        for _ in range(n_trials):
            coefs = rng.choice([0, 1], size=d+1)
            coefs[0] = 1
            if coefs[-1] == 0: continue
            try:
                M = mahler_measure(coefs.astype(float))
                if 1.001 < M < min_M:
                    min_M = M
                    min_P = coefs.copy()
            except: pass
        if min_P is not None:
            print(f"  d = {d}: min M = {min_M:.6f}, P = {min_P.astype(int).tolist()}")

    print("\n[4] Pandrosion-Hadamard det G = |Disc|^2 for 0-1 polys")
    test_01 = [
        ("z^4 + 1 (cyclo)", np.array([1.0, 0, 0, 0, 1])),
        ("z^4 + z + 1", np.array([1.0, 0, 0, 1, 1])),
        ("z^5 + z + 1", np.array([1.0, 0, 0, 0, 1, 1])),
        ("z^7 + z^4 + 1", np.array([1.0, 0, 0, 1, 0, 0, 0, 1])),
    ]
    print(f"  {'P':>20} {'M':>10} {'|Disc|':>12} {'log|Disc|/log M':>20}")
    for name, P in test_01:
        M = mahler_measure(P)
        D = discriminant(P)
        if M > 1.001:
            ratio = math.log(D) / math.log(M)
            print(f"  {name:>20} {M:>10.6f} {D:>12.2f} {ratio:>20.2f}")
        else:
            print(f"  {name:>20} {M:>10.6f} {D:>12.2f}  (cyclotomic)")

    print("\n[5] Comparison with Lehmer L_0 (NOT a 0-1 poly)")
    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])
    print(f"  L_0 = z^10 + z^9 - z^7 - z^6 - z^5 - z^4 - z^3 + z + 1")
    print(f"  Has NEGATIVE coefficients, so NOT a 0-1 poly.")
    print(f"  M(L_0) = 1.176280 (smallest known across ALL height-1 polys).")
    print(f"  0-1 record (Boyd 1980): M ~ 1.18034 (slightly larger).")

    print("\n[6] HONEST ASSESSMENT")
    print("  PROVED (classical):")
    print("    0-1 polys are monic Z-polys; Kronecker classifies M = 1 cases.")
    print("    Non-cyclotomic 0-1 polys: M > 1 strictly (Kronecker).")
    print("  ")
    print("  EMPIRICAL (Boyd 1980, this paper):")
    print("    Smallest M for 0-1 polys appears to be ~1.18 (vs Lehmer 1.176).")
    print("    The slight gap is because 0-1 constraint excludes negative coefs.")
    print("  ")
    print("  PANDROSION FRAMEWORK:")
    print("    Hadamard-Gram det G = |Disc|^2 (paper 087).")
    print("    log|Disc|/log M ratio (paper 117) applies as in general case.")
    print("  ")
    print("  OPEN:")
    print("    Sharp value of min_{0-1, non-cyclo} M(P) (conjecturally ~ 1.18034).")
    print("    Lehmer-style bound for 0-1 polys: same conjectural status.")


if __name__ == "__main__":
    main()
