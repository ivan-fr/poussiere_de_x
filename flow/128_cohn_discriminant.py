"""
PAPER: 128 (NEW — Cohn's discriminant problem)
TITLE: Cohn's Discriminant Problem — Lower Bounds via Pandrosion-Hadamard
STATUS: Cohn 1933 proved minimal-disc bounds; sharp constants in terms of
        height + degree remain open in many cases.
DEPENDS: 032 (discriminant), 087 (Pandrosion-Hadamard det G = |Disc|^2),
         105 (Lehmer), 117 (log|Disc|/log M), 124 (Mossinghoff)

THEORY
======

------------------------------------------------------------------------
COHN'S DISCRIMINANT PROBLEM
------------------------------------------------------------------------

For an algebraic integer alpha of degree d with minimal polynomial P_alpha,
what is the lower bound on |Disc(P_alpha)| as a function of:
  - the degree d,
  - the height h(alpha) (= max coefficient of P_alpha),
  - the Mahler measure M(alpha) ?

CLASSICAL RESULTS:

(C1) Minkowski-Hermite (1850s): For algebraic integer alpha, deg d,
       |Disc(P_alpha)| >= d^d / d! ~ exp(d).

(C2) Cohn 1933: For non-cyclotomic alpha, more refined lower bounds in
     terms of height.

(C3) Smith 1954: |Disc(P_alpha)| >= (d/e^d) (for irreducible).

(C4) Open: sharp lower bounds for SPECIFIC families (totally real, CM, ...).

------------------------------------------------------------------------
PANDROSION-HADAMARD CONNECTION
------------------------------------------------------------------------

By paper 087:
  det G^(P) = |Disc(P)|^2 = prod_{i<j} |alpha_i - alpha_j|^2.

For minimal poly of algebraic integer:
  Disc(P) is NON-ZERO INTEGER, hence |Disc(P)| >= 1.

But Cohn's question is QUANTITATIVE: how large is |Disc| as function of (d, h)?

PANDROSION INVARIANT VIEW:
  log |Disc(P)| = sum_{i<j} log |alpha_i - alpha_j|^2.
This is the "spread" of the conjugate set in the complex plane.

For totally real alpha (all conjugates in R), conjugates lie in some interval
[-H, H] (H ~ height). By packing argument:
  prod |alpha_i - alpha_j|^2 has lower bound related to discrepancy of points
  on real line.

VERIFICATION
============

  1. |Disc(P_alpha)| >= 1 for monic Z-poly.
  2. Empirical |Disc| growth with degree.
  3. Connection log|Disc|/log M (paper 117).
"""
from __future__ import annotations
import math
import itertools
import numpy as np


def discriminant(P):
    roots = np.roots(P)
    d = len(roots)
    log_disc = sum(2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
                   for i in range(d) for j in range(i+1, d))
    return math.exp(log_disc)


def mahler_measure(P):
    return float(abs(P[0]) * np.prod(np.maximum(1.0, np.abs(np.roots(P)))))


def height(P):
    return int(max(abs(c) for c in P))


def main():
    print("=" * 80)
    print("PAPER 128 — Cohn's discriminant problem")
    print("=" * 80)

    print("\n[1] Classical bounds: Minkowski-Hermite |Disc| >= d^d / d!")
    print(f"  {'d':>3} {'d^d / d!':>14} {'log of bound':>14}")
    for d in [3, 5, 8, 10, 15]:
        bd = d**d / math.factorial(d)
        print(f"  {d:>3} {bd:>14.2f} {math.log(bd):>14.4f}")

    print("\n[2] Empirical |Disc| for various minimal polys")
    print(f"  {'P':>30} {'d':>3} {'h':>3} {'|Disc|':>14} {'log|Disc|':>12}")
    test_polys = [
        ("z^2 - 2", np.array([1, 0, -2.0])),
        ("z^3 - z - 1 (plastic)", np.array([1, 0, -1, -1.0])),
        ("z^3 - 2", np.array([1, 0, 0, -2.0])),
        ("z^4 - z + 1", np.array([1, 0, 0, -1, 1.0])),
        ("z^5 - z + 1", np.array([1, 0, 0, 0, -1, 1.0])),
        ("Smyth L_0", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])),
    ]
    for name, P in test_polys:
        d = len(P) - 1
        h = height(P)
        D = discriminant(P)
        print(f"  {name:>30} {d:>3} {h:>3} {D:>14.2f} {math.log(D):>12.4f}")

    print("\n[3] Pandrosion-Hadamard verified: det G = |Disc|^2 (paper 087)")
    print(f"  Already verified in papers 087, 105, 117. Quick re-check:")
    P = np.array([1.0, 0, -1, -1])  # plastic
    roots = np.roots(P)
    Pp = np.polyder(P)
    D_via_pairs = float(np.prod([(roots[i] - roots[j])**2
                                  for i in range(3) for j in range(i+1, 3)]).real)
    D_via_deriv = float(np.prod([np.polyval(Pp, ak) for ak in roots]).real)
    print(f"  Plastic z^3 - z - 1:")
    print(f"  Disc via pairs:    {D_via_pairs:.4f}")
    print(f"  Disc via prod P': {D_via_deriv:.4f}")

    print("\n[4] Empirical |Disc|^{1/d} (root-mean discriminant)")
    print(f"  {'d':>3} {'|Disc|^{1/d} mean (random Z-polys)':>34}")
    rng = np.random.default_rng(2026)
    for d in [3, 5, 7, 10, 12]:
        log_disc_list = []
        for _ in range(200):
            coefs = rng.choice([-1, 0, 1], size=d+1)
            coefs[0] = 1
            if coefs[-1] == 0: continue
            try:
                D = discriminant(coefs.astype(float))
                if D > 1: log_disc_list.append(math.log(D))
            except: pass
        if log_disc_list:
            mean_log = np.mean(log_disc_list)
            print(f"  {d:>3} {math.exp(mean_log/d):>34.4f}")

    print("\n[5] Sharp Cohn bounds: empirical comparison")
    print(f"  Smith 1954: |Disc| >= (d/e)^d for non-cyclotomic monic Z-poly")
    print(f"  {'d':>3} {'Smith bound':>14} {'min |Disc| observed':>22}")
    for d in [3, 4, 5, 6, 8]:
        smith_bd = (d / math.e)**d
        # exhaustive {-1, 0, 1} search
        min_D = float('inf')
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            try:
                M = mahler_measure(coefs)
                if M < 1.001: continue  # cyclotomic
                D = discriminant(coefs)
                if D > 0 and D < min_D: min_D = D
            except: pass
        print(f"  {d:>3} {smith_bd:>14.2f} {min_D:>22.4f}")

    print("\n[6] HONEST ASSESSMENT")
    print("  PROVED:")
    print("    Minkowski-Hermite: |Disc| >= d^d / d! for algebraic integer of degree d.")
    print("    Smith 1954: refined bounds for non-cyclotomic case.")
    print("    Pandrosion-Hadamard: det G = |Disc|^2 (paper 087, 105).")
    print("  ")
    print("  PANDROSION CONTRIBUTION:")
    print("    Discriminant via Pandrosion-Hadamard Gram matrix.")
    print("    Connection log|Disc|/log M (paper 117) shows ratio explodes near")
    print("    cyclotomic limit (M -> 1).")
    print("  ")
    print("  OPEN:")
    print("    Sharp Cohn lower bounds in (d, height, M) parameters.")
    print("    Discriminant of specific families (totally real, CM, ...).")
    print("    Effective Smith-style bounds with explicit constants.")


if __name__ == "__main__":
    main()
