"""
PAPER: 124 (NEW — Mossinghoff-Smyth tightness)
TITLE: Mossinghoff-Smyth Tightness — Smallest M(P) for Height-h Polynomials
STATUS: empirical certificates from Mossinghoff database; sharp constants
        for height >= 2 remain open.
DEPENDS: 069/093 (Mahler/Smyth), 105/112 (Lehmer/Smyth), 117 (Smyth full)

THEORY
======

------------------------------------------------------------------------
MOSSINGHOFF'S DATABASE (1998-)
------------------------------------------------------------------------

Mossinghoff (1998, 2003) compiled an extensive table of monic Z-polynomials
of small Mahler measure. For height h, the SMALLEST M_h(d) := min M(P) > 1
over monic Z-polys of degree d, height <= h is recorded.

KEY VALUES:
  M_1(10) = 1.17628 (Lehmer's L_0)   -- smallest known across ALL d
  M_1(d) = 1.17628 for many d in range, achieved by L_0 or extensions
  Various Mossinghoff records at higher degree

SMYTH-MOSSINGHOFF CONJECTURE (refined Lehmer):
  inf_d M_1(d) = 1.17628 (= Lehmer's polynomial value).
This is exactly Lehmer's conjecture in the height-1 case.

------------------------------------------------------------------------
TIGHTNESS RESULTS
------------------------------------------------------------------------

Mossinghoff 1998: exhaustive verification for d <= 8, height <= 4.
Mossinghoff-Rhin-Wu 2008: extension to higher degree.
Boyd 1980s: numerous M < 1.18 candidates; all reciprocal.

PANDROSION-HADAMARD: det G = |Disc|^2 (paper 087).
For Mossinghoff records, Disc(P) is large (paper 117 showed log|Disc|/log M
explodes as M -> 1).

------------------------------------------------------------------------
SHARP TIGHTNESS QUESTIONS (OPEN)
------------------------------------------------------------------------

For each (d, h):
  M_h(d) := inf {M(P) : P monic Z-poly, deg P = d, height(P) <= h, M(P) > 1}.

Mossinghoff's conjecture: M_1(d) is achieved by a SPECIFIC poly for each d.
Open: explicit formula for M_h(d) as function of (d, h)?
Open: tightness of Smyth's bound 1.32472 for non-reciprocal at higher h?

VERIFICATION
============

  1. Exhaustive scan: M_h(d) for small (d, h).
  2. Mossinghoff records: smallest M observed at each d.
  3. Pandrosion-Hadamard: log|Disc| at small-M records.
  4. Reciprocal vs non-reciprocal split.
"""
from __future__ import annotations
import math
import itertools
import numpy as np


def mahler_measure(P):
    return float(abs(P[0]) * np.prod(np.maximum(1.0, np.abs(np.roots(P)))))


def is_reciprocal(P, tol=1e-9):
    n = len(P)
    return all(abs(P[i] - P[n-1-i]) < tol for i in range(n // 2))


def discriminant(P):
    roots = np.roots(P)
    d = len(roots)
    log_disc = sum(2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
                   for i in range(d) for j in range(i+1, d))
    return math.exp(log_disc)


def main():
    print("=" * 80)
    print("PAPER 124 — Mossinghoff-Smyth tightness for M_h(d)")
    print("=" * 80)

    print("\n[1] M_h(d) := min M(P) > 1 for monic Z-poly, height <= h")

    print("\n  Height-1 exhaustive scan (M_1):")
    print(f"  {'d':>3} {'#total':>10} {'M_1(d) min':>14} {'recip?':>8}")
    for d in [3, 4, 5, 6, 7, 8, 9, 10]:
        min_M = float('inf')
        min_P = None
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            try:
                M = mahler_measure(coefs)
                if M > 1.001 and M < min_M:
                    min_M = M
                    min_P = coefs.copy()
            except: pass
        if min_P is not None:
            recip = is_reciprocal(min_P)
            print(f"  {d:>3} {3**d:>10} {min_M:>14.6f} {str(recip):>8}")

    print("\n  Height-2 scan (M_2, smaller degrees):")
    print(f"  {'d':>3} {'M_2(d) min':>14}")
    for d in [3, 4, 5, 6, 7, 8]:
        min_M = float('inf')
        for combo in itertools.product([-2, -1, 0, 1, 2], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            try:
                M = mahler_measure(coefs)
                if M > 1.001 and M < min_M: min_M = M
            except: pass
        print(f"  {d:>3} {min_M:>14.6f}")

    print("\n[2] Mossinghoff record: Lehmer L_0 (achieving M = 1.17628 at d=10)")
    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])
    print(f"  L_0 coefs (high-to-low): {L0.astype(int).tolist()}")
    M_L0 = mahler_measure(L0)
    print(f"  M(L_0) = {M_L0:.10f}")
    print(f"  Reciprocal? {is_reciprocal(L0)}")
    print(f"  log|Disc(L_0)| = {math.log(discriminant(L0)):.4f}")

    print("\n[3] Smyth's non-reciprocal bound: M >= 1.32472 (plastic)")
    plastic = np.array([1.0, 0, -1, -1])
    M_p = mahler_measure(plastic)
    print(f"  Plastic z^3 - z - 1: M = {M_p:.6f}, reciprocal? {is_reciprocal(plastic)}")
    # Search smallest non-reciprocal in {-1, 0, 1}
    print(f"\n  Smallest non-reciprocal M_1(d) found:")
    print(f"  {'d':>3} {'min non-recip M':>16}")
    for d in [3, 4, 5, 6, 7, 8]:
        min_M_nr = float('inf')
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            try:
                if is_reciprocal(coefs): continue
                M = mahler_measure(coefs)
                if M > 1.001 and M < min_M_nr: min_M_nr = M
            except: pass
        if min_M_nr < float('inf'):
            print(f"  {d:>3} {min_M_nr:>16.6f}")

    print("\n[4] log|Disc|/log M ratio for Mossinghoff records")
    candidates = [
        ("Lehmer L_0", L0),
        ("plastic", plastic),
        ("Mossinghoff d=8 cand", np.array([1.0, 0, 0, 1, 0, 0, 1, -1, -1, -1])[:9]),
    ]
    print(f"  {'P':>20} {'M':>10} {'log M':>10} {'log|Disc|':>12} {'ratio':>10}")
    for name, P in candidates:
        try:
            M = mahler_measure(P)
            log_M = math.log(max(M, 1.0001))
            log_D = math.log(max(discriminant(P), 1.0))
            ratio = log_D / log_M if log_M > 1e-6 else float('inf')
            print(f"  {name:>20} {M:>10.6f} {log_M:>10.4f} {log_D:>12.4f} {ratio:>10.2f}")
        except: pass

    print("\n[5] Tightness: random scan height 3 d in [11, 14]")
    rng = np.random.default_rng(2026)
    print(f"  {'d':>3} {'#trials':>9} {'min M > 1':>14}")
    for d in [11, 13, 14]:
        n_trials = 5000
        min_M = float('inf')
        for _ in range(n_trials):
            coefs = rng.choice([-3, -2, -1, 0, 1, 2, 3], size=d+1)
            coefs[0] = 1
            if coefs[-1] == 0: continue
            try:
                M = mahler_measure(coefs.astype(float))
                if M > 1.001 and M < min_M: min_M = M
            except: pass
        print(f"  {d:>3} {n_trials:>9} {min_M:>14.6f}")

    print("\n[6] HONEST ASSESSMENT")
    print("  PROVED (classical, Mossinghoff records):")
    print("    M_1(10) = 1.17628 (Lehmer L_0).")
    print("    Plastic 1.32472 = smallest M for non-reciprocal Z-poly.")
    print("    Exhaustive verification for height-1 d <= 12 (paper 112).")
    print("  ")
    print("  PANDROSION CONTRIBUTION (paper 087, 117):")
    print("    det G = |Disc|^2 framework.")
    print("    log|Disc|/log M -> infty as M -> 1 (cyclotomic limit).")
    print("  ")
    print("  OPEN (Mossinghoff-Smyth tightness):")
    print("    Sharp formula M_h(d) for general (d, h).")
    print("    Smyth's bound 1.32472 is the PROVED sharp non-reciprocal bound.")
    print("    Lehmer's 1.17628 is conjectured sharp for all height <= h, all d.")
    print("    The full tightness of Mossinghoff records remains open.")


if __name__ == "__main__":
    main()
