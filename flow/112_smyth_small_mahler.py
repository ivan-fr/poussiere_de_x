"""
PAPER: 112 (canonical: 112_smyth_small_mahler.pdf)
TITLE: Smyth's Small-Mahler Problem
STATUS: empirical certificate (exhaustive {-1,0,1} to d=12)
DEPENDS: 020, 069, 087, 093

THEORY
======

SMYTH'S CONJECTURE: For monic Z-poly P, if 1 < M(P) < L_0 = 1.17628,
then M(P) = 1 (P cyclotomic). I.e., L_0 is the conjectural infimum of
non-cyclotomic Mahler measures.

This is STRONGER than Lehmer (which only asserts existence of L > 1).

PANDROSION-HADAMARD lens: |Disc(P)| in Z, >= 1, gives discriminant lower
bound. log|Disc|/log M ratio grows as M -> 1 (cyclotomic limit).

VERIFICATION
============

  1. Exhaustive scan {-1,0,1} for d <= 12: smallest M > 1.
  2. log|Disc|/log M ratio.
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def mahler_measure(P):
    return float(abs(P[0]) * np.prod(np.maximum(1.0, np.abs(np.roots(P)))))


def discriminant(P):
    roots = np.roots(P)
    d = len(roots)
    log_disc = sum(2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
                   for i in range(d) for j in range(i+1, d))
    return math.exp(log_disc)


def main():
    print("=" * 80)
    print("PAPER 112 — Smyth's small-Mahler problem")
    print("=" * 80)

    L_0 = 1.17628081

    print(f"\n[1] Smyth's conjecture: no M(P) in (1, {L_0}) for non-cyclotomic Z-poly P")

    print("\n[2] Exhaustive scan {-1, 0, 1} for d <= 12")
    print(f"  {'d':>3} {'#total':>10} {'min M > 1':>14} {'#in (1, L_0)':>15}")
    for d in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
        min_M = float('inf')
        n_total = 0
        n_below_smyth = 0
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            n_total += 1
            try:
                roots = np.roots(coefs)
                M = float(np.prod(np.maximum(1.0, np.abs(roots))))
                if M > 1.001:
                    if M < min_M: min_M = M
                    if M < L_0 - 1e-6: n_below_smyth += 1
            except: pass
        print(f"  {d:>3} {n_total:>10} {min_M:>14.6f} {n_below_smyth:>15}")

    print("\n[3] log|Disc|/log M ratio (grows as M -> 1)")
    polys = [
        ("z^3 - z - 1 (plastic)", np.array([1, 0, -1, -1.0])),
        ("z^4 + z + 1", np.array([1, 0, 0, 1, 1.0])),
        ("Smyth L_0", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])),
    ]
    print(f"  {'P':>30} {'M':>10} {'log M':>10} {'log|Disc|':>12} {'ratio':>10}")
    for name, P in polys:
        M = mahler_measure(P)
        D = discriminant(P)
        log_M = math.log(M)
        log_D = math.log(D)
        ratio = log_D / log_M if log_M > 1e-6 else float('inf')
        print(f"  {name:>30} {M:>10.6f} {log_M:>10.4f} {log_D:>12.4f} {ratio:>10.2f}")


if __name__ == "__main__":
    main()
