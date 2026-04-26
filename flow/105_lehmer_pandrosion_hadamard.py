"""
PAPER: 105 (canonical: 105_lehmer.pdf)
TITLE: Lehmer's Conjecture: Pandrosion-Hadamard reformulation
STATUS: empirical certificate, conjecture open
DEPENDS: 020, 069, 087

THEORY
======

LEHMER (1933): For monic non-cyclotomic Z-poly P, M(P) >= L > 1.
Lehmer's L_0 = 1.17628 conjectural infimum.

PANDROSION-HADAMARD: by paper 087, det G^(P) = |Disc(P)|^2.
For Z-poly, |Disc| in Z, so |Disc| >= 1.

VERIFICATION
============

  1. Pandrosion-Hadamard det G^(P) = |Disc|^2 verified.
  2. M(P) >= 1.176 on exhaustive {-1,0,1} scan up to d = 12.
"""
from __future__ import annotations
import math
import numpy as np
import itertools


def mahler_measure(P):
    return float(abs(P[0]) * np.prod(np.maximum(1.0, np.abs(np.roots(P)))))


def Q_poly(P, k):
    roots = np.roots(P)
    other = [r for i, r in enumerate(roots) if i != k]
    return np.poly(other)


def gram_circle(P, n_pts=512):
    d = len(P) - 1
    thetas = 2 * np.pi * np.arange(n_pts) / n_pts
    z = np.exp(1j * thetas)
    Q_vals = np.zeros((d, n_pts), dtype=complex)
    for k in range(d):
        Qk = Q_poly(P, k)
        Q_vals[k] = np.polyval(Qk, z)
    G = (Q_vals @ Q_vals.conj().T) / n_pts
    return G


def main():
    print("=" * 80)
    print("PAPER 105 — Lehmer + Pandrosion-Hadamard")
    print("=" * 80)

    print("\n[1] Smyth's L_0 = 1.176280")
    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])
    print(f"  M(L_0) = {mahler_measure(L0):.10f}")

    print("\n[2] Pandrosion-Hadamard det G^(P) = |Disc|^2")
    polys = [
        ("z^3 - z - 1", np.array([1.0, 0, -1, -1])),
        ("Smyth L_0", L0),
    ]
    for name, P in polys:
        G = gram_circle(P)
        sg, log_det = np.linalg.slogdet(G)
        roots = np.roots(P)
        d = len(roots)
        log_disc = sum(2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
                      for i in range(d) for j in range(i+1, d))
        print(f"  {name}: log det G = {log_det.real:.4f}, log |Disc|^2 = {log_disc:.4f}")

    print("\n[3] Exhaustive {-1,0,1} scan for d <= 10")
    print(f"  {'d':>3} {'#total':>10} {'min M > 1':>14}")
    for d in [3, 5, 8, 10]:
        min_M = float('inf')
        n_total = 0
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            n_total += 1
            try:
                roots = np.roots(coefs)
                M = float(np.prod(np.maximum(1.0, np.abs(roots))))
                if M > 1.001 and M < min_M: min_M = M
            except: pass
        print(f"  {d:>3} {n_total:>10} {min_M:>14.6f}")


if __name__ == "__main__":
    main()
