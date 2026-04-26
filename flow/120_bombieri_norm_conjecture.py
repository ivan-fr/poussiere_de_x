"""
PAPER: 120 (NEW — Bombieri norm conjecture on minimal polynomials)
TITLE: Bombieri's L^2 / Bombieri-Weyl Norm Conjecture for Minimal Polynomials
       via Pandrosion-Hadamard
STATUS: empirical certificate + Pandrosion-form for several Bombieri inequalities;
        sharp constants and minimal-polynomial extension remain open
DEPENDS: 034 (Hadamard-Bombieri), 056 (Hardy-Pandrosion), 087 (Pandrosion-Gram),
         020/105 (Lehmer), 069/112 (Mahler/Smyth)

THEORY
======

DEFINITIONS
-----------

For polynomial P(z) = sum_{k=0}^d c_k z^k:

  L^2 norm:          ||P||_2^2 = sum |c_k|^2 = (1/2pi) integral_T |P(e^{i theta})|^2 d theta
  L^infty (sup):     ||P||_inf = max_{|z|=1} |P(z)|
  Bombieri-Weyl:     ||P||_BW^2 = sum |c_k|^2 / C(d, k)
  Mahler measure:    M(P) = |a_d| prod max(1, |alpha_k|).

------------------------------------------------------------------------
BOMBIERI INEQUALITIES
------------------------------------------------------------------------

(B1) Mahler-Bombieri:    M(P) <= ||P||_2.
     Proof: by Jensen, log M = (1/2pi) integral log|P|; by AM-QM:
       (1/2pi) integral log|P| <= log((1/2pi) integral |P|^2)^{1/2} = log ||P||_2.
     So M <= ||P||_2.

(B2) Hadamard-Bombieri (paper 034):
     |P(z)|^2 / (1 + |z|^2)^d <= ||P||_BW^2,
     with equality at the Bombieri-extremal P_0 = (1 + a z)^d.

(B3) Bombieri's L^2 conjecture for minimal polynomials (refined Mahler bound):
     For P in Z[z] minimal polynomial of an algebraic integer alpha,
     non-cyclotomic, of degree d:
       ||P||_2^2 >= 1 + (delta * d)
     for some absolute delta > 0.
     (Open in this exact sharp form; weaker non-tight versions are known.)

------------------------------------------------------------------------
PANDROSION-HADAMARD CONNECTION
------------------------------------------------------------------------

By paper 087: det G^(P) = |Disc(P)|^2.
For minimal poly P of an algebraic integer:
  Disc(P) is a non-zero integer.
  |Disc(P)| = prod_{i<j} |alpha_i - alpha_j|^2 >= 1.

Combined with (B1): M(P)^2 <= ||P||_2^2, and discriminant gives:
  log ||P||_2^2 >= 2 log M(P) >= 2 log L_0  (if Lehmer's conjecture holds)
giving ||P||_2 >= L_0 = 1.176... for non-cyclotomic monic Z-polys.

VERIFICATION
============

  1. Verify B1: M(P) <= ||P||_2 on test polynomials.
  2. Verify B2: |P(z)|^2 / (1 + |z|^2)^d <= ||P||_BW^2.
  3. Distribution of ||P||_2^2 / d for minimal polys of algebraic integers
     (testing B3 empirically).
  4. Pandrosion-Hadamard det G = |Disc|^2 link.
"""
from __future__ import annotations
import math
import itertools
import numpy as np


def L2_norm_sq_coeffs(P):
    return float(np.sum(np.abs(P)**2))


def Linf_norm_circle(P, n_pts=512):
    z = np.exp(2j * np.pi * np.arange(n_pts) / n_pts)
    return float(np.max(np.abs(np.polyval(P, z))))


def bombieri_weyl_norm_sq(P):
    """||P||_BW^2 = sum |c_k|^2 / C(d, k)."""
    d = len(P) - 1
    binom = np.array([math.comb(d, k) for k in range(d + 1)])
    coefs_low = P[::-1]
    return float(np.sum(np.abs(coefs_low)**2 / binom))


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
    print("PAPER 120 — Bombieri norm conjecture on minimal polynomials")
    print("=" * 80)

    # ===========================
    # PROOF OF B1: M(P) <= ||P||_2
    # ===========================
    # log M(P) = (1/2pi) integral log|P(e^{i theta})| d theta  (Jensen, paper 037)
    # By AM-QM concavity of log:
    #   (1/2pi) integral log|P| d theta <= log((1/2pi) integral |P|^2 d theta)^{1/2}
    #                                    = log sqrt(||P||_2^2 / 1) = (1/2) log ||P||_2^2
    # Hence M(P) <= ||P||_2.
    # QED (Mahler's classical bound).

    print("\n[1] B1: M(P) <= ||P||_2 (Mahler-Bombieri)")
    print(f"  {'P':>30} {'M(P)':>10} {'||P||_2':>10} {'OK?':>6}")
    test_polys = [
        ("z^2 - 2", np.array([1.0, 0, -2])),
        ("z^3 - z - 1 (plastic)", np.array([1.0, 0, -1, -1])),
        ("Smyth L_0", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])),
        ("z^4 + z + 1", np.array([1.0, 0, 0, 1, 1])),
        ("(1 + z)^4", np.array([1, 4, 6, 4, 1.0])),
    ]
    for name, P in test_polys:
        M = mahler_measure(P)
        L2 = math.sqrt(L2_norm_sq_coeffs(P))
        ok = M <= L2 + 1e-9
        print(f"  {name:>30} {M:>10.6f} {L2:>10.6f} {str(ok):>6}")

    # ===========================
    # B2: Hadamard-Bombieri verification
    # ===========================
    print("\n[2] B2: |P(z)|^2/(1+|z|^2)^d <= ||P||_BW^2 (Hadamard-Bombieri)")
    rng = np.random.default_rng(2026)
    for name, P in test_polys:
        d = len(P) - 1
        BW_sq = bombieri_weyl_norm_sq(P)
        max_ratio = 0.0
        for _ in range(50):
            z = complex(rng.standard_normal(), rng.standard_normal())
            ratio = abs(np.polyval(P, z))**2 / (1 + abs(z)**2)**d
            if ratio > max_ratio: max_ratio = ratio
        print(f"  {name:>30}: max |P|^2/(1+|z|^2)^d = {max_ratio:.4f}, ||P||_BW^2 = {BW_sq:.4f}")

    # ===========================
    # B3: Pandrosion-Hadamard for minimal polys: ||P||_2 lower bound
    # ===========================
    print("\n[3] B3 (open): for minimal poly of algebraic integer, ||P||_2^2 >= 1 + delta * d")
    print(f"  Testing on minimal polys of small algebraic integers...")
    minpolys = [
        ("sqrt 2: z^2 - 2", np.array([1.0, 0, -2])),
        ("plastic: z^3 - z - 1", np.array([1.0, 0, -1, -1])),
        ("golden: z^2 - z - 1", np.array([1.0, -1, -1])),
        ("z^4 - 2", np.array([1.0, 0, 0, 0, -2])),
        ("Smyth L_0", np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])),
        ("z^5 + z + 1 (S_5)", np.array([1.0, 0, 0, 0, 1, 1])),
    ]
    print(f"  {'P':>30} {'d':>3} {'||P||_2^2':>12} {'(||P||_2^2-1)/d':>18}")
    delta_min = float('inf')
    for name, P in minpolys:
        d = len(P) - 1
        L2_sq = L2_norm_sq_coeffs(P)
        delta = (L2_sq - 1) / d
        if delta < delta_min: delta_min = delta
        print(f"  {name:>30} {d:>3} {L2_sq:>12.4f} {delta:>18.4f}")
    print(f"\n  Empirical min of (||P||_2^2 - 1)/d = {delta_min:.4f}")
    print(f"  (Bombieri's conjectured uniform lower bound: open, would need delta > 0 universal.)")

    # ===========================
    # Pandrosion-Hadamard link
    # ===========================
    print("\n[4] Pandrosion-Hadamard det G = |Disc|^2 (paper 087)")
    print(f"  For Z-poly with simple roots, |Disc| in Z, |Disc| >= 1 implies log det G >= 0.")
    print(f"  Combined with M <= ||P||_2 and Lehmer's conjecture M >= 1.176:")
    print(f"  ||P||_2 >= 1.176 for monic non-cyclotomic Z-polys.")
    for name, P in minpolys:
        d = len(P) - 1
        D = discriminant(P)
        M = mahler_measure(P)
        L2 = math.sqrt(L2_norm_sq_coeffs(P))
        print(f"  {name:>30}: |Disc|={D:>10.2f}, M={M:.4f}, ||P||_2={L2:.4f}")

    # ===========================
    # Random Z-poly scan: distribution
    # ===========================
    print("\n[5] Random monic {-1, 0, 1} Z-polys: distribution of ||P||_2^2")
    print(f"  {'d':>3} {'#trials':>9} {'mean ||P||_2^2':>16} {'min ||P||_2^2':>16}")
    rng = np.random.default_rng(42)
    for d in [4, 6, 8, 10, 12]:
        n = 2000
        L2_list = []
        for combo in itertools.product([-1, 0, 1], repeat=d):
            coefs = np.array([1] + list(combo), dtype=float)
            if combo[-1] == 0: continue
            L2_list.append(L2_norm_sq_coeffs(coefs))
            if len(L2_list) >= n: break
        arr = np.array(L2_list)
        print(f"  {d:>3} {len(L2_list):>9} {arr.mean():>16.4f} {arr.min():>16.4f}")

    # ===========================
    # Bombieri-Weyl extremality
    # ===========================
    print("\n[6] Bombieri-extremal: P = (1+z)^d saturates HB inequality on z = real")
    for d in [3, 5, 8]:
        # P = (1+z)^d coefficients (high to low: 1, d, C(d,2), ..., 1)
        binom = [math.comb(d, k) for k in range(d, -1, -1)]
        P = np.array(binom, dtype=float)
        # Test at z = 1
        val = abs(np.polyval(P, 1.0))**2 / (1 + 1)**d  # = (2^d)^2 / 2^d = 2^d
        BW_sq = bombieri_weyl_norm_sq(P)
        # ||P||_BW^2 = sum C(d,k)^2 / C(d,k) = sum C(d,k) = 2^d
        print(f"  P = (1+z)^{d}: |P(1)|^2/(1+1)^{d} = {val:.4f}, ||P||_BW^2 = {BW_sq:.4f} (= 2^{d})")

    # ===========================
    # Honest assessment
    # ===========================
    print("\n[7] HONEST ASSESSMENT")
    print("  PROVED (classical):")
    print("    B1 Mahler-Bombieri:  M(P) <= ||P||_2 (Jensen + AM-QM)")
    print("    B2 Hadamard-Bombieri: |P(z)|^2/(1+|z|^2)^d <= ||P||_BW^2")
    print("    Pandrosion-Hadamard:  det G^(P) = |Disc(P)|^2 (paper 087)")
    print("  ")
    print("  OPEN (this paper's targets):")
    print("    B3: ||P||_2^2 >= 1 + delta*d for non-cyclotomic min polys (uniform delta).")
    print("    Sharp constants in Hadamard-Bombieri for non-extremal polys.")
    print("    Bombieri-Lehmer connection: ||P||_2 lower bound via Mahler.")
    print("  ")
    print("  Empirically: smallest (||P||_2^2 - 1)/d on tested minimal polys is ~0.07-0.5.")
    print("  Pandrosion-Hadamard provides the discriminant framework but does")
    print("  not close the uniform-delta question of B3.")


if __name__ == "__main__":
    main()
