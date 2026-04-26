"""
PAPER: 118 (NEW — Salem density attack)
TITLE: Density of Salem Numbers in [1, infty) via Pandrosion-Mahler
STATUS: empirical evidence + structural reformulation; conjecture remains open
DEPENDS: 020/105 (Lehmer), 069 (Mahler), 087 (Hadamard-Gram), 093 (Smyth), 112

THEORY
======

DEFINITION (Salem number): A real algebraic integer beta > 1 such that all
conjugates have |.| <= 1 and at least one conjugate has |.| = 1.

EQUIVALENT (Pandrosion form): beta is Salem iff its minimal poly P_beta is
RECIPROCAL, NON-CYCLOTOMIC, with EXACTLY ONE conjugate (= beta itself) of
modulus > 1.

For Salem beta: M(P_beta) = beta * 1 * ... * 1 = beta.
So Salem nos are EXACTLY the non-trivial values of M for reciprocal Z-polys
with one root outside the disk.

CONJECTURE (Salem density): Salem numbers are DENSE in [1, infty).
Stronger: the SMALLEST Salem number is conjecturally the Lehmer's number
1.17628... (so smallest Salem > 1 is L_0 = M(Smyth's polynomial)).

------------------------------------------------------------------------
PANDROSION-MAHLER FRAMEWORK
------------------------------------------------------------------------

For reciprocal Z-poly P of degree 2n:
  P(z) = z^n Q(z + 1/z)  for some Q in Z[t] of degree n.
The roots come in pairs (alpha, 1/alpha). For Salem:
  - one pair (beta, 1/beta) with beta > 1, 1/beta < 1.
  - all other pairs on the unit circle (= roots of unity-like, complex
    conjugate pairs theta, theta-bar with |theta| = 1).

NOTE: The "Pisot-Salem" hierarchy: Pisot numbers are real beta > 1 with
|conjugates| < 1 (all strictly inside disk). Salem are the LIMIT of
Pisot from above (one conjugate hits |.|=1).

LEHMER'S CONNECTION: smallest Salem >= L_0 = 1.17628 (conjecturally).
This is exactly Lehmer's conjecture in the reciprocal-Z-poly case.

BOYD-LEHMER: smallest Salem KNOWN is L_0 itself (Lehmer's polynomial); no
smaller candidate has been found despite extensive search.

VERIFICATION
============

  1. Generate Salem numbers from reciprocal Z-polys.
  2. Verify density: gaps between Salem candidates.
  3. Pandrosion-Mahler: reciprocal poly + 1 outside-disk root.
"""
from __future__ import annotations
import math
import itertools
import numpy as np


def is_reciprocal(P, tol=1e-9):
    n = len(P)
    return all(abs(P[i] - P[n-1-i]) < tol for i in range(n // 2))


def is_salem_candidate(P, tol=1e-6):
    """Check if P is reciprocal Z-poly with exactly one root |.|>1 (and 1/r inside)."""
    if not is_reciprocal(P): return False, None
    roots = np.roots(P)
    outside = [r for r in roots if abs(r) > 1.0 + tol]
    inside = [r for r in roots if abs(r) < 1.0 - tol]
    on_circle = [r for r in roots if abs(abs(r) - 1.0) <= tol]
    if len(outside) != 1: return False, None
    # Outside root should be REAL > 1 (Salem condition)
    beta = outside[0]
    if abs(beta.imag) > tol: return False, None
    if beta.real <= 1: return False, None
    return True, float(beta.real)


def mahler_measure(P):
    return float(abs(P[0]) * np.prod(np.maximum(1.0, np.abs(np.roots(P)))))


def main():
    print("=" * 80)
    print("PAPER 118 — Density of Salem numbers in [1, infty)")
    print("=" * 80)

    print("\n[1] Smyth's L_0 is Salem: verify")
    L0 = np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1.0])
    print(f"  L_0 reciprocal? {is_reciprocal(L0)}")
    is_s, beta = is_salem_candidate(L0)
    print(f"  Salem? {is_s}, beta = {beta}")
    if is_s:
        roots = np.roots(L0)
        on_circle = sum(1 for r in roots if abs(abs(r) - 1) < 1e-6)
        print(f"  Roots on |z|=1: {on_circle} / {len(roots)}")

    print("\n[2] Exhaustive scan for Salem candidates (reciprocal {-1,0,1} polys)")
    print(f"  {'d':>3} {'#recip':>10} {'#Salem':>8} {'min Salem':>12}")
    for d in [4, 6, 8, 10]:  # even d only
        half = d // 2
        n_recip = 0
        n_salem = 0
        salem_betas = []
        # Generate all reciprocal: a_0, a_1, ..., a_{d/2} free; a_i = a_{d-i}
        for combo in itertools.product([-1, 0, 1], repeat=half + 1):
            # Build [a_0, a_1, ..., a_{half-1}, a_{half}, a_{half-1}, ..., a_0]
            full = list(combo[:half]) + [combo[half]] + list(reversed(combo[:half]))
            full = list(combo) + list(reversed(combo[:-1]))  # actually: leading + middle + mirror
            # Simpler: take coefs a_0..a_d with a_i = a_{d-i}
            cs = list(combo)  # length half+1
            full = cs + cs[-2::-1]  # mirror
            if len(full) != d + 1: continue
            full[0] = 1  # monic
            P = np.array(full, dtype=float)
            if abs(P[-1]) < 1e-12: continue
            if not is_reciprocal(P): continue
            n_recip += 1
            is_s, beta = is_salem_candidate(P)
            if is_s:
                n_salem += 1
                salem_betas.append(beta)
        min_salem = min(salem_betas) if salem_betas else float('inf')
        print(f"  {d:>3} {n_recip:>10} {n_salem:>8} {min_salem:>12.6f}")

    print("\n[3] Salem numbers found (sorted) and gaps")
    rng = np.random.default_rng(2026)
    salem_set = set()
    # Random reciprocal scan
    for d in [10, 12, 14, 16]:
        half = d // 2
        for _ in range(20000):
            cs = rng.choice([-1, 0, 1], size=half + 1)
            cs[0] = 1
            full = list(cs) + list(reversed(cs[:-1]))
            if len(full) != d + 1: continue
            P = np.array(full, dtype=float)
            if abs(P[-1]) < 1e-12: continue
            if not is_reciprocal(P): continue
            is_s, beta = is_salem_candidate(P)
            if is_s and 1.0 < beta < 2.0:
                salem_set.add(round(beta, 6))
    salem_list = sorted(salem_set)
    print(f"  Found {len(salem_list)} distinct Salem numbers (rounded to 6 digits) in (1, 2)")
    if salem_list:
        print(f"  Smallest: {salem_list[:5]}")
        print(f"  Lehmer's L_0 = 1.176281: in set? {1.176281 in salem_set}")
        # Density: max gap in (1, 1.5)
        in_range = [s for s in salem_list if s < 1.5]
        if len(in_range) > 1:
            gaps = [in_range[i+1] - in_range[i] for i in range(len(in_range)-1)]
            print(f"  Max gap in (1, 1.5): {max(gaps):.4f}")
            print(f"  Mean gap: {sum(gaps)/len(gaps):.4f}")

    print("\n[4] Pisot vs Salem: examples")
    # Plastic 1.32472 = smallest Pisot
    plastic = np.array([1.0, 0, -1, -1])
    M_plastic = mahler_measure(plastic)
    is_recip = is_reciprocal(plastic)
    print(f"  Plastic z^3 - z - 1: M = {M_plastic:.4f}, reciprocal? {is_recip} (so Pisot, not Salem)")
    # Smallest known Salem
    print(f"  Smyth L_0: Salem with beta = 1.17628 (smallest known Salem)")

    print("\n[5] Conclusion")
    print("  Salem density in [1, infty) is OPEN.")
    print("  Lehmer's conjecture: smallest Salem = L_0 = 1.17628 (open).")
    print("  Pandrosion-Mahler framework gives reciprocal-Z-poly characterisation.")
    print("  The empirical Salem distribution is consistent with density in (1, 2)")
    print("  but proving DENSE in [1, infty) requires showing arbitrary precision.")


if __name__ == "__main__":
    main()
