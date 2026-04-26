"""
PAPER: 041 (canonical: 31pandrosion_sturm.pdf)
TITLE: Sturm Sequences via Pandrosion
STATUS: proved (Sturm 1829, Pandrosion-form)
DEPENDS: 011

THEORY
======

STURM'S THEOREM: For real polynomial P with simple real roots, the Sturm
sequence (P_0 = P, P_1 = P', P_{k+1} = -rem(P_{k-1}, P_k)) gives the count
of real roots in (a, b] as V(a) - V(b), where V(x) is the number of sign
variations of (P_0(x), ..., P_m(x)).

PANDROSION CONNECTION: P_1 = P' is the Pandrosion-anchor of P at 'self';
each P_k is a Pandrosion remainder.

VERIFICATION
============

  1. Sturm sequence of test polynomial.
  2. Real-root count via sign variations.
"""
from __future__ import annotations
import numpy as np


def sturm_sequence(P):
    """Generate the Sturm sequence."""
    seq = [P.copy(), np.polyder(P).copy()]
    while len(seq[-1]) > 1:
        rem = np.polydiv(seq[-2], seq[-1])[1]
        if len(rem) == 0 or all(abs(c) < 1e-12 for c in rem): break
        seq.append(-rem)
    return seq


def sign_variations(values):
    nonzero = [v for v in values if abs(v) > 1e-12]
    sgns = [1 if v > 0 else -1 for v in nonzero]
    return sum(1 for i in range(len(sgns)-1) if sgns[i] * sgns[i+1] < 0)


def real_roots_in(P, a, b):
    """Count real roots in (a, b] via Sturm."""
    seq = sturm_sequence(P)
    Va = sign_variations([float(np.polyval(p, a)) for p in seq])
    Vb = sign_variations([float(np.polyval(p, b)) for p in seq])
    return Va - Vb


def main():
    print("=" * 80)
    print("PAPER 41 — Sturm sequences (Pandrosion-form)")
    print("=" * 80)

    polys = [
        ("(z-1)(z-2)(z-3)", np.poly([1, 2, 3.0])),
        ("z^3 - 1", np.array([1.0, 0, 0, -1])),
        ("z^4 - 5 z^2 + 4 = (z-1)(z+1)(z-2)(z+2)", np.array([1.0, 0, -5, 0, 4])),
    ]
    print("\n[1] Real-root count via Sturm")
    for name, P in polys:
        seq = sturm_sequence(P)
        n_real = real_roots_in(P, -10, 10)
        true_n = sum(1 for r in np.roots(P) if abs(r.imag) < 1e-9)
        print(f"  {name}: Sturm count = {n_real}, true = {true_n}")

    print("\n[2] Localizing roots via interval bisection (Sturm)")
    P = np.poly([1, 2, 3, 4.0])
    print(f"  P = (z-1)(z-2)(z-3)(z-4)")
    for a, b in [(0, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 5)]:
        n = real_roots_in(P, a, b)
        print(f"  ({a}, {b}]: {n} root(s)")


if __name__ == "__main__":
    main()
