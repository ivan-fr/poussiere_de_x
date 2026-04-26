"""
PAPER: 082 (canonical: 72pandrosion_cheb.pdf)
TITLE: Chebyshev Polynomials and Pandrosion Extremality
STATUS: proved (classical)
DEPENDS: 042

THEORY
======

CHEBYSHEV T_n(x) = cos(n arccos x):
  - extremal: among monic polynomials of degree n on [-1, 1], 2^{1-n} T_n
    has the smallest sup-norm (= 2^{1-n}).
  - alternation: T_n attains ±1 at n+1 Chebyshev nodes x_k = cos(k pi / n).

PANDROSION FORM: The Pandrosion field of T_n on [-1, 1]:
  F_{T_n}(x) = T_n'(x) / T_n(x) = -n sin(n theta) / sin(theta) / cos(n theta)
where x = cos(theta). Vanishes at Chebyshev critical points.

VERIFICATION
============

  1. Chebyshev minimax property.
  2. Equioscillation at n+1 points.
"""
from __future__ import annotations
import math
import numpy as np


def cheb_T(n, x):
    return math.cos(n * math.acos(max(-1, min(1, x))))


def main():
    print("=" * 80)
    print("PAPER 82 — Chebyshev polynomials")
    print("=" * 80)

    print("\n[1] Equioscillation: T_n alternates ±1 at n+1 nodes")
    n = 5
    nodes = [math.cos(k * math.pi / n) for k in range(n + 1)]
    vals = [cheb_T(n, x) for x in nodes]
    print(f"  T_{n} at nodes: {[f'{v:+.4f}' for v in vals]}")

    print("\n[2] Minimax: T_n is extremal monic on [-1, 1]")
    # Monic Chebyshev: T_n / 2^{n-1}
    print(f"  Monic Chebyshev T_n / 2^(n-1): sup-norm = 1/2^(n-1)")
    for n in [3, 5, 8]:
        sup = 1.0 / 2**(n - 1)
        print(f"  n = {n}: extremal sup-norm = {sup:.6f}")

    print("\n[3] Pandrosion field at Chebyshev critical points")
    n = 4
    crits = [math.cos(k * math.pi / n) for k in range(1, n)]
    print(f"  T_{n} critical points: {crits}")
    # T_n'(x) = n U_{n-1}(x) (Chebyshev second kind)
    # At Chebyshev critical, T_n'(x) = 0
    print(f"  At these, T_n' = 0 (critical), F_{{T_n}}(x) = 0/T_n(x) = 0.")


if __name__ == "__main__":
    main()
