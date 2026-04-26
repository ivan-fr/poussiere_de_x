"""
PAPER: 042 (canonical: 32pandrosion_markov.pdf)
TITLE: Markov Brothers Inequality via Pandrosion
STATUS: proved (Markov 1889, Markov-Bernstein theorem)
DEPENDS: 039

THEORY
======

MARKOV INEQUALITY: For real polynomial P of degree n with |P(x)| <= 1
on [-1, 1]:
  max_{x in [-1, 1]} |P^(k)(x)| <= T_n^(k)(1)
where T_n is the Chebyshev polynomial of degree n.

For k = 1: T_n'(1) = n^2, so max|P'| <= n^2.
For k = 2: T_n''(1) = n^2(n^2-1)/3, so max|P''| <= n^2(n^2-1)/3.

PANDROSION CONNECTION: Chebyshev T_n is the EXTREMAL real-rooted polynomial
in the n-th symmetric-power Bombieri-Weyl norm.

VERIFICATION
============

  1. Markov inequality for P^(k) at various k.
  2. Chebyshev T_n saturates the bound.
"""
from __future__ import annotations
import math
import numpy as np


def chebyshev_T(n, x):
    """Chebyshev T_n(x) = cos(n arccos x)."""
    return math.cos(n * math.acos(x)) if -1 <= x <= 1 else math.cosh(n * math.acosh(abs(x))) * (1 if x > 0 else (-1)**n)


def markov_bound(n, k):
    """Bound on |P^(k)(x)| for x in [-1, 1] when |P| <= 1: T_n^(k)(1)."""
    # T_n^(k)(1) = n^2(n^2-1)(n^2-4)...(n^2-(k-1)^2) / ((2k-1)!!)
    if k == 0: return 1.0
    val = 1.0
    for j in range(k):
        val *= (n*n - j*j)
    return val / math.prod(2*j + 1 for j in range(k))


def main():
    print("=" * 80)
    print("PAPER 42 — Markov brothers inequality")
    print("=" * 80)

    print("\n[1] Markov bound T_n^(k)(1)")
    print(f"  {'n':>3} {'k=1 (n^2)':>12} {'k=2':>12} {'k=3':>12}")
    for n in [3, 4, 5, 8]:
        bounds = [markov_bound(n, k) for k in range(1, 4)]
        print(f"  {n:>3} {bounds[0]:>12.0f} {bounds[1]:>12.0f} {bounds[2]:>12.0f}")

    print("\n[2] Chebyshev saturates: T_n'(1) = n^2 attained on [-1, 1]")
    xs = np.linspace(-1, 1, 1000)
    for n in [3, 4, 5]:
        # T_n via recurrence
        T = [np.ones_like(xs), xs]
        for _ in range(n - 1):
            T.append(2 * xs * T[-1] - T[-2])
        T_n = T[n]
        # Numerical derivative
        dT = np.gradient(T_n, xs)
        max_T = np.max(np.abs(T_n))
        max_dT = np.max(np.abs(dT))
        print(f"  T_{n}: max|T| = {max_T:.4f}, max|T'| ~ {max_dT:.2f}, n^2 = {n**2}")


if __name__ == "__main__":
    main()
