"""
PAPER: 039 (canonical: 29pandrosion_bernstein.pdf)
TITLE: Bernstein Inequalities via Pandrosion
STATUS: proved (Bernstein 1912, Pandrosion-form derivative bounds)
DEPENDS: 011

THEORY
======

BERNSTEIN'S INEQUALITY: For polynomial P of degree n with |P(z)| <= M on |z| = 1:
  max_{|z|=1} |P'(z)| <= n M.
Equality iff P(z) = M e^{i theta} z^n.

PANDROSION FORM: F_P = P'/P, so |F_P| = |P'|/|P|. On |z| = 1 with |P| <= M,
  max|P'| / max|P| <= n.

MARKOV BROTHERS' INEQUALITY: For real P with |P(x)| <= 1 on [-1, 1]:
  max |P'(x)| <= n^2 on [-1, 1] (factor n^2 instead of n).

VERIFICATION
============

  1. Bernstein on unit circle: max|P'| <= n max|P|.
  2. Equality at P(z) = z^n.
  3. Markov on [-1, 1]: max|P'| <= n^2 max|P|.
"""
from __future__ import annotations
import math
import numpy as np


def max_on_circle(P, n_pts=512):
    z = np.exp(2j * np.pi * np.arange(n_pts) / n_pts)
    return float(np.max(np.abs(np.polyval(P, z))))


def main():
    print("=" * 80)
    print("PAPER 39 — Bernstein and Markov inequalities")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Bernstein on unit circle: max|P'| <= n max|P|")
    label_pp = "max|Pp|"
    print(f"  {'P':>30} {'max|P|':>10} {label_pp:>12} {'ratio':>10} {'n':>4}")
    test_polys = [
        ("z^4", np.array([1.0, 0, 0, 0, 0])),  # extremal
        ("z^4 + z^2", np.array([1.0, 0, 1, 0, 0])),
        ("z^5 - 1", np.array([1.0, 0, 0, 0, 0, -1])),
        ("random d=8", rng.standard_normal(9)),
    ]
    for name, P in test_polys:
        Pp = np.polyder(P)
        n = len(P) - 1
        mP = max_on_circle(P)
        mPp = max_on_circle(Pp)
        ratio = mPp / mP if mP > 0 else 0
        ok = "OK" if ratio <= n + 1e-6 else "VIOLATES"
        print(f"  {name:>30} {mP:>10.4f} {mPp:>12.4f} {ratio:>10.4f} {n:>4} {ok}")

    print("\n[2] Markov on [-1, 1]: max|P'| <= n^2 max|P|")
    xs = np.linspace(-1, 1, 500)
    label_pp = "max|Pp|"
    print(f"  {'P':>30} {'max|P|':>10} {label_pp:>12} {'ratio':>10} {'n^2':>6}")
    for name, P in [
        ("Cheb T_3 = 4x^3-3x", np.array([4, 0, -3, 0])),
        ("Cheb T_4 = 8x^4-8x^2+1", np.array([8, 0, -8, 0, 1])),
        ("z^3 - z", np.array([1, 0, -1, 0.0])),
    ]:
        Pp = np.polyder(P)
        n = len(P) - 1
        mP = max(abs(np.polyval(P, x)) for x in xs)
        mPp = max(abs(np.polyval(Pp, x)) for x in xs)
        ratio = mPp / mP if mP > 0 else 0
        ok = "OK" if ratio <= n**2 + 1e-6 else "VIOLATES"
        print(f"  {name:>30} {mP:>10.4f} {mPp:>12.4f} {ratio:>10.4f} {n**2:>6} {ok}")


if __name__ == "__main__":
    main()
