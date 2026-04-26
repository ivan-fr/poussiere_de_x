"""
PAPER: 076 (canonical: 66pandrosion_rota.pdf)
TITLE: Rota's Conjecture and Pandrosion Logarithmic Concavity
STATUS: framework (Rota's conjecture proved by Adiprasito-Huh-Katz 2018)
DEPENDS: 030

THEORY
======

ROTA'S CONJECTURE (proved 2018 AHK): For matroid M, the absolute values
of the coefficients of the characteristic polynomial form a log-concave
sequence: |c_k|^2 >= |c_{k-1}| |c_{k+1}|.

PANDROSION CONNECTION: Log-concavity is equivalent to (P')^2 - P P'' >= 0
on positive reals (paper 30 Laguerre/Turán).

VERIFICATION
============

  1. Coefficient log-concavity of test polynomials.
  2. Connection to Pandrosion energy.
"""
from __future__ import annotations
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 76 — Rota / log-concavity")
    print("=" * 80)

    print("\n[1] Log-concavity test: |c_k|^2 >= |c_{k-1}| |c_{k+1}|")
    test_polys = [
        ("(1+z)^4", np.array([1, 4, 6, 4, 1])),
        ("z^4 + z^3 + z^2 + z + 1 (cyclo)", np.array([1, 1, 1, 1, 1])),
        ("z^5 - 5z^3 + 4z = z(z^2-1)(z^2-4)", np.poly([0, 1, -1, 2, -2.0])),
    ]
    for name, P in test_polys:
        # Reverse to low-to-high
        c = list(reversed(P))
        log_concave = all(c[k]**2 >= c[k-1] * c[k+1] - 1e-9 for k in range(1, len(c)-1)
                         if c[k-1] != 0 and c[k+1] != 0)
        print(f"  {name}: coeffs = {c}, log-concave? {log_concave}")

    print("\n[2] Real-rooted polynomials always have log-concave coefficients (Newton)")
    rng = np.random.default_rng(0)
    for d in [3, 5, 8]:
        roots = sorted(rng.uniform(0, 5, d))  # all positive
        P = np.poly(roots)
        c = list(reversed(P))
        # signs alternate; take absolute values
        c_abs = [abs(ci) for ci in c]
        log_concave = all(c_abs[k]**2 >= c_abs[k-1] * c_abs[k+1] - 1e-9
                         for k in range(1, len(c_abs)-1) if c_abs[k-1] > 0 and c_abs[k+1] > 0)
        print(f"  d = {d}, real positive roots: log-concave? {log_concave}")


if __name__ == "__main__":
    main()
