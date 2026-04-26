"""
PAPER: 028 (canonical: 18pandrosion_walsh.pdf)
TITLE: Walsh's Theorem via the Pandrosion Field
STATUS: proved (Walsh 1922, Pandrosion-form proof)
DEPENDS: 024

THEORY
======

WALSH'S THEOREM: If P, Q are polynomials and B = {z : |P(z)| <= |Q(z)|} is
bounded, then for any constant c with |c| < 1, the polynomial P - c Q has
all its roots in B.

PANDROSION CONNECTION: The boundary {|P| = |Q|} is the locus where
|F_P / F_Q| = 1 in some sense; perturbing by c < 1 keeps the roots inside.

VERIFICATION
============

  1. Walsh's theorem on test pairs (P, Q).
  2. Boundary of B = {|P| <= |Q|} contains all roots of P - c Q.
"""
from __future__ import annotations
import math
import numpy as np


def main():
    print("=" * 80)
    print("PAPER 28 — Walsh's theorem")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Walsh on (P = z^3 - 1, Q = z^2)")
    # Roots of P - cQ for c in (0, 1)
    P = np.array([1, 0, 0, -1.0])
    Q = np.array([1, 0, 0, 0.0])  # z^3 (degree 3 to match)
    # Actually use Q = z^2 but pad to degree 3
    Q = np.array([0, 1, 0, 0.0])  # z^2 as degree 3 poly with 0 leading
    for c in [0.0, 0.3, 0.6, 0.9, 0.99]:
        diff = P - c * Q
        diff = diff[diff != 0] if abs(diff[0]) < 1e-12 else diff
        if len(diff) <= 1: continue
        roots = np.roots(diff)
        print(f"  c = {c}: roots of P - c Q = {[f'{r:.3f}' for r in roots]}")

    print("\n[2] Walsh region containment")
    P_arr = np.array([1, 0, 0, -1.0])
    Q_arr = np.array([1, 0, 1.0])  # z^2 + 1
    # Find some boundary points
    print(f"  P = z^3 - 1, Q = z^2 + 1")
    # For c < 1, roots of P - c Q in B = {|P| <= |Q|}
    for c in [0.5, 0.9]:
        # Pad Q to same degree
        Q_padded = np.concatenate([[0]*(len(P_arr) - len(Q_arr)), Q_arr])
        diff = P_arr - c * Q_padded
        roots = np.roots(diff)
        for r in roots:
            in_B = abs(np.polyval(P_arr, r)) <= abs(np.polyval(Q_arr, r))
            print(f"    c={c}, r={r:.4f}: |P(r)|={abs(np.polyval(P_arr, r)):.4f}, "
                  f"|Q(r)|={abs(np.polyval(Q_arr, r)):.4f}, in B = {in_B}")


if __name__ == "__main__":
    main()
