"""
PAPER: 048 (canonical: 38pandrosion_knots.pdf)
TITLE: Pandrosion Quotients of Alexander Polynomials
STATUS: framework
DEPENDS: 047

THEORY
======

For a knot K, the Alexander polynomial Delta_K(t) in Z[t] is a knot
invariant. The Pandrosion Q-chain at t = 1 produces:
  Q_n(1, 1) = Delta_K^(n)(1) / n!.
These are related to Vassiliev (finite-type) invariants.

  Q_0 = Delta_K(1) = ±1
  Q_1 = Delta_K'(1)
  Q_2 = Delta_K''(1)/2  =  v_2 (second Vassiliev invariant)

VERIFICATION
============

  1. Compute Q_n for known knot Alexander polynomials.
  2. Compare with Vassiliev invariants.
"""
from __future__ import annotations
import numpy as np


def alexander_3_1():
    """Trefoil 3_1: Delta = t^2 - t + 1."""
    return np.array([1, -1, 1])


def alexander_4_1():
    """Figure-8 4_1: Delta = -t^2 + 3t - 1, normalized: t^2 - 3t + 1."""
    return np.array([1, -3, 1])


def alexander_5_1():
    """5_1: Delta = t^4 - t^3 + t^2 - t + 1."""
    return np.array([1, -1, 1, -1, 1])


def Q_chain_at_1(P, n_max=4):
    """Q_n(1, 1) = P^(n)(1) / n!."""
    import math
    Q = []
    Pk = P.copy()
    for n in range(n_max):
        Q.append(np.polyval(Pk, 1.0) / math.factorial(n))
        Pk = np.polyder(Pk)
    return Q


def main():
    print("=" * 80)
    print("PAPER 48 — Vassiliev knot invariants via Pandrosion Q-chain")
    print("=" * 80)

    knots = [
        ("3_1 (trefoil)", alexander_3_1()),
        ("4_1 (figure-8)", alexander_4_1()),
        ("5_1", alexander_5_1()),
    ]

    print("\n[1] Pandrosion Q_n(1, 1) = Delta^(n)(1) / n!")
    print(f"  {'knot':>20} {'Q_0':>6} {'Q_1':>6} {'Q_2':>6} {'Q_3':>6} {'Q_4':>6}")
    for name, P in knots:
        Q = Q_chain_at_1(P, 5)
        Q_str = " ".join(f"{int(round(q.real)):>6}" for q in Q[:5])
        print(f"  {name:>20} {Q_str}")

    print("\n[2] Vassiliev v_2 = Q_2 (second invariant)")
    for name, P in knots:
        Q = Q_chain_at_1(P, 3)
        v2 = Q[2]
        print(f"  {name}: v_2 = Q_2 = {v2.real:.4f}")


if __name__ == "__main__":
    main()
