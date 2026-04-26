"""
PAPER: 035 (canonical: 25pandrosion_arg.pdf)
TITLE: Argument Control via Pandrosion Energy Dominance
STATUS: framework

THEORY
======

ARGUMENT CONTROL: The Pandrosion step z -> a - P(a)/Q(a, z) has bounded
argument when the anchor field Q(a, z) "dominates" the tidal field
T(z) = E_P(z)/P(z)^2 = sum |Q(alpha_k, z)|^2 / |P(z)|^2.

ENERGY DOMINANCE: |Q(a, z)|^2 > E_P(z) ⟹ |arg(C)| < pi for the contraction
ratio C (paper 23). This is "the anchor's force is greater than the
combined tidal force from all roots".

VERIFICATION
============

  1. Energy dominance check on test polynomials.
  2. Argument bound verification.
"""
from __future__ import annotations
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def main():
    print("=" * 80)
    print("PAPER 35 — Argument control via energy dominance")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Energy dominance: |Q(a,z)|^2 > E_P(z)?")
    for d in [3, 5, 8]:
        roots_real = rng.uniform(-1, 1, d)
        P = np.poly(roots_real)
        z = 0.3 + 0.4j
        E_z = sum(abs(Q(P, ak, z))**2 for ak in roots_real)
        # Try anchor near a root
        a_good = roots_real[0] + 0.05  # near root: |Q(a, z)| big
        a_bad = 5.0  # far: |Q(a, z)| small (but in some sense)
        Q_good = abs(Q(P, a_good, z))**2
        Q_bad = abs(Q(P, a_bad, z))**2
        print(f"  d={d}: E_P(z) = {E_z:.4e}")
        print(f"    near-root anchor: |Q|^2 = {Q_good:.4e}, dominance = {Q_good > E_z}")
        print(f"    far anchor:       |Q|^2 = {Q_bad:.4e}, dominance = {Q_bad > E_z}")

    print("\n[2] Argument of contraction ratio C with energy dominance")
    P = np.array([1, 0, 0, -1.0])  # z^3 - 1
    roots = np.roots(P)
    z = 0.5 + 0.2j
    for a_test in [1.0, 0.7 + 0.2j, 5.0]:
        Q_a = Q(P, a_test, z)
        E_z = sum(abs(Q(P, r, z))**2 for r in roots)
        rs = [Q(P, r, z) / Q_a for r in roots]
        C = np.prod([1 - r for r in rs])
        dominance = abs(Q_a)**2 > E_z
        print(f"  anchor a={a_test:.3f}: dominance = {dominance}, "
              f"|arg C|/pi = {abs(np.angle(C))/np.pi:.4f}")


if __name__ == "__main__":
    main()
