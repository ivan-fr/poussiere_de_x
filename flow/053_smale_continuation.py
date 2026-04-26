"""
PAPER: 053 (canonical: 43pandrosion_smale.pdf)
TITLE: Smale-related Pandrosion identities (continuation)
STATUS: framework (continuation of papers 1, 17)
DEPENDS: 017

THEORY
======

Continuation of Pandrosion-Smale work: refinements of basin radius,
gamma-theory bounds, and convergence results.

Key formula (Smale 1986): If alpha(P, z) < alpha* = (13 - 3 sqrt 17)/4,
then Newton iteration starting at z converges quadratically to a root.
Pandrosion-T_n inherits this: in the lock-in regime, |error_{n+1}| <=
gamma * |error_n|^2.

VERIFICATION
============

  1. alpha* = (13 - 3 sqrt 17)/4 numerical value.
  2. Quadratic convergence in lock-in.
"""
from __future__ import annotations
import math
import numpy as np


def smale_alpha(P, z):
    Pp_z = np.polyval(np.polyder(P), z)
    if abs(Pp_z) < 1e-15: return float('inf')
    beta = abs(np.polyval(P, z) / Pp_z)
    d = len(P) - 1
    gammas = []
    for k in range(2, d + 1):
        Pk = np.polyval(np.polyder(P, k), z)
        v = abs(Pk / (math.factorial(k) * Pp_z))
        if v > 0: gammas.append(v**(1.0 / (k - 1)))
    return beta * (max(gammas) if gammas else 0)


def main():
    print("=" * 80)
    print("PAPER 53 — Smale-related Pandrosion identities (continuation)")
    print("=" * 80)

    alpha_star = (13 - 3 * math.sqrt(17)) / 4
    print(f"\n[1] Smale's alpha* = (13 - 3 sqrt 17)/4 = {alpha_star:.10f}")

    print("\n[2] Quadratic convergence in lock-in")
    P = np.array([1.0, 0, 0, 0, -1])  # z^4 - 1
    root = 1.0
    for delta in [0.05, 0.1, 0.2, 0.3]:
        z = root + delta
        a = smale_alpha(P, z)
        in_lockin = a < alpha_star
        # Newton step
        Pp = np.polyval(np.polyder(P), z)
        z_new = z - np.polyval(P, z) / Pp
        ratio = abs(z_new - root) / (abs(z - root)**2)
        print(f"  z = root + {delta}: alpha = {a:.4f} ({'lock-in' if in_lockin else 'OUT'}), "
              f"e_new/e^2 = {ratio:.4f}")

    print("\n[3] Smale's basin radius: 1/(4 gamma)")
    P = np.array([1.0, 0, 0, -1])  # z^3 - 1
    root = 1.0
    z = root + 0.0001
    Pp_z = np.polyval(np.polyder(P), z)
    gammas = []
    for k in range(2, 4):
        Pk = np.polyval(np.polyder(P, k), z)
        v = abs(Pk / (math.factorial(k) * Pp_z))
        if v > 0: gammas.append(v**(1.0 / (k - 1)))
    gamma = max(gammas) if gammas else 0
    print(f"  P = z^3 - 1, near root 1: gamma = {gamma:.4f}, basin = 1/(4 gamma) = {1/(4*gamma):.4f}")


if __name__ == "__main__":
    main()
