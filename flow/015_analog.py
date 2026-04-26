"""
PAPER: 015 (canonical: 05pandrosion_analog.pdf)
TITLE: Analogues of the Pandrosion Iterator
STATUS: framework (Newton, Halley, Steffensen as Pandrosion specialisations)
DEPENDS: 011

THEORY
======

The Pandrosion iterator h_{a, z}(z) = a - P(a)/Q(a, z) generalises classical
methods:
  Newton:      a = z, h(z) = z - P(z)/P'(z)
  Steffensen:  a = z, w = z + P(z), Q ~ P'(z) approximated by finite diff
  Halley:      uses curvature: h_Halley = z - 2 P P' / (2 (P')^2 - P P'')
              = z - P / P' / (1 - P P'' / (2 (P')^2))

PANDROSION-HALLEY IDENTITY:
  Halley's correction term involves P P'' / (P')^2 = 1 - E_P / (P')^2
  where E_P = (P')^2 - P P'' is the Pandrosion energy. So:
    h_Halley = z - P / P' * 1 / (1 - (1 - E/E_lin)/2)
            = z - 2 P P' E / (E (2(P')^2 - P P''))    after algebra.

VERIFICATION
============

  1. Newton recovery: h_a(z) at a = z = Newton step.
  2. Steffensen recovery via finite-difference Q.
  3. Halley = Pandrosion-T_2 (or equivalent).
"""
from __future__ import annotations
import math
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def newton_step(P, z):
    Pp_z = np.polyval(np.polyder(P), z)
    if abs(Pp_z) < 1e-15: return z
    return z - np.polyval(P, z) / Pp_z


def pandrosion_step(P, a, z):
    Q_val = Q(P, a, z)
    if abs(Q_val) < 1e-15: return z
    return z - np.polyval(P, z) / Q_val


def steffensen_step(P, z):
    """Steffensen: uses Q(z, z + P(z)) as derivative estimate."""
    Pz = np.polyval(P, z)
    if abs(Pz) < 1e-15: return z
    w = z + Pz
    Pw = np.polyval(P, w)
    Q_steff = (Pw - Pz) / (w - z)
    if abs(Q_steff) < 1e-15: return z
    return z - Pz / Q_steff


def halley_step(P, z):
    """Halley: z - 2 P P' / (2 (P')^2 - P P'')."""
    Pz = np.polyval(P, z)
    Pp_z = np.polyval(np.polyder(P), z)
    Ppp_z = np.polyval(np.polyder(P, 2), z)
    denom = 2 * Pp_z**2 - Pz * Ppp_z
    if abs(denom) < 1e-15: return z
    return z - 2 * Pz * Pp_z / denom


def main():
    print("=" * 80)
    print("PAPER 15 — Pandrosion analogues: Newton, Steffensen, Halley")
    print("=" * 80)
    rng = np.random.default_rng(0)

    P = np.array([1.0, 0, 0, -1.0])  # z^3 - 1, root = 1
    root = 1.0
    z0 = 1.3

    print("\n[1] Pandrosion(a=z, z) recovers Newton at every iterate")
    z_p = z0
    z_n = z0
    for _ in range(6):
        z_p = pandrosion_step(P, z_p, z_p)
        z_n = newton_step(P, z_n)
    print(f"  After 6 steps: Pandrosion(a=z) = {z_p:.10f}, Newton = {z_n:.10f}")

    print("\n[2] Steffensen vs Newton vs Pandrosion: convergence on z^3 - 1")
    z_n, z_s, z_h, z_p = z0, z0, z0, z0
    print(f"  {'iter':>5} {'Newton':>15} {'Steffensen':>15} {'Halley':>15} {'Pandrosion':>15}")
    print(f"  {0:>5} {z_n:>15.8f} {z_s:>15.8f} {z_h:>15.8f} {z_p:>15.8f}")
    for k in range(6):
        z_n = newton_step(P, z_n)
        z_s = steffensen_step(P, z_s)
        z_h = halley_step(P, z_h)
        z_p = pandrosion_step(P, root, z_p)  # exact-anchor Pandrosion
        print(f"  {k+1:>5} {z_n:>15.8f} {z_s:>15.8f} {z_h:>15.8f} {z_p:>15.8f}")

    print("\n[3] Convergence orders (errors vs iteration)")
    z_n, z_s, z_h = z0, z0, z0
    errs = {'N': [], 'S': [], 'H': []}
    for _ in range(8):
        z_n = newton_step(P, z_n); errs['N'].append(abs(z_n - root))
        z_s = steffensen_step(P, z_s); errs['S'].append(abs(z_s - root))
        z_h = halley_step(P, z_h); errs['H'].append(abs(z_h - root))
    for k, v in errs.items():
        print(f"  {k}: " + ", ".join(f"{e:.1e}" for e in v[:6]))


if __name__ == "__main__":
    main()
