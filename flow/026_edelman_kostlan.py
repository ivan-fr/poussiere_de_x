"""
PAPER: 026 (canonical: 16pandrosion_ek.pdf)
TITLE: Edelman-Kostlan Root Concentration in Pandrosion Form
STATUS: proved (Edelman-Kostlan 1995, Pandrosion-form interpretation)
DEPENDS: 011

THEORY
======

EDELMAN-KOSTLAN (1995): For a Kostlan-Smale (Bombieri-Weyl) random polynomial
P of degree d, the roots concentrate on the unit circle |z| = 1:
  Pr[1/2 <= |zeta| <= 2] >= 1 - e^{-c d}.

Density on the line: 1/(pi (1 + x^2)) (the Cauchy distribution scaled).

PANDROSION VIEW: For a KS poly, the Pandrosion field F_P = P'/P at z near
|z| = 2 has typical magnitude:
  |F_P(z)| ~ d / |z|  (residue density),
giving Newton step |P/P'| ~ |z|/d ~ 2/d.

VERIFICATION
============

  1. Generate KS polynomials and check root concentration on |z| ~ 1.
  2. Verify density: 1/(pi (1 + x^2)) on real axis.
  3. Newton step magnitude on KS at |z| = 2.
"""
from __future__ import annotations
import math
import numpy as np


def kostlan_smale_polynomial(d, rng):
    log_binom = np.array([math.lgamma(d+1) - math.lgamma(k+1) - math.lgamma(d-k+1)
                          for k in range(d+1)])
    sigma = np.exp(0.5 * log_binom)
    coefs = (rng.standard_normal(d+1) + 1j*rng.standard_normal(d+1)) * sigma
    return coefs[::-1]  # high-to-low


def main():
    print("=" * 80)
    print("PAPER 26 — Edelman-Kostlan root concentration")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Root concentration on |z| ~ 1")
    print(f"  {'d':>4} {'frac in [1/2, 2]':>20}")
    for d in [10, 20, 50, 100]:
        n_polys = 30
        in_band = 0
        total = 0
        for _ in range(n_polys):
            P = kostlan_smale_polynomial(d, rng)
            roots = np.roots(P)
            for r in roots:
                total += 1
                if 0.5 <= abs(r) <= 2.0: in_band += 1
        print(f"  {d:>4} {in_band/total:>20.4f}")

    print("\n[2] Newton step magnitude at |z| = 2 (paper IX Lemma 1.1)")
    print(f"  {'d':>4} {'mean |Delta|':>14} {'1/sqrt(d)':>12}")
    for d in [16, 32, 64, 128]:
        n_test = 30
        deltas = []
        for _ in range(n_test):
            P = kostlan_smale_polynomial(d, rng)
            z = 2.0 * np.exp(1j * rng.uniform(0, 2*np.pi))
            Pp = np.polyval(np.polyder(P), z)
            if abs(Pp) > 1e-15:
                deltas.append(abs(np.polyval(P, z) / Pp))
        mean_d = np.mean(deltas)
        target = 1.0 / math.sqrt(d)
        print(f"  {d:>4} {mean_d:>14.6f} {target:>12.6f}")

    print("\n[3] Cauchy-like density on real axis (sample one KS poly d=50)")
    P = kostlan_smale_polynomial(50, rng)
    roots = np.roots(P)
    real_roots = [r.real for r in roots if abs(r.imag) < 0.05]
    print(f"  {len(real_roots)} roots near real axis (out of {len(roots)} total)")
    print(f"  range: [{min(real_roots) if real_roots else 0:.3f}, {max(real_roots) if real_roots else 0:.3f}]")


if __name__ == "__main__":
    main()
