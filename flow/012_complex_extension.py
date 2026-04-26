"""
PAPER: 012 (canonical: 02pandrosion_complex.pdf)
TITLE: Complex Extension of the Pandrosion Iterator
STATUS: proved (universal multi-start convergence theorem)
DEPENDS: 000, 011

THEORY
======

The univariate Pandrosion iterator h_{P, z_0}(z) = z_0 - P(z_0)/Q(z_0, z)
extends naturally to C. Its fixed points are the roots of P (since
P(zeta) = 0 implies h(zeta) = z_0).

MULTI-START SCHEME (paper 0 Part II):
seed_{a, eps}(z) = a + eps * (z - a)
places any z in C inside the quadratic-contraction disk of anchor a, with
eps small. Universal convergence in N_eff ~ log_2 log_2(1/eps) iterations.

VERIFICATION
============

  1. Pandrosion iterator on z^2 - 1 converges to ±1.
  2. Multi-start coverage: random starts hit all roots.
  3. Quadratic convergence rate near a root.
"""
from __future__ import annotations
import math
import numpy as np


def pandrosion_iter_complex(P, a, z, max_iter=100, tol=1e-13):
    orbit = [z]
    for _ in range(max_iter):
        Pz0 = np.polyval(P, a)
        if abs(z - a) < 1e-15:
            Q = np.polyval(np.polyder(P), a)
        else:
            Q = (np.polyval(P, z) - Pz0) / (z - a)
        if abs(Q) < 1e-15: break
        z_new = a - Pz0 / Q
        orbit.append(z_new)
        if abs(z_new - z) < tol: break
        z = z_new
    return orbit


def main():
    print("=" * 80)
    print("PAPER 12 — Complex Pandrosion iterator")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] z^2 - 1: convergence to ±1 from random complex starts")
    P = np.array([1, 0, -1.0])
    targets = [1.0, -1.0]
    hits = {0: 0, 1: 0}
    for _ in range(50):
        a = complex(rng.uniform(-2, 2), rng.uniform(-2, 2))
        z = complex(rng.uniform(-2, 2), rng.uniform(-2, 2))
        orbit = pandrosion_iter_complex(P, a, z)
        z_final = orbit[-1]
        d_to_t = [abs(z_final - t) for t in targets]
        if min(d_to_t) < 1e-3:
            hits[d_to_t.index(min(d_to_t))] += 1
    print(f"  hits per root: {hits}")

    print("\n[2] Quadratic convergence rate near a root")
    P = np.array([1, 0, 0, -1.0])  # z^3 - 1
    root = 1.0
    z0 = root + 0.3
    a = root  # exact-anchor case: quadratic
    orbit = pandrosion_iter_complex(P, a, z0, max_iter=10)
    errs = [abs(z - root) for z in orbit]
    print(f"  z^3 - 1, anchor at root, z_0 = 1.3:")
    print(f"  errors: {[f'{e:.2e}' for e in errs[:6]]}")
    # quadratic ratios e_{n+1} / e_n^2
    quads = [errs[i+1]/errs[i]**2 if errs[i] > 1e-12 else 0 for i in range(len(errs)-1)]
    print(f"  quadratic ratios: {[f'{r:.3f}' for r in quads[:4]]}")

    print("\n[3] Multi-start on z^4 - 1 (4 roots ±1, ±i)")
    P = np.array([1, 0, 0, 0, -1.0])
    targets = [1, -1, 1j, -1j]
    hits = [0, 0, 0, 0]
    for _ in range(80):
        a = complex(rng.uniform(-1.5, 1.5), rng.uniform(-1.5, 1.5))
        z = complex(rng.uniform(-1.5, 1.5), rng.uniform(-1.5, 1.5))
        orbit = pandrosion_iter_complex(P, a, z)
        zf = orbit[-1]
        ds = [abs(zf - t) for t in targets]
        if min(ds) < 1e-3:
            hits[ds.index(min(ds))] += 1
    print(f"  hits per root: {hits}, total = {sum(hits)}/80")


if __name__ == "__main__":
    main()
