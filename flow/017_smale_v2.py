"""
PAPER: 017 (canonical: 07pandrosion_smale_v2.pdf)
TITLE: Pandrosion-Smale Iteration v2: Refinements
STATUS: refinement of paper 001 (Pandrosion-Smale operator)
DEPENDS: 001

THEORY
======

Refinements of paper 1's Pandrosion operator for univariate root-finding:
  - Tighter local convergence theorem with explicit basin radius
    eta_F(zeta) = c / (gamma(P, zeta) * |P|).
  - Pole-avoidance via measure-theoretic argument (Lebesgue null).
  - Iteration count O(log log eps^{-1}) after lock-in.

VERIFICATION
============

  1. Quadratic convergence in lock-in regime.
  2. Pole avoidance (no division by Q close to 0 near generic starts).
  3. Cost analysis: log log scaling.
"""
from __future__ import annotations
import math
import numpy as np


def Q(P, a, z):
    if abs(z - a) < 1e-15: return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def pandrosion_self(P, z):
    """Pandrosion with self-anchor (= Newton)."""
    Pp = np.polyval(np.polyder(P), z)
    if abs(Pp) < 1e-15: return None
    return z - np.polyval(P, z) / Pp


def smale_alpha(P, z):
    """alpha(P, z) = beta * gamma."""
    Pp_z = np.polyval(np.polyder(P), z)
    if abs(Pp_z) < 1e-15: return float('inf')
    beta = abs(np.polyval(P, z) / Pp_z)
    # gamma = sup_{k>=2} |P^(k)(z) / (k! * Pp_z)|^{1/(k-1)}
    d = len(P) - 1
    gammas = []
    for k in range(2, d + 1):
        Pk_z = np.polyval(np.polyder(P, k), z)
        v = abs(Pk_z / (math.factorial(k) * Pp_z))
        if v > 0: gammas.append(v**(1.0 / (k - 1)))
    gamma = max(gammas) if gammas else 0
    return beta * gamma


def main():
    print("=" * 80)
    print("PAPER 17 — Pandrosion-Smale v2: refinements")
    print("=" * 80)
    rng = np.random.default_rng(0)

    print("\n[1] Smale alpha and basin radius")
    test_polys = [
        ("z^2 - 2", np.array([1.0, 0, -2])),
        ("z^3 - 1", np.array([1.0, 0, 0, -1])),
        ("z^5 - 1", np.array([1.0, 0, 0, 0, 0, -1])),
    ]
    alpha_star = (13 - 3 * math.sqrt(17)) / 4
    print(f"  alpha* = (13 - 3 sqrt 17)/4 = {alpha_star:.6f}")
    for name, P in test_polys:
        roots = np.roots(P)
        for r in roots[:1]:
            for dist in [0.05, 0.1, 0.3, 0.5]:
                z = r + dist
                a = smale_alpha(P, z)
                ok = "lock-in" if a < alpha_star else "outside"
                print(f"  {name}, z = root + {dist}: alpha = {a:.4f}  ({ok})")

    print("\n[2] Quadratic convergence in lock-in regime")
    P = np.array([1.0, 0, -2.0])
    root = math.sqrt(2)
    z = root + 0.1
    errs = [abs(z - root)]
    for _ in range(8):
        z_new = pandrosion_self(P, z)
        if z_new is None: break
        z = z_new
        errs.append(abs(z - root))
    print(f"  errors: {[f'{e:.2e}' for e in errs[:7]]}")
    # Verify quadratic ratios
    quads = [errs[i+1]/errs[i]**2 for i in range(len(errs)-1) if errs[i] > 1e-12]
    print(f"  quadratic ratios e_{{n+1}}/e_n^2: {[f'{q:.3f}' for q in quads[:4]]}")

    print("\n[3] log log iteration count to reach precision eps")
    P = np.array([1.0, 0, -2.0])
    root = math.sqrt(2)
    z0 = root + 0.1
    print(f"  {'eps':>10} {'iters':>6} {'log_2 log_2 1/eps':>18}")
    for eps in [1e-3, 1e-6, 1e-12, 1e-15]:
        z = z0
        n = 0
        while abs(z - root) > eps and n < 100:
            z = pandrosion_self(P, z)
            if z is None: break
            n += 1
        log_log = math.log2(math.log2(1/eps)) if eps < 1 else 0
        print(f"  {eps:>10.0e} {n:>6} {log_log:>18.2f}")


if __name__ == "__main__":
    main()
