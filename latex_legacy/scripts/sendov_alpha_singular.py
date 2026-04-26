"""Sendov on the alpha-singular subset.
After paper 104, Sendov is verified empirically to d=1024 and proved analytically
on the alpha-regular subset (Proposition 4.1). The gap: alpha-singular polynomials.

This script:
  1. Identifies alpha-singular polynomials (those with small alpha(P, z) at critical points)
  2. Tests Sendov on this restricted subset
  3. Measures fraction alpha-singular vs alpha-regular as d grows
"""
from __future__ import annotations
import math, time
import numpy as np


def smale_alpha(coefs, z):
    """Compute alpha(P, z) = beta * gamma."""
    d = len(coefs) - 1
    if d < 1: return float('inf')
    poly = coefs.astype(complex)
    cs = []
    fact = 1.0
    for k in range(d + 1):
        try:
            cs.append(complex(np.polyval(poly[::-1], z)))
        except (OverflowError, ZeroDivisionError):
            cs.append(complex('inf'))
        if len(poly) <= 1: break
        idx = np.arange(1, len(poly))
        poly = poly[1:] * idx / (k + 1)
    if len(cs) < 2 or abs(cs[1]) < 1e-30 or not np.isfinite(cs[1]):
        return float('inf')
    beta = abs(cs[0]) / abs(cs[1])
    gamma = 0.0
    for k in range(2, len(cs)):
        if not np.isfinite(cs[k]): continue
        ratio = abs(cs[k]) / abs(cs[1])
        if ratio <= 0: continue
        try:
            val = ratio ** (1.0 / (k - 1))
        except (OverflowError, ValueError):
            val = float('inf')
        if val > gamma: gamma = val
    return beta * gamma


ALPHA_STAR = (13.0 - 3.0 * math.sqrt(17.0)) / 4.0


def sendov_check(roots):
    """Returns (V, max_alpha) where V is Sendov violation and max_alpha is the
    maximum alpha(P, zeta_j) over all roots (alpha-singularity measure)."""
    d = len(roots)
    P = np.poly(roots)  # descending
    P_asc = P[::-1]      # ascending
    deriv_asc = P_asc[1:] * np.arange(1, len(P_asc))
    crits = np.roots(deriv_asc[::-1])
    if len(crits) == 0: return float('inf'), float('inf')
    # Sendov violation
    D = np.abs(roots[:, None] - crits[None, :])
    min_dists = D.min(axis=1)
    V = float(min_dists.max()) - 1.0
    # Max alpha at critical points
    max_alpha = 0.0
    for xi in crits:
        a = smale_alpha(P_asc, xi)
        if np.isfinite(a) and a > max_alpha:
            max_alpha = a
    return V, max_alpha


def sample_disk_BW_like(d, rng):
    r = rng.beta(d, 1.0, size=d)
    theta = rng.uniform(0, 2*np.pi, size=d)
    return r * np.exp(1j * theta)


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 95, flush=True)
    print("Sendov scan with alpha-singularity tracking", flush=True)
    print(f"  alpha_star = {ALPHA_STAR:.4f}", flush=True)
    print("=" * 95, flush=True)
    print(f"{'d':>5} {'#polys':>7} {'#alpha-sing':>11} {'%sing':>6} "
          f"{'max V (sing)':>13} {'max V (reg)':>12}", flush=True)
    print("-" * 95, flush=True)
    for d in [9, 16, 32, 64, 128, 256]:
        rng = np.random.default_rng(20260427 + d)
        n_polys = max(50, 5000 // d)
        n_sing = 0
        max_V_sing = -float('inf')
        max_V_reg = -float('inf')
        for _ in range(n_polys):
            roots = sample_disk_BW_like(d, rng)
            try:
                V, ma = sendov_check(roots)
            except (np.linalg.LinAlgError, ValueError):
                continue
            if not math.isfinite(V): continue
            if ma >= ALPHA_STAR:  # alpha-singular
                n_sing += 1
                if V > max_V_sing: max_V_sing = V
            else:
                if V > max_V_reg: max_V_reg = V
        pct_sing = 100.0 * n_sing / n_polys
        print(f"{d:>5} {n_polys:>7} {n_sing:>11} {pct_sing:>5.1f}% "
              f"{max_V_sing:>13.4f} {max_V_reg:>12.4f}", flush=True)


if __name__ == "__main__":
    main()
