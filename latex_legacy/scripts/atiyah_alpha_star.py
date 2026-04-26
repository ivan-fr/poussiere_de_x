"""Atiyah-Sutcliffe: measuring alpha_star_AS for the Atiyah-Berry system.

Paper 103 proved |D| >= (alpha_star)^n on alpha-regular configurations.
Conjecture: alpha_star_AS = 1 (would give |D| >= 1 uniform).

This script estimates alpha_star_AS empirically by computing the Smale invariant
of the Atiyah-Berry polynomial system at random S^2 configurations.
If alpha < 1 always, the conjecture-bound holds.
"""
from __future__ import annotations
import math, time
import numpy as np


def stereo(v):
    z = v[2]
    if z >= 1.0 - 1e-13: return complex('inf')
    return complex(v[0], v[1]) / (1.0 - z)


def random_S2(rng, n):
    v = rng.standard_normal((n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def atiyah_polynomials(points):
    """Build the n Atiyah-Berry polynomials, return as list of coefs (ascending)."""
    n = len(points)
    polys = []
    for i in range(n):
        finite = []
        n_inf = 0
        for j in range(n):
            if i == j: continue
            v = points[j] - points[i]
            v = v / np.linalg.norm(v)
            a = stereo(v)
            if math.isinf(a.real) or math.isinf(a.imag):
                n_inf += 1
            else:
                finite.append(a)
        c = np.array([1.0 + 0j])
        for a in finite:
            c = np.convolve(c, np.array([-a, 1.0 + 0j]))
        if n_inf > 0:
            cc = np.zeros(len(c) + n_inf, dtype=complex)
            cc[n_inf:] = c
            c = cc
        polys.append(c)
    return polys


def smale_alpha_univariate(coefs, z):
    d = len(coefs) - 1
    poly = coefs.astype(complex)
    cs = []
    for k in range(d + 1):
        try:
            cs.append(complex(np.polyval(poly[::-1], z)))
        except (OverflowError, ZeroDivisionError):
            cs.append(complex('inf'))
        if len(poly) <= 1: break
        idx = np.arange(1, len(poly))
        poly = poly[1:] * idx / (k + 1)
    if len(cs) < 2 or abs(cs[1]) < 1e-30:
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


def estimate_alpha_AS(points):
    """Estimate alpha_AS = max_i alpha(p_i, z) over the Atiyah-Berry polynomials
    evaluated at z = 0 (the natural Pandrosion basepoint)."""
    polys = atiyah_polynomials(points)
    max_alpha = 0.0
    for p in polys:
        a = smale_alpha_univariate(p, 0.0+0j)
        if math.isfinite(a) and a > max_alpha:
            max_alpha = a
    return max_alpha


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 80, flush=True)
    print("Atiyah-Sutcliffe alpha_star_AS estimation", flush=True)
    print("=" * 80, flush=True)
    print(f"{'n':>4} {'#configs':>9} {'min alpha':>10} {'median':>10} "
          f"{'max alpha':>10} {'#alpha>1':>10}", flush=True)
    print("-" * 80, flush=True)
    for n in [5, 10, 15, 20, 30, 50]:
        rng = np.random.default_rng(20260428 + n)
        n_cfg = max(50, 1000 // (n // 2))
        alphas = []
        for _ in range(n_cfg):
            pts = random_S2(rng, n)
            a = estimate_alpha_AS(pts)
            if math.isfinite(a):
                alphas.append(a)
        if not alphas:
            print(f"{n:>4} no valid", flush=True)
            continue
        arr = np.array(alphas)
        n_above_1 = int((arr > 1.0).sum())
        print(f"{n:>4} {n_cfg:>9} {arr.min():>10.4f} {float(np.median(arr)):>10.4f} "
              f"{arr.max():>10.4f} {n_above_1:>10}", flush=True)


if __name__ == "__main__":
    main()
