"""Push n to large values and check if ratio -log det G / L_n converges.

If beta_n := sup -log det G_norm(x_1,...,x_n) / L_n converges to some
beta_infty < 1 as n -> infty, Atiyah-Sutcliffe holds asymptotically with
margin (1 - beta_infty) * L_n -> infty.

Empirically beta_n <= 0.30 for n <= 80. Push to n = 200.
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


def antipodal_split_S2(rng, n):
    v = rng.standard_normal((n, 3)) * 0.1
    half = n // 2
    v[:half, 2] += 1.0
    v[half:, 2] -= 1.0
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def equispaced_circle_S2(n, eps):
    pts = np.zeros((n, 3))
    for k in range(n):
        theta = 2 * np.pi * k / n
        pts[k, 0] = math.cos(theta) * math.sqrt(1 - eps**2)
        pts[k, 1] = math.sin(theta) * math.sqrt(1 - eps**2)
        pts[k, 2] = eps
    return pts


def two_cluster_S2(rng, n, sep):
    """Two clusters separated by spherical distance sep."""
    half = n // 2
    v = np.zeros((n, 3))
    for i in range(half):
        u = rng.standard_normal(3) * 0.05
        u[2] += math.cos(sep / 2)
        v[i] = u / np.linalg.norm(u)
    for i in range(half, n):
        u = rng.standard_normal(3) * 0.05
        u[2] -= math.cos(sep / 2)
        v[i] = u / np.linalg.norm(u)
    return v


def atiyah_polynomials(points):
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


def atiyah_matrix(points):
    n = len(points)
    polys = atiyah_polynomials(points)
    M = np.zeros((n, n), dtype=complex)
    for i, p in enumerate(polys):
        if len(p) < n:
            pp = np.zeros(n, dtype=complex)
            pp[:len(p)] = p
            p = pp
        elif len(p) > n:
            p = p[:n]
        M[i] = p
    return M


def gram_norm(points):
    n = len(points)
    M = atiyah_matrix(points)
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    W = np.diag(1.0 / binom)
    G = M @ W @ M.conj().T
    d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
    return G / np.outer(d, d)


def log_det_G_norm(points):
    Gn = gram_norm(points)
    sg, ld = np.linalg.slogdet(Gn)
    return float(ld) if np.isfinite(ld) else None


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 105, flush=True)
    print("CONVERGENCE TEST: does -log det G_norm / L_n converge as n -> infty?",
          flush=True)
    print("=" * 105, flush=True)
    print(f"\n{'n':>5} {'L_n':>10} {'random':>10} {'antipod':>10} "
          f"{'eq-eq01':>10} {'2cluster':>10} {'WORST':>10} {'WORST/L_n':>12}",
          flush=True)
    print("-" * 105, flush=True)
    for n in [10, 20, 30, 50, 70, 100, 130, 160, 200]:
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        results = {}
        # random
        rng = np.random.default_rng(20260512 + n)
        worst_random = 0.0
        for _ in range(max(5, 50 // (n // 10))):
            pts = random_S2(rng, n)
            ld = log_det_G_norm(pts)
            if ld is not None and -ld > worst_random: worst_random = -ld
        results['random'] = worst_random
        # antipodal
        rng = np.random.default_rng(20260513 + n)
        worst_anti = 0.0
        for _ in range(max(5, 50 // (n // 10))):
            pts = antipodal_split_S2(rng, n)
            ld = log_det_G_norm(pts)
            if ld is not None and -ld > worst_anti: worst_anti = -ld
        results['antipod'] = worst_anti
        # equispaced near equator
        pts = equispaced_circle_S2(n, 0.01)
        ld = log_det_G_norm(pts)
        results['eq-eq01'] = -ld if ld is not None else 0.0
        # 2-cluster with various separations
        rng = np.random.default_rng(20260514 + n)
        worst_2c = 0.0
        for sep in [0.5, 1.0, 1.5, 2.0]:
            for _ in range(3):
                pts = two_cluster_S2(rng, n, sep)
                ld = log_det_G_norm(pts)
                if ld is not None and -ld > worst_2c: worst_2c = -ld
        results['2cluster'] = worst_2c
        worst = max(results.values())
        ratio = worst / L_n
        t0 = time.time()
        print(f"{n:>5} {L_n:>10.2f} {results['random']:>10.2f} "
              f"{results['antipod']:>10.2f} {results['eq-eq01']:>10.2f} "
              f"{results['2cluster']:>10.2f} {worst:>10.2f} {ratio:>12.4f}",
              flush=True)


if __name__ == "__main__":
    main()
