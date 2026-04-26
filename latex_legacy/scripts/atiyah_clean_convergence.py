"""Clean convergence test using log|D| from slogdet(M), avoiding eigenvalue
clipping artifacts.

The identity 2 log|D| = log det G_norm + L_n is exact, so we compute
   -log det G_norm = L_n - 2 log|D|
where log|D| is computed via slogdet(M) - sum log ||p_i||_SU(2). slogdet on M
is numerically MUCH more stable than eigvalsh on G_norm because the singular
values of M are spread O(n) decades vs O(n^2) for G_norm.
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


def neg_log_det_clean(points):
    """Compute -log det G_norm = L_n - 2 log|D| via slogdet(M)."""
    n = len(points)
    M = atiyah_matrix(points)
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    log_norms_2 = np.log(np.maximum(
        np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
    sg, logabs = np.linalg.slogdet(M)
    if not np.isfinite(logabs):
        return None
    log_D = logabs - 0.5 * log_norms_2.sum()
    L_n = float(np.sum(np.log(binom)))
    return L_n - 2.0 * log_D, log_D, L_n


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 110, flush=True)
    print("CLEAN convergence test using log|D| via slogdet(M), no eigenvalue clipping",
          flush=True)
    print("=" * 110, flush=True)
    print(f"\n{'n':>5} {'L_n':>11} {'random':>10} {'antipod':>10} "
          f"{'eq-eq01':>10} {'2cluster':>10} {'WORST -log det':>15} "
          f"{'beta_n':>9} {'log|D|>=0?':>11}", flush=True)
    print("-" * 110, flush=True)
    for n in [10, 30, 50, 100, 150, 200, 300, 400]:
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        all_log_Ds = []  # to verify Atiyah-Sutcliffe
        # random
        rng = np.random.default_rng(20260530 + n)
        worst_random = 0.0
        nc = 5 if n < 200 else 3
        for _ in range(nc):
            pts = random_S2(rng, n)
            res = neg_log_det_clean(pts)
            if res is not None:
                neg_ld, log_D, _ = res
                all_log_Ds.append(log_D)
                if neg_ld > worst_random: worst_random = neg_ld
        # antipodal
        rng = np.random.default_rng(20260531 + n)
        worst_anti = 0.0
        for _ in range(nc):
            pts = antipodal_split_S2(rng, n)
            res = neg_log_det_clean(pts)
            if res is not None:
                neg_ld, log_D, _ = res
                all_log_Ds.append(log_D)
                if neg_ld > worst_anti: worst_anti = neg_ld
        # equispaced
        pts = equispaced_circle_S2(n, 0.01)
        res = neg_log_det_clean(pts)
        if res is not None:
            neg_ld_eq, log_D_eq, _ = res
            all_log_Ds.append(log_D_eq)
        else:
            neg_ld_eq = 0.0
        # 2-cluster
        rng = np.random.default_rng(20260532 + n)
        worst_2c = 0.0
        for sep in [1.0, 2.0]:
            for _ in range(2):
                pts = two_cluster_S2(rng, n, sep)
                res = neg_log_det_clean(pts)
                if res is not None:
                    neg_ld, log_D, _ = res
                    all_log_Ds.append(log_D)
                    if neg_ld > worst_2c: worst_2c = neg_ld
        worst = max(worst_random, worst_anti, neg_ld_eq, worst_2c)
        ratio = worst / L_n
        min_log_D = min(all_log_Ds) if all_log_Ds else float('nan')
        AS_holds = "YES" if min_log_D >= 0 else "NO"
        print(f"{n:>5} {L_n:>11.2f} {worst_random:>10.2f} "
              f"{worst_anti:>10.2f} {neg_ld_eq:>10.2f} {worst_2c:>10.2f} "
              f"{worst:>15.2f} {ratio:>9.4f} {AS_holds:>11}",
              flush=True)
        print(f"      min log|D| over all tested cfgs: {min_log_D:.4f}",
              flush=True)


if __name__ == "__main__":
    main()
