"""Push to n = 300+ to confirm beta_n convergence.

Use eigvalsh + sum log to avoid slogdet underflow.
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


def gram_norm(points):
    n = len(points)
    M = atiyah_matrix(points)
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    W = np.diag(1.0 / binom)
    G = M @ W @ M.conj().T
    d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
    return G / np.outer(d, d)


def neg_log_det_via_eigh(Gn):
    """Use eigvalsh; sum log(eig) avoiding zero by clipping at machine eps."""
    w = np.real(np.linalg.eigvalsh(Gn))
    # clip negative due to floating point
    w = np.maximum(w, 1e-300)
    return -float(np.sum(np.log(w)))


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 105, flush=True)
    print("PUSH n LARGE: confirm convergence beta_n -> ?", flush=True)
    print("=" * 105, flush=True)
    print(f"\n{'n':>5} {'L_n':>11} {'random':>10} {'antipod':>10} "
          f"{'eq-eq01':>10} {'2cluster':>10} {'WORST':>11} {'WORST/L_n':>12}",
          flush=True)
    print("-" * 105, flush=True)
    for n in [50, 100, 150, 200, 250, 300, 400]:
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        # random
        rng = np.random.default_rng(20260520 + n)
        worst_random = 0.0
        nc = 5 if n < 200 else 3
        for _ in range(nc):
            pts = random_S2(rng, n)
            Gn = gram_norm(pts)
            ld = neg_log_det_via_eigh(Gn)
            if ld > worst_random: worst_random = ld
        # antipodal
        rng = np.random.default_rng(20260521 + n)
        worst_anti = 0.0
        for _ in range(nc):
            pts = antipodal_split_S2(rng, n)
            Gn = gram_norm(pts)
            ld = neg_log_det_via_eigh(Gn)
            if ld > worst_anti: worst_anti = ld
        # equispaced
        pts = equispaced_circle_S2(n, 0.01)
        Gn = gram_norm(pts)
        ld_eq = neg_log_det_via_eigh(Gn)
        # 2-cluster
        rng = np.random.default_rng(20260522 + n)
        worst_2c = 0.0
        for sep in [1.0, 2.0]:
            for _ in range(2):
                pts = two_cluster_S2(rng, n, sep)
                Gn = gram_norm(pts)
                ld = neg_log_det_via_eigh(Gn)
                if ld > worst_2c: worst_2c = ld
        worst = max(worst_random, worst_anti, ld_eq, worst_2c)
        ratio = worst / L_n
        print(f"{n:>5} {L_n:>11.2f} {worst_random:>10.2f} "
              f"{worst_anti:>10.2f} {ld_eq:>10.2f} {worst_2c:>10.2f} "
              f"{worst:>11.2f} {ratio:>12.4f}",
              flush=True)


if __name__ == "__main__":
    main()
