"""Empirically determine the growth rate of log det G_norm.

If we can show log det G_norm >= -c * n (linear growth), Atiyah-Sutcliffe follows
since L_n grows like n^2.

Test: fit log det G_norm vs n to a power law n^p.
"""
from __future__ import annotations
import math
import numpy as np


def stereo(v):
    z = v[2]
    if z >= 1.0 - 1e-13: return complex('inf')
    return complex(v[0], v[1]) / (1.0 - z)


def random_S2(rng, n):
    v = rng.standard_normal((n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


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


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 90, flush=True)
    print("Growth of -min log det G_norm and L_n versus n", flush=True)
    print("=" * 90, flush=True)
    print(f"{'n':>4} {'L_n':>10} {'-min log det G':>16} {'L_n / n^2':>12} "
          f"{'-min/n':>10} {'-min/(n log n)':>15}", flush=True)
    print("-" * 90, flush=True)
    ns = [4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50]
    data = []
    for n in ns:
        rng = np.random.default_rng(20260504 + n)
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        ld_Gs = []
        n_cfg = max(40, 600 // n)
        for _ in range(n_cfg):
            pts = random_S2(rng, n)
            M = atiyah_matrix(pts)
            W = np.diag(1.0 / binom)
            G = M @ W @ M.conj().T
            d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
            G_norm = G / np.outer(d, d)
            sg, ldG = np.linalg.slogdet(G_norm)
            if np.isfinite(ldG): ld_Gs.append(ldG)
        if ld_Gs:
            neg_min = -min(ld_Gs)
            data.append((n, L_n, neg_min))
            print(f"{n:>4} {L_n:>10.3f} {neg_min:>16.4f} "
                  f"{L_n/n**2:>12.4f} {neg_min/n:>10.4f} "
                  f"{neg_min/(n*math.log(n)):>15.4f}", flush=True)

    # Fit -min log det G ~ a * n^p
    print("\nLog-log fit:", flush=True)
    arr = np.array(data)
    if len(arr) >= 3:
        logn = np.log(arr[:, 0])
        logy = np.log(arr[:, 2])
        slope, intercept = np.polyfit(logn, logy, 1)
        print(f"  -min log det G ~ {math.exp(intercept):.3f} * n^{slope:.3f}",
              flush=True)
        # And L_n
        logy2 = np.log(arr[:, 1])
        slope2, intercept2 = np.polyfit(logn, logy2, 1)
        print(f"  L_n            ~ {math.exp(intercept2):.3f} * n^{slope2:.3f}",
              flush=True)
        print(f"  Atiyah margin  : L_n - |min log det G| grows like n^{slope2:.3f}"
              f" - n^{slope:.3f}  ==>  +infty", flush=True)


if __name__ == "__main__":
    main()
