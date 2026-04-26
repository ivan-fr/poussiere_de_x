"""Spectral analysis of G_norm: can we prove -log det G_norm <= beta * L_n?

Key facts:
  - G_norm hermitian PSD, diag = 1, so trace = n.
  - lambda_i(G) in [0, max], sum = n.
  - log det G = sum log lambda_i.
  - L_n = sum_k log C(n-1,k) ~ n^2 / 2.

Strategy: bound the smallest eigenvalue lambda_min(G) from below and bound
sum_i log lambda_i in terms of lambda_min and trace.

We have: log det G >= n * log(lambda_min) (trivial, weak)
       : log det G >= sum log lambda_i >= log(prod lambda) (=log det G, circular)

Useful: for PSD with trace = n, log det <= 0 (AM-GM). For lower bound,
need spectral gap.

Concrete test:
  1. Plot eigenvalue distribution of G_norm for several configurations.
  2. Test: is min eigenvalue bounded below by exp(-c * n) or polynomial?
  3. If lambda_min(G) >= 1/poly(n), then -log det <= n log poly(n) = O(n log n).
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


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("Spectral analysis: smallest eigenvalue of G_norm vs n", flush=True)
    print("=" * 100, flush=True)
    print(f"{'family':>10} {'n':>4} {'min lambda_min':>15} {'max lambda_max':>16} "
          f"{'-log det':>10} {'L_n':>10} {'ratio':>8}", flush=True)
    print("-" * 100, flush=True)
    fams = [
        ("random",   lambda rng, n: random_S2(rng, n)),
        ("antipod",  lambda rng, n: antipodal_split_S2(rng, n)),
        ("eq-eq",    lambda rng, n: equispaced_circle_S2(n, 0.01)),
    ]
    for fname, fgen in fams:
        for n in [10, 20, 30, 50, 80]:
            rng = np.random.default_rng(20260506 + n + hash(fname) % 1000)
            n_cfg = 30 if "eq-" not in fname else 1
            min_lams = []
            max_lams = []
            log_dets = []
            for _ in range(n_cfg):
                if "eq-" in fname:
                    pts = fgen(rng, n)
                else:
                    pts = fgen(rng, n)
                Gn = gram_norm(pts)
                w = np.linalg.eigvalsh(Gn)
                w = np.real(w)
                min_lams.append(float(w.min()))
                max_lams.append(float(w.max()))
                sg, ld = np.linalg.slogdet(Gn)
                if np.isfinite(ld):
                    log_dets.append(-float(ld))
            if log_dets:
                ml = min(min_lams)
                Mx = max(max_lams)
                neg_ld = max(log_dets)
                binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
                L_n = float(np.sum(np.log(binom)))
                ratio = neg_ld / L_n
                print(f"{fname:>10} {n:>4} {ml:>15.4e} {Mx:>16.4f} "
                      f"{neg_ld:>10.4f} {L_n:>10.4f} {ratio:>8.4f}", flush=True)

    # Test: is min lambda_min >= c * (1/L_n) ?
    print("\n" + "=" * 100, flush=True)
    print("Test: log min(lambda_min) vs log L_n  -- is lambda_min >= c/poly(n)?",
          flush=True)
    print("=" * 100, flush=True)
    data = []
    for n in [4, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70]:
        rng = np.random.default_rng(20260507 + n)
        worst_lm = 1.0
        for _ in range(40):
            pts = antipodal_split_S2(rng, n)
            Gn = gram_norm(pts)
            w = np.real(np.linalg.eigvalsh(Gn))
            if w.min() < worst_lm:
                worst_lm = float(w.min())
        binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        data.append((n, worst_lm, L_n))
        print(f"  n={n:>3}: min lambda_min = {worst_lm:.3e},  L_n = {L_n:.2f},  "
              f"-log lambda_min / L_n = {-math.log(max(worst_lm, 1e-300))/L_n:.4f}",
              flush=True)
    arr = np.array(data)
    if len(arr) >= 3:
        logn = np.log(arr[:, 0])
        log_lm = -np.log(np.maximum(arr[:, 1], 1e-300))
        slope, b = np.polyfit(logn, log_lm, 1)
        print(f"\n  -log lambda_min ~ {math.exp(b):.3f} * n^{slope:.3f}", flush=True)


if __name__ == "__main__":
    main()
