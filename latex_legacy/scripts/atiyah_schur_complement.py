"""Schur-complement attack on the Atiyah Gram deficit.

For PSD G with G_ii = 1, by iterated Schur complement:
    det(G) = prod_i (1 - r_i^2)
where r_i^2 = g_i^* G_{<i}^{-1} g_i is the squared multiple correlation of the
i-th polynomial on the first i-1 polynomials.

If we can show r_i^2 <= rho < 1 uniformly (rho independent of n), then
    |log det G_norm| <= n * log(1 / (1 - rho)) = O(n)
which closes Atiyah-Sutcliffe since L_n ~ n^2 / 2.

Test:
  1. Compute r_i^2 for adversarial configurations (random, antipodal,
     equispaced) at various n.
  2. Report max r_i^2 and how it scales with n.
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


def schur_residuals(G):
    """Compute r_i^2 = 1 - (det leading minor i / det leading minor (i-1))
    using Schur complement: r_i^2 = g_i^* G_<i^{-1} g_i.
    """
    n = G.shape[0]
    r2 = np.zeros(n)
    for i in range(n):
        if i == 0:
            r2[i] = 0.0  # by convention
            continue
        G_lt = G[:i, :i]
        g_i = G[:i, i]
        try:
            sol = np.linalg.solve(G_lt, g_i)
            r2[i] = float(np.real(np.dot(g_i.conj(), sol)))
        except np.linalg.LinAlgError:
            r2[i] = 1.0  # degenerate
    return r2


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("Schur-complement Attack: max r_i^2 across configurations", flush=True)
    print("=" * 100, flush=True)
    print(f"\n{'family':>10} {'n':>4} {'#cfg':>5} {'max r_i^2':>11} "
          f"{'1 - max r^2':>13} {'-log(1-r^2)':>13} {'sum -log(1-r^2)':>17} "
          f"{'-log det':>10} {'L_n':>10}",
          flush=True)
    print("-" * 100, flush=True)
    fams = [
        ("random",   lambda rng, n: random_S2(rng, n)),
        ("antipod",  lambda rng, n: antipodal_split_S2(rng, n)),
        ("eq-eq",    lambda rng, n: equispaced_circle_S2(n, 0.01)),
        ("eq-eq-tiny", lambda rng, n: equispaced_circle_S2(n, 1e-5)),
    ]
    for fname, fgen in fams:
        for n in [6, 10, 15, 20, 30, 50]:
            rng = np.random.default_rng(20260510 + n + hash(fname) % 1000)
            n_cfg = 20 if "eq-" not in fname else 1
            max_r2_overall = 0.0
            sum_neg_log_overall = 0.0
            neg_log_det_overall = 0.0
            for _ in range(n_cfg):
                pts = fgen(rng, n)
                Gn = gram_norm(pts)
                r2 = schur_residuals(Gn)
                # Numerical stability: clip
                r2 = np.clip(r2, 0.0, 1.0 - 1e-15)
                neg_log_resid = -np.log(1 - r2)
                if r2.max() > max_r2_overall: max_r2_overall = r2.max()
                if neg_log_resid.sum() > sum_neg_log_overall:
                    sum_neg_log_overall = neg_log_resid.sum()
                sg, ld = np.linalg.slogdet(Gn)
                if np.isfinite(ld) and -ld > neg_log_det_overall:
                    neg_log_det_overall = -ld
            binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
            L_n = float(np.sum(np.log(binom)))
            print(f"{fname:>10} {n:>4} {n_cfg:>5} {float(max_r2_overall):>11.6f} "
                  f"{float(1 - max_r2_overall):>13.4e} "
                  f"{float(-math.log(max(1 - max_r2_overall, 1e-300))):>13.4f} "
                  f"{float(sum_neg_log_overall):>17.4f} "
                  f"{float(neg_log_det_overall):>10.4f} {L_n:>10.2f}",
                  flush=True)

    # Per-i analysis: what's the worst r_i^2 by index?
    print("\n" + "=" * 100, flush=True)
    print("Per-index r_i^2 spectrum (n=15, antipodal split, worst case)", flush=True)
    print("=" * 100, flush=True)
    n = 15
    rng = np.random.default_rng(424242)
    worst_r2 = None
    worst_max = 0
    for _ in range(50):
        pts = antipodal_split_S2(rng, n)
        Gn = gram_norm(pts)
        r2 = schur_residuals(Gn)
        if r2.max() > worst_max:
            worst_max = r2.max()
            worst_r2 = r2
    print(f"Worst r_i^2 sequence (max = {worst_max:.6f}):")
    for i, v in enumerate(worst_r2):
        print(f"  r_{i}^2 = {v:.6f},  1 - r_i^2 = {max(1-v, 1e-15):.4e}",
              flush=True)


if __name__ == "__main__":
    main()
