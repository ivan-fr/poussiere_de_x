"""Greedy Schur-complement: pick at each step the index minimizing r_i^2.

If under greedy ordering, r_i^2 is bounded by some rho_n -> rho < 1, then
   sum -log(1-r_i^2) <= n log(1/(1-rho))
which gives O(n) deficit bound.

Note: det G is INVARIANT under permutation, so sum -log(1-r_i^2) is the same.
But the MAX r_i^2 depends on order — greedy may avoid the bad pivots.

Greedy: at step i, choose argmin_{j not yet chosen} (1 - g_j^* G_<^{-1} g_j)?
Actually, for det = prod (1-r^2), we want to MAXIMIZE 1-r^2 at each step,
i.e., MINIMIZE r^2 at each step.

But to BOUND the worst case, we'd want max_i r_i^2 small. Greedy minimizes
max_i r_i^2 by choosing the most-orthogonal vector first.

Wait actually the SUM is fixed (=-log det) by Schur. So greedy doesn't change
the sum, only the distribution. What changes? The DISTRIBUTION of r_i^2.
If greedy makes ALL r_i^2 below some threshold, we get a NON-TRIVIAL bound.
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


def greedy_schur(G):
    """Greedy Schur: at each step choose the column with smallest residual.

    Maintains:
      perm = list of chosen indices (in chosen order)
      remaining = set of unchosen indices
      r2_seq = the r_i^2 sequence under the greedy ordering
    """
    n = G.shape[0]
    perm = []
    remaining = set(range(n))
    r2_seq = []
    # Start: choose the index with maximum norm (all 1, so any)
    # but for stability let's pick index 0
    perm.append(0)
    remaining.remove(0)
    r2_seq.append(0.0)
    while remaining:
        best_j = None
        best_r2 = float('inf')
        # Submatrix at chosen indices
        idx = perm
        G_chosen = G[np.ix_(idx, idx)]
        try:
            G_chosen_inv = np.linalg.inv(G_chosen)
        except np.linalg.LinAlgError:
            G_chosen_inv = np.linalg.pinv(G_chosen)
        for j in remaining:
            g_j = G[idx, j]
            r2 = float(np.real(g_j.conj() @ G_chosen_inv @ g_j))
            if r2 < best_r2:
                best_r2 = r2
                best_j = j
        perm.append(best_j)
        remaining.remove(best_j)
        r2_seq.append(best_r2)
    return perm, r2_seq


def antigreedy_schur(G):
    """Anti-greedy: at each step choose the WORST (largest r2) — to find
    the worst-case max r_i^2 under any ordering."""
    n = G.shape[0]
    perm = [0]
    remaining = set(range(1, n))
    r2_seq = [0.0]
    while remaining:
        worst_j = None
        worst_r2 = -1.0
        idx = perm
        G_chosen = G[np.ix_(idx, idx)]
        try:
            G_chosen_inv = np.linalg.inv(G_chosen)
        except np.linalg.LinAlgError:
            G_chosen_inv = np.linalg.pinv(G_chosen)
        for j in remaining:
            g_j = G[idx, j]
            r2 = float(np.real(g_j.conj() @ G_chosen_inv @ g_j))
            if r2 > worst_r2:
                worst_r2 = r2
                worst_j = j
        perm.append(worst_j)
        remaining.remove(worst_j)
        r2_seq.append(worst_r2)
    return perm, r2_seq


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("Greedy Schur-complement: max r_i^2 under best ordering", flush=True)
    print("=" * 100, flush=True)
    print(f"\n{'family':>10} {'n':>4} {'#cfg':>5} {'naive max r^2':>15} "
          f"{'greedy max r^2':>16} {'1 - max':>13} {'-log(1-max)':>13}",
          flush=True)
    print("-" * 100, flush=True)
    fams = [
        ("random",   lambda rng, n: random_S2(rng, n)),
        ("antipod",  lambda rng, n: antipodal_split_S2(rng, n)),
        ("eq-eq",    lambda rng, n: equispaced_circle_S2(n, 0.01)),
    ]
    for fname, fgen in fams:
        for n in [6, 10, 15, 20, 30]:
            rng = np.random.default_rng(20260511 + n + hash(fname) % 1000)
            n_cfg = 10 if "eq-" not in fname else 1
            naive_max = 0.0
            greedy_max = 0.0
            for _ in range(n_cfg):
                pts = fgen(rng, n)
                Gn = gram_norm(pts)
                # Naive (lex) ordering
                _, r2_naive = greedy_schur_with_order(Gn, list(range(n)))
                # Greedy
                _, r2_greedy = greedy_schur(Gn)
                if max(r2_naive) > naive_max: naive_max = max(r2_naive)
                if max(r2_greedy) > greedy_max: greedy_max = max(r2_greedy)
            print(f"{fname:>10} {n:>4} {n_cfg:>5} {naive_max:>15.6f} "
                  f"{greedy_max:>16.6f} {1-greedy_max:>13.4e} "
                  f"{-math.log(max(1-greedy_max, 1e-300)):>13.4f}",
                  flush=True)


def greedy_schur_with_order(G, order):
    """Apply Schur in the given fixed order, return r_i^2 sequence."""
    n = len(order)
    r2_seq = [0.0]
    for i in range(1, n):
        idx = order[:i]
        j = order[i]
        G_chosen = G[np.ix_(idx, idx)]
        try:
            G_chosen_inv = np.linalg.inv(G_chosen)
        except np.linalg.LinAlgError:
            G_chosen_inv = np.linalg.pinv(G_chosen)
        g_j = G[idx, j]
        r2 = float(np.real(g_j.conj() @ G_chosen_inv @ g_j))
        r2_seq.append(r2)
    return order, r2_seq


if __name__ == "__main__":
    main()
