"""Attack the deficit bound: study off-diagonal structure of G_norm.

Goal: find a uniform bound |log det G_norm| <= C * n^p with p < 2.

Strategy: study off-diagonal magnitudes of G_norm, the operator/Frobenius norms,
and adversarial configurations (clustered, antipodal, equispaced, ...).

Bound to test: ||G_norm - I||_op < 1 implies det G_norm > 0; and Bhatia gives
  -log det G_norm <= ||G_norm - I||_F^2 / (1 - ||G_norm - I||_op)
when ||G_norm - I||_op < 1.

If we can show max off-diagonal |G_norm_ij| stays bounded by some rho_max < 1
uniformly, with rho_max independent of n, then ||G_norm-I||_op <= (n-1)*rho_max
which doesn't help. We need ||G_norm-I||_F^2 bound that grows like n^p.

By Hadamard-Fischer: log det G_norm = sum_i log(G_ii principal minor ratio).
For PSD matrices this is bounded.

This script does:
  1. Sample G_norm over many adversarial families and many n.
  2. Report max |off-diagonal|, ||G_norm-I||_F, ||G_norm-I||_op, log det.
  3. See if any family pushes ||G_norm-I||_op close to 1 (bad case).
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


def clustered_S2(rng, n, eps):
    """n points clustered in spherical cap of radius eps."""
    v = rng.standard_normal((n, 3))
    v[:, 2] = np.abs(v[:, 2]) + 1.0 / eps  # bias toward +z pole
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def antipodal_split_S2(rng, n):
    """Half points clustered near north pole, half near south."""
    v = rng.standard_normal((n, 3)) * 0.1
    half = n // 2
    v[:half, 2] += 1.0
    v[half:, 2] -= 1.0
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def equispaced_circle_S2(n, eps):
    """n points on a circle near the equator."""
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


def stats(G):
    n = G.shape[0]
    R = G - np.eye(n)
    off = np.abs(R)
    return dict(
        max_off=float(off.max()),
        fro=float(np.linalg.norm(R, 'fro')),
        op=float(np.linalg.norm(R, 2)),
    )


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("Off-diagonal structure of G_norm by configuration family", flush=True)
    print("=" * 100, flush=True)

    families = [
        ("random",   lambda rng, n: random_S2(rng, n)),
        ("eq-eps01", lambda rng, n: equispaced_circle_S2(n, 0.1)),
        ("eq-eps001",lambda rng, n: equispaced_circle_S2(n, 0.01)),
        ("antipod",  lambda rng, n: antipodal_split_S2(rng, n)),
        ("cluster1", lambda rng, n: clustered_S2(rng, n, 1.0)),
        ("cluster01",lambda rng, n: clustered_S2(rng, n, 0.1)),
    ]

    print(f"\n{'family':>10} {'n':>4} {'#cfg':>5} {'max|off|':>10} "
          f"{'||R||_F':>10} {'||R||_op':>10} {'-log det':>10} "
          f"{'L_n':>10} {'margin':>10}", flush=True)
    print("-" * 100, flush=True)
    for fname, fgen in families:
        for n in [6, 10, 15, 20, 30, 50]:
            rng = np.random.default_rng(20260505 + n + hash(fname) % 1000)
            n_cfg = 30 if "eq-" not in fname else 1  # equispaced is deterministic
            stats_list = []
            for _ in range(n_cfg):
                pts = fgen(rng, n) if "eq-" not in fname else fgen(None, n)
                if len(pts) != n: continue
                Gn = gram_norm(pts)
                s = stats(Gn)
                sg, ld = np.linalg.slogdet(Gn)
                if not np.isfinite(ld): continue
                s['nld'] = -ld
                stats_list.append(s)
            if not stats_list: continue
            binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
            L_n = float(np.sum(np.log(binom)))
            mx = max(s['max_off'] for s in stats_list)
            mxF = max(s['fro'] for s in stats_list)
            mxOp = max(s['op'] for s in stats_list)
            mxNld = max(s['nld'] for s in stats_list)
            margin = L_n - mxNld
            print(f"{fname:>10} {n:>4} {len(stats_list):>5} "
                  f"{mx:>10.4f} {mxF:>10.4f} {mxOp:>10.4f} "
                  f"{mxNld:>10.4f} {L_n:>10.4f} {margin:>10.4f}", flush=True)

    # Search for adversarial: vary cluster radius
    print("\n" + "=" * 100, flush=True)
    print("Adversarial sweep: cluster radius vs deficit (n=15)", flush=True)
    print("=" * 100, flush=True)
    n = 15
    binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
    L_n = float(np.sum(np.log(binom)))
    print(f"L_n at n=15 = {L_n:.3f}", flush=True)
    print(f"\n{'eps':>10} {'max|off|':>10} {'||R||_op':>10} {'-log det':>12}",
          flush=True)
    for eps in [1.0, 0.5, 0.3, 0.1, 0.05, 0.01, 0.001, 1e-4, 1e-5]:
        rng = np.random.default_rng(424242)
        worst_nld = 0.0
        worst_op = 0.0
        worst_off = 0.0
        for _ in range(20):
            pts = clustered_S2(rng, n, eps)
            Gn = gram_norm(pts)
            sg, ld = np.linalg.slogdet(Gn)
            if not np.isfinite(ld): continue
            s = stats(Gn)
            if -ld > worst_nld: worst_nld = -ld
            if s['op'] > worst_op: worst_op = s['op']
            if s['max_off'] > worst_off: worst_off = s['max_off']
        print(f"{eps:>10.0e} {worst_off:>10.4f} {worst_op:>10.4f} {worst_nld:>12.4f}",
              flush=True)


if __name__ == "__main__":
    main()
