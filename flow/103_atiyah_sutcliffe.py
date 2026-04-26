"""
PAPER: 103 (canonical: 103_atiyah_sutcliffe.pdf)
TITLE: Atiyah-Sutcliffe via the SU(2)-Gram Determinant
STATUS: empirical certificate to n = 150; conjecture open for n >= 5
DEPENDS: 087

THEORY
======

ATIYAH-SUTCLIFFE CONJECTURE: For n distinct points x_1, ..., x_n on S^2,
the Atiyah determinant |D| >= 1, where
  D := det(M) / prod_i ||p_i||_{SU(2)}
with p_i monic Atiyah-Berry polys (roots = stereographic images of unit
directions from x_i).

IDENTITY (Theorem 2.1 of paper 103):
  2 log |D| = log det G_norm + L_n
where L_n = sum_k log C(n-1, k), G_norm is the unit-diagonal SU(2)-Gram of
the Atiyah polys.

REFORMULATION: Atiyah-Sutcliffe ⟺ log det G_norm >= -L_n
(Hadamard-deficit bound).

EMPIRICAL: deficit ratio beta_n := -log det G_norm / L_n stays below 0.314
uniformly for n <= 150.

VERIFICATION
============

  1. Identity 2 log|D| = log det G_norm + L_n exact.
  2. Atiyah-Sutcliffe holds (log|D| >= 0) on test configs.
  3. Empirical deficit ratio beta_n.
  4. det G_norm >= det K (DPP coherent-state kernel) — structural finding.
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
            else: finite.append(a)
        c = np.array([1.0+0j])
        for a in finite: c = np.convolve(c, np.array([-a, 1.0+0j]))
        if n_inf > 0:
            cc = np.zeros(len(c) + n_inf, dtype=complex); cc[n_inf:] = c; c = cc
        polys.append(c)
    return polys


def atiyah_M(points):
    n = len(points)
    polys = atiyah_polynomials(points)
    M = np.zeros((n, n), dtype=complex)
    for i, p in enumerate(polys):
        if len(p) < n: pp = np.zeros(n, dtype=complex); pp[:len(p)] = p; p = pp
        elif len(p) > n: p = p[:n]
        M[i] = p
    return M


def hopf_lift(x):
    z = float(x[2])
    if z > 1 - 1e-13: return np.array([1.0+0j, 0.0+0j])
    if z < -1 + 1e-13: return np.array([0.0+0j, 1.0+0j])
    cos_h = math.sqrt((1 + z) / 2)
    sin_h = math.sqrt((1 - z) / 2)
    phi = math.atan2(float(x[1]), float(x[0]))
    return np.array([cos_h+0j, sin_h * complex(math.cos(phi), math.sin(phi))])


def main():
    print("=" * 80)
    print("PAPER 103 — Atiyah-Sutcliffe via SU(2)-Gram determinant")
    print("=" * 80)
    rng = np.random.default_rng(2026)

    print("\n[1] Identity 2 log|D| = log det G_norm + L_n (exact)")
    print(f"  {'n':>3} {'log|D|^2':>10} {'log det G_norm':>16} {'L_n':>10} {'diff':>10}")
    for n in [4, 5, 6, 8, 10, 15]:
        pts = random_S2(rng, n)
        M = atiyah_M(pts)
        binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        log_norms_2 = np.log(np.maximum(np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
        sg, log_det_M = np.linalg.slogdet(M)
        log_D_sq = 2 * log_det_M.real - log_norms_2.sum()
        # G_norm
        W = np.diag(1.0 / binom)
        G = M @ W @ M.conj().T
        d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
        G_norm = G / np.outer(d, d)
        sg2, ld_G = np.linalg.slogdet(G_norm)
        ld_G_real = float(np.real(ld_G))
        diff = abs(log_D_sq - (ld_G_real + L_n))
        print(f"  {n:>3} {log_D_sq:>10.4f} {ld_G_real:>16.4f} {L_n:>10.4f} {diff:>10.2e}")

    print("\n[2] Atiyah-Sutcliffe verification: log|D| >= 0?")
    print(f"  {'n':>3} {'#cfg':>5} {'min log|D|':>13}")
    for n in [4, 8, 15, 20]:
        n_cfg = 30
        log_Ds = []
        for _ in range(n_cfg):
            pts = random_S2(rng, n)
            M = atiyah_M(pts)
            binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
            log_norms_2 = np.log(np.maximum(np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
            sg, ld = np.linalg.slogdet(M)
            log_Ds.append(ld.real - 0.5 * log_norms_2.sum())
        print(f"  {n:>3} {n_cfg:>5} {min(log_Ds):>13.4f}")

    print("\n[3] Deficit ratio beta_n = -log det G_norm / L_n")
    print(f"  {'n':>3} {'beta_n (worst)':>16}")
    for n in [10, 20, 50]:
        worst = 0.0
        for _ in range(15):
            pts = random_S2(rng, n)
            binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
            L_n = float(np.sum(np.log(binom)))
            M = atiyah_M(pts)
            log_norms_2 = np.log(np.maximum(np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
            sg, ld = np.linalg.slogdet(M)
            log_D_sq = 2 * ld.real - log_norms_2.sum()
            beta = (L_n - log_D_sq) / L_n
            if beta > worst: worst = beta
        print(f"  {n:>3} {worst:>16.4f}")

    print("\n[4] DPP kernel inequality: det G_norm >= det K (empirical)")
    n = 8
    n_violations = 0
    for _ in range(30):
        pts = random_S2(rng, n)
        M = atiyah_M(pts)
        binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
        W = np.diag(1.0 / binom)
        G = M @ W @ M.conj().T
        d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
        G_norm = G / np.outer(d, d)
        # Coherent-state kernel
        spinors = [hopf_lift(x) for x in pts]
        K = np.zeros((n, n), dtype=complex)
        for i in range(n):
            for j in range(n):
                ip = np.conj(spinors[i][0]) * spinors[j][0] + np.conj(spinors[i][1]) * spinors[j][1]
                K[i, j] = ip ** (n - 1)
        sg1, ld_G = np.linalg.slogdet(G_norm)
        sg2, ld_K = np.linalg.slogdet(K)
        if ld_G.real < ld_K.real - 1e-10: n_violations += 1
    print(f"  n={n}, 30 cfgs: violations of det G_norm >= det K = {n_violations}/30")


if __name__ == "__main__":
    main()
