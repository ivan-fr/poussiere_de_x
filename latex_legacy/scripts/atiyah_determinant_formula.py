"""Test: does det G_Atiyah have a closed form in terms of pairwise distances?

Inspired by paper 77 (Pandrosion-Hadamard): for the Gram matrix of Pandrosion
quotients of P on |z|=1, det G = |Delta(P)|^2 = prod |alpha_i - alpha_j|^2.

For Atiyah-Berry, the polynomials p_i have DIFFERENT root sets (a_{ij} =
stereographic projection of (x_j - x_i) from x_i). So Paper 77 doesn't directly
apply.

BUT: maybe there's an analogous formula in terms of |x_i - x_j| (chord distances)?

Test:
  log det G_norm vs. various candidate forms:
  (a) sum_{i<j} log |x_i - x_j|^p for various p
  (b) sum_{i<j} log (1 - <x_i, x_j>)^p
  (c) sum_{i<j} log sin^p(angle(x_i, x_j) / 2)

Look for clean linear fit.
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


def gram_norm(points):
    n = len(points)
    M = atiyah_matrix(points)
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    W = np.diag(1.0 / binom)
    G = M @ W @ M.conj().T
    d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
    return G / np.outer(d, d)


def log_chord_sums(pts):
    """Returns various symmetric statistics of pairwise chord distances."""
    n = len(pts)
    sums = {
        'log_chord': 0.0,
        'log_chord_sq': 0.0,
        'log_1_minus_dot': 0.0,  # log(1 - <x_i, x_j>) = log(|x_i - x_j|^2 / 2)
        'sum_log_pairs': 0.0,
        'count': 0,
    }
    for i in range(n):
        for j in range(i+1, n):
            d = float(np.linalg.norm(pts[i] - pts[j]))
            dot = float(np.dot(pts[i], pts[j]))
            sums['log_chord'] += math.log(max(d, 1e-300))
            sums['log_chord_sq'] += math.log(max(d, 1e-300)) * 2
            sums['log_1_minus_dot'] += math.log(max(1.0 - dot, 1e-300))
            sums['count'] += 1
    return sums


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("Hunt: does log|D|^2 admit a closed-form decomposition in terms of pairwise distances?",
          flush=True)
    print("=" * 100, flush=True)

    for n in [4, 5, 6, 8, 10]:
        rng = np.random.default_rng(20260508 + n)
        binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        Y_logD = []
        X_chord = []  # sum log |x_i - x_j|^2
        X_dot = []  # sum log(1 - <x_i, x_j>)
        for _ in range(80):
            pts = random_S2(rng, n)
            M = atiyah_matrix(pts)
            log_norms_2 = np.log(np.maximum(
                np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
            sign, logabs = np.linalg.slogdet(M)
            if not np.isfinite(logabs): continue
            log_D = logabs - 0.5 * log_norms_2.sum()
            sums = log_chord_sums(pts)
            Y_logD.append(log_D)
            X_chord.append(sums['log_chord_sq'])
            X_dot.append(sums['log_1_minus_dot'])
        Y = np.array(Y_logD)
        Xc = np.array(X_chord)
        Xd = np.array(X_dot)
        if len(Y) >= 3:
            # Linear fit: log|D| = a + b * X_chord
            slope_c, b_c = np.polyfit(Xc, Y, 1)
            slope_d, b_d = np.polyfit(Xd, Y, 1)
            corr_c = float(np.corrcoef(Xc, Y)[0, 1])
            corr_d = float(np.corrcoef(Xd, Y)[0, 1])
            print(f"\nn = {n}:")
            print(f"  log|D| vs sum log|x_i-x_j|^2: slope = {slope_c:+.4f},  intercept = {b_c:+.3f},  corr = {corr_c:.4f}")
            print(f"  log|D| vs sum log(1-<x_i,x_j>): slope = {slope_d:+.4f},  intercept = {b_d:+.3f},  corr = {corr_d:.4f}")

    # If correlation is ~1, the formula is essentially exact.
    # Try also: log|D| = c_n + sum_{i<j} f(|x_i - x_j|)
    # where f might be (n-2) * log(...) or similar.

    print("\n" + "=" * 100, flush=True)
    print("Is there an EXACT formula of the form log|D| = a_n + b_n * sum_{i<j} log|x_i-x_j|^2?",
          flush=True)
    print("=" * 100, flush=True)

    print(f"\n{'n':>4} {'#cfg':>5} {'mean log|D|':>14} {'mean sum log|d_ij|^2':>22} "
          f"{'fit slope':>11} {'fit intercept':>14} {'corr':>8} {'residual std':>14}",
          flush=True)
    for n in [4, 5, 6, 8, 10, 15]:
        rng = np.random.default_rng(20260509 + n)
        Y = []; X = []
        for _ in range(150):
            pts = random_S2(rng, n)
            M = atiyah_matrix(pts)
            binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
            log_norms_2 = np.log(np.maximum(
                np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
            sign, logabs = np.linalg.slogdet(M)
            if not np.isfinite(logabs): continue
            log_D = logabs - 0.5 * log_norms_2.sum()
            sumlog = sum(math.log(max(np.linalg.norm(pts[i]-pts[j]), 1e-300))
                        for i in range(n) for j in range(i+1, n))
            Y.append(log_D); X.append(2.0 * sumlog)
        Y = np.array(Y); X = np.array(X)
        if len(Y) >= 3:
            slope, b = np.polyfit(X, Y, 1)
            corr = float(np.corrcoef(X, Y)[0, 1])
            resid = Y - (slope * X + b)
            print(f"{n:>4} {len(Y):>5} {float(Y.mean()):>14.4f} "
                  f"{float(X.mean()):>22.4f} {slope:>11.4f} {b:>14.4f} "
                  f"{corr:>8.5f} {float(resid.std()):>14.4f}", flush=True)


if __name__ == "__main__":
    main()
