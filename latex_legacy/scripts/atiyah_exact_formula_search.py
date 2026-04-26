"""Search for an EXACT formula linking Atiyah-Berry |D|^2 to spinor invariants.

We have:
  |det V|^2 = prod_{i<j} |[u_i, u_j]|^2 = (1/4)^{n(n-1)/2} prod |x_i-x_j|^2  [EXACT]
  |D_Lag|^2 = |det M_L|^2 / prod ||P_i^Lag||^2  with M_L = V^{-T}
  |D_AB|^2 = |det M_AB|^2 / prod ||p_i^AB||^2

Test possible exact identities:
  Hypothesis 1: log|D_AB|^2 = log|D_Lag|^2 + a + b * sum log |x_i - x_j|^2
  Hypothesis 2: log|D_AB|^2 = -log|det V|^2 + c_n  (i.e., 1 / discriminant up to const)
  Hypothesis 3: log|D_AB|^2 = log(something) where something = symmetric pol in |[u_i, u_j]|^2

Specifically test: log|D_AB|^2 + log prod_{i<j} |x_i - x_j|^2 = constant of n only?
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


def atiyah_polys_AB(points):
    n = len(points)
    polys = []
    for i in range(n):
        finite = []; n_inf = 0
        for j in range(n):
            if i == j: continue
            v = points[j] - points[i]; v = v / np.linalg.norm(v)
            a = stereo(v)
            if math.isinf(a.real) or math.isinf(a.imag): n_inf += 1
            else: finite.append(a)
        c = np.array([1.0+0j])
        for a in finite: c = np.convolve(c, np.array([-a, 1.0+0j]))
        if n_inf > 0:
            cc = np.zeros(len(c)+n_inf, dtype=complex); cc[n_inf:] = c; c = cc
        polys.append(c)
    return polys


def atiyah_M_AB(points):
    n = len(points)
    polys = atiyah_polys_AB(points)
    M = np.zeros((n, n), dtype=complex)
    for i, p in enumerate(polys):
        if len(p) < n: pp = np.zeros(n, dtype=complex); pp[:len(p)] = p; p = pp
        elif len(p) > n: p = p[:n]
        M[i] = p
    return M


def log_D_AB_sq(points):
    n = len(points)
    M = atiyah_M_AB(points)
    binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
    log_norms_2 = np.log(np.maximum(np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
    sg, ld = np.linalg.slogdet(M)
    return 2*ld - log_norms_2.sum()


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("Hypothesis: log|D_AB|^2 + s*sum log|x_i-x_j|^2 = c_n (constant of n)", flush=True)
    print("Search for s that makes the LHS constant.", flush=True)
    print("=" * 100, flush=True)

    for n in [4, 5, 6, 8, 10]:
        rng = np.random.default_rng(20260602 + n)
        Y = []  # log|D_AB|^2
        X = []  # sum log |x_i - x_j|^2
        for _ in range(150):
            pts = random_S2(rng, n)
            ld_AB_2 = log_D_AB_sq(pts)
            sumlog_dsq = sum(2 * math.log(max(np.linalg.norm(pts[i] - pts[j]), 1e-300))
                            for i in range(n) for j in range(i+1, n))
            Y.append(ld_AB_2); X.append(sumlog_dsq)
        Y = np.array(Y); X = np.array(X)
        slope, intercept = np.polyfit(X, Y, 1)
        residuals = Y - (slope * X + intercept)
        std_resid = float(residuals.std())
        corr = float(np.corrcoef(X, Y)[0, 1])
        print(f"\nn = {n}:")
        print(f"  slope = {slope:+.6f},  intercept = {intercept:+.4f}")
        print(f"  correlation = {corr:.6f}")
        print(f"  residual std = {std_resid:.4f}")
        # Predicted slope value: -1 (i.e., log|D_AB|^2 + sum log|x_i-x_j|^2 = constant)?
        print(f"  Test slope = -1 (perfect cancellation): residual_std under fit slope = {std_resid:.4f}")
        # Try: log|D_AB|^2 = c_n - sum log|x_i - x_j|^2
        Y_alt = Y + X  # = log|D_AB|^2 + sum log|d|^2
        print(f"  log|D_AB|^2 + sum log|d_ij|^2:  mean = {Y_alt.mean():.4f}, std = {Y_alt.std():.4f}")
        # Try other slopes: 0.5, 1, 2
        for s in [-2, -1, -0.5, 0.5, 1, 2]:
            Y_test = Y - s * X
            print(f"    slope = {s:+5.2f}: log|D_AB|^2 - {s:+.2f} * sum log|d|^2: "
                  f"mean = {Y_test.mean():+.4f}, std = {Y_test.std():.4f}")

    # Try: maybe it's log(prod |[u_i, u_j]|^2) NOT log(prod |x_i - x_j|^2)?
    # These differ by constant 4^{-n(n-1)/2}, so same effective dependence.

    # Try: maybe slope DIFFERS by n? Check fit slopes systematically.
    print("\n" + "=" * 100, flush=True)
    print("Slope of log|D_AB|^2 vs sum log|x_i-x_j|^2 across n:", flush=True)
    print("=" * 100, flush=True)
    for n in [4, 5, 6, 8, 10, 12, 15]:
        rng = np.random.default_rng(20260603 + n)
        Y = []; X = []
        for _ in range(150):
            pts = random_S2(rng, n)
            ld_AB_2 = log_D_AB_sq(pts)
            sumlog_dsq = sum(2 * math.log(max(np.linalg.norm(pts[i] - pts[j]), 1e-300))
                            for i in range(n) for j in range(i+1, n))
            Y.append(ld_AB_2); X.append(sumlog_dsq)
        Y = np.array(Y); X = np.array(X)
        slope, intercept = np.polyfit(X, Y, 1)
        corr = float(np.corrcoef(X, Y)[0, 1])
        # Possible expected slope: number of pairs n(n-1)/2 weighted somehow
        # Or: -(n-2)? Let's see.
        print(f"  n = {n:>3}: slope = {slope:+.4f}, corr = {corr:.4f},  "
              f"  -(n-2)/n = {-(n-2)/n:+.4f},  -1+1/n = {-1+1/n:+.4f}")


if __name__ == "__main__":
    main()
