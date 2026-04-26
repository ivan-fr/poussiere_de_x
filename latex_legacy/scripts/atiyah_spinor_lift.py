"""Spinor lift of Atiyah-Sutcliffe.

For each x_i in S^2, lift to a unit spinor u_i in C^2 via the Hopf map:
    x_i = (2 Re(a_i bar(b_i)), 2 Im(a_i bar(b_i)), |a_i|^2 - |b_i|^2)
The lift is well-defined up to U(1) phase per i.

Define:
  V_{k,j} = (u_k^1)^j * (u_k^2)^{n-1-j}  -- spinor Vandermonde
  P_i^Lag(z, w) := prod_{j != i} (a_j w - b_j z) / (a_j b_i - b_j a_i)
                 = prod_{j != i} [z, u_j] / [u_i, u_j]    (normalized to P_i(u_i) = 1)

Classical identity:
  det V = prod_{i<j} [u_i, u_j]      (spinor Vandermonde)
  |det V|^2 = prod_{i<j} |[u_i, u_j]|^2 = 4^{-n(n-1)/2} prod |x_i - x_j|^2

Compare with Atiyah-Berry: do they give the same |D|? Are they related?
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


def hopf_lift(x):
    """Lift x in S^2 to unit spinor u in C^2 (up to U(1) phase).

    Choice: if x = (sin t cos p, sin t sin p, cos t), then
       u = (cos(t/2), sin(t/2) e^{i p})
    """
    z = x[2]
    if z > 1 - 1e-13:
        return np.array([1.0 + 0j, 0.0 + 0j])
    if z < -1 + 1e-13:
        return np.array([0.0 + 0j, 1.0 + 0j])
    cos_half = math.sqrt((1 + z) / 2)
    sin_half = math.sqrt((1 - z) / 2)
    # phi from (x[0], x[1]) = sin(theta) (cos phi, sin phi)
    phi = math.atan2(x[1], x[0])
    return np.array([cos_half + 0j, sin_half * complex(math.cos(phi), math.sin(phi))])


def spinor_skew(u, v):
    """[u, v] = u^1 v^2 - u^2 v^1 (skew product, SU(2)-invariant)."""
    return u[0] * v[1] - u[1] * v[0]


def spinor_vandermonde(spinors):
    """V_{k, j} = (u_k^1)^j * (u_k^2)^{n-1-j}."""
    n = len(spinors)
    V = np.zeros((n, n), dtype=complex)
    for k, u in enumerate(spinors):
        for j in range(n):
            V[k, j] = u[0]**j * u[1]**(n - 1 - j)
    return V


def lagrange_atiyah_polys(spinors):
    """P_i^{Lag}(z, w) = prod_{j != i} (u_j^1 w - u_j^2 z),
    a homogeneous polynomial of degree n-1 in (z, w).

    Returns coefficients in the basis z^k w^{n-1-k}, indexed by k.
    """
    n = len(spinors)
    polys = []
    for i in range(n):
        c = np.array([1.0 + 0j])  # constant 1 = w^0... start as 0-degree
        # For each j != i, multiply by linear factor (-u_j^2 z + u_j^1 w)
        # = (u_j^1) w + (-u_j^2) z
        for j in range(n):
            if j == i: continue
            u = spinors[j]
            # linear factor: -u[1] z + u[0] w
            # In our convention coefficients in basis z^k w^{n-1-k},
            # current poly has degree d_now, we multiply by (-u[1] z + u[0] w),
            # which has coeffs (in z^k w^{1-k}): [u[0], -u[1]] for k=0, k=1.
            new = np.zeros(len(c) + 1, dtype=complex)
            for k_old, coef in enumerate(c):
                new[k_old] += coef * u[0]      # multiply by w (k unchanged, deg+=1)
                new[k_old + 1] += coef * (-u[1])  # multiply by z (k+=1)
            c = new
        polys.append(c)
    # Convert to "monomial basis" matrix: row i = coefficients of P_i in basis z^0 w^{n-1}, ..., z^{n-1} w^0
    M_L = np.array(polys, dtype=complex)
    # Each row should have length n
    return M_L


def su2_norm_sq(coeffs, n):
    """SU(2)-invariant norm squared of poly with given coefficients.
    Coeffs are in basis z^k w^{n-1-k} for k = 0, ..., n-1.
    Norm: sum |c_k|^2 / C(n-1, k).
    """
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    return float(np.sum(np.abs(coeffs)**2 / binom))


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 110, flush=True)
    print("SPINOR LIFT TEST: relate Atiyah-Berry |D| to spinor Vandermonde", flush=True)
    print("=" * 110, flush=True)

    print(f"\n{'n':>3} {'|det V|^2':>14} {'pred 4^(-n(n-1)/2) prod|x-x|^2':>32} "
          f"{'ratio':>10} {'log|D_Lag|':>12} {'log|D_AB|':>12} {'diff':>8}",
          flush=True)
    for n in [3, 4, 5, 6, 8, 10, 12]:
        rng = np.random.default_rng(20260601 + n)
        diffs = []
        log_D_Ls = []
        log_D_ABs = []
        for trial in range(20):
            pts = random_S2(rng, n)
            spinors = [hopf_lift(x) for x in pts]
            # Spinor Vandermonde
            V = spinor_vandermonde(spinors)
            sg, log_det_V_abs = np.linalg.slogdet(V)
            log_det_V_sq = 2 * log_det_V_abs
            # Predicted: |det V|^2 = prod |[u_i, u_j]|^2
            pred_log_det_V_sq = 0.0
            for i in range(n):
                for j in range(i+1, n):
                    s = spinor_skew(spinors[i], spinors[j])
                    pred_log_det_V_sq += 2 * math.log(max(abs(s), 1e-300))
            ratio = math.exp(log_det_V_sq - pred_log_det_V_sq)

            # Lagrange |D_L|^2 = |det M_L|^2 / prod ||P_i^Lag||^2 = ?
            M_L = lagrange_atiyah_polys(spinors)
            sg2, ld2 = np.linalg.slogdet(M_L)
            log_norms_2 = []
            for i in range(n):
                lns = math.log(max(su2_norm_sq(M_L[i], n), 1e-300))
                log_norms_2.append(lns)
            log_D_L = ld2 - 0.5 * sum(log_norms_2)
            log_D_Ls.append(log_D_L)

            # Atiyah-Berry |D|^2 (paper 110 method)
            n_ = n
            polys_AB = []
            for i in range(n):
                finite = []
                n_inf = 0
                for j in range(n):
                    if i == j: continue
                    v = pts[j] - pts[i]
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
                polys_AB.append(c)
            M_AB = np.zeros((n, n), dtype=complex)
            for i, p in enumerate(polys_AB):
                if len(p) < n:
                    pp = np.zeros(n, dtype=complex)
                    pp[:len(p)] = p
                    p = pp
                elif len(p) > n:
                    p = p[:n]
                M_AB[i] = p
            log_norms_AB = []
            for i in range(n):
                log_norms_AB.append(math.log(max(su2_norm_sq(M_AB[i], n), 1e-300)))
            sg3, ld3 = np.linalg.slogdet(M_AB)
            log_D_AB = ld3 - 0.5 * sum(log_norms_AB)
            log_D_ABs.append(log_D_AB)

            diffs.append(abs(ratio - 1.0))

        # Test claim
        print(f"{n:>3} {math.exp(log_det_V_sq):>14.4e} {math.exp(pred_log_det_V_sq):>32.4e} "
              f"{ratio:>10.6f} {np.mean(log_D_Ls):>12.4f} {np.mean(log_D_ABs):>12.4f} "
              f"{np.mean(np.array(log_D_Ls) - np.array(log_D_ABs)):>8.4f}",
              flush=True)

    print("\n" + "=" * 110, flush=True)
    print("Specific test: is log|D_AB| - log|D_Lag| a function of pairwise distances only?",
          flush=True)
    print("=" * 110, flush=True)
    n = 6
    rng = np.random.default_rng(424242)
    rows = []
    for _ in range(50):
        pts = random_S2(rng, n)
        spinors = [hopf_lift(x) for x in pts]
        # log_D_L
        M_L = lagrange_atiyah_polys(spinors)
        sg, ld = np.linalg.slogdet(M_L)
        log_norms_L = sum(math.log(max(su2_norm_sq(M_L[i], n), 1e-300)) for i in range(n))
        log_D_L = ld - 0.5 * log_norms_L

        # log_D_AB
        polys_AB = []
        for i in range(n):
            finite = []; n_inf = 0
            for j in range(n):
                if i == j: continue
                v = pts[j] - pts[i]; v = v / np.linalg.norm(v)
                a = stereo(v)
                if math.isinf(a.real) or math.isinf(a.imag): n_inf += 1
                else: finite.append(a)
            c = np.array([1.0+0j])
            for a in finite: c = np.convolve(c, np.array([-a, 1.0+0j]))
            if n_inf > 0:
                cc = np.zeros(len(c)+n_inf, dtype=complex); cc[n_inf:] = c; c = cc
            polys_AB.append(c)
        M_AB = np.zeros((n, n), dtype=complex)
        for i, p in enumerate(polys_AB):
            if len(p) < n: pp = np.zeros(n, dtype=complex); pp[:len(p)] = p; p = pp
            elif len(p) > n: p = p[:n]
            M_AB[i] = p
        sg3, ld3 = np.linalg.slogdet(M_AB)
        log_norms_AB = sum(math.log(max(su2_norm_sq(M_AB[i], n), 1e-300)) for i in range(n))
        log_D_AB = ld3 - 0.5 * log_norms_AB

        diff = log_D_AB - log_D_L
        sumlog_dist = sum(math.log(np.linalg.norm(pts[i] - pts[j]))
                         for i in range(n) for j in range(i+1, n))
        rows.append((log_D_L, log_D_AB, diff, sumlog_dist))
    rows = np.array(rows)
    print(f"\n  log|D_AB| - log|D_Lag| vs sum log|x_i - x_j|:")
    print(f"  range of diff: [{rows[:, 2].min():+.4f}, {rows[:, 2].max():+.4f}]")
    print(f"  std of diff: {rows[:, 2].std():.4f}")
    print(f"  correlation diff vs sumlog: {np.corrcoef(rows[:, 2], rows[:, 3])[0,1]:.4f}")
    if rows[:, 2].std() < 0.01:
        print("  ** DIFFERENCE IS CONSTANT → Atiyah-Berry = Lagrange + constant! **")


if __name__ == "__main__":
    main()
