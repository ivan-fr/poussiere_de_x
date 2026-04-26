"""Atiyah-Sutcliffe via the Pandrosion-Vandermonde route.

The alpha-route was refuted (paper 109/103 correction). We try a different
Pandrosion-native invariant: the singular-value gap of the Schmidt slope
matrix Q_F associated to the Atiyah-Berry system.

Specifically, paper 103 Theorem 2.1 reformulates D as
    |D| = |det Q_F(0, zeta_star)| / prod_i ||p_i||_{SU(2)}.
The alpha-route bounds sigma_min(Q_F) via Smale beta. The Pandrosion-native
route bounds sigma_min(Q_F) via the cross-ratio invariants of the Atiyah-Berry
configuration -- specifically, the SU(2)-equivariant Cayley-Menger determinant
of the inter-point distances.

Numerical test: compute sigma_min(Q_F) directly via SVD, compare to log|D|,
and check if there's a clean relationship that bypasses Smale alpha.
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


def atiyah_polynomials(points):
    """Build n Atiyah-Berry polynomials, return as list of coefs (ascending)."""
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


def atiyah_matrix_monomial(points):
    """The Atiyah matrix M with rows = monomial coefficients of p_i."""
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


def normalize_su2(M, n):
    """Divide each row of M by its SU(2)-invariant norm."""
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    norms = np.zeros(n)
    M_norm = M.copy()
    for i in range(n):
        nrm2 = float(np.sum(np.abs(M[i])**2 / binom))
        if nrm2 < 1e-300:
            return M_norm, norms, False
        norms[i] = math.sqrt(nrm2)
        M_norm[i] = M[i] / norms[i]
    return M_norm, norms, True


def compute_atiyah_invariants(points):
    """Compute multiple Atiyah invariants for comparison.

    Returns dict with:
      log_D: Atiyah determinant in log scale
      sigma_min_M: smallest singular value of M (monomial basis, unnormalized)
      sigma_min_M_norm: smallest singular value of M (SU(2)-normalized rows)
      cond_M: condition number of M
      log_norms: log of SU(2) norms
    """
    n = len(points)
    M = atiyah_matrix_monomial(points)

    # SU(2) norms
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    log_norms = np.zeros(n)
    for i in range(n):
        nrm2 = float(np.sum(np.abs(M[i])**2 / binom))
        if nrm2 < 1e-300:
            return None
        log_norms[i] = 0.5 * math.log(nrm2)

    # log|D| = log|det M| - sum log_norms
    sign, logabs = np.linalg.slogdet(M)
    if not np.isfinite(logabs):
        return None
    log_D = logabs - log_norms.sum()

    # Singular values
    s = np.linalg.svd(M, compute_uv=False)
    sigma_min = float(s[-1])
    sigma_max = float(s[0])

    # Normalized version
    M_norm, _, ok = normalize_su2(M, n)
    if not ok:
        return None
    s_norm = np.linalg.svd(M_norm, compute_uv=False)
    sigma_min_norm = float(s_norm[-1])
    sigma_max_norm = float(s_norm[0])

    return dict(
        log_D=log_D,
        sigma_min=sigma_min,
        sigma_max=sigma_max,
        sigma_min_norm=sigma_min_norm,
        sigma_max_norm=sigma_max_norm,
        log_norms=log_norms,
    )


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 90, flush=True)
    print("Atiyah via Pandrosion-Vandermonde route", flush=True)
    print("Test: relate log|D| to sigma_min(M_norm)?", flush=True)
    print("=" * 90, flush=True)

    print(f"\n{'n':>4} {'#configs':>9} {'min log|D|':>11} "
          f"{'min sigma_min_norm':>19} {'predicted':>12} {'gap':>8}",
          flush=True)
    print("-" * 90, flush=True)
    for n in [5, 10, 15, 20, 30]:
        rng = np.random.default_rng(20260428 + n)
        n_cfg = max(50, 500 // n)
        log_Ds = []
        sigma_norms = []
        for _ in range(n_cfg):
            pts = random_S2(rng, n)
            inv = compute_atiyah_invariants(pts)
            if inv is None: continue
            log_Ds.append(inv['log_D'])
            sigma_norms.append(inv['sigma_min_norm'])
        log_Ds = np.array(log_Ds)
        sigma_norms = np.array(sigma_norms)
        # In ON-Hadamard: log|D| = sum log s_i, so n * log(sigma_min_norm) is a lower bound
        predicted = n * np.log(sigma_norms)
        gap = log_Ds - predicted
        print(f"{n:>4} {n_cfg:>9} {float(log_Ds.min()):>11.3f} "
              f"{float(sigma_norms.min()):>19.6f} {float(predicted.min()):>12.3f} "
              f"{float(gap.min()):>8.3f}",
              flush=True)

    # Test 2: relationship sigma_min_norm and Cayley-Menger of points
    print("\n" + "=" * 90, flush=True)
    print("Test 2: Cayley-Menger style — does sigma_min_norm relate to "
          "min_{i,j} |x_i - x_j|?", flush=True)
    print("=" * 90, flush=True)
    print(f"{'n':>4} {'#cfg':>5} {'corr(sigma, min_dist)':>22}", flush=True)
    for n in [5, 10, 15, 20]:
        rng = np.random.default_rng(20260429 + n)
        sigmas = []
        min_dists = []
        for _ in range(100):
            pts = random_S2(rng, n)
            inv = compute_atiyah_invariants(pts)
            if inv is None: continue
            # min pairwise distance
            md = float('inf')
            for i in range(n):
                for j in range(i+1, n):
                    d = float(np.linalg.norm(pts[i] - pts[j]))
                    if d < md: md = d
            sigmas.append(inv['sigma_min_norm'])
            min_dists.append(md)
        sigmas = np.array(sigmas); min_dists = np.array(min_dists)
        if len(sigmas) > 5:
            corr = float(np.corrcoef(sigmas, min_dists)[0, 1])
            print(f"{n:>4} {len(sigmas):>5} {corr:>22.4f}", flush=True)

    # Test 3: a Pandrosion-Vandermonde candidate invariant
    # The minimum singular value of the un-normalized M (monomial basis)
    print("\n" + "=" * 90, flush=True)
    print("Test 3: bound |D| via sigma_min(M) / prod(SU(2) norms) directly?",
          flush=True)
    print("=" * 90, flush=True)
    print(f"{'n':>4} {'#cfg':>5} {'min log|D|':>12} "
          f"{'log sigma_min':>14} {'-sum log norms':>15}", flush=True)
    for n in [5, 10, 15, 20]:
        rng = np.random.default_rng(20260430 + n)
        records = []
        for _ in range(100):
            pts = random_S2(rng, n)
            inv = compute_atiyah_invariants(pts)
            if inv is None: continue
            records.append((inv['log_D'],
                           math.log(inv['sigma_min']),
                           -inv['log_norms'].sum()))
        records = np.array(records)
        i_worst = int(np.argmin(records[:, 0]))
        log_D = records[i_worst, 0]
        log_smin = records[i_worst, 1]
        neg_log_n = records[i_worst, 2]
        print(f"{n:>4} {len(records):>5} {log_D:>12.4f} "
              f"{log_smin:>14.4f} {neg_log_n:>15.4f}", flush=True)


if __name__ == "__main__":
    main()
