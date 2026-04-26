"""Atiyah-Sutcliffe via the FULL singular spectrum.

Lessons from paper 109:
  - sigma_min alone fails (one term in sum_i log sigma_i)
  - alpha-route fails (per-polynomial)
  - D is a JOINT invariant: log|D| = sum_i log sigma_i(M) - sum_i log ||p_i||_SU(2)

New strategy: instead of bounding sigma_min, bound the *Frobenius* sum directly.
By AM-GM and the trace identity:
  sum_i log sigma_i^2 = log det(M^* M)
  sum_i sigma_i^2 = ||M||_F^2 = sum_{i,j} |M_{ij}|^2

For SU(2)-normalized M_norm:
  ||M_norm||_F^2 = sum_i ||p_i||_2^2 / ||p_i||_SU(2)^2
                 = sum_i (sum_k |M_ik|^2) / (sum_k |M_ik|^2 / C(n-1,k))

So row i contributes a ratio between Euclidean and SU(2)-weighted norms.
The MGM-AM inequality gives:
  log|det M_norm|^2 = sum log sigma_i^2 <= n log( ||M_norm||_F^2 / n )
NO wait that's an UPPER bound. For the LOWER we need different inequality.

Right idea: Hadamard's inequality says |det M|^2 <= prod_i ||row_i||_2^2.
For the LOWER bound, we use that the MIN over configurations is bounded:

  log|D|^2 = log det(M_norm M_norm^*) = sum log sigma_i^2

If we can show sum log sigma_i^2 >= 0, equivalently det(M_norm M_norm^*) >= 1,
that's exactly Atiyah-Sutcliffe (squared).

KEY CANDIDATE INVARIANT: the Gram matrix G = M_norm M_norm^* has 1's on
diagonal (rows are SU(2)-unit, but in SU(2)-norm not Euclidean!). Wait,
SU(2)-norm uses binomial weights so even diagonal is NOT 1 in Euclidean.

Let G_SU2 := the matrix whose (i,j) entry is the SU(2) inner product
<p_i, p_j>_SU(2) / (||p_i|| ||p_j||). This IS hermitian with diagonal 1.

NEW IDEA: |D| = |det G_SU2|^{1/2}. We test Atiyah-Sutcliffe via the
SU(2)-Gram matrix directly. The diagonal is 1, off-diagonal entries are
the SU(2)-cosines between the Atiyah polynomials — which depend ONLY on
the cross-ratio data of the configuration.

The conjecture becomes: det(G_SU2) >= 1 for all configurations.
But det of hermitian PSD matrix with 1 on diagonal satisfies det <= 1
by Hadamard. So we'd need det = 1 exactly, i.e. orthogonal. False.

I'm confused on the convention. Let me just COMPUTE both
det(M_norm) and det(G_SU2) and see what's >= 1.
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


def su2_inner_product_matrix(M, n):
    """Compute the n x n SU(2)-inner-product Gram matrix of the rows of M.

    SU(2)-inner: <p, q>_SU(2) = sum_k p_k bar(q_k) / C(n-1, k).
    """
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    W = np.diag(1.0 / binom)
    return M @ W @ M.conj().T


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 90, flush=True)
    print("Atiyah via the FULL singular spectrum + SU(2) Gram matrix", flush=True)
    print("=" * 90, flush=True)
    print(f"\n{'n':>4} {'#cfg':>5} {'log|D|_min':>12} {'log|D|_med':>12} "
          f"{'sum_log_sigma_min':>18} {'gram_det_min':>13}", flush=True)
    print("-" * 90, flush=True)
    for n in [4, 5, 6, 8, 10, 12, 15, 20]:
        rng = np.random.default_rng(20260428 + n)
        n_cfg = max(50, 800 // n)
        log_Ds = []
        sum_log_s = []
        gram_dets = []
        for _ in range(n_cfg):
            pts = random_S2(rng, n)
            M = atiyah_matrix(pts)
            binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
            log_norms = 0.5 * np.log(np.maximum(
                np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
            sign, logabs = np.linalg.slogdet(M)
            if not np.isfinite(logabs): continue
            log_D = logabs - log_norms.sum()
            log_Ds.append(log_D)
            # Full singular spectrum: sum log sigma_i = log|det M|
            s = np.linalg.svd(M, compute_uv=False)
            sum_log_s.append(np.sum(np.log(np.maximum(s, 1e-300))))
            # SU(2) Gram matrix
            G = su2_inner_product_matrix(M, n)
            # Normalize to unit-diagonal
            d = np.sqrt(np.real(np.diag(G)))
            G_norm = G / np.outer(d, d)
            sign_g, log_det_g = np.linalg.slogdet(G_norm)
            if np.isfinite(log_det_g):
                gram_dets.append(log_det_g)
        if not log_Ds: continue
        a = np.array(log_Ds); ss = np.array(sum_log_s); gd = np.array(gram_dets)
        print(f"{n:>4} {len(a):>5} {float(a.min()):>12.4f} "
              f"{float(np.median(a)):>12.4f} {float(ss.min()):>18.4f} "
              f"{float(gd.min()):>13.4f}", flush=True)

    # Now the KEY test: is det(G_norm)^{-1/2} >= 1 i.e., log det(G_norm) <= 0?
    # If YES: |D| = product of poly norms / |det M| ... wait reorganize.
    print("\n" + "=" * 90, flush=True)
    print("KEY: log|D| = - (1/2) log det(G_SU2_normalized)?", flush=True)
    print("=" * 90, flush=True)
    print(f"{'n':>4} {'#cfg':>5} {'log|D|_min':>12} {'-1/2 log det G':>18}",
          flush=True)
    for n in [4, 6, 8, 10, 12, 15]:
        rng = np.random.default_rng(20260501 + n)
        diffs = []
        log_Ds = []
        half_neg_log_det_Gs = []
        for _ in range(50):
            pts = random_S2(rng, n)
            M = atiyah_matrix(pts)
            binom = np.array([math.comb(n-1, k) for k in range(n)], dtype=float)
            log_norms = 0.5 * np.log(np.maximum(
                np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
            sign, logabs = np.linalg.slogdet(M)
            if not np.isfinite(logabs): continue
            log_D = logabs - log_norms.sum()
            G = su2_inner_product_matrix(M, n)
            d = np.sqrt(np.real(np.diag(G)))
            G_norm = G / np.outer(d, d)
            sg, log_det_G = np.linalg.slogdet(G_norm)
            if not np.isfinite(log_det_G): continue
            log_Ds.append(log_D)
            half_neg_log_det_Gs.append(-0.5 * log_det_G)
            diffs.append(log_D - (-0.5 * log_det_G))
        if log_Ds:
            print(f"{n:>4} {len(log_Ds):>5} {float(min(log_Ds)):>12.4f} "
                  f"{float(min(half_neg_log_det_Gs)):>18.4f}", flush=True)
            print(f"      diff range: [{min(diffs):+.6f}, {max(diffs):+.6f}]",
                  flush=True)


if __name__ == "__main__":
    main()
