"""Test the IDENTITY: |D|^2 = det(G_norm) * prod_k C(n-1, k)

Derivation:
  G = M W M^* with W = diag(1/C(n-1,k))   [SU(2) inner-product Gram of rows of M]
  det G = |det M|^2 * det W = |det M|^2 / prod_k C(n-1, k)
  G_norm = D^{-1/2} G D^{-1/2} with D = diag(||p_i||^2_SU(2))
  det G_norm = det G / prod_i ||p_i||^2 = |det M|^2 / (prod C(n-1,k) * prod ||p_i||^2)
  And |D| = |det M| / prod ||p_i||, so |D|^2 = |det M|^2 / prod ||p_i||^2.

Therefore:  det G_norm = |D|^2 / prod_k C(n-1, k).
Equivalently:  log det G_norm + sum_k log C(n-1,k) = 2 log|D|.

If this identity HOLDS exactly, Atiyah-Sutcliffe (log|D| >= 0) becomes:
  log det G_norm >= -sum_k log C(n-1, k) =: -L_n

and Hadamard gives log det G_norm <= 0 (since G_norm is hermitian PSD with unit diagonal).

So Atiyah-Sutcliffe becomes a HADAMARD-DEFICIT bound on the Gram matrix.
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


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 95, flush=True)
    print("IDENTITY TEST: 2*log|D| = log det G_norm + sum_k log C(n-1,k)?",
          flush=True)
    print("=" * 95, flush=True)
    print(f"{'n':>4} {'#cfg':>5} {'L_n=sum log C':>15} {'max |LHS-RHS|':>15}", flush=True)
    print("-" * 95, flush=True)
    for n in [4, 5, 6, 8, 10, 12, 15, 20]:
        rng = np.random.default_rng(20260502 + n)
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        max_err = 0.0
        for _ in range(80):
            pts = random_S2(rng, n)
            M = atiyah_matrix(pts)
            log_norms_2 = np.log(np.maximum(
                np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
            sign, logabs = np.linalg.slogdet(M)
            if not np.isfinite(logabs): continue
            log_D = logabs - 0.5 * log_norms_2.sum()
            # G = M W M^*
            W = np.diag(1.0 / binom)
            G = M @ W @ M.conj().T
            d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
            G_norm = G / np.outer(d, d)
            sg, ldG = np.linalg.slogdet(G_norm)
            if not np.isfinite(ldG): continue
            lhs = 2.0 * log_D
            rhs = ldG + L_n
            err = abs(lhs - rhs)
            if err > max_err: max_err = err
        print(f"{n:>4} {80:>5} {L_n:>15.4f} {max_err:>15.2e}", flush=True)

    # Now test: is the Hadamard deficit bound -L_n actually attained?
    # We compare: log det G_norm vs -L_n (the Atiyah threshold).
    print("\n" + "=" * 95, flush=True)
    print("MARGIN: log det G_norm >= -L_n? (Atiyah-Sutcliffe reformulated)",
          flush=True)
    print("=" * 95, flush=True)
    print(f"{'n':>4} {'-L_n (threshold)':>18} {'min log det G_norm':>22} "
          f"{'margin (>= 0?)':>16}", flush=True)
    for n in [4, 5, 6, 8, 10, 12, 15, 20, 25, 30]:
        rng = np.random.default_rng(20260503 + n)
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        ld_Gs = []
        for _ in range(100):
            pts = random_S2(rng, n)
            M = atiyah_matrix(pts)
            W = np.diag(1.0 / binom)
            G = M @ W @ M.conj().T
            d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
            G_norm = G / np.outer(d, d)
            sg, ldG = np.linalg.slogdet(G_norm)
            if np.isfinite(ldG): ld_Gs.append(ldG)
        if ld_Gs:
            min_ld = float(min(ld_Gs))
            margin = min_ld - (-L_n)
            print(f"{n:>4} {-L_n:>18.4f} {min_ld:>22.4f} {margin:>16.4f}",
                  flush=True)


if __name__ == "__main__":
    main()
