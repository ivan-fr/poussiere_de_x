"""SU(2)-equivariant barrier method attempt for Atiyah-Sutcliffe.

Inspired by Marcus-Spielman-Srivastava (paper 79):
  - Define an "Atiyah barrier" Phi(t; x_1, ..., x_n) such that bounding the
    smallest root of Phi gives bound on smallest eigenvalue of G_norm.
  - Find a randomization sigma over which the EXPECTED barrier is real-stable.
  - Use real-stability + interlacing to get worst-case bound.

Key obstacle for Atiyah:
  - G_norm is SU(2)-invariant — global SU(2) rotation doesn't change it.
  - So we need a DIFFERENT randomization.

Three candidates tested here:
  (R1) Random pair-swaps: sigma in S_n acts by permuting points.
       (G transforms by P sigma G P sigma^T, det invariant.)
  (R2) Random reflections: epsilon_i in {+1, -1} acts by x_i -> epsilon_i x_i.
       (NOT SU(2) action, but moves us through configuration space.)
  (R3) Random small perturbations: x_i -> R_i x_i with R_i random.
       (Does change G but stays nearby.)

Test for each: is the EXPECTED characteristic polynomial mu(t) = E[det(tI-G)]
real-rooted? If yes, MSS-style argument applies.
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


def antipodal_split_S2(rng, n):
    v = rng.standard_normal((n, 3)) * 0.1
    half = n // 2
    v[:half, 2] += 1.0
    v[half:, 2] -= 1.0
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def atiyah_polys(points):
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


def atiyah_M(points):
    n = len(points)
    polys = atiyah_polys(points)
    M = np.zeros((n, n), dtype=complex)
    for i, p in enumerate(polys):
        if len(p) < n: pp = np.zeros(n, dtype=complex); pp[:len(p)] = p; p = pp
        elif len(p) > n: p = p[:n]
        M[i] = p
    return M


def gram_norm(points):
    n = len(points)
    M = atiyah_M(points)
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    W = np.diag(1.0 / binom)
    G = M @ W @ M.conj().T
    d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
    return G / np.outer(d, d)


def char_poly(G):
    """Return coefficients of det(tI - G), high to low order."""
    n = G.shape[0]
    # Compute via eigvals (PSD so real); then expand poly
    w = np.real(np.linalg.eigvalsh(G))
    # poly = prod (t - w_i)
    coeffs = np.poly(w)  # numpy convention: high to low
    return coeffs, w


def is_real_rooted(coeffs, tol=1e-8):
    """Check if poly with given coeffs has all real roots."""
    roots = np.roots(coeffs)
    return all(abs(r.imag) < tol * (1 + abs(r.real)) for r in roots), roots


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("SU(2) BARRIER METHOD: test if E[det(tI-G)] is real-rooted under various randomizations",
          flush=True)
    print("=" * 100, flush=True)

    n = 6
    rng = np.random.default_rng(20260701)
    base_pts = antipodal_split_S2(rng, n)
    G_base = gram_norm(base_pts)
    base_coeffs, base_eigs = char_poly(G_base)

    print(f"\nBase configuration (n=6 antipodal):")
    print(f"  eigenvalues of G_base: {[f'{e:.3e}' for e in base_eigs]}")
    print(f"  -log det G_base = {-np.sum(np.log(np.maximum(base_eigs, 1e-300))):.4f}")

    # Approach R2: epsilon-randomization (x_i -> epsilon_i * x_i? or x_i -> -x_i?)
    print(f"\nApproach R2 (point-by-point reflection):")
    n_samples = 2 ** n
    avg_coeffs = np.zeros(n + 1)
    for s in range(n_samples):
        eps = np.array([1 if (s >> i) & 1 == 0 else -1 for i in range(n)])
        pts = base_pts * eps[:, None]
        G = gram_norm(pts)
        c, _ = char_poly(G)
        avg_coeffs += np.real(c)
    avg_coeffs /= n_samples
    is_real, roots = is_real_rooted(avg_coeffs)
    print(f"  avg char poly coeffs: {avg_coeffs}")
    print(f"  roots: {sorted([(r.real, r.imag) for r in roots])}")
    print(f"  all real-rooted? {is_real}")
    if is_real:
        real_roots = sorted([r.real for r in roots])
        print(f"  smallest root: {real_roots[0]:.4e}")
        print(f"  -log of product: {-sum(math.log(max(r, 1e-300)) for r in real_roots):.4f}")

    # Approach R3: small random perturbations of all points
    print(f"\nApproach R3 (small i.i.d. perturbations):")
    rng2 = np.random.default_rng(20260702)
    eigs_list = []
    for trial in range(100):
        perturb = 0.1 * rng2.standard_normal((n, 3))
        new_pts = base_pts + perturb
        new_pts = new_pts / np.linalg.norm(new_pts, axis=1, keepdims=True)
        G = gram_norm(new_pts)
        _, w = char_poly(G)
        eigs_list.append(np.sort(w))
    eigs_arr = np.array(eigs_list)
    avg_eigs = np.mean(eigs_arr, axis=0)
    print(f"  Average eigenvalues over 100 random perturbations:")
    print(f"  {[f'{e:.3e}' for e in avg_eigs]}")
    print(f"  Avg -log det = {-np.sum(np.log(np.maximum(avg_eigs, 1e-300))):.4f}")

    # MSS-style: pick a random "matching" — pair up points and see if Gram has structure
    print(f"\nApproach R4 (random pairing of antipodal points):")
    half = n // 2
    rng3 = np.random.default_rng(20260703)
    # Simulate matching: random orthogonal rotation per pair
    eigs_pair = []
    for trial in range(20):
        pts = antipodal_split_S2(rng3, n)
        G = gram_norm(pts)
        _, w = char_poly(G)
        eigs_pair.append(np.sort(np.real(w)))
    eigs_pair = np.array(eigs_pair)
    print(f"  smallest eigenvalues range: [{eigs_pair[:, 0].min():.3e}, {eigs_pair[:, 0].max():.3e}]")
    print(f"  geometric mean of smallest: {math.exp(np.mean(np.log(np.maximum(eigs_pair[:, 0], 1e-300)))):.3e}")

    # Check key MSS criterion: Hyperbolicity / interlacing of pairs
    print("\n" + "=" * 100, flush=True)
    print("MSS interlacing test: do char polys of nearby configs interlace?", flush=True)
    print("=" * 100, flush=True)
    pts1 = random_S2(np.random.default_rng(1), 6)
    pts2 = random_S2(np.random.default_rng(2), 6)
    # Path between pts1 and pts2
    interlace_found = True
    prev_roots = None
    for t in [0, 0.25, 0.5, 0.75, 1.0]:
        pts_t = (1-t) * pts1 + t * pts2
        pts_t = pts_t / np.linalg.norm(pts_t, axis=1, keepdims=True)
        G = gram_norm(pts_t)
        _, w = char_poly(G)
        w_sorted = np.sort(np.real(w))
        if prev_roots is not None:
            # Interlace: w_i^t should be between w_i and w_{i+1} of prev
            interlaces = all(prev_roots[i] <= w_sorted[i] <= prev_roots[i+1] if i+1 < len(prev_roots) else True
                           for i in range(len(w_sorted)))
        prev_roots = w_sorted
        print(f"  t = {t}: eigs = {[f'{e:.3e}' for e in w_sorted]}")


if __name__ == "__main__":
    main()
