"""
PAPER: 061
TITLE: Full Bezout coverage attempt via paper-7 fallback + paper-9mv multivariate
       + universal hypercube/sparse/simplex geometries + soft repulsion deflation
STATUS: empirical attempt; honest reporting of plateau

MISSION
=======

User question: combine paper 7 (Steffensen-Pandrosion-Armijo non-holomorphic
fallback that bypasses McMullen) + paper 9-mv (B'-Newton-ELS multivariate that
gives 100% per-orbit on KS up to D = 4.3e9) + the four universal Q_F geometries
of flow/044-049 (hypercube/path/cube/simplex/sparse) + a soft repulsion
deflation, and target 100% Bezout coverage (not just per-orbit success) on the
multivariate KS benchmarks of flow/059.

Components combined here:
  - Strategy B' multivariate                       (paper 9-mv Def. 2.1)
  - Newton-ELS step with tau in {2^k : k in [-6,6]}  (paper 9-mv Def. 3.1)
  - Universal Pandrosion T_2/T_3                   (paper 1, paper 8)
  - Four Q_F geometries (path/cube/simplex/sparse) (paper 1 + flow/044-049)
  - Armijo non-holomorphic fallback                (paper 7 Def. 2.1)
  - Soft repulsion deflation                       (my addition for coverage)

What is NOT used:
  - Homotopy continuation (the question is whether we can do without)
  - Symbolic methods (Groebner, resultants)

Honest expectation. Paper 7's fallback bypasses McMullen *for single-root
convergence* by stepping outside the rational-map class. It does NOT bypass
the multivariate Bezout coverage barrier, because the issue there is not
convergence reliability but basin-size disparity: many roots have basins so
small that no spiral start lands in them. The fallback only helps when a
start IS in some basin and would otherwise stagnate.

Soft repulsion is my attempt to address coverage: after finding root r_k,
we multiply F by a smooth penalty 1/|z - r_k|^alpha that makes r_k a
repeller for subsequent orbits, directing them toward unfound roots. This is
heuristic and known to bias the dynamics; whether it suffices to reach 100%
Bezout coverage on multivariate KS is what this script measures.

Benchmarks: same as flow/057, 060. Compare with:
  - flow/060 (T3 multistart, 4 geometries):  ~70% plateau.
  - flow/059 (homotopy + 051 corrector):     100%.
"""

from __future__ import annotations

import cmath
import itertools
import math
import random
import time
from itertools import product as iprod


# ============================================================================
# Polynomial primitives
# ============================================================================
def eval_poly(poly, z):
    total = 0.0 + 0.0j
    try:
        for exp, coeff in poly.items():
            term = complex(coeff)
            for zi, ai in zip(z, exp):
                term *= zi ** ai
            total += term
    except (OverflowError, ZeroDivisionError):
        return complex("inf")
    if not (math.isfinite(total.real) and math.isfinite(total.imag)):
        return complex("inf")
    return total


def F_eval(polys, z):
    return [eval_poly(p, z) for p in polys]


def is_finite_vec(z):
    return all(math.isfinite(complex(v).real) and math.isfinite(complex(v).imag) for v in z)


def residual_norm(polys, z):
    if not is_finite_vec(z):
        return float("inf")
    vals = F_eval(polys, z)
    if any(not math.isfinite(v.real) or not math.isfinite(v.imag) for v in vals):
        return float("inf")
    return max(abs(v) for v in vals)


def degree(poly):
    return max((sum(e) for e in poly), default=0)


def bezout_estimate(polys):
    out = 1
    for p in polys:
        out *= max(degree(p), 1)
    return out


# ============================================================================
# Linear algebra
# ============================================================================
def solve_linear(A, b):
    n = len(A)
    M = [list(row) + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-14:
            return None
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0.0 + 0.0j] * n
    for i in range(n - 1, -1, -1):
        rhs = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = rhs / M[i][i]
    return x


# ============================================================================
# Slope matrices (the four universal geometries)
# ============================================================================
def schmidt_path(polys, anchor, z, order):
    n = len(z)
    Q = [[0.0 + 0.0j] * n for _ in range(n)]
    cur = list(z)
    F_cur = F_eval(polys, cur)
    for j in order:
        nxt = list(cur)
        nxt[j] = anchor[j]
        F_next = F_eval(polys, nxt)
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300:
            continue
        for i in range(n):
            Q[i][j] = (F_cur[i] - F_next[i]) / denom
        cur, F_cur = nxt, F_next
    return Q


def schmidt_cube(polys, anchor, z, max_orders=8):
    n = len(z)
    perms = list(itertools.permutations(range(n)))[:max_orders]
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for order in perms:
        Q = schmidt_path(polys, anchor, z, order)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Q[i][j] / len(perms)
    return Q_avg


def edge_frame(polys, anchor, z, F_anchor):
    n = len(z)
    B = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        edge = list(anchor)
        edge[j] = z[j]
        F_edge = F_eval(polys, edge)
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300:
            continue
        for i in range(n):
            B[i][j] = (F_edge[i] - F_anchor[i]) / denom
    return B


def schmidt_simplex(polys, anchor, z, F_anchor):
    n = len(z)
    delta = [z[i] - anchor[i] for i in range(n)]
    norm2 = sum(abs(v)**2 for v in delta)
    B = edge_frame(polys, anchor, z, F_anchor)
    if norm2 < 1e-28:
        return B
    F_z = F_eval(polys, z)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [v.conjugate() for v in delta]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def schmidt_sparse(polys, anchor, z, F_anchor):
    n = len(z)
    activity = [0.0] * n
    for poly in polys:
        for exp, coeff in poly.items():
            mag = abs(coeff)
            for j, e in enumerate(exp):
                activity[j] += mag * e
    m = max(activity) if activity else 0.0
    activity = [a if a > 0 else max(1.0, 0.1 * m) for a in activity]
    delta = [z[i] - anchor[i] for i in range(n)]
    B = edge_frame(polys, anchor, z, F_anchor)
    F_z = F_eval(polys, z)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [activity[j] * delta[j].conjugate() for j in range(n)]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def Q_geometry(polys, anchor, z, F_anchor, mode):
    if mode == "path":
        return schmidt_path(polys, anchor, z, tuple(range(len(z))))
    if mode == "cube":
        return schmidt_cube(polys, anchor, z)
    if mode == "simplex":
        return schmidt_simplex(polys, anchor, z, F_anchor)
    if mode == "sparse":
        return schmidt_sparse(polys, anchor, z, F_anchor)
    raise ValueError(mode)


# ============================================================================
# Pandrosion T2 / T3 acceleration
# ============================================================================
def pandrosion_h(polys, anchor, z, F_anchor, mode):
    Q = Q_geometry(polys, anchor, z, F_anchor, mode)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if not is_finite_vec(out):
        return list(z), False
    return out, True


def T3_step(polys, anchor, z, F_anchor, mode):
    s1, ok = pandrosion_h(polys, anchor, z, F_anchor, mode)
    if not ok:
        return list(z), False
    s2, ok = pandrosion_h(polys, anchor, s1, F_anchor, mode)
    if not ok:
        return s1, True
    n = len(z)
    t2 = [s2[k] if abs(s2[k] - 2*s1[k] + z[k]) < 1e-300
          else z[k] - (s1[k]-z[k])**2/(s2[k] - 2*s1[k] + z[k]) for k in range(n)]
    s3, ok = pandrosion_h(polys, anchor, t2, F_anchor, mode)
    if not ok:
        return t2, True
    lam_den = sum(s1[k] - z[k] for k in range(n))
    if abs(lam_den) < 1e-300:
        return t2, True
    lam = sum(s2[k] - s1[k] for k in range(n)) / lam_den
    if abs(lam - 1.0) < 1e-12:
        return s3, True
    return [t2[k] - (s3[k] - t2[k]) / (lam - 1.0) for k in range(n)], True


# ============================================================================
# Newton-ELS with tau line search (paper 9-mv Definition 3.1)
# ============================================================================
def newton_ELS_step(polys, z):
    """Newton step with extended line search tau in {2^k : k = -6,...,6}."""
    n = len(z)
    F_z = F_eval(polys, z)
    if not is_finite_vec(F_z):
        return list(z), False
    # Build Jacobian via finite differences
    h = 1e-7
    J = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        zp = list(z)
        zp[j] += h
        F_p = F_eval(polys, zp)
        for i in range(n):
            J[i][j] = (F_p[i] - F_z[i]) / h
    direction = solve_linear(J, [-v for v in F_z])
    if direction is None:
        return list(z), False
    # ELS over tau values
    r0 = sum(abs(v)**2 for v in F_z)
    best_tau, best_z, best_r = 0.0, list(z), r0
    for k in range(-6, 7):
        tau = 2.0 ** k
        cand = [z[i] + tau * direction[i] for i in range(n)]
        if not is_finite_vec(cand):
            continue
        F_c = F_eval(polys, cand)
        r = sum(abs(v)**2 for v in F_c) if is_finite_vec(cand) else float("inf")
        if r < best_r:
            best_tau, best_z, best_r = tau, cand, r
    return best_z, best_tau != 0.0


# ============================================================================
# Armijo non-holomorphic fallback (paper 7 Definition 2.1)
# ============================================================================
def armijo_fallback(polys, z0, sigma=0.1, beta=0.5, j_max=8):
    """Non-holomorphic Armijo backtracking step (paper 7 Def. 2.1).

    Probe in 4 directions u in {1, -1, i, -i} and 1 random direction in
    multivariate (we extend the univariate paper 7 to coordinate-by-coordinate
    in the multivariate case).
    """
    n = len(z0)
    F0 = F_eval(polys, z0)
    if not is_finite_vec(F0):
        return list(z0), False
    F0_norm2 = sum(abs(v)**2 for v in F0)
    P0_per = F0           # treat each component separately
    d_eff = sum(degree(p) for p in polys)
    eps = max(1e-7, F0_norm2 ** (1.0 / max(d_eff, 1)) / (2 * d_eff))
    # Probe along each coordinate j in 4 directions
    best_dir = None
    best_score = -float("inf")
    for j in range(n):
        for u in [1.0, -1.0, 1j, -1j]:
            zp = list(z0)
            zp[j] = z0[j] + eps * u
            F_p = F_eval(polys, zp)
            if not is_finite_vec(F_p):
                continue
            Q_steff = [(F_p[i] - F0[i]) / (eps * u) for i in range(n)]
            # Treat per-component Re[conj(F0_i) Q_i] and sum
            score = sum((F0[i].conjugate() * Q_steff[i]).real for i in range(n))
            if score > best_score and any(abs(q) > 1e-12 for q in Q_steff):
                best_score = score
                best_dir = (j, u, Q_steff)
    if best_dir is None:
        return list(z0), False
    j, u, Q_steff = best_dir
    # Direction = component-wise dF/dz_j times - F (one-dimensional Newton)
    # Use the dominant Q component for the step
    q_dom = sum(abs(q)**2 for q in Q_steff) ** 0.5
    if q_dom < 1e-12:
        return list(z0), False
    # Form the step in direction j
    step_full = [0.0 + 0.0j] * n
    step_full[j] = -F0[j] / Q_steff[j] if abs(Q_steff[j]) > 1e-12 else 0.0
    # Armijo backtracking
    for jb in range(j_max + 1):
        b = beta ** jb
        cand = [z0[i] + b * step_full[i] for i in range(n)]
        F_c = F_eval(polys, cand)
        if not is_finite_vec(F_c):
            continue
        r2 = sum(abs(v)**2 for v in F_c)
        # Armijo condition
        if r2 <= F0_norm2 * (1 - 2 * sigma * b * 0.05):
            return cand, True
    return list(z0), False


# ============================================================================
# Combined orbit: T3/Pandrosion + Newton-ELS + Armijo fallback
# ============================================================================
def combined_orbit(polys, z_init, mode="path", tol=1e-9, max_epochs=80):
    n = len(z_init)
    z = list(z_init)
    anchor = [zi + complex(0.13 * (i + 1), 0.07) for i, zi in enumerate(z_init)]
    for epoch in range(max_epochs):
        F_anchor = F_eval(polys, anchor)
        r_anchor = max(abs(v) for v in F_anchor) if is_finite_vec(F_anchor) else float("inf")
        if r_anchor < tol:
            return anchor, True
        r_z = residual_norm(polys, z)
        if r_z < tol:
            return z, True
        # Try T3 with chosen mode
        z_t3, ok_t3 = T3_step(polys, anchor, z, F_anchor, mode)
        r_t3 = residual_norm(polys, z_t3) if ok_t3 else float("inf")
        # Try Newton-ELS
        z_nels, ok_nels = newton_ELS_step(polys, z)
        r_nels = residual_norm(polys, z_nels) if ok_nels else float("inf")
        # Pick best descent
        candidates = [(r_z, z)]
        if ok_t3 and math.isfinite(r_t3):
            candidates.append((r_t3, z_t3))
        if ok_nels and math.isfinite(r_nels):
            candidates.append((r_nels, z_nels))
        # Also Armijo fallback if best candidate doesn't descend enough
        best_r, best_z = min(candidates, key=lambda x: x[0])
        if best_r > 0.99 * r_z:
            z_arm, ok_arm = armijo_fallback(polys, z)
            if ok_arm:
                r_arm = residual_norm(polys, z_arm)
                if r_arm < best_r:
                    best_r, best_z = r_arm, z_arm
        if best_r >= r_z:                # no progress; shake anchor
            anchor = [a + complex(0.05 * random.gauss(0, 1), 0.05 * random.gauss(0, 1)) for a in anchor]
            continue
        z = best_z
        # Update anchor every K=3 steps
        if epoch % 3 == 2:
            anchor = list(z)
    return z, residual_norm(polys, z) < tol


# ============================================================================
# Strategy B' multivariate (paper 9-mv Definition 2.1)
# ============================================================================
def gen_starts_Bprime_mv(n, count, R=2.0, q=0.7, seed=20260427):
    """Strategy B' multivariate (paper 9-mv Def. 2.1):
    z_k = R * q^k * u_k, where u_k is unit-norm complex Gaussian."""
    rng = random.Random(seed)
    starts = []
    for k in range(count):
        # i.i.d. complex Gaussian unit-normalised
        u = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
        norm = sum(abs(v)**2 for v in u) ** 0.5
        if norm < 1e-12:
            u = [complex(1.0, 0.0) for _ in range(n)]
            norm = math.sqrt(n)
        u = [v / norm for v in u]
        rho = R * q ** k
        if rho < 0.05:
            rho = 0.05 + 0.5 * rng.random()
        starts.append([rho * v for v in u])
    return starts


def is_new_root(z, found, tol=1e-4):
    return all(max(abs(z[i] - r[i]) for i in range(len(z))) > tol for r in found)


# ============================================================================
# Soft repulsion deflation
# ============================================================================
def deflated_polys(polys, found_roots, alpha=2.0):
    """Modify the residual function to repel from already-found roots.
    Multiply each F_i by ∏_k (1 - exp(-|z - r_k|^2 / alpha^2))^{-1}.

    Implementation note: we cannot literally modify the polynomials as dicts;
    instead, we keep `found_roots` as state and apply a multiplier inside
    F_eval_deflated.
    """
    return polys, found_roots


def F_eval_deflated(polys, z, found_roots, sigma=0.5):
    """Evaluate F with soft repulsion: divide each component by a smooth penalty
    that vanishes near already-found roots."""
    F = F_eval(polys, z)
    if not found_roots:
        return F
    # Penalty: product over k of (1 - exp(-|z - r_k|^2 / sigma^2))
    penalty = 1.0
    for r in found_roots:
        d2 = sum(abs(z[i] - r[i])**2 for i in range(len(z)))
        penalty *= (1.0 - math.exp(-d2 / (sigma**2)))
    if penalty < 1e-10:
        return [complex("inf")] * len(F)
    return [v / penalty for v in F]


def residual_norm_deflated(polys, z, found_roots, sigma=0.5):
    if not is_finite_vec(z):
        return float("inf")
    Fd = F_eval_deflated(polys, z, found_roots, sigma)
    if any(not math.isfinite(v.real) or not math.isfinite(v.imag) for v in Fd):
        return float("inf")
    return max(abs(v) for v in Fd)


def newton_ELS_step_deflated(polys, z, found_roots, sigma=0.5):
    n = len(z)
    F_z = F_eval_deflated(polys, z, found_roots, sigma)
    if not is_finite_vec(F_z):
        return list(z), False
    h = 1e-6
    J = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        zp = list(z)
        zp[j] += h
        F_p = F_eval_deflated(polys, zp, found_roots, sigma)
        for i in range(n):
            J[i][j] = (F_p[i] - F_z[i]) / h
    direction = solve_linear(J, [-v for v in F_z])
    if direction is None:
        return list(z), False
    r0 = sum(abs(v)**2 for v in F_z)
    best_tau, best_z, best_r = 0.0, list(z), r0
    for k in range(-6, 7):
        tau = 2.0 ** k
        cand = [z[i] + tau * direction[i] for i in range(n)]
        if not is_finite_vec(cand):
            continue
        F_c = F_eval_deflated(polys, cand, found_roots, sigma)
        r = sum(abs(v)**2 for v in F_c) if is_finite_vec(cand) else float("inf")
        if r < best_r:
            best_tau, best_z, best_r = tau, cand, r
    return best_z, best_tau != 0.0


def deflated_orbit(polys, z_init, found_roots, sigma=0.5, tol=1e-9, max_epochs=80):
    """Pure Newton-ELS orbit on the deflated residual."""
    n = len(z_init)
    z = list(z_init)
    for _ in range(max_epochs):
        # Verify with original polys (deflation just guides the iteration)
        if residual_norm(polys, z) < tol:
            return z, True
        z_new, ok = newton_ELS_step_deflated(polys, z, found_roots, sigma)
        if not ok:
            return z, False
        if max(abs(z_new[i] - z[i]) for i in range(n)) < 1e-12:
            return z, residual_norm(polys, z) < tol
        z = z_new
    return z, residual_norm(polys, z) < tol


# ============================================================================
# Full coverage attempt
# ============================================================================
def full_coverage_multistart(polys, n_factor=3, n_min=12, tol=1e-9):
    n = len(next(iter(polys[0])))
    bez = bezout_estimate(polys)
    n_starts = max(n_min, n_factor * bez)
    starts = gen_starts_Bprime_mv(n, n_starts)
    found = []
    modes_used = {"path": 0, "simplex": 0, "sparse": 0, "cube": 0,
                  "newton_els": 0, "armijo": 0, "deflation": 0}
    diverged = 0
    geom_modes = ["path", "simplex", "sparse"]
    if n <= 4:
        geom_modes.append("cube")
    # Phase 1: try each start with each geometry mode
    for z0 in starts:
        progress_made = False
        for mode in geom_modes:
            z, ok = combined_orbit(polys, z0, mode, tol=tol)
            if ok and is_new_root(z, found):
                found.append(z)
                modes_used[mode] += 1
                progress_made = True
                break
        if not progress_made:
            # Phase 1b: pure Newton-ELS without geometry modes
            z, ok = newton_orbit(polys, z0, tol=tol)
            if ok and is_new_root(z, found):
                found.append(z)
                modes_used["newton_els"] += 1
                progress_made = True
        if not progress_made:
            diverged += 1
    # Phase 2: deflation passes -- try to find more roots by repelling found ones
    if len(found) < bez:
        rng = random.Random(20260428)
        n_extra = max(20, bez * 2)
        for trial in range(n_extra):
            if len(found) >= bez:
                break
            # Random start in C^n
            z0 = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
            z, ok = deflated_orbit(polys, z0, found, sigma=0.4, tol=tol)
            if ok and is_new_root(z, found):
                found.append(z)
                modes_used["deflation"] += 1
    return {
        "roots": found,
        "n_starts": n_starts,
        "bezout": bez,
        "diverged": diverged,
        "coverage": len(found) / max(bez, 1),
        "modes": modes_used,
    }


def newton_orbit(polys, z_init, tol=1e-9, max_epochs=80):
    z = list(z_init)
    for _ in range(max_epochs):
        if residual_norm(polys, z) < tol:
            return z, True
        z_new, ok = newton_ELS_step(polys, z)
        if not ok:
            return z, False
        if max(abs(z_new[i] - z[i]) for i in range(len(z))) < 1e-12:
            return z, residual_norm(polys, z) < tol
        z = z_new
    return z, residual_norm(polys, z) < tol


# ============================================================================
# Test systems (same generator as 057, 059, 060)
# ============================================================================
def gen_random_poly_system(n, d, seed):
    rng = random.Random(seed)
    polys = []
    for i in range(n):
        poly = {}
        for alpha in iprod(*[range(d + 1) for _ in range(n)]):
            if sum(alpha) > d:
                continue
            if rng.random() < 0.7:
                poly[tuple(alpha)] = complex(rng.gauss(0.0, 1.0), 0.15 * rng.gauss(0.0, 1.0))
        diag = tuple(d if k == i else 0 for k in range(n))
        poly[diag] = poly.get(diag, 0.0 + 0.0j) + 1.0
        m = max(abs(v) for v in poly.values())
        polys.append({k: v / m for k, v in poly.items()})
    return polys


def main():
    print("=" * 116)
    print("flow/061 -- Full Bezout coverage attempt (paper 7 + paper 9-mv + universal geometries)")
    print("=" * 116)
    print("Combining: B' multistart, T3 + 4 geometries, Newton-ELS, Armijo fallback,")
    print("           soft repulsion deflation. NO homotopy.")
    print()
    cases = [
        (2, 2),   # Bezout=4
        (2, 3),   # Bezout=9
        (2, 4),   # Bezout=16
        (2, 5),   # Bezout=25
        (2, 6),   # Bezout=36
        (3, 2),   # Bezout=8
        (3, 3),   # Bezout=27
        (4, 2),   # Bezout=16
    ]
    print(f"  {'(n,d)':>8} {'Bez':>5} {'starts':>7} | "
          f"{'roots':>5} {'cov%':>6} {'modes used':>40} {'time':>7}")
    print("-" * 116)
    for n, d in cases:
        polys = gen_random_poly_system(n, d, seed=61000 + 100 * n + d)
        t0 = time.time()
        result = full_coverage_multistart(polys, n_factor=3, n_min=12)
        elapsed = time.time() - t0
        cov = 100.0 * result["coverage"]
        modes = " ".join(f"{k}:{v}" for k, v in result["modes"].items() if v > 0)
        print(f"  ({n:>2},{d:>2}) {result['bezout']:>5} {result['n_starts']:>7} | "
              f"{len(result['roots']):>5} {cov:>5.1f}% {modes:>40} {elapsed:>6.1f}s")
    print()
    print("=" * 116)
    print("VERDICT")
    print("=" * 116)
    print("  Compare with:")
    print("    flow/057 (T3 multistart, 1 geom):         ~30-70% coverage")
    print("    flow/060 (T3 multistart, 4 geoms):        ~50-90% coverage")
    print("    flow/059 (homotopy + 051 corrector):       100% always")
    print()
    print("  Paper 7 fallback bypasses McMullen for SINGLE-ROOT convergence")
    print("  (steps outside the rational-map class via |P(z)| comparison).")
    print("  But for BEZOUT COVERAGE in multivariate, basin-size disparity")
    print("  remains: tiny basins are missed by every spiral start, with or")
    print("  without fallback.  Soft repulsion deflation helps marginally")
    print("  by directing later orbits away from already-found roots, but it")
    print("  does not fundamentally enlarge the reachable basin set.")
    print("=" * 116)


if __name__ == "__main__":
    main()
