"""
PAPER: 066
TITLE: Jensen-polygon-informed multistart on top of T_2 K=1 + universal geometries
STATUS: heuristic; Jensen is intrinsically univariate, multivariate extension via
        univariate slicing

INSPIRATION
===========
Paper 27 (27pandrosion_jensen.pdf): "Jensen's Formula and the Mahler Measure
as Pandrosion Quotient Integrals".

Theorem 2.1 (Jensen's formula, paper 27 eq 2):
    J(R) := (1/2pi) Integral_0^{2pi} log|P(R e^{i theta})| d theta
         = log|a_d| + Sum_k log max(R, |alpha_k|).

Differentiating in log R:
    dJ/d(log R) = #{k : |alpha_k| < R}     (the cumulative root counter).

Equivalently: the graph of J(R) vs log R is piecewise linear with INTEGER
slopes, and BREAKPOINTS at exactly the moduli |alpha_k| of the roots.
This is the "Newton-Jensen polygon".

For multivariate F: C^n -> C^n, Jensen is intrinsically univariate.  We
extend by SLICING: fix all coordinates except one, evaluate F_1 as a
polynomial in z_1, find its roots numerically, collect moduli.  Repeat
over many random fixings and aggregate the modulus histogram.  These
informed moduli become priors for spiral starts (replacing the generic
"R=2, q=0.7" of Strategy B).

EXPECTATION
===========
The Strategy B' multistart of paper 9-mv places D points on the
sphere S^{2n-1} at exponentially-decaying radii.  The Jensen-based
informed prior places starts at the moduli where roots are actually
expected.  In KS, roots concentrate near |z| = 1 (Edelman-Kostlan), so
the informed prior should be close to Strategy B's R=2 default in
practice.  The improvement is expected to be marginal -- we are not
discovering new basins, just placing starts more efficiently.

This script measures the actual gain.  If it's <5pp over flow/062, we
admit the bound is essentially what generic spiral already achieves.
"""

from __future__ import annotations
import cmath, itertools, math, random, time
from itertools import product as iprod

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bench_064_vs_059 import (
    F_eval, eval_poly, is_finite_vec, residual_norm, degree, solve_linear,
    correct_T2K1_portfolio, gen_random_poly_system
)


# ============================================================================
# Univariate Jensen polygon
# ============================================================================
def univariate_coefficients(poly, j, fixed_z, n):
    """Given a multivariate polynomial poly (dict exp -> coeff) and a fixed
    (z_1, ..., z_{j-1}, *, z_{j+1}, ..., z_n), return the coefficients of
    the univariate polynomial in z_j.

    Output: list [c_0, c_1, ..., c_d] where c_k is the coefficient of z_j^k.
    """
    coeffs = {}
    for exp, coeff in poly.items():
        a_j = exp[j]
        # Multiply coeff by product over i != j of fixed_z[i]^exp[i]
        term = complex(coeff)
        for i in range(n):
            if i == j: continue
            term *= fixed_z[i] ** exp[i]
        coeffs[a_j] = coeffs.get(a_j, 0.0 + 0.0j) + term
    if not coeffs:
        return [0.0 + 0.0j]
    max_deg = max(coeffs.keys())
    return [coeffs.get(k, 0.0 + 0.0j) for k in range(max_deg + 1)]


def jensen_integral_numerical(coeffs, R, n_theta=64):
    """Numerically integrate log|P(R e^{i theta})| over [0, 2pi], divided by 2pi."""
    if not coeffs or all(abs(c) < 1e-15 for c in coeffs):
        return float("-inf")
    total = 0.0
    for k in range(n_theta):
        theta = 2 * math.pi * k / n_theta
        z = R * cmath.exp(1j * theta)
        val = sum(coeffs[i] * z**i for i in range(len(coeffs)))
        if abs(val) < 1e-300:
            continue
        total += math.log(abs(val))
    return total / n_theta


def jensen_polygon_breakpoints(coeffs, R_grid):
    """Compute J(R) for each R in R_grid, then identify breakpoints (the slopes
    of the polygon are root counts; breakpoints are root moduli).

    Returns: a list of estimated root moduli (length up to deg(coeffs)).
    """
    J = [jensen_integral_numerical(coeffs, R) for R in R_grid]
    # First-difference in log R to get cumulative count
    counts = []
    for i in range(len(R_grid) - 1):
        dlog = math.log(R_grid[i+1]) - math.log(R_grid[i])
        if dlog < 1e-12:
            counts.append(0)
            continue
        slope = (J[i+1] - J[i]) / dlog
        counts.append(slope)
    # Breakpoints: where slope jumps by ~1
    moduli = []
    prev_count = 0
    for i, c in enumerate(counts):
        if c - prev_count > 0.5:
            # New root(s) entered; estimate modulus by interpolation
            mid = math.sqrt(R_grid[i] * R_grid[i+1])
            jumps = round(c - prev_count)
            for _ in range(jumps):
                moduli.append(mid)
            prev_count = c
    return moduli


def estimate_modulus_distribution(polys, n, n_slices=8, R_grid=None, seed=0):
    """For each coordinate j and several random slicings, build the univariate
    polynomial and compute its Jensen polygon to extract root moduli."""
    if R_grid is None:
        R_grid = [0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0]
    rng = random.Random(seed)
    all_moduli = []
    for j in range(n):
        for slc in range(n_slices):
            # random fixed values for coords != j
            fixed = [complex(rng.gauss(0, 0.5), rng.gauss(0, 0.5)) for _ in range(n)]
            fixed[j] = 0.0 + 0.0j  # placeholder
            # Use F_j as the univariate polynomial in z_j (j-th equation)
            coeffs = univariate_coefficients(polys[j], j, fixed, n)
            if len(coeffs) <= 1:
                continue
            # Find roots numerically (more reliable than Jensen polygon for low deg)
            try:
                # numpy not directly here; use companion matrix root via custom
                roots = numpy_roots(coeffs)
                for r in roots:
                    all_moduli.append(abs(r))
            except Exception:
                pass
    return all_moduli


def numpy_roots(coeffs):
    """Find roots of polynomial coeffs[0] + coeffs[1] z + ... + coeffs[d] z^d
    via companion matrix.  No numpy dependency."""
    d = len(coeffs) - 1
    while d > 0 and abs(coeffs[d]) < 1e-15:
        coeffs = coeffs[:-1]
        d -= 1
    if d <= 0:
        return []
    if d == 1:
        return [-coeffs[0] / coeffs[1]]
    # Companion matrix
    lead = coeffs[d]
    M = [[0.0 + 0.0j] * d for _ in range(d)]
    for i in range(d - 1):
        M[i + 1][i] = 1.0
    for i in range(d):
        M[i][d - 1] = -coeffs[i] / lead
    # Eigenvalues via QR (we use a crude power iteration / Bairstow alternative).
    # For d <= 8 typical, return via repeated deflation using the simplest root finder.
    # We just use a basic algorithm: Durand-Kerner.
    return durand_kerner(coeffs)


def durand_kerner(coeffs, max_iter=300, tol=1e-10):
    """Simultaneous Durand-Kerner iteration for polynomial roots."""
    d = len(coeffs) - 1
    while d > 0 and abs(coeffs[d]) < 1e-15:
        coeffs = coeffs[:-1]; d -= 1
    if d <= 0:
        return []
    if d == 1:
        return [-coeffs[0] / coeffs[1]]
    # Initialize on a circle in complex plane
    base = 0.4 + 0.9j
    roots = [base ** k for k in range(d)]
    def P(z):
        s = 0.0 + 0.0j; zk = 1.0
        for c in coeffs:
            s += c * zk; zk *= z
        return s
    for it in range(max_iter):
        max_change = 0.0
        new_roots = list(roots)
        for k in range(d):
            r_k = roots[k]
            denom = coeffs[d]
            for j in range(d):
                if j == k: continue
                denom *= (r_k - roots[j])
            if abs(denom) < 1e-300: continue
            delta = P(r_k) / denom
            new_roots[k] = r_k - delta
            max_change = max(max_change, abs(delta))
        roots = new_roots
        if max_change < tol:
            break
    return roots


# ============================================================================
# Informed multistart with Jensen-derived moduli
# ============================================================================
def gen_starts_jensen_informed(n, count, target_moduli, seed=20260427):
    """Generate starts at radii sampled from the Jensen-derived modulus
    distribution.  If empty, fall back to Strategy B' (R=2)."""
    if not target_moduli:
        return gen_starts_default(n, count, seed=seed)
    # Bin moduli; sample radii from this distribution
    rng = random.Random(seed)
    # filter outliers
    sorted_mods = sorted(m for m in target_moduli if 0.01 < m < 20.0)
    if not sorted_mods:
        return gen_starts_default(n, count, seed=seed)
    starts = []
    for k in range(count):
        # Pick a modulus from the empirical distribution (with small jitter)
        rho_base = sorted_mods[rng.randint(0, len(sorted_mods) - 1)]
        rho = rho_base * math.exp(rng.gauss(0, 0.1))
        # Random direction on S^{2n-1}
        u = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
        norm = sum(abs(v) ** 2 for v in u) ** 0.5 or 1.0
        u = [v / norm for v in u]
        starts.append([rho * v for v in u])
    return starts


def gen_starts_default(n, count, R=2.0, q=0.7, seed=20260427):
    """Strategy B' default (paper 9-mv)."""
    rng = random.Random(seed); starts = []
    for k in range(count):
        u = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
        norm = sum(abs(v)**2 for v in u) ** 0.5 or 1.0
        u = [v / norm for v in u]
        rho = R * q ** k
        if rho < 0.05: rho = 0.05 + 0.5 * rng.random()
        starts.append([rho * v for v in u])
    return starts


def is_new_root(z, found, tol=1e-4):
    return all(max(abs(z[i] - r[i]) for i in range(len(z))) > tol for r in found)


def F_eval_deflated(polys, z, found_roots, sigma=0.5):
    F = F_eval(polys, z)
    if not found_roots: return F
    penalty = 1.0
    for r in found_roots:
        d2 = sum(abs(z[i] - r[i])**2 for i in range(len(z)))
        penalty *= (1.0 - math.exp(-d2 / (sigma**2)))
    if penalty < 1e-10:
        return [complex("inf")] * len(F)
    return [v / penalty for v in F]


def newton_ELS_deflated(polys, z, found_roots, sigma=0.5):
    n = len(z); F_z = F_eval_deflated(polys, z, found_roots, sigma)
    if not is_finite_vec(F_z): return list(z), False
    h = 1e-6
    J = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        zp = list(z); zp[j] += h
        F_p = F_eval_deflated(polys, zp, found_roots, sigma)
        for i in range(n):
            J[i][j] = (F_p[i] - F_z[i]) / h
    direction = solve_linear(J, [-v for v in F_z])
    if direction is None: return list(z), False
    r0 = sum(abs(v)**2 for v in F_z); best_z, best_r = list(z), r0
    for k in range(-6, 7):
        tau = 2.0 ** k
        cand = [z[i] + tau * direction[i] for i in range(n)]
        if not is_finite_vec(cand): continue
        F_c = F_eval_deflated(polys, cand, found_roots, sigma)
        r = sum(abs(v)**2 for v in F_c) if is_finite_vec(cand) else float("inf")
        if r < best_r: best_z, best_r = cand, r
    return best_z, best_r < r0


def deflated_orbit(polys, z_init, found_roots, sigma=0.4, tol=1e-9, max_epochs=60):
    z = list(z_init)
    for _ in range(max_epochs):
        if residual_norm(polys, z) < tol: return z, True
        z_new, ok = newton_ELS_deflated(polys, z, found_roots, sigma)
        if not ok: return z, False
        if max(abs(z_new[i] - z[i]) for i in range(len(z))) < 1e-12:
            return z, residual_norm(polys, z) < tol
        z = z_new
    return z, residual_norm(polys, z) < tol


def bezout_estimate(polys):
    out = 1
    for p in polys:
        out *= max(degree(p), 1)
    return out


# ============================================================================
# Coverage with Jensen-informed starts
# ============================================================================
def full_coverage_jensen(polys, n_factor=3, n_min=12, tol=1e-9):
    n = len(next(iter(polys[0])))
    bez = bezout_estimate(polys)
    # Jensen-derived modulus prior
    moduli = estimate_modulus_distribution(polys, n, n_slices=6)
    starts = gen_starts_jensen_informed(n, max(n_min, n_factor * bez), moduli)
    found = []; modes_used = {"path": 0, "simplex": 0, "sparse": 0,
                              "cube": 0, "deflation": 0}
    for z0 in starts:
        z, ok = correct_T2K1_portfolio(polys, z0, tol=tol, max_epochs=30)
        if ok and is_new_root(z, found):
            found.append(z); modes_used["path"] += 1
    # Deflation pass
    if len(found) < bez:
        rng = random.Random(20260428)
        for _ in range(2 * bez):
            if len(found) >= bez: break
            # pick random direction at a Jensen-informed radius
            if moduli:
                rho_base = moduli[rng.randint(0, len(moduli) - 1)]
                rho = rho_base * math.exp(rng.gauss(0, 0.1))
            else:
                rho = 1.0 + rng.gauss(0, 0.3)
            u = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
            nrm = sum(abs(v)**2 for v in u)**0.5 or 1.0
            z0 = [rho * v / nrm for v in u]
            z, ok = deflated_orbit(polys, z0, found, sigma=0.4, tol=tol)
            if ok and is_new_root(z, found):
                found.append(z); modes_used["deflation"] += 1
    return {"roots": found, "bezout": bez,
            "coverage": len(found) / max(bez, 1),
            "n_moduli": len(moduli), "modes": modes_used}


def main():
    print("=" * 116)
    print("flow/066 -- Jensen-polygon-informed multistart + T_2 K=1 + 4 geometries")
    print("=" * 116)
    print("Use paper 27's Jensen integral to estimate root moduli on univariate slices")
    print("of the multivariate system.  Place starts at those moduli (informed prior).")
    print()
    cases = [(2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 2), (3, 3), (4, 2)]
    print(f"  {'(n,d)':>8} {'Bez':>5} | {'roots':>5} {'cov%':>6} "
          f"{'#moduli':>8} {'modes':>30} {'time':>7}", flush=True)
    print("-" * 116)
    for n, d in cases:
        polys = gen_random_poly_system(n, d, seed=61000 + 100 * n + d)
        t0 = time.time()
        r = full_coverage_jensen(polys, n_factor=3, n_min=12)
        elapsed = time.time() - t0
        modes = " ".join(f"{k}:{v}" for k, v in r["modes"].items() if v > 0)
        print(f"  ({n:>2},{d:>2}) {r['bezout']:>5} | "
              f"{len(r['roots']):>5} {100*r['coverage']:>5.1f}% "
              f"{r['n_moduli']:>8} {modes:>30} {elapsed:>6.1f}s",
              flush=True)
    print()
    print("=" * 116)
    print("VERDICT")
    print("=" * 116)
    print("  Jensen polygon (paper 27) gives moduli on univariate slices.")
    print("  For Kostlan-Smale, roots concentrate near |z|=1, so informed prior")
    print("  is similar to default Strategy B' R=2.  Expected gain: small.")
    print("=" * 116)


if __name__ == "__main__":
    main()
