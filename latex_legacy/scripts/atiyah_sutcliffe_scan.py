"""Compute the Atiyah-Sutcliffe determinant |D| for random point configurations
on S^2 in R^3, and scan for counterexamples (|D| < 1).

Construction:
  For n distinct points x_1, ..., x_n in R^3:
  - For each i, compute n-1 unit vectors u_{ij} = (x_j - x_i)/|x_j - x_i|.
  - Stereographic project u_{ij} to alpha_{ij} in C cup {infty}.
  - Build polynomial p_i(z) = z^{n_inf} prod_{finite} (z - alpha_{ij}), of degree n-1.
  - Normalize by the SU(2)-invariant norm:
        ||p||^2 = sum_k |c_k|^2 / binom(n-1, k).
  - Form n x n matrix M with row i = normalized coefficient vector of p_i.
  - D = |det M|.

Sendov / Atiyah-Sutcliffe conjecture: |D| >= 1 always, with equality iff collinear.
"""
from __future__ import annotations
import math, time
import numpy as np


def stereo(v):
    """Stereographic projection from north pole. v in S^2 -> C cup {inf}."""
    z = v[2]
    if z >= 1.0 - 1e-13:
        return complex('inf')
    return complex(v[0], v[1]) / (1.0 - z)


def atiyah_sutcliffe_logD(points):
    """Compute log|D| for n points using Atiyah's original normalization
    (monomial basis):
        D := |det M_monomial| / prod_i ||p_i||_{SU(2)}
    where M_monomial has rows = coefficient vectors of p_i in the standard
    monomial basis (1, t, t^2, ..., t^{n-1}), and ||p||^2_{SU(2)} = sum_k |c_k|^2 / binom(n-1, k).

    The Atiyah-Sutcliffe conjecture asserts log|D| >= 0, with equality iff collinear.

    Returns log|D| (in nats).  Use log-scale to handle very large or very small values.
    """
    n = points.shape[0]
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)

    M = np.zeros((n, n), dtype=complex)
    log_norms = np.zeros(n)
    for i in range(n):
        finite_roots = []
        n_inf = 0
        for j in range(n):
            if i == j:
                continue
            v = points[j] - points[i]
            nrm = np.linalg.norm(v)
            if nrm < 1e-15:
                u = v / max(nrm, 1e-15)
            else:
                u = v / nrm
            a = stereo(u)
            if math.isinf(a.real) or math.isinf(a.imag) or math.isnan(a.real):
                n_inf += 1
            else:
                finite_roots.append(a)
        c = np.array([1.0 + 0j])
        for a in finite_roots:
            c = np.convolve(c, np.array([-a, 1.0 + 0j]))
        if n_inf > 0:
            cc = np.zeros(len(c) + n_inf, dtype=complex)
            cc[n_inf:] = c
            c = cc
        if len(c) < n:
            cc = np.zeros(n, dtype=complex)
            cc[:len(c)] = c
            c = cc
        elif len(c) > n:
            c = c[:n]
        # SU(2)-invariant norm
        norm2_su2 = float(np.sum(np.abs(c)**2 / binom))
        if norm2_su2 < 1e-300:
            return -float('inf')
        log_norms[i] = 0.5 * math.log(norm2_su2)
        # Use monomial-basis row directly (unnormalized)
        M[i] = c

    # log|det M| via SLOG
    sign, logabsdet = np.linalg.slogdet(M)
    if not np.isfinite(logabsdet):
        return -float('inf')
    log_D = logabsdet - log_norms.sum()
    return log_D


def atiyah_sutcliffe_det(points):
    """Wrapper returning |D| from log|D|."""
    log_D = atiyah_sutcliffe_logD(points)
    if not math.isfinite(log_D):
        return 0.0
    return math.exp(log_D)


def random_S2(rng, n):
    """Sample n points uniformly on S^2."""
    v = rng.standard_normal((n, 3))
    nrms = np.linalg.norm(v, axis=1, keepdims=True)
    return v / nrms


def scan(n, n_configs, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260427 + n)
    log_Ds = []
    for _ in range(n_configs):
        pts = random_S2(rng, n)
        log_D = atiyah_sutcliffe_logD(pts)
        log_Ds.append(log_D)
    return np.array(log_Ds)


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    CASES = [
        (5,  5000),
        (10, 5000),
        (20, 2000),
        (30, 1000),
        (40,  500),
        (50,  200),
    ]
    t0 = time.perf_counter()
    print("=" * 105, flush=True)
    print("Atiyah-Sutcliffe log|D| scan on uniform S^2 configurations  (Atiyah's normalisation)",
          flush=True)
    print("Conjecture: log|D| >= 0 always (i.e., |D| >= 1).", flush=True)
    print("=" * 105, flush=True)
    print(f"{'n':>4} {'#configs':>10} {'min log|D|':>11} {'median':>11} "
          f"{'max':>11} {'#below 0':>9} {'time':>10}", flush=True)
    print("-" * 105, flush=True)
    summary = {}
    for n, ncfg in CASES:
        t1 = time.perf_counter()
        log_Ds = scan(n, ncfg)
        elapsed = time.perf_counter() - t1
        valid = log_Ds[np.isfinite(log_Ds)]
        if len(valid) == 0:
            print(f"{n:>4} {ncfg:>10}  no valid", flush=True)
            continue
        min_lD = float(valid.min())
        med_lD = float(np.median(valid))
        max_lD = float(valid.max())
        n_below = int((valid < 0.0).sum())
        flag = "  <-- VIOLATION!" if n_below > 0 else ""
        print(f"{n:>4} {ncfg:>10} {min_lD:>11.4f} {med_lD:>11.4f} "
              f"{max_lD:>11.4f} {n_below:>9} {elapsed:>9.1f}s{flag}",
              flush=True)
        summary[n] = (min_lD, med_lD, max_lD, n_below, ncfg)

    # Adversarial test 1: near-collinear (perturbations of points on a line in R^3)
    print("\n" + "=" * 95, flush=True)
    print("Adversarial 1: TRULY collinear (points on z-axis at heights 1, 2, ..., n)",
          flush=True)
    print("=" * 95, flush=True)
    for n in [5, 10, 20, 50]:
        # Points on the line {(0, 0, t) : t in R}, well-separated
        pts = np.column_stack([np.zeros(n), np.zeros(n), np.arange(1, n+1, dtype=float)])
        try:
            log_D = atiyah_sutcliffe_logD(pts)
        except Exception as e:
            log_D = float('nan')
        print(f"  n={n:>3}, collinear: log|D| = {log_D:+.6f}  (conjectured = 0)",
              flush=True)

    # Adversarial test 2: near-collinear (perturbations of a line)
    print("\n" + "=" * 95, flush=True)
    print("Adversarial 2: line-perturbed (points on z-axis + delta noise)",
          flush=True)
    print("=" * 95, flush=True)
    rng = np.random.default_rng(424242)
    for n in [5, 10, 20, 50]:
        for delta in [1e-1, 1e-2, 1e-3, 1e-4]:
            pts = np.column_stack([np.zeros(n), np.zeros(n), np.arange(1, n+1, dtype=float)])
            pts += delta * rng.standard_normal((n, 3))
            try:
                log_D = atiyah_sutcliffe_logD(pts)
            except Exception as e:
                log_D = float('nan')
            print(f"  n={n:>3}, delta={delta:>5.0e}: log|D| = {log_D:+.4f}", flush=True)

    # Adversarial test 3: equispaced great circle (planar but not collinear)
    print("\n" + "=" * 95, flush=True)
    print("Adversarial 3: equispaced great circle (coplanar)", flush=True)
    print("=" * 95, flush=True)
    for n in [5, 10, 20, 50]:
        angles = np.linspace(0, 2*np.pi, n, endpoint=False)
        pts = np.column_stack([np.cos(angles), np.sin(angles), np.zeros(n)])
        try:
            log_D = atiyah_sutcliffe_logD(pts)
        except Exception as e:
            log_D = float('nan')
        print(f"  n={n:>3}: log|D| = {log_D:+.4f}", flush=True)

    # LaTeX-ready table
    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table for paper 103 §4:", flush=True)
    print("=" * 70, flush=True)
    for n_cfg in CASES:
        n = n_cfg[0]
        if n in summary:
            min_lD, med_lD, max_lD, n_below, ncfg = summary[n]
            status = ("\\textbf{COUNTEREXAMPLE}" if n_below > 0 else "$\\log|D| \\ge 0$ everywhere")
            print(f"  ${n}$  & ${ncfg:,}$ & ${min_lD:+.3f}$ & ${med_lD:+.3f}$ & "
                  f"${max_lD:+.3f}$ & {status} \\\\",
                  flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
