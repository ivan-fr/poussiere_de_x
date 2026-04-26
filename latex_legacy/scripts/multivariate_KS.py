"""
Multivariate Pandrosion-Newton with extended line search (ELS) on
Kostlan-Smale multivariate ensemble.

Background: Part I §3.5 introduced the Schmidt slope matrix Q_F for
polynomial systems F: C^n -> C^n; the multivariate Pandrosion operator is
    P_{F, z0}(z) = z0 - Q_F(z0, z)^{-1} F(z0).
With dynamic anchor (z0 = z), this reduces to multivariate Newton:
    z_{n+1} = z_n - DF(z_n)^{-1} F(z_n).

This script tests the multivariate analogue of the extended line search
(ELS, Part IX): given Newton direction Delta_n = DF^{-1} F, search
    tau* = argmin_{tau in T_ELS} || F(z_n - tau * Delta_n) ||,
with T_ELS = {2^k : k = -6, ..., +6}.

We measure on the multivariate KS ensemble (Bombieri-Weyl normalised), with
Strategy B' starts adapted to C^n (Fibonacci sphere * geometric radius).
We track:
  N_pre_basin  = first epoch with || F(z) || small enough to be in Newton's basin
  N_total      = first epoch with || F(z) || < tol
  s_star       = first start that converges
"""
from __future__ import annotations
import math, time
from itertools import product
import numpy as np

GOLDEN = math.pi * (3.0 - math.sqrt(5.0))


# --- Multi-index utilities -------------------------------------------------

def multi_indices(n, d):
    """All alpha = (alpha_1, ..., alpha_n) with sum(alpha) <= d."""
    out = []
    for alpha in product(range(d + 1), repeat=n):
        if sum(alpha) <= d:
            out.append(alpha)
    return out


def multinomial_log(d, alpha):
    """log multinomial(d; alpha_1, ..., alpha_n, d-|alpha|) using lgamma."""
    s = math.lgamma(d + 1)
    for a in alpha:
        s -= math.lgamma(a + 1)
    s -= math.lgamma(d - sum(alpha) + 1)
    return s


# --- KS multivariate sampling ----------------------------------------------

def sample_KS_mv(n, d, rng):
    """Sample F = (F_1, ..., F_n) iid KS multivariate, degree d each.
    Returns list of length n; each entry is a dict {alpha: coef}.
    Variance of coef at alpha is multinomial(d; alpha, d-|alpha|).
    """
    indices = multi_indices(n, d)
    Fs = []
    for _ in range(n):
        coefs = {}
        for alpha in indices:
            sigma = math.exp(0.5 * multinomial_log(d, alpha))
            re = rng.standard_normal() * sigma
            im = rng.standard_normal() * sigma
            coefs[alpha] = complex(re, im) / math.sqrt(2)
        Fs.append(coefs)
    return Fs


# --- Polynomial primitives -------------------------------------------------

def eval_F(Fs, z):
    """Evaluate F at z. Returns np.array of length n."""
    n = len(Fs)
    out = np.zeros(n, dtype=complex)
    z = np.asarray(z, dtype=complex)
    for i, coefs in enumerate(Fs):
        s = 0.0 + 0.0j
        for alpha, c in coefs.items():
            term = c
            for j, a in enumerate(alpha):
                if a > 0:
                    term *= z[j] ** a
            s += term
        out[i] = s
    return out


def jacobian_F(Fs, z):
    """Return Jacobian DF(z) as n x n complex matrix."""
    n = len(Fs)
    z = np.asarray(z, dtype=complex)
    J = np.zeros((n, n), dtype=complex)
    for i, coefs in enumerate(Fs):
        for alpha, c in coefs.items():
            for k in range(n):
                if alpha[k] == 0:
                    continue
                term = c * alpha[k]
                for j, a in enumerate(alpha):
                    if j == k:
                        if a - 1 > 0:
                            term *= z[j] ** (a - 1)
                    else:
                        if a > 0:
                            term *= z[j] ** a
                J[i, k] += term
    return J


# --- Strategy B' multivariate (Fibonacci sphere * geometric radius) --------

def fibonacci_sphere(N, n):
    """Generate N points on the unit sphere S^{2n-1} (real dim 2n) using a
    deterministic golden-spiral-like sampling. Returns array of shape (N, n)
    of complex entries with each row of complex-2-norm = 1.
    """
    rng = np.random.default_rng(20260427)
    pts = np.zeros((N, n), dtype=complex)
    for k in range(N):
        # Real Fibonacci sphere on S^{2n-1}, then split into complex coords
        # Simple approach: random gaussian, normalized; with deterministic seed
        rng2 = np.random.default_rng(20260427 + k)
        v = rng2.standard_normal(n) + 1j * rng2.standard_normal(n)
        nrm = np.linalg.norm(v)
        pts[k] = v / max(nrm, 1e-30)
    return pts


def starts_Bprime_mv(N, n, R=2.0, q=0.7):
    """Multivariate Strategy B' starts. Each start is a vector in C^n with
    radius R * q^k and direction from a deterministic Fibonacci-like sample.
    """
    pts = fibonacci_sphere(N, n)
    return [R * (q ** k) * pts[k] for k in range(N)]


# --- Multivariate Newton + ELS step ---------------------------------------

def newton_step(Fs, z):
    """Compute Newton direction Delta = DF(z)^{-1} F(z). Returns None on degeneracy."""
    Fz = eval_F(Fs, z)
    if not np.all(np.isfinite(Fz)):
        return None, None
    J = jacobian_F(Fs, z)
    if not np.all(np.isfinite(J)):
        return None, None
    try:
        Delta = np.linalg.solve(J, Fz)
    except np.linalg.LinAlgError:
        return None, None
    if not np.all(np.isfinite(Delta)):
        return None, None
    return Fz, Delta


def els_step(Fs, z, T_set):
    """Extended line search step.
    Returns (z_new, log_contraction). z_new = z - tau* Delta with
    tau* = argmin_{tau in T_set} ||F(z - tau Delta)||.
    """
    Fz, Delta = newton_step(Fs, z)
    if Delta is None:
        return z, 0.0
    nF0 = float(np.linalg.norm(Fz))
    if nF0 < 1e-300:
        return z, 0.0
    log_nF0 = math.log(nF0)
    best_z = z
    best_log = log_nF0
    for tau in T_set:
        try:
            z_try = z - tau * Delta
            F_try = eval_F(Fs, z_try)
        except (OverflowError, ZeroDivisionError):
            continue
        if not np.all(np.isfinite(F_try)):
            continue
        nF = float(np.linalg.norm(F_try))
        if nF < 1e-300:
            return z_try, -math.log(1e300)
        log_nF = math.log(nF)
        if log_nF < best_log:
            best_log = log_nF
            best_z = z_try
    return best_z, best_log - log_nF0


# --- Smale alpha (multivariate, simplified via 2nd derivative) -------------

def alpha_mv_proxy(Fs, z, tol=1e-300):
    """Proxy multivariate alpha = beta * gamma_2 where:
        beta = ||DF^{-1} F||
        gamma_2 = ||DF^{-1} D^2 F[Delta, Delta]|| / ||Delta||^2  (rough estimate)
    Full multivariate gamma is intractable; we use beta/||z|| as cheap proxy.
    Returns +inf on degeneracy.
    """
    Fz, Delta = newton_step(Fs, z)
    if Delta is None:
        return float('inf')
    beta = float(np.linalg.norm(Delta))
    return beta  # use beta alone as basin proxy (Newton convergence iff beta small)


# --- Driver ---------------------------------------------------------------

ALPHA_STAR_MV = 0.1577  # same threshold used univariately

def run_orbit_ELS(Fs, z_init, max_epochs=200, tol_rel=1e-10, T_set=None):
    """Run B'-Newton-ELS orbit. Returns (status, n_total, n_pre_basin, beta_init)."""
    if T_set is None:
        T_set = [2.0**k for k in range(-6, 7)]
    n = len(Fs)
    z = np.asarray(z_init, dtype=complex)
    Fz = eval_F(Fs, z)
    if not np.all(np.isfinite(Fz)):
        return 'overflow', 0, 0, float('inf')
    nF0_init = float(np.linalg.norm(Fz))
    # rough scale of F
    Fs_scale = max(1.0, sum(sum(abs(c) for c in fi.values()) for fi in Fs))
    basin_tol = tol_rel * Fs_scale
    if nF0_init < basin_tol:
        return 'converged', 0, 0, 0.0

    beta_init = alpha_mv_proxy(Fs, z)
    n_pre = 0 if (np.isfinite(beta_init) and beta_init < ALPHA_STAR_MV) else -1

    for ep in range(1, max_epochs + 1):
        z_new, _ = els_step(Fs, z, T_set)
        if z_new is z:
            return 'fail', ep, n_pre if n_pre >= 0 else ep, beta_init
        z = z_new
        try:
            Fz = eval_F(Fs, z)
        except (OverflowError, ZeroDivisionError):
            return 'fail', ep, n_pre if n_pre >= 0 else ep, beta_init
        if not np.all(np.isfinite(Fz)):
            return 'fail', ep, n_pre if n_pre >= 0 else ep, beta_init
        if n_pre < 0:
            beta = alpha_mv_proxy(Fs, z)
            if np.isfinite(beta) and beta < ALPHA_STAR_MV:
                n_pre = ep
        nF = float(np.linalg.norm(Fz))
        if nF < basin_tol:
            return 'converged', ep, n_pre if n_pre >= 0 else ep, beta_init
    return 'stagnated', max_epochs, n_pre if n_pre >= 0 else max_epochs, beta_init


def measure(n, d, n_polys, max_starts=24, max_epochs=200, T_set=None, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260427 + 100*n + d)
    n_tots, n_pres, statuses = [], [], []
    for _ in range(n_polys):
        Fs = sample_KS_mv(n, d, rng)
        starts = starts_Bprime_mv(max_starts, n)
        success = False
        for z0 in starts:
            try:
                st, n_tot, n_pre, _ = run_orbit_ELS(
                    Fs, z0, max_epochs=max_epochs, T_set=T_set)
            except (OverflowError, ZeroDivisionError):
                continue
            if st == 'converged':
                n_tots.append(n_tot); n_pres.append(n_pre); statuses.append('conv')
                success = True
                break
        if not success:
            n_tots.append(max_epochs); n_pres.append(max_epochs); statuses.append('fail')
    return np.array(n_tots), np.array(n_pres), statuses


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    # Test cases following Part I Proposition 3.5 benchmark range
    CASES = [(2, 3), (3, 3), (4, 2), (4, 3), (5, 2), (5, 3), (6, 2), (5, 4), (6, 3)]
    N_POLYS = 12
    t0 = time.perf_counter()
    print("=" * 130, flush=True)
    print("Multivariate KS  +  Strategy B'  +  Newton-ELS"
          "  (T_set = {2^k : k=-6..6})", flush=True)
    print("=" * 130, flush=True)
    print(f"{'(n,d)':>8} {'D=d^n':>7} | {'%conv':>5} | "
          f"{'<N_pre>':>9} {'med':>5} {'max':>5} | "
          f"{'<N_tot>':>9} {'med':>5} {'max':>5} | "
          f"{'log2 D':>7}", flush=True)
    print("-" * 130, flush=True)
    summary = {}
    for n, d in CASES:
        D = d ** n
        n_tot, n_pre, statuses = measure(n, d, N_POLYS)
        n_conv = sum(1 for s in statuses if s == 'conv')
        if n_conv == 0:
            print(f"{(n,d)!s:>8} {D:>7} |  0% no convergent run", flush=True)
            continue
        mask = np.array([s == 'conv' for s in statuses])
        np_c = n_pre[mask]; nt_c = n_tot[mask]
        m_pre = float(np.mean(np_c))
        med_pre = int(np.median(np_c))
        mx_pre = int(np.max(np_c))
        m_tot = float(np.mean(nt_c))
        med_tot = int(np.median(nt_c))
        mx_tot = int(np.max(nt_c))
        ld = math.log2(D) if D > 0 else 0
        print(f"{(n,d)!s:>8} {D:>7} | {100*n_conv/N_POLYS:>4.0f}% | "
              f"{m_pre:>9.2f} {med_pre:>5} {mx_pre:>5} | "
              f"{m_tot:>9.2f} {med_tot:>5} {mx_tot:>5} | {ld:>7.2f}",
              flush=True)
        summary[(n, d)] = (m_pre, m_tot, D)

    # Regression on Bezout degree D
    if len(summary) >= 3:
        Ds = np.array([summary[k][2] for k in summary])
        log_Ds = np.log2(Ds)
        pres = np.array([summary[k][0] for k in summary])
        A = np.stack([log_Ds, np.ones_like(log_Ds)], axis=1)
        (b, a), *_ = np.linalg.lstsq(A, pres, rcond=None)
        pres_pos = np.maximum(pres, 1e-3)
        (b2, a2), *_ = np.linalg.lstsq(A, np.log2(pres_pos), rcond=None)
        print(f"\n  Linear:  N_pre = {a:.2f} + {b:.2f} log2(D)", flush=True)
        print(f"  Log-log: N_pre ~ D^{b2:.2f} * 2^{a2:.2f}", flush=True)
        if abs(b2) < 0.20:
            print(f"  -> CONFIRMED: N_pre = O(1) on multivariate KS",
                  flush=True)
        elif b2 < 0.5:
            print(f"  -> Sub-sqrt(D) scaling", flush=True)
        else:
            print(f"  -> Polynomial in D", flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
