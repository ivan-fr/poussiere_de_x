"""
Vectorized multivariate KS with Newton-ELS, scaled to large Bezout degrees.
Target: D up to 10^6.

Optimizations:
  - Multi-indices stored as (n_mons, n) int array
  - Coefficients as (n, n_mons) complex array
  - Vectorized eval via z**alpha_arr broadcasting + np.prod
  - Vectorized Jacobian
"""
from __future__ import annotations
import math, time
from itertools import product
import numpy as np

GOLDEN = math.pi * (3.0 - math.sqrt(5.0))
ALPHA_STAR_MV = 0.1577


def multi_indices_array(n, d):
    """Return (n_mons, n) array of all alpha with sum(alpha) <= d."""
    out = []
    for alpha in product(range(d + 1), repeat=n):
        if sum(alpha) <= d:
            out.append(alpha)
    return np.array(out, dtype=int)


def multinomial_log_arr(d, alpha_arr):
    """Vectorized log multinomial(d; alpha_1, ..., alpha_n, d-|alpha|)."""
    n_mons = alpha_arr.shape[0]
    s = np.full(n_mons, math.lgamma(d + 1))
    for j in range(alpha_arr.shape[1]):
        s -= np.array([math.lgamma(int(a) + 1) for a in alpha_arr[:, j]])
    sums = alpha_arr.sum(axis=1)
    s -= np.array([math.lgamma(int(d - s_i) + 1) for s_i in sums])
    return s


def sample_KS_mv_fast(n, d, alpha_arr, sigmas, rng):
    """Sample F = (F_1, ..., F_n) iid KS multivariate. Returns coefs_arr (n, n_mons)."""
    n_mons = alpha_arr.shape[0]
    re = rng.standard_normal((n, n_mons)) * sigmas
    im = rng.standard_normal((n, n_mons)) * sigmas
    return (re + 1j * im) / math.sqrt(2)


def eval_F_fast(coefs_arr, alpha_arr, z):
    """Fast eval F(z) using vectorized numpy."""
    z = np.asarray(z, dtype=complex)
    # alpha_arr (n_mons, n), z (n,) -> z**alpha_arr broadcasts to (n_mons, n)
    # Avoid 0**0 = 1 by careful handling: z[j]**0 = 1 always (numpy convention)
    with np.errstate(invalid='ignore', divide='ignore'):
        powers = z ** alpha_arr  # (n_mons, n)
    monomials = np.prod(powers, axis=1)  # (n_mons,)
    return coefs_arr @ monomials  # (n,)


def jacobian_F_fast(coefs_arr, alpha_arr, z):
    """Fast Jacobian using vectorized numpy."""
    z = np.asarray(z, dtype=complex)
    n = len(z)
    n_mons = alpha_arr.shape[0]
    with np.errstate(invalid='ignore', divide='ignore'):
        base_powers = z ** alpha_arr  # (n_mons, n)
    J = np.zeros((n, n), dtype=complex)
    for k in range(n):
        # d/dz_k: replace col k by alpha[:, k] * z[k]^(alpha[:, k]-1)
        ak = alpha_arr[:, k]
        mask = ak > 0
        zpow_minus1 = np.zeros(n_mons, dtype=complex)
        if abs(z[k]) > 1e-300:
            zpow_minus1[mask] = z[k] ** (ak[mask] - 1)
        else:
            mask1 = ak == 1
            zpow_minus1[mask1] = 1.0
        deriv_col = ak * zpow_minus1  # (n_mons,)
        derived = base_powers.copy()
        derived[:, k] = deriv_col
        derived_monomials = np.prod(derived, axis=1)
        J[:, k] = coefs_arr @ derived_monomials
    return J


def newton_step_fast(coefs_arr, alpha_arr, z):
    Fz = eval_F_fast(coefs_arr, alpha_arr, z)
    if not np.all(np.isfinite(Fz)):
        return None, None
    J = jacobian_F_fast(coefs_arr, alpha_arr, z)
    if not np.all(np.isfinite(J)):
        return None, None
    try:
        Delta = np.linalg.solve(J, Fz)
    except np.linalg.LinAlgError:
        return None, None
    if not np.all(np.isfinite(Delta)):
        return None, None
    return Fz, Delta


def els_step_fast(coefs_arr, alpha_arr, z, T_set):
    Fz, Delta = newton_step_fast(coefs_arr, alpha_arr, z)
    if Delta is None:
        return z, 0.0
    nF0 = float(np.linalg.norm(Fz))
    if nF0 < 1e-300:
        return z, 0.0
    log_nF0 = math.log(nF0)
    best_z = z
    best_log = log_nF0
    for tau in T_set:
        z_try = z - tau * Delta
        try:
            F_try = eval_F_fast(coefs_arr, alpha_arr, z_try)
        except (OverflowError, ZeroDivisionError, FloatingPointError):
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


def starts_Bprime_mv_fast(N, n, R=2.0, q=0.7):
    out = []
    for k in range(N):
        rng = np.random.default_rng(20260427 + k)
        v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        v = v / max(np.linalg.norm(v), 1e-30)
        out.append(R * (q ** k) * v)
    return out


def alpha_mv_proxy_fast(coefs_arr, alpha_arr, z):
    Fz, Delta = newton_step_fast(coefs_arr, alpha_arr, z)
    if Delta is None:
        return float('inf')
    return float(np.linalg.norm(Delta))


def run_orbit_ELS_fast(coefs_arr, alpha_arr, z_init, max_epochs=200,
                       tol_rel=1e-8, T_set=None):
    if T_set is None:
        T_set = [2.0**k for k in range(-6, 7)]
    z = np.asarray(z_init, dtype=complex)
    Fz = eval_F_fast(coefs_arr, alpha_arr, z)
    if not np.all(np.isfinite(Fz)):
        return 'overflow', 0, 0
    Fs_scale = float(np.max(np.abs(coefs_arr)))
    if Fs_scale < 1e-30:
        Fs_scale = 1.0
    basin_tol = tol_rel * Fs_scale
    if float(np.linalg.norm(Fz)) < basin_tol:
        return 'converged', 0, 0
    beta_init = alpha_mv_proxy_fast(coefs_arr, alpha_arr, z)
    n_pre = 0 if (np.isfinite(beta_init) and beta_init < ALPHA_STAR_MV) else -1
    for ep in range(1, max_epochs + 1):
        z_new, _ = els_step_fast(coefs_arr, alpha_arr, z, T_set)
        if z_new is z:
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        z = z_new
        try:
            Fz = eval_F_fast(coefs_arr, alpha_arr, z)
        except (OverflowError, ZeroDivisionError, FloatingPointError):
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        if not np.all(np.isfinite(Fz)):
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        if n_pre < 0:
            beta = alpha_mv_proxy_fast(coefs_arr, alpha_arr, z)
            if np.isfinite(beta) and beta < ALPHA_STAR_MV:
                n_pre = ep
        if float(np.linalg.norm(Fz)) < basin_tol:
            return 'converged', ep, n_pre if n_pre >= 0 else ep
    return 'stagnated', max_epochs, n_pre if n_pre >= 0 else max_epochs


def measure(n, d, n_polys, max_starts=24, max_epochs=200, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260427 + 100*n + d)
    alpha_arr = multi_indices_array(n, d)
    sigmas = np.exp(0.5 * multinomial_log_arr(d, alpha_arr))
    n_pres, n_tots, statuses = [], [], []
    for _ in range(n_polys):
        coefs_arr = sample_KS_mv_fast(n, d, alpha_arr, sigmas, rng)
        starts = starts_Bprime_mv_fast(max_starts, n)
        success = False
        for z0 in starts:
            try:
                st, n_tot, n_pre = run_orbit_ELS_fast(
                    coefs_arr, alpha_arr, z0, max_epochs=max_epochs)
            except (OverflowError, ZeroDivisionError, FloatingPointError):
                continue
            if st == 'converged':
                n_pres.append(n_pre); n_tots.append(n_tot); statuses.append('conv')
                success = True; break
        if not success:
            n_pres.append(max_epochs); n_tots.append(max_epochs); statuses.append('fail')
    return np.array(n_pres), np.array(n_tots), statuses


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    # Cases ranging from D = 9 up to D ~ 10^6
    CASES = [
        (2, 3),    # D = 9
        (3, 3),    # D = 27
        (4, 2),    # D = 16
        (4, 3),    # D = 81
        (5, 3),    # D = 243
        (5, 4),    # D = 1024
        (6, 3),    # D = 729
        (6, 4),    # D = 4096
        (4, 8),    # D = 4096
        (4, 16),   # D = 65536
        (5, 6),    # D = 7776
        (6, 6),    # D = 46656
        (4, 32),   # D = 1048576 ~ 10^6
        (5, 16),   # D = 1048576 ~ 10^6
        (6, 10),   # D = 1000000 ~ 10^6
    ]
    N_POLYS = 8
    t0 = time.perf_counter()
    print("=" * 130, flush=True)
    print("LARGE-SCALE multivariate KS  +  B'-Newton-ELS"
          "  (scaled to D ~ 10^6)", flush=True)
    print("=" * 130, flush=True)
    print(f"{'(n,d)':>8} {'D':>10} {'#mons':>6} | {'%conv':>5} | "
          f"{'<N_pre>':>9} {'med':>5} {'max':>5} | "
          f"{'<N_tot>':>9} {'med':>5} {'max':>5} | "
          f"{'log2 D':>7} {'time/poly':>10}", flush=True)
    print("-" * 130, flush=True)
    summary = {}
    for n, d in CASES:
        D = d ** n
        # Number of monomials
        # = C(d+n, n)
        nmons = math.comb(d + n, n)
        if nmons > 50000:
            print(f"{(n,d)!s:>8} {D:>10} {nmons:>6} -- skipped (too many monomials)",
                  flush=True)
            continue
        t1 = time.perf_counter()
        n_pre, n_tot, statuses = measure(n, d, N_POLYS)
        elapsed = time.perf_counter() - t1
        n_conv = sum(1 for s in statuses if s == 'conv')
        if n_conv == 0:
            print(f"{(n,d)!s:>8} {D:>10} {nmons:>6} |  0% no conv "
                  f"({elapsed/N_POLYS:.1f}s/poly)", flush=True)
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
        print(f"{(n,d)!s:>8} {D:>10} {nmons:>6} | {100*n_conv/N_POLYS:>4.0f}% | "
              f"{m_pre:>9.2f} {med_pre:>5} {mx_pre:>5} | "
              f"{m_tot:>9.2f} {med_tot:>5} {mx_tot:>5} | {ld:>7.2f} "
              f"{elapsed/N_POLYS:>9.2f}s",
              flush=True)
        summary[(n, d)] = (m_pre, m_tot, D)

    # Regression
    if len(summary) >= 3:
        Ds = np.array([summary[k][2] for k in summary])
        log_Ds = np.log2(Ds)
        pres = np.array([summary[k][0] for k in summary])
        A = np.stack([log_Ds, np.ones_like(log_Ds)], axis=1)
        (b, a), *_ = np.linalg.lstsq(A, pres, rcond=None)
        pres_pos = np.maximum(pres, 1e-3)
        (b2, a2), *_ = np.linalg.lstsq(A, np.log2(pres_pos), rcond=None)
        print(f"\n  Linear:  N_pre = {a:.2f} + {b:.2f} log2(D)", flush=True)
        print(f"  Log-log: N_pre ~ D^{b2:.3f} * 2^{a2:.2f}", flush=True)
        if abs(b2) < 0.20:
            print(f"  -> N_pre = O(1), holds at D ~ 10^6", flush=True)
        elif b2 < 0.5:
            print(f"  -> Sub-sqrt(D), still good", flush=True)
        else:
            print(f"  -> Polynomial in D", flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
