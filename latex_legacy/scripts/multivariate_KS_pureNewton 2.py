"""Test pure Newton (tau=1, no ELS) on multivariate KS to confirm:
the ELS mechanism is not the load-bearing fix multivariately."""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from multivariate_KS import (
    sample_KS_mv, eval_F, jacobian_F, starts_Bprime_mv, newton_step,
    alpha_mv_proxy, ALPHA_STAR_MV,
)


def run_orbit_pureNewton(Fs, z_init, max_epochs=200, tol_rel=1e-10):
    n = len(Fs)
    z = np.asarray(z_init, dtype=complex)
    Fz = eval_F(Fs, z)
    if not np.all(np.isfinite(Fz)):
        return 'overflow', 0, 0
    Fs_scale = max(1.0, sum(sum(abs(c) for c in fi.values()) for fi in Fs))
    basin_tol = tol_rel * Fs_scale
    if float(np.linalg.norm(Fz)) < basin_tol:
        return 'converged', 0, 0
    beta_init = alpha_mv_proxy(Fs, z)
    n_pre = 0 if (np.isfinite(beta_init) and beta_init < ALPHA_STAR_MV) else -1
    for ep in range(1, max_epochs + 1):
        Fz, Delta = newton_step(Fs, z)
        if Delta is None:
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        z = z - Delta  # Pure Newton (tau = 1)
        try:
            Fz_new = eval_F(Fs, z)
        except (OverflowError, ZeroDivisionError):
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        if not np.all(np.isfinite(Fz_new)):
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        if n_pre < 0:
            beta = alpha_mv_proxy(Fs, z)
            if np.isfinite(beta) and beta < ALPHA_STAR_MV:
                n_pre = ep
        if float(np.linalg.norm(Fz_new)) < basin_tol:
            return 'converged', ep, n_pre if n_pre >= 0 else ep
    return 'stagnated', max_epochs, n_pre if n_pre >= 0 else max_epochs


def measure(n, d, n_polys, max_starts=24, max_epochs=200, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260427 + 100*n + d)
    n_tots, n_pres, statuses = [], [], []
    for _ in range(n_polys):
        Fs = sample_KS_mv(n, d, rng)
        starts = starts_Bprime_mv(max_starts, n)
        success = False
        for z0 in starts:
            try:
                st, n_tot, n_pre = run_orbit_pureNewton(Fs, z0, max_epochs=max_epochs)
            except (OverflowError, ZeroDivisionError):
                continue
            if st == 'converged':
                n_tots.append(n_tot); n_pres.append(n_pre); statuses.append('conv')
                success = True; break
        if not success:
            n_tots.append(max_epochs); n_pres.append(max_epochs); statuses.append('fail')
    return np.array(n_tots), np.array(n_pres), statuses


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    CASES = [(2, 3), (3, 3), (4, 2), (4, 3), (5, 2), (5, 3), (6, 2), (5, 4), (6, 3)]
    N_POLYS = 12
    print("=" * 110, flush=True)
    print("Pure Newton (tau=1, no ELS) on mvKS — control for Part IX-mv §5.4", flush=True)
    print("=" * 110, flush=True)
    print(f"{'(n,d)':>8} {'D':>5} | {'%conv':>5} | "
          f"{'<N_pre>':>9} {'med':>5} {'max':>5} | {'<N_tot>':>9} {'med':>5} {'max':>5}",
          flush=True)
    print("-" * 110, flush=True)
    for n, d in CASES:
        D = d ** n
        n_tot, n_pre, statuses = measure(n, d, N_POLYS)
        n_conv = sum(1 for s in statuses if s == 'conv')
        if n_conv == 0:
            print(f"{(n,d)!s:>8} {D:>5} |  0% no conv", flush=True)
            continue
        mask = np.array([s == 'conv' for s in statuses])
        np_c = n_pre[mask]; nt_c = n_tot[mask]
        m_pre = float(np.mean(np_c))
        m_tot = float(np.mean(nt_c))
        print(f"{(n,d)!s:>8} {D:>5} | {100*n_conv/N_POLYS:>4.0f}% | "
              f"{m_pre:>9.2f} {int(np.median(np_c)):>5} {int(np.max(np_c)):>5} | "
              f"{m_tot:>9.2f} {int(np.median(nt_c)):>5} {int(np.max(nt_c)):>5}",
              flush=True)


if __name__ == "__main__":
    main()
