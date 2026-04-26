"""Push to D ~ 10^9 with optimized numpy."""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from multivariate_KS_fast import (
    multi_indices_array, multinomial_log_arr, sample_KS_mv_fast,
    starts_Bprime_mv_fast, run_orbit_ELS_fast,
)


def measure_huge(n, d, n_polys, max_starts=24, max_epochs=200, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260427 + 100*n + d)
    alpha_arr = multi_indices_array(n, d)
    sigmas = np.exp(0.5 * multinomial_log_arr(d, alpha_arr))
    n_pres, n_tots, statuses = [], [], []
    t0 = time.perf_counter()
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
    elapsed = time.perf_counter() - t0
    return np.array(n_pres), np.array(n_tots), statuses, elapsed


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    # Configurations spanning D = 10^4 to 10^9
    CASES = [
        (4, 10),    # D = 1e4
        (5, 10),    # D = 1e5
        (6, 10),    # D = 1e6
        (8, 10),    # D = 1e8
        (15, 4),    # D = 1.07e9     (3876 monomials, fast)
        (20, 3),    # D = 3.49e9     (1771 monomials, fast)
        (12, 6),    # D = 2.18e9     (18564 monomials)
        (8, 8),     # D = 1.68e7     (12870 monomials)
        (10, 6),    # D = 6.05e7     (8008 monomials)
        (10, 8),    # D = 1.07e9     (43758 monomials)
    ]
    N_POLYS = 4
    t0 = time.perf_counter()
    print("=" * 140, flush=True)
    print("HUGE-SCALE multivariate KS  +  B'-Newton-ELS  (D up to 10^9-10^10)",
          flush=True)
    print("=" * 140, flush=True)
    print(f"{'(n,d)':>8} {'D':>12} {'#mons':>8} | {'%conv':>5} | "
          f"{'<N_pre>':>9} {'med':>5} {'max':>5} | "
          f"{'<N_tot>':>9} {'med':>5} {'max':>5} | "
          f"{'log2 D':>7} {'time/poly':>11}", flush=True)
    print("-" * 140, flush=True)
    summary = {}
    for n, d in CASES:
        D = d ** n
        nmons = math.comb(d + n, n)
        if nmons > 250000:
            print(f"{(n,d)!s:>8} {D:>12} {nmons:>8} -- skipped (too many monomials)",
                  flush=True)
            continue
        n_pre, n_tot, statuses, elapsed = measure_huge(n, d, N_POLYS)
        n_conv = sum(1 for s in statuses if s == 'conv')
        if n_conv == 0:
            print(f"{(n,d)!s:>8} {D:>12} {nmons:>8} |  0% no conv "
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
        ld = math.log2(D)
        print(f"{(n,d)!s:>8} {D:>12} {nmons:>8} | {100*n_conv/N_POLYS:>4.0f}% | "
              f"{m_pre:>9.2f} {med_pre:>5} {mx_pre:>5} | "
              f"{m_tot:>9.2f} {med_tot:>5} {mx_tot:>5} | {ld:>7.2f} "
              f"{elapsed/N_POLYS:>10.2f}s",
              flush=True)
        summary[(n, d)] = (m_pre, m_tot, D)

    # Regression
    if len(summary) >= 3:
        Ds = np.array([summary[k][2] for k in summary], dtype=float)
        log_Ds = np.log2(Ds)
        pres = np.array([summary[k][0] for k in summary])
        tots = np.array([summary[k][1] for k in summary])
        A = np.stack([log_Ds, np.ones_like(log_Ds)], axis=1)
        for label, ys in [("N_pre", pres), ("N_tot", tots)]:
            (b, a), *_ = np.linalg.lstsq(A, ys, rcond=None)
            ys_pos = np.maximum(ys, 1e-3)
            (b2, a2), *_ = np.linalg.lstsq(A, np.log2(ys_pos), rcond=None)
            print(f"\n  {label}: linear y = {a:.2f} + {b:.3f} log2 D  | "
                  f"log-log y ~ D^{b2:.4f} * 2^{a2:.2f}", flush=True)
        # log log log
        loglog_Ds = np.log2(np.log2(Ds + 1.0))
        A2 = np.stack([loglog_Ds, np.ones_like(loglog_Ds)], axis=1)
        (b3, a3), *_ = np.linalg.lstsq(A2, tots, rcond=None)
        print(f"\n  N_tot vs log log D: y = {a3:.2f} + {b3:.2f} log2 log2 D", flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
