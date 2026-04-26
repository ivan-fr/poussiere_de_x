"""Robust large-scale measurement on mvKS, with per-case wall-clock budget."""
from __future__ import annotations
import math, time, sys, os, signal
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from multivariate_KS_fast import (
    multi_indices_array, multinomial_log_arr, sample_KS_mv_fast,
    starts_Bprime_mv_fast, run_orbit_ELS_fast,
)


def measure_with_budget(n, d, n_polys, wall_budget_sec,
                        max_starts=6, max_epochs=80, rng=None,
                        per_orbit_sec=20.0):
    """Like measure_huge, but stops after wall_budget_sec elapsed.
    Each orbit is also bounded by per_orbit_sec (wall-clock timeout via signal.alarm)."""
    if rng is None:
        rng = np.random.default_rng(20260427 + 100*n + d)
    alpha_arr = multi_indices_array(n, d)
    sigmas = np.exp(0.5 * multinomial_log_arr(d, alpha_arr))
    n_pres, n_tots, statuses = [], [], []
    t0 = time.perf_counter()

    class OrbitTimeout(Exception):
        pass

    def _alarm_handler(signum, frame):
        raise OrbitTimeout()

    old_handler = signal.signal(signal.SIGALRM, _alarm_handler)

    try:
        for k in range(n_polys):
            if time.perf_counter() - t0 > wall_budget_sec:
                break
            coefs_arr = sample_KS_mv_fast(n, d, alpha_arr, sigmas, rng)
            starts = starts_Bprime_mv_fast(max_starts, n)
            success = False
            for z0 in starts:
                signal.setitimer(signal.ITIMER_REAL, per_orbit_sec)
                try:
                    st, n_tot, n_pre = run_orbit_ELS_fast(
                        coefs_arr, alpha_arr, z0, max_epochs=max_epochs)
                except OrbitTimeout:
                    st, n_tot, n_pre = 'timeout', max_epochs, max_epochs
                except (OverflowError, ZeroDivisionError, FloatingPointError):
                    st, n_tot, n_pre = 'fail', max_epochs, max_epochs
                finally:
                    signal.setitimer(signal.ITIMER_REAL, 0)  # cancel alarm
                if st == 'converged':
                    n_pres.append(n_pre); n_tots.append(n_tot); statuses.append('conv')
                    success = True; break
            if not success:
                n_pres.append(max_epochs); n_tots.append(max_epochs)
                statuses.append('fail')
    finally:
        signal.signal(signal.SIGALRM, old_handler)
    elapsed = time.perf_counter() - t0
    n_done = len(statuses)
    return np.array(n_pres), np.array(n_tots), statuses, elapsed, n_done


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    # Configurations spanning D from 10^4 to ~10^9
    # Format: (n, d, wall_budget_sec)
    CASES = [
        (4, 10,   30),     # D = 1e4
        (5, 10,   30),     # D = 1e5
        (6, 10,   60),     # D = 1e6
        (8, 10,  120),     # D = 1e8
        (15, 4,   60),     # D = 1.07e9
        (20, 3,   60),     # D = 3.49e9
        (12, 6,   90),     # D = 2.18e9
        (10, 8,  120),     # D = 1.07e9
        (16, 4,  120),     # D = 4.29e9
    ]
    N_POLYS = 3
    t0 = time.perf_counter()
    print("=" * 140, flush=True)
    print("HUGE-SCALE mvKS  +  B'-Newton-ELS  v2 (per-case wall budget)",
          flush=True)
    print("=" * 140, flush=True)
    print(f"{'(n,d)':>8} {'D':>13} {'#mons':>8} | {'budget':>7} | {'done':>5} | "
          f"{'<N_pre>':>9} {'max':>5} | "
          f"{'<N_tot>':>9} {'max':>5} | "
          f"{'log2 D':>7} {'time/poly':>11}", flush=True)
    print("-" * 140, flush=True)
    summary = {}
    for n, d, budget in CASES:
        D = d ** n
        nmons = math.comb(d + n, n)
        if nmons > 250000:
            print(f"{(n,d)!s:>8} {D:>13} {nmons:>8} -- skipped (too many monomials)",
                  flush=True)
            continue
        n_pre, n_tot, statuses, elapsed, n_done = measure_with_budget(
            n, d, N_POLYS, wall_budget_sec=budget)
        if n_done == 0:
            print(f"{(n,d)!s:>8} {D:>13} {nmons:>8} | {budget:>5}s | "
                  f"{0:>5}/{N_POLYS} | timeout before any poly", flush=True)
            continue
        n_conv = sum(1 for s in statuses if s == 'conv')
        if n_conv == 0:
            print(f"{(n,d)!s:>8} {D:>13} {nmons:>8} | {budget:>5}s | "
                  f"{n_done:>5}/{N_POLYS} |  no conv ({elapsed/n_done:.1f}s/poly)",
                  flush=True)
            continue
        mask = np.array([s == 'conv' for s in statuses])
        np_c = n_pre[mask]; nt_c = n_tot[mask]
        m_pre = float(np.mean(np_c))
        mx_pre = int(np.max(np_c))
        m_tot = float(np.mean(nt_c))
        mx_tot = int(np.max(nt_c))
        ld = math.log2(D)
        print(f"{(n,d)!s:>8} {D:>13} {nmons:>8} | {budget:>5}s | "
              f"{n_done:>2}/{N_POLYS}({n_conv}c) | "
              f"{m_pre:>9.2f} {mx_pre:>5} | "
              f"{m_tot:>9.2f} {mx_tot:>5} | {ld:>7.2f} "
              f"{elapsed/n_done:>10.2f}s",
              flush=True)
        summary[(n, d)] = (m_pre, m_tot, D)

    # Regression
    if len(summary) >= 3:
        Ds = np.array([summary[k][2] for k in summary], dtype=float)
        log_Ds = np.log2(Ds)
        loglog_Ds = np.log2(np.log2(Ds + 1.0))
        pres = np.array([summary[k][0] for k in summary])
        tots = np.array([summary[k][1] for k in summary])
        A = np.stack([log_Ds, np.ones_like(log_Ds)], axis=1)
        A2 = np.stack([loglog_Ds, np.ones_like(loglog_Ds)], axis=1)
        for label, ys in [("N_pre", pres), ("N_tot", tots)]:
            (b, a), *_ = np.linalg.lstsq(A, ys, rcond=None)
            ys_pos = np.maximum(ys, 1e-3)
            (b2, a2), *_ = np.linalg.lstsq(A, np.log2(ys_pos), rcond=None)
            (b3, a3), *_ = np.linalg.lstsq(A2, ys, rcond=None)
            print(f"\n  {label}:", flush=True)
            print(f"    linear:  y = {a:.2f} + {b:.3f} log2 D",       flush=True)
            print(f"    log-log: y ~ D^{b2:.4f} * 2^{a2:.2f}",         flush=True)
            print(f"    loglog:  y = {a3:.2f} + {b3:.2f} log2 log2 D", flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
