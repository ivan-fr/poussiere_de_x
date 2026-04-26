"""Measure observed optimal tau* on multivariate KS, for table in Part IX-mv §5.4."""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from multivariate_KS import (
    sample_KS_mv, eval_F, jacobian_F, starts_Bprime_mv, newton_step,
)


def first_step_tau(Fs, z, T_set):
    Fz, Delta = newton_step(Fs, z)
    if Delta is None:
        return None
    nF0 = float(np.linalg.norm(Fz))
    if nF0 < 1e-300:
        return 1.0
    log_nF0 = math.log(nF0)
    best_tau = None
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
            return tau
        log_nF = math.log(nF)
        if log_nF < best_log:
            best_log = log_nF
            best_tau = tau
    return best_tau


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    CASES = [(2, 3), (3, 3), (4, 2), (4, 3), (5, 2), (5, 3), (6, 2), (5, 4), (6, 3)]
    N_POLYS = 20
    T_SET = [2.0**k for k in range(-6, 7)]
    rng = np.random.default_rng(20260427)
    print("=" * 90, flush=True)
    print(f"Observed optimal tau* at first epoch on mvKS (T_ELS = {{2^k : k=-6..6}})",
          flush=True)
    print("=" * 90, flush=True)
    print(f"{'(n,d)':>8} {'D':>5} {'sqrt(D)':>8} | "
          f"{'mean tau*':>10} {'med tau*':>10} {'mode tau*':>10}",
          flush=True)
    print("-" * 90, flush=True)
    for n, d in CASES:
        D = d ** n
        taus = []
        for _ in range(N_POLYS):
            Fs = sample_KS_mv(n, d, rng)
            starts = starts_Bprime_mv(8, n)
            for z0 in starts:
                tau = first_step_tau(Fs, z0, T_SET)
                if tau is not None:
                    taus.append(tau)
                    break
        if not taus:
            print(f"{(n,d)!s:>8} {D:>5} {math.sqrt(D):>8.1f} | (no data)", flush=True)
            continue
        arr = np.array(taus)
        m = float(np.mean(arr))
        med = float(np.median(arr))
        from collections import Counter
        cnt = Counter(taus)
        mode = max(cnt, key=cnt.get)
        print(f"{(n,d)!s:>8} {D:>5} {math.sqrt(D):>8.1f} | "
              f"{m:>10.2f} {med:>10.2f} {mode:>10.2f}",
              flush=True)


if __name__ == "__main__":
    main()
