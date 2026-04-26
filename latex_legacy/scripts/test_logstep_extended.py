"""Confirm LOG-STEP scaling on KS at higher d with more polynomials."""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from test_super_linear_prebasin import run_logstep, starts_Bprime
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from verify_armijo_KS import sample_KS


def measure(d, n_polys, max_starts=24, max_epochs=400, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260427 + d)
    n_pres, n_tots, statuses = [], [], []
    for _ in range(n_polys):
        coefs = sample_KS(d, rng)
        starts = starts_Bprime(max_starts)
        success = False
        for z0 in starts:
            try:
                st, n_tot, n_pre = run_logstep(coefs, z0, max_epochs=max_epochs)
            except (OverflowError, ZeroDivisionError):
                continue
            if st == 'converged':
                n_pres.append(n_pre)
                n_tots.append(n_tot)
                statuses.append('conv')
                success = True
                break
        if not success:
            n_pres.append(max_epochs)
            n_tots.append(max_epochs)
            statuses.append('fail')
    return np.array(n_pres), np.array(n_tots), statuses


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    DEGREES = [16, 32, 64, 128, 256]
    N_POLYS = 30
    t0 = time.perf_counter()
    print("=" * 110, flush=True)
    print(f"LOG-STEP extended verification on KS  (N_POLYS = {N_POLYS})", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'%conv':>5} | {'<N_pre>':>9} {'med_pre':>8} {'p90_pre':>8} "
          f"{'max_pre':>8} | {'<N_tot>':>9} {'med_tot':>8} {'max_tot':>8} | log2 d",
          flush=True)
    print("-" * 110, flush=True)
    summary = {}
    for d in DEGREES:
        n_pre, n_tot, statuses = measure(d, N_POLYS)
        nc = sum(1 for s in statuses if s == 'conv')
        if nc == 0:
            print(f"{d:>4} |  0%   no convergent runs", flush=True)
            continue
        mask = np.array([s == 'conv' for s in statuses])
        np_c = n_pre[mask]; nt_c = n_tot[mask]
        m_pre = float(np.mean(np_c))
        med_pre = int(np.median(np_c))
        p90 = int(np.percentile(np_c, 90))
        mx_pre = int(np.max(np_c))
        m_tot = float(np.mean(nt_c))
        med_tot = int(np.median(nt_c))
        mx_tot = int(np.max(nt_c))
        ld = math.log2(d)
        print(f"{d:>4} | {100*nc/N_POLYS:>4.0f}% | "
              f"{m_pre:>9.2f} {med_pre:>8} {p90:>8} {mx_pre:>8} | "
              f"{m_tot:>9.2f} {med_tot:>8} {mx_tot:>8} | {ld:>5.2f}",
              flush=True)
        summary[d] = (m_pre, m_tot)

    if len(summary) >= 3:
        ds = list(summary.keys())
        log_ds = np.array([math.log2(d) for d in ds])
        pres = np.array([summary[d][0] for d in ds])
        A = np.stack([log_ds, np.ones_like(log_ds)], axis=1)
        (b, a), *_ = np.linalg.lstsq(A, pres, rcond=None)
        pres_pos = np.maximum(pres, 1e-3)
        (b2, a2), *_ = np.linalg.lstsq(A, np.log2(pres_pos), rcond=None)
        print(f"\n  Linear:  N_pre = {a:.2f} + {b:.2f} log2(d)", flush=True)
        print(f"  Log-log: N_pre ~ d^{b2:.2f} * 2^{a2:.2f}", flush=True)
        if abs(b2) < 0.20:
            print(f"  -> CONFIRMED: N_pre = O(1) on KS, ABD rescued", flush=True)
        elif b2 < 0.5:
            print(f"  -> Sub-sqrt(d) scaling", flush=True)
        else:
            print(f"  -> Polynomial scaling, ABD still refuted", flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
