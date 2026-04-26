"""
Same epoch-count measurement as verify_epochs_KS.py, but using Strategy A
(equispaced Cauchy ring) instead of Strategy B'.

Strategy A from Part VIII Sec 4 (Strategy A baseline):
    z_k = R * exp(2 pi i k / d),     k = 0, 1, ..., d-1.

We use R = 2 (KS-adaptive radius) so the comparison with verify_epochs_KS.py
isolates the start-pattern effect and avoids the Cauchy-bound overflow on KS
documented in paper 101bis.
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import (
    eval_poly, T2_step, armijo_fallback,
)
from verify_armijo_KS import sample_KS
from verify_epochs_KS import run_T2_armijo_count_epochs


def starts_A(N, R=2.0):
    return [R * np.exp(2j * np.pi * k / N) for k in range(N)]


def measure(d, n_polys, max_starts, max_epochs=400, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260426 + d)
    epoch_counts = []
    for _ in range(n_polys):
        coefs = sample_KS(d, rng)
        starts = starts_A(max_starts)
        for k, z0 in enumerate(starts):
            try:
                st, ep = run_T2_armijo_count_epochs(coefs, z0, max_epochs=max_epochs)
            except (OverflowError, ZeroDivisionError):
                continue
            if st == 'converged':
                epoch_counts.append(ep)
                break
        else:
            epoch_counts.append(max_epochs)
    return np.array(epoch_counts, dtype=float)


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    DEGREES = [16, 32, 64, 128, 256, 512]
    N_POLYS = 60
    t0 = time.perf_counter()
    print("=" * 110, flush=True)
    print("Strategy A (equispaced ring, R=2) — KS epoch-count measurement", flush=True)
    print(f"  N_POLYS = {N_POLYS}", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'mean ep':>9} {'med':>5} {'p90':>5} {'p99':>5} {'max':>5} | "
          f"{'log2(d)':>8} {'mean/log2':>10}", flush=True)
    print("-" * 110, flush=True)
    summary = {}
    for d in DEGREES:
        # Strategy A uses d starts (one per root) per Part VIII
        eps = measure(d, N_POLYS, max_starts=d)
        mean = float(np.mean(eps))
        med = int(np.median(eps))
        p90 = int(np.percentile(eps, 90))
        p99 = int(np.percentile(eps, 99))
        mx = int(np.max(eps))
        ld = math.log2(d)
        ratio = mean / ld if ld > 0 else float('nan')
        print(f"{d:>4} | {mean:>9.2f} {med:>5} {p90:>5} {p99:>5} {mx:>5} | "
              f"{ld:>8.2f} {ratio:>10.3f}",
              flush=True)
        summary[d] = (mean, med, p90, p99, mx, eps)

    ds = list(summary.keys())
    means = [summary[d][0] for d in ds]
    log_ds = [math.log2(d) for d in ds]
    A = np.stack([np.array(log_ds), np.ones(len(ds))], axis=1)
    (b, a), *_ = np.linalg.lstsq(A, np.array(means), rcond=None)
    print(f"\nLinear fit:  E[N_epochs] = {a:.2f} + {b:.2f} * log2(d)", flush=True)
    log_means = [math.log2(m) for m in means]
    (b2, a2), *_ = np.linalg.lstsq(A, np.array(log_means), rcond=None)
    print(f"Log-log fit: log2(E[N_epochs]) = {a2:.2f} + {b2:.2f} * log2(d)", flush=True)
    print(f"  -> E[N_epochs] ~ d^{b2:.2f}", flush=True)

    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table (Strategy A, R=2, KS):", flush=True)
    print("=" * 70, flush=True)
    for d, (mean, med, p90, p99, mx, _) in summary.items():
        print(f"  ${d}$ & ${mean:.2f}$ & ${med}$ & ${p90}$ & ${p99}$ & ${mx}$ & "
              f"${math.log2(d):.2f}$ \\\\",
              flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
