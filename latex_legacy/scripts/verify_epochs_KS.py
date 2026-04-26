"""
Numerical verification for paper 101 (ABD conjecture):
how many Pandrosion-T2/K=1 epochs are needed per orbit to reach the basin
and then machine precision, on KS-distributed monic polynomials, with
Strategy B' starts (radius=2)?

Hypothesis: E[N_epochs per orbit] = O(log d) on KS.
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import (
    eval_poly, T2_step, armijo_fallback,
)
from verify_armijo_KS import sample_KS

GOLDEN = math.pi * (3.0 - math.sqrt(5.0))


def run_T2_armijo_count_epochs(coefs, z_init, max_epochs=300, eta=0.95, rel_tol=1e-12):
    """Modified run_T2_armijo that returns the number of epochs to convergence."""
    d = len(coefs) - 1
    coef_scale = float(np.max(np.abs(coefs)))
    if coef_scale < 1e-30:
        coef_scale = 1.0
    basin_tol = rel_tol * coef_scale

    z0 = complex(z_init); z = complex(z_init)
    P_prev = abs(eval_poly(coefs, z0))
    epochs_to_converge = 0
    if P_prev < basin_tol:
        return 'converged', 0
    for ep in range(1, max_epochs + 1):
        try:
            z_cand = T2_step(coefs, z0, z)
        except (OverflowError, ZeroDivisionError):
            z_cand = None
        if z_cand is None or not np.isfinite(z_cand):
            try:
                z_new, j_acc = armijo_fallback(coefs, z0, d=d)
            except (OverflowError, ZeroDivisionError):
                return 'fail', ep
            if z_new is None or not np.isfinite(z_new):
                return 'fail', ep
            z = z_new
        else:
            try:
                pcand = abs(eval_poly(coefs, z_cand))
            except (OverflowError, ZeroDivisionError):
                pcand = float('inf')
            if pcand > eta * P_prev:
                try:
                    z_new, j_acc = armijo_fallback(coefs, z0, d=d)
                except (OverflowError, ZeroDivisionError):
                    return 'fail', ep
                if z_new is None or not np.isfinite(z_new):
                    return 'fail', ep
                z = z_new
            else:
                z = z_cand
        try:
            P_curr = abs(eval_poly(coefs, z))
        except (OverflowError, ZeroDivisionError):
            return 'fail', ep
        if P_curr < basin_tol:
            return 'converged', ep
        z0 = z; P_prev = P_curr
    return 'stagnated', max_epochs


def starts_Bprime(N, R=2.0, q=0.7):
    return [R * (q ** k) * np.exp(1j * k * GOLDEN) for k in range(N)]


def measure(d, n_polys, max_starts=24, max_epochs=400, rng=None):
    """Returns list of epoch counts for the first successful orbit per poly."""
    if rng is None:
        rng = np.random.default_rng(20260425 + d)
    epoch_counts = []
    for _ in range(n_polys):
        coefs = sample_KS(d, rng)
        starts = starts_Bprime(max_starts)
        for k, z0 in enumerate(starts):
            try:
                st, ep = run_T2_armijo_count_epochs(coefs, z0, max_epochs=max_epochs)
            except (OverflowError, ZeroDivisionError):
                continue
            if st == 'converged':
                epoch_counts.append(ep)
                break
        else:
            epoch_counts.append(max_epochs)  # never converged
    return np.array(epoch_counts, dtype=float)


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    DEGREES = [16, 32, 64, 128, 256, 512]
    N_POLYS = 60
    t0 = time.perf_counter()
    print("=" * 110, flush=True)
    print("Epoch-count measurement on KS  (paper 101, ABD conjecture)", flush=True)
    print(f"  N_POLYS = {N_POLYS}", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'mean ep':>9} {'med':>5} {'p90':>5} {'p99':>5} {'max':>5} | "
          f"{'log2(d)':>8} {'mean/log2':>10}", flush=True)
    print("-" * 110, flush=True)
    summary = {}
    for d in DEGREES:
        eps = measure(d, N_POLYS)
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

    # Test: is E[N_epochs] = O(log d)? Fit mean = a + b * log(d)
    ds = list(summary.keys())
    means = [summary[d][0] for d in ds]
    log_ds = [math.log2(d) for d in ds]
    A = np.stack([np.array(log_ds), np.ones(len(ds))], axis=1)
    (b, a), *_ = np.linalg.lstsq(A, np.array(means), rcond=None)
    print(f"\nLinear fit:  E[N_epochs] = {a:.2f} + {b:.2f} * log2(d)", flush=True)
    print(f"  -> E[N_epochs] = O(log d) iff b > 0 and bounded.", flush=True)

    # Also fit log-log: is it polynomial in d?
    log_means = [math.log2(m) for m in means]
    (b2, a2), *_ = np.linalg.lstsq(
        np.stack([np.array(log_ds), np.ones(len(ds))], axis=1),
        np.array(log_means), rcond=None)
    print(f"Log-log fit: log2(E[N_epochs]) = {a2:.2f} + {b2:.2f} * log2(d)", flush=True)
    print(f"  -> E[N_epochs] ~ d^{b2:.2f}", flush=True)

    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table:", flush=True)
    print("=" * 70, flush=True)
    for d, (mean, med, p90, p99, mx, _) in summary.items():
        print(f"  ${d}$ & ${mean:.2f}$ & ${med}$ & ${p90}$ & ${p99}$ & ${mx}$ & "
              f"${math.log2(d):.2f}$ \\\\",
              flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
