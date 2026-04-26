"""
Same epoch-count measurement as verify_epochs_KS.py, but using the
Pandrosion-T_3 / K=3 scheme of Part VII (the configuration before
the Part-VIII T_2/K=1 simplification).

T_n hierarchy (Part VIII, Sec 3):
    s_1 = h(z), s_2 = h(s_1), s_3 = h(s_2)
    T_2(z) = z - (s_1 - z)^2 / (s_2 - 2 s_1 + z)
    T_3(z) = T_2(z) - (s_3 - T_2(z)) / (lambda_hat - 1)
where  lambda_hat = (s_2 - s_1) / (s_1 - z).

K=3 anchor schedule: keep z_0 fixed for 3 successive epochs, then re-anchor.
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import (
    eval_poly, h_step, armijo_fallback,
)
from verify_armijo_KS import sample_KS

GOLDEN = math.pi * (3.0 - math.sqrt(5.0))


def T3_step(coefs, z0, z):
    """T_3 extrapolation built from three h-iterates anchored at z0."""
    P0 = eval_poly(coefs, z0)
    if abs(P0) < 1e-50:
        return z0
    s1 = h_step(coefs, z0, P0, z)
    if s1 is None or not np.isfinite(s1):
        return None
    s2 = h_step(coefs, z0, P0, s1)
    if s2 is None or not np.isfinite(s2):
        return s1
    s3 = h_step(coefs, z0, P0, s2)
    den2 = s2 - 2 * s1 + z
    if abs(den2) < 1e-50 or not np.isfinite(den2):
        return s1
    T2 = z - (s1 - z) ** 2 / den2
    if not np.isfinite(T2):
        return s1
    if s3 is None or not np.isfinite(s3):
        return T2
    den_l = s1 - z
    if abs(den_l) < 1e-50 or not np.isfinite(den_l):
        return T2
    lam = (s2 - s1) / den_l
    if abs(lam - 1) < 1e-50 or not np.isfinite(lam):
        return T2
    T3 = T2 - (s3 - T2) / (lam - 1)
    if not np.isfinite(T3):
        return T2
    return T3


def run_T3_K3_armijo_count_epochs(coefs, z_init, max_epochs=400, eta=0.95, rel_tol=1e-12):
    """Pandrosion-T_3 with anchor period K=3, Armijo fallback. Returns epoch count."""
    d = len(coefs) - 1
    coef_scale = float(np.max(np.abs(coefs)))
    if coef_scale < 1e-30:
        coef_scale = 1.0
    basin_tol = rel_tol * coef_scale
    K = 3

    z0 = complex(z_init)
    z = complex(z_init)
    P_prev = abs(eval_poly(coefs, z0))
    if P_prev < basin_tol:
        return 'converged', 0
    epoch_in_anchor = 0
    for ep in range(1, max_epochs + 1):
        # T_3 step from current anchor z0 and current iterate z
        try:
            z_cand = T3_step(coefs, z0, z)
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
        # K=3 anchor schedule: re-anchor every 3 epochs
        epoch_in_anchor += 1
        if epoch_in_anchor >= K:
            z0 = z
            P_prev = P_curr
            epoch_in_anchor = 0
        else:
            P_prev = P_curr
    return 'stagnated', max_epochs


def starts_Bprime(N, R=2.0, q=0.7):
    return [R * (q ** k) * np.exp(1j * k * GOLDEN) for k in range(N)]


def measure(d, n_polys, max_starts=24, max_epochs=400, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260426 + d)
    epoch_counts = []
    for _ in range(n_polys):
        coefs = sample_KS(d, rng)
        starts = starts_Bprime(max_starts)
        for k, z0 in enumerate(starts):
            try:
                st, ep = run_T3_K3_armijo_count_epochs(coefs, z0, max_epochs=max_epochs)
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
    print("Pandrosion-T_3 / K=3 — KS epoch-count measurement (Part VII config)",
          flush=True)
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
    print("LaTeX-ready table (T_3, K=3, B' starts, KS):", flush=True)
    print("=" * 70, flush=True)
    for d, (mean, med, p90, p99, mx, _) in summary.items():
        print(f"  ${d}$ & ${mean:.2f}$ & ${med}$ & ${p90}$ & ${p99}$ & ${mx}$ & "
              f"${math.log2(d):.2f}$ \\\\",
              flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
