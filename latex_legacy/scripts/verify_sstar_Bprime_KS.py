"""
Strategy B' (KS-adaptive radius) verification.

Observation: under KS, coefficients have std sqrt(binom(d, k)), so the
Cauchy bound rho = 1 + max|a_k/a_d| scales as 2^{d/2}, exponentially
loose for KS polynomials whose roots concentrate near |z| = 1
(Edelman-Kostlan).

Strategy B' fix: replace the Cauchy radius rho by R = 2, the natural
KS root radius. The spiral becomes z_k = 2 * q^k * exp(i k phi).
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import run_T2_armijo
from verify_armijo_KS import sample_KS


GOLDEN = math.pi * (3.0 - math.sqrt(5.0))


def starts_Bprime(N, R=2.0, q=0.7):
    return [R * (q ** k) * np.exp(1j * k * GOLDEN) for k in range(N)]


def first_success(coefs, max_starts, max_epochs):
    starts = starts_Bprime(max_starts)
    for k, z0 in enumerate(starts):
        try:
            st, ep, _ = run_T2_armijo(coefs, z0, max_epochs=max_epochs)
        except (OverflowError, ZeroDivisionError, FloatingPointError):
            continue
        if st == 'converged':
            return k
    return max_starts


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    rng = np.random.default_rng(20260425)
    DEGREES = [16, 32, 64, 128, 256, 512]
    N_POLYS = 60
    MAX_STARTS = 24
    MAX_EPOCHS = 300

    t0 = time.perf_counter()
    print("=" * 110, flush=True)
    print("Strategy B' (R=2, KS-adaptive) verification — paper 101bis Theorem 4.1", flush=True)
    print(f"  N_POLYS = {N_POLYS}, MAX_STARTS = {MAX_STARTS}, MAX_EPOCHS = {MAX_EPOCHS}",
          flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'mean s*':>9} {'med':>5} {'p90':>5} {'p99':>5} {'max':>5} | "
          f"{'Pr[s*=0]':>10} {'Pr[s*>=1]':>10} {'Pr[s*>=4]':>10} {'rho_emp':>9}", flush=True)
    print("-" * 110, flush=True)
    summary = {}
    for d in DEGREES:
        s_list = []
        for _ in range(N_POLYS):
            coefs = sample_KS(d, rng)
            s_list.append(first_success(coefs, MAX_STARTS, MAX_EPOCHS))
        s_arr = np.array(s_list, dtype=float)
        mean = float(np.mean(s_arr))
        med = int(np.median(s_arr))
        p90 = int(np.percentile(s_arr, 90))
        p99 = int(np.percentile(s_arr, 99))
        mx = int(np.max(s_arr))
        N = len(s_arr)
        pr0 = float(np.sum(s_arr == 0)) / N
        pr_ge1 = float(np.sum(s_arr >= 1)) / N
        pr_ge4 = float(np.sum(s_arr >= 4)) / N
        rhos = []
        prev = 1.0
        for t in range(1, mx + 2):
            p = float(np.sum(s_arr >= t)) / N
            if p > 0.5/N and prev > 0.5/N:
                rhos.append(p / prev)
            prev = p
        rho_emp = float(np.mean(rhos)) if rhos else float('nan')
        print(f"{d:>4} | {mean:>9.3f} {med:>5} {p90:>5} {p99:>5} {mx:>5} | "
              f"{pr0:>10.3f} {pr_ge1:>10.3f} {pr_ge4:>10.3f} {rho_emp:>9.3f}",
              flush=True)
        summary[d] = (mean, mx, pr0, pr_ge1, pr_ge4, rho_emp, s_arr)

    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table:", flush=True)
    print("=" * 70, flush=True)
    for d, (mean, mx, pr0, pr1, pr4, rho, _) in summary.items():
        print(f"  ${d}$ & ${mean:.2f}$ & ${mx}$ & ${pr0*100:.1f}\\%$ & ${pr1*100:.1f}\\%$ & "
              f"${pr4*100:.1f}\\%$ & ${rho:.3f}$ \\\\",
              flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
