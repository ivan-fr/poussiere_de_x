"""
Numerical verification of the phase-transition-elimination claim
(paper 101bis, Theorem 4.1).

Statement under test: under the Kostlan-Smale ensemble, the Strategy-B
first-success index s_star := min{k : z_k^(B) -> some root within the budget}
satisfies a geometric tail Pr[s_star > t] <= C * rho^t with rho in (0, 1)
universal, hence E[s_star] = O(1).

We measure on KS (NOT UNI) for d in {16, 32, 64, 128, 256, 512}.
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import (
    cauchy_root_bound, run_T2_armijo,
)
from verify_armijo_KS import sample_KS


def starts_B_full(R_root, N, q=0.7, seed=0.0):
    golden = math.pi * (3.0 - math.sqrt(5.0))
    return [R_root * (q ** k) * np.exp(1j * (k * golden + seed)) for k in range(N)]


def first_success(coefs, max_starts, max_epochs):
    """Return the smallest k such that Strategy-B start k converges; else max_starts."""
    R = cauchy_root_bound(coefs)
    starts = starts_B_full(R, max_starts)
    for k, z0 in enumerate(starts):
        try:
            st, ep, _ = run_T2_armijo(coefs, z0, max_epochs=max_epochs)
        except (OverflowError, ZeroDivisionError, FloatingPointError):
            continue  # treat as non-convergent for this start
        if st == 'converged':
            return k
    return max_starts


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    rng = np.random.default_rng(20260425)
    DEGREES = [16, 32, 64, 128, 256, 512]
    N_POLYS = 60
    MAX_STARTS = 32
    MAX_EPOCHS = 300

    t0 = time.perf_counter()
    print("=" * 110, flush=True)
    print("Phase-transition test (paper 101bis, Theorem 4.1):  s_star_B under KS", flush=True)
    print(f"  N_POLYS = {N_POLYS}, MAX_STARTS = {MAX_STARTS}, MAX_EPOCHS = {MAX_EPOCHS}",
          flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'mean s*':>9} {'med':>5} {'p90':>5} {'p99':>5} {'max':>5} | "
          f"{'Pr[s*=0]':>10} {'Pr[s*>=1]':>10} {'Pr[s*>=2]':>10} "
          f"{'rho_emp':>9}",
          flush=True)
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
        pr_ge2 = float(np.sum(s_arr >= 2)) / N
        # Geometric tail rate: log Pr[s*>=t]/log Pr[s*>=t-1]
        rhos = []
        prev = 1.0
        for t in range(1, mx + 1):
            p = float(np.sum(s_arr >= t)) / N
            if p > 0.5/N and prev > 0.5/N:
                rhos.append(p / prev)
            prev = p
        rho_emp = float(np.mean(rhos)) if rhos else float('nan')
        print(f"{d:>4} | {mean:>9.3f} {med:>5} {p90:>5} {p99:>5} {mx:>5} | "
              f"{pr0:>10.3f} {pr_ge1:>10.3f} {pr_ge2:>10.3f} "
              f"{rho_emp:>9.3f}",
              flush=True)
        summary[d] = (mean, mx, pr0, pr_ge1, pr_ge2, rho_emp, s_arr)

    # LaTeX-ready row block
    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table:", flush=True)
    print("=" * 70, flush=True)
    for d, (mean, mx, pr0, pr1, pr2, rho, _) in summary.items():
        print(f"  ${d}$ & ${mean:.2f}$ & ${mx}$ & ${pr0*100:.1f}\\%$ & ${pr1*100:.1f}\\%$ & "
              f"${pr2*100:.1f}\\%$ & ${rho:.3f}$ \\\\",
              flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
