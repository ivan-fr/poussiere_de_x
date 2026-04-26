"""
Numerical verification of Theorem 4.1 of paper 101ter:
on the Kostlan-Smale ensemble, Pr[j_accept >= t] <= 2 C_2 * 2^{-t}.

We measure on the KS ensemble (NOT UNI) the geometric decay rate c_emp and
report it alongside the theoretical lower bound c_thm = 1 (Lemma 3.3).
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

# Reuse the algorithmic primitives of Part VIII verbatim
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import (
    eval_poly, T2_step, armijo_fallback, run_T2_armijo,
    cauchy_root_bound, starts_B,
)


# --- Kostlan-Smale ensemble ------------------------------------------------

def sample_KS(d, rng):
    """Bombieri-Weyl-normalized monic polynomial of degree d.

    Coefficients a_k for k = 0..d are i.i.d. complex Gaussian with variance
    binomial(d, k); we then divide by the leading coefficient to get monic.
    """
    # Use lgamma to avoid integer overflow for binomial coefficients at large d
    sigmas = np.exp(0.5 * np.array(
        [math.lgamma(d + 1) - math.lgamma(k + 1) - math.lgamma(d - k + 1)
         for k in range(d + 1)]))
    re = rng.standard_normal(d + 1) * sigmas
    im = rng.standard_normal(d + 1) * sigmas
    coefs = (re + 1j * im) / math.sqrt(2)
    # Monic normalisation: divide by leading coefficient
    if abs(coefs[-1]) > 1e-30:
        coefs = coefs / coefs[-1]
    return coefs


# --- Driver ----------------------------------------------------------------

def run(degrees, n_polys=80, max_epochs=300, max_starts=12):
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    rng = np.random.default_rng(20260425)
    summary = {}
    print("=" * 100, flush=True)
    print(f"{'d':>4} | {'N':>7} {'Pr[j=0]':>9} {'mean j':>8} {'max j':>6} | "
          f"{'c_emp':>7} {'c_thm':>7} {'verdict':>15}", flush=True)
    print("-" * 100, flush=True)
    for d in degrees:
        all_j = []
        n_converged = 0
        for _ in range(n_polys):
            coefs = sample_KS(d, rng)
            R = cauchy_root_bound(coefs)
            starts = starts_B(R, max_starts)
            # Take first start that converges; record j_accept logs from that orbit
            for z0 in starts:
                st, ep, j_log = run_T2_armijo(coefs, z0, max_epochs=max_epochs)
                if st == 'converged':
                    all_j.extend(j_log)
                    n_converged += 1
                    break
        if not all_j:
            print(f"{d:>4} | no convergent runs", flush=True)
            continue
        j_arr = np.array(all_j, dtype=float)
        N = len(j_arr)
        pr0 = float(np.sum(j_arr == 0)) / N
        mean_j = float(np.mean(j_arr))
        max_j = int(np.max(j_arr))
        # Geometric-tail regression on ts where Pr[j>=t] >= 0.5/N (i.e. at least 1 obs)
        ts, logps = [], []
        thresh = 0.5 / N
        for t in range(0, max_j + 2):
            p = float(np.sum(j_arr >= t)) / N
            if p > thresh:
                ts.append(t); logps.append(math.log2(p))
        c_emp = float('nan')
        if len(ts) >= 2:
            A = np.stack([np.array(ts, float), np.ones(len(ts))], axis=1)
            (slope, _), *_ = np.linalg.lstsq(A, np.array(logps), rcond=None)
            c_emp = -slope
        verdict = "OK (c_emp >= 1)" if c_emp >= 1.0 else "loose"
        print(f"{d:>4} | {N:>7} {pr0:>9.3f} {mean_j:>8.4f} {max_j:>6} | "
              f"{c_emp:>7.2f} {1.0:>7.2f} {verdict:>15}  [{n_converged}/{n_polys} converged]",
              flush=True)
        summary[d] = (pr0, c_emp, mean_j, max_j, j_arr)
    return summary


def main():
    t0 = time.perf_counter()
    print("KS-ensemble verification of Theorem 4.1 (paper 101ter)", flush=True)
    print(f"Predicted geometric decay rate c_thm = 1 (Lemma 3.3)\n", flush=True)
    summary = run([16, 32, 64, 128, 256], n_polys=80)
    # Print the LaTeX-ready row block for §5 of the paper
    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table for §5:", flush=True)
    print("=" * 70, flush=True)
    for d, (pr0, c_emp, mean_j, max_j, _) in summary.items():
        print(f"  ${d}$  & ${pr0*100:.1f}\\%$ & ${c_emp:.2f}$ & $1$ & ${mean_j:.3f}$ \\\\",
              flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
