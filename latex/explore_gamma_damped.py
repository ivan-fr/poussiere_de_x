"""
Option B: replace Armijo backtracking with a gamma-damped step computed analytically
from Smale alpha-theory.

Idea:
    With 3 Steffensen probes at z0, z0+eps, z0+2*eps, we compute:
        Q1   ~~  P'(z0)            from (P(z0+eps) - P0) / eps
        Q2   ~~  P''(z0)/2         from (P(z0+2eps) - 2*P(z0+eps) + P0) / eps^2 / 2
        gamma_local := |Q2/Q1|     a 1-term lower bound on Smale gamma
        beta := |P0/Q1|
        alpha := beta * gamma_local
    Then choose the damping:
        if alpha < alpha0 (~ 0.157, Smale's alpha-test): lambda = 1 (full Newton)
        else:                                            lambda = c / alpha  (e.g. c = 0.05)
    Apply z_new = z0 - lambda * P0 / Q1.

Theorem (Smale 1986, restated): with this damping the descent factor is universal:
    |P(z_new)| <= (1 - lambda/2) * |P0|
in O(1) probes (no backtracking). Total per-epoch cost: O(d), giving worst-case O(d^2 log d).

We empirically verify: does the damped step always descend by at least factor (1 - c)?
"""
from __future__ import annotations
import math, time
import numpy as np

RNG = np.random.default_rng(20260424)

ALPHA0 = 0.1576       # Smale's alpha test constant, (3 - 2 sqrt 2)/2
C_DAMP = 0.05         # constant in the damping formula


def eval_poly(coefs, z):
    return complex(np.polyval(coefs[::-1], z))


def gamma_damped_step(coefs, z0, d=None):
    """Compute one gamma-damped Newton-Steffensen step from anchor z0.
    Returns (z_new, descent_factor, alpha_est, used_damping_bool)
    where descent_factor = |P(z_new)| / |P(z0)| (smaller is better).
    Cost: 3 probe evaluations of P (constant)."""
    if d is None: d = len(coefs) - 1
    P0 = eval_poly(coefs, z0)
    aP0 = abs(P0)
    if aP0 < 1e-50 or not np.isfinite(aP0):
        return z0, 0.0, 0.0, False
    eps = max(1e-6, abs(z0) * 1e-3)
    P1 = eval_poly(coefs, z0 + eps)
    P2 = eval_poly(coefs, z0 + 2*eps)
    if not (np.isfinite(P1) and np.isfinite(P2)):
        return None, float('inf'), float('inf'), False
    # Q1 ~~ P'(z0)
    Q1 = (P1 - P0) / eps
    if abs(Q1) < 1e-50 or not np.isfinite(Q1):
        return None, float('inf'), float('inf'), False
    # Q2 ~~ P''(z0)/2
    Q2 = (P2 - 2*P1 + P0) / (2 * eps * eps)
    gamma_loc = abs(Q2) / abs(Q1)
    beta = aP0 / abs(Q1)
    alpha = beta * gamma_loc
    # Choose damping
    if alpha < ALPHA0:
        lam = 1.0
        used = False
    else:
        lam = C_DAMP / max(alpha, 1e-50)
        used = True
    # Step
    try:
        z_new = z0 - lam * P0 / Q1
    except (OverflowError, ZeroDivisionError):
        return None, float('inf'), alpha, used
    if not np.isfinite(z_new):
        return None, float('inf'), alpha, used
    Pn = eval_poly(coefs, z_new)
    aPn = abs(Pn)
    if not np.isfinite(aPn):
        return None, float('inf'), alpha, used
    return z_new, aPn / aP0, alpha, used


def cauchy_root_bound(coefs):
    a = np.abs(coefs)
    leading = a[-1]
    if leading < 1e-30: return 10.0
    d = len(coefs) - 1
    if d == 0: return 1.0
    ratios = a[:-1] / leading
    cauchy = 1.0 + float(np.max(ratios))
    expos = 1.0 / np.arange(d, 0, -1)
    safe = np.where(ratios > 1e-30, ratios, 1e-30)
    fuj = 2.0 * float(np.max(safe ** expos))
    return min(cauchy, fuj)


def starts_B(R_root, N, q=0.7):
    golden = math.pi * (3.0 - math.sqrt(5.0))
    return [R_root * (q ** k) * np.exp(1j * (k * golden)) for k in range(N)]


def sample_UNI(d, rng):
    return (rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)) / math.sqrt(2)


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    t0 = time.perf_counter()

    DEGREES = [16, 32, 64, 128]
    N_POLYS = 100
    N_PROBE_POINTS = 8       # number of Cauchy-circle anchors per polynomial to sample

    print("=" * 110, flush=True)
    print("Gamma-damped step empirical test  --  does the analytic damping give descent in 1 probe-set?", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'N samples':>10} {'descent %':>10} "
          f"{'mean log10 ratio':>17} {'median log10 ratio':>19} | "
          f"{'Pr[alpha<a0]':>13} {'Pr[lam=1]':>10}", flush=True)
    print("-" * 110, flush=True)

    rng = np.random.default_rng(20260424)
    summary = {}
    for d in DEGREES:
        ratios = []
        alphas = []
        used_count = 0
        descent_count = 0
        polys = [sample_UNI(d, rng) for _ in range(N_POLYS)]
        for coefs in polys:
            R_root = cauchy_root_bound(coefs)
            for z0 in starts_B(R_root, N_PROBE_POINTS):
                z_new, ratio, alpha, used = gamma_damped_step(coefs, z0, d=d)
                if not np.isfinite(ratio): continue
                ratios.append(ratio)
                alphas.append(alpha)
                if not used: used_count += 1
                if ratio < 0.99: descent_count += 1
        if not ratios:
            print(f"{d:>4} | no valid samples", flush=True)
            continue
        N = len(ratios)
        ratios_arr = np.array(ratios)
        alphas_arr = np.array(alphas)
        descent_pct = 100 * descent_count / N
        mean_log = float(np.mean(np.log10(np.maximum(ratios_arr, 1e-300))))
        med_log = float(np.median(np.log10(np.maximum(ratios_arr, 1e-300))))
        pr_smale = 100 * float(np.sum(alphas_arr < ALPHA0)) / N
        pr_full = 100 * used_count / N
        print(f"{d:>4} | {N:>10} {descent_pct:>9.1f}% "
              f"{mean_log:>17.3f} {med_log:>19.3f} | "
              f"{pr_smale:>12.1f}% {pr_full:>9.1f}%", flush=True)
        summary[d] = (descent_pct, ratios_arr, alphas_arr)

    # Distribution of descent factors
    print("\n" + "=" * 110, flush=True)
    print("Distribution of descent ratio  |P(z_new)| / |P(z0)|  (smaller = better; 1.0 = stagnation)", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'< 0.1':>7} {'< 0.5':>7} {'< 0.9':>7} {'< 1.0':>7} {'< 2.0':>7} {'>= 2':>7} {'>= 100':>8}",
          flush=True)
    for d, (_, ratios_arr, _) in summary.items():
        N = len(ratios_arr)
        thr = [0.1, 0.5, 0.9, 1.0, 2.0]
        cnts = [100 * float(np.sum(ratios_arr < t)) / N for t in thr]
        ge2 = 100 * float(np.sum(ratios_arr >= 2.0)) / N
        ge100 = 100 * float(np.sum(ratios_arr >= 100.0)) / N
        row = f"{d:>4} | "
        for c in cnts: row += f"{c:>6.1f}% "
        row += f"{ge2:>6.1f}% {ge100:>7.1f}%"
        print(row, flush=True)

    # Verdict
    print("\n" + "=" * 70, flush=True)
    print("Verdict for Option B  (worst-case O(d^2 log d) without backtracking)", flush=True)
    print("=" * 70, flush=True)
    for d, (descent_pct, _, _) in summary.items():
        if descent_pct >= 95:
            verdict = "+++  Damping reliably descends (rare backtrack needed)"
        elif descent_pct >= 70:
            verdict = "+    Damping descends most of the time, hybrid needed"
        else:
            verdict = "-    Damping insufficient; backtracking still required"
        print(f"  d = {d:>3}: descent in {descent_pct:5.1f}% of samples   {verdict}", flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
