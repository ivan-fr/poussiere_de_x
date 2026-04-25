"""
Iterated gamma-damped step: measure full convergence cost (no backtracking).
If the cost stays O(d log d) BSS ops per orbit, Option B succeeds.
"""
from __future__ import annotations
import math, time
import numpy as np

ALPHA0 = 0.1576
C_DAMP = 0.05


def eval_poly(coefs, z):
    return complex(np.polyval(coefs[::-1], z))


def gamma_damped_step(coefs, z0, d):
    """3 oracle probes -> one descending step (no backtracking)."""
    P0 = eval_poly(coefs, z0)
    aP0 = abs(P0)
    if aP0 < 1e-50 or not np.isfinite(aP0):
        return z0, aP0, 0
    eps = max(1e-6, abs(z0) * 1e-3)
    P1 = eval_poly(coefs, z0 + eps)
    P2 = eval_poly(coefs, z0 + 2*eps)
    if not (np.isfinite(P1) and np.isfinite(P2)):
        return None, float('inf'), 3
    Q1 = (P1 - P0) / eps
    if abs(Q1) < 1e-50 or not np.isfinite(Q1):
        return None, float('inf'), 3
    Q2 = (P2 - 2*P1 + P0) / (2 * eps * eps)
    gamma_loc = abs(Q2) / abs(Q1)
    beta = aP0 / abs(Q1)
    alpha = beta * gamma_loc
    if alpha < ALPHA0:
        lam = 1.0
    else:
        lam = C_DAMP / max(alpha, 1e-50)
    try:
        z_new = z0 - lam * P0 / Q1
    except (OverflowError, ZeroDivisionError):
        return None, float('inf'), 3
    if not np.isfinite(z_new):
        return None, float('inf'), 3
    Pn = eval_poly(coefs, z_new)
    aPn = abs(Pn)
    if not np.isfinite(aPn):
        return None, float('inf'), 4
    return z_new, aPn, 4   # 3 probes + 1 verification


def run_orbit(coefs, z_init, d, max_epochs=500, rel_tol=1e-12):
    coef_scale = float(np.max(np.abs(coefs)))
    if coef_scale < 1e-30: coef_scale = 1.0
    basin_tol = rel_tol * coef_scale
    z = complex(z_init)
    aP = abs(eval_poly(coefs, z))
    if aP < basin_tol: return 'converged', 0, 1, aP
    n_oracles = 1
    for ep in range(1, max_epochs + 1):
        z_new, aP_new, cost = gamma_damped_step(coefs, z, d)
        n_oracles += cost
        if z_new is None or not np.isfinite(z_new):
            return 'overflow', ep, n_oracles, aP
        # Accept any descent; if aP_new >= aP, declare stagnation
        if aP_new >= aP * 0.999:
            return 'stagnated', ep, n_oracles, aP_new
        z = z_new; aP = aP_new
        if aP < basin_tol:
            return 'converged', ep, n_oracles, aP
    return 'budget', max_epochs, n_oracles, aP


def cauchy_root_bound(coefs):
    a = np.abs(coefs); leading = a[-1]
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
    N_POLYS = 50

    print("=" * 110, flush=True)
    print("Iterated gamma-damped step (NO backtracking) on UNI -- worst-case cost per orbit", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'N orbits':>8} {'conv %':>7} {'stag %':>7} {'budget %':>9} | "
          f"{'mean orac':>10} {'med orac':>9} {'max orac':>9} | "
          f"{'/ d log d':>10} {'/ d^2 log d':>13}", flush=True)
    print("-" * 110, flush=True)

    rng = np.random.default_rng(20260424)
    summary = {}
    for d in DEGREES:
        polys = [sample_UNI(d, rng) for _ in range(N_POLYS)]
        n_conv = n_stag = n_budget = 0
        oracle_counts = []
        for coefs in polys:
            R_root = cauchy_root_bound(coefs)
            starts = starts_B(R_root, max(4, d))
            best_cost = None
            for z0 in starts[:4]:    # try up to 4 starts to find a converged orbit
                st, ep, n_ora, aP = run_orbit(coefs, z0, d,
                                              max_epochs=20 * d)
                if st == 'converged':
                    if best_cost is None or n_ora < best_cost:
                        best_cost = n_ora
                    break    # one start enough
            if best_cost is None:
                n_stag += 1
            else:
                n_conv += 1
                oracle_counts.append(best_cost)
        N = N_POLYS
        if oracle_counts:
            mean_o = float(np.mean(oracle_counts))
            med_o = int(np.median(oracle_counts))
            max_o = int(np.max(oracle_counts))
            logd = math.log(d)
            ratio_dlogd = mean_o / (d * logd)
            ratio_d2logd = mean_o / (d * d * logd)
        else:
            mean_o = med_o = max_o = -1
            ratio_dlogd = ratio_d2logd = float('nan')
        print(f"{d:>4} | {N:>8} {100*n_conv/N:>6.1f}% {100*n_stag/N:>6.1f}% "
              f"{100*n_budget/N:>8.1f}% | "
              f"{mean_o:>10.1f} {med_o:>9} {max_o:>9} | "
              f"{ratio_dlogd:>10.3f} {ratio_d2logd:>13.5f}", flush=True)
        summary[d] = (mean_o, ratio_dlogd, ratio_d2logd)

    # Fit: log(cost) ~ alpha * log(d)
    print("\nFit cost(d) ~ C * d^alpha (log d)^beta", flush=True)
    if len(summary) >= 3:
        ds = np.array(sorted(summary.keys()), dtype=float)
        cs = np.array([summary[int(d)][0] for d in ds])
        log_d = np.log(ds); log_log_d = np.log(np.log(ds))
        log_c = np.log(np.maximum(cs, 1.0))
        A = np.stack([log_d, log_log_d, np.ones_like(log_d)], axis=1)
        coef, *_ = np.linalg.lstsq(A, log_c, rcond=None)
        alpha, beta, C = coef
        print(f"  cost ~ exp({C:.2f}) * d^{alpha:.3f} * (log d)^{beta:.3f}", flush=True)
        if alpha < 1.5 and beta < 1.5:
            print(f"  ==> Empirical complexity sub-O(d^2 log d):  YES", flush=True)
        else:
            print(f"  ==> Empirical complexity exceeds O(d^2 log d): exponent {alpha:.2f}", flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
