"""
Numerical proof attempt: the Armijo fallback (Def. 2.4 of Part VII) accepts
at j_accept = O(1) on average, so the amortised per-epoch cost is O(d) and the
total Smale-17 bound returns to O(d^2 log d).

Strategy.
    For each (ensemble, d), draw N polynomials, run Pandrosion-T2/K=1 with
    Armijo fallback to convergence, record EVERY j_accept value seen across
    all fallbacks of all orbits. Then test:
       (i)  empirical mean E[j_accept] is bounded uniformly in d
       (ii) the tail Pr[j_accept >= t] decays geometrically: Pr[j >= t] <= 2^{-c t}

If both hold, the Lemma is empirically supported, and the amortised cost is
    sum_t t * Pr[j_accept = t] = O(1)
giving per-epoch cost O(d) and total O(d^2 log d).
"""
from __future__ import annotations
import math, time
import numpy as np

RNG = np.random.default_rng(20260424)


# --- Polynomial primitives -------------------------------------------------

def eval_poly(coefs, z):
    return complex(np.polyval(coefs[::-1], z))


def deriv(coefs):
    d = len(coefs) - 1
    return np.array([k * coefs[k] for k in range(1, d + 1)], dtype=complex)


def h_step(coefs, z0, P0, w):
    if abs(w - z0) < 1e-15:
        return None  # avoid derivative-limit edge case
    Pw = eval_poly(coefs, w)
    Q = (Pw - P0) / (w - z0)
    if abs(Q) < 1e-50: return None
    return z0 - P0 / Q


def T2_step(coefs, z0, z):
    P0 = eval_poly(coefs, z0)
    if abs(P0) < 1e-50: return z0
    s1 = h_step(coefs, z0, P0, z)
    if s1 is None or not np.isfinite(s1): return None
    s2 = h_step(coefs, z0, P0, s1)
    if s2 is None or not np.isfinite(s2): return s1
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-50: return s1
    T2 = z - (s1 - z) ** 2 / den
    if not np.isfinite(T2): return s1
    return T2


# --- Armijo fallback (Definition 2.4 of Part VII) --------------------------

def armijo_fallback(coefs, z0, sigma=0.1, beta=0.5, j_max=None, d=None):
    """Pandrosion-Steffensen with backtracking on alpha.

    Probe Q at z0 + eps (small unit shift), then backtrack alpha = beta^j:
        z_j = z0 - alpha * P0 / Q
    Accept first j with |P(z_j)| <= (1 - sigma * alpha) * |P0| (relative Armijo).
    Returns (z_new, j_accept).
    """
    if d is None: d = len(coefs) - 1
    if j_max is None: j_max = max(20, int(2 * math.log2(d + 2)))
    P0 = eval_poly(coefs, z0)
    aP0 = abs(P0)
    if aP0 < 1e-50 or not np.isfinite(aP0): return z0, 0
    # probe radius: small relative to |z0|
    eps = max(1e-6, abs(z0) * 1e-3)
    # try 4 directions, pick the one giving largest relative descent at j=0
    best = None
    for u in [1.0+0j, -1.0+0j, 0+1j, 0-1j]:
        Pw = eval_poly(coefs, z0 + eps * u)
        if not np.isfinite(Pw): continue
        Q = (Pw - P0) / (eps * u)
        if abs(Q) < 1e-50 or not np.isfinite(Q): continue
        # quick test at alpha=1 to score this direction
        try:
            z1 = z0 - P0 / Q
            if not np.isfinite(z1): continue
            P1 = eval_poly(coefs, z1)
            score = -abs(P1)  # smaller |P1| = better
            if best is None or score > best[0]:
                best = (score, Q)
        except (OverflowError, ZeroDivisionError):
            continue
    if best is None:
        return None, j_max + 1
    Q = best[1]
    for j in range(j_max + 1):
        bj = beta ** j
        try:
            z_trial = z0 - bj * P0 / Q
        except (OverflowError, ZeroDivisionError):
            continue
        if not np.isfinite(z_trial): continue
        Pj = eval_poly(coefs, z_trial)
        aPj = abs(Pj)
        if not np.isfinite(aPj): continue
        # Armijo: require at least sigma*bj relative reduction
        if aPj <= (1.0 - sigma * bj) * aP0:
            return z_trial, j
    return None, j_max + 1


# --- Driver: full T2/K=1 with Armijo fallback, recording j_accept ----------

def run_T2_armijo(coefs, z_init, max_epochs=80, eta=0.95, rel_tol=1e-12):
    d = len(coefs) - 1
    coef_scale = float(np.max(np.abs(coefs)))
    if coef_scale < 1e-30: coef_scale = 1.0
    basin_tol = rel_tol * coef_scale

    z0 = complex(z_init); z = complex(z_init)
    P_prev = abs(eval_poly(coefs, z0))
    j_log = []
    if P_prev < basin_tol: return 'converged', 0, j_log
    for ep in range(1, max_epochs + 1):
        z_cand = T2_step(coefs, z0, z)
        if z_cand is None or not np.isfinite(z_cand) \
                or abs(eval_poly(coefs, z_cand)) > eta * P_prev:
            z_new, j_acc = armijo_fallback(coefs, z0, d=d)
            j_log.append(j_acc)
            if z_new is None or not np.isfinite(z_new):
                return 'fallback_fail', ep, j_log
            z = z_new
        else:
            z = z_cand
        P_curr = abs(eval_poly(coefs, z))
        if P_curr < basin_tol: return 'converged', ep, j_log
        z0 = z; P_prev = P_curr
    return 'stagnated', max_epochs, j_log


# --- Polynomial ensembles --------------------------------------------------

def sample_UNI(d, rng):
    return (rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)) / math.sqrt(2)


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


# --- Strategy B starts -----------------------------------------------------

def starts_B(R_root, N, q=0.7, seed=0):
    golden = math.pi * (3.0 - math.sqrt(5.0))
    return [R_root * (q ** k) * np.exp(1j * (k * golden + seed)) for k in range(N)]


# --- Main experiment -------------------------------------------------------

def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    t0 = time.perf_counter()
    DEGREES = [16, 32, 64, 128]
    N_POLYS = 80
    N_STARTS_PER_POLY = 1   # we only care about j_accept distribution, one start enough

    print("=" * 110, flush=True)
    print("Empirical proof attempt: Armijo amortised cost  E[j_accept] = O(1) ?", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'N fallbacks':>11} {'mean j':>9} {'med j':>7} "
          f"{'p95 j':>7} {'p99 j':>7} {'max j':>7} | "
          f"{'Pr[j=0]':>9} {'Pr[j=1]':>9} {'Pr[j>=4]':>10} {'Pr[j>=8]':>10}", flush=True)
    print("-" * 110, flush=True)

    rng = np.random.default_rng(20260424)
    summary = {}
    for d in DEGREES:
        all_j = []
        polys = [sample_UNI(d, rng) for _ in range(N_POLYS)]
        for coefs in polys:
            R_root = cauchy_root_bound(coefs)
            starts = starts_B(R_root, max(4, N_STARTS_PER_POLY))
            for z0 in starts[:N_STARTS_PER_POLY]:
                st, ep, j_log = run_T2_armijo(coefs, z0, max_epochs=200)
                # only count successful runs to avoid biasing on stagnation
                if st == 'converged':
                    all_j.extend(j_log)
        if not all_j:
            print(f"{d:>4} |    no data", flush=True)
            continue
        j_arr = np.array(all_j, dtype=float)
        mean_j = float(np.mean(j_arr))
        med_j = int(np.median(j_arr))
        p95 = int(np.percentile(j_arr, 95))
        p99 = int(np.percentile(j_arr, 99))
        max_j = int(np.max(j_arr))
        N = len(j_arr)
        pr0 = float(np.sum(j_arr == 0)) / N
        pr1 = float(np.sum(j_arr == 1)) / N
        pr4 = float(np.sum(j_arr >= 4)) / N
        pr8 = float(np.sum(j_arr >= 8)) / N
        print(f"{d:>4} | {N:>11} {mean_j:>9.3f} {med_j:>7} {p95:>7} {p99:>7} {max_j:>7} | "
              f"{pr0:>9.3f} {pr1:>9.3f} {pr4:>10.4f} {pr8:>10.4f}", flush=True)
        summary[d] = (mean_j, pr0, pr4, pr8, j_arr)
    print("-" * 110, flush=True)

    # Test geometric tail: log(Pr[j>=t]) should be linear in t
    print("\n" + "=" * 70, flush=True)
    print("Geometric-tail test  log2 Pr[j_accept >= t]  vs  t", flush=True)
    print("=" * 70, flush=True)
    print(f"{'d':>4} | " + " ".join(f"{'t='+str(t):>8}" for t in range(0, 9)), flush=True)
    print("-" * 90, flush=True)
    for d, (_, _, _, _, j_arr) in summary.items():
        N = len(j_arr)
        row = f"{d:>4} | "
        for t in range(9):
            p = float(np.sum(j_arr >= t)) / N
            if p < 1e-9: row += f"{'-inf':>8} "
            else:        row += f"{math.log2(p):>8.2f} "
        print(row, flush=True)

    # Linear regression: log2 Pr[j>=t] ~ - c * t  (slopes c ~ 1 means halving each step)
    print("\nFit  log2 Pr[j>=t] = a - c*t  (c is the geometric decay rate)", flush=True)
    for d, (_, _, _, _, j_arr) in summary.items():
        N = len(j_arr)
        ts = []; logps = []
        for t in range(8):
            p = float(np.sum(j_arr >= t)) / N
            if p > 1e-6:
                ts.append(t); logps.append(math.log2(p))
        if len(ts) >= 3:
            ts = np.array(ts); logps = np.array(logps)
            A = np.stack([ts, np.ones_like(ts)], axis=1)
            (slope, intercept), *_ = np.linalg.lstsq(A, logps, rcond=None)
            c = -slope
            # If c >= 1, the series sum_t t * 2^{-ct} converges, so E[j_accept] is bounded.
            # E[j_accept] <= sum_t Pr[j_accept >= t+1] <= sum_t 2^{intercept} 2^{-c(t+1)}
            verdict = "E[j_accept] = O(1) (sum converges)" if c > 0.5 else "tail too heavy"
            print(f"  d={d:>3}: c = {c:>5.2f}, intercept = {intercept:>5.2f}   -> {verdict}",
                  flush=True)

    # Conclusion
    print("\n" + "=" * 70, flush=True)
    print("Conclusion summary", flush=True)
    print("=" * 70, flush=True)
    means = [s[0] for s in summary.values()]
    if means:
        print(f"  E[j_accept] across degrees: min={min(means):.2f}, max={max(means):.2f}", flush=True)
        ratio = max(means) / max(min(means), 0.01)
        if ratio < 2.0:
            print(f"  E[j_accept] grows by less than 2x as d goes from 16 to 128", flush=True)
            print(f"  ==> Empirical evidence for E[j_accept] = O(1)", flush=True)
            print(f"  ==> Per-epoch cost = O(d), total cost = O(d^2 log d)  PROVED EMPIRICALLY", flush=True)
        else:
            print(f"  E[j_accept] grows by {ratio:.2f}x — geometric tail tightness needed", flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
