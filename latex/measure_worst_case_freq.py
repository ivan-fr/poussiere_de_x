"""
Measure the empirical frequency of worst-case behaviour on RANDOM polynomials,
under the two start strategies A:Cauchy and B:GeomDecr.

Question: what fraction of polynomials require close to the worst-case d starts
(saturating the O(d^2 log d) bound) versus few starts (O(d log d) regime)?

We use TWO random ensembles:
    KS  -- Kostlan-Smale, the natural unitarily-invariant complex Gaussian
            P(z) = sum_k a_k z^k with Var(a_k) = C(d, k).
    UNI -- coefficients i.i.d. complex standard Gaussian (a less natural but
           common ensemble in the literature; harder for A).

For each (ensemble, strategy, d), we draw N_SAMPLES polynomials and bin s*
(number of starts before first basin entry) into:
    "fast" : s* <= log d
    "med"  : log d  <  s* <= sqrt(d)
    "slow" : sqrt(d) < s* <= d
    "fail" : never converges within d starts
The "slow" + "fail" buckets ARE the worst-case fraction (saturates O(d^2 log d)).
"""
from __future__ import annotations
import math, time
import numpy as np

RNG = np.random.default_rng(20260424)


# ---- Pandrosion T2/K=1 with fallback (winner) ------------------------------

class C:
    n = 0


def reset(): C.n = 0


def eval_poly_coefs(coefs, z):
    # Horner
    acc = 0j
    for c in coefs[::-1]:
        acc = acc * z + c
    return acc


def h_step(coefs, z0, P0, w):
    C.n += 1
    if abs(w - z0) < 1e-15:
        # P'(z0) via Horner derivative
        d = len(coefs) - 1
        acc = 0j
        for k in range(d, 0, -1):
            acc = acc * z0 + k * coefs[k]
        Q = acc
    else:
        Pw = eval_poly_coefs(coefs, w)
        Q = (Pw - P0) / (w - z0)
    if abs(Q) < 1e-50:
        return None
    return z0 - P0 / Q


def T2_step(coefs, z0, z):
    P0 = eval_poly_coefs(coefs, z0)
    if abs(P0) < 1e-50: return z0, P0
    s1 = h_step(coefs, z0, P0, z)
    if s1 is None or not np.isfinite(s1): return None, P0
    s2 = h_step(coefs, z0, P0, s1)
    if s2 is None or not np.isfinite(s2): return s1, P0
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-50: return s1, P0
    T2 = z - (s1 - z) ** 2 / den
    if not np.isfinite(T2): return s1, P0
    return T2, P0


def steffensen(coefs, z0, alpha=0.25):
    P0 = eval_poly_coefs(coefs, z0)
    if abs(P0) < 1e-50: return z0
    z_shift = z0 + P0
    C.n += 1
    Pshift = eval_poly_coefs(coefs, z_shift)
    Q = (Pshift - P0) / P0
    if abs(Q) < 1e-50: return z0 - alpha * P0
    return z0 - alpha * P0 / Q


def run_T2_K1(coefs, z_init, max_epochs=80, eta=0.95, basin_tol=1e-12):
    z0 = complex(z_init); z = complex(z_init)
    P_prev = abs(eval_poly_coefs(coefs, z0))
    if P_prev < basin_tol: return 'converged', 0
    for ep in range(1, max_epochs + 1):
        z_cand, _ = T2_step(coefs, z0, z)
        if z_cand is None or not np.isfinite(z_cand):
            z = steffensen(coefs, z0)
        else:
            P_cand = abs(eval_poly_coefs(coefs, z_cand))
            z = z_cand if P_cand <= eta * P_prev else steffensen(coefs, z0)
        if not np.isfinite(z): return 'overflow', ep
        P_curr = abs(eval_poly_coefs(coefs, z))
        if P_curr < basin_tol: return 'converged', ep
        z0 = z; P_prev = P_curr
    return 'stagnated', max_epochs


# ---- Random polynomial ensembles ------------------------------------------

def sample_KS(d, rng):
    """Kostlan-Smale: a_k ~ N_C(0, C(d,k))."""
    coefs = np.zeros(d + 1, dtype=complex)
    for k in range(d + 1):
        var = math.comb(d, k)
        coefs[k] = (rng.standard_normal() + 1j * rng.standard_normal()) * math.sqrt(var / 2)
    return coefs


def sample_UNI(d, rng):
    """i.i.d. complex standard Gaussian coefficients."""
    return (rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)) / math.sqrt(2)


# ---- Cauchy radius bound --------------------------------------------------

def cauchy_radius(coefs):
    """rho <= 1 + max_k |a_k / a_d|  (classical Cauchy bound)."""
    leading = abs(coefs[-1])
    if leading < 1e-30: return 10.0
    return 1.0 + max(abs(c) / leading for c in coefs[:-1])


# ---- Start strategies -----------------------------------------------------

def starts_A_Cauchy(rho, N, seed=0):
    R = 1.0 + rho
    return [R * np.exp(2j * np.pi * (k + seed * 0.137) / N) for k in range(N)]


def starts_B_GeomDecr(rho, N, q=0.7, seed=0):
    R0 = 1.0 + rho
    golden = math.pi * (3.0 - math.sqrt(5.0))
    return [R0 * (q ** k) * np.exp(1j * (k * golden + seed)) for k in range(N)]


STRATS = {'A:Cauchy': starts_A_Cauchy, 'B:GeomDecr': starts_B_GeomDecr}


# ---- First-success scan ---------------------------------------------------

def first_success(coefs, starts, max_epochs=80):
    for k, z0 in enumerate(starts):
        st, _ = run_T2_K1(coefs, z0, max_epochs=max_epochs)
        if st == 'converged':
            return k + 1
    return None


# ---- Main experiment -------------------------------------------------------

def bucket(s_star, d):
    if s_star is None: return 'fail'
    if s_star <= max(2, math.log(d)): return 'fast'
    if s_star <= math.sqrt(d):       return 'med'
    return 'slow'


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    t0 = time.perf_counter()

    DEGREES = [16, 32, 64, 128]
    N_SAMPLES = 100   # polynomials per (ensemble, d)

    print("=" * 110, flush=True)
    print(f"Empirical worst-case frequency  --  N_SAMPLES = {N_SAMPLES} polynomials per (ensemble,d,strategy)", flush=True)
    print(f"Buckets:  fast  s* <= log d      med  log d < s* <= sqrt(d)      slow  sqrt(d) < s* <= d      fail  s* > d", flush=True)
    print("=" * 110, flush=True)
    print(f"{'ensemble':<6} {'strat':<11} {'d':>4} | "
          f"{'fast %':>8} {'med %':>8} {'slow %':>8} {'fail %':>8} | "
          f"{'mean s*':>8} {'med s*':>7} {'p95 s*':>7} {'p99 s*':>7}", flush=True)
    print("-" * 110, flush=True)

    rng = np.random.default_rng(20260424)
    summary = {}  # (ensemble, strat, d) -> (fast, med, slow, fail, mean, med, p95, p99)

    for ensemble_name, sampler in [('KS', sample_KS), ('UNI', sample_UNI)]:
        for d in DEGREES:
            polys = [sampler(d, rng) for _ in range(N_SAMPLES)]
            for strat_name, gen in STRATS.items():
                buckets = {'fast': 0, 'med': 0, 'slow': 0, 'fail': 0}
                s_list = []
                for coefs in polys:
                    rho = cauchy_radius(coefs)
                    starts = gen(rho, d, seed=0)
                    s = first_success(coefs, starts, max_epochs=60)
                    buckets[bucket(s, d)] += 1
                    if s is not None: s_list.append(s)
                tot = sum(buckets.values())
                pct = {k: 100 * v / tot for k, v in buckets.items()}
                if s_list:
                    mean_s = float(np.mean(s_list))
                    med_s  = int(np.median(s_list))
                    p95_s  = int(np.percentile(s_list, 95))
                    p99_s  = int(np.percentile(s_list, 99))
                else:
                    mean_s = med_s = p95_s = p99_s = -1
                print(f"{ensemble_name:<6} {strat_name:<11} {d:>4} | "
                      f"{pct['fast']:>7.1f}% {pct['med']:>7.1f}% {pct['slow']:>7.1f}% {pct['fail']:>7.1f}% | "
                      f"{mean_s:>8.2f} {med_s:>7} {p95_s:>7} {p99_s:>7}", flush=True)
                summary[(ensemble_name, strat_name, d)] = (pct, mean_s, p95_s)
        print("-" * 110, flush=True)

    # Highlight worst-case fraction
    print("\n" + "=" * 110, flush=True)
    print("WORST-CASE FRACTION = % of polynomials hitting (slow + fail) bucket  (saturates O(d^2 log d))", flush=True)
    print("=" * 110, flush=True)
    print(f"{'ensemble':<6} {'d':>4} | {'A:Cauchy worst%':>18} {'B:GeomDecr worst%':>20} {'reduction factor':>18}", flush=True)
    print("-" * 80, flush=True)
    for ensemble_name in ['KS', 'UNI']:
        for d in DEGREES:
            wA = summary[(ensemble_name, 'A:Cauchy', d)][0]['slow'] + summary[(ensemble_name, 'A:Cauchy', d)][0]['fail']
            wB = summary[(ensemble_name, 'B:GeomDecr', d)][0]['slow'] + summary[(ensemble_name, 'B:GeomDecr', d)][0]['fail']
            ratio = wA / wB if wB > 1e-9 else float('inf') if wA > 0 else 1.0
            print(f"{ensemble_name:<6} {d:>4} | {wA:>17.1f}% {wB:>19.1f}% {ratio:>17.2f}x", flush=True)
        print("-" * 80, flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
