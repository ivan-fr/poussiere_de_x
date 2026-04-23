"""
Attack on Conjecture 5.1 (Amortised Block Descent), Paper 1 [pandrosion_smale.pdf].

Conjecture: there exist universal K_0 = O(1), c > 0 such that for every
normalised P in C[z] of degree d with simple roots, any reanchoring period
K <= K_0, and every block of K consecutive adaptive Pandrosion-T3 steps:

    |P(z_K)| <= exp(-c) |P(z_0)|,

equivalently the "sum of logarithms"

    Sigma_log(z0, z) := sum_k log |Q_k(z0, z) / Q(z0, z)| <= -c.

Vectorised numpy implementation for speed: sigma_log is O(d^2) vectorised;
a full orbit at d=500 runs in a few seconds.
"""
from __future__ import annotations
import numpy as np
import time
import math

RNG = np.random.default_rng(20260424)


def eval_poly(roots, z):
    """P(z) = prod_k (z - root_k), monic."""
    return np.prod(z - roots)


def Q_diff(roots, z0, z):
    """Q(z0, z) = (P(z) - P(z0))/(z - z0), with P'(z0) at z = z0."""
    if abs(z - z0) < 1e-18:
        # P'(z0) = sum_k prod_{j != k} (z0 - root_j)
        mat = z0 - roots
        prod_tot = np.prod(mat)
        # divide out each entry (with zero-handling)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(np.where(mat != 0, prod_tot / mat, 0))
    return (eval_poly(roots, z) - eval_poly(roots, z0)) / (z - z0)


def sigma_log_and_Q(roots, z0, z):
    """Return (Sigma_log, Q_value).  Q_k computed via product-ratio trick:
    R_k(w) = P(w)/(w - zeta_k) = prod_{j != k}(w - zeta_j).
    We compute R_k(z) and R_k(z0) as two vectors, Q_k = (R_k(z)-R_k(z0))/(z-z0),
    then sum log|Q_k/Q|.
    """
    d = len(roots)
    dz = z - roots   # length d, dz[k] = z - root_k
    dz0 = z0 - roots # length d, dz0[k] = z0 - root_k
    Pz = np.prod(dz)
    Pz0 = np.prod(dz0)
    # R_k(z) = Pz / dz[k]  (handle zero entries)
    if np.any(np.abs(dz) < 1e-18) or np.any(np.abs(dz0) < 1e-18):
        # rare; fall back to explicit loop
        Rk_z = np.array([np.prod(np.delete(dz, k)) for k in range(d)])
        Rk_z0 = np.array([np.prod(np.delete(dz0, k)) for k in range(d)])
    else:
        Rk_z = Pz / dz
        Rk_z0 = Pz0 / dz0
    if abs(z - z0) < 1e-18:
        # Qk = R_k'(z0) = sum_{m != k} prod_{j != k, m}(z0 - root_j)
        # Vectorise: for each k, sum over m != k of Rk_z0 / dz0[m] ... not
        # Use explicit loop as it's rarely hit.
        Qks = np.zeros(d, dtype=complex)
        for k in range(d):
            for m in range(d):
                if m == k:
                    continue
                t = 1.0 + 0j
                for j in range(d):
                    if j == k or j == m:
                        continue
                    t *= (z0 - roots[j])
                Qks[k] += t
    else:
        Qks = (Rk_z - Rk_z0) / (z - z0)
    Q = (Pz - Pz0) / (z - z0) if abs(z - z0) > 1e-18 else Q_diff(roots, z0, z)
    if abs(Q) < 1e-80:
        return None, Q
    absQks = np.abs(Qks)
    absQ = abs(Q)
    # skip zeros
    mask = absQks > 1e-80
    if not mask.any():
        return None, Q
    logs = np.log(absQks[mask] / absQ)
    return float(np.sum(logs)), Q


def pandrosion_T3_step(roots, z0, z):
    """One T3 (Aitken on base map) from anchor z0, iterate z."""
    P0 = eval_poly(roots, z0)
    Q = Q_diff(roots, z0, z)
    if abs(Q) < 1e-60:
        return None
    s1 = z0 - P0 / Q
    Q2 = Q_diff(roots, z0, s1)
    if abs(Q2) < 1e-60:
        return s1
    s2 = z0 - P0 / Q2
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-60:
        return s1
    return z - (s1 - z) ** 2 / den


# ---------------------------------------------------------------------------
# Adversarial families
# ---------------------------------------------------------------------------

def family_roots_unity(d):
    return np.exp(2j * np.pi * np.arange(d) / d)


def family_random_circle(d):
    return np.exp(2j * np.pi * RNG.random(d))


def family_clustered_line(d):
    return np.linspace(-1, 1, d).astype(complex)


def family_wilkinson(d):
    return (np.arange(1, d + 1) / d).astype(complex)


def family_mignotte(d):
    eps = 2.0 ** (-d / 2)
    base = np.exp(2j * np.pi * np.arange(d - 1) / (d - 1))
    return np.concatenate([base, [base[0] * (1 + eps)]])[:d]


def family_chebyshev(d):
    return np.cos((2 * np.arange(d) + 1) * np.pi / (2 * d)).astype(complex)


def family_random_gaussian(d):
    raw = RNG.standard_normal(d) + 1j * RNG.standard_normal(d)
    return raw / np.max(np.abs(raw))


def family_multiscale(d):
    scales = np.logspace(-3, 0, d)
    return scales * np.exp(2j * np.pi * np.arange(d) / d)


FAMILIES = {
    "roots_unity": family_roots_unity,
    "random_circle": family_random_circle,
    "clustered_line": family_clustered_line,
    "wilkinson": family_wilkinson,
    "mignotte": family_mignotte,
    "chebyshev": family_chebyshev,
    "gaussian": family_random_gaussian,
    "multiscale": family_multiscale,
}


# ---------------------------------------------------------------------------
# Orbit with sigma_log tracking
# ---------------------------------------------------------------------------

def run_orbit(roots, z_init, K, max_epochs):
    """Run adaptive Pandrosion-T3 (period K) for max_epochs epochs.
    Record (Sigma_log_min, log|P_K/P_0|) per epoch.  Stop on convergence."""
    z0 = z_init
    z = z_init
    records = []
    for ep in range(max_epochs):
        P0 = eval_poly(roots, z0)
        if abs(P0) < 1e-35:
            break
        sigmas = []
        ok = True
        for _step in range(K):
            sigma, _Q = sigma_log_and_Q(roots, z0, z)
            if sigma is None:
                ok = False; break
            sigmas.append(sigma)
            z_new = pandrosion_T3_step(roots, z0, z)
            if z_new is None or not np.isfinite(z_new):
                ok = False; break
            z = z_new
        if not ok or not sigmas:
            break
        PK = eval_poly(roots, z)
        if abs(PK) < 1e-35 or abs(P0) < 1e-35:
            records.append((min(sigmas), -float('inf')))
            break
        log_ratio = math.log(abs(PK) / abs(P0))
        records.append((min(sigmas), log_ratio))
        # reanchor
        z0 = z
    return records


def test_config(family_name, d, K, n_starts):
    roots = FAMILIES[family_name](d)
    rho = float(np.max(np.abs(roots)))
    R = 1 + rho
    sigma_mins = []
    log_ratios = []
    max_ep = max(10, int(2 * d * math.log(max(d, 2)) / math.pi))
    max_ep = min(max_ep, 40)
    for s in range(n_starts):
        theta = 2 * np.pi * s / n_starts + np.pi / (2 * max(d, 3))
        z_init = R * np.exp(1j * theta)
        for sigma_min, log_r in run_orbit(roots, z_init, K, max_ep):
            sigma_mins.append(sigma_min)
            if log_r > -float('inf'):
                log_ratios.append(log_r)
    if not sigma_mins:
        return None
    return {
        'n_blocks': len(sigma_mins),
        'sigma_min': min(sigma_mins),
        'sigma_median': float(np.median(sigma_mins)),
        'log_median': float(np.median(log_ratios)) if log_ratios else 0.0,
        'log_max': float(np.max(log_ratios)) if log_ratios else 0.0,
        'violations': sum(1 for s in sigma_mins if s > 1e-10),
    }


def main():
    t0 = time.perf_counter()
    print("=" * 105, flush=True)
    print("Conjecture 5.1 attack: Sigma_log <= -c over adaptive Pandrosion-T3 orbits", flush=True)
    print("=" * 105, flush=True)
    print(f"{'Family':<18} {'d':>5} {'K':>2} {'blocks':>7} {'min Sigma':>11} {'med Sigma':>11} {'med log':>10} {'max log':>10} {'viol':>5}",
          flush=True)
    print("-" * 105, flush=True)

    degrees = [10, 30, 100, 200]
    K_values = [1, 3, 5]
    total_blocks = 0
    total_violations = 0
    all_log_ratios = []

    for family in FAMILIES:
        for d in degrees:
            # skip expensive combinations
            if family == "mignotte" and d > 80:
                continue
            if family == "wilkinson" and d > 150:
                continue
            for K in K_values:
                n_starts = 8 if d <= 100 else 4
                t_start = time.perf_counter()
                r = test_config(family, d, K, n_starts)
                dt = time.perf_counter() - t_start
                if r is None:
                    continue
                total_blocks += r['n_blocks']
                total_violations += r['violations']
                flag = " <<< VIOL" if r['violations'] > 0 else ""
                print(f"{family:<18} {d:>5} {K:>2} {r['n_blocks']:>7} "
                      f"{r['sigma_min']:>11.3f} {r['sigma_median']:>11.3f} "
                      f"{r['log_median']:>10.3f} {r['log_max']:>10.3f} "
                      f"{r['violations']:>5}{flag} [{dt:.1f}s]",
                      flush=True)

    print("-" * 105, flush=True)
    print(f"\nTotal blocks tested: {total_blocks}", flush=True)
    print(f"Violations of Sigma_log <= 0: {total_violations}", flush=True)
    print(f"Elapsed: {time.perf_counter()-t0:.1f}s", flush=True)

    # Aggregate: effective c on one representative configuration
    print("\n" + "=" * 105, flush=True)
    print("Effective universal descent constant c (d = 50, K = 3, across all families)", flush=True)
    print("=" * 105, flush=True)
    all_log_ratios = []
    all_sigmas = []
    for family in FAMILIES:
        if family in ("mignotte", "wilkinson"):
            continue
        roots = FAMILIES[family](50)
        rho = float(np.max(np.abs(roots)))
        R = 1 + rho
        for s in range(10):
            theta = 2 * np.pi * s / 10
            z_init = R * np.exp(1j * theta)
            for sigma_min, log_r in run_orbit(roots, z_init, 3, 20):
                if log_r > -float('inf'):
                    all_log_ratios.append(log_r)
                all_sigmas.append(sigma_min)
    c_eff = -float(np.median(all_log_ratios))
    print(f"  Blocks analysed:         {len(all_log_ratios)}", flush=True)
    print(f"  Median Sigma_log:        {float(np.median(all_sigmas)):.3f}", flush=True)
    print(f"  Median log|P_K/P_0|:     {float(np.median(all_log_ratios)):.3f}", flush=True)
    print(f"  Effective c:             {c_eff:.3f}", flush=True)
    print(f"  Theoretical -pi/2:       {-math.pi/2:.3f}  (exp = {math.exp(-math.pi/2):.3f})", flush=True)
    print(f"  Descent predicted as:    c = pi/2 + O(1/d) = {math.pi/2:.3f} asymptotically", flush=True)


if __name__ == "__main__":
    main()
