"""
Exploration: K=4 (anchor reset period) and T4 (Shanks-4 acceleration).

Background.
    Part VII fixes K=3 and uses T3 = Aitken Delta^2 on three base-map iterates
    (z, s1, s2) where s_{i+1} = z0 - P(z0)/Q(z0, s_i).  T3 yields one
    extrapolated value at cost ~3 oracle calls per epoch.

    T4 = Shanks transform of order 2 on four iterates (z, s1, s2, s3): a
    cubic-convergence extrapolation that costs ~4 oracle calls per epoch.

    K is the period after which we reset the anchor z0 <- z.  K=3 is the
    setting in Part VII.  K=4 amortises one more step before paying the cost
    of re-anchoring (which itself costs 1 oracle call to evaluate the new P0).

Hypothesis.
    Per-epoch cost: T3 ~ 3 oracle calls; T4 ~ 4 oracle calls.  If T4 cuts the
    number of epochs by more than 4/3, it is a net win.  Same logic for K.

We measure (cost, success) on the same adversarial families as benchmark_no_fallback,
in BSS oracle-call units (one call to P or Q is one unit; we ignore O(1) ops).
"""
from __future__ import annotations
import math
import time
import numpy as np

RNG = np.random.default_rng(20260424)


# ---------------------------------------------------------------------------
# Polynomial primitives — each Q_diff or eval_poly counts as 1 oracle call.
# ---------------------------------------------------------------------------

class OracleCounter:
    def __init__(self):
        self.n_eval = 0
        self.n_q = 0

    def reset(self):
        self.n_eval = 0; self.n_q = 0

    def total(self):
        return self.n_eval + self.n_q


CTR = OracleCounter()


def eval_poly(roots, z):
    CTR.n_eval += 1
    return complex(np.prod(z - roots))


def Q_diff(roots, z0, z):
    CTR.n_q += 1
    if abs(z - z0) < 1e-15:
        dz0 = z0 - roots
        if np.any(np.abs(dz0) < 1e-15):
            return complex(sum(np.prod(np.delete(dz0, k)) for k in range(len(roots))))
        prod_tot = np.prod(dz0)
        return complex(np.sum(prod_tot / dz0))
    # Compute via two evaluations but charge as 1 Q oracle (the divided-difference oracle)
    return (complex(np.prod(z - roots)) - complex(np.prod(z0 - roots))) / (z - z0)


# ---------------------------------------------------------------------------
# Base map and accelerators
# ---------------------------------------------------------------------------

def base_step(roots, z0, z, P0):
    Q = Q_diff(roots, z0, z)
    if abs(Q) < 1e-50:
        return None
    return z0 - P0 / Q


def aitken_delta2(z, s1, s2):
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-50:
        return s1
    return z - (s1 - z) ** 2 / den


def T3_step(roots, z0, z):
    """Part-VII T3: two base steps then Aitken Delta^2 on (z, s1, s2)."""
    P0 = eval_poly(roots, z0)
    s1 = base_step(roots, z0, z, P0)
    if s1 is None: return None
    s2 = base_step(roots, z0, s1, P0)
    if s2 is None: return s1
    return aitken_delta2(z, s1, s2)


def T4_step(roots, z0, z):
    """T4: three base steps then Shanks/Wynn epsilon on (z, s1, s2, s3),
    yielding a cubic-convergence extrapolation. We use the iterated-Aitken
    formula: t1 = A(z, s1, s2), t2 = A(s1, s2, s3), then return
    A(z, t1, t2) — a Wynn epsilon-2 transform.
    """
    P0 = eval_poly(roots, z0)
    s1 = base_step(roots, z0, z, P0)
    if s1 is None: return None
    s2 = base_step(roots, z0, s1, P0)
    if s2 is None: return s1
    s3 = base_step(roots, z0, s2, P0)
    if s3 is None: return aitken_delta2(z, s1, s2)
    t1 = aitken_delta2(z, s1, s2)
    t2 = aitken_delta2(s1, s2, s3)
    return aitken_delta2(z, t1, t2)


def T2_step(roots, z0, z):
    """T2: a single base step (no acceleration). Just the Pandrosion base map."""
    P0 = eval_poly(roots, z0)
    return base_step(roots, z0, z, P0)


# ---------------------------------------------------------------------------
# Steffensen fallback (Part VII)
# ---------------------------------------------------------------------------

def steffensen_fallback(roots, z0, alpha=0.25):
    P0 = eval_poly(roots, z0)
    if abs(P0) < 1e-50: return z0
    z_shift = z0 + P0
    Pshift = eval_poly(roots, z_shift)
    Q_st = (Pshift - P0) / P0
    if abs(Q_st) < 1e-50:
        return z0 - alpha * P0
    return z0 - alpha * P0 / Q_st


# ---------------------------------------------------------------------------
# Generic Pandrosion-Tk with anchor period K and fallback
# ---------------------------------------------------------------------------

def run_Tk_K(roots, z_init, accel_step, K, eta=0.95, max_epochs=120,
             basin_tol=None, alpha=0.25):
    d = len(roots)
    if basin_tol is None:
        basin_tol = (1.0 / max(d, 2)) ** d
    z0 = complex(z_init); z = complex(z_init)
    P_prev = abs(eval_poly(roots, z0))
    if P_prev < basin_tol:
        return {'status': 'converged', 'epochs': 0, 'residual': P_prev}
    for ep in range(1, max_epochs + 1):
        z_cand = z; ok = True
        for _ in range(K):
            zn = accel_step(roots, z0, z_cand)
            if zn is None or not np.isfinite(zn):
                ok = False; break
            z_cand = zn
        if ok and abs(eval_poly(roots, z_cand)) <= eta * P_prev:
            z = z_cand
        else:
            z = steffensen_fallback(roots, z0, alpha)
            if not np.isfinite(z):
                return {'status': 'overflow', 'epochs': ep, 'residual': P_prev}
        P_curr = abs(eval_poly(roots, z))
        if P_curr < basin_tol:
            return {'status': 'converged', 'epochs': ep, 'residual': P_curr}
        z0 = z; P_prev = P_curr
    return {'status': 'stagnated', 'epochs': max_epochs, 'residual': P_prev}


# ---------------------------------------------------------------------------
# Adversarial root families (same as benchmark_no_fallback.py)
# ---------------------------------------------------------------------------

def fam_unity(d):       return np.exp(2j * np.pi * np.arange(d) / d)
def fam_chebyshev(d):   return np.cos((2 * np.arange(d) + 1) * np.pi / (2 * d)).astype(complex)
def fam_clustered(d):   return np.linspace(-1, 1, d).astype(complex)
def fam_mignotte(d):
    eps = 2.0 ** (-d / 3)
    base = np.exp(2j * np.pi * np.arange(d - 1) / (d - 1))
    return np.concatenate([base, [base[0] * (1 + eps)]])[:d]
def fam_multiscale(d):
    s = np.logspace(-2, 0, d)
    return s * np.exp(2j * np.pi * np.arange(d) / d)
def fam_random(d):
    raw = RNG.standard_normal(d) + 1j * RNG.standard_normal(d)
    return raw / np.max(np.abs(raw))

FAMILIES = {
    'unity':      fam_unity,
    'chebyshev':  fam_chebyshev,
    'clustered':  fam_clustered,
    'mignotte':   fam_mignotte,
    'multiscale': fam_multiscale,
    'random':     fam_random,
}


# ---------------------------------------------------------------------------
# Driver: for each (family, d), find min single-start oracle cost across
# 12 starts on the Cauchy circle. Compare T2/T3/T4 x K in {2,3,4,5}.
# ---------------------------------------------------------------------------

ACCEL = {'T2': T2_step, 'T3': T3_step, 'T4': T4_step}


def best_cost(roots, accel_name, K, n_starts=12, max_epochs=80):
    rho = float(np.max(np.abs(roots))); R = 1.0 + rho
    accel = ACCEL[accel_name]
    best = None
    best_eps = None
    for k in range(n_starts):
        theta = 2 * math.pi * k / n_starts
        z0 = R * complex(math.cos(theta), math.sin(theta))
        CTR.reset()
        r = run_Tk_K(roots, z0, accel, K, max_epochs=max_epochs)
        if r['status'] == 'converged':
            cost = CTR.total()
            if best is None or cost < best:
                best = cost; best_eps = r['epochs']
    return best, best_eps


def main():
    t0 = time.perf_counter()
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 110, flush=True)
    print("Sweep over (Tk, K): oracle-call cost on first successful start (smallest of 12 Cauchy starts)", flush=True)
    print("=" * 110, flush=True)
    header = f"{'family':<11} {'d':>4} | "
    pairs = [('T2', 1), ('T2', 3), ('T3', 3), ('T3', 4), ('T3', 5),
             ('T4', 2), ('T4', 3), ('T4', 4)]
    for name, K in pairs:
        header += f"{name+'/K='+str(K):>10} "
    print(header, flush=True)
    print("-" * len(header), flush=True)

    raw = {p: [] for p in pairs}
    for family, builder in FAMILIES.items():
        for d in [8, 16, 32, 64, 128]:
            roots = builder(d)
            row = f"{family:<11} {d:>4} | "
            for name, K in pairs:
                cost, eps = best_cost(roots, name, K,
                                       n_starts=12, max_epochs=80)
                if cost is None:
                    row += f"{'  --  ':>10} "
                else:
                    row += f"{cost:>10} "
                    raw[(name, K)].append((d, cost))
            print(row, flush=True)
        print("-" * len(header), flush=True)

    # Summarise: for each (Tk, K), fit cost(d) = C * d^alpha * (log d)^beta
    print("\n" + "=" * 110, flush=True)
    print("Empirical fits  log(cost) ~ alpha * log d + beta * log log d  (only successful runs)", flush=True)
    print("=" * 110, flush=True)
    print(f"{'config':<10} {'n':>4} {'alpha':>7} {'beta':>7} {'R^2':>7} "
          f"{'mean cost / (d log d)':>22}", flush=True)
    print("-" * 80, flush=True)
    for p in pairs:
        lst = raw[p]
        if len(lst) < 4: continue
        ds = np.array([d for d, _ in lst], dtype=float)
        cs = np.array([c for _, c in lst], dtype=float)
        log_d = np.log(ds); log_log_d = np.log(np.log(ds)); log_c = np.log(cs)
        A = np.stack([log_d, log_log_d, np.ones_like(log_d)], axis=1)
        coef, *_ = np.linalg.lstsq(A, log_c, rcond=None)
        alpha, beta, _ = coef
        pred = A @ coef
        ss_res = float(np.sum((log_c - pred) ** 2))
        ss_tot = float(np.sum((log_c - np.mean(log_c)) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else float('nan')
        ratio = float(np.mean(cs / (ds * np.log(ds))))
        name = f"{p[0]}/K={p[1]}"
        print(f"{name:<10} {len(lst):>4} {alpha:>7.3f} {beta:>7.3f} {r2:>7.3f} "
              f"{ratio:>22.2f}", flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
