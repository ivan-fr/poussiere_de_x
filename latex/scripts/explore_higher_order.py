"""
Creative exploration of the Pandrosion Tn hierarchy (Def. 4.2 of Paper 1) toward
the most aggressive complexity achievable.

Per Paper 1 §4.4 + §5.3:
    * Each Tn uses n base-map evaluations h(.) and reaches convergence order n
      via Richardson extrapolation (NOT iterated Aitken!):
          T2(z) = z - (s1-z)^2 / (s2 - 2 s1 + z),    lambda_hat = (s2-s1)/(s1-z)
          T3(z) = T2(z) - (s3 - T2(z)) / (lambda_hat - 1)
          T4(z) = T3(z) - (s4 - T3(z)) / (lambda_hat - 1)
          ...
      where s1 = h(z), s_{k+1} = h(T_k(z)) for k >= 2 and s2 = h(s1).
    * After lock-in to a quadratic basin, residual decay is O(log log eps^-1).

Theoretical optimum.  In oracle units, Tn costs n calls per epoch and gains
log_n(1/eps) bits per epoch after lock-in. The ratio n/log n is minimised at
n = e ~ 2.7  =>  T3 is theoretically optimal (factor 2.73 versus 2.88 for T2).

Three creative strategies tested below:
    (S1) Static high-order Tn for n in {2, 3, 4, 6, 8, 12}.
    (S2) "Doubling"  -- gear-shift from T2 (cheap, robust, pre-lock-in) to
         T_{ceil(log d)} (fast, post-lock-in) once |P(z)| < threshold.
    (S3) Dynamic-n: at each epoch, choose n that minimises the predicted cost
         given the current empirical contraction ratio  lambda_hat.

Target.  The user requests O(log(d log d)) per-orbit cost.  In BSS, every oracle
call costs O(d) -- so a per-orbit count of O(log(d log d)) BASE-MAP calls would
yield O(d log d log log d) total ops, sub-Smale-17.  We test if this is reachable.
"""
from __future__ import annotations
import math, time
import numpy as np

RNG = np.random.default_rng(20260424)


# ---------------------------------------------------------------------------
# Oracle counter (BSS).  Each h evaluation = 1 oracle.
# ---------------------------------------------------------------------------
class Counter:
    def __init__(self): self.n = 0
    def reset(self): self.n = 0

CTR = Counter()


def eval_poly(roots, z):
    return complex(np.prod(z - roots))


def h_step(roots, z0, P0, w):
    """Generalized Pandrosion base map evaluated at w from anchor z0:
       h(w) = z0 - P(z0)/Q(z0, w)   with  Q(z0, w) = (P(w) - P(z0))/(w - z0).
    Costs 1 oracle call (one P(w) evaluation; P(z0) is reused)."""
    CTR.n += 1
    if abs(w - z0) < 1e-15:
        # P'(z0) limit
        dz0 = z0 - roots
        if np.any(np.abs(dz0) < 1e-15):
            Q = complex(sum(np.prod(np.delete(dz0, k)) for k in range(len(roots))))
        else:
            Q = complex(np.sum(np.prod(dz0) / dz0))
    else:
        Q = (eval_poly(roots, w) - P0) / (w - z0)
    if abs(Q) < 1e-50:
        return None
    return z0 - P0 / Q


# ---------------------------------------------------------------------------
# Tn hierarchy (Definition 4.2, Paper 1) -- proper Richardson extrapolation
# ---------------------------------------------------------------------------

def T_n_step(roots, z0, z, n):
    """Compute T_n(z) using n evaluations of h.

    Returns (T_n_value, lambda_hat) or (None, None) on failure.
    """
    P0 = eval_poly(roots, z0)
    if abs(P0) < 1e-50:
        return z0, 1.0

    s1 = h_step(roots, z0, P0, z)
    if s1 is None or not np.isfinite(s1):
        return None, None
    if n == 1:
        return s1, None

    s2 = h_step(roots, z0, P0, s1)
    if s2 is None or not np.isfinite(s2):
        return s1, None

    den = s2 - 2 * s1 + z
    if abs(den) < 1e-50:
        return s1, None
    T2 = z - (s1 - z) ** 2 / den
    if not np.isfinite(T2):
        return s1, None
    if n == 2:
        return T2, (s2 - s1) / (s1 - z) if abs(s1 - z) > 1e-50 else None

    if abs(s1 - z) < 1e-50:
        return T2, None
    lam = (s2 - s1) / (s1 - z)

    Tk = T2
    for k in range(3, n + 1):
        sk = h_step(roots, z0, P0, Tk)
        if sk is None or not np.isfinite(sk):
            return Tk, lam
        if abs(lam - 1) < 1e-50:
            return Tk, lam
        Tk_new = Tk - (sk - Tk) / (lam - 1)
        if not np.isfinite(Tk_new):
            return Tk, lam
        Tk = Tk_new
    return Tk, lam


# ---------------------------------------------------------------------------
# Steffensen fallback (Part VII) -- 2 oracle calls
# ---------------------------------------------------------------------------

def steffensen_fallback(roots, z0, alpha=0.25):
    P0 = eval_poly(roots, z0)
    if abs(P0) < 1e-50:
        return z0
    z_shift = z0 + P0
    CTR.n += 1
    Pshift = eval_poly(roots, z_shift)
    Q = (Pshift - P0) / P0
    if abs(Q) < 1e-50:
        return z0 - alpha * P0
    return z0 - alpha * P0 / Q


# ---------------------------------------------------------------------------
# Driver: adaptive-anchor Pandrosion-Tn with K=1 (Newton-like anchor) +
# Steffensen fallback.  Returns (status, epochs, residual, oracle_calls).
# ---------------------------------------------------------------------------

def run_orbit(roots, z_init, accel_callable, max_epochs, eta=0.95,
              basin_tol=None):
    """accel_callable(roots, z0, z) -> (z_new, lam) or (None, None)."""
    d = len(roots)
    if basin_tol is None:
        # Smale-17 basin: |P(z)|^(1/d) <= 1/(2d) i.e. |P(z)| <= (2d)^(-d).
        # Cap below the float underflow threshold so high d remains testable.
        basin_tol = max((1.0 / (2.0 * max(d, 2))) ** d, 1e-300)
    z0 = complex(z_init)
    z = complex(z_init)
    P_prev = abs(eval_poly(roots, z0))
    if P_prev < basin_tol:
        return 'converged', 0, P_prev
    for ep in range(1, max_epochs + 1):
        z_cand, lam = accel_callable(roots, z0, z)
        if z_cand is None or not np.isfinite(z_cand):
            z = steffensen_fallback(roots, z0)
        else:
            P_cand = abs(eval_poly(roots, z_cand))
            if P_cand <= eta * P_prev:
                z = z_cand
            else:
                z = steffensen_fallback(roots, z0)
        if not np.isfinite(z):
            return 'overflow', ep, P_prev
        P_curr = abs(eval_poly(roots, z))
        if P_curr < basin_tol:
            return 'converged', ep, P_curr
        z0 = z; P_prev = P_curr
    return 'stagnated', max_epochs, P_prev


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

def make_static_Tn(n):
    return lambda roots, z0, z: T_n_step(roots, z0, z, n)


def make_doubling(d):
    """S2: T2 until pre-lock-in, then T_{ceil(log2 d)} after.
    We approximate lock-in detection by tracking |P(z)|.
    Since the driver doesn't expose P(z) here, we use a closure with a counter
    that ramps n linearly with epoch index."""
    n_late = max(2, int(math.ceil(math.log2(max(d, 2)))))
    state = {'ep': 0}

    def step(roots, z0, z):
        state['ep'] += 1
        n = 2 if state['ep'] <= 4 else n_late
        return T_n_step(roots, z0, z, n)

    return step


def make_dynamic(d):
    """S3: pick n* = argmin_n (n / max(log(1/|lam|), 1e-3)) -- the cost of one
    additional Tn step relative to the gain in correct digits."""
    state = {'lam': 0.5, 'ep': 0}
    n_max = max(3, int(math.ceil(math.log2(max(d, 2)))))

    def step(roots, z0, z):
        state['ep'] += 1
        # Translate empirical contraction lam into an order:
        # decay per Tn step ~ |lam|^n, so we want n where n / log(1/|lam|) is min,
        # but n must be >= 2 for stability.  Use a quick rule: n = ceil(2 - log|lam|).
        absl = max(min(abs(state['lam']), 0.99), 1e-3)
        n_pick = max(2, min(n_max, int(round(2 - math.log(absl)))))
        z_new, lam = T_n_step(roots, z0, z, n_pick)
        if lam is not None and abs(lam) > 1e-6:
            state['lam'] = complex(lam)
        return z_new, lam

    return step


# ---------------------------------------------------------------------------
# Adversarial families
# ---------------------------------------------------------------------------

def fam_unity(d):       return np.exp(2j * np.pi * np.arange(d) / d)
def fam_chebyshev(d):   return np.cos((2 * np.arange(d) + 1) * np.pi / (2 * d)).astype(complex)
def fam_clustered(d):   return np.linspace(-1, 1, d).astype(complex)
def fam_wilkinson(d):   return (np.arange(1, d + 1) / d).astype(complex)
def fam_mignotte(d):
    eps = 2.0 ** (-d / 3)
    base = np.exp(2j * np.pi * np.arange(d - 1) / (d - 1))
    return np.concatenate([base, [base[0] * (1 + eps)]])[:d]
def fam_random(d):
    raw = RNG.standard_normal(d) + 1j * RNG.standard_normal(d)
    return raw / np.max(np.abs(raw))

FAMILIES = {
    'unity':      fam_unity,
    'chebyshev':  fam_chebyshev,
    'clustered':  fam_clustered,
    'wilkinson':  fam_wilkinson,
    'mignotte':   fam_mignotte,
    'random':     fam_random,
}


# ---------------------------------------------------------------------------
# Per-config experiment
# ---------------------------------------------------------------------------

def best_oracle_cost(roots, accel_factory, n_starts=12, max_epochs=80):
    """Try `n_starts` Cauchy-circle starts, return (oracle_calls, epochs) of
    the cheapest successful run. accel_factory is either a fixed callable or a
    factory of the form lambda d: callable (rebuilt for each start so dynamic
    state resets)."""
    rho = float(np.max(np.abs(roots))); R = 1.0 + rho
    best_cost = None; best_eps = None
    for k in range(n_starts):
        theta = 2 * math.pi * k / n_starts
        z0 = R * complex(math.cos(theta), math.sin(theta))
        if callable(accel_factory) and not isinstance(accel_factory, type(make_static_Tn(2))):
            accel = accel_factory  # static
        else:
            accel = accel_factory  # already a callable
        # Rebuild fresh dynamic state if the accel was returned from a factory:
        if hasattr(accel, '__name__') is False and callable(accel):
            pass  # plain callables (T_n closures) have no per-orbit state to reset
        CTR.reset()
        status, ep, res = run_orbit(roots, z0, accel, max_epochs)
        if status == 'converged':
            if best_cost is None or CTR.n < best_cost:
                best_cost = CTR.n; best_eps = ep
    return best_cost, best_eps


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    t0 = time.perf_counter()

    print("=" * 110, flush=True)
    print("S1 -- Static Tn (Richardson, Def. 4.2 of Paper 1)", flush=True)
    print("=" * 110, flush=True)
    static_orders = [2, 3, 4, 6, 8, 12]
    header = f"{'family':<11} {'d':>4} | "
    for n in static_orders: header += f"{'T'+str(n):>8} "
    header += "  S2:double   S3:dyn"
    print(header, flush=True)
    print("-" * len(header), flush=True)

    raw = {n: [] for n in static_orders}
    raw['double'] = []
    raw['dyn'] = []

    DEGREES = [8, 16, 32, 64, 128, 256, 512]
    for family, builder in FAMILIES.items():
        for d in DEGREES:
            roots = builder(d)
            row = f"{family:<11} {d:>4} | "
            for n in static_orders:
                accel = make_static_Tn(n)
                c, e = best_oracle_cost(roots, accel, n_starts=12, max_epochs=120)
                if c is None:
                    row += f"{'  --  ':>8} "
                else:
                    row += f"{c:>8} "
                    raw[n].append((d, c))

            # S2 doubling -- factory rebuilds per start
            best_c2 = None
            for k in range(12):
                theta = 2 * math.pi * k / 12
                rho = float(np.max(np.abs(roots))); R = 1.0 + rho
                z0 = R * complex(math.cos(theta), math.sin(theta))
                CTR.reset()
                accel = make_doubling(d)
                st, ep, res = run_orbit(roots, z0, accel, 120)
                if st == 'converged':
                    if best_c2 is None or CTR.n < best_c2:
                        best_c2 = CTR.n
            row += (f" {best_c2:>9} " if best_c2 else f" {'--':>9} ")
            if best_c2 is not None:
                raw['double'].append((d, best_c2))

            # S3 dynamic
            best_c3 = None
            for k in range(12):
                theta = 2 * math.pi * k / 12
                rho = float(np.max(np.abs(roots))); R = 1.0 + rho
                z0 = R * complex(math.cos(theta), math.sin(theta))
                CTR.reset()
                accel = make_dynamic(d)
                st, ep, res = run_orbit(roots, z0, accel, 120)
                if st == 'converged':
                    if best_c3 is None or CTR.n < best_c3:
                        best_c3 = CTR.n
            row += (f" {best_c3:>7}" if best_c3 else f" {'--':>7}")
            if best_c3 is not None:
                raw['dyn'].append((d, best_c3))

            print(row, flush=True)
        print("-" * len(header), flush=True)

    # Summarise
    print("\n" + "=" * 110, flush=True)
    print("Empirical fit  cost ~ C * d^alpha * (log d)^beta  -- want alpha=1, beta close to 0", flush=True)
    print("=" * 110, flush=True)
    print(f"{'config':<10} {'n':>4} {'alpha':>7} {'beta':>7} {'R^2':>7} "
          f"{'mean cost / (d log d)':>22}  per_orbit_at_d=256", flush=True)
    print("-" * 92, flush=True)
    for key in static_orders + ['double', 'dyn']:
        lst = raw.get(key, [])
        if len(lst) < 6:
            continue
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
        d256 = [c for d, c in lst if d == 256]
        d256_mean = np.mean(d256) if d256 else float('nan')
        name = f"T{key}" if isinstance(key, int) else key
        print(f"{name:<10} {len(lst):>4} {alpha:>7.3f} {beta:>7.3f} {r2:>7.3f} "
              f"{ratio:>22.3f}  {d256_mean:>10.0f}", flush=True)

    # Pure base-map count after lock-in: each Tn step gives log_n(d^d) digits
    # so post-lock-in epochs ~ d log d / (n log n). Total oracle ~ d log d / log n
    # which is minimal at n -> infinity in theory but pre-lock-in dominates.
    print("\n" + "=" * 110, flush=True)
    print("Theoretical lower bound (oracle units) per orbit assuming clean lock-in", flush=True)
    print("=" * 110, flush=True)
    print(f"{'n':>3} {'n / log n':>10} {'epochs to reach |P|=d^-d':>24} {'oracle cost':>14}", flush=True)
    print("-" * 60, flush=True)
    for n in [2, 3, 4, 6, 8, 12, 24]:
        if n < 2: continue
        ratio = n / math.log(n)
        # post-lock-in: convergence order n => log_n(target/eps); target = d, here eps = d^-d
        # epochs ~ log(d log d) / log n, total oracle ~ n * epochs.
        for d in [128]:
            ep = math.log(d * math.log(d)) / math.log(n)
            tot = n * ep
            print(f"{n:>3} {ratio:>10.2f} {ep:>24.2f} {tot:>14.2f}  (d={d})", flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
