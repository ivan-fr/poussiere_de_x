"""
Creative exploration: start-point strategies beyond the equispaced Cauchy circle.

Baseline (Part VII):
    s_k = R * exp(2*pi*i*k/N)   with  R = 1+rho,  rho = max |root|.

Hypothesis: the equispaced Cauchy circle wastes starts on regions where every
root's basin is far away. By varying *both* the angle AND the radius -- a
geometric spiral that probes the interior of the root convex hull -- we should
hit a basin in fewer trials, lowering the empirical exponent alpha that
appeared in the previous regression (T2/K=1 fit  cost ~ 0.34 * d^0.91).

Strategies tested
-----------------
A) Cauchy             : R fixed = 1+rho, equispaced angles (baseline)
B) GeomDecreasing     : R_k = (1+rho) * q^k   (q in (0,1)), equispaced angles
C) LogSpiral          : R(theta) = (1+rho) * exp(-c * theta / (2 pi))
D) HubbardSchleicher  : multi-circle, 3 radii { (1+rho), 0.5(1+rho), 1/(1+rho) }
E) Disk uniform       : i.i.d. uniform in disk of radius 1+rho
F) ConvHullScan       : starts placed at the geometric mean of root pairs
                        (an "inside" probe directly tied to the root constellation,
                         useful only as an upper-bound diagnostic since it peeks
                         at the roots; we use it as oracle of best-possible)

Reported metric
---------------
For each (strategy, family, d): the number  s*  = first start index that
converges within max_epochs, AND the oracle cost of that single orbit.
We average over 5 cyclic shifts of the start sequence to remove offset bias.
"""
from __future__ import annotations
import math, time
import numpy as np

RNG = np.random.default_rng(20260424)


# ---------------------------------------------------------------------------
# Pandrosion T2/K=1 with fallback (winner of previous experiment)
# ---------------------------------------------------------------------------

class C:
    n = 0

def reset_ctr(): C.n = 0


def eval_poly(roots, z):
    return complex(np.prod(z - roots))


def h(roots, z0, P0, w):
    C.n += 1
    if abs(w - z0) < 1e-15:
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


def T2_step(roots, z0, z):
    P0 = eval_poly(roots, z0)
    if abs(P0) < 1e-50: return z0, P0
    s1 = h(roots, z0, P0, z)
    if s1 is None or not np.isfinite(s1): return None, P0
    s2 = h(roots, z0, P0, s1)
    if s2 is None or not np.isfinite(s2): return s1, P0
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-50: return s1, P0
    T2 = z - (s1 - z) ** 2 / den
    if not np.isfinite(T2): return s1, P0
    return T2, P0


def steffensen(roots, z0, alpha=0.25):
    P0 = eval_poly(roots, z0)
    if abs(P0) < 1e-50: return z0
    z_shift = z0 + P0
    C.n += 1
    Pshift = eval_poly(roots, z_shift)
    Q = (Pshift - P0) / P0
    if abs(Q) < 1e-50: return z0 - alpha * P0
    return z0 - alpha * P0 / Q


def run_T2_K1(roots, z_init, max_epochs=80, eta=0.95):
    d = len(roots)
    basin_tol = max((1.0 / (2.0 * max(d, 2))) ** d, 1e-300)
    z0 = complex(z_init); z = complex(z_init)
    P_prev = abs(eval_poly(roots, z0))
    if P_prev < basin_tol: return 'converged', 0, P_prev
    for ep in range(1, max_epochs + 1):
        z_cand, _ = T2_step(roots, z0, z)
        if z_cand is None or not np.isfinite(z_cand):
            z = steffensen(roots, z0)
        else:
            P_cand = abs(eval_poly(roots, z_cand))
            z = z_cand if P_cand <= eta * P_prev else steffensen(roots, z0)
        if not np.isfinite(z): return 'overflow', ep, P_prev
        P_curr = abs(eval_poly(roots, z))
        if P_curr < basin_tol: return 'converged', ep, P_curr
        z0 = z; P_prev = P_curr
    return 'stagnated', max_epochs, P_prev


# ---------------------------------------------------------------------------
# Start-point generators
# ---------------------------------------------------------------------------

def starts_cauchy(rho, N, seed=0):
    R = 1.0 + rho
    return [R * np.exp(2j * np.pi * (k + seed * 0.137) / N) for k in range(N)]


def starts_geom_decreasing(rho, N, q=0.7, seed=0):
    """R_k = (1+rho) * q^k, equispaced angles (golden-angle to avoid resonance)."""
    R0 = 1.0 + rho
    golden = math.pi * (3.0 - math.sqrt(5.0))
    return [R0 * (q ** k) * np.exp(1j * (k * golden + seed)) for k in range(N)]


def starts_log_spiral(rho, N, c=0.3, seed=0):
    """Logarithmic spiral: R(theta) = (1+rho) exp(-c theta / (2 pi))."""
    R0 = 1.0 + rho
    out = []
    for k in range(N):
        theta = 2 * np.pi * k / max(1, int(math.sqrt(N))) + seed
        R = R0 * math.exp(-c * theta / (2 * np.pi))
        if R < 1e-3: R = 1e-3
        out.append(R * np.exp(1j * theta))
    return out


def starts_hubbard_schleicher(rho, N, seed=0):
    """Multi-circle grid: 3 radii x N/3 angles each.
    Inspired by Hubbard-Schleicher-Sutherland's universal initial set for Newton."""
    R = 1.0 + rho
    radii = [R, 0.5 * R, 1.0 / max(R, 1e-3)]
    out = []
    per = max(1, N // len(radii))
    for r in radii:
        for k in range(per):
            theta = 2 * np.pi * (k + seed * 0.137) / per
            out.append(r * np.exp(1j * theta))
    while len(out) < N:
        out.append(R * np.exp(2j * np.pi * len(out) / N))
    return out[:N]


def starts_disk_uniform(rho, N, seed=0):
    rng = np.random.default_rng(20260424 + seed)
    R = 1.0 + rho
    pts = []
    while len(pts) < N:
        x = rng.uniform(-1, 1); y = rng.uniform(-1, 1)
        if x * x + y * y <= 1:
            pts.append(R * (x + 1j * y))
    return pts


def starts_annular_geom(rho, N, q=0.5, seed=0):
    """Concentric circles with geometrically decreasing radius. Each ring has
    N_ring = ceil(R_k / R_0 * sqrt(N)) angles, so the outer ring is densest
    (where most basins of high-degree polynomials live, by Cauchy bound)."""
    R0 = 1.0 + rho
    out = []
    ring = 0
    while len(out) < N:
        Rk = R0 * (q ** ring)
        if Rk < 1e-3: Rk = 1e-3
        n_angles = max(1, int(round(math.sqrt(N) * (q ** ring))))
        for k in range(n_angles):
            if len(out) >= N: break
            theta = 2 * np.pi * (k + seed * 0.137 + 0.5 * ring) / n_angles
            out.append(Rk * np.exp(1j * theta))
        ring += 1
    return out[:N]


STRATEGIES = {
    'A:Cauchy':       starts_cauchy,
    'B:GeomDecr':     starts_geom_decreasing,
    'C:LogSpiral':    starts_log_spiral,
    'D:HubSchle':     starts_hubbard_schleicher,
    'E:DiskUniform':  starts_disk_uniform,
    'F:AnnularGeom':  starts_annular_geom,
}


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
# Driver: scan starts in order, record (s*, oracle_cost_of_winner)
# ---------------------------------------------------------------------------

def first_success(roots, starts_seq, max_epochs=80):
    """Return (s_index, oracle_cost) of the FIRST start that converges, or
    (None, None) if all fail."""
    for k, z0 in enumerate(starts_seq):
        reset_ctr()
        st, ep, res = run_T2_K1(roots, z0, max_epochs)
        if st == 'converged':
            return k + 1, C.n      # +1 = 1-indexed for human reading
    return None, None


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    t0 = time.perf_counter()

    DEGREES = [16, 32, 64, 128, 256]
    N_BUDGET = lambda d: max(8, d)        # at most d starts (matches Part VII)
    N_SHIFTS = 5                           # average over 5 cyclic shifts

    print("=" * 130, flush=True)
    print("Start-point strategies for Pandrosion-T2/K=1 -- median (s*, cost) over 5 shifts of the start sequence", flush=True)
    print("=" * 130, flush=True)
    head = f"{'family':<11} {'d':>4} | "
    for name in STRATEGIES: head += f"{name+' s*/cost':>20}"
    print(head, flush=True)
    print("-" * len(head), flush=True)

    # raw[strategy] -> list of (d, cost_winner)
    raw_cost = {name: [] for name in STRATEGIES}
    raw_sstar = {name: [] for name in STRATEGIES}

    for family, builder in FAMILIES.items():
        for d in DEGREES:
            roots = builder(d)
            rho = float(np.max(np.abs(roots)))
            row = f"{family:<11} {d:>4} | "
            for name, gen in STRATEGIES.items():
                s_list = []; c_list = []
                for shift in range(N_SHIFTS):
                    starts = gen(rho, N_BUDGET(d), seed=shift)
                    s, c = first_success(roots, starts, max_epochs=80)
                    if s is not None:
                        s_list.append(s); c_list.append(c)
                if s_list:
                    sm = int(np.median(s_list)); cm = int(np.median(c_list))
                    row += f" {sm:>4}/{cm:<14}"
                    raw_cost[name].append((d, cm))
                    raw_sstar[name].append((d, sm))
                else:
                    row += f" {'--':>4}/{'--':<14}"
            print(row, flush=True)
        print("-" * len(head), flush=True)

    # Regressions
    print("\n" + "=" * 130, flush=True)
    print("Empirical fits across (family,d): cost ~ d^alpha (log d)^beta  ;  s* ~ d^gamma", flush=True)
    print("=" * 130, flush=True)
    print(f"{'strategy':<14} {'n':>4} | {'cost alpha':>10} {'cost beta':>10} {'cost R^2':>9} "
          f"{'mean cost / (d log d)':>22} | {'s* gamma':>9} {'s* mean@d=128':>14}", flush=True)
    print("-" * 110, flush=True)
    for name in STRATEGIES:
        lst = raw_cost[name]; ls = raw_sstar[name]
        if len(lst) < 6: continue
        ds = np.array([d for d, _ in lst], dtype=float)
        cs = np.array([c for _, c in lst], dtype=float)
        log_d = np.log(ds); log_log_d = np.log(np.log(ds)); log_c = np.log(cs)
        A = np.stack([log_d, log_log_d, np.ones_like(log_d)], axis=1)
        coef, *_ = np.linalg.lstsq(A, log_c, rcond=None)
        alpha, beta, _ = coef
        pred = A @ coef
        ss_res = float(np.sum((log_c - pred) ** 2)); ss_tot = float(np.sum((log_c - np.mean(log_c)) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else float('nan')
        ratio = float(np.mean(cs / (ds * np.log(ds))))

        # s* ~ d^gamma
        ss_arr = np.array([s for _, s in ls], dtype=float)
        ds2 = np.array([d for d, _ in ls], dtype=float)
        if len(ss_arr) > 3:
            B = np.stack([np.log(ds2), np.ones_like(ds2)], axis=1)
            coefs, *_ = np.linalg.lstsq(B, np.log(np.maximum(ss_arr, 1.0)), rcond=None)
            gamma = coefs[0]
        else:
            gamma = float('nan')
        s128 = [s for d, s in ls if d == 128]
        s128m = float(np.mean(s128)) if s128 else float('nan')
        print(f"{name:<14} {len(lst):>4} | {alpha:>10.3f} {beta:>10.3f} {r2:>9.3f} "
              f"{ratio:>22.3f} | {gamma:>9.3f} {s128m:>14.2f}", flush=True)

    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
