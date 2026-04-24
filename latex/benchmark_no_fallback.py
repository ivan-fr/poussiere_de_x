"""
Benchmark: Pandrosion-T3 WITHOUT fallback + search for complexity better than O(d^2 log d).

Background (Universitas Pandrosion VII):
    The article proves that Pandrosion-T3 *with* a non-holomorphic fallback achieves
    O(d^2 log d) BSS operations on any normalised polynomial of degree d.
    Decomposition:    d  starts  x  O(d)  per oracle  x  O(log d)  epochs.
    Without the fallback the algorithm hits McMullen's topological barrier and is
    expected to stagnate on a positive-measure set of root constellations for d >= 4.

This script measures, on adversarial root families (Mignotte clusters, roots of
unity, Chebyshev nodes, multiscale, etc.):

    (A) Failure rate of Pandrosion-T3 *without* fallback  ->  empirical witness
        of the McMullen barrier.

    (B) The empirical scaling of #epochs(d) and #successful_starts(d), so that
        we can fit  cost(d) ~ d^alpha (log d)^beta and ask whether (alpha, beta)
        can be pushed below (2, 1).

The complexity of one epoch from a single start is O(d) BSS operations (one P
oracle call dominates). With s concurrent starts and E epochs to reach the basin,
the total cost is  s * E * O(d).  Part VII chooses s = d, E = O(log d), giving
O(d^2 log d). We empirically scan s and E.
"""
from __future__ import annotations
import math
import time
import numpy as np

RNG = np.random.default_rng(20260424)

# ---------------------------------------------------------------------------
# Polynomial primitives in the BSS-style "oracle = P(z)" model.
# Each call to eval_poly counts as O(d) BSS ops.
# ---------------------------------------------------------------------------

def eval_poly(roots: np.ndarray, z: complex) -> complex:
    return complex(np.prod(z - roots))


def Q_diff(roots: np.ndarray, z0: complex, z: complex) -> complex:
    """Pandrosion divided difference Q(z0, z) = (P(z) - P(z0))/(z - z0)."""
    if abs(z - z0) < 1e-15:
        # P'(z0) via product-ratio
        dz0 = z0 - roots
        if np.any(np.abs(dz0) < 1e-15):
            return complex(sum(np.prod(np.delete(dz0, k)) for k in range(len(roots))))
        prod_tot = np.prod(dz0)
        return complex(np.sum(prod_tot / dz0))
    return (eval_poly(roots, z) - eval_poly(roots, z0)) / (z - z0)


def pandrosion_T3_step(roots, z0, z):
    """One T3-accelerated Pandrosion step from anchor z0, current iterate z.
    T3 = Aitken acceleration on two consecutive base-map iterates."""
    P0 = eval_poly(roots, z0)
    Q1 = Q_diff(roots, z0, z)
    if abs(Q1) < 1e-50:
        return None
    s1 = z0 - P0 / Q1
    Q2 = Q_diff(roots, z0, s1)
    if abs(Q2) < 1e-50:
        return s1
    s2 = z0 - P0 / Q2
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-50:
        return s1
    return z - (s1 - z) ** 2 / den


# ---------------------------------------------------------------------------
# Pandrosion-T3 WITHOUT fallback: pure rational iteration.
# We DO NOT execute the Steffensen escape step, so McMullen's theorem (1987)
# predicts there must exist root constellations + anchors yielding stagnation.
# ---------------------------------------------------------------------------

def run_pandrosion_no_fallback(roots, z_init, K=3, eta=0.95,
                               max_epochs=200, basin_tol=None):
    """Run pure Pandrosion-T3 until basin entry or stagnation.

    Convention: epoch = K base-map steps, then anchor reset z0 <- z.
    A "stagnation trap" is detected when |P(z_K)| > eta * |P(z_0)|, exactly the
    test that the fallback in Part VII would react to. WITHOUT the fallback we
    simply give up on this start.

    Returns a dict with 'status' in {'converged', 'stagnated', 'overflow'},
    'epochs', and the residual modulus reached.
    """
    d = len(roots)
    if basin_tol is None:
        basin_tol = (1.0 / max(d, 2)) ** d  # |P(z)|^(1/d) <= 1/d
    z0 = complex(z_init)
    z = complex(z_init)
    P_prev = abs(eval_poly(roots, z0))
    if P_prev < basin_tol:
        return {'status': 'converged', 'epochs': 0, 'residual': P_prev}

    for ep in range(1, max_epochs + 1):
        ok = True
        for _ in range(K):
            z_new = pandrosion_T3_step(roots, z0, z)
            if z_new is None or not np.isfinite(z_new):
                ok = False
                break
            z = z_new
        if not ok:
            return {'status': 'overflow', 'epochs': ep, 'residual': P_prev}

        P_curr = abs(eval_poly(roots, z))
        if P_curr < basin_tol:
            return {'status': 'converged', 'epochs': ep, 'residual': P_curr}
        # McMullen-style stagnation test (no fallback to repair it)
        if P_curr > eta * P_prev:
            return {'status': 'stagnated', 'epochs': ep, 'residual': P_curr}

        z0 = z          # anchor reset (T3 with adaptive anchor)
        P_prev = P_curr

    return {'status': 'stagnated', 'epochs': max_epochs, 'residual': P_prev}


# ---------------------------------------------------------------------------
# Multistart driver: s concurrent starts on the Cauchy circle.
# Returns the FIRST start that reaches the basin (epochs spent on that start)
# and the total number of starts attempted. This is the natural way to
# parallelise: in BSS we count cost per orbit but report the minimum.
# ---------------------------------------------------------------------------

def multistart(roots, n_starts, K=3, eta=0.95, max_epochs=200, basin_tol=None,
               with_fallback=False, alpha=0.25):
    rho = float(np.max(np.abs(roots)))
    R = 1.0 + rho
    successes = []
    for k in range(n_starts):
        theta = 2 * math.pi * k / n_starts + math.pi / (2 * max(n_starts, 3))
        z_init = R * complex(math.cos(theta), math.sin(theta))
        if with_fallback:
            res = run_pandrosion_with_fallback(roots, z_init, K, eta,
                                               max_epochs, basin_tol, alpha)
        else:
            res = run_pandrosion_no_fallback(roots, z_init, K, eta,
                                             max_epochs, basin_tol)
        if res['status'] == 'converged':
            successes.append((k, res['epochs'], res['residual']))
    return successes


# ---------------------------------------------------------------------------
# Reference: Pandrosion-T3 WITH fallback (Definition 2.1, Part VII).
# We add the derivative-free Steffensen step when |P(z_cand)| > eta |P(z0)|.
# ---------------------------------------------------------------------------

def steffensen_fallback_step(roots, z0, alpha=0.25):
    P0 = eval_poly(roots, z0)
    if abs(P0) < 1e-50:
        return z0
    z_shift = z0 + P0
    Pshift = eval_poly(roots, z_shift)
    Q_steff = (Pshift - P0) / P0  # since z_shift - z0 = P0
    if abs(Q_steff) < 1e-50:
        return z0 - alpha * P0  # safe drift
    return z0 - alpha * P0 / Q_steff


def run_pandrosion_with_fallback(roots, z_init, K=3, eta=0.95,
                                 max_epochs=200, basin_tol=None, alpha=0.25):
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
            zn = pandrosion_T3_step(roots, z0, z_cand)
            if zn is None or not np.isfinite(zn):
                ok = False; break
            z_cand = zn
        if ok and abs(eval_poly(roots, z_cand)) <= eta * P_prev:
            z = z_cand
        else:
            z = steffensen_fallback_step(roots, z0, alpha)
            if not np.isfinite(z):
                return {'status': 'overflow', 'epochs': ep, 'residual': P_prev}
        P_curr = abs(eval_poly(roots, z))
        if P_curr < basin_tol:
            return {'status': 'converged', 'epochs': ep, 'residual': P_curr}
        z0 = z; P_prev = P_curr
    return {'status': 'stagnated', 'epochs': max_epochs, 'residual': P_prev}


# ---------------------------------------------------------------------------
# Adversarial root families
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
def fam_anti(d):
    """Adversarial McMullen-style: roots placed so the reverse polynomial has
    a saddle of |P| inside the convex hull aligned with d equispaced anchors."""
    base = np.exp(2j * np.pi * np.arange(d) / d)
    return base * (1 + 0.05 * np.exp(1j * np.pi * np.arange(d) / d))

FAMILIES = {
    'unity':       fam_unity,
    'chebyshev':   fam_chebyshev,
    'clustered':   fam_clustered,
    'mignotte':    fam_mignotte,
    'multiscale':  fam_multiscale,
    'random':      fam_random,
    'anti_mcm':    fam_anti,
}


# ---------------------------------------------------------------------------
# Experiment A: stagnation rate WITHOUT fallback
# ---------------------------------------------------------------------------

def experiment_no_fallback():
    print("=" * 95, flush=True)
    print("EXPERIMENT A — Pandrosion-T3 WITHOUT fallback : McMullen barrier in action", flush=True)
    print("=" * 95, flush=True)
    print(f"{'family':<12} {'d':>4} {'starts':>6} {'conv':>5} {'stag':>5} {'over':>5} "
          f"{'minE':>5} {'medE':>6}  worst_residual", flush=True)
    print("-" * 95, flush=True)
    rows = []
    for family, builder in FAMILIES.items():
        for d in [4, 8, 16, 32, 64, 128]:
            roots = builder(d)
            rho = float(np.max(np.abs(roots)))
            R = 1.0 + rho
            n_starts = max(8, d)         # Part VII uses d concurrent starts
            n_conv = n_stag = n_over = 0
            epochs_ok = []
            worst_res = 0.0
            for k in range(n_starts):
                theta = 2 * math.pi * k / n_starts
                z0 = R * complex(math.cos(theta), math.sin(theta))
                r = run_pandrosion_no_fallback(roots, z0, K=3,
                                               max_epochs=4 * max(int(math.log2(d + 1)), 2))
                if r['status'] == 'converged':
                    n_conv += 1; epochs_ok.append(r['epochs'])
                elif r['status'] == 'stagnated':
                    n_stag += 1; worst_res = max(worst_res, r['residual'])
                else:
                    n_over += 1
            minE = min(epochs_ok) if epochs_ok else -1
            medE = float(np.median(epochs_ok)) if epochs_ok else -1
            print(f"{family:<12} {d:>4} {n_starts:>6} {n_conv:>5} {n_stag:>5} "
                  f"{n_over:>5} {minE:>5} {medE:>6}  {worst_res:.3e}", flush=True)
            rows.append((family, d, n_starts, n_conv, n_stag, n_over))
    return rows


# ---------------------------------------------------------------------------
# Experiment B: empirical scaling of (#starts needed, #epochs needed)
# When at least ONE of s starts converges in E epochs, total cost is s*E*O(d)
# We minimise (s, E) for each (family, d) and fit cost(d).
# ---------------------------------------------------------------------------

def first_success_epochs(roots, n_starts, K=3, max_epochs=80, with_fallback=True):
    rho = float(np.max(np.abs(roots))); R = 1.0 + rho
    runner = (run_pandrosion_with_fallback if with_fallback else
              run_pandrosion_no_fallback)
    best_E = None
    success_idx = None
    for k in range(n_starts):
        theta = 2 * math.pi * k / n_starts
        z0 = R * complex(math.cos(theta), math.sin(theta))
        r = runner(roots, z0, K=K, max_epochs=max_epochs)
        if r['status'] == 'converged':
            if best_E is None or r['epochs'] < best_E:
                best_E = r['epochs']; success_idx = k
    return best_E, success_idx


def experiment_scaling():
    print("\n" + "=" * 95, flush=True)
    print("EXPERIMENT B — Empirical cost scan (with fallback to guarantee descent)", flush=True)
    print("=" * 95, flush=True)
    print(f"{'family':<12} {'d':>4} | {'min starts':>10} {'min epochs':>10} | "
          f"{'cost / d':>10} {'/(d log d)':>10} {'/(d^2 log d)':>13}", flush=True)
    print("-" * 95, flush=True)

    rows = []
    for family, builder in FAMILIES.items():
        for d in [4, 8, 16, 32, 64, 128]:
            roots = builder(d)
            # Search the smallest (s, E) such that some start succeeds in <= E epochs.
            # Scan starts in geometric ladder and epochs implicitly via the run.
            min_starts_needed = None
            best_E_at_min = None
            for s in [1, 2, 4, max(4, int(math.log2(d + 1))), 2 * int(math.log2(d + 1)),
                      int(math.sqrt(d)), d]:
                E, _ = first_success_epochs(roots, s, K=3,
                                            max_epochs=4 * max(int(math.log2(d + 1)), 2))
                if E is not None:
                    min_starts_needed = s
                    best_E_at_min = E
                    break
            if min_starts_needed is None:
                # all probed s failed within epoch budget; force s=d with a generous budget
                E, _ = first_success_epochs(roots, d, K=3,
                                            max_epochs=8 * max(int(math.log2(d + 1)), 2))
                min_starts_needed = d
                best_E_at_min = E if E is not None else -1

            # Per-orbit cost = O(d) per epoch. Total cost = s * E * d (BSS units).
            if best_E_at_min and best_E_at_min > 0:
                cost = min_starts_needed * best_E_at_min * d
                norm_d = cost / d if d else float('nan')
                logd = math.log(max(d, 2))
                norm_dlogd = cost / (d * logd)
                norm_d2logd = cost / (d * d * logd)
            else:
                cost = norm_d = norm_dlogd = norm_d2logd = float('nan')
            print(f"{family:<12} {d:>4} | {min_starts_needed:>10} {best_E_at_min:>10} | "
                  f"{norm_d:>10.2f} {norm_dlogd:>10.3f} {norm_d2logd:>13.4f}", flush=True)
            rows.append((family, d, min_starts_needed, best_E_at_min, cost))
    return rows


# ---------------------------------------------------------------------------
# Experiment C: regression  log(cost) = alpha * log(d) + beta * log(log d) + C
# across degrees, using the WITH-FALLBACK results from Experiment B.
# ---------------------------------------------------------------------------

def experiment_regression(rows_b):
    print("\n" + "=" * 95, flush=True)
    print("EXPERIMENT C — Regression  log cost  ~  alpha * log d  +  beta * log log d", flush=True)
    print("=" * 95, flush=True)
    print(f"{'family':<12} {'alpha':>8} {'beta':>8} {'R^2':>8}   interpretation", flush=True)
    print("-" * 95, flush=True)
    by_family = {}
    for family, d, s, E, cost in rows_b:
        if cost != cost or cost is None or d < 4 or E < 0:  # NaN / invalid
            continue
        by_family.setdefault(family, []).append((d, cost))
    for family, lst in by_family.items():
        if len(lst) < 3:
            continue
        ds = np.array([d for d, _ in lst], dtype=float)
        cs = np.array([c for _, c in lst], dtype=float)
        log_d = np.log(ds)
        log_log_d = np.log(np.log(ds))
        log_c = np.log(cs)
        # Solve [log d, log log d, 1] @ [alpha, beta, C] = log c
        A = np.stack([log_d, log_log_d, np.ones_like(log_d)], axis=1)
        coef, *_ = np.linalg.lstsq(A, log_c, rcond=None)
        alpha, beta, C = coef
        pred = A @ coef
        ss_res = float(np.sum((log_c - pred) ** 2))
        ss_tot = float(np.sum((log_c - np.mean(log_c)) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else float('nan')
        if alpha < 1.85:
            verdict = f"<<< sub-O(d^2 log d)?  alpha={alpha:.2f}"
        elif alpha < 2.05 and beta < 0.6:
            verdict = "matches O(d^2 log d)"
        else:
            verdict = "no improvement"
        print(f"{family:<12} {alpha:>8.3f} {beta:>8.3f} {r2:>8.3f}   {verdict}", flush=True)


def main():
    t0 = time.perf_counter()
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    experiment_no_fallback()
    rows_b = experiment_scaling()
    experiment_regression(rows_b)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
