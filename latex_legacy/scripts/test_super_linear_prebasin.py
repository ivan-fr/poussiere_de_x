"""
Test two candidate super-linear pre-basin contraction mechanisms on KS,
to see if either rescues the ABD conjecture (N_pre_basin = O(log d)).

Mechanism 1 (LOG-STEP): Newton direction with aggressive line-search on log|P|.
    direction d = P(z) / P'(z)
    z' = z - tau * d, where tau is chosen by line search over {2^k : k = -4..4}
    that minimises log|P(z - tau d)|.

Mechanism 2 (MOBIUS-SCALED): scaled-Steffensen with adaptive probe radius
    s_n = |P(z_n)|^{1/d}.  This is the "natural scale" at which P has unit
    modulus near z_n.  In scaled coordinates w = (z - z_n)/s_n, the polynomial
    P_n(w) = P(z_n + s_n w) has effective gamma(P_n, 0) = s_n * gamma(P, z_n).

Both run with Strategy B' starts (R=2, golden spiral) on KS.
We track alpha-split (pre-basin / basin) via verify_epochs_KS_alpha_split.
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import eval_poly
from verify_armijo_KS import sample_KS
from verify_epochs_KS_alpha_split import smale_alpha, ALPHA_STAR, taylor_coefs

GOLDEN = math.pi * (3.0 - math.sqrt(5.0))


# --- Mechanism 1: log-step (aggressive line-search on log|P|) --------------

def log_step(coefs, z):
    """One log-step: damped Newton with aggressive line search on log|P|.

    Direction = Newton's. Step size tau searched over {2^k : k = -4..4}.
    Returns (z_new, log_contraction = log|P(z_new)| - log|P(z)|).
    """
    cs = taylor_coefs(coefs, z)  # c_0 = P(z), c_1 = P'(z), ...
    if len(cs) < 2:
        return z, 0.0
    P0, P1 = cs[0], cs[1]
    if abs(P1) < 1e-300 or not np.isfinite(P1):
        return z, 0.0
    direction = P0 / P1
    aP0 = abs(P0)
    if aP0 < 1e-300:
        return z, 0.0
    log_aP0 = math.log(aP0 + 1e-300)
    best_z = z
    best_log = log_aP0
    # Aggressive line search: try a wide range of tau
    taus = [2.0**k for k in range(-6, 7)]  # 2^-6 = 0.0156 ... 2^6 = 64
    for tau in taus:
        z_try = z - tau * direction
        try:
            P_try = eval_poly(coefs, z_try)
        except (OverflowError, ZeroDivisionError):
            continue
        if not np.isfinite(P_try):
            continue
        a = abs(P_try)
        if a < 1e-300:
            continue
        log_a = math.log(a + 1e-300)
        if log_a < best_log:
            best_log = log_a
            best_z = z_try
    return best_z, (best_log - log_aP0)


def run_logstep(coefs, z_init, max_epochs=400, rel_tol=1e-12):
    d = len(coefs) - 1
    coef_scale = float(np.max(np.abs(coefs)))
    if coef_scale < 1e-30:
        coef_scale = 1.0
    basin_tol = rel_tol * coef_scale
    z = complex(z_init)
    try:
        P_prev = abs(eval_poly(coefs, z))
    except (OverflowError, ZeroDivisionError):
        return 'overflow', 0, 0
    if P_prev < basin_tol:
        return 'converged', 0, 0

    a_init, _, _ = smale_alpha(coefs, z)
    n_pre = -1
    if np.isfinite(a_init) and a_init < ALPHA_STAR:
        n_pre = 0

    for ep in range(1, max_epochs + 1):
        z_new, log_contr = log_step(coefs, z)
        if z_new == z and log_contr == 0.0:
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        z = z_new
        if n_pre < 0:
            try:
                a_n, _, _ = smale_alpha(coefs, z)
            except (OverflowError, ZeroDivisionError):
                a_n = float('inf')
            if np.isfinite(a_n) and a_n < ALPHA_STAR:
                n_pre = ep
        try:
            P_curr = abs(eval_poly(coefs, z))
        except (OverflowError, ZeroDivisionError):
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        if P_curr < basin_tol:
            return 'converged', ep, n_pre if n_pre >= 0 else ep
    return 'stagnated', max_epochs, n_pre if n_pre >= 0 else max_epochs


# --- Mechanism 2: Mobius scaled-Steffensen ---------------------------------

def mobius_scaled_step(coefs, z):
    """One Mobius scaled-Steffensen step:
        s_n = |P(z)|^{1/d} (natural scale)
        anchor at z, probe at z + s_n
        w_new = -P(z) / Q(z, z+s_n)
        z_new = z + (line-searched tau) * w_new
    """
    d = len(coefs) - 1
    try:
        P0 = eval_poly(coefs, z)
    except (OverflowError, ZeroDivisionError):
        return z, 0.0
    aP0 = abs(P0)
    if aP0 < 1e-300:
        return z, 0.0
    # Natural scale
    try:
        s = aP0 ** (1.0 / d) if d > 0 else 1.0
    except (OverflowError, ValueError):
        s = 1.0
    if s < 1e-12 or not np.isfinite(s):
        s = 1e-3
    # Probe: try 4 directions, pick the one giving best descent
    best_z = z
    best_log = math.log(aP0 + 1e-300)
    for u in [1.0+0j, -1.0+0j, 0+1j, 0-1j]:
        z_probe = z + s * u
        try:
            P_probe = eval_poly(coefs, z_probe)
        except (OverflowError, ZeroDivisionError):
            continue
        if not np.isfinite(P_probe):
            continue
        # Q = (P_probe - P0) / (s * u)
        Q = (P_probe - P0) / (s * u)
        if abs(Q) < 1e-300 or not np.isfinite(Q):
            continue
        # Newton-like step in scaled coords
        full_step = -P0 / Q  # This is a complex displacement in z
        # Line search on log|P(z + tau * full_step)|
        for tau in [1.0, 0.5, 0.25, 0.125, 0.0625, 2.0, 4.0]:
            z_try = z + tau * full_step
            try:
                P_try = eval_poly(coefs, z_try)
            except (OverflowError, ZeroDivisionError):
                continue
            if not np.isfinite(P_try):
                continue
            a = abs(P_try)
            if a < 1e-300:
                a = 1e-300
            log_a = math.log(a)
            if log_a < best_log:
                best_log = log_a
                best_z = z_try
    return best_z, (best_log - math.log(aP0 + 1e-300))


def run_mobius(coefs, z_init, max_epochs=400, rel_tol=1e-12):
    d = len(coefs) - 1
    coef_scale = float(np.max(np.abs(coefs)))
    if coef_scale < 1e-30:
        coef_scale = 1.0
    basin_tol = rel_tol * coef_scale
    z = complex(z_init)
    try:
        P_prev = abs(eval_poly(coefs, z))
    except (OverflowError, ZeroDivisionError):
        return 'overflow', 0, 0
    if P_prev < basin_tol:
        return 'converged', 0, 0

    a_init, _, _ = smale_alpha(coefs, z)
    n_pre = -1
    if np.isfinite(a_init) and a_init < ALPHA_STAR:
        n_pre = 0

    for ep in range(1, max_epochs + 1):
        z_new, log_contr = mobius_scaled_step(coefs, z)
        if z_new == z and log_contr == 0.0:
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        z = z_new
        if n_pre < 0:
            try:
                a_n, _, _ = smale_alpha(coefs, z)
            except (OverflowError, ZeroDivisionError):
                a_n = float('inf')
            if np.isfinite(a_n) and a_n < ALPHA_STAR:
                n_pre = ep
        try:
            P_curr = abs(eval_poly(coefs, z))
        except (OverflowError, ZeroDivisionError):
            return 'fail', ep, n_pre if n_pre >= 0 else ep
        if P_curr < basin_tol:
            return 'converged', ep, n_pre if n_pre >= 0 else ep
    return 'stagnated', max_epochs, n_pre if n_pre >= 0 else max_epochs


# --- Strategy B' starts ----------------------------------------------------

def starts_Bprime(N, R=2.0, q=0.7):
    return [R * (q ** k) * np.exp(1j * k * GOLDEN) for k in range(N)]


def measure(d, n_polys, run_fn, max_starts=24, max_epochs=400, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260426 + d)
    n_totals = []
    n_pres = []
    for _ in range(n_polys):
        coefs = sample_KS(d, rng)
        starts = starts_Bprime(max_starts)
        for z0 in starts:
            try:
                st, n_tot, n_pre = run_fn(coefs, z0, max_epochs=max_epochs)
            except (OverflowError, ZeroDivisionError):
                continue
            if st == 'converged':
                n_totals.append(n_tot)
                n_pres.append(n_pre)
                break
        else:
            n_totals.append(max_epochs)
            n_pres.append(max_epochs)
    return np.array(n_totals, dtype=float), np.array(n_pres, dtype=float)


def report(name, run_fn, degrees, n_polys):
    print("=" * 110, flush=True)
    print(f"Mechanism: {name}   (Strategy B' starts, KS)", flush=True)
    print("=" * 110, flush=True)
    print(f"{'d':>4} | {'<N_pre>':>9} {'med_pre':>8} {'max_pre':>8} | "
          f"{'<N_tot>':>9} {'med_tot':>8} | {'log2 d':>7}", flush=True)
    print("-" * 110, flush=True)
    summary = {}
    for d in degrees:
        n_tot, n_pre = measure(d, n_polys, run_fn)
        m_pre = float(np.mean(n_pre))
        med_pre = int(np.median(n_pre))
        mx_pre = int(np.max(n_pre))
        m_tot = float(np.mean(n_tot))
        med_tot = int(np.median(n_tot))
        ld = math.log2(d)
        print(f"{d:>4} | {m_pre:>9.2f} {med_pre:>8} {mx_pre:>8} | "
              f"{m_tot:>9.2f} {med_tot:>8} | {ld:>7.2f}", flush=True)
        summary[d] = (m_pre, m_tot)
    # Regression
    ds = list(summary.keys())
    log_ds = np.array([math.log2(d) for d in ds])
    pres = np.array([summary[d][0] for d in ds])
    A = np.stack([log_ds, np.ones_like(log_ds)], axis=1)
    (b, a), *_ = np.linalg.lstsq(A, pres, rcond=None)
    pres_pos = np.maximum(pres, 1e-3)
    (b2, a2), *_ = np.linalg.lstsq(A, np.log2(pres_pos), rcond=None)
    print(f"\n  Linear: N_pre = {a:.2f} + {b:.2f} log2(d)", flush=True)
    print(f"  Log-log: N_pre ~ d^{b2:.2f} * 2^{a2:.2f}", flush=True)
    if abs(b2) < 0.20:
        print(f"  -> N_pre = O(log d): ABD-RESCUING", flush=True)
    elif b2 < 0.5:
        print(f"  -> N_pre = sub-sqrt(d): improvement over Pandrosion-T2", flush=True)
    else:
        print(f"  -> N_pre = polynomial in d: ABD still refuted", flush=True)
    return summary


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    DEGREES = [16, 32, 64, 128]
    N_POLYS = 15
    t0 = time.perf_counter()
    print("Comparing super-linear pre-basin candidates on KS\n", flush=True)
    s_log = report("LOG-STEP (aggressive line-search Newton on log|P|)",
                   run_logstep, DEGREES, N_POLYS)
    print()
    s_mob = report("MOBIUS-SCALED (s_n = |P(z)|^{1/d}, scaled Steffensen)",
                   run_mobius, DEGREES, N_POLYS)

    # Comparison vs T2/K=1 baseline (paper 101 v2 alpha-split numbers)
    print("\n" + "=" * 110, flush=True)
    print("Comparison vs T2/K=1 baseline (paper 101 v2)", flush=True)
    print("=" * 110, flush=True)
    baseline = {16: 11.72, 32: 22.18, 64: 38.10, 128: 73.28, 256: 118.12}
    print(f"{'d':>4} | {'T2/K=1':>10} {'LOG-STEP':>10} {'MOBIUS':>10} | "
          f"{'LOG/T2':>8} {'MOB/T2':>8}", flush=True)
    print("-" * 110, flush=True)
    for d in DEGREES:
        b = baseline[d]
        l = s_log[d][0]
        m = s_mob[d][0]
        print(f"{d:>4} | {b:>10.2f} {l:>10.2f} {m:>10.2f} | "
              f"{l/b:>8.2f} {m/b:>8.2f}",
              flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
