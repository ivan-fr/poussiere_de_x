"""
PROPERLY measure the basin-entry epoch count vs total epoch count
on KS for Pandrosion-T2/K=1 with Strategy B' starts.

Tracks alpha(P, z_n) = beta * gamma at each epoch:
  beta(P, z)  = |P(z) / P'(z)|
  gamma(P, z) = max_{k>=2} |c_k / c_1|^{1/(k-1)},  c_k = P^(k)(z)/k!
  alpha*      = (13 - 3 sqrt(17))/4 ~ 0.1577 (Smale)

Definitions (paper 101 v2):
  N_pre_basin = first epoch n such that alpha(P, z_n) < alpha*
  N_basin     = (final convergence epoch) - N_pre_basin
  N_total     = N_pre_basin + N_basin

ABD claim: N_pre_basin = O(log d).
This script tests that claim directly, no longer conflated with N_total.
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import eval_poly, T2_step, armijo_fallback
from verify_armijo_KS import sample_KS

GOLDEN = math.pi * (3.0 - math.sqrt(5.0))
ALPHA_STAR = (13.0 - 3.0 * math.sqrt(17.0)) / 4.0   # ~ 0.15767


# --- Smale alpha theory at a point -----------------------------------------

def taylor_coefs(coefs, z):
    """Return Taylor coefficients c_k = P^(k)(z)/k! for k = 0..d.  Vectorized."""
    poly = np.asarray(coefs, dtype=complex).copy()
    n = len(poly)
    cs = np.zeros(n, dtype=complex)
    for k in range(n):
        try:
            cs[k] = complex(np.polyval(poly[::-1], z))
        except (OverflowError, ZeroDivisionError):
            cs[k] = complex('inf')
        if len(poly) <= 1:
            break
        # Vectorized: new_poly[j-1] = poly[j] * j / (k+1) for j=1..len-1
        idx = np.arange(1, len(poly))
        poly = poly[1:] * idx / (k + 1)
    return cs


def smale_alpha(coefs, z):
    """Return (alpha, beta, gamma) at z. Returns (inf, ..., ...) on degeneracy."""
    cs = taylor_coefs(coefs, z)
    if len(cs) < 2:
        return float('inf'), float('inf'), float('inf')
    c0, c1 = cs[0], cs[1]
    if abs(c1) < 1e-300 or not np.isfinite(c1):
        return float('inf'), float('inf'), float('inf')
    beta = abs(c0) / abs(c1)
    gamma = 0.0
    for k in range(2, len(cs)):
        ck = cs[k]
        if not np.isfinite(ck):
            continue
        ratio = abs(ck) / abs(c1)
        if ratio <= 0.0:
            continue
        # |c_k / c_1|^{1/(k-1)}
        try:
            val = ratio ** (1.0 / (k - 1))
        except (OverflowError, ValueError):
            val = float('inf')
        if val > gamma:
            gamma = val
    if not np.isfinite(gamma):
        return float('inf'), beta, gamma
    return beta * gamma, beta, gamma


# --- T2/K=1 + Armijo, recording alpha at each epoch ------------------------

def run_T2_K1_armijo_with_alpha(coefs, z_init, max_epochs=400, eta=0.95, rel_tol=1e-12):
    """Run Pandrosion-T2/K=1 with Armijo. Return:
        status, n_total, n_pre_basin, alpha_initial.
    n_pre_basin is the smallest n>=0 with alpha(P, z_n) < alpha_star;
    set to n_total if basin entry never observed.
    """
    d = len(coefs) - 1
    coef_scale = float(np.max(np.abs(coefs)))
    if coef_scale < 1e-30:
        coef_scale = 1.0
    basin_tol = rel_tol * coef_scale

    z0 = complex(z_init); z = complex(z_init)
    P_prev = abs(eval_poly(coefs, z0))
    if P_prev < basin_tol:
        return 'converged', 0, 0, 0.0

    # alpha at start (epoch 0)
    alpha0, _, _ = smale_alpha(coefs, z0)
    alpha_init = alpha0
    n_pre_basin = -1
    if alpha0 < ALPHA_STAR:
        n_pre_basin = 0

    for ep in range(1, max_epochs + 1):
        try:
            z_cand = T2_step(coefs, z0, z)
        except (OverflowError, ZeroDivisionError):
            z_cand = None
        if z_cand is None or not np.isfinite(z_cand):
            try:
                z_new, _ = armijo_fallback(coefs, z0, d=d)
            except (OverflowError, ZeroDivisionError):
                return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, alpha_init
            if z_new is None or not np.isfinite(z_new):
                return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, alpha_init
            z = z_new
        else:
            try:
                pcand = abs(eval_poly(coefs, z_cand))
            except (OverflowError, ZeroDivisionError):
                pcand = float('inf')
            if pcand > eta * P_prev:
                try:
                    z_new, _ = armijo_fallback(coefs, z0, d=d)
                except (OverflowError, ZeroDivisionError):
                    return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, alpha_init
                if z_new is None or not np.isfinite(z_new):
                    return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, alpha_init
                z = z_new
            else:
                z = z_cand

        # Probe alpha at the new iterate (only if we have not yet entered basin)
        if n_pre_basin < 0:
            alpha_n, _, _ = smale_alpha(coefs, z)
            if np.isfinite(alpha_n) and alpha_n < ALPHA_STAR:
                n_pre_basin = ep

        try:
            P_curr = abs(eval_poly(coefs, z))
        except (OverflowError, ZeroDivisionError):
            return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, alpha_init
        if P_curr < basin_tol:
            return 'converged', ep, n_pre_basin if n_pre_basin >= 0 else ep, alpha_init
        z0 = z
        P_prev = P_curr

    # Stagnated
    return 'stagnated', max_epochs, (n_pre_basin if n_pre_basin >= 0 else max_epochs), alpha_init


def starts_Bprime(N, R=2.0, q=0.7):
    return [R * (q ** k) * np.exp(1j * k * GOLDEN) for k in range(N)]


def measure(d, n_polys, max_starts=24, max_epochs=400, rng=None):
    """Returns lists (N_total, N_pre_basin, alpha_initial) for first converging orbit per poly."""
    if rng is None:
        rng = np.random.default_rng(20260426 + d)
    n_totals = []
    n_pres = []
    alphas0 = []
    for _ in range(n_polys):
        coefs = sample_KS(d, rng)
        starts = starts_Bprime(max_starts)
        for k, z0 in enumerate(starts):
            try:
                st, n_tot, n_pre, a0 = run_T2_K1_armijo_with_alpha(
                    coefs, z0, max_epochs=max_epochs)
            except (OverflowError, ZeroDivisionError):
                continue
            if st == 'converged':
                n_totals.append(n_tot)
                n_pres.append(n_pre)
                alphas0.append(a0)
                break
        else:
            n_totals.append(max_epochs)
            n_pres.append(max_epochs)
            alphas0.append(float('inf'))
    return (np.array(n_totals, dtype=float),
            np.array(n_pres, dtype=float),
            np.array(alphas0, dtype=float))


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    DEGREES = [16, 32, 64, 128, 256]   # d=512 dropped: taylor_coefs O(d^2) too slow there
    N_POLYS = 40
    t0 = time.perf_counter()
    print("=" * 120, flush=True)
    print("alpha-split epoch measurement on KS  (paper 101 v2, ABD properly tested)",
          flush=True)
    print(f"  alpha_star = {ALPHA_STAR:.5f}, N_POLYS = {N_POLYS}", flush=True)
    print("=" * 120, flush=True)
    header = (f"{'d':>4} | {'<N_pre>':>9} {'med_pre':>8} {'max_pre':>8} | "
              f"{'<N_tot>':>9} {'med_tot':>8} {'max_tot':>8} | "
              f"{'<N_basin>':>10} | {'pre_pct':>8} {'log2 d':>7}")
    print(header, flush=True)
    print("-" * 120, flush=True)
    summary = {}
    for d in DEGREES:
        n_tot, n_pre, a0 = measure(d, N_POLYS)
        m_pre = float(np.mean(n_pre))
        med_pre = int(np.median(n_pre))
        mx_pre = int(np.max(n_pre))
        m_tot = float(np.mean(n_tot))
        med_tot = int(np.median(n_tot))
        mx_tot = int(np.max(n_tot))
        m_basin = float(np.mean(n_tot - n_pre))
        pre_pct = 100.0 * m_pre / m_tot if m_tot > 0 else 0.0
        ld = math.log2(d)
        print(f"{d:>4} | {m_pre:>9.2f} {med_pre:>8} {mx_pre:>8} | "
              f"{m_tot:>9.2f} {med_tot:>8} {mx_tot:>8} | "
              f"{m_basin:>10.2f} | {pre_pct:>7.1f}% {ld:>7.2f}",
              flush=True)
        summary[d] = (m_pre, m_tot, m_basin, n_pre, n_tot, a0)
    print("-" * 120, flush=True)

    # Regressions
    ds = list(summary.keys())
    log_ds = np.array([math.log2(d) for d in ds])
    pres = np.array([summary[d][0] for d in ds])
    tots = np.array([summary[d][1] for d in ds])
    basins = np.array([summary[d][2] for d in ds])

    print("\nRegressions  (decisive test of ABD: N_pre_basin = O(log d) ?)", flush=True)
    for label, ys in [("N_pre_basin", pres), ("N_total", tots), ("N_basin", basins)]:
        # Linear fit: y = a + b * log2(d)
        A = np.stack([log_ds, np.ones_like(log_ds)], axis=1)
        (b, a), *_ = np.linalg.lstsq(A, ys, rcond=None)
        # Log-log fit: log y = a' + b' log d
        ys_pos = np.maximum(ys, 1e-3)
        log_ys = np.log2(ys_pos)
        (b2, a2), *_ = np.linalg.lstsq(A, log_ys, rcond=None)
        verdict = ""
        if label == "N_pre_basin":
            if abs(b2) < 0.15:
                verdict = "  -> consistent with O(log d): ABD holds"
            else:
                verdict = "  -> grows polynomially: ABD refuted"
        print(f"  {label:>13}:  linear  y = {a:>7.2f} + {b:>5.2f} log2 d   "
              f"|  log-log  y ~ d^{b2:>5.2f} * 2^{a2:>5.2f}   {verdict}",
              flush=True)

    # Sanity check: are starts already in basin frequently?
    print("\nFraction of orbits already in basin at start (n_pre_basin = 0):", flush=True)
    for d in DEGREES:
        _, _, _, n_pre_arr, _, _ = summary[d]
        frac = float(np.sum(n_pre_arr == 0)) / len(n_pre_arr)
        print(f"  d = {d}:  {frac*100:.1f}%", flush=True)

    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table:", flush=True)
    print("=" * 70, flush=True)
    for d in DEGREES:
        m_pre, m_tot, m_basin, _, _, _ = summary[d]
        print(f"  ${d}$ & ${m_pre:.2f}$ & ${m_basin:.2f}$ & ${m_tot:.2f}$ & "
              f"${100*m_pre/max(m_tot, 1e-9):.1f}\\%$ & ${math.log2(d):.2f}$ \\\\",
              flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
