"""
Pandrosion-T_3/K=3 + Strategy A (equispaced ring at the **Cauchy bound**)
on KS, with full alpha-split measurement (paper 101 v2 methodology).

This is the original Part-VII configuration. We expect:
  (a) Strong overflow at large d because rho ~ 2^{d/2} on KS;
  (b) When it does converge, polynomial scaling of N_pre_basin with d.
"""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "latex", "scripts"))
from prove_armijo_O1 import eval_poly, h_step, armijo_fallback, cauchy_root_bound
from verify_armijo_KS import sample_KS
from verify_epochs_KS_T3K3 import T3_step
from verify_epochs_KS_alpha_split import smale_alpha, ALPHA_STAR


def run_T3_K3_armijo_with_alpha(coefs, z_init, max_epochs=400, eta=0.95, rel_tol=1e-12):
    d = len(coefs) - 1
    coef_scale = float(np.max(np.abs(coefs)))
    if coef_scale < 1e-30:
        coef_scale = 1.0
    basin_tol = rel_tol * coef_scale
    K = 3

    z0 = complex(z_init); z = complex(z_init)
    try:
        P_prev = abs(eval_poly(coefs, z0))
    except (OverflowError, ZeroDivisionError):
        return 'overflow', 0, 0, float('inf')
    if not np.isfinite(P_prev):
        return 'overflow', 0, 0, float('inf')
    if P_prev < basin_tol:
        return 'converged', 0, 0, 0.0

    try:
        a_init, _, _ = smale_alpha(coefs, z0)
    except (OverflowError, ZeroDivisionError):
        a_init = float('inf')
    n_pre_basin = -1
    if np.isfinite(a_init) and a_init < ALPHA_STAR:
        n_pre_basin = 0

    epoch_in_anchor = 0
    for ep in range(1, max_epochs + 1):
        try:
            z_cand = T3_step(coefs, z0, z)
        except (OverflowError, ZeroDivisionError):
            z_cand = None
        if z_cand is None or not np.isfinite(z_cand):
            try:
                z_new, _ = armijo_fallback(coefs, z0, d=d)
            except (OverflowError, ZeroDivisionError):
                return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, a_init
            if z_new is None or not np.isfinite(z_new):
                return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, a_init
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
                    return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, a_init
                if z_new is None or not np.isfinite(z_new):
                    return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, a_init
                z = z_new
            else:
                z = z_cand

        if n_pre_basin < 0:
            try:
                a_n, _, _ = smale_alpha(coefs, z)
            except (OverflowError, ZeroDivisionError):
                a_n = float('inf')
            if np.isfinite(a_n) and a_n < ALPHA_STAR:
                n_pre_basin = ep

        try:
            P_curr = abs(eval_poly(coefs, z))
        except (OverflowError, ZeroDivisionError):
            return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, a_init
        if not np.isfinite(P_curr):
            return 'fail', ep, n_pre_basin if n_pre_basin >= 0 else ep, a_init
        if P_curr < basin_tol:
            return 'converged', ep, n_pre_basin if n_pre_basin >= 0 else ep, a_init
        epoch_in_anchor += 1
        if epoch_in_anchor >= K:
            z0 = z
            P_prev = P_curr
            epoch_in_anchor = 0
        else:
            P_prev = P_curr
    return 'stagnated', max_epochs, (n_pre_basin if n_pre_basin >= 0 else max_epochs), a_init


def starts_A_cauchy(coefs, N):
    """Strategy A at the actual Cauchy bound rho = 1 + max|a_k/a_d|."""
    R = cauchy_root_bound(coefs)
    return [R * np.exp(2j * np.pi * k / N) for k in range(N)], R


def measure(d, n_polys, max_starts=None, max_epochs=400, rng=None):
    if rng is None:
        rng = np.random.default_rng(20260426 + d)
    if max_starts is None:
        max_starts = d
    n_totals, n_pres, statuses, rhos = [], [], [], []
    for _ in range(n_polys):
        coefs = sample_KS(d, rng)
        starts, R = starts_A_cauchy(coefs, max_starts)
        rhos.append(R)
        used_status = 'no_start_finite'
        for k, z0 in enumerate(starts):
            try:
                st, n_tot, n_pre, _ = run_T3_K3_armijo_with_alpha(
                    coefs, z0, max_epochs=max_epochs)
            except (OverflowError, ZeroDivisionError):
                continue
            if st == 'converged':
                n_totals.append(n_tot)
                n_pres.append(n_pre)
                statuses.append('converged')
                used_status = 'converged'
                break
            if st == 'overflow':
                continue  # try next start
        else:
            # No start converged -> record as failure
            n_totals.append(max_epochs)
            n_pres.append(max_epochs)
            statuses.append(used_status)
    return (np.array(n_totals, dtype=float),
            np.array(n_pres, dtype=float),
            statuses, np.array(rhos, dtype=float))


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    DEGREES = [16, 32, 64, 128, 256]
    N_POLYS = 30
    t0 = time.perf_counter()
    print("=" * 130, flush=True)
    print("T_3/K=3 + Strategy A (Cauchy bound radius) on KS", flush=True)
    print(f"  alpha_star = {ALPHA_STAR:.5f}, N_POLYS = {N_POLYS}", flush=True)
    print("=" * 130, flush=True)
    print(f"{'d':>4} | {'<rho>':>14} | {'%conv':>6} | "
          f"{'<N_pre>':>9} {'med_pre':>8} {'max_pre':>8} | "
          f"{'<N_tot>':>9} {'med_tot':>8} {'max_tot':>8} | "
          f"{'<N_basin>':>10} {'log2 d':>7}",
          flush=True)
    print("-" * 130, flush=True)
    summary = {}
    for d in DEGREES:
        n_tot, n_pre, statuses, rhos = measure(d, N_POLYS)
        n_conv = sum(1 for s in statuses if s == 'converged')
        if n_conv == 0:
            print(f"{d:>4} | {np.mean(rhos):>14.2e} | "
                  f"{0.0:>5.0f}% | no convergent run", flush=True)
            summary[d] = (None, None, None, None, None, np.mean(rhos), 0)
            continue
        # Restrict statistics to converged runs
        mask = np.array([s == 'converged' for s in statuses])
        n_pre_c = n_pre[mask]
        n_tot_c = n_tot[mask]
        m_pre = float(np.mean(n_pre_c))
        med_pre = int(np.median(n_pre_c))
        mx_pre = int(np.max(n_pre_c))
        m_tot = float(np.mean(n_tot_c))
        med_tot = int(np.median(n_tot_c))
        mx_tot = int(np.max(n_tot_c))
        m_basin = float(np.mean(n_tot_c - n_pre_c))
        rho_mean = float(np.mean(rhos))
        ld = math.log2(d)
        pct = 100.0 * n_conv / N_POLYS
        print(f"{d:>4} | {rho_mean:>14.2e} | "
              f"{pct:>5.0f}% | "
              f"{m_pre:>9.2f} {med_pre:>8} {mx_pre:>8} | "
              f"{m_tot:>9.2f} {med_tot:>8} {mx_tot:>8} | "
              f"{m_basin:>10.2f} {ld:>7.2f}", flush=True)
        summary[d] = (m_pre, m_tot, m_basin, n_pre_c, n_tot_c, rho_mean, pct)

    # Regressions on degrees that converged enough
    valid = [d for d in DEGREES if summary[d][0] is not None]
    if len(valid) >= 3:
        log_ds = np.array([math.log2(d) for d in valid])
        pres = np.array([summary[d][0] for d in valid])
        tots = np.array([summary[d][1] for d in valid])
        basins = np.array([summary[d][2] for d in valid])
        print("\nRegressions over converged degrees:", flush=True)
        for label, ys in [("N_pre_basin", pres), ("N_total", tots), ("N_basin", basins)]:
            A = np.stack([log_ds, np.ones_like(log_ds)], axis=1)
            (b, a), *_ = np.linalg.lstsq(A, ys, rcond=None)
            ys_pos = np.maximum(ys, 1e-3)
            (b2, a2), *_ = np.linalg.lstsq(A, np.log2(ys_pos), rcond=None)
            verdict = ""
            if label == "N_pre_basin":
                verdict = ("  -> O(log d): ABD holds" if abs(b2) < 0.15
                           else "  -> polynomial: ABD refuted")
            print(f"  {label:>13}: y = {a:>7.2f} + {b:>5.2f} log2 d  |  "
                  f"y ~ d^{b2:>5.2f} * 2^{a2:>5.2f}  {verdict}", flush=True)

    print("\n" + "=" * 70, flush=True)
    print("LaTeX-ready table (T_3/K=3 + A Cauchy):", flush=True)
    print("=" * 70, flush=True)
    for d in DEGREES:
        m_pre, m_tot, m_basin, _, _, rho, pct = summary[d]
        if m_pre is None:
            print(f"  ${d}$ & ${rho:.1e}$ & ${pct:.0f}\\%$ & --- & --- & --- & "
                  f"${math.log2(d):.2f}$ \\\\", flush=True)
        else:
            print(f"  ${d}$ & ${rho:.1e}$ & ${pct:.0f}\\%$ & ${m_pre:.2f}$ & "
                  f"${m_basin:.2f}$ & ${m_tot:.2f}$ & ${math.log2(d):.2f}$ \\\\",
                  flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
