"""Diagnose what happens at (n=15, d=4) D=10^9."""
from __future__ import annotations
import math, time, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from multivariate_KS_fast import (
    multi_indices_array, multinomial_log_arr, sample_KS_mv_fast,
    starts_Bprime_mv_fast, eval_F_fast, jacobian_F_fast, newton_step_fast,
    els_step_fast,
)


def diag(n, d, n_polys=2):
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print(f"=== Diag for (n={n}, d={d}) D={d**n}", flush=True)
    rng = np.random.default_rng(20260427 + 100*n + d)
    alpha_arr = multi_indices_array(n, d)
    sigmas = np.exp(0.5 * multinomial_log_arr(d, alpha_arr))
    print(f"  n_mons={alpha_arr.shape[0]}, alpha_arr range = [{alpha_arr.min()}, {alpha_arr.max()}]",
          flush=True)
    print(f"  sigmas range = [{sigmas.min():.2e}, {sigmas.max():.2e}]", flush=True)

    for k in range(n_polys):
        print(f"\n  -- Poly {k} --", flush=True)
        coefs_arr = sample_KS_mv_fast(n, d, alpha_arr, sigmas, rng)
        starts = starts_Bprime_mv_fast(6, n)
        for s_idx, z0 in enumerate(starts):
            t0 = time.perf_counter()
            try:
                Fz = eval_F_fast(coefs_arr, alpha_arr, z0)
                t_eval = time.perf_counter() - t0
                norm_F = float(np.linalg.norm(Fz))
                print(f"    start {s_idx}: ||z||={float(np.linalg.norm(z0)):.3f}, "
                      f"||F(z)||={norm_F:.2e}, eval={t_eval*1000:.1f}ms", flush=True)

                t0 = time.perf_counter()
                Fz_, Delta = newton_step_fast(coefs_arr, alpha_arr, z0)
                t_jac = time.perf_counter() - t0
                if Delta is None:
                    print(f"      Newton step FAILED (Delta=None)", flush=True)
                    continue
                norm_Delta = float(np.linalg.norm(Delta))
                print(f"      ||Delta||={norm_Delta:.2e}, "
                      f"newton_step={t_jac*1000:.1f}ms", flush=True)

                # Try one ELS step
                t0 = time.perf_counter()
                T_set = [2.0**kk for kk in range(-6, 7)]
                z_new, log_contr = els_step_fast(coefs_arr, alpha_arr, z0, T_set)
                t_els = time.perf_counter() - t0
                Fz_new = eval_F_fast(coefs_arr, alpha_arr, z_new)
                norm_F_new = float(np.linalg.norm(Fz_new))
                print(f"      after 1 ELS: ||F||={norm_F_new:.2e} (contraction {norm_F_new/norm_F:.2e}), "
                      f"els={t_els*1000:.1f}ms", flush=True)

                # Run more epochs to see what happens
                z = z_new
                for ep in range(1, 30):
                    t0 = time.perf_counter()
                    z, _ = els_step_fast(coefs_arr, alpha_arr, z, T_set)
                    t_step = time.perf_counter() - t0
                    Fz_now = eval_F_fast(coefs_arr, alpha_arr, z)
                    norm = float(np.linalg.norm(Fz_now))
                    if not np.isfinite(norm):
                        print(f"      ep {ep+1}: NON-FINITE", flush=True)
                        break
                    if ep <= 5 or ep % 5 == 0:
                        print(f"      ep {ep+1}: ||F||={norm:.2e}, step={t_step*1000:.1f}ms",
                              flush=True)
                    if norm < 1e-6:
                        print(f"      CONVERGED at ep {ep+1}", flush=True)
                        break
                break  # only test first start per poly
            except Exception as e:
                print(f"    EXCEPTION: {e}", flush=True)
                continue


if __name__ == "__main__":
    diag(15, 4, n_polys=2)
