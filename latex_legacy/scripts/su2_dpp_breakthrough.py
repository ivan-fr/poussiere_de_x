"""
Test the HYPOTHESIS:  log|D_AB|^2 - log det K  ~  L_n  asymptotically?
And:  log det G_norm  ~  log det K  for large n?

If TRUE, Atiyah-Sutcliffe reduces to a question about the
SU(2)-coherent-state DPP kernel K_{ij} = <u_i, u_j>^{n-1}, which has
KNOWN positivity / concentration properties from the DPP literature.
"""
from __future__ import annotations
import math
import numpy as np
from su2_atiyah_toolkit import (
    hopf_lift, atiyah_D_squared, gram_normalized,
    reproducing_kernel_matrix, random_S2_points, antipodal_split,
    spinor_inner,
)


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("HYPOTHESIS: log|D|^2 = L_n + log det K asymptotically?", flush=True)
    print("=" * 100, flush=True)

    # Random configurations
    print(f"\n{'family':>10} {'n':>4} {'log|D|^2':>10} {'log det K':>11} "
          f"{'log det G':>11} {'L_n':>10} "
          f"{'log|D|^2 - log det K':>22} {'... - L_n':>11}",
          flush=True)
    for family in ['random', 'antipod']:
        for n in [4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40]:
            rng = np.random.default_rng(2026 + n + (1000 if family == 'antipod' else 0))
            n_cfg = 30 if n < 20 else 10
            log_D_sqs = []
            log_det_Ks = []
            log_det_Gs = []
            for _ in range(n_cfg):
                if family == 'random':
                    pts = random_S2_points(rng, n)
                else:
                    pts = antipodal_split(rng, n)
                logD_sq, L_n = atiyah_D_squared(pts)
                spinors = [hopf_lift(x) for x in pts]
                K = reproducing_kernel_matrix(spinors, n)
                K_h = (K + K.conj().T) / 2
                sg, ld_K = np.linalg.slogdet(K_h)
                Gn = gram_normalized(pts)
                sg2, ld_G = np.linalg.slogdet(Gn)
                if not np.isfinite(logD_sq) or not np.isfinite(ld_K) or not np.isfinite(ld_G):
                    continue
                log_D_sqs.append(logD_sq)
                log_det_Ks.append(float(np.real(ld_K)))
                log_det_Gs.append(float(ld_G))
            if not log_D_sqs:
                continue
            arr_D = np.array(log_D_sqs)
            arr_K = np.array(log_det_Ks)
            arr_G = np.array(log_det_Gs)
            binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
            L_n = float(np.sum(np.log(binom)))
            mean_D = arr_D.mean()
            mean_K = arr_K.mean()
            mean_G = arr_G.mean()
            diff = arr_D - arr_K
            print(f"{family:>10} {n:>4} {mean_D:>10.4f} {mean_K:>11.4f} "
                  f"{mean_G:>11.4f} {L_n:>10.2f} "
                  f"{diff.mean():>22.4f} {(diff.mean() - L_n):>11.4f}",
                  flush=True)

    # Sharper test: is log det G_norm = log det K + correction?
    print("\n" + "=" * 100, flush=True)
    print("Sharper: is log det G_norm = log det K + alpha_n? (alpha_n = const of n)",
          flush=True)
    print("=" * 100, flush=True)
    print(f"\n{'n':>4} {'mean log G':>12} {'mean log K':>12} {'diff (log G - log K)':>22} "
          f"{'std diff':>10} {'pcorr':>8}", flush=True)
    for n in [4, 5, 6, 8, 10, 12, 15, 20]:
        rng = np.random.default_rng(20260901 + n)
        log_Gs = []; log_Ks = []
        for _ in range(60):
            pts = random_S2_points(rng, n)
            spinors = [hopf_lift(x) for x in pts]
            K = reproducing_kernel_matrix(spinors, n)
            K_h = (K + K.conj().T) / 2
            sg1, ld_K = np.linalg.slogdet(K_h)
            Gn = gram_normalized(pts)
            sg2, ld_G = np.linalg.slogdet(Gn)
            if np.isfinite(ld_K) and np.isfinite(ld_G):
                log_Ks.append(float(np.real(ld_K)))
                log_Gs.append(float(ld_G))
        if not log_Ks:
            continue
        arr_K = np.array(log_Ks); arr_G = np.array(log_Gs)
        diff = arr_G - arr_K
        corr = float(np.corrcoef(arr_G, arr_K)[0, 1]) if len(arr_K) >= 3 else 0.0
        print(f"{n:>4} {arr_G.mean():>12.4f} {arr_K.mean():>12.4f} "
              f"{diff.mean():>22.4f} {diff.std():>10.4f} {corr:>8.4f}",
              flush=True)

    # Test: does the DPP kernel K bound G_norm in the sense det K <= det G_norm?
    # Or det G_norm <= det K?
    print("\n" + "=" * 100, flush=True)
    print("Test: which is bigger, det G_norm or det K (normalized)?", flush=True)
    print("=" * 100, flush=True)
    # K has diagonal K_ii = <u_i, u_i>^{n-1} = 1 (since u_i unit).
    # So K has unit diagonal, just like G_norm! Both PSD with unit diag.
    # Hadamard: det K <= 1 and det G_norm <= 1.
    print(f"\n{'n':>4} {'log det K (avg)':>16} {'log det G (avg)':>17} "
          f"{'K bigger?':>12}",
          flush=True)
    for n in [4, 5, 6, 8, 10, 12, 15, 20]:
        rng = np.random.default_rng(20261001 + n)
        log_Ks, log_Gs = [], []
        for _ in range(40):
            pts = random_S2_points(rng, n)
            spinors = [hopf_lift(x) for x in pts]
            K = reproducing_kernel_matrix(spinors, n)
            K_h = (K + K.conj().T) / 2
            # Verify K_ii = 1
            assert np.max(np.abs(np.real(np.diag(K)) - 1.0)) < 1e-12
            sg1, ld_K = np.linalg.slogdet(K_h)
            Gn = gram_normalized(pts)
            sg2, ld_G = np.linalg.slogdet(Gn)
            if np.isfinite(ld_K) and np.isfinite(ld_G):
                log_Ks.append(float(np.real(ld_K)))
                log_Gs.append(float(ld_G))
        arr_K = np.array(log_Ks); arr_G = np.array(log_Gs)
        bigger = "K" if arr_K.mean() > arr_G.mean() else "G_norm"
        print(f"{n:>4} {arr_K.mean():>16.4f} {arr_G.mean():>17.4f} {bigger:>12}",
              flush=True)


if __name__ == "__main__":
    main()
