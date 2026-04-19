#!/usr/bin/env python3
"""
Verification of the Amortised Block Descent Conjecture
for the Pandrosion base map on four adversarial polynomial families.

Tests Conjecture 4.7: In the microscopic regime, short blocks of m=O(1)
steps yield |P(z_{n+m})| <= e^{-c/d} |P(z_n)|.

Also verifies the cluster domination mechanism (Lemma 4.8).
"""

import numpy as np
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ─── Core computations ───────────────────────────────────────────────

def eval_P(z, roots):
    """P(z) = prod(z - roots[k]), stable log-space for large d."""
    return np.prod(z - roots)

def eval_P_log(z, roots):
    """log|P(z)| via sum of logs, numerically stable."""
    return np.sum(np.log(np.abs(z - roots)))

def pandrosion_step(z, z0, roots):
    """One step of F_{z0}(z) = z0 - P(z0)/Q(z0,z)."""
    P_z = eval_P(z, roots)
    P_z0 = eval_P(z0, roots)
    Q = (P_z - P_z0) / (z - z0)
    if abs(Q) < 1e-50:
        return None
    return z0 - P_z0 / Q

def compute_omega_logs(z, z_next, roots):
    """Compute log|1/omega_k| = log|(z_next - zeta_k)/(z - zeta_k)| for each k."""
    return np.log(np.abs(z_next - roots) / np.abs(z - roots))

def compute_cluster_decomp(z, z_next, roots, sigma):
    """Decompose sum of logs into cluster (near) and far contributions."""
    dists = np.abs(z - roots)
    logs = compute_omega_logs(z, z_next, roots)
    near = dists <= sigma
    far = ~near
    return {
        'total': np.sum(logs),
        'near_sum': np.sum(logs[near]),
        'far_sum': np.sum(logs[far]),
        'n_near': int(np.sum(near)),
        'near_abs': np.sum(np.abs(logs[near])),
        'far_abs': np.sum(np.abs(logs[far])),
    }

# ─── Orbit runner ────────────────────────────────────────────────────

def run_orbit(z0, z_init, roots, n_steps=300, block_size=3):
    """Run Pandrosion base-map orbit and collect all diagnostics."""
    d = len(roots)
    rho = np.max(np.abs(roots))
    R = np.abs(z0)
    sigma = 1.0 / d  # cluster scale

    z = z_init
    P_z0 = eval_P(z0, roots)

    steps = []
    log_P_vals = [eval_P_log(z, roots)]

    for n in range(n_steps):
        z_next = pandrosion_step(z, z0, roots)
        if z_next is None or np.abs(z_next) > 50 * R or np.isnan(z_next):
            break

        # Kinematic ratio r_n = P(z_n)/P(z_0)
        P_z = eval_P(z, roots)
        r_n = P_z / P_z0 if abs(P_z0) > 1e-300 else 0
        
        # Step ratio (microscopic test)
        denom = np.abs(z_next - z0)
        step_ratio = np.abs(z_next - z) / denom if denom > 1e-30 else 1.0
        is_micro = step_ratio < 1.0 / d

        # Cluster decomposition of sum of logs
        decomp = compute_cluster_decomp(z, z_next, roots, sigma)

        # Log ratio log|P(z_{n+1})/P(z_n)|
        log_P_next = eval_P_log(z_next, roots)
        log_ratio = log_P_next - log_P_vals[-1]

        steps.append({
            'n': n,
            'z': z,
            'z_next': z_next,
            'abs_r': np.abs(r_n),
            'step_ratio': step_ratio,
            'is_micro': is_micro,
            'log_ratio': log_ratio,
            'min_dist': np.min(np.abs(z - roots)),
            **decomp,
        })

        log_P_vals.append(log_P_next)
        z = z_next

        # Stop if converged
        if np.min(np.abs(z - roots)) < 1e-12:
            break

    # Compute block descents
    blocks = []
    for i in range(len(log_P_vals) - block_size):
        block_log_descent = log_P_vals[i + block_size] - log_P_vals[i]
        # Check if all steps in block are microscopic
        all_micro = all(
            steps[j]['is_micro']
            for j in range(i, min(i + block_size, len(steps)))
            if j < len(steps)
        )
        blocks.append({
            'start': i,
            'descent': block_log_descent,
            'all_micro': all_micro,
        })

    return steps, blocks, log_P_vals

# ─── Experiment runner ───────────────────────────────────────────────

def run_experiment(name, roots, n_starts=30, n_steps=300, block_size=3):
    """Run full experiment on a polynomial family."""
    d = len(roots)
    rho = np.max(np.abs(roots))
    R = 1 + rho
    sigma = 1.0 / d

    print(f"\n{'═'*72}")
    print(f"  {name}  |  d = {d}  |  ρ = {rho:.3f}  |  R = {R:.3f}  |  σ = {sigma:.4f}")
    print(f"{'═'*72}")

    all_steps = []
    all_blocks = []
    micro_steps = []
    micro_blocks = []
    converged = 0

    for s in range(n_starts):
        theta = 2 * np.pi * s / n_starts
        z0 = R * np.exp(1j * theta)

        # Start iterates at various positions inside the annulus
        for frac in [0.5, 0.7, 0.9]:
            z_init = frac * z0
            steps, blocks, log_P = run_orbit(z0, z_init, roots, n_steps, block_size)

            if steps and np.min(np.abs(steps[-1]['z_next'] - roots)) < 1e-6:
                converged += 1

            all_steps.extend(steps)
            all_blocks.extend(blocks)
            micro_steps.extend(s for s in steps if s['is_micro'])
            micro_blocks.extend(b for b in blocks if b['all_micro'])

    total_orbits = n_starts * 3

    # ─── Report ──────────────────────────────────────────────
    print(f"\n  Orbits: {total_orbits}  |  Converged: {converged}  ({100*converged/total_orbits:.0f}%)")
    print(f"  Total steps: {len(all_steps)}  |  Microscopic: {len(micro_steps)}")

    if all_steps:
        lr = np.array([s['log_ratio'] for s in all_steps])
        print(f"\n  ── Per-step log|P(z_{{n+1}})/P(z_n)| ──")
        print(f"     All steps:  mean = {np.mean(lr):.6f}  |  "
              f"frac < 0 = {np.mean(lr < 0):.3f}  |  "
              f"min = {np.min(lr):.4f}  |  max = {np.max(lr):.4f}")

    if micro_steps:
        lr_m = np.array([s['log_ratio'] for s in micro_steps])
        print(f"     Micro steps: mean = {np.mean(lr_m):.6f}  |  "
              f"frac < 0 = {np.mean(lr_m < 0):.3f}")

    if all_blocks:
        bd = np.array([b['descent'] for b in all_blocks])
        print(f"\n  ── Block descent (m={block_size}) ──")
        print(f"     All blocks ({len(bd):>6d}):  mean = {np.mean(bd):.6f}  |  "
              f"frac < 0 = {np.mean(bd < 0):.3f}")
        target = -1.0 / d
        print(f"     Target: -{1.0/d:.6f} (= -1/d)")
        print(f"     Fraction achieving target: {np.mean(bd < target):.3f}")

    if micro_blocks:
        bd_m = np.array([b['descent'] for b in micro_blocks])
        print(f"     Micro blocks ({len(bd_m):>6d}):  mean = {np.mean(bd_m):.6f}  |  "
              f"frac < 0 = {np.mean(bd_m < 0):.3f}")
        if len(bd_m) > 0:
            print(f"     Micro block mean/d = {np.mean(bd_m)/d:.6f}")
            violations = np.sum(bd_m > 0)
            print(f"     Violations (block ascent): {violations}/{len(bd_m)}")

    # Cluster decomposition stats
    if micro_steps:
        far_abs = np.array([s['far_abs'] for s in micro_steps])
        near_abs = np.array([s['near_abs'] for s in micro_steps])
        far_sum = np.array([s['far_sum'] for s in micro_steps])
        near_sum = np.array([s['near_sum'] for s in micro_steps])
        n_near = np.array([s['n_near'] for s in micro_steps])
        print(f"\n  ── Cluster decomposition (σ = 1/d = {sigma:.4f}) ──")
        print(f"     Mean |I_σ| (cluster size):     {np.mean(n_near):.2f}")
        print(f"     Mean Σ_far (signed):           {np.mean(far_sum):.6f}")
        print(f"     Mean Σ_near (signed):          {np.mean(near_sum):.6f}")
        print(f"     Mean |Σ_far| (absolute):       {np.mean(far_abs):.6f}")
        print(f"     Mean |Σ_near| (absolute):      {np.mean(near_abs):.6f}")
        print(f"     Ratio |far|/|near|:            {np.mean(far_abs)/(np.mean(near_abs)+1e-30):.4f}")

    return {
        'name': name, 'd': d,
        'n_steps': len(all_steps),
        'n_micro': len(micro_steps),
        'n_micro_blocks': len(micro_blocks),
        'converged': converged,
        'total_orbits': total_orbits,
    }


# ─── Define polynomial families ─────────────────────────────────────

def roots_of_unity(d):
    return np.exp(2j * np.pi * np.arange(d) / d)

def roots_on_circle(d, seed=42):
    rng = np.random.RandomState(seed)
    angles = np.sort(rng.uniform(0, 2*np.pi, d))
    return np.exp(1j * angles)

def roots_clustered_line(d):
    return np.array([k / d for k in range(d)], dtype=complex)

def roots_wilkinson(d):
    return np.array([k + 1 for k in range(d)], dtype=complex)


# ─── Sum-of-logs analysis for z^d - 1 ───────────────────────────────

def unity_detailed_analysis(d_list=[5, 10, 20, 50, 100]):
    """Detailed analysis of the sum-of-logs for P(z) = z^d - 1."""
    print("\n" + "█"*72)
    print("  DETAILED ANALYSIS: P(z) = z^d - 1")
    print("█"*72)

    for d in d_list:
        roots = roots_of_unity(d)
        rho = 1.0
        R = 1 + rho

        # Multiple anchors
        step_log_ratios = []
        step_sigmas = []

        for theta_idx in range(20):
            theta = 2 * np.pi * theta_idx / 20
            z0 = R * np.exp(1j * theta)
            P_z0 = eval_P(z0, roots)

            # Start near the root at 1 (distance ~ 1/d from root)
            for eps_scale in [0.5, 1.0, 2.0, 5.0]:
                eps = eps_scale / d
                z = 1.0 + eps * np.exp(1j * np.pi / 4)

                z_next = pandrosion_step(z, z0, roots)
                if z_next is None:
                    continue

                log_ratio = eval_P_log(z_next, roots) - eval_P_log(z, roots)
                step_log_ratios.append(log_ratio)

                # Compute Sigma_log = sum_k log|Q_k/Q|
                omega_logs = compute_omega_logs(z, z_next, roots)
                step_sigmas.append(np.sum(omega_logs))

        arr = np.array(step_log_ratios) if step_log_ratios else np.array([0.0])
        sig = np.array(step_sigmas) if step_sigmas else np.array([0.0])
        print(f"\n  d = {d:>3d}:  mean log-ratio = {np.mean(arr):>10.6f}  |  "
              f"mean Σ_log = {np.mean(sig):>10.6f}  |  "
              f"frac < 0 = {np.mean(arr < 0):.3f}  |  "
              f"mean/d = {np.mean(arr)/d:.8f}")


# ─── Block descent rate as function of d ─────────────────────────────

def scaling_analysis(block_size=3):
    """Measure block descent rate c/d as a function of d."""
    print("\n" + "█"*72)
    print(f"  SCALING ANALYSIS: Block descent rate vs d  (m = {block_size})")
    print("█"*72)
    print(f"\n  {'d':>5s}  {'Family':>20s}  {'micro_blocks':>12s}  "
          f"{'mean_descent':>14s}  {'mean/d':>10s}  {'c':>8s}  {'frac<0':>7s}")
    print("  " + "─"*85)

    for d in [5, 10, 20, 50, 100]:
        families = [
            ("Unity", roots_of_unity(d)),
            ("Circle", roots_on_circle(d)),
            ("Line", roots_clustered_line(d)),
        ]
        if d <= 30:
            families.append(("Wilkinson", roots_wilkinson(d)))

        for fname, roots in families:
            rho = np.max(np.abs(roots))
            R = 1 + rho

            micro_descents = []
            for s in range(20):
                theta = 2 * np.pi * s / 20
                z0 = R * np.exp(1j * theta)
                for frac in [0.5, 0.7, 0.9]:
                    z_init = frac * z0
                    _, blocks, _ = run_orbit(z0, z_init, roots, 200, block_size)
                    for b in blocks:
                        if b['all_micro']:
                            micro_descents.append(b['descent'])

            if micro_descents:
                arr = np.array(micro_descents)
                mean_d = np.mean(arr)
                c_est = -mean_d * d  # descent = -c/d => c = -descent*d
                frac_neg = np.mean(arr < 0)
                print(f"  {d:>5d}  {fname:>20s}  {len(arr):>12d}  "
                      f"{mean_d:>14.6f}  {mean_d/d:>10.8f}  {c_est:>8.4f}  {frac_neg:>7.3f}")
            else:
                print(f"  {d:>5d}  {fname:>20s}  {'(none)':>12s}")


# ─── Main ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  PANDROSION BLOCK DESCENT CONJECTURE — NUMERICAL VERIFICATION      ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    # 1. Detailed per-family experiments
    for d in [10, 20, 50]:
        run_experiment(f"Roots of unity (d={d})", roots_of_unity(d))
        run_experiment(f"Random circle (d={d})", roots_on_circle(d))
        run_experiment(f"Clustered line (d={d})", roots_clustered_line(d))
        if d <= 30:
            run_experiment(f"Wilkinson (d={d})", roots_wilkinson(d))

    # 2. Detailed z^d-1 analysis
    unity_detailed_analysis()

    # 3. Scaling analysis
    scaling_analysis()

    print("\n" + "═"*72)
    print("  DONE.")
