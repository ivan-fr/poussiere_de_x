#!/usr/bin/env python3
"""
Phase 2: Refined verification with:
1. Adaptive anchor (K=1 = Newton) vs fixed anchor
2. Anchor at intermediate radius (not just Cauchy circle)
3. Precise scaling of the descent constant c as a function of d
4. Full 1/omega_k decomposition
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def eval_P(z, roots):
    return np.prod(z - roots)

def eval_P_log(z, roots):
    return np.sum(np.log(np.abs(z - roots)))

def pandrosion_step_fixed(z, z0, roots):
    """Base map with FIXED anchor z0."""
    P_z = eval_P(z, roots)
    P_z0 = eval_P(z0, roots)
    Q = (P_z - P_z0) / (z - z0)
    if abs(Q) < 1e-50:
        return None
    return z0 - P_z0 / Q

def newton_step(z, roots):
    """Newton step = Pandrosion with z0 = z (K=1)."""
    d = len(roots)
    P_z = eval_P(z, roots)
    P_prime = P_z * np.sum(1.0 / (z - roots))
    if abs(P_prime) < 1e-50:
        return None
    return z - P_z / P_prime

def run_adaptive_orbit(z_init, roots, n_steps=200):
    """Pandrosion T3 with K=1 (= Newton with Steffensen acceleration)."""
    z = z_init
    d = len(roots)
    log_Ps = [eval_P_log(z, roots)]
    
    for n in range(n_steps):
        # Base map step (Newton)
        z_next = newton_step(z, roots)
        if z_next is None or np.abs(z_next) > 100 or np.isnan(z_next):
            break
        log_Ps.append(eval_P_log(z_next, roots))
        z = z_next
        if np.min(np.abs(z - roots)) < 1e-12:
            break
    return log_Ps

def run_base_map_orbit(z_init, z0, roots, n_steps=500):
    """Base map with fixed anchor z0."""
    z = z_init
    d = len(roots)
    log_Ps = [eval_P_log(z, roots)]
    
    for n in range(n_steps):
        z_next = pandrosion_step_fixed(z, z0, roots)
        if z_next is None or np.abs(z_next) > 100 or np.isnan(z_next):
            break
        log_Ps.append(eval_P_log(z_next, roots))
        z = z_next
        if np.min(np.abs(z - roots)) < 1e-12:
            break
    return log_Ps


# ═══════════════════════════════════════════════════════════════════
# Experiment 1: Fixed anchor at R = 1+ρ vs closer anchors
# ═══════════════════════════════════════════════════════════════════

def exp1_anchor_radius_effect():
    print("\n" + "█"*72)
    print("  EXP 1: Effect of anchor radius on per-step descent")
    print("█"*72)
    print(f"\n  {'d':>4s}  {'R_factor':>10s}  {'R':>6s}  {'mean_log_ratio':>16s}  "
          f"{'c_eff':>8s}  {'frac<0':>7s}  {'n_steps':>8s}")
    print("  " + "─"*70)

    for d in [10, 20, 50, 100]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        rho = 1.0

        for r_factor in [2.0, 1.5, 1.2, 1.05, 1.01]:
            R = r_factor  # since rho = 1, R = 1+rho would be 2
            log_ratios = []

            for s in range(20):
                theta = 2 * np.pi * s / 20
                z0 = R * np.exp(1j * theta)
                    
                for frac in [0.5, 0.8]:
                    z_init = frac * np.exp(1j * (theta + 0.3))
                    log_Ps = run_base_map_orbit(z_init, z0, roots, 500)
                    
                    for i in range(len(log_Ps)-1):
                        lr = log_Ps[i+1] - log_Ps[i]
                        if np.isfinite(lr):
                            log_ratios.append(lr)

            if log_ratios:
                arr = np.array(log_ratios)
                arr = arr[np.isfinite(arr)]
                if len(arr) > 0:
                    mean_lr = np.mean(arr)
                    c_eff = -mean_lr * d
                    frac_neg = np.mean(arr < 0)
                    print(f"  {d:>4d}  {r_factor:>10.2f}  {R:>6.2f}  {mean_lr:>16.8f}  "
                          f"{c_eff:>8.4f}  {frac_neg:>7.3f}  {len(arr):>8d}")


# ═══════════════════════════════════════════════════════════════════
# Experiment 2: Newton (adaptive K=1) descent rate
# ═══════════════════════════════════════════════════════════════════

def exp2_newton_descent():
    print("\n" + "█"*72)
    print("  EXP 2: Newton (K=1 adaptive) descent rate vs d")
    print("█"*72)
    print(f"\n  {'d':>4s}  {'Family':>15s}  {'mean_descent':>14s}  {'c_eff':>8s}  "
          f"{'frac<0':>7s}  {'conv':>6s}  {'n_steps':>8s}")
    print("  " + "─"*70)

    for d in [5, 10, 20, 50, 100, 200]:
        families = [
            ("Unity", np.exp(2j * np.pi * np.arange(d) / d)),
            ("Line", np.array([k/d for k in range(d)], dtype=complex)),
        ]
        if d <= 50:
            families.append(("Wilkinson", np.array([k+1 for k in range(d)], dtype=complex)))

        for fname, roots in families:
            rho = np.max(np.abs(roots))
            R = 1 + rho
            
            all_lr = []
            conv = 0
            total = 0

            for s in range(30):
                theta = 2 * np.pi * s / 30
                z_init = R * np.exp(1j * theta)
                total += 1

                log_Ps = run_adaptive_orbit(z_init, roots, 300)
                
                if len(log_Ps) > 1:
                    if eval_P_log(1.0, roots) > -20:  # roughly means z is near a root
                        pass
                    
                    for i in range(min(len(log_Ps)-1, 100)):  # pre-lock-in steps
                        lr = log_Ps[i+1] - log_Ps[i]
                        if np.isfinite(lr):
                            all_lr.append(lr)
                    
                    # Check convergence
                    # Find if last iterate is near a root
                    # Approximate: if log|P| dropped a lot
                    if log_Ps[-1] < -20:
                        conv += 1

            if all_lr:
                arr = np.array(all_lr)
                arr = arr[np.isfinite(arr)]
                mean_lr = np.mean(arr)
                c_eff = -mean_lr * d
                frac_neg = np.mean(arr < 0)
                print(f"  {d:>4d}  {fname:>15s}  {mean_lr:>14.6f}  {c_eff:>8.4f}  "
                      f"{frac_neg:>7.3f}  {conv:>3d}/{total:<2d}  {len(arr):>8d}")


# ═══════════════════════════════════════════════════════════════════
# Experiment 3: Precise 1/omega_k analysis for z^d - 1
# ═══════════════════════════════════════════════════════════════════

def exp3_omega_analysis():
    print("\n" + "█"*72)
    print("  EXP 3: Detailed 1/ω_k structure for P(z) = z^d - 1")
    print("█"*72)

    for d in [10, 20, 50]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        rho = 1.0

        print(f"\n  ── d = {d} ──")
        
        # Pick specific anchor and iterate near root
        z0 = 2.0 + 0j  # anchor on positive real axis
        
        for eps in [0.1, 0.01, 0.001]:
            z = 1.0 + eps * np.exp(1j * np.pi/4)
            
            z_next = pandrosion_step_fixed(z, z0, roots)
            if z_next is None:
                continue
            
            # Compute 1/omega_k
            inv_omega = (z_next - roots) / (z - roots)
            log_inv_omega = np.log(np.abs(inv_omega))
            
            # Kinematic ratio
            P_z = eval_P(z, roots)
            P_z0 = eval_P(z0, roots)
            r = P_z / P_z0
            
            # Delta_k = (z_next - z)/(z - zeta_k)
            delta_k = (z_next - z) / (z - roots)
            
            total_log = np.sum(log_inv_omega)
            
            # Signed decomposition: near root 0 (= root at 1)
            dists = np.abs(z - roots)
            sorted_idx = np.argsort(dists)
            
            print(f"\n     ε = {eps:.4f}  |  |r| = {abs(r):.2e}  |  "
                  f"|step| = {abs(z_next - z):.2e}  |  "
                  f"Σ log|1/ω| = {total_log:.8f}")
            print(f"     Top 5 |Δ_k|: {np.sort(np.abs(delta_k))[-5:][::-1]}")
            print(f"     Top 5 |log(1/ω_k)|: {np.sort(np.abs(log_inv_omega))[-5:][::-1]}")
            
            # Check: does sum equal log|P(z_next)/P(z)|?
            check = eval_P_log(z_next, roots) - eval_P_log(z, roots)
            print(f"     log|P(z')/P(z)| = {check:.8f}  (should equal Σ = {total_log:.8f})")


# ═══════════════════════════════════════════════════════════════════
# Experiment 4: What scaling of c actually works?
# ═══════════════════════════════════════════════════════════════════

def exp4_scaling():
    print("\n" + "█"*72)
    print("  EXP 4: Precise scaling of descent rate for BASE MAP")
    print("█"*72)
    print(f"\n  Base map descent: log|P(z_{{n+1}})/P(z_n)| ~ -c/d^α")
    print(f"\n  {'d':>4s}  {'mean(-lr)':>12s}  {'c (if α=1)':>12s}  {'c (if α=2)':>12s}  "
          f"{'log_d':>6s}  {'c/log_d':>10s}")
    print("  " + "─"*65)

    ds = [5, 10, 20, 50, 100, 200]
    means = []

    for d in ds:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        
        all_lr = []
        for s in range(30):
            theta = 2 * np.pi * s / 30
            z0 = R * np.exp(1j * theta)
            for frac in [0.3, 0.5, 0.7, 0.9]:
                z_init = frac * np.exp(1j * (theta + 0.5))
                log_Ps = run_base_map_orbit(z_init, z0, roots, 500)
                for i in range(len(log_Ps) - 1):
                    lr = log_Ps[i+1] - log_Ps[i]
                    if np.isfinite(lr) and abs(lr) < 10:
                        all_lr.append(lr)
        
        if all_lr:
            arr = np.array(all_lr)
            mean_neg = -np.mean(arr)
            c1 = mean_neg * d
            c2 = mean_neg * d**2
            log_d = np.log(d)
            means.append((d, mean_neg))
            print(f"  {d:>4d}  {mean_neg:>12.8f}  {c1:>12.6f}  {c2:>12.4f}  "
                  f"{log_d:>6.2f}  {c1/log_d:>10.6f}")

    # Power-law fit
    if len(means) >= 3:
        ds_arr = np.array([m[0] for m in means])
        mn_arr = np.array([m[1] for m in means])
        valid = mn_arr > 0
        if np.sum(valid) >= 3:
            log_d = np.log(ds_arr[valid])
            log_m = np.log(mn_arr[valid])
            # Linear regression: log(mean) = a + b*log(d) => mean ~ d^b
            coeffs = np.polyfit(log_d, log_m, 1)
            print(f"\n  Power-law fit: descent ~ d^{{{coeffs[0]:.3f}}}")
            print(f"  Equivalently: descent ~ 1/d^{{{-coeffs[0]:.3f}}}")


# ═══════════════════════════════════════════════════════════════════
# Experiment 5: Amortized descent with REANCHORING
# ═══════════════════════════════════════════════════════════════════

def exp5_reanchoring():
    """Test block descent with anchor updated every K steps."""
    print("\n" + "█"*72)
    print("  EXP 5: Block descent with reanchoring (K = 1, 3, 10, ∞)")  
    print("█"*72)
    print(f"\n  {'d':>4s}  {'K':>4s}  {'mean_block_desc':>16s}  {'c (×d)':>10s}  "
          f"{'frac<0':>7s}  {'conv%':>6s}")
    print("  " + "─"*55)

    for d in [10, 20, 50]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        rho = 1.0
        R = 2.0
        block_size = 3

        for K in [1, 3, 10, 1000]:  # K=1000 ≈ infinity (fixed anchor)
            all_block_desc = []
            conv = 0
            total = 0

            for s in range(30):
                theta = 2 * np.pi * s / 30
                z0 = R * np.exp(1j * theta)
                
                for frac in [0.5, 0.7, 0.9]:
                    z = frac * np.exp(1j * (theta + 0.3))
                    total += 1
                    
                    log_Ps = [eval_P_log(z, roots)]
                    current_z0 = z0
                    
                    for step in range(300):
                        if K == 1:
                            z_next = newton_step(z, roots)
                        else:
                            z_next = pandrosion_step_fixed(z, current_z0, roots)
                        
                        if z_next is None or np.abs(z_next) > 100 or np.isnan(z_next):
                            break
                        
                        log_Ps.append(eval_P_log(z_next, roots))
                        z = z_next
                        
                        # Reanchor every K steps
                        if K < 1000 and (step + 1) % K == 0:
                            current_z0 = z
                        
                        if np.min(np.abs(z - roots)) < 1e-10:
                            conv += 1
                            break
                    
                    # Block descents
                    for i in range(len(log_Ps) - block_size):
                        bd = log_Ps[i + block_size] - log_Ps[i]
                        if np.isfinite(bd) and abs(bd) < 100:
                            all_block_desc.append(bd)

            if all_block_desc:
                arr = np.array(all_block_desc)
                mean_bd = np.mean(arr)
                c_eff = -mean_bd * d
                frac_neg = np.mean(arr < 0)
                conv_pct = f"{100*conv/total:.0f}%"
                K_str = "∞" if K == 1000 else str(K)
                print(f"  {d:>4d}  {K_str:>4s}  {mean_bd:>16.6f}  {c_eff:>10.4f}  "
                      f"{frac_neg:>7.3f}  {conv_pct:>6s}")


if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  PANDROSION BLOCK DESCENT — REFINED EXPERIMENTS                    ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    
    exp1_anchor_radius_effect()
    exp2_newton_descent()
    exp3_omega_analysis()
    exp4_scaling()
    exp5_reanchoring()
    
    print("\n" + "═"*72)
    print("  ALL EXPERIMENTS DONE.")
