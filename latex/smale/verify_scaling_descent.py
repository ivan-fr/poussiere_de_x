#!/usr/bin/env python3
"""
Phase 3: Scaling-aware Pandrosion-T3 descent verification.

Key insight from pandrosion_en_improved.tex: the contraction ratio λ_{p,x}
is monotone in x, and goes to 0 as x' → 1. The scaling principle
x^{1/p} = A^{1/p} · (x/A)^{1/p} reduces the problem to x' ≈ 1.

For general polynomials, the ADAPTIVE anchor z_0 = z_n is the analogue:
it keeps the divided difference Q(z_0,z) ≈ P'(z), making the effective
ratio P(z)/P(z_0) ≈ 1 (the "easy" regime).

This script verifies block descent for:
1. Pandrosion-T3 with K=1 (Newton = maximally adaptive)
2. Pandrosion-T3 with K=3 (anchor updated every 3 steps)
3. Fixed anchor at various radii (showing the scaling obstruction)
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def eval_P(z, roots):
    return np.prod(z - roots)

def eval_P_log(z, roots):
    return np.sum(np.log(np.abs(z - roots)))

def pandrosion_step(z, z0, roots):
    """F_{z0}(z) = z0 - P(z0)/Q(z0,z)"""
    P_z = eval_P(z, roots)
    P_z0 = eval_P(z0, roots)
    if abs(z - z0) < 1e-30:
        # Degenerate: use Newton instead
        P_prime = P_z * np.sum(1.0 / (z - roots))
        if abs(P_prime) < 1e-50:
            return None
        return z - P_z / P_prime
    Q = (P_z - P_z0) / (z - z0)
    if abs(Q) < 1e-50:
        return None
    return z0 - P_z0 / Q

def steffensen_step(z, z0, roots):
    """Steffensen (T3) acceleration: three-point Aitken on Pandrosion."""
    z1 = pandrosion_step(z, z0, roots)
    if z1 is None: return None
    z2 = pandrosion_step(z1, z0, roots)
    if z2 is None: return z1
    denom = z2 - 2*z1 + z
    if abs(denom) < 1e-50:
        return z2
    return z - (z1 - z)**2 / denom

def run_pandrosion_T3(z_init, roots, K=3, n_steps=300, use_steffensen=True):
    """
    Pandrosion-T3 with anchor update every K steps.
    K=1: Newton-like (anchor = iterate)
    K=3: base Pandrosion-T3
    K=∞: fixed anchor
    """
    d = len(roots)
    rho = np.max(np.abs(roots))
    R = 1 + rho
    
    z = z_init
    z0 = z_init  # initial anchor = starting point
    log_Ps = [eval_P_log(z, roots)]
    step_count = 0
    
    for epoch in range(n_steps):
        # Within epoch: run K Pandrosion base-map steps with fixed z0
        for k_step in range(K):
            if use_steffensen and K >= 3:
                z_next = steffensen_step(z, z0, roots)
            else:
                z_next = pandrosion_step(z, z0, roots)
            
            if z_next is None or np.abs(z_next) > 100 * R or np.isnan(z_next):
                return log_Ps, step_count, True
            
            z = z_next
            step_count += 1
            log_Ps.append(eval_P_log(z, roots))
            
            if np.min(np.abs(z - roots)) < 1e-12:
                return log_Ps, step_count, True
        
        # Update anchor
        z0 = z
    
    return log_Ps, step_count, False


def run_fixed_anchor(z_init, z0_fixed, roots, n_steps=500):
    """Base map with fixed anchor (no reanchoring)."""
    z = z_init
    d = len(roots)
    rho = np.max(np.abs(roots))
    R = 1 + rho
    log_Ps = [eval_P_log(z, roots)]
    
    for n in range(n_steps):
        z_next = pandrosion_step(z, z0_fixed, roots)
        if z_next is None or np.abs(z_next) > 100 * R or np.isnan(z_next):
            break
        z = z_next
        log_Ps.append(eval_P_log(z, roots))
        if np.min(np.abs(z - roots)) < 1e-12:
            break
    
    return log_Ps


# ═══════════════════════════════════════════════════════════════════
# Main experiment: compare K=1, K=3, K=∞ for all families
# ═══════════════════════════════════════════════════════════════════

def main_experiment():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  PANDROSION-T3 SCALING-AWARE DESCENT VERIFICATION                  ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    
    block_size = 3
    
    for d in [10, 20, 50, 100, 200, 500]:
        families = [
            ("Unity", np.exp(2j * np.pi * np.arange(d) / d)),
        ]
        if d <= 200:
            families.append(("Line", np.array([k/d for k in range(d)], dtype=complex)))
        if d <= 50:
            families.append(("Wilkinson", np.array([k+1 for k in range(d)], dtype=complex)))
        # Random circle
        rng = np.random.RandomState(42)
        families.append(("Circle", np.exp(1j * np.sort(rng.uniform(0, 2*np.pi, d)))))
        
        print(f"\n{'═'*72}")
        print(f"  d = {d}")
        print(f"{'═'*72}")
        print(f"  {'Family':>12s}  {'K':>5s}  {'T3':>4s}  {'steps':>7s}  "
              f"{'conv%':>6s}  {'mean_blk':>10s}  {'c (×d)':>8s}  "
              f"{'frac<0':>7s}  {'viol':>5s}")
        print("  " + "─"*78)
        
        for fname, roots in families:
            rho = np.max(np.abs(roots))
            R = 1 + rho
            
            for K, use_T3, label in [(1, False, "1"), (3, True, "3"), (3, False, "3"), (10, False, "10")]:
                # Skip T3 for K=1 (Newton is already order 2)
                if K == 1 and use_T3:
                    continue
                if K == 3 and not use_T3:
                    K_label = "3"
                    T3_label = "no"
                elif K == 3 and use_T3:
                    K_label = "3"
                    T3_label = "T3"
                else:
                    K_label = str(K)
                    T3_label = "no"
                
                if K == 1:
                    T3_label = "N"  # Newton
                
                all_block_desc = []
                conv = 0
                total_steps = 0
                total_orbits = 0
                
                for s in range(min(30, 2*d)):
                    theta = 2 * np.pi * s / min(30, 2*d)
                    z_init = R * np.exp(1j * theta)
                    total_orbits += 1
                    
                    log_Ps, steps, converged = run_pandrosion_T3(
                        z_init, roots, K=K, n_steps=300, use_steffensen=use_T3)
                    
                    if converged:
                        conv += 1
                    total_steps += steps
                    
                    # Block descents
                    for i in range(len(log_Ps) - block_size):
                        bd = log_Ps[i + block_size] - log_Ps[i]
                        if np.isfinite(bd) and abs(bd) < 1000:
                            all_block_desc.append(bd)
                
                if all_block_desc:
                    arr = np.array(all_block_desc)
                    mean_bd = np.mean(arr)
                    c_eff = -mean_bd * d
                    frac_neg = np.mean(arr < 0)
                    violations = np.sum(arr > 0)
                    conv_pct = f"{100*conv/total_orbits:.0f}%"
                    
                    print(f"  {fname:>12s}  {K_label:>5s}  {T3_label:>4s}  "
                          f"{total_steps:>7d}  {conv_pct:>6s}  "
                          f"{mean_bd:>10.4f}  {c_eff:>8.2f}  "
                          f"{frac_neg:>7.3f}  {violations:>5d}")
                else:
                    print(f"  {fname:>12s}  {K_label:>5s}  {T3_label:>4s}  "
                          f"{'---':>7s}")

    # ═══════════════════════════════════════════════════════════════════
    # Scaling analysis: c as function of d for K=1 (Newton)
    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'█'*72}")
    print(f"  SCALING ANALYSIS: c vs d for Newton (K=1)")
    print(f"{'█'*72}")
    print(f"\n  {'d':>5s}  {'Family':>12s}  {'mean_desc':>12s}  {'c':>8s}  "
          f"{'c/d':>10s}  {'c/log(d)':>10s}  {'steps':>7s}  {'conv%':>6s}")
    print("  " + "─"*80)
    
    results = []
    for d in [5, 10, 20, 50, 100, 200, 500]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        rho = 1.0
        R = 2.0
        
        all_lr = []
        conv = 0
        total = 0
        total_steps = 0
        
        for s in range(min(30, 2*d)):
            theta = 2 * np.pi * s / min(30, 2*d)
            z_init = R * np.exp(1j * theta)
            total += 1
            
            log_Ps, steps, converged = run_pandrosion_T3(
                z_init, roots, K=1, n_steps=500, use_steffensen=False)
            
            if converged:
                conv += 1
            total_steps += steps
            
            for i in range(len(log_Ps) - 1):
                lr = log_Ps[i+1] - log_Ps[i]
                if np.isfinite(lr) and abs(lr) < 1000:
                    all_lr.append(lr)
        
        if all_lr:
            arr = np.array(all_lr)
            mean_desc = np.mean(arr)
            c_eff = -mean_desc * d
            c_over_d = c_eff / d
            c_over_logd = c_eff / np.log(d)
            conv_pct = f"{100*conv/total:.0f}%"
            results.append((d, c_eff))
            
            print(f"  {d:>5d}  {'Unity':>12s}  {mean_desc:>12.6f}  {c_eff:>8.2f}  "
                  f"{c_over_d:>10.4f}  {c_over_logd:>10.4f}  "
                  f"{total_steps:>7d}  {conv_pct:>6s}")
    
    # Power-law fit
    if len(results) >= 3:
        ds = np.array([r[0] for r in results])
        cs = np.array([r[1] for r in results])
        valid = cs > 0
        if np.sum(valid) >= 3:
            log_d = np.log(ds[valid])
            log_c = np.log(cs[valid])
            coeffs = np.polyfit(log_d, log_c, 1)
            print(f"\n  Power-law fit: c ~ d^{{{coeffs[0]:.3f}}}")
            print(f"  => Per-step descent ~ -c/d ~ -d^{{{coeffs[0]-1:.3f}}}")


    # ═══════════════════════════════════════════════════════════════════
    # The scaling analogy: anchor proximity = scaling reduction
    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'█'*72}")
    print(f"  THE SCALING ANALOGY: anchor proximity = Pandrosion scaling")
    print(f"{'█'*72}")
    print(f"\n  For P(z) = z^d - 1, the 'effective ratio' is x' = P(z)/P(z0).")
    print(f"  Closer anchor → smaller x' → smaller λ → faster convergence.")
    xp_hdr = "|x'|"
    print(f"\n  {'d':>4s}  {'|z-z0|':>8s}  {xp_hdr:>10s}  {'mean_desc':>12s}  {'c':>8s}")
    print("  " + "─"*50)
    
    for d in [10, 50]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        
        # Start near root 1
        z_near = 1.0 + 0.1/d
        
        for anchor_dist in [0.01, 0.1, 0.5, 1.0, 2.0]:
            z0 = z_near + anchor_dist * np.exp(1j * np.pi/3)
            
            P_z = eval_P(z_near, roots)
            P_z0 = eval_P(z0, roots)
            x_prime = abs(P_z / P_z0)
            
            # Run a few base-map steps
            z = z_near
            log_ratios = []
            for n in range(20):
                z_next = pandrosion_step(z, z0, roots)
                if z_next is None or np.isnan(z_next) or abs(z_next) > 100:
                    break
                lr = eval_P_log(z_next, roots) - eval_P_log(z, roots)
                if np.isfinite(lr):
                    log_ratios.append(lr)
                z = z_next
            
            if log_ratios:
                mean_lr = np.mean(log_ratios)
                c_eff = -mean_lr * d
                print(f"  {d:>4d}  {anchor_dist:>8.3f}  {x_prime:>10.2e}  "
                      f"{mean_lr:>12.6f}  {c_eff:>8.4f}")


if __name__ == "__main__":
    main_experiment()
    print(f"\n{'═'*72}")
    print("  DONE.")
