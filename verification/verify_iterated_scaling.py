#!/usr/bin/env python3
"""
Iterated Scaling Verification.

Key idea from pandrosion_en_improved.tex:
  1. Optimal starting point: s_0 = h(1) — one step from naive → best residual
  2. Scaling: x^{1/p} = A^{1/p} · (x/A)^{1/p}, with A redefined at each step

For polynomial root-finding, the adaptive reanchoring IS iterated scaling:
  - Epoch n: anchor z_0^(n) plays role of A (scaling factor)
  - Pandrosion steps compute on "reduced ratio" x' = P(z)/P(z_0) ≈ 1
  - At epoch end: A ← z_K (redefine A = best current approximation)
  - Contraction λ(x') → 0 as x' → 1

This script traces the epoch-by-epoch structure showing:
  - |P(z_0)| (the "A" at each epoch)
  - |P(z_K)/P(z_0)| (the effective ratio x')
  - The resulting contraction
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
        # L'Hôpital: Q(z0,z0) = P'(z0) → Newton step
        P_prime = P_z * np.sum(1.0 / (z - roots))
        if abs(P_prime) < 1e-50:
            return None
        return z - P_z / P_prime
    Q = (P_z - P_z0) / (z - z0)
    if abs(Q) < 1e-50:
        return None
    return z0 - P_z0 / Q

def steffensen_T3(z0, roots):
    """One Steffensen (T3) epoch from anchor z0.
    
    z1 = F_{z0}(z0) = Newton step (= optimal starting point h(1))
    z2 = F_{z0}(z1)
    z3 = Aitken acceleration
    
    Returns (z3, [z0, z1, z2, z3])
    """
    # Step 1: F_{z0}(z0) = Newton (first step = optimal starting point)
    z1 = pandrosion_step(z0, z0, roots)
    if z1 is None:
        return None, [z0]
    
    # Step 2: F_{z0}(z1) with same anchor
    z2 = pandrosion_step(z1, z0, roots)
    if z2 is None:
        return z1, [z0, z1]
    
    # Step 3: Aitken Δ² acceleration
    denom = z2 - 2*z1 + z0
    if abs(denom) < 1e-50:
        return z2, [z0, z1, z2]
    
    z3 = z0 - (z1 - z0)**2 / denom
    return z3, [z0, z1, z2, z3]


def trace_epochs(z_init, roots, max_epochs=30):
    """Trace the iterated scaling epoch by epoch."""
    d = len(roots)
    z0 = z_init
    
    epoch_data = []
    
    for epoch in range(max_epochs):
        log_P_z0 = eval_P_log(z0, roots)
        
        # Run one T3 epoch
        z_new, trajectory = steffensen_T3(z0, roots)
        
        if z_new is None:
            break
        
        log_P_new = eval_P_log(z_new, roots)
        
        # Effective ratio x' = |P(z_K)/P(z_0)|
        log_ratio = log_P_new - log_P_z0
        x_prime = np.exp(log_ratio)
        
        # Distance to nearest root
        dist_root = np.min(np.abs(z_new - roots))
        
        # Kinematic ratio r = P(z1)/P(z0) after first step
        z1 = trajectory[1] if len(trajectory) > 1 else z0
        r_first = abs(eval_P(z1, roots) / eval_P(z0, roots)) if abs(eval_P(z0, roots)) > 1e-300 else 0
        
        epoch_data.append({
            'epoch': epoch,
            'log_P_z0': log_P_z0,
            'log_P_new': log_P_new,
            'log_ratio': log_ratio,
            'x_prime': x_prime,
            'dist_root': dist_root,
            'r_first': r_first,
            'z0': z0,
            'z_new': z_new,
            'n_steps': len(trajectory) - 1,
        })
        
        # Update anchor (= redefine A)
        z0 = z_new
        
        if dist_root < 1e-14:
            break
    
    return epoch_data


def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  ITERATED SCALING: Adaptive Reanchoring = Redéfinir A à chaque pas ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    
    # ═══════════════════════════════════════════════════════════════════
    # Experiment 1: Epoch-by-epoch trace for specific examples
    # ═══════════════════════════════════════════════════════════════════
    
    for d in [10, 50, 100]:
        families = [
            ("Unity", np.exp(2j * np.pi * np.arange(d) / d)),
        ]
        if d <= 50:
            rng = np.random.RandomState(123)
            families.append(("Circle", np.exp(1j * np.sort(rng.uniform(0, 2*np.pi, d)))))
        
        for fname, roots in families:
            rho = np.max(np.abs(roots))
            R = 1 + rho
            
            # Start from Cauchy circle
            z_init = R * np.exp(1j * 0.3)
            
            print(f"\n{'═'*72}")
            print(f"  {fname} (d = {d})  |  z_init on Cauchy circle (R = {R:.2f})")
            print(f"{'═'*72}")
            xp_h = "|x'|"
            print(f"  {'Epoch':>5s}  {'log|P(z0)|':>12s}  {'log|P(zK)|':>12s}  "
                  f"{'Descent':>8s}  {xp_h:>10s}  {'|r1|':>10s}  {'d(root)':>10s}")
            print("  " + "─"*75)
            
            epochs = trace_epochs(z_init, roots, max_epochs=30)
            
            for e in epochs:
                xp = f"{e['x_prime']:.2e}" if e['x_prime'] < 1e6 else "large"
                print(f"  {e['epoch']:>5d}  {e['log_P_z0']:>12.4f}  {e['log_P_new']:>12.4f}  "
                      f"{e['log_ratio']:>8.2f}  {xp:>10s}  {e['r_first']:>10.2e}  "
                      f"{e['dist_root']:>10.2e}")
    
    # ═══════════════════════════════════════════════════════════════════
    # Experiment 2: Connect to Pandrosion contraction ratio formula
    # ═══════════════════════════════════════════════════════════════════
    
    print(f"\n{'█'*72}")
    print(f"  CONTRACTION RATIO ANALOGY")
    print(f"  For P(z) = z^p - 1, the Pandrosion iteration on the reduced")
    print(f"  problem x' has contraction λ(x') ≈ (p-1)(α'-1)/2 as x'→1.")
    print(f"  In the polynomial case, x' = P(z_K)/P(z_0) plays the same role.")
    print(f"{'█'*72}")
    
    xp_h2 = "|x'| = |r|"
    lam_h = "lam_predict"
    print(f"\n  {'d':>4s}  {'Epoch':>5s}  {xp_h2:>14s}  "
          f"{lam_h:>14s}  {'Total desc.':>12s}  {'Steps':>10s}")
    print("  " + "─"*65)
    
    for d in [10, 50, 100, 200]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        z_init = R * np.exp(1j * 0.3)
        
        epochs = trace_epochs(z_init, roots, max_epochs=50)
        
        total_desc = 0
        for e in epochs[:5]:  # Show first 5 epochs
            r = e['x_prime']
            # Pandrosion-like contraction prediction
            if r > 0 and r < 10:
                alpha_prime = r ** (1.0/d) if r > 0 else 0
                # λ ≈ (d-1)(α'-1)/2, but adjusted for polynomial case
                lambda_pred = (d-1) * abs(alpha_prime - 1) / 2 if abs(alpha_prime - 1) < 1 else 1.0
            else:
                lambda_pred = float('nan')
            total_desc += e['log_ratio']
            print(f"  {d:>4d}  {e['epoch']:>5d}  {r:>14.6e}  "
                  f"{lambda_pred:>14.6f}  {total_desc:>12.4f}  "
                  f"{e['dist_root']:>10.2e}")
        
        if epochs:
            n_epochs = len(epochs)
            total_steps = sum(e['n_steps'] for e in epochs)
            final_dist = epochs[-1]['dist_root']
            print(f"  {d:>4d}  {'Total':>5s}  {'':>14s}  {'':>14s}  "
                  f"{sum(e['log_ratio'] for e in epochs):>12.4f}  "
                  f"{'conv' if final_dist < 1e-10 else 'no':>10s}  "
                  f"({n_epochs} epochs, {total_steps} steps)")
    
    # ═══════════════════════════════════════════════════════════════════
    # Experiment 3: Per-epoch descent is O(1) — the key claim
    # ═══════════════════════════════════════════════════════════════════
    
    print(f"\n{'█'*72}")
    print(f"  PER-EPOCH DESCENT STATISTICS (T3 with iterated scaling)")
    print(f"{'█'*72}")
    print(f"\n  {'d':>5s}  {'Family':>10s}  {'#orbits':>8s}  {'#epochs':>8s}  "
          f"{'mean desc':>10s}  {'median':>10s}  {'min':>10s}  {'frac<0':>7s}")
    print("  " + "─"*75)
    
    for d in [10, 20, 50, 100, 200, 500]:
        families = [
            ("Unity", np.exp(2j * np.pi * np.arange(d) / d)),
        ]
        rng = np.random.RandomState(42)
        families.append(("Circle", np.exp(1j * np.sort(rng.uniform(0, 2*np.pi, d)))))
        if d <= 100:
            families.append(("Line", np.array([k/d for k in range(d)], dtype=complex)))
        
        for fname, roots in families:
            rho = np.max(np.abs(roots))
            R = 1 + rho
            
            all_descents = []
            n_orbits = 0
            total_epochs = 0
            
            for s in range(min(30, 2*d)):
                theta = 2 * np.pi * s / min(30, 2*d)
                z_init = R * np.exp(1j * theta)
                n_orbits += 1
                
                epochs = trace_epochs(z_init, roots, max_epochs=100)
                total_epochs += len(epochs)
                
                for e in epochs:
                    if np.isfinite(e['log_ratio']) and abs(e['log_ratio']) < 1000:
                        all_descents.append(e['log_ratio'])
            
            if all_descents:
                arr = np.array(all_descents)
                print(f"  {d:>5d}  {fname:>10s}  {n_orbits:>8d}  {total_epochs:>8d}  "
                      f"{np.mean(arr):>10.4f}  {np.median(arr):>10.4f}  "
                      f"{np.min(arr):>10.4f}  {np.mean(arr < 0):>7.3f}")
    
    # ═══════════════════════════════════════════════════════════════════
    # Experiment 4: "Optimal starting point" effect
    # Show that the FIRST step (z0 → Newton(z0)) is the analogue of h(1)
    # ═══════════════════════════════════════════════════════════════════
    
    print(f"\n{'█'*72}")
    print(f"  OPTIMAL STARTING POINT: First step = Newton = h(1)")
    print(f"  The first step of each epoch is F_{{z0}}(z0) = N_P(z0)")
    print(f"  This is the analogue of s_0^opt = h(1) — best initial residual.")
    print(f"{'█'*72}")
    
    print(f"\n  {'d':>4s}  {'Step':>6s}  {'log|P| drop':>14s}  {'Fraction of total':>18s}")
    print("  " + "─"*50)
    
    for d in [10, 50, 100]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        
        step1_drops = []
        step2_drops = []
        step3_drops = []
        total_drops = []
        
        for s in range(30):
            theta = 2 * np.pi * s / 30
            z_init = R * np.exp(1j * theta)
            
            z0 = z_init
            log_P0 = eval_P_log(z0, roots)
            
            # Step 1: Newton (= h(1), optimal starting point)
            z1 = pandrosion_step(z0, z0, roots)
            if z1 is None:
                continue
            log_P1 = eval_P_log(z1, roots)
            
            # Step 2: F_{z0}(z1)
            z2 = pandrosion_step(z1, z0, roots)
            if z2 is None:
                continue
            log_P2 = eval_P_log(z2, roots)
            
            # Step 3: F_{z0}(z2)
            z3 = pandrosion_step(z2, z0, roots)
            if z3 is None:
                continue
            log_P3 = eval_P_log(z3, roots)
            
            d1 = log_P1 - log_P0
            d2 = log_P2 - log_P1
            d3 = log_P3 - log_P2
            dt = log_P3 - log_P0
            
            if all(np.isfinite([d1, d2, d3, dt])):
                step1_drops.append(d1)
                step2_drops.append(d2)
                step3_drops.append(d3)
                total_drops.append(dt)
        
        if step1_drops:
            m1 = np.mean(step1_drops)
            m2 = np.mean(step2_drops)
            m3 = np.mean(step3_drops)
            mt = np.mean(total_drops)
            
            print(f"  {d:>4d}  {'1(N)':>6s}  {m1:>14.4f}  {m1/mt*100:>17.1f}%")
            print(f"  {d:>4d}  {'2(F)':>6s}  {m2:>14.4f}  {m2/mt*100:>17.1f}%")
            print(f"  {d:>4d}  {'3(F)':>6s}  {m3:>14.4f}  {m3/mt*100:>17.1f}%")
            print(f"  {d:>4d}  {'Total':>6s}  {mt:>14.4f}  {'100.0':>17s}%")
            print()


if __name__ == "__main__":
    main()
    print(f"\n{'═'*72}")
    print("  DONE.")
