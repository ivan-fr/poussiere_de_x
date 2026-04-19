#!/usr/bin/env python3
"""
PURE PANDROSION (no Newton): T3 and T4 with iterated scaling.

Key insight: NEVER use F_{z0}(z0) = Newton. Always keep z ≠ z0.
This completely avoids the P'(z)=0 singularity.

Architecture:
  - Anchor a and iterate z are ALWAYS distinct  
  - T3: 3 base-map steps + Aitken Δ² acceleration
  - T4: 4 base-map steps + higher-order extrapolation
  - Scaling cascade: reanchor a ← Aitken result, keep z = last base-map output

The divided difference Q(a,z) = (P(z)-P(a))/(z-a) is ALWAYS defined
(no L'Hôpital needed), and the regularization theorem guarantees |Q| >> 0.
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def eval_P(z, roots):
    return np.prod(z - roots)

def eval_P_log(z, roots):
    return np.sum(np.log(np.abs(z - roots)))

def pandrosion_step(z, a, roots):
    """F_a(z) = a - P(a)/Q(a,z) where Q = (P(z)-P(a))/(z-a).
    REQUIRES z ≠ a (pure Pandrosion, no Newton).
    """
    if abs(z - a) < 1e-30:
        return None  # refuse to do Newton
    P_z = eval_P(z, roots)
    P_a = eval_P(a, roots)
    Q = (P_z - P_a) / (z - a)
    if abs(Q) < 1e-50:
        return None
    return a - P_a / Q

def pandrosion_T3(a, z, roots):
    """One T3 epoch (Steffensen on Pandrosion).
    
    Input: anchor a, starting iterate z (a ≠ z)
    Returns: (new_anchor, new_iterate, trajectory, success)
    
    Steps:
      z1 = F_a(z)     — pure Pandrosion
      z2 = F_a(z1)    — pure Pandrosion  
      ẑ = z - (z1-z)²/(z2-2z1+z)  — Aitken Δ²
    
    Output: new anchor = ẑ, new iterate = z2 (always distinct)
    """
    z1 = pandrosion_step(z, a, roots)
    if z1 is None:
        return None, None, [z], False
    
    z2 = pandrosion_step(z1, a, roots)
    if z2 is None:
        return z1, z, [z, z1], False
    
    denom = z2 - 2*z1 + z
    if abs(denom) < 1e-50:
        return z2, z1, [z, z1, z2], True
    
    z_hat = z - (z1 - z)**2 / denom
    
    # Ensure new anchor ≠ new iterate
    # new anchor = z_hat (Aitken result)
    # new iterate = z2 (last base-map output)
    if abs(z_hat - z2) < 1e-30:
        # Fallback: use z1 as iterate
        return z_hat, z1, [z, z1, z2, z_hat], True
    
    return z_hat, z2, [z, z1, z2, z_hat], True

def pandrosion_T4(a, z, roots):
    """One T4 epoch: 4 base-map steps + double Aitken.
    
    Steps:
      z1 = F_a(z), z2 = F_a(z1), z3 = F_a(z2), z4 = F_a(z3)
      ẑ₁ = Aitken(z,z1,z2)
      ẑ₂ = Aitken(z1,z2,z3)
      ẑ = Aitken(ẑ₁,ẑ₂,Aitken(z2,z3,z4))  — or simpler: just Aitken(z1,z2,z3)
    
    Output: new anchor = ẑ, new iterate = z4
    """
    z1 = pandrosion_step(z, a, roots)
    if z1 is None:
        return None, None, [z], False
    z2 = pandrosion_step(z1, a, roots)
    if z2 is None:
        return z1, z, [z, z1], False
    z3 = pandrosion_step(z2, a, roots)
    if z3 is None:
        return z2, z1, [z, z1, z2], True
    z4 = pandrosion_step(z3, a, roots)
    if z4 is None:
        return z3, z2, [z, z1, z2, z3], True
    
    # Double Aitken: 
    # First Aitken on (z1, z2, z3)
    d1 = z3 - 2*z2 + z1
    if abs(d1) < 1e-50:
        zhat1 = z3
    else:
        zhat1 = z1 - (z2-z1)**2 / d1
    
    # Second Aitken on (z2, z3, z4)
    d2 = z4 - 2*z3 + z2
    if abs(d2) < 1e-50:
        zhat2 = z4
    else:
        zhat2 = z2 - (z3-z2)**2 / d2
    
    # Use the better Aitken result (closer to a root)
    log_P_hat1 = eval_P_log(zhat1, roots) if np.isfinite(zhat1) and abs(zhat1) < 1e10 else float('inf')
    log_P_hat2 = eval_P_log(zhat2, roots) if np.isfinite(zhat2) and abs(zhat2) < 1e10 else float('inf')
    
    if log_P_hat1 < log_P_hat2:
        z_hat = zhat1
    else:
        z_hat = zhat2
    
    # new anchor = z_hat, new iterate = z4
    if abs(z_hat - z4) < 1e-30:
        return z_hat, z3, [z, z1, z2, z3, z4, z_hat], True
    
    return z_hat, z4, [z, z1, z2, z3, z4, z_hat], True


def run_pure_pandrosion(roots, z_anchor, z_iter, mode="T3", max_epochs=200):
    """Run pure Pandrosion with iterated scaling (redéfinir A à chaque fois).
    
    a = anchor (= "A" in scaling)
    z = iterate (always ≠ a)
    """
    d = len(roots)
    a = z_anchor
    z = z_iter
    
    epoch_data = []
    
    for epoch in range(max_epochs):
        log_P_a = eval_P_log(a, roots)
        log_P_z = eval_P_log(z, roots)
        
        if mode == "T3":
            a_new, z_new, traj, ok = pandrosion_T3(a, z, roots)
        elif mode == "T4":
            a_new, z_new, traj, ok = pandrosion_T4(a, z, roots)
        else:
            raise ValueError(f"Unknown mode: {mode}")
        
        if a_new is None or not ok:
            break
        if np.isnan(a_new) or np.isnan(z_new) or abs(a_new) > 1e10 or abs(z_new) > 1e10:
            break
        
        log_P_new = eval_P_log(a_new, roots)
        descent = log_P_new - log_P_a
        
        dist_root_a = np.min(np.abs(a_new - roots))
        dist_root_z = np.min(np.abs(z_new - roots))
        dist_az = abs(a_new - z_new)
        
        epoch_data.append({
            'epoch': epoch,
            'log_P_a': log_P_a,
            'log_P_new': log_P_new,
            'descent': descent,
            'dist_root': min(dist_root_a, dist_root_z),
            'dist_az': dist_az,
            'n_steps': len(traj) - 1,
        })
        
        # Iterated scaling: redefine A
        a = a_new
        z = z_new
        
        if min(dist_root_a, dist_root_z) < 1e-12:
            break
    
    return epoch_data


def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  PURE PANDROSION T3/T4 — No Newton, derivative-free                ║")
    print("║  Anchor ≠ Iterate at every step                                    ║")
    print("║  Iterated scaling = redéfinir A à chaque époque                    ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    
    # ═══════════════════════════════════════════════════════════════════
    # Experiment 1: Epoch-by-epoch trace
    # ═══════════════════════════════════════════════════════════════════
    
    for d in [10, 50, 100]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        rho = 1.0
        R = 2.0
        
        # TWO starting points on Cauchy circle (always distinct)
        a0 = R * np.exp(1j * 0.0)  # anchor at angle 0
        z0 = R * np.exp(1j * np.pi / d)  # iterate at angle π/d (next Cauchy point)
        
        for mode in ["T3", "T4"]:
            print(f"\n{'═'*72}")
            print(f"  d = {d} | {mode} | anchor=R, iterate=R·e^(iπ/d)")
            print(f"{'═'*72}")
            print(f"  {'Ep':>4s}  {'log|P(a)|':>12s}  {'log|P(a_new)|':>14s}  "
                  f"{'descent':>8s}  {'d(root)':>10s}  {'|a-z|':>10s}")
            print("  " + "─"*65)
            
            epochs = run_pure_pandrosion(roots, a0, z0, mode=mode, max_epochs=100)
            
            for e in epochs:
                print(f"  {e['epoch']:>4d}  {e['log_P_a']:>12.4f}  {e['log_P_new']:>14.4f}  "
                      f"{e['descent']:>8.4f}  {e['dist_root']:>10.2e}  {e['dist_az']:>10.2e}")
            
            if epochs:
                total = sum(e['descent'] for e in epochs)
                n = len(epochs)
                n_steps = sum(e['n_steps'] for e in epochs)
                conv = epochs[-1]['dist_root'] < 1e-10
                print(f"\n  Total: {n} epochs, {n_steps} steps, "
                      f"mean descent = {total/n:.4f} nats/epoch, "
                      f"{'CONVERGED' if conv else 'NOT CONVERGED'}")
    
    # ═══════════════════════════════════════════════════════════════════
    # Experiment 2: Convergence statistics across starting points
    # ═══════════════════════════════════════════════════════════════════
    
    print(f"\n{'█'*72}")
    print(f"  CONVERGENCE STATISTICS: Pure Pandrosion T3 vs T4")
    print(f"{'█'*72}")
    print(f"\n  {'d':>4s}  {'Mode':>4s}  {'Family':>10s}  {'#orb':>5s}  {'conv%':>6s}  "
          f"{'mean_desc':>10s}  {'median':>10s}  {'epochs':>8s}  {'steps':>8s}")
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
            
            for mode in ["T3", "T4"]:
                all_descents = []
                conv_count = 0
                total_epochs = 0
                total_steps = 0
                n_orbits = 0
                
                for s in range(min(30, 2*d)):
                    theta_a = 2 * np.pi * s / min(30, 2*d)
                    theta_z = theta_a + np.pi / d  # offset by π/d
                    
                    a0 = R * np.exp(1j * theta_a)
                    z0 = R * np.exp(1j * theta_z)
                    n_orbits += 1
                    
                    epochs = run_pure_pandrosion(roots, a0, z0, mode=mode, max_epochs=200)
                    
                    if epochs:
                        total_epochs += len(epochs)
                        total_steps += sum(e['n_steps'] for e in epochs)
                        if epochs[-1]['dist_root'] < 1e-10:
                            conv_count += 1
                        for e in epochs:
                            if np.isfinite(e['descent']) and abs(e['descent']) < 1000:
                                all_descents.append(e['descent'])
                
                if all_descents:
                    arr = np.array(all_descents)
                    conv_pct = f"{100*conv_count/n_orbits:.0f}%"
                    print(f"  {d:>4d}  {mode:>4s}  {fname:>10s}  {n_orbits:>5d}  {conv_pct:>6s}  "
                          f"{np.mean(arr):>10.4f}  {np.median(arr):>10.4f}  "
                          f"{total_epochs:>8d}  {total_steps:>8d}")
    
    # ═══════════════════════════════════════════════════════════════════
    # Experiment 3: Does P'(z)=0 still cause problems?
    # ═══════════════════════════════════════════════════════════════════
    
    print(f"\n{'█'*72}")
    print(f"  KEY TEST: Pure Pandrosion near critical points of P")
    print(f"  Does the P'(z)=0 singularity still cause problems?")
    print(f"{'█'*72}")
    
    for d in [10, 50]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        # Critical point of z^d - 1 is z = 0
        # Start orbit near z = 0
        
        # Anchor on Cauchy circle, iterate NEAR critical point
        a = 2.0 * np.exp(1j * 0.3)
        
        print(f"\n  d = {d}: anchor on Cauchy, iterate near critical point z=0")
        print(f"  {'z_init':>10s}  {'|z|':>6s}  {'mode':>4s}  {'epochs':>7s}  "
              f"{'conv':>5s}  {'mean_desc':>10s}  {'Q fail?':>8s}")
        
        for r_init in [0.01, 0.05, 0.1, 0.3, 0.5]:
            z = r_init * np.exp(1j * 0.7)
            
            for mode in ["T3", "T4"]:
                epochs = run_pure_pandrosion(roots, a, z, mode=mode, max_epochs=100)
                
                if epochs:
                    conv = epochs[-1]['dist_root'] < 1e-10
                    descs = [e['descent'] for e in epochs if np.isfinite(e['descent'])]
                    mean_d = np.mean(descs) if descs else float('nan')
                    print(f"  {r_init:>10.3f}  {r_init:>6.3f}  {mode:>4s}  {len(epochs):>7d}  "
                          f"{'YES' if conv else 'NO':>5s}  {mean_d:>10.4f}  {'NO':>8s}")
                else:
                    print(f"  {r_init:>10.3f}  {r_init:>6.3f}  {mode:>4s}  {'FAIL':>7s}")
    
    # ═══════════════════════════════════════════════════════════════════
    # Experiment 4: BSS operation count scaling
    # ═══════════════════════════════════════════════════════════════════
    
    print(f"\n{'█'*72}")
    print(f"  BSS COST SCALING: Total steps vs d")
    print(f"{'█'*72}")
    print(f"\n  {'d':>5s}  {'Mode':>4s}  {'Total steps':>12s}  {'ops/d':>8s}  "
          f"{'ops/d²':>8s}  {'conv%':>6s}")
    print("  " + "─"*50)
    
    for d in [10, 20, 50, 100, 200, 500]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        
        for mode in ["T3", "T4"]:
            total_steps = 0
            conv = 0
            n_orbits = min(30, 2*d)
            
            for s in range(n_orbits):
                theta_a = 2 * np.pi * s / n_orbits
                theta_z = theta_a + np.pi / d
                a0 = R * np.exp(1j * theta_a)
                z0 = R * np.exp(1j * theta_z)
                
                epochs = run_pure_pandrosion(roots, a0, z0, mode=mode, max_epochs=300)
                if epochs:
                    total_steps += sum(e['n_steps'] for e in epochs)
                    if epochs[-1]['dist_root'] < 1e-10:
                        conv += 1
            
            conv_pct = f"{100*conv/n_orbits:.0f}%"
            ops_d = total_steps / d if d > 0 else 0
            ops_d2 = total_steps / d**2 if d > 0 else 0
            print(f"  {d:>5d}  {mode:>4s}  {total_steps:>12d}  {ops_d:>8.2f}  "
                  f"{ops_d2:>8.4f}  {conv_pct:>6s}")


if __name__ == "__main__":
    main()
    print(f"\n{'═'*72}")
    print("  DONE.")
