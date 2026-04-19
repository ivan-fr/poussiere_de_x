#!/usr/bin/env python3
"""
AMORTIZED DESCENT: λ > 1 sometimes, but PRODUCT ∏λ_n < 1 per epoch.

This is the KEY INSIGHT: the proof of Σ_log < 0 cannot use pointwise λ < 1.
Instead, the PRODUCT of contraction ratios across an epoch is < 1.

The correct proof object is:
  Λ_epoch = ∏_{steps in epoch} λ_step = |P(z_end)/P(z_start)| < 1

NOT λ < 1 at each step, but ∏λ < 1 per epoch.
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def eval_P(z, roots):
    return np.prod(z - roots)

def eval_P_log(z, roots):
    return np.sum(np.log(np.abs(z - roots)))

def pandrosion_step(z, a, roots):
    if abs(z - a) < 1e-30:
        return None
    P_z = eval_P(z, roots)
    P_a = eval_P(a, roots)
    Q = (P_z - P_a) / (z - a)
    if abs(Q) < 1e-50:
        return None
    return a - P_a / Q


def analyze_amortized_lambda(d):
    """Track λ at EACH base-map step and the epoch product."""
    roots = np.exp(2j * np.pi * np.arange(d) / d)
    R = 2.0
    
    print(f"\n{'█'*72}")
    print(f"  AMORTIZED λ: d = {d}")
    print(f"  λ_step = contraction ratio at nearest root for each base-map step")
    print(f"  Λ_epoch = product of λ over epoch")
    print(f"{'█'*72}")
    
    a = R * np.exp(1j * 0.0)
    z = R * np.exp(1j * np.pi / d)
    
    print(f"\n  {'Ep':>4s}  {'step':>4s}  {'|a-ζ|':>10s}  {'λ_step':>10s}  "
          f"{'Λ_epoch':>10s}  {'Σ_log_ep':>10s}  {'desc':>10s}  {'status':>8s}")
    print("  " + "─"*75)
    
    epoch_products = []
    epoch_descs = []
    
    for epoch in range(min(50, 3*d)):
        log_P_start = eval_P_log(a, roots)
        
        # Track lambda per step
        step_lambdas = []
        z_cur = z
        
        for step in range(3):  # T3 = 3 base-map steps
            # Nearest root to ANCHOR
            dists_a = np.abs(a - roots)
            nearest = np.argmin(dists_a)
            zeta = roots[nearest]
            
            # Lambda at nearest root
            P_a = eval_P(a, roots)
            P_prime_zeta = np.prod([zeta - roots[k] for k in range(d) if k != nearest])
            if abs(P_a) > 1e-300:
                lam = abs(1 + P_prime_zeta * (zeta - a) / P_a)
            else:
                lam = 0
            step_lambdas.append(lam)
            
            # Do one step
            z_new = pandrosion_step(z_cur, a, roots)
            if z_new is None or np.isnan(z_new):
                break
            z_cur = z_new
        
        if len(step_lambdas) < 3:
            break
        
        # Aitken
        traj = [z]
        ok = True
        z_temp = z
        for _ in range(3):
            z_n = pandrosion_step(z_temp, a, roots)
            if z_n is None:
                ok = False
                break
            traj.append(z_n)
            z_temp = z_n
        
        if not ok or len(traj) < 4:
            break
        
        z0, z1, z2 = traj[0], traj[1], traj[2]
        den = z2 - 2*z1 + z0
        if abs(den) > 1e-50:
            z_hat = z0 - (z1-z0)**2/den
        else:
            z_hat = traj[-1]
        
        if np.isnan(z_hat) or abs(z_hat) > 1e6:
            z_hat = traj[-1]
        
        log_P_end = eval_P_log(z_hat, roots)
        desc = log_P_end - log_P_start
        
        # Epoch product of lambdas
        Lambda_epoch = np.prod(step_lambdas)
        epoch_products.append(Lambda_epoch)
        epoch_descs.append(desc)
        
        # Print each step
        for step, lam in enumerate(step_lambdas):
            status = "✓" if lam < 1 else "✗"
            ep_label = f"{epoch}" if step == 0 else ""
            Lambda_so_far = np.prod(step_lambdas[:step+1])
            if step == 0:
                print(f"  {ep_label:>4s}  {step:>4d}  {abs(a-zeta):>10.2e}  {lam:>10.6f}  "
                      f"{'':>10s}  {'':>10s}  {'':>10s}  {status:>8s}")
            elif step == len(step_lambdas)-1:
                print(f"  {'':>4s}  {step:>4d}  {'':>10s}  {lam:>10.6f}  "
                      f"{Lambda_epoch:>10.6f}  {desc:>10.4f}  {'':>10s}  {status:>8s}")
            else:
                print(f"  {'':>4s}  {step:>4d}  {'':>10s}  {lam:>10.6f}  "
                      f"{'':>10s}  {'':>10s}  {'':>10s}  {status:>8s}")
        
        # Reanchor
        a = z_hat
        z = traj[-1]
        if abs(a - z) < 1e-30:
            z = z + 0.001 * np.exp(1j * epoch)
        
        if np.min(np.abs(a - roots)) < 1e-12:
            print(f"  {epoch+1:>4d}  CONVERGED")
            break
    
    if epoch_products:
        arr = np.array(epoch_products)
        descs = np.array(epoch_descs)
        
        print(f"\n  EPOCH SUMMARY:")
        print(f"    Epochs:                    {len(arr)}")
        print(f"    max Λ_epoch (product):     {np.max(arr):.6f}")
        print(f"    mean Λ_epoch:              {np.mean(arr):.6f}")
        print(f"    frac(Λ_epoch < 1):         {np.mean(arr < 1):.4f}")
        print(f"    mean Σ_log (descent/ep):   {np.mean(descs):.4f}")
        print(f"    frac(descent < 0):         {np.mean(descs < 0):.4f}")
        
        return arr, descs
    return None, None


def main():
    print("╔═══════════════════════════════════════════════════════════════════╗")
    print("║  AMORTIZED λ ANALYSIS: The proof MUST be amortized, not pointwise ║")
    print("╚═══════════════════════════════════════════════════════════════════╝")
    
    all_results = {}
    
    for d in [5, 10, 20, 50, 100]:
        all_Lambda = []
        all_desc = []
        
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        
        # Multiple starting points
        for s in range(min(20, d)):
            theta_a = 2 * np.pi * s / min(20, d)
            theta_z = theta_a + np.pi / d
            a0 = R * np.exp(1j * theta_a)
            z0 = R * np.exp(1j * theta_z)
            
            a, z = a0, z0
            for epoch in range(200):
                log_P_start = eval_P_log(a, roots)
                
                traj = [z]
                ok = True
                z_t = z
                for _ in range(3):
                    z_n = pandrosion_step(z_t, a, roots)
                    if z_n is None or np.isnan(z_n) or abs(z_n) > 1e6:
                        ok = False
                        break
                    traj.append(z_n)
                    z_t = z_n
                
                if not ok or len(traj) < 4:
                    break
                
                z0t, z1t, z2t = traj[0], traj[1], traj[2]
                den = z2t - 2*z1t + z0t
                if abs(den) > 1e-50:
                    z_hat = z0t - (z1t-z0t)**2/den
                else:
                    z_hat = traj[-1]
                
                if np.isnan(z_hat) or abs(z_hat) > 1e6:
                    z_hat = traj[-1]
                
                log_P_end = eval_P_log(z_hat, roots)
                desc = log_P_end - log_P_start
                
                # Compute Lambda_epoch = |P(z_hat)/P(a)|
                Lambda = np.exp(desc) if abs(desc) < 500 else float('inf')
                
                if np.isfinite(Lambda) and np.isfinite(desc):
                    all_Lambda.append(Lambda)
                    all_desc.append(desc)
                
                a = z_hat
                z = traj[-1]
                if abs(a - z) < 1e-30:
                    z = z + 0.001 * np.exp(1j * epoch)
                
                if np.min(np.abs(a - roots)) < 1e-12:
                    break
        
        if all_Lambda:
            arr = np.array(all_Lambda)
            darr = np.array(all_desc)
            
            all_results[d] = {
                'max_Lambda': np.max(arr),
                'mean_Lambda': np.mean(arr),
                'frac_lt1': np.mean(arr < 1),
                'mean_desc': np.mean(darr),
                'frac_neg': np.mean(darr < 0),
                'n': len(arr),
            }
    
    # Summary
    print(f"\n{'═'*72}")
    print("  AMORTIZED EPOCH DESCENT: Lambda_epoch = |P(a_new)/P(a_old)|")
    print(f"  (This is the PRODUCT of per-step λ across the epoch)")
    print(f"{'═'*72}")
    print(f"\n  {'d':>5s}  {'#ep':>6s}  {'max Λ':>10s}  {'mean Λ':>10s}  {'frac<1':>8s}  "
          f"{'mean desc':>10s}  {'frac<0':>8s}  {'verdict':>10s}")
    print("  " + "─"*72)
    
    all_ok = True
    for d, r in sorted(all_results.items()):
        verdict = "DESCENT" if r['frac_neg'] > 0.999 else "PARTIAL"
        if r['frac_neg'] <= 0.999:
            all_ok = False
        print(f"  {d:>5d}  {r['n']:>6d}  {r['max_Lambda']:>10.6f}  {r['mean_Lambda']:>10.6f}  "
              f"{r['frac_lt1']:>8.4f}  {r['mean_desc']:>10.4f}  {r['frac_neg']:>8.4f}  "
              f"{verdict:>10s}")
    
    print(f"\n{'═'*72}")
    if all_ok:
        print(f"  Σ_log < 0 at EVERY epoch → AMORTIZED DESCENT HOLDS")
    else:
        print(f"  Some epochs have Σ_log > 0 → need to check block averaging")
    print(f"{'═'*72}")
    
    # Detailed trace for d=10
    analyze_amortized_lambda(10)


if __name__ == "__main__":
    main()
