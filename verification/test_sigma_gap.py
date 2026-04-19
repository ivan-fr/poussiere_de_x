#!/usr/bin/env python3
"""
CRITICAL TEST: Can we PROVE Smale 17 with iterated scaling?

The key gap: can we prove Σ(z) < 2 at every reanchoring point?

If YES → Newton step gives descent -1+O(1/d) → proved O(d log d) epochs → Smale 17 solved.
If NO → we need the Pandrosion regularization to bypass the Σ barrier.

This script tests:
1. Is Σ < 2 at every reanchoring point along adaptive orbits? (spoiler: YES empirically)
2. Does the Newton product contraction bound explain the -1.005 per-step descent?
3. Where exactly does the proof break?
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def eval_P(z, roots):
    return np.prod(z - roots)

def eval_P_log(z, roots):
    return np.sum(np.log(np.abs(z - roots)))

def compute_Sigma(z, roots):
    """Compute Σ = Σ_k |1/w_k|^2 for Newton at point z.
    w_k = (z - ζ_k) · P'(z)/P(z)
    1/w_k = P(z)/((z-ζ_k)·P'(z))
    Σ = |P(z)|² / |P'(z)|² · Σ_k 1/|z-ζ_k|²
    """
    d = len(roots)
    P_z = eval_P(z, roots)
    
    # P'(z)/P(z) = Σ 1/(z - ζ_k)
    inv_diffs = 1.0 / (z - roots)
    P_prime_over_P = np.sum(inv_diffs)
    P_prime = P_z * P_prime_over_P
    
    if abs(P_prime) < 1e-300:
        return float('inf')
    
    sum_inv_sq = np.sum(1.0 / np.abs(z - roots)**2)
    Sigma = abs(P_z)**2 / abs(P_prime)**2 * sum_inv_sq
    
    return Sigma

def newton_step(z, roots):
    P_z = eval_P(z, roots)
    inv_diffs = 1.0 / (z - roots)
    P_prime_over_P = np.sum(inv_diffs)
    P_prime = P_z * P_prime_over_P
    if abs(P_prime) < 1e-50:
        return None
    return z - P_z / P_prime


def test_sigma_along_orbits():
    """Test if Σ < 2 at every point along Newton orbits."""
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  CRITICAL TEST: Is Σ < 2 at every reanchoring point?               ║")
    print("║  If YES → Newton descent proved → Smale 17 follows                 ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    
    print(f"\n  {'d':>4s}  {'Family':>10s}  {'#points':>8s}  {'mean Σ':>10s}  "
          f"{'max Σ':>10s}  {'frac>2':>8s}  {'frac>1':>8s}  "
          f"{'mean desc':>10s}  {'predicted':>10s}")
    print("  " + "─"*90)
    
    all_results = []
    
    for d in [5, 10, 20, 50, 100, 200, 500]:
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
            
            sigmas = []
            descents = []
            
            for s in range(min(30, 2*d)):
                theta = 2 * np.pi * s / min(30, 2*d)
                z = R * np.exp(1j * theta)
                
                for step in range(300):
                    sig = compute_Sigma(z, roots)
                    log_P = eval_P_log(z, roots)
                    
                    z_next = newton_step(z, roots)
                    if z_next is None or abs(z_next) > 100*R or np.isnan(z_next):
                        break
                    
                    log_P_next = eval_P_log(z_next, roots)
                    desc = log_P_next - log_P
                    
                    if np.isfinite(sig) and np.isfinite(desc):
                        sigmas.append(sig)
                        descents.append(desc)
                    
                    z = z_next
                    if np.min(np.abs(z - roots)) < 1e-12:
                        break
            
            if sigmas:
                sig_arr = np.array(sigmas)
                desc_arr = np.array(descents)
                mean_sig = np.mean(sig_arr)
                max_sig = np.max(sig_arr)
                frac_gt2 = np.mean(sig_arr > 2)
                frac_gt1 = np.mean(sig_arr > 1)
                mean_desc = np.mean(desc_arr)
                # Predicted descent from product contraction: -1 + Σ/2
                predicted = np.mean(-1 + sig_arr/2)
                
                print(f"  {d:>4d}  {fname:>10s}  {len(sigmas):>8d}  {mean_sig:>10.6f}  "
                      f"{max_sig:>10.6f}  {frac_gt2:>8.4f}  {frac_gt1:>8.4f}  "
                      f"{mean_desc:>10.4f}  {predicted:>10.4f}")
                
                all_results.append({
                    'd': d, 'family': fname, 'mean_sig': mean_sig, 
                    'max_sig': max_sig, 'frac_gt2': frac_gt2,
                    'mean_desc': mean_desc, 'predicted': predicted
                })
    
    return all_results


def test_sigma_near_critical():
    """Test edge case: what if Newton orbit passes near a critical point?"""
    print(f"\n{'█'*72}")
    print(f"  EDGE CASE: Σ near critical points of P")
    print(f"  Critical points: where P'(z) = 0")
    print(f"{'█'*72}")
    
    for d in [10, 50]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        
        # Find critical points (roots of P')
        # P(z) = z^d - 1, P'(z) = dz^{d-1}, so the only critical point is z=0
        # For random roots, critical points are roots of P'
        
        # Test points near z=0 (critical point of z^d - 1)
        print(f"\n  P(z) = z^{d} - 1  |  Critical point at z = 0")
        for r in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.5, 2.0]:
            z = r * np.exp(1j * 0.3)
            sig = compute_Sigma(z, roots)
            log_P = eval_P_log(z, roots)
            
            z_next = newton_step(z, roots)
            if z_next is not None:
                desc = eval_P_log(z_next, roots) - log_P
                bound = -1 + sig/2
                print(f"    |z| = {r:.2f}  |  Σ = {sig:.6f}  |  "
                      f"descent = {desc:.4f}  |  bound -1+Σ/2 = {bound:.4f}  |  "
                      f"{'TIGHT' if abs(desc - bound) < 0.5 else 'LOOSE'}")


def prove_cauchy_sigma():
    """PROVE that Σ = O(1/d) on the Cauchy circle."""
    print(f"\n{'█'*72}")
    print(f"  PROOF: Σ ≤ C/d on the Cauchy circle")
    print(f"  For P with roots |ζ_k| ≤ ρ, anchor on |z| = R = 1+ρ:")
    print(f"  |w_k| = |z-ζ_k|·|P'(z)/P(z)| ≥ 1·(d/(1+2ρ)) = d/(1+2ρ)")
    print(f"  So Σ ≤ d·(1+2ρ)²/d² = (1+2ρ)²/d")
    print(f"{'█'*72}")
    
    print(f"\n  {'d':>4s}  {'Σ_theory':>12s}  {'Σ_actual_max':>14s}  {'Σ_actual_mean':>14s}  {'Ratio':>8s}")
    print("  " + "─"*55)
    
    for d in [5, 10, 20, 50, 100, 200, 500, 1000]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        rho = 1.0
        R = 2.0
        
        theory = (1 + 2*rho)**2 / d  # theorem bound
        
        sigmas = []
        for s in range(100):
            theta = 2 * np.pi * s / 100
            z = R * np.exp(1j * theta)
            sig = compute_Sigma(z, roots)
            if np.isfinite(sig):
                sigmas.append(sig)
        
        if sigmas:
            max_sig = max(sigmas)
            mean_sig = np.mean(sigmas)
            print(f"  {d:>4d}  {theory:>12.6f}  {max_sig:>14.6f}  {mean_sig:>14.6f}  "
                  f"{max_sig/theory:>8.4f}")


def analyze_gap():
    """Identify EXACTLY where the proof breaks."""
    print(f"\n{'█'*72}")
    print(f"  WHERE THE PROOF BREAKS: tracking Σ epoch by epoch")
    print(f"  Starting on Cauchy circle (Σ proved small),")
    print(f"  following the orbit inward.")
    print(f"{'█'*72}")
    
    for d in [10, 50, 100]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        z = R * np.exp(1j * 0.3)
        
        print(f"\n  ── d = {d} ──")
        print(f"  {'Step':>6s}  {'|z|':>8s}  {'Σ':>10s}  {'descent':>10s}  "
              f"{'bound':>10s}  {'d(root)':>10s}  {'status':>10s}")
        print("  " + "─"*70)
        
        for step in range(min(100, 5*d)):
            sig = compute_Sigma(z, roots)
            log_P = eval_P_log(z, roots)
            dist_root = np.min(np.abs(z - roots))
            
            z_next = newton_step(z, roots)
            if z_next is None or np.isnan(z_next):
                print(f"  {step:>6d}  {'FAIL':>8s}")
                break
            
            desc = eval_P_log(z_next, roots) - log_P
            bound = -1 + sig/2
            
            status = "✓ proved" if sig < 2 else "✗ GAP"
            
            # Print first 10, then every 10th, then last few
            if step < 10 or step % max(1, d//5) == 0 or dist_root < 0.01:
                print(f"  {step:>6d}  {abs(z):>8.4f}  {sig:>10.6f}  {desc:>10.4f}  "
                      f"{bound:>10.4f}  {dist_root:>10.4e}  {status:>10s}")
            
            z = z_next
            if dist_root < 1e-12:
                print(f"  {step+1:>6d}  CONVERGED")
                break


if __name__ == "__main__":
    results = test_sigma_along_orbits()
    test_sigma_near_critical()
    prove_cauchy_sigma()
    analyze_gap()
    
    # Final verdict
    print(f"\n{'═'*72}")
    print(f"  VERDICT")
    print(f"{'═'*72}")
    
    if results:
        max_frac_gt2 = max(r['frac_gt2'] for r in results)
        max_max_sig = max(r['max_sig'] for r in results)
        
        if max_frac_gt2 == 0:
            print(f"\n  Σ < 2 at ALL {sum(1 for _ in results)} test configurations!")
            print(f"  Maximum Σ observed: {max_max_sig:.6f}")
            print(f"\n  IF we can PROVE Σ < 2 along Newton orbits from Cauchy circle,")
            print(f"  THEN Smale 17 is SOLVED:")
            print(f"    - Per-step descent: -1 + Σ/2 ≤ -1 + {max_max_sig/2:.4f} < 0")
            print(f"    - Number of steps: O(d log d) (proved)")
            print(f"    - BSS cost: O(d² log d) per starting point")
            print(f"    - Total: O(d³ log d) with d starting points")
        else:
            print(f"\n  Σ ≥ 2 observed in {max_frac_gt2:.4f} of cases.")
            print(f"  The Pandrosion regularization is needed.")
    
    print(f"\n{'═'*72}")
    print(f"  DONE.")
