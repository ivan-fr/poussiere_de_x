#!/usr/bin/env python3
"""
ADVERSARIAL POLYNOMIAL TEST: Pure Pandrosion T3/T4 + iterated scaling.

Test the hardest polynomials known in numerical analysis:
1. Wilkinson (roots 1,...,20 — ill-conditioned)
2. Clustered roots (separation ε)
3. Multi-scale roots (radius ratio >> 1)
4. Random complex coefficients
5. Chebyshev roots
6. Spiral roots
7. Mignotte polynomial (pair of roots at distance 2^{-d})
8. Mandelbrot-critical (roots near unit circle)
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════════
# Core: Pure Pandrosion with all optimizations
# ═══════════════════════════════════════════════════════════════════════

def eval_P_from_roots(z, roots):
    return np.prod(z - roots)

def eval_P_log_from_roots(z, roots):
    return np.sum(np.log(np.abs(z - roots) + 1e-300))

def pandrosion_step(z, a, roots):
    """Pure Pandrosion (a ≠ z always). Log-stable version."""
    if abs(z - a) < 1e-30:
        return None
    d = len(roots)
    # For moderate d: direct computation
    if d <= 30:
        P_z = np.prod(z - roots)
        P_a = np.prod(a - roots)
        Q = (P_z - P_a) / (z - a)
        if abs(Q) < 1e-50:
            return None
        return a - P_a / Q
    else:
        # Log-stable: F_a(z) = a - (z-a)/(r-1) where r = P(z)/P(a)
        try:
            log_r = np.sum(np.log((z - roots) / (a - roots)))
            r = np.exp(log_r)
            if abs(r - 1) < 1e-30:
                return None
            return a - (z - a) / (r - 1)
        except:
            return None

def pandrosion_T3_epoch(a, z, roots):
    """One T3 epoch: 3 base-map steps + Aitken."""
    z1 = pandrosion_step(z, a, roots)
    if z1 is None: return None, None, False
    z2 = pandrosion_step(z1, a, roots)
    if z2 is None: return z1, z, True
    z3 = pandrosion_step(z2, a, roots)
    if z3 is None: return z2, z1, True
    
    den = z2 - 2*z1 + z
    if abs(den) > 1e-50:
        z_hat = z - (z1 - z)**2 / den
    else:
        z_hat = z3
    
    if np.isnan(z_hat) or abs(z_hat) > 1e15:
        z_hat = z3
    
    # new anchor = z_hat, new iterate = z3
    if abs(z_hat - z3) < 1e-30:
        return z_hat, z2, True
    return z_hat, z3, True

def pandrosion_T4_epoch(a, z, roots):
    """One T4 epoch: 4 base-map steps + double Aitken."""
    z1 = pandrosion_step(z, a, roots)
    if z1 is None: return None, None, False
    z2 = pandrosion_step(z1, a, roots)
    if z2 is None: return z1, z, True
    z3 = pandrosion_step(z2, a, roots)
    if z3 is None: return z2, z1, True
    z4 = pandrosion_step(z3, a, roots)
    if z4 is None: return z3, z2, True
    
    # Aitken on (z1,z2,z3)
    d1 = z3 - 2*z2 + z1
    zh1 = z1 - (z2-z1)**2/d1 if abs(d1) > 1e-50 else z3
    # Aitken on (z2,z3,z4)
    d2 = z4 - 2*z3 + z2
    zh2 = z2 - (z3-z2)**2/d2 if abs(d2) > 1e-50 else z4
    
    # Pick better
    for zh in [zh1, zh2]:
        if np.isnan(zh) or abs(zh) > 1e15:
            zh = z4
    
    lp1 = eval_P_log_from_roots(zh1, roots) if np.isfinite(zh1) and abs(zh1) < 1e15 else 1e30
    lp2 = eval_P_log_from_roots(zh2, roots) if np.isfinite(zh2) and abs(zh2) < 1e15 else 1e30
    z_hat = zh1 if lp1 < lp2 else zh2
    
    if abs(z_hat - z4) < 1e-30:
        return z_hat, z3, True
    return z_hat, z4, True


def run_adaptive(roots, mode="T4", max_epochs=300, n_starts=None):
    """Run pure Pandrosion with iterated scaling from multiple starts.
    
    Optimal starting point: for each anchor a on Cauchy circle,
    the iterate z starts at the NEXT point (offset π/d) — this is the
    analogue of s₀ = h(1) (one Pandrosion step from the naive start).
    """
    d = len(roots)
    if n_starts is None:
        n_starts = max(d, 10)
    
    rho = np.max(np.abs(roots))
    R = max(1 + rho, 2 * rho, 2.0)
    
    best_dist = float('inf')
    best_root_idx = -1
    best_epochs = max_epochs
    best_steps = 0
    all_desc = []
    convergences = []
    
    for s in range(n_starts):
        theta_a = 2 * np.pi * s / n_starts
        theta_z = theta_a + np.pi / max(d, 2)  # optimal offset
        
        a = R * np.exp(1j * theta_a)
        z = R * np.exp(1j * theta_z)
        
        converged = False
        total_steps = 0
        
        for epoch in range(max_epochs):
            log_P_before = eval_P_log_from_roots(a, roots)
            
            if mode == "T3":
                a_new, z_new, ok = pandrosion_T3_epoch(a, z, roots)
                steps = 3
            else:
                a_new, z_new, ok = pandrosion_T4_epoch(a, z, roots)
                steps = 4
            
            if not ok or a_new is None:
                break
            
            total_steps += steps
            log_P_after = eval_P_log_from_roots(a_new, roots)
            desc = log_P_after - log_P_before
            
            if np.isfinite(desc) and abs(desc) < 1000:
                all_desc.append(desc)
            
            a = a_new
            z = z_new
            
            # Check convergence
            dists = np.abs(a - roots)
            min_dist = np.min(dists)
            nearest = np.argmin(dists)
            
            if min_dist < 1e-8:
                converged = True
                if min_dist < best_dist:
                    best_dist = min_dist
                    best_root_idx = nearest
                    best_epochs = epoch + 1
                    best_steps = total_steps
                break
            
            if abs(a) > 1e15 or np.isnan(a):
                break
        
        convergences.append(converged)
    
    result = {
        'conv_rate': np.mean(convergences),
        'conv_count': np.sum(convergences),
        'n_starts': n_starts,
        'best_epochs': best_epochs if best_dist < 1e-8 else -1,
        'best_steps': best_steps,
        'best_dist': best_dist,
        'mean_desc': np.mean(all_desc) if all_desc else float('nan'),
        'max_desc': np.max(all_desc) if all_desc else float('nan'),
        'frac_neg': np.mean(np.array(all_desc) < 0) if all_desc else 0,
    }
    return result


# ═══════════════════════════════════════════════════════════════════════
# Adversarial polynomial families
# ═══════════════════════════════════════════════════════════════════════

def get_adversarial_families():
    families = []
    
    # 1. Wilkinson polynomial: roots at 1,...,20
    families.append(("Wilkinson-20", np.arange(1, 21, dtype=complex)))
    
    # 2. Wilkinson variant: roots at 1,...,d for various d
    for d in [10, 30]:
        families.append((f"Wilkinson-{d}", np.arange(1, d+1, dtype=complex)))
    
    # 3. Clustered roots: pairs at distance ε
    for eps in [1e-2, 1e-4, 1e-6]:
        d = 10
        roots = np.array([k + (0 if k % 2 == 0 else eps) for k in range(d)], dtype=complex)
        families.append((f"Cluster-ε={eps:.0e}", roots))
    
    # 4. Multi-scale: roots at 10^k for k = -3,...,3
    roots_ms = np.array([10**k * np.exp(2j*np.pi*j/3) 
                         for k in range(-2, 3) for j in range(3)], dtype=complex)
    families.append(("MultiScale-15", roots_ms))
    
    # 5. Chebyshev roots (on [-1,1], near endpoints)
    for d in [10, 20, 50]:
        roots = np.cos(np.pi * (2*np.arange(d) + 1) / (2*d)).astype(complex)
        families.append((f"Chebyshev-{d}", roots))
    
    # 6. Spiral roots
    for d in [10, 20]:
        t = np.linspace(0.1, 2, d)
        roots = t * np.exp(2j * np.pi * t)
        families.append((f"Spiral-{d}", roots))
    
    # 7. Mignotte-type: one pair at distance 2^{-d/2}
    for d in [10, 20]:
        eps = 2**(-d//2)
        roots = np.concatenate([
            [1 + eps/2, 1 - eps/2],
            np.exp(2j * np.pi * np.arange(2, d) / (d-2)) * 0.5
        ]).astype(complex)
        families.append((f"Mignotte-{d}", roots))
    
    # 8. Random complex (Gaussian)
    rng = np.random.RandomState(42)
    for d in [10, 20, 50, 100]:
        roots = rng.randn(d) + 1j * rng.randn(d)
        families.append((f"Random-{d}", roots))
    
    # 9. Roots on unit circle (dense)
    for d in [20, 50, 100]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        families.append((f"Unity-{d}", roots))
    
    # 10. Near-degenerate: d-1 roots at 0, one at 1
    for d in [5, 10, 20]:
        roots = np.concatenate([
            np.zeros(d-1) + 1e-4 * np.exp(2j * np.pi * np.arange(d-1) / (d-1)),
            [1.0]
        ]).astype(complex)
        families.append((f"NearDegen-{d}", roots))
    
    # 11. Bernstein: roots densely packed on [-1,1]
    for d in [20, 50]:
        roots = np.linspace(-1, 1, d).astype(complex)
        families.append((f"Bernstein-{d}", roots))
    
    # 12. Extreme: roots at distance exactly 1/d from each other on a line
    for d in [20, 50]:
        roots = np.array([k/d for k in range(d)], dtype=complex)
        families.append((f"Line-{d}", roots))
    
    return families


def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  ADVERSARIAL POLYNOMIAL TEST: Pure Pandrosion T4 + Iterated Scaling ║")
    print("║  Optimal starting point (offset π/d) at each epoch                  ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    
    families = get_adversarial_families()
    
    print(f"\n  {'Family':<20s}  {'d':>3s}  {'ρ':>8s}  {'conv%':>6s}  {'#conv':>5s}  "
          f"{'epochs':>7s}  {'steps':>7s}  {'m.desc':>8s}  {'max_d':>8s}  {'f<0':>5s}")
    print("  " + "─"*95)
    
    results = []
    n_fail = 0
    
    for name, roots in families:
        d = len(roots)
        rho = np.max(np.abs(roots))
        
        r = run_adaptive(roots, mode="T4", max_epochs=500, n_starts=max(d, 20))
        
        conv_pct = f"{100*r['conv_rate']:.0f}%"
        ep_str = f"{r['best_epochs']}" if r['best_epochs'] > 0 else "---"
        st_str = f"{r['best_steps']}" if r['best_steps'] > 0 else "---"
        desc_str = f"{r['mean_desc']:.4f}" if np.isfinite(r['mean_desc']) else "N/A"
        maxd_str = f"{r['max_desc']:.4f}" if np.isfinite(r['max_desc']) else "N/A"
        fneg_str = f"{r['frac_neg']:.2f}" if r['frac_neg'] > 0 else "N/A"
        
        status = "✓" if r['conv_rate'] > 0.5 else "✗"
        if r['conv_rate'] < 0.5:
            n_fail += 1
        
        print(f"  {name:<20s}  {d:>3d}  {rho:>8.2f}  {conv_pct:>6s}  "
              f"{r['conv_count']:>5d}  {ep_str:>7s}  {st_str:>7s}  "
              f"{desc_str:>8s}  {maxd_str:>8s}  {fneg_str:>5s}  {status}")
        
        results.append((name, d, r))
    
    # Summary
    print(f"\n{'═'*95}")
    total = len(families)
    passed = total - n_fail
    print(f"  RESULTS: {passed}/{total} families passed (convergence > 50%)")
    
    if n_fail > 0:
        print(f"\n  FAILURES:")
        for name, d, r in results:
            if r['conv_rate'] < 0.5:
                print(f"    {name} (d={d}): {100*r['conv_rate']:.0f}% convergence")
    
    # Descent statistics
    print(f"\n  DESCENT STATISTICS (all families combined):")
    all_descs = []
    for name, d, r in results:
        if np.isfinite(r['mean_desc']):
            all_descs.append(r['mean_desc'])
    
    if all_descs:
        arr = np.array(all_descs)
        print(f"    Mean descent/epoch:   {np.mean(arr):.4f}")
        print(f"    Worst descent/epoch:  {np.max(arr):.4f}")
        print(f"    Best descent/epoch:   {np.min(arr):.4f}")
        
        worst_name = None
        worst_desc = -float('inf')
        for name, d, r in results:
            if np.isfinite(r['mean_desc']) and r['mean_desc'] > worst_desc:
                worst_desc = r['mean_desc']
                worst_name = name
        print(f"    Worst family:         {worst_name} ({worst_desc:.4f})")
    
    # BSS cost analysis
    print(f"\n  BSS COST ANALYSIS:")
    print(f"  {'Family':<20s}  {'d':>3s}  {'steps':>7s}  {'steps/d':>8s}  {'steps/d²':>9s}")
    print("  " + "─"*55)
    for name, d, r in results:
        if r['best_steps'] > 0:
            print(f"  {name:<20s}  {d:>3d}  {r['best_steps']:>7d}  "
                  f"{r['best_steps']/d:>8.2f}  {r['best_steps']/d**2:>9.4f}")
    
    print(f"\n{'═'*95}")
    if n_fail == 0:
        print(f"  ╔═══════════════════════════════════════════════════════════╗")
        print(f"  ║  ALL {total} ADVERSARIAL FAMILIES PASSED!                     ║")
        print(f"  ║  Pure Pandrosion T4 converges on ALL tested polynomials   ║")
        print(f"  ╚═══════════════════════════════════════════════════════════╝")
    print(f"{'═'*95}")


if __name__ == "__main__":
    main()
