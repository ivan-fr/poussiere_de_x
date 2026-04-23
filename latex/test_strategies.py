import numpy as np
import time
import math

RNG = np.random.default_rng(20260424)

# ---------------------------------------------------------------------------
# Polynomial and Base Math
# ---------------------------------------------------------------------------

def eval_poly(roots, z):
    return np.prod(z - roots)

def Q_diff(roots, z0, z):
    if abs(z - z0) < 1e-18:
        mat = z0 - roots
        prod_tot = np.prod(mat)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(np.where(mat != 0, prod_tot / mat, 0))
    return (eval_poly(roots, z) - eval_poly(roots, z0)) / (z - z0)

def pandrosion_T3_step(roots, z0, z):
    P0 = eval_poly(roots, z0)
    Q = Q_diff(roots, z0, z)
    if abs(Q) < 1e-60: return None
    s1 = z0 - P0 / Q
    Q2 = Q_diff(roots, z0, s1)
    if abs(Q2) < 1e-60: return s1
    s2 = z0 - P0 / Q2
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-60: return s1
    return z - (s1 - z) ** 2 / den

# ---------------------------------------------------------------------------
# Adversarial families
# ---------------------------------------------------------------------------

def family_clustered_line(d):
    return np.linspace(-1, 1, d).astype(complex)

def family_roots_unity(d):
    return np.exp(2j * np.pi * np.arange(d) / d)

def family_counter_example(d=3):
    # The specific counter-example trap we found earlier
    return np.array([-0.96731984-0.06459744j,  0.22987706-0.45925592j,  0.06489314-0.22842954j])

FAMILIES = {
    "roots_unity": family_roots_unity,
    "clustered_line": family_clustered_line,
    "counter_trap": family_counter_example
}

# ---------------------------------------------------------------------------
# Strategy 1: Fallback Non-Holomorphe
# ---------------------------------------------------------------------------

def run_strategy_1(roots, z_init, max_epochs, K=3):
    z = z_init
    ops = 0
    d = len(roots)
    fallback_count = 0
    for ep in range(max_epochs):
        z0 = z
        P0 = eval_poly(roots, z0)
        ops += d
        if abs(P0) < 1e-35:
            return True, ops, fallback_count
        
        # K steps
        for _ in range(K):
            z_new = pandrosion_T3_step(roots, z0, z)
            ops += 3 * d # roughly 3 evaluations per T3 step
            if z_new is None or not np.isfinite(z_new):
                return True, ops, fallback_count
            z = z_new
            
        PK = eval_poly(roots, z)
        ops += d
        
        # Descent check (Fallback condition)
        if abs(PK) > 0.99 * abs(P0):
            # Fallback: Fractional Newton step
            fallback_count += 1
            Pprime0 = Q_diff(roots, z0, z0)
            ops += d
            if abs(Pprime0) > 1e-60:
                z = z0 - 0.25 * P0 / Pprime0
            else:
                return False, ops, fallback_count # completely stuck
                
    return False, ops, fallback_count

# ---------------------------------------------------------------------------
# Strategy 2: Multi-Start with Pruning
# ---------------------------------------------------------------------------

def run_strategy_2(roots, max_epochs, K=3):
    d = len(roots)
    rho = float(np.max(np.abs(roots)))
    R = 1 + rho # Cauchy circle
    
    # Initialize d starting points
    active_paths = []
    for s in range(d):
        theta = 2 * np.pi * s / d + np.pi / (2 * max(d, 3))
        active_paths.append(R * np.exp(1j * theta))
        
    ops = 0
    epochs_taken = 0
    for ep in range(max_epochs):
        epochs_taken += 1
        next_paths = []
        for z_init in active_paths:
            z0 = z_init
            z = z0
            P0 = eval_poly(roots, z0)
            ops += d
            if abs(P0) < 1e-35:
                return True, ops, epochs_taken, len(active_paths)
                
            # K steps
            ok = True
            for _ in range(K):
                z_new = pandrosion_T3_step(roots, z0, z)
                ops += 3 * d
                if z_new is None or not np.isfinite(z_new):
                    ok = False
                    break
                z = z_new
                
            if not ok:
                return True, ops, epochs_taken, len(active_paths)
                
            PK = eval_poly(roots, z)
            ops += d
            
            # Pruning condition
            if abs(PK) <= 0.99 * abs(P0):
                next_paths.append(z)
                
        active_paths = next_paths
        if not active_paths:
            return False, ops, epochs_taken, 0 # all paths died!
            
    return False, ops, epochs_taken, len(active_paths)

# ---------------------------------------------------------------------------
# Test Runner
# ---------------------------------------------------------------------------

def main():
    degrees = [3, 10, 50, 100]
    
    print("=" * 80)
    print("Strategy 1: Single Path with Fractional Newton Fallback (z0 = R)")
    print("=" * 80)
    print(f"{'Family':<15} {'d':>4} {'Success':>8} {'Total Ops':>12} {'Fallbacks':>10}")
    for family in FAMILIES:
        for d in degrees:
            if family == "counter_trap" and d != 3: continue
            roots = FAMILIES[family](d)
            if family == "counter_trap":
                # Start exactly AT the trap for the counter_trap to stress test the fallback
                z_init = -0.6496697856287815-0.28032195127130505j
            else:
                rho = float(np.max(np.abs(roots)))
                z_init = (1 + rho) * np.exp(1j * np.pi / 7)
                
            success, ops, fallbacks = run_strategy_1(roots, z_init, max_epochs=200)
            print(f"{family:<15} {d:>4} {str(success):>8} {ops:>12} {fallbacks:>10}")
            
    print("\n" + "=" * 80)
    print("Strategy 2: Multi-Start with Pruning (d paths on Cauchy Circle)")
    print("=" * 80)
    print(f"{'Family':<15} {'d':>4} {'Success':>8} {'Total Ops':>12} {'Surviving':>10}")
    for family in FAMILIES:
        for d in degrees:
            if family == "counter_trap" and d != 3: continue
            roots = FAMILIES[family](d)
            success, ops, eps, survived = run_strategy_2(roots, max_epochs=200)
            print(f"{family:<15} {d:>4} {str(success):>8} {ops:>12} {survived:>10}")

if __name__ == "__main__":
    main()
