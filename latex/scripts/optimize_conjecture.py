import numpy as np
from scipy.optimize import minimize, basinhopping
import time
import math

def eval_poly(roots, z):
    """P(z) = prod_k (z - root_k), monic."""
    return np.prod(z - roots)

def Q_diff(roots, z0, z):
    """Q(z0, z) = (P(z) - P(z0))/(z - z0), with P'(z0) at z = z0."""
    if abs(z - z0) < 1e-18:
        mat = z0 - roots
        prod_tot = np.prod(mat)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(np.where(mat != 0, prod_tot / mat, 0))
    return (eval_poly(roots, z) - eval_poly(roots, z0)) / (z - z0)

def sigma_log_and_Q(roots, z0, z):
    """Return (Sigma_log, Q_value)."""
    d = len(roots)
    dz = z - roots
    dz0 = z0 - roots
    Pz = np.prod(dz)
    Pz0 = np.prod(dz0)
    
    if np.any(np.abs(dz) < 1e-18) or np.any(np.abs(dz0) < 1e-18):
        Rk_z = np.array([np.prod(np.delete(dz, k)) for k in range(d)])
        Rk_z0 = np.array([np.prod(np.delete(dz0, k)) for k in range(d)])
    else:
        Rk_z = Pz / dz
        Rk_z0 = Pz0 / dz0
        
    if abs(z - z0) < 1e-18:
        Qks = np.zeros(d, dtype=complex)
        for k in range(d):
            for m in range(d):
                if m == k: continue
                t = 1.0 + 0j
                for j in range(d):
                    if j == k or j == m: continue
                    t *= (z0 - roots[j])
                Qks[k] += t
    else:
        Qks = (Rk_z - Rk_z0) / (z - z0)
        
    Q = (Pz - Pz0) / (z - z0) if abs(z - z0) > 1e-18 else Q_diff(roots, z0, z)
    if abs(Q) < 1e-80:
        return None, Q
    absQks = np.abs(Qks)
    absQ = abs(Q)
    mask = absQks > 1e-80
    if not mask.any():
        return None, Q
    logs = np.log(absQks[mask] / absQ)
    return float(np.sum(logs)), Q

def pandrosion_T3_step(roots, z0, z):
    """One T3 (Aitken on base map) from anchor z0, iterate z."""
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

def objective_function(x, d, K):
    """
    x contains: 
    - 2*d parameters for roots (real, imag)
    - 2 parameters for z0 (real, imag)
    Total parameters = 2*d + 2
    
    We want to MAXIMIZE Sigma_log, so we return -Sigma_log.
    """
    roots = x[:2*d].reshape((d, 2))
    roots_c = roots[:, 0] + 1j * roots[:, 1]
    
    z0_c = x[-2] + 1j * x[-1]
    z_c = z0_c
    
    # We want to keep roots somewhat bounded to avoid trivial infinity scaling
    # Let's add a penalty if max|root| > 1 or |z0| > 2
    penalty = 0.0
    max_root_abs = np.max(np.abs(roots_c))
    if max_root_abs > 1.0:
        penalty += 100 * (max_root_abs - 1.0)**2
        
    if abs(z0_c) > 2.0:
        penalty += 100 * (abs(z0_c) - 2.0)**2
        
    for step in range(K):
        z_new = pandrosion_T3_step(roots_c, z0_c, z_c)
        if z_new is None or not np.isfinite(z_new):
            return 1e6 # Invalid step, large penalty
        z_c = z_new
        
    sigma, _ = sigma_log_and_Q(roots_c, z0_c, z_c)
    if sigma is None:
        return 1e6
        
    # Return negative since we minimize
    return -sigma + penalty

def search_counter_example(d=5, K=1, n_trials=5):
    print(f"=== Searching for counter-examples: d={d}, K={K} ===")
    
    best_sigma = -float('inf')
    best_x = None
    
    for trial in range(n_trials):
        # Random initialization
        roots_init = np.random.randn(d, 2) * 0.5
        z0_init = np.random.randn(2)
        x0 = np.concatenate([roots_init.flatten(), z0_init])
        
        # Optimize using L-BFGS-B
        res = minimize(objective_function, x0, args=(d, K), method='L-BFGS-B', 
                       options={'maxfun': 50000, 'maxiter': 5000})
        
        # Evaluate without penalty
        roots_c = res.x[:2*d].reshape((d, 2))
        roots_c = roots_c[:, 0] + 1j * roots_c[:, 1]
        z0_c = res.x[-2] + 1j * res.x[-1]
        
        z_c = z0_c
        for _ in range(K):
            z_new = pandrosion_T3_step(roots_c, z0_c, z_c)
            if z_new is not None:
                z_c = z_new
                
        sigma, _ = sigma_log_and_Q(roots_c, z0_c, z_c)
        
        print(f"Trial {trial+1}/{n_trials} - Max Sigma found: {sigma:.6f}")
        
        if sigma > best_sigma:
            best_sigma = sigma
            best_x = res.x
            
        if sigma > 0:
            print(">>> COUNTER-EXAMPLE FOUND! Sigma > 0 <<<")
            print("Roots:", roots_c)
            print("z0:", z0_c)
            return roots_c, z0_c, sigma
            
    print(f"Best Sigma found overall: {best_sigma:.6f}")
    return None, None, best_sigma

if __name__ == "__main__":
    # Test for different degrees and steps
    for d in [3, 5, 8]:
        for K in [1, 2, 3]:
            search_counter_example(d=d, K=K, n_trials=3)
            print("-" * 50)
