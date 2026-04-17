import numpy as np

def P_eval(roots, z):
    r = 1.0+0j
    for rk in roots: r *= (z-rk)
    return r

def pandrosion_F(roots, z0, z):
    Pz = P_eval(roots, z); Pz0 = P_eval(roots, z0)
    denom = Pz - Pz0
    if abs(denom) < 1e-30: return None
    return (z0 * Pz - z * Pz0) / denom

def test_dichotomy(trials=1000):
    np.random.seed(42)
    min_c = float('inf')
    violations = 0
    total_steps = 0
    
    for _ in range(trials):
        d = np.random.randint(10, 50)
        roots = np.random.randn(d) + 1j*np.random.randn(d)
        
        rho = max(abs(r_k) for r_k in roots)
        R = 1 + rho
        
        theta_0 = np.random.uniform(0, 2*np.pi)
        z0 = R * np.exp(1j*theta_0)
        
        theta = np.random.uniform(0, 2*np.pi)
        z = R * np.exp(1j*theta)
        
        for step in range(100):
            Pz = P_eval(roots, z)
            F = pandrosion_F(roots, z0, z)
            if F is None: break
            
            # Check if we locked in
            dist = min(abs(F - rk) for rk in roots)
            if dist < 1e-3:
                break
                
            PF = P_eval(roots, F)
            if abs(PF) == 0: break
            
            alg_progress = np.log(abs(Pz) / abs(PF))
            geom_progress = d * abs(F - z) / abs(F - z0)
            
            current_c = max(alg_progress, geom_progress)
            
            if current_c <= 0 and violations < 5:
                print(f"Violation details! alg={alg_progress:.6f}, geom={geom_progress:.6f}")
                print(f"|Pz| = {abs(Pz):.2e}, |PF| = {abs(PF):.2e}")
                print(f"|z-z0| = {abs(z-z0):.2e}, |F-z0| = {abs(F-z0):.2e}, |F-z| = {abs(F-z):.2e}")
                print(f"r_n = {abs(Pz)/abs(P_eval(roots, z0)):.6f}")
            
            if current_c < min_c:
                min_c = current_c
            
            if current_c <= 0:
                violations += 1
            
            total_steps += 1
            z = F
            
    print(f"Violations: {violations} out of {total_steps} steps")

test_dichotomy()
