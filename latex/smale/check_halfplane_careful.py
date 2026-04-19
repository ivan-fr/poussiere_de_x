#!/usr/bin/env python3
"""
CAREFUL CHECK: |arg(C)| with roots STRICTLY in unit disk.
Is the threshold π/2 or π ?
"""
import numpy as np
from scipy.optimize import minimize, differential_evolution

def compute_C_and_r(roots, s, d, R=2.0):
    """Compute C and r_s with roots strictly in unit disk."""
    theta_s = 2*np.pi*s/d
    a_s = R * np.exp(1j*theta_s)
    z_s = R * np.exp(1j*(theta_s + np.pi/d))
    
    r_s = np.prod([(z_s - zk)/(a_s - zk) for zk in roots])
    
    rot = np.exp(-1j*np.pi/d)
    C = 1.0+0j
    for zk in roots:
        w_k = zk / (R * np.exp(1j*theta_s))
        C *= (1 - w_k*rot) / (1 - w_k)
    
    return C, r_s

def project_to_disk(roots):
    """Project roots to closed unit disk."""
    for k in range(len(roots)):
        if abs(roots[k]) > 1.0:
            roots[k] = roots[k] / abs(roots[k]) * 0.999
    return roots

# ═══════════════════════════════════════════════════════════════
print("="*70)
print("  CAREFUL CHECK: |arg(C)| with |ζ_k| ≤ 1 strictly enforced")
print("="*70)

# TEST 1: Random roots in unit disk — what is max|arg(C)|?
np.random.seed(42)

for d in [3, 5, 7, 10, 15, 20, 50]:
    max_arg_C = 0
    max_re_r = -1e10  # max Re(r_s) — should be < 0
    
    for trial in range(5000):
        # Random roots STRICTLY in unit disk
        r = np.sqrt(np.random.random(d))  # uniform in disk
        theta = 2*np.pi*np.random.random(d)
        roots = r * np.exp(1j*theta)
        
        for s in range(d):
            C, r_s = compute_C_and_r(roots, s, d)
            arg_C = abs(np.angle(C))
            if arg_C > max_arg_C:
                max_arg_C = arg_C
            if r_s.real > max_re_r:
                max_re_r = r_s.real
    
    print(f"  d={d:3d}: max|arg(C)| = {max_arg_C:.4f} = {max_arg_C/np.pi:.3f}π  "
          f"max Re(r_s) = {max_re_r:.6f}  "
          f"{'Re<0 ✓' if max_re_r < 0 else 'Re≥0 ✗!'}")

# TEST 2: Extremal optimization WITH ENFORCEMENT  
print(f"\n{'='*70}")
print("  EXTREMAL: Maximize |arg(C)| with |ζ_k| ≤ 1 ENFORCED")
print("="*70)

def neg_arg_C_constrained(params, d, R=2.0):
    """Maximize |arg(C)| over roots in unit disk (enforced)."""
    n = d
    roots = np.zeros(n, dtype=complex)
    for k in range(n):
        r_k = min(abs(params[2*k]), 1.0)  # enforce |root| ≤ 1
        phi_k = params[2*k+1]
        roots[k] = r_k * np.exp(1j*phi_k)
    
    max_arg = 0
    for s in range(d):
        C, _ = compute_C_and_r(roots, s, d, R)
        arg_C = abs(np.angle(C))
        if arg_C > max_arg:
            max_arg = arg_C
    
    return -max_arg

def neg_re_r_constrained(params, d, R=2.0):
    """Maximize Re(r_s) — DIRECT test of conjecture violation."""
    n = d
    roots = np.zeros(n, dtype=complex)
    for k in range(n):
        r_k = min(abs(params[2*k]), 1.0)
        phi_k = params[2*k+1]
        roots[k] = r_k * np.exp(1j*phi_k)
    
    max_re = -1e10
    for s in range(d):
        _, r_s = compute_C_and_r(roots, s, d, R)
        if r_s.real > max_re:
            max_re = r_s.real
    
    return -max_re  # minimize → maximize Re(r_s)

for d in [3, 5, 7, 10]:
    best_arg = 0
    best_re = -1e10
    
    for trial in range(1000):
        x0 = np.random.randn(2*d)
        x0[::2] = np.abs(x0[::2]) * 0.5  # start with |root| ~ 0.5
        
        # Maximize |arg(C)|
        try:
            res = minimize(neg_arg_C_constrained, x0, args=(d,),
                          method='Nelder-Mead',
                          options={'maxiter': 10000, 'xatol': 1e-14})
            arg_val = -res.fun
            if arg_val > best_arg:
                best_arg = arg_val
        except:
            pass
        
        # Maximize Re(r_s) — DIRECT violation test
        try:
            res = minimize(neg_re_r_constrained, x0, args=(d,),
                          method='Nelder-Mead',
                          options={'maxiter': 10000, 'xatol': 1e-14})
            re_val = -res.fun
            if re_val > best_re:
                best_re = re_val
        except:
            pass
    
    margin_half = np.pi/2 - best_arg
    margin_pi = np.pi - best_arg
    
    print(f"  d={d:3d}: max|arg(C)| = {best_arg:.6f} = {best_arg/np.pi:.4f}π  "
          f"< π/2? {'✓' if best_arg < np.pi/2 else '✗'}  "
          f"< π? {'✓' if best_arg < np.pi else '✗'}  "
          f"max Re(r_s) = {best_re:.8f}  "
          f"{'Re<0 ✓' if best_re < 0 else 'Re≥0 ✗!'}")

# TEST 3: Focus on d=3 with polar coordinates (most likely to violate)
print(f"\n{'='*70}")
print("  FOCUSED: d=3, roots on unit circle (boundary)")
print("="*70)

d = 3
best_arg = 0
best_re = -1e10
best_roots = None

for trial in range(50000):
    # Roots on unit circle (extreme case)
    phi = 2*np.pi*np.random.random(d)
    roots = np.exp(1j*phi)
    
    for s in range(d):
        C, r_s = compute_C_and_r(roots, s, d)
        arg_C = abs(np.angle(C))
        
        if arg_C > best_arg:
            best_arg = arg_C
            best_roots = (roots.copy(), s, C, r_s)
        
        if r_s.real > best_re:
            best_re = r_s.real

print(f"  d={d}: max|arg(C)| = {best_arg:.6f} = {best_arg/np.pi:.4f}π")
print(f"  d={d}: max Re(r_s)  = {best_re:.8f}")

if best_roots:
    rr, ss, CC, rr_s = best_roots
    print(f"  Extremal config: s={ss}")
    print(f"    roots = {[f'{z:.4f}' for z in rr]}")
    print(f"    C = {CC:.6f}  |C| = {abs(CC):.6f}  arg(C) = {np.angle(CC):.6f}")
    print(f"    r_s = {rr_s:.6f}  Re(r_s) = {rr_s.real:.8f}")

# TEST 4: THE KEY QUESTION
print(f"\n{'='*70}")
print("  THE KEY QUESTION: What is the supremum of |arg(C)|?")
print("="*70)

# Use differential_evolution (global optimizer) for d=3
for d in [3, 5]:
    bounds = [(0, 1), (-np.pi, np.pi)] * d
    
    result = differential_evolution(neg_arg_C_constrained, bounds, args=(d,),
                                    maxiter=2000, seed=42, tol=1e-15,
                                    polish=True, popsize=50)
    
    best_arg = -result.fun
    
    # Also try to directly maximize Re(r_s)
    result2 = differential_evolution(neg_re_r_constrained, bounds, args=(d,),
                                     maxiter=2000, seed=42, tol=1e-15,
                                     polish=True, popsize=50)
    
    best_re = -result2.fun
    
    print(f"  d={d}: GLOBAL OPT max|arg(C)| = {best_arg:.8f} = {best_arg/np.pi:.6f}π")
    print(f"         GLOBAL OPT max Re(r_s) = {best_re:.10f}  "
          f"{'✓ Re<0: CONJECTURE HOLDS!' if best_re < 0 else '✗ Re≥0: COUNTEREXAMPLE!'}")
    print(f"         < π/2={np.pi/2:.6f}? {'YES ✓' if best_arg < np.pi/2 else 'NO'}")
    print(f"         < π={np.pi:.6f}?   {'YES ✓' if best_arg < np.pi else 'NO'}")
    print()

print("="*70)
print("  SUMMARY")
print("="*70)
print("""
  With roots STRICTLY in the unit disk |ζ_k| ≤ 1:
  
  - Previous "violations" were due to roots ESCAPING the unit disk
    during optimization (no proper constraint enforcement)
  
  - With proper enforcement, the maximum |arg(C)| and Re(r_s) 
    values determine whether the conjecture holds.
""")
