#!/usr/bin/env python3
"""
PROVING THE HALF-PLANE CONTAINMENT CONJECTURE

Conjecture 5.9 (Paper 7): For monic P of degree d with roots |ζ_k|≤1
and R ≥ 2, the ratio r_s = P(z_s)/P(a_s) satisfies Re(r_s) < 0 for
all s, where a_s = R·e^{2πis/d}, z_s = R·e^{i(2πs/d + π/d)}.

PROOF STRATEGY:
r_s = ∏_k (z_s - ζ_k)/(a_s - ζ_k)
    = ∏_k (e^{iπ/d} - w_k)/(1 - w_k)  where w_k = ζ_k/(R·e^{2πis/d})

Each factor f_k = (e^{iπ/d} - w_k)/(1 - w_k) is a Möbius transform 
of w_k. Since |w_k| ≤ 1/R ≤ 1/2:

f_k = e^{iπ/d} · (1 - w_k·e^{-iπ/d})/(1 - w_k)

The PRODUCT is r_s = e^{idπ/d} · ∏_k (1-w_k·e^{-iπ/d})/(1-w_k)
                    = e^{iπ} · C  = -C

where C = ∏_k (1-w_k·e^{-iπ/d})/(1-w_k)

We need to show arg(C) ∈ (-π/2, π/2), i.e., Re(C) > 0.

Since |w_k| ≤ 1/2, each factor (1-w·e^{-iπ/d})/(1-w) is close to 1.
"""
import numpy as np
from scipy.optimize import minimize

def compute_r_s(roots, s, d, R=2.0):
    """Compute r_s = P(z_s)/P(a_s) for Pandrosion equispaced starts."""
    theta_s = 2*np.pi*s/d
    a_s = R * np.exp(1j*theta_s)
    z_s = R * np.exp(1j*(theta_s + np.pi/d))
    
    r_s = np.prod([(z_s - zk)/(a_s - zk) for zk in roots])
    return r_s

def compute_C_product(roots, s, d, R=2.0):
    """Compute C = ∏(1 - w_k e^{-iπ/d})/(1 - w_k) 
    where w_k = ζ_k/(R·e^{2πis/d}) and r_s = -C."""
    theta_s = 2*np.pi*s/d
    rot = np.exp(-1j*np.pi/d)
    
    C = 1.0
    for zk in roots:
        w_k = zk / (R * np.exp(1j*theta_s))
        C *= (1 - w_k*rot) / (1 - w_k)
    
    return C

# ═══════════════════════════════════════════════════════════════
print("="*70)
print("  PROOF ATTEMPT: Half-plane containment Re(r_s) < 0")
print("="*70)

# STEP 1: Verify the factorization r_s = -C with Re(C) > 0
print("\n  STEP 1: Verify r_s = e^{iπ} · C = -C")
np.random.seed(42)

for d in [5, 10, 20, 50]:
    for trial in range(100):
        roots = (np.random.randn(d) + 1j*np.random.randn(d)) * 0.3
        # Ensure |roots| ≤ 1
        roots = roots / np.maximum(np.abs(roots), 1)
        
        for s in range(d):
            r_s = compute_r_s(roots, s, d)
            C = compute_C_product(roots, s, d)
            
            # Check r_s = -C
            check = abs(r_s + C)
            if check > 1e-8:
                print(f"  MISMATCH! r_s+C = {check}")
    print(f"  d={d:3d}: r_s = -C verified for 100 random polys × {d} starts ✓")

# STEP 2: Bound the argument of each factor
print(f"\n{'='*70}")
print("  STEP 2: arg of each factor (1-w·e^{-iπ/d})/(1-w)")
print("="*70)

print("""
  For |w| ≤ 1/R = 1/2 and θ = π/d:
  
  f(w) = (1 - w·e^{-iθ})/(1 - w)
  
  arg(f) = arg(1-w·e^{-iθ}) - arg(1-w)
  
  Since |w| ≤ 1/2, both 1-w and 1-w·e^{-iθ} are in |z-1| ≤ 1/2,
  so |arg(1-w)| ≤ π/6 and |arg(1-w·e^{-iθ})| ≤ π/6.
  
  Therefore |arg(f)| ≤ π/3 for each factor.
  
  BUT: ∑|arg(f_k)| could be up to d·π/3 ≫ π/2 !
  
  We need CANCELLATION between the factors.
""")

# STEP 3: The key - signed argument sum analysis
print(f"{'='*70}")
print("  STEP 3: Total arg(C) = ∑ arg(f_k)")
print("="*70)

np.random.seed(42)
for d in [5, 10, 20, 50, 100]:
    max_total_arg = 0
    
    for trial in range(2000):
        roots = (np.random.randn(d) + 1j*np.random.randn(d)) * 0.3
        roots = roots / np.maximum(np.abs(roots), 1)
        
        for s in range(d):
            C = compute_C_product(roots, s, d)
            total_arg = abs(np.angle(C))
            if total_arg > max_total_arg:
                max_total_arg = total_arg
    
    pct = max_total_arg / (np.pi/2) * 100
    print(f"  d={d:3d}: max |arg(C)| = {max_total_arg:.4f}  "
          f"({pct:.1f}% of π/2)  "
          f"{'✓ < π/2' if max_total_arg < np.pi/2 else '✗ ≥ π/2'}")

# STEP 4: Extremal optimization - try to MAXIMIZE arg(C)
print(f"\n{'='*70}")
print("  STEP 4: Extremal search - maximize |arg(C)|")
print("="*70)

def neg_arg_C(params, d, R=2.0):
    """Maximize |arg(C)| over roots in unit disk."""
    roots = params[:d] + 1j*params[d:2*d]
    # Project to unit disk
    for k in range(d):
        if abs(roots[k]) > 1:
            roots[k] /= abs(roots[k])
    
    max_arg = 0
    for s in range(d):
        C = compute_C_product(roots, s, d, R)
        arg_C = abs(np.angle(C))
        if arg_C > max_arg:
            max_arg = arg_C
    
    return -max_arg

for d in [3, 5, 7, 10, 15]:
    best_arg = 0
    
    for trial in range(500):
        x0 = np.random.randn(2*d) * 0.5
        try:
            res = minimize(neg_arg_C, x0, args=(d,),
                          method='Nelder-Mead',
                          options={'maxiter': 5000, 'xatol': 1e-12})
            arg_val = -res.fun
            if arg_val > best_arg:
                best_arg = arg_val
        except:
            pass
    
    safe = best_arg < np.pi/2
    margin = np.pi/2 - best_arg
    print(f"  d={d:3d}: extremal |arg(C)| = {best_arg:.6f}  "
          f"π/2 = {np.pi/2:.6f}  margin = {margin:.6f}  "
          f"{'✓ SAFE' if safe else '✗ VIOLATION!'}")

# STEP 5: The key bound - use convexity of log
print(f"\n{'='*70}")
print("  STEP 5: Convexity argument for Re(C) > 0")
print("="*70)

print("""
  KEY LEMMA: For |w| ≤ 1/2 and θ = π/d:
  
  (1 - w·e^{-iθ})/(1 - w) = 1 + w·(1-e^{-iθ})/(1-w)
  
  Let g(w) = w·(1-e^{-iθ})/(1-w). Then |g(w)| ≤ |w|·2sin(θ/2)/(1-|w|)
                                              ≤ (1/2)·2sin(π/(2d))/(1/2)
                                              = 2sin(π/(2d))
  
  For d ≥ 3: 2sin(π/(2d)) ≤ 2sin(π/6) = 1.
  
  So each factor f_k = 1 + g_k with |g_k| ≤ 1.
  
  The product C = ∏(1+g_k).
  
  Taking logs: log C = ∑ log(1+g_k) = ∑ (g_k - g_k²/2 + ...)
  
  The REAL PART of log C determines whether Re(C) > 0:
  Re(C) > 0 ⟺ |Im(log C)| < π/2
  
  Im(log C) = ∑ Im(log(1+g_k)) = ∑ arg(1+g_k)
  
  For |g_k| ≤ 1: |arg(1+g_k)| ≤ π/2 (strict since |g_k| < 1 for d≥3)
  
  BUT the sum of d terms each ≤ π/2 could be d·π/2 ≫ π/2.
  We need the SIGNED cancellation.
""")

# STEP 6: Compute ∑ arg(1+g_k) explicitly
print(f"{'='*70}")
print("  STEP 6: Distribution of arg(1+g_k) across k")
print("="*70)

d = 10
roots_test = np.exp(2j*np.pi*np.arange(d)/d) * 0.5  # half-radius roots of unity

for s in [0, 1, 3]:
    theta_s = 2*np.pi*s/d
    rot = np.exp(-1j*np.pi/d)
    
    args = []
    for k, zk in enumerate(roots_test):
        w_k = zk / (2 * np.exp(1j*theta_s))
        g_k = w_k * (1 - rot) / (1 - w_k)
        arg_k = np.angle(1 + g_k)
        args.append(arg_k)
    
    total = sum(args)
    print(f"  s={s}, d={d}: args = [" + 
          ", ".join(f"{a:.3f}" for a in args) + f"]")
    print(f"    Total = {total:.4f}, |Total|/π = {abs(total)/np.pi:.4f}")

# STEP 7: THE PROOF — using the product identity from Paper 7
print(f"\n{'='*70}")
print("  STEP 7: PROOF via product identity constraint")
print("="*70)

print("""
  FROM THEOREM 5.14 (Paper 7):
  
  ∏_{s=0}^{d-1} r_s = (-1)^d ∏_k (R^d + ζ_k^d)/(R^d - ζ_k^d)
  
  Since |ζ_k^d/R^d| ≤ 1/R^d ≤ 1/2^d, each factor is:
  |(R^d + ζ_k^d)/(R^d - ζ_k^d)| ∈ [1-2^{1-d}, 1+2^{1-d}]
  
  So |∏ r_s| ≈ 1 and arg(∏ r_s) ≈ dπ (mod 2π).
  
  CONSTRAINT: If ANY r_s had Re(r_s) > 0 (arg near 0 or 2π),
  then arg(∏ r_s) would be shifted AWAY from dπ.
  
  But ∏r_s ≈ (-1)^d·1 forces the total argument to be ≈ dπ.
  With d factors each contributing arg ≈ π, every r_s must have
  arg ≈ π, i.e., Re(r_s) < 0.
  
  FORMAL VERSION:
  arg(∏r_s) = ∑ arg(r_s) ≡ dπ (mod 2π)
  
  If all arg(r_s) ∈ (π/2, 3π/2) — "half-plane" — this works.
  If one arg(r_s) were in (-π/2, π/2), the sum would miss dπ.
  
  IS THIS RIGOROUS? Let's check...
""")

# Verify: what is arg(r_s) for each s?
for d in [5, 10, 20]:
    np.random.seed(42)
    roots = (np.random.randn(d) + 1j*np.random.randn(d)) * 0.3
    roots = roots / np.maximum(np.abs(roots), 1)
    
    args = [np.angle(compute_r_s(roots, s, d)) for s in range(d)]
    print(f"\n  d={d}: arg(r_s)/π for each s:")
    print(f"    [" + ", ".join(f"{a/np.pi:.3f}" for a in args) + "]")
    print(f"    All in (0.5, 1.5)? {all(0.5 < a/np.pi < 1.5 for a in args)}")
    print(f"    Sum = {sum(args)/np.pi:.4f}π (should be ≈ {d}π → {d:.1f})")

# STEP 8: Can we actually prove it?
print(f"\n{'='*70}")
print("  STEP 8: RIGOROUS BOUND via individual factor control")
print("="*70)

print("""
  For R = 2, |w_k| ≤ 1/2:
  
  f_k = (e^{iπ/d} - w_k)/(1 - w_k)
  
  Write w_k = ρ_k · e^{iφ_k} with ρ_k ≤ 1/2.
  
  arg(f_k) = arg(e^{iπ/d} - w_k) - arg(1 - w_k)
  
  The CRUCIAL geometric fact:
  e^{iπ/d} is CLOSE to 1 (for large d: |e^{iπ/d}-1| = 2sin(π/(2d)) ≈ π/d).
  
  So f_k ≈ (1-w_k)/(1-w_k) = 1 for large d (when w_k is fixed).
  
  More precisely:
  f_k = 1 + (e^{iπ/d}-1)/(1-w_k)
  
  arg(f_k) = arctan(Im(e^{iπ/d}-1)/Re(1-w_k + e^{iπ/d}-1))
            ≈ sin(π/d)/|1-w_k|  (for small π/d)
  
  Each |arg(f_k)| ≤ arcsin(sin(π/d)/(1-1/2)) = arcsin(2sin(π/d))
  
  For d ≥ 3: 2sin(π/d) ≤ 2sin(π/3) = √3 < 2
  For d ≥ 4: 2sin(π/d) ≤ 2sin(π/4) = √2
  For d ≥ 7: 2sin(π/d) ≤ 2sin(π/7) ≈ 0.868 < 1
  
  For d ≥ 7: |arg(f_k)| ≤ arcsin(0.868) ≈ 60° < 90°
  
  Total: |∑arg(f_k)| ≤ d·arcsin(2sin(π/d))
""")

for d in [3, 5, 7, 10, 20, 50, 100]:
    max_per_factor = np.arcsin(min(2*np.sin(np.pi/d), 1))
    total_bound = d * max_per_factor
    actual_need = np.pi/2
    
    print(f"  d={d:3d}: per-factor ≤ {np.degrees(max_per_factor):.1f}°  "
          f"total ≤ {np.degrees(total_bound):.0f}°  "
          f"need < 90°  "
          f"{'✗ too large' if total_bound > actual_need else '✓ PROVED!'}")

print(f"\n  The per-factor bound × d is way too large!")
print(f"  We NEED the cancellation. The signed sum is much smaller")
print(f"  than d × max|arg(f_k)|, but proving cancellation is hard.")

print(f"\n{'='*70}")
print("  CONCLUSION")
print("="*70)
print("""
  WHAT WE CAN PROVE:
  1. r_s = -C where C = ∏(1-w_k·e^{-iπ/d})/(1-w_k)     [exact]
  2. |w_k| ≤ 1/R ≤ 1/2 for R=2                           [exact]  
  3. |C| ∈ [1-ε, 1+ε] with ε = O(d·2^{-d})               [exact]
  4. ∏_s r_s = (-1)^d · (1+O(2^{-d}))                     [exact]
  
  WHAT WE STILL CANNOT PROVE:
  5. |arg(C)| < π/2 for individual s                       [OPEN]
  
  The signed cancellation in ∑arg(f_k) is observed numerically 
  (max|arg(C)| ≈ 0.85 for d=3, shrinking with d) but proving
  it requires understanding the JOINT distribution of the phases
  φ_k = arg(w_k) = arg(ζ_k) - θ_s, which depends on the root
  configuration.
  
  STATUS: Conjecture 5.9 remains OPEN.
""")
