"""
Find the CORRECT theorem: existential half-plane containment.

Key discovery: Re(r_s) < 0 fails for ALL s (universal), 
but ALWAYS holds for SOME s (existential).

Moreover, the product identity ∏_s r_s = (-1)^d + O(2^{-d})
constrains the problem algebraically: the product of all r_s
has modulus ≈ 1 and argument ≈ dπ (mod 2π).

THEOREM TO PROVE: There exists s ∈ {0,...,d-1} such that Re(r_s) < 0.
Equivalently: ∑_s 1_{Re(r_s) > 0} < d.

Actually, the product identity gives something STRONGER.
"""

import numpy as np
from itertools import product as iprod

def compute_all_rs(roots, R):
    d = len(roots)
    thetas = 2 * np.pi * np.arange(d) / d
    a_s = R * np.exp(1j * thetas)
    z_s = R * np.exp(1j * (thetas + np.pi/d))
    
    rs = np.array([np.prod(z_s[s] - roots) / np.prod(a_s[s] - roots) for s in range(d)])
    return rs

print("=" * 70)
print("THE PRODUCT IDENTITY: ∏_s r_s = (-1)^d · ∏_k (R^d+ζ^d)/(R^d-ζ^d)")
print("=" * 70)

for d in [5, 10, 20]:
    roots = 0.99 * np.ones(d) + 0.01j * np.linspace(-1, 1, d)
    R = 2.0
    rs = compute_all_rs(roots, R)
    
    prod_r = np.prod(rs)
    expected = (-1)**d
    for z in roots:
        expected *= (R**d + z**d) / (R**d - z**d)
    
    print(f"\nd={d}: ∏r_s = {prod_r:.6f}")
    print(f"  Expected = {expected:.6f}")
    print(f"  |∏r_s| = {abs(prod_r):.6f}")
    print(f"  arg(∏r_s) = {np.degrees(np.angle(prod_r)):.2f}°")
    print(f"  Sum arg(r_s) = {np.degrees(sum(np.angle(rs))):.2f}°")
    
    n_pos = sum(np.real(rs) > 0)
    print(f"  Starts with Re > 0: {n_pos}/{d}")
    print(f"  Starts with Re < 0: {d - n_pos}/{d}")

print("\n" + "=" * 70)
print("EXISTENTIAL THEOREM: ∃ s with Re(r_s) < 0")
print("Proof strategy via the product identity")
print("=" * 70)

print("""
The product identity gives: ∏_{s=0}^{d-1} r_s = (-1)^d · C_d
where C_d = ∏_k (R^d + ζ_k^d)/(R^d - ζ_k^d) and |C_d| → 1.

Key argument:
- arg(∏ r_s) = d·π + arg(C_d) ≡ d·π + small  (mod 2π)
- The sum of arguments: ∑ arg(r_s) = d·π + arg(C_d) + 2kπ

If ALL r_s had Re(r_s) > 0, then arg(r_s) ∈ (-π/2, π/2) for all s.
So ∑ arg(r_s) ∈ (-dπ/2, dπ/2).

But we need ∑ arg(r_s) = dπ + arg(C_d) + 2kπ.
For d odd: dπ is NOT in (-dπ/2, dπ/2) since dπ > dπ/2. CONTRADICTION!
For d even: dπ ≡ 0 (mod 2π) could be in the range. Need more care.

Actually, let me be more precise...
""")

# Verify the argument constraint
print("Verifying the argument constraint:")
for d in [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50, 100]:
    # Test many random root configurations
    violations = 0
    total = 1000
    
    for trial in range(total):
        roots = np.random.rand(d) * np.exp(2j * np.pi * np.random.rand(d))
        roots /= np.maximum(abs(roots), 1)
        R = 2.0
        
        rs = compute_all_rs(roots, R)
        
        # Check: exists s with Re(r_s) < 0?
        if not any(np.real(rs) < 0):
            violations += 1
    
    print(f"d={d:>3}: {violations}/{total} cases where NO s has Re(r_s) < 0")

print("\n" + "=" * 70)
print("PROOF SKETCH for d ODD")
print("=" * 70)

print("""
THEOREM: For d odd, there exists s with Re(r_s) < 0.

PROOF:
  ∏ r_s = (-1)^d · C_d = -C_d  (since d odd, (-1)^d = -1).
  
  Since |C_d| > 0 and C_d is close to 1: ∏ r_s ≈ -1.
  So arg(∏ r_s) ≈ π (modulo 2π).
  
  If Re(r_s) > 0 for ALL s, then arg(r_s) ∈ (-π/2, π/2).
  Sum of arguments ∈ (-dπ/2, dπ/2).
  
  We need: ∑ arg(r_s) ≡ π (mod 2π).
  
  The smallest positive value of ∑ arg ≡ π (mod 2π) in (-dπ/2, dπ/2) is just π.
  Since π ∈ (-dπ/2, dπ/2) for d ≥ 1, we cannot yet derive a contradiction
  from the argument alone!
  
  BUT: we also have the MODULUS constraint.
  If all Re(r_s) > 0 with |∏ r_s| ≈ 1, the r_s are constrained to a narrow
  cone. For the product to have arg ≈ π, at least one factor must contribute
  significant negative imaginary part...
  
  Actually, the clean argument is SIMPLER:
""")

print("=" * 70)
print("CLEAN ARGUMENT (using the product being real-negative)")
print("=" * 70)

print("""
For d ODD and R large enough (so |C_d - 1| < ε):
  ∏ r_s ≈ -1   (a negative real number)

If ALL r_s had Re(r_s) > 0, write r_s = |r_s|·e^{iφ_s} with |φ_s| < π/2.
Then ∏ r_s = (∏|r_s|) · e^{i·∑φ_s}.

For ∏r_s to be close to -1, we need:
  ∑ φ_s ≡ π  (mod 2π)

But each |φ_s| < π/2, and there are d terms.
Question: can ∑ φ_s = π with d odd and |φ_s| < π/2?

Yes! For example, d=3: φ = (π/3, π/3, π/3) gives ∑ = π. POSSIBLE!

So the clean argument FAILS for the universal statement. BUT:
we need ∏|r_s| ≈ 1 simultaneously. This extra constraint is severe.
""")

# Let's check: for NearDegen, what's the product?
print("=" * 70)
print("PRODUCT ANALYSIS for NearDegen")
print("=" * 70)

for d in [5, 7, 10, 20]:
    roots = 0.99 * np.ones(d) + 0.01j * np.linspace(-1, 1, d)
    R = 2.0
    rs = compute_all_rs(roots, R)
    
    prod_mod = abs(np.prod(rs))
    prod_arg = np.degrees(np.angle(np.prod(rs)))
    
    mods = abs(rs)
    args = np.degrees(np.angle(rs))
    
    print(f"\nd={d}:")
    print(f"  |∏r_s| = {prod_mod:.6f}")
    print(f"  arg(∏r_s) = {prod_arg:.2f}°")
    print(f"  Individual |r_s|: min={min(mods):.4f}, max={max(mods):.4f}, mean={np.mean(mods):.4f}")
    print(f"  Individual arg(r_s): [{min(args):.1f}°, {max(args):.1f}°]")
    print(f"  Positive-Re starts: Re(r) = {', '.join(f'{x:.3f}' for x in np.real(rs) if x > 0)}")

print("\n" + "=" * 70)
print("REFINED SEARCH: What is the TIGHTEST correct statement?")
print("=" * 70)

# Test: "at least d/2 starts have Re < 0"
# Test: "at least 1 start has Re < 0"  
# Test: "the best start has Re(r_s) < -c for some universal c"

np.random.seed(123)
for d in [3, 5, 7, 10, 20, 50, 100]:
    min_neg_fraction = 1.0
    min_best_re = 0.0
    
    for trial in range(2000):
        roots = np.random.rand(d) * np.exp(2j * np.pi * np.random.rand(d))
        roots /= np.maximum(abs(roots), 1)
        R = 2.0
        
        rs = compute_all_rs(roots, R)
        
        neg_frac = np.sum(np.real(rs) < 0) / d
        best_re = np.min(np.real(rs))
        
        min_neg_fraction = min(min_neg_fraction, neg_frac)
        min_best_re = min(min_best_re, best_re)
    
    # Also test the adversarial NearDegen
    for spread in [0.01, 0.05, 0.1, 0.5]:
        for center_angle in np.linspace(0, 2*np.pi, 20):
            center = 0.99 * np.exp(1j * center_angle)
            roots_adv = center + spread * (np.random.randn(d) + 1j*np.random.randn(d)) / np.sqrt(d)
            roots_adv /= np.maximum(abs(roots_adv), 1)
            
            rs = compute_all_rs(roots_adv, R)
            
            neg_frac = np.sum(np.real(rs) < 0) / d
            best_re = np.min(np.real(rs))
            
            min_neg_fraction = min(min_neg_fraction, neg_frac)
    
    print(f"d={d:>3}: worst fraction with Re<0 = {min_neg_fraction:.2f}, "
          f"best min(Re(r_s)) = {min_best_re:.4f}")
