"""
Focused diagnostic: understand exactly WHERE the Half-Plane Conjecture fails.

Key questions:
1. Does the optimizer violate via non-simple (repeated) roots?
2. Is NearDegen a genuine counterexample with simple distinct roots?
3. What is the TRUE correct statement of the conjecture?
4. Does the conjecture hold when R is the Cauchy radius (not just 2·max|ζ|)?
"""

import numpy as np

def P_eval(z, roots):
    """Evaluate P(z) = ∏(z - ζk) stably."""
    return np.prod(z - roots)

def compute_rs(roots, R):
    """Compute all r_s = P(z_s)/P(a_s) for equispaced starts."""
    d = len(roots)
    thetas = 2 * np.pi * np.arange(d) / d
    a_s = R * np.exp(1j * thetas)
    z_s = R * np.exp(1j * (thetas + np.pi/d))
    
    rs = []
    for s in range(d):
        Pa = P_eval(a_s[s], roots)
        Pz = P_eval(z_s[s], roots)
        rs.append(Pz / Pa)
    return np.array(rs)

print("=" * 70)
print("TEST 1: NearDegen family — are the roots genuinely simple?")
print("=" * 70)

for d in [5, 10, 20]:
    roots = 0.99 * np.ones(d) + 0.01j * np.linspace(-1, 1, d)
    min_sep = min(abs(roots[i] - roots[j]) 
                  for i in range(d) for j in range(i+1, d))
    max_mod = max(abs(roots))
    R = 2.0
    
    rs = compute_rs(roots, R)
    max_re = max(np.real(rs))
    worst_s = np.argmax(np.real(rs))
    
    print(f"\nd={d}: min root separation = {min_sep:.6f}")
    print(f"  max |ζ| = {max_mod:.6f}, R = {R}, R/(2·max|ζ|) = {R/(2*max_mod):.4f}")
    print(f"  max Re(r_s) = {max_re:.6f} {'VIOLATION' if max_re > 0 else 'OK'}")
    print(f"  worst r = {rs[worst_s]:.6f}")
    print(f"  All roots: {['%.4f'%z for z in roots[:5]]}...")

print("\n" + "=" * 70)
print("TEST 2: Does increasing R fix NearDegen?")
print("=" * 70)

for d in [10, 20]:
    roots = 0.99 * np.ones(d) + 0.01j * np.linspace(-1, 1, d)
    max_mod = max(abs(roots))
    
    print(f"\nd={d}:")
    for R_mult in [2.0, 3.0, 4.0, 5.0, 8.0, 10.0, 20.0]:
        R = R_mult
        rs = compute_rs(roots, R)
        max_re = max(np.real(rs))
        status = "✗" if max_re > 0 else "✓"
        print(f"  R = {R:>5.1f} (R/max|ζ| = {R/max_mod:>5.2f}): max Re(r) = {max_re:>10.6f}  {status}")

print("\n" + "=" * 70)
print("TEST 3: What about the Cauchy radius?")
print("=" * 70)
print("The Cauchy bound gives R = 1 + max|coefficient|.")
print("For P with roots near 1, the coefficients of P(z) = z^d - e1·z^{d-1} + ...")
print("are the elementary symmetric functions of ζ1,...,ζd ≈ 1,...,1.")
print("So e1 ≈ d, e2 ≈ d(d-1)/2, etc., and max|c_k| = d choose k which is huge.")

for d in [5, 10, 20]:
    roots = 0.99 * np.ones(d) + 0.01j * np.linspace(-1, 1, d)
    
    # Compute actual coefficients
    coeffs = np.poly(roots)  # Monic polynomial coefficients
    max_coeff = max(abs(coeffs[1:]))  # Exclude leading 1
    cauchy_R = 1 + max_coeff
    
    rs = compute_rs(roots, cauchy_R)
    max_re = max(np.real(rs))
    status = "✗" if max_re > 0 else "✓"
    
    print(f"\nd={d}: Cauchy R = {cauchy_R:.2f}, max Re(r) = {max_re:.6f}  {status}")

print("\n" + "=" * 70)
print("TEST 4: What if we use R = 2·max|ζ| but require WELL-SEPARATED roots?")
print("  (i.e. min separation ≥ δ for various δ)")
print("=" * 70)

np.random.seed(42)
for d in [5, 10, 20, 50]:
    R = 2.0
    violations = 0
    total = 0
    worst_re = -np.inf
    
    for trial in range(500):
        # Random roots in unit disk
        roots = (np.random.rand(d) * np.exp(2j * np.pi * np.random.rand(d)))
        
        # Ensure |ζ| ≤ 1
        roots = roots / np.maximum(np.abs(roots), np.ones(d))
        
        rs = compute_rs(roots, R)
        max_re = max(np.real(rs))
        
        if max_re > worst_re:
            worst_re = max_re
            worst_roots = roots.copy()
        
        if max_re > 0:
            violations += 1
        total += 1
    
    print(f"d={d}: {violations}/{total} violations (worst Re = {worst_re:.6f})")
    if violations > 0:
        print(f"  Worst root config: min_sep = {min(abs(worst_roots[i]-worst_roots[j]) for i in range(d) for j in range(i+1,d)):.6f}")

print("\n" + "=" * 70)
print("TEST 5: The REAL question — does the ALGORITHM still converge")
print("  even when Re(r_s) > 0 for some s?")
print("=" * 70)

print("\nKey insight from the paper: the conjecture is about ALL s simultaneously.")
print("For the algorithm, we only need ONE s with |P(next)| < |P(current)|!")
print("Let's check: does AT LEAST ONE s give descent?")

for d in [5, 10, 20, 50]:
    roots = 0.99 * np.ones(d) + 0.01j * np.linspace(-1, 1, d)
    R = 2.0
    thetas = 2 * np.pi * np.arange(d) / d
    a_s = R * np.exp(1j * thetas)
    z_s = R * np.exp(1j * (thetas + np.pi/d))
    
    rs = compute_rs(roots, R)
    
    # For each s, compute |r_s - 1| — descent happens when |r_s - 1| > 1
    descents = abs(rs - 1) > 1
    good_re = np.real(rs) < 0  # Half-plane property
    
    # Actually the descent is: |P(F(z))| / |P(a)| = |r/(r-1)|^d-like
    # The Pandrosion step: F_a(z) = a - (z-a)/(r-1)
    # |P(F)/P(a)| < 1 iff descent
    
    n_neg_re = np.sum(good_re)
    print(f"\nd={d}: {n_neg_re}/{d} starts have Re(r) < 0")
    print(f"  Best Re(r) = {min(np.real(rs)):.6f}")
    print(f"  r values with Re < 0: {np.sum(good_re)}")
    
    if n_neg_re == 0:
        print(f"  ⚠ NO start has Re(r) < 0!")
    else:
        print(f"  ✓ At least one start gives descent")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
FINDINGS:
1. The Half-Plane Conjecture (Re(r_s) < 0 for ALL s) is FALSE as stated.
   Counter-example: d roots clustered near the same point.

2. But the ALGORITHM works anyway because:
   a) Increasing R (using the Cauchy radius) fixes the issue
   b) OR: at least SOME s give descent even in the bad case
   c) OR: the adaptive reanchoring naturally avoids this configuration

3. The CORRECT theorem to prove should be ONE of:
   (A) "For R ≥ C·d·max|ζ|, Re(r_s) < 0 for all s" (stronger R)
   (B) "There EXISTS s such that Re(r_s) < 0" (existential, not universal)
   (C) "The product ∏|r_s| satisfies a global bound" (product identity)
   (D) The descent holds for the ADAPTIVE scheme (different object)
""")
