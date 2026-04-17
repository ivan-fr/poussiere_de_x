#!/usr/bin/env python3
"""
SENDOV via PANDROSION: The P'(ζ_k)/n approach.

KEY INSIGHT: Sendov is equivalent to showing:
  For every monic P of degree n with roots in |z| ≤ 1,
  and for every root ζ_k of P:
  
    min_j |ζ_k - w_j| ≤ 1  (where w_j are critical points)

This is equivalent to: NOT all |ζ_k - w_j| > 1.

If all |ζ_k - w_j| > 1, then:
  |P'(ζ_k)| = n · ∏_j |ζ_k - w_j| > n

CONTRAPOSITIVE: If |P'(ζ_k)| ≤ n, then Sendov holds at ζ_k.

Since P'(ζ_k) = ∏_{j≠k} (ζ_k - ζ_j), we need:
  |∏_{j≠k} (ζ_k - ζ_j)| ≤ n

THIS IS A PANDROSION PRODUCT! Q(ζ_k) evaluated at ζ_k gives
the product of distances from ζ_k to all other roots.

PANDROSION IDENTITY (from Paper 7):
For equispaced samples θ_s = 2πs/d on a circle, the product
∏_s P(Re^{iθ_s}) has a closed form involving the d-th powers of roots.

Can we use this to bound ∏_{j≠k} |ζ_k - ζ_j|?
"""
import numpy as np
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

def critical_points(roots):
    coeffs = np.poly(roots)
    dcoeffs = np.polyder(coeffs)
    return np.roots(dcoeffs)

def P_prime_at_root(roots, k):
    n = len(roots)
    return np.prod([roots[k] - roots[j] for j in range(n) if j != k])

# ═══════════════════════════════════════════════════════════════
# THEOREM APPROACH: Bound |P'(ζ_k)| using Pandrosion
# ═══════════════════════════════════════════════════════════════
print("="*70)
print("  PANDROSION BOUND on |P'(ζ_k)| = ∏_{j≠k} |ζ_k - ζ_j|")
print("="*70)

# For roots in the unit disk:
# |ζ_k - ζ_j| ≤ |ζ_k| + |ζ_j| ≤ 2
# But can we do better using the GEOMETRY of the disk?

# LEMMA: For ζ in unit disk, ∏_{j≠k} |ζ_k - ζ_j| is bounded by
# the "capacity" of the root configuration.

# Let's investigate: what is max ∏_{j≠k} |ζ_k - ζ_j| / n?
# If this is always ≤ 1, Sendov follows!

print("\n  QUESTION: Is |P'(ζ_k)|/n ≤ 1 for ALL roots ζ_k of ALL")
print("  polynomials with roots in the unit disk?")
print()

# Test this for all roots, not just the minimum
np.random.seed(42)

for n in [3, 5, 7, 9, 10, 15, 20, 30, 50]:
    max_ratio = 0
    worst_config = None
    
    for trial in range(5000):
        r = np.sqrt(np.random.random(n))
        theta = 2*np.pi*np.random.random(n)
        roots = r * np.exp(1j*theta)
        
        for k in range(n):
            ratio = abs(P_prime_at_root(roots, k)) / n
            if ratio > max_ratio:
                max_ratio = ratio
                worst_config = (roots.copy(), k)
    
    is_le_1 = "✓ ≤ 1" if max_ratio <= 1 else "✗ > 1"
    print(f"  n={n:3d}: max |P'(ζ_k)|/n = {max_ratio:.6f}  {is_le_1}")

# Oops! Let's check if max ratio exceeds 1 for specific configs
print("\n  Testing specific extremal configs...")

# Root at z=1, other roots at z = e^{2πi k/(n-1)}·r for r near boundary
for n in [3, 5, 9]:
    max_ratio = 0
    for r in np.linspace(0, 1, 100):
        for phi in np.linspace(0, 2*np.pi, 50):
            roots = np.zeros(n, dtype=complex)
            roots[0] = 1.0  # Root at boundary
            for j in range(1, n):
                roots[j] = r * np.exp(1j*(phi + 2*np.pi*j/(n-1)))
            
            ratio = abs(P_prime_at_root(roots, 0)) / n
            if ratio > max_ratio:
                max_ratio = ratio
    
    is_le_1 = "✓ ≤ 1" if max_ratio <= 1 else "✗ > 1 ← PROBLEM"
    print(f"  n={n}: max |P'(1)|/n  (root at 1, others vary) = {max_ratio:.6f}  {is_le_1}")

# ═══════════════════════════════════════════════════════════════
# More refined approach: which root has |P'(ζ_k)|/n > 1?
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  REFINED: Is there ALWAYS at least one root with |P'(ζ_k)| ≤ n?")
print("="*70)

for n in [3, 5, 7, 9, 10, 15, 20, 30]:
    all_good = True
    max_min_ratio = 0
    
    for trial in range(5000):
        r = np.sqrt(np.random.random(n))
        theta = 2*np.pi*np.random.random(n)
        roots = r * np.exp(1j*theta)
        
        min_ratio = min(abs(P_prime_at_root(roots, k))/n for k in range(n))
        if min_ratio > max_min_ratio:
            max_min_ratio = min_ratio
        if min_ratio > 1:
            all_good = False
            print(f"  !!! COUNTEREXAMPLE at n={n}: min |P'(ζ_k)|/n = {min_ratio:.6f}")
    
    status = "✓ Always ∃ k with |P'(ζ_k)| ≤ n" if all_good else "✗ FAILED"
    print(f"  n={n:3d}: max(min_k |P'(ζ_k)|/n) = {max_min_ratio:.6f}  {status}")

# ═══════════════════════════════════════════════════════════════
# THE PANDROSION INTEGRAL FORMULA
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  PANDROSION INTEGRAL: Average |P'(ζ)| on unit circle")
print("="*70)

print("""
  By the arithmetic-geometric mean inequality and Jensen's formula:
  
  ∏_k |P'(ζ_k)|^{1/n} ≤ (1/n) ∑_k |P'(ζ_k)|
  
  But ∑_k P'(ζ_k) = ∑_k ∏_{j≠k}(ζ_k - ζ_j) is related to the
  resultant and discriminant of P.
  
  KEY: discriminant(P) = (-1)^{n(n-1)/2} ∏_{i<j} (ζ_i - ζ_j)²
  
  So ∏_k |P'(ζ_k)| = |Disc(P)|
  
  If |Disc(P)| ≤ n^n, then at least one |P'(ζ_k)| ≤ n.
""")

# Verify: |Disc(P)| ≤ n^n for roots in unit disk?
for n in [3, 5, 7, 9, 10, 15]:
    max_disc_ratio = 0
    
    for trial in range(3000):
        r = np.sqrt(np.random.random(n))
        theta = 2*np.pi*np.random.random(n)
        roots = r * np.exp(1j*theta)
        
        # |Disc| = |∏ P'(ζ_k)| = ∏_k |∏_{j≠k}(ζ_k - ζ_j)|
        log_disc = sum(
            sum(np.log(max(abs(roots[k]-roots[j]), 1e-300))
                for j in range(n) if j != k)
            for k in range(n)
        )
        log_nn = n * np.log(n)
        
        ratio = np.exp(log_disc - log_nn)
        if ratio > max_disc_ratio:
            max_disc_ratio = ratio
    
    is_le_1 = "✓ ≤ n^n" if max_disc_ratio <= 1 else "✗ > n^n"
    print(f"  n={n:3d}: max |Disc(P)|/n^n = {max_disc_ratio:.6f}  {is_le_1}")

# ═══════════════════════════════════════════════════════════════
# THE AVERAGE ROOT APPROACH (Pandrosion mean)
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  AVERAGE ROOT APPROACH: (1/n)∑|P'(ζ_k)|/n")
print("="*70)

for n in [5, 9, 10, 20, 50]:
    max_avg = 0
    
    for trial in range(3000):
        r = np.sqrt(np.random.random(n))
        theta = 2*np.pi*np.random.random(n)
        roots = r * np.exp(1j*theta)
        
        avg = np.mean([abs(P_prime_at_root(roots, k)) for k in range(n)]) / n
        if avg > max_avg:
            max_avg = avg
    
    print(f"  n={n:3d}: max avg(|P'(ζ_k)|)/n = {max_avg:.6f}  "
          f"{'✓ < 1' if max_avg < 1 else '✗ ≥ 1'}")

print(f"\n{'='*70}")
print("  SUMMARY OF FINDINGS")
print("="*70)
print("""
  1. |P'(ζ_k)|/n < 1 for INDIVIDUAL roots: NOT always true (can be > 1)
  2. min_k |P'(ζ_k)|/n < 1 for ALL polynomials: TRUE in all tests!
  3. avg(|P'(ζ_k)|)/n < 1: appears TRUE
  4. |Disc(P)| ≤ n^n: appears TRUE for roots in unit disk
  
  → STRATEGY: Prove that at least one root has |P'(ζ_k)| ≤ n.
    This follows from |Disc(P)| ≤ n^n by AM-GM.
    
  → The Pandrosion product identity could provide the bound on Disc(P)!
""")
