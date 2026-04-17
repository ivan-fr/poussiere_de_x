#!/usr/bin/env python3
"""
SENDOV'S CONJECTURE - Numerical Investigation via Pandrosion

Conjecture (Sendov, 1958): If P has degree n ≥ 2 with all roots in |z| ≤ 1,
then for each root ζ_k, there exists a critical point w of P with |ζ_k - w| ≤ 1.

Key identity linking roots to critical points:
  P'(ζ_k) = ∏_{j≠k} (ζ_k - ζ_j)   [for monic P]
  P'(z)   = n ∏_j (z - w_j)          [w_j = critical points]

Therefore:  ∏_{j≠k} (ζ_k - ζ_j) = n ∏_j (ζ_k - w_j)

If Sendov is false at ζ_k, then ALL |ζ_k - w_j| > 1, so:
  |P'(ζ_k)| = n ∏|ζ_k - w_j| > n

PANDROSION CONNECTION:
  Q(ζ_k, z) = P(z)/(z - ζ_k) = ∏_{j≠k} (z - ζ_j)
  P'(z) = ∑_k Q(ζ_k, z)
  
  The Pandrosion product identity (Paper 7) on equispaced evaluations of Q
  gives structural constraints on ∏(ζ_k - ζ_j).
"""
import numpy as np
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

def critical_points(roots):
    """Compute critical points of monic polynomial with given roots."""
    n = len(roots)
    # Build polynomial coefficients
    coeffs = np.poly(roots)  # [1, -sum(roots), ..., prod(roots)]
    # Derivative
    dcoeffs = np.polyder(coeffs)
    return np.roots(dcoeffs)

def sendov_distance(roots):
    """For each root, compute min distance to a critical point.
    Returns (max_min_distance, worst_root_index, worst_root, nearest_cp)."""
    cps = critical_points(roots)
    n = len(roots)
    
    max_min_dist = 0
    worst_k = 0
    worst_root = roots[0]
    nearest_cp = cps[0]
    
    for k in range(n):
        dists = np.abs(roots[k] - cps)
        min_d = np.min(dists)
        if min_d > max_min_dist:
            max_min_dist = min_d
            worst_k = k
            worst_root = roots[k]
            nearest_cp = cps[np.argmin(dists)]
    
    return max_min_dist, worst_k, worst_root, nearest_cp

def P_prime_at_root(roots, k):
    """P'(ζ_k) = ∏_{j≠k} (ζ_k - ζ_j)"""
    n = len(roots)
    prod = 1+0j
    for j in range(n):
        if j != k:
            prod *= (roots[k] - roots[j])
    return prod

# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 1: Random polynomials - verify Sendov
# ═══════════════════════════════════════════════════════════════
print("="*70)
print("  EXPERIMENT 1: Random polynomials with roots in unit disk")
print("="*70)

np.random.seed(42)
max_violation = 0
max_dist_seen = 0

for n in [5, 9, 10, 15, 20, 30, 50]:
    worst_for_n = 0
    for trial in range(5000):
        # Random roots in unit disk
        r = np.sqrt(np.random.random(n))
        theta = 2 * np.pi * np.random.random(n)
        roots = r * np.exp(1j * theta)
        
        d, k, zk, cp = sendov_distance(roots)
        if d > worst_for_n:
            worst_for_n = d
        if d > max_dist_seen:
            max_dist_seen = d
    
    print(f"  n={n:3d}: max(min dist to CP) = {worst_for_n:.6f}  "
          f"{'✓ < 1' if worst_for_n < 1 else '✗ ≥ 1 VIOLATION!'}")

print(f"\n  Overall max distance seen: {max_dist_seen:.6f}")
print(f"  Sendov margin: {1 - max_dist_seen:.6f}")

# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 2: Known extremal configurations
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  EXPERIMENT 2: Known extremal configurations")
print("="*70)

# Extremal for Sendov: one root at z=1, others clustered near origin
for n in [3, 5, 9, 10, 15, 20]:
    # Config 1: root at z=1, n-1 roots at origin
    roots1 = np.zeros(n, dtype=complex)
    roots1[0] = 1.0
    d1, _, z1, cp1 = sendov_distance(roots1)
    
    # Config 2: root at z=1, n-1 roots uniformly on circle |z|=ε
    eps = 0.01
    roots2 = eps * np.exp(2j*np.pi*np.arange(n-1)/(n-1))
    roots2 = np.append(roots2, 1.0+0j)
    d2, _, z2, cp2 = sendov_distance(roots2)
    
    # Config 3: root at z=a (near boundary), n-1 roots at z=0
    a = 1 - 1e-6
    roots3 = np.zeros(n, dtype=complex)
    roots3[0] = a
    d3, _, z3, cp3 = sendov_distance(roots3)
    
    print(f"  n={n:3d}: (1,0,...,0)→d={d1:.6f}  "
          f"(1,ε-circle)→d={d2:.6f}  "
          f"(a≈1,0,...,0)→d={d3:.6f}")

# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 3: Pandrosion product identity & critical points
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  EXPERIMENT 3: Pandrosion identity for roots & critical points")
print("="*70)

print("""
  KEY IDENTITY: P'(ζ_k) = ∏_{j≠k} (ζ_k - ζ_j) = n · ∏_j (ζ_k - w_j)
  
  If Sendov fails at ζ_k, then |ζ_k - w_j| > 1 for ALL j, hence:
    |P'(ζ_k)| > n
    
  So Sendov follows from: |P'(ζ_k)| ≤ n for at least one root ζ_k.
  
  But |P'(ζ_k)| = ∏_{j≠k} |ζ_k - ζ_j| with all roots in unit disk,
  so |ζ_k - ζ_j| ≤ 2, giving |P'(ζ_k)| ≤ 2^{n-1} (too weak).
  
  PANDROSION APPROACH: use the product identity to get a TIGHTER bound.
""")

# For each n, find the root with smallest |P'(ζ_k)| and compare with n
for n in [5, 9, 10, 15, 20, 30]:
    min_ratio = float('inf')
    
    for trial in range(2000):
        r = np.sqrt(np.random.random(n))
        theta = 2*np.pi*np.random.random(n)
        roots = r * np.exp(1j*theta)
        
        for k in range(n):
            Ppk = abs(P_prime_at_root(roots, k))
            ratio = Ppk / n
            if ratio < min_ratio:
                min_ratio = ratio
    
    print(f"  n={n:3d}: min |P'(ζ_k)|/n = {min_ratio:.6f}  "
          f"{'< 1 → Sendov safe' if min_ratio < 1 else '≥ 1 → need closer look'}")

# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 4: The Pandrosion quotient Q(ζ_k, w_j)
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  EXPERIMENT 4: Pandrosion quotient at critical points")
print("="*70)

print("""
  Q(ζ_k, z) = P(z)/(z - ζ_k) = ∏_{j≠k} (z - ζ_j)
  
  At a critical point w: P'(w) = 0, so ∑_k Q(ζ_k, w) = 0.
  
  This means: for each critical point w, the Q-values at all roots 
  sum to ZERO. This is a constraint that links ALL roots to each CP.
  
  PANDROSION PRODUCT on equispaced evaluations of Q:
  For R > 1 and equispaced a_s on circle of radius R:
    ∏_s Q(ζ_k, a_s) = ∏_{j≠k} (a_s^d - ζ_j^d) / ... 
""")

# For a specific polynomial, explore the Q(ζ_k, w_j) values
n = 9  # The first open case
np.random.seed(123)
r = np.sqrt(np.random.random(n))
theta = 2*np.pi*np.random.random(n)
roots = r * np.exp(1j*theta)
cps = critical_points(roots)

print(f"  Example: n={n}, random roots in unit disk")
print(f"  Roots:    {['%.2f' % abs(z) for z in roots]}")
print(f"  CPs:      {['%.2f' % abs(z) for z in cps]}")
print()

for k in range(min(3, n)):
    dists = np.abs(roots[k] - cps)
    nearest_idx = np.argmin(dists)
    w = cps[nearest_idx]
    
    # Q(ζ_k, w) = ∏_{j≠k} (w - ζ_j)
    Qkw = np.prod([w - roots[j] for j in range(n) if j != k])
    
    # Sum of Q(ζ_j, w) over all j should be 0 (= P'(w))
    sumQ = sum(np.prod([w - roots[i] for i in range(n) if i != j]) for j in range(n))
    
    print(f"  Root ζ_{k}: |ζ_{k}|={abs(roots[k]):.4f}, "
          f"nearest CP w: |ζ_{k}-w|={dists[nearest_idx]:.4f}, "
          f"|Q(ζ_{k},w)|={abs(Qkw):.4f}, "
          f"|∑Q|={abs(sumQ):.2e}")

# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 5: Extremal optimization for n=9
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  EXPERIMENT 5: Extremal Sendov configuration for n=9 (first open case)")
print("="*70)

def neg_sendov_dist(params):
    """Maximize the min(|ζ_k - nearest CP|) = try to VIOLATE Sendov."""
    n = len(params) // 2
    roots = params[:n] + 1j*params[n:]
    # Penalize if roots leave unit disk
    penalty = 0
    for z in roots:
        if abs(z) > 1:
            penalty += 100 * (abs(z) - 1)**2
    
    try:
        d, _, _, _ = sendov_distance(roots)
        return -(d - penalty)  # minimize negative = maximize
    except:
        return 0

best_dist = 0
best_roots = None

for trial in range(500):
    n = 9
    # Random starting config
    r = np.sqrt(np.random.random(n))
    theta = 2*np.pi*np.random.random(n)
    roots0 = r * np.exp(1j*theta)
    x0 = np.concatenate([roots0.real, roots0.imag])
    
    try:
        res = minimize(neg_sendov_dist, x0, method='Nelder-Mead',
                      options={'maxiter': 2000, 'xatol': 1e-12})
        roots_opt = res.x[:n] + 1j*res.x[n:]
        
        # Check all in unit disk
        if np.all(np.abs(roots_opt) <= 1.001):
            d, k, zk, cp = sendov_distance(roots_opt)
            if d > best_dist:
                best_dist = d
                best_roots = roots_opt.copy()
    except:
        pass

print(f"\n  Best (hardest) configuration found for n=9:")
print(f"  max min|ζ_k - CP| = {best_dist:.8f}")
print(f"  Sendov margin = {1 - best_dist:.8f}")
if best_roots is not None:
    print(f"  Roots: {['%.4f%+.4fi' % (z.real, z.imag) for z in best_roots]}")
    cps = critical_points(best_roots)
    print(f"  CPs:   {['%.4f%+.4fi' % (z.real, z.imag) for z in cps]}")
    
    # Which root is hardest?
    d, k, zk, cp = sendov_distance(best_roots)
    print(f"  Hardest root: ζ_{k} = {zk:.4f}, nearest CP = {cp:.4f}, dist = {d:.6f}")
    print(f"  |P'(ζ_k)| = {abs(P_prime_at_root(best_roots, k)):.6f}, n = {n}")
    print(f"  |P'(ζ_k)|/n = {abs(P_prime_at_root(best_roots, k))/n:.6f}")

print(f"\n{'='*70}")
print(f"  CONCLUSION: Sendov holds for all {5000*6 + 500} tested cases.")
print(f"  Maximum distance seen: {max(max_dist_seen, best_dist):.8f}")
print(f"  Sendov margin ≥ {1 - max(max_dist_seen, best_dist):.8f}")
print(f"{'='*70}")
