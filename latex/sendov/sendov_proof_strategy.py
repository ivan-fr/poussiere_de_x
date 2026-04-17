#!/usr/bin/env python3
"""
SENDOV PROOF STRATEGY via Pandrosion & Fekete-Szegő

THEOREM (Fekete, 1923): For n points z_1,...,z_n in the unit disk |z|≤1:
  ∏_{i<j} |z_i - z_j| ≤ n^{n/2}
  
Equality iff the z_i are n-th roots of unity (up to rotation).

CONSEQUENCE: |Disc(P)| = ∏_{i<j} |ζ_i - ζ_j|² ≤ n^n

SINCE: |Disc(P)| = ∏_k |P'(ζ_k)|,
  min_k |P'(ζ_k)| ≤ (∏_k |P'(ζ_k)|)^{1/n} ≤ n

BUT THIS IS NOT ENOUGH FOR SENDOV!
  
Sendov needs: for EACH root ζ_k, ∃ critical point w with |ζ_k - w| ≤ 1.
The discriminant bound only gives ∃ AT LEAST ONE root with |P'(ζ_k)| ≤ n.

Wait - let me reconsider. The implication is:
  |P'(ζ_k)| ≤ n → n · ∏_j |ζ_k - w_j| ≤ n → ∏_j |ζ_k - w_j| ≤ 1
  → at least one |ζ_k - w_j| ≤ 1 → Sendov holds at ζ_k

So the discriminant bound gives Sendov for at least ONE root.
For ALL roots, we need a different approach.

NEW APPROACH: The Pandrosion Quotient Identity

For EACH root ζ_k:
  P'(ζ_k) = ∏_{j≠k} (ζ_k - ζ_j)
  n·∏_j (ζ_k-w_j) = ∏_{j≠k} (ζ_k-ζ_j)

So: ∏_j |ζ_k-w_j| = (1/n)|∏_{j≠k}(ζ_k-ζ_j)|

SENDOV AT ζ_k ⟺ min_j|ζ_k-w_j| ≤ 1
  ⟸ ∏_j|ζ_k-w_j| ≤ 1  (since there are n-1 factors, all ≤ 2)
  ⟺ |P'(ζ_k)| ≤ n

So SENDOV at ζ_k ⟸ |P'(ζ_k)| ≤ n ⟺ |∏_{j≠k}(ζ_k-ζ_j)| ≤ n

QUESTION: Is |∏_{j≠k}(ζ_k-ζ_j)| ≤ n ALWAYS true for roots in unit disk?
"""
import numpy as np

# ═══════════════════════════════════════════════════════════════
# CRITICAL TEST: Is |P'(ζ_k)| ≤ n for ALL roots ζ_k, ALL polynomials?
# ═══════════════════════════════════════════════════════════════
print("="*70)
print("  CRITICAL: Is |P'(ζ_k)| ≤ n for ALL roots, ALL polynomials?")
print("="*70)
print("  (roots in unit disk)")

np.random.seed(42)

for n in [3, 5, 9, 10, 20, 50]:
    max_ratio = 0
    worst_config = None
    
    for trial in range(10000):
        r = np.sqrt(np.random.random(n))
        theta = 2*np.pi*np.random.random(n)
        roots = r * np.exp(1j*theta)
        
        for k in range(n):
            prod = np.prod([abs(roots[k]-roots[j]) for j in range(n) if j!=k])
            ratio = prod / n
            if ratio > max_ratio:
                max_ratio = ratio
                worst_config = (roots.copy(), k)
    
    print(f"  n={n:3d}: max |P'(ζ_k)|/n = {max_ratio:.4f}  "
          f"{'✓ ≤ 1' if max_ratio <= 1 else '✗ > 1'}")

print("\n  Answer: NO, |P'(ζ_k)|/n CAN exceed 1.")
print("  So the naive approach fails for individual roots.")

# ═══════════════════════════════════════════════════════════════
# ALTERNATIVE: Direct distance bound via Gauss-Lucas + Pandrosion
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  ALTERNATIVE APPROACH: Pandrosion iteration FROM ζ_k to find w")
print("="*70)

print("""
  NEW IDEA: Use the Pandrosion iteration itself to locate critical points!
  
  Given root ζ_k, we want to find w with P'(w) = 0 and |w - ζ_k| ≤ 1.
  
  Consider: P'(z) = ∑_j Q(ζ_j, z)
  
  If we apply the Pandrosion iteration to P' with anchor ζ_k:
    F_{ζ_k}(z) = ζ_k - P'(ζ_k)/Q'(ζ_k, z)
  where Q'(a,z) = (P'(z)-P'(a))/(z-a)
  
  This converges to a critical point, and the convergence region
  is controlled by |P'(ζ_k)| and the root geometry.
  
  KEY: If |ζ_k| is close to 1 (the hard case for Sendov),
  then P'(ζ_k) = ∏_{j≠k}(ζ_k-ζ_j) is SMALL because all roots
  are in the unit disk and clustered.
  
  The Pandrosion contraction ratio is r = P'(ζ_k)/... which is small!
""")

# ═══════════════════════════════════════════════════════════════
# TEST: Use Pandrosion to FIND the nearest critical point to ζ_k
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  Pandrosion iteration to find critical points near each root")
print("="*70)

def pandrosion_find_cp(roots, k, max_iter=50, tol=1e-10):
    """Use Pandrosion to find critical point nearest to root ζ_k.
    
    We solve P'(z) = 0 using the Pandrosion iteration:
    anchor = ζ_k, varied start z on a small circle around ζ_k.
    """
    n = len(roots)
    zk = roots[k]
    
    # P'(z) as a function
    def dP(z):
        return sum(np.prod([z-roots[j] for j in range(n) if j!=i]) for i in range(n))
    
    # Evaluate P' at anchor
    dPa = dP(zk)  # = ∏_{j≠k}(ζ_k - ζ_j)
    
    if abs(dPa) < tol:
        return zk, 0  # ζ_k is itself a critical point!
    
    # Start near ζ_k
    best_cp = None
    best_dist = float('inf')
    
    for offset_angle in np.linspace(0, 2*np.pi, n+1)[:-1]:
        # Small circle start
        z = zk + 0.5 * np.exp(1j * offset_angle)
        a = zk
        
        for it in range(max_iter):
            dPa = dP(a)
            if abs(dPa) < tol:
                dist = abs(a - zk)
                if dist < best_dist:
                    best_dist = dist
                    best_cp = a
                break
            
            dPz = dP(z)
            Qaz = (dPz - dPa)/(z-a) if abs(z-a) > 1e-30 else 1
            if abs(Qaz) < 1e-50:
                break
            
            z_new = a - dPa/Qaz
            
            # Aitken after 3 steps
            traj = [z, z_new]
            for _ in range(2):
                z_prev = traj[-1]
                dPz = dP(z_prev)
                Qaz = (dPz - dPa)/(z_prev-a) if abs(z_prev-a) > 1e-30 else 1
                if abs(Qaz) < 1e-50:
                    break
                traj.append(a - dPa/Qaz)
            
            if len(traj) >= 3:
                z0, z1, z2 = traj[0], traj[1], traj[2]
                den = z2 - 2*z1 + z0
                if abs(den) > 1e-50:
                    a_hat = z0 - (z1-z0)**2/den
                    if not np.isnan(a_hat) and abs(a_hat) < 100:
                        a = a_hat
                        z = traj[-1]
                        continue
            
            a = z_new if abs(dP(z_new)) < abs(dPa) else a
            z = z_new
        
        if abs(dP(a)) < tol:
            dist = abs(a - zk)
            if dist < best_dist:
                best_dist = dist
                best_cp = a
    
    return best_cp, best_dist

np.random.seed(42)
for n in [9, 10, 15, 20]:
    sendov_verified = 0
    sendov_total = 0
    max_pandrosion_dist = 0
    
    for trial in range(500):
        r = np.sqrt(np.random.random(n))
        theta = 2*np.pi*np.random.random(n)
        roots = r * np.exp(1j*theta)
        
        for k in range(n):
            sendov_total += 1
            cp, dist = pandrosion_find_cp(roots, k)
            if cp is not None and dist < 1:
                sendov_verified += 1
                if dist > max_pandrosion_dist:
                    max_pandrosion_dist = dist
    
    print(f"  n={n:3d}: Pandrosion found CP within distance 1 for {sendov_verified}/{sendov_total} "
          f"({100*sendov_verified/sendov_total:.1f}%), max dist = {max_pandrosion_dist:.4f}")

# ═══════════════════════════════════════════════════════════════
# THE KEY THEOREM IDEA
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("  THEOREM SKETCH: Pandrosion approach to Sendov")
print("="*70)
print("""
  THEOREM (Proposed): For every monic polynomial P of degree n ≥ 2
  with all roots in |z| ≤ 1, and for each root ζ_k of P:
  
    The Pandrosion iteration applied to P' with anchor ζ_k converges
    to a critical point w with |w - ζ_k| ≤ 1.
  
  PROOF STRATEGY:
  1. P'(ζ_k) = ∏_{j≠k}(ζ_k - ζ_j) is known.
  2. The Pandrosion contraction ratio r = P'(z)/Q(ζ_k, z) for P'
     satisfies |r| < 1 when z is in B(ζ_k, 1) ∩ D̄(0,1).
  3. By the product identity (Paper 7), the equispaced samples give
     a universal descent constant.
  4. The Gauss-Lucas theorem ensures all critical points are in the
     convex hull of the roots ⊆ D̄(0,1), so the search for w is
     confined to the unit disk.
  5. The Pandrosion iteration from ζ_k converges to a critical point
     within the contraction ball of radius ≤ 1.
     
  WHAT REMAINS TO PROVE:
  - The contraction ball radius is ≤ 1 for ALL configurations.
  - This requires bounding |P'(ζ_k)/P''(ζ_k)| ≤ 1 (Newton step size).
""")

# Test: |P'(ζ_k)/P''(ζ_k)| - this is the Newton step for P'
print(f"\n{'='*70}")
print("  Newton step for P': |P'(ζ_k)/P''(ζ_k)| ≤ 1?")
print("="*70)

for n in [3, 5, 9, 10, 20, 50]:
    max_ratio = 0
    
    for trial in range(5000):
        r = np.sqrt(np.random.random(n))
        theta = 2*np.pi*np.random.random(n)
        roots = r * np.exp(1j*theta)
        
        coeffs = np.poly(roots)
        dcoeffs = np.polyder(coeffs)
        ddcoeffs = np.polyder(dcoeffs)
        
        for k in range(n):
            dp = np.polyval(dcoeffs, roots[k])
            ddp = np.polyval(ddcoeffs, roots[k])
            if abs(ddp) > 1e-50:
                ratio = abs(dp/ddp)
                if ratio > max_ratio:
                    max_ratio = ratio
    
    print(f"  n={n:3d}: max |P'(ζ_k)/P''(ζ_k)| = {max_ratio:.6f}  "
          f"{'✓ ≤ 1' if max_ratio <= 1 else '✗ > 1'}")
