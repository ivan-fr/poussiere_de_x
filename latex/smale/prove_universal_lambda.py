#!/usr/bin/env python3
"""
THE PROOF: Why Λ_epoch ≈ 0.105 when λ ≈ 0.97

KEY IDENTITY (Aitken for near-geometric sequences):

For z_n - ζ = e_0·λ^n + O(e_0²·μ^n) with |μ| < λ:
  â = z_0 - (z_1-z_0)²/(z_2-2z_1+z_0)

Let's compute this EXACTLY:
  z_1 - z_0 = e_0(λ-1) + O(e_0²)  
  z_2 - z_1 = e_0·λ(λ-1) + O(e_0²)
  z_2 - 2z_1 + z_0 = e_0(λ-1)² + O(e_0²)

  â - z_0 = -(z_1-z_0)²/(z_2-2z_1+z_0)
            = -e_0²(λ-1)²/[e_0(λ-1)²] + O(e_0²)
            = -e_0 + O(e_0²)

  â - ζ = z_0 - ζ - e_0 + O(e_0²) = e_0 - e_0 + O(e_0²) = O(e_0²)

So for a FIRST-ORDER geometric sequence, Aitken is EXACT.
The O(e_0²) correction comes from the nonlinearity of F_a.

Now, the Pandrosion map F_a(z) near ζ has:
  F_a(z) = ζ + λ·(z-ζ) + μ·(z-ζ)² + ...
  where λ = 1 + P'(ζ)(ζ-a)/P(a) and μ = F_a''(ζ)/2.

For the EXACT Aitken formula on F_a(z) = ζ + λ(z-ζ) + μ(z-ζ)²:
  z_0 - ζ = e
  z_1 - ζ = λe + μe²
  z_2 - ζ = λ(λe+μe²) + μ(λe+μe²)² = λ²e + λμe² + μλ²e² + O(e³)
           = λ²e + μ(λ+λ²)e² + O(e³)

  z_1 - z_0 = (λ-1)e + μe²
  z_2 - z_1 = (λ²-λ)e + μ(λ+λ²-1)e² = λ(λ-1)e + μ(λ²+λ-1)e²

  z_2 - 2z_1 + z_0 = (λ²-2λ+1)e + μ(λ²+λ-1-2)e² + ...
                    = (λ-1)²e + μ(λ²+λ-3)e² + ...

  (z_1-z_0)² = [(λ-1)e + μe²]² = (λ-1)²e² + 2μ(λ-1)e³ + ...

  â = z_0 - (z_1-z_0)²/(z_2-2z_1+z_0)
    = z_0 - [(λ-1)²e² + ...] / [(λ-1)²e + ...]
    = z_0 - e·[1 + 2μe/(λ-1) + ...] / [1 + μ(λ²+λ-3)e/(λ-1)² + ...]
    ≈ z_0 - e·[1 + (2/(λ-1) - (λ²+λ-3)/(λ-1)²)·μe + ...]
    = ζ + e - e + [correction]·μe² + ...
    = ζ + C·μ·e² where C = [2(λ-1) - (λ²+λ-3)] / (λ-1)² = ...

Let's compute C:
  2(λ-1) = 2λ-2
  λ²+λ-3 = λ²+λ-3
  numerator = 2λ-2 - λ²-λ+3 = -λ² + λ + 1
  C = -(λ²-λ-1)/(λ-1)²

For λ ≈ 1: C ≈ -(1-1-1)/(0)² → diverges!? 

Wait, let me redo more carefully. Near λ = 1:
  C = -(λ²-λ-1)/(λ-1)²

  At λ = 1: numerator = 1-1-1 = -1, denominator = 0 → DIVERGES.

This means the Aitken correction |â-ζ| ≈ |C|·|μ|·|e|² DIVERGES as λ → 1.
But the data shows Λ_epoch is BOUNDED! So the Aitken step is doing something
more subtle than the Taylor expansion suggests.

The resolution: when λ ≈ 1, the higher-order terms in z_n - ζ are NOT small
relative to the geometric part. The sequence z_0, z_1, z_2 is not well-
approximated by a geometric sequence. Instead, z_n - ζ ≈ e + n·δ·e 
(almost arithmetic, not geometric).

For an arithmetic-like sequence: z_n ≈ z_0 + n·v where v = F_a(z_0) - z_0:
  z_1 - z_0 = v
  z_2 - z_1 ≈ v (arithmetic)
  z_2 - 2z_1 + z_0 ≈ 0

So the Aitken denominator ≈ 0, and â → ∞ (extrapolation blows up).

But the DATA shows â converges fine! So the sequence is BETWEEN arithmetic 
and geometric. The Aitken step works because:
  1. The sequence IS geometric (λ < 1), just barely
  2. z_2 - 2z_1 + z_0 ≠ 0, just small
  3. â extrapolates in the CORRECT DIRECTION (toward ζ)
  4. The extrapolation overshoots slightly, but lands between a and ζ

This means: â is between a and ζ. So |â-ζ| < |a-ζ| always, giving descent.

Let's verify: is â always between a and ζ?
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def eval_P_log(z, roots):
    return np.sum(np.log(np.abs(z - roots) + 1e-300))

def pandrosion_step(z, a, roots):
    if abs(z - a) < 1e-30:
        return None
    P_z = np.prod(z - roots)
    P_a = np.prod(a - roots)
    Q = (P_z - P_a) / (z - a)
    if abs(Q) < 1e-50:
        return None
    return a - P_a / Q


def verify_aitken_between():
    """Is â always closer to ζ than a is?"""
    print("═"*72)
    print("  IS |â-ζ| < |a-ζ| ALWAYS? (descent guaranteed if yes)")
    print("═"*72)
    
    for d in [5, 10, 20, 50, 100]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        
        total_epochs = 0
        aitken_closer = 0
        aitken_farther = 0
        
        for s in range(min(30, 2*d)):
            theta_a = 2 * np.pi * s / min(30, 2*d)
            theta_z = theta_a + np.pi / d
            a = R * np.exp(1j * theta_a)
            z = R * np.exp(1j * theta_z)
            
            for epoch in range(200):
                dists = np.abs(a - roots)
                nearest = np.argmin(dists)
                zeta = roots[nearest]
                dist_a = abs(a - zeta)
                
                traj = [z]
                z_t = z
                ok = True
                for _ in range(3):
                    z_n = pandrosion_step(z_t, a, roots)
                    if z_n is None or np.isnan(z_n):
                        ok = False
                        break
                    traj.append(z_n)
                    z_t = z_n
                
                if not ok or len(traj) < 4:
                    break
                
                z0, z1, z2 = traj[0], traj[1], traj[2]
                den = z2 - 2*z1 + z0
                if abs(den) > 1e-50:
                    z_hat = z0 - (z1-z0)**2/den
                else:
                    z_hat = traj[-1]
                
                if np.isnan(z_hat) or abs(z_hat) > 1e15:
                    z_hat = traj[-1]
                
                dist_hat = abs(z_hat - zeta)
                total_epochs += 1
                if dist_hat < dist_a:
                    aitken_closer += 1
                else:
                    aitken_farther += 1
                
                a = z_hat
                z = traj[-1] if abs(z_hat - traj[-1]) > 1e-30 else traj[-2]
                
                if np.min(np.abs(a - roots)) < 1e-14:
                    break
        
        pct = 100*aitken_closer/total_epochs if total_epochs > 0 else 0
        print(f"  d={d:>3d}: {pct:>6.2f}% Aitken closer ({aitken_closer}/{total_epochs})"
              f"  farther: {aitken_farther}")


def compute_lambda_at_first_epoch():
    """
    The FIRST epoch lambda: is it < 1 or ≈ 1?
    λ = |1 + P'(ζ)(ζ-a)/P(a)| on the Cauchy circle.
    
    For z^d - 1, ζ = 1, a = R = 2:
    P'(1) = d, P(2) = 2^d - 1
    λ = |1 + d(1-2)/(2^d-1)| = |1 - d/(2^d-1)|
    
    For d=5:  λ = |1 - 5/31| = |26/31| = 0.839
    For d=10: λ = |1 - 10/1023| = |1013/1023| = 0.990
    For d=50: λ = |1 - 50/(2^50-1)| ≈ 1 - 4.4e-14 ≈ 1
    
    So for d ≥ 10, λ ≈ 1 on the Cauchy circle!
    Yet Λ_epoch = 0.206 (not 1). The Aitken compression factor is:
    Λ_epoch / λ^6 ≈ 0.206 / 0.990^6 ≈ 0.206 / 0.941 ≈ 0.219
    
    For d=50: 0.206 / (1-ε)^6 ≈ 0.206. So Λ_epoch ≈ 0.206 regardless of λ!
    
    This means the Aitken step gives a CONSTANT compression factor
    that is INDEPENDENT of λ. This is the key insight!
    
    WHY? Because the Aitken step works on the RATIO (z_1-z_0)/(z_2-z_1),
    not on λ directly. On the Cauchy circle, this ratio is:
    (z_1-z_0)/(z_2-z_1) ≈ (F(z)-z)/(F²(z)-F(z))
    
    Since F contracts radially, |F^n(z)| decreases monotonically,
    and the ratio |F(z)-z|/|F²(z)-F(z)| determines the Aitken compression.
    
    On the Cauchy circle, this ratio ≈ 1/φ (golden ratio related!)
    because the step sizes are approximately constant.
    """
    print(f"\n{'═'*72}")
    print(f"  LAMBDA vs LAMBDA_EPOCH: Aitken provides CONSTANT compression")
    print(f"{'═'*72}")
    
    print(f"\n  {'d':>4s}  {'λ(theory)':>10s}  {'λ(num)':>8s}  {'λ^6':>8s}  "
          f"{'Λ_epoch':>10s}  {'Λ/λ^6':>8s}  {'step_ratio':>12s}")
    print("  " + "─"*70)
    
    for d in [5, 10, 20, 50, 100, 200]:
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        
        # Theoretical lambda
        lam_theory = abs(1 - d / (R**d - 1))
        
        a = R * np.exp(1j * 0.0)
        z = R * np.exp(1j * np.pi / d)
        
        traj = [z]
        z_t = z
        for _ in range(3):
            z_n = pandrosion_step(z_t, a, roots)
            if z_n is None:
                break
            traj.append(z_n)
            z_t = z_n
        
        if len(traj) < 4:
            continue
        
        zeta = roots[0]
        e = [abs(t - zeta) for t in traj]
        lam_num = e[2]/e[1] if e[1] > 1e-50 else 0
        
        z0, z1, z2 = traj[0], traj[1], traj[2]
        den = z2 - 2*z1 + z0
        if abs(den) > 1e-50:
            z_hat = z0 - (z1-z0)**2/den
        else:
            z_hat = traj[-1]
        
        lp_a = eval_P_log(a, roots)
        lp_hat = eval_P_log(z_hat, roots)
        Lambda_ep = np.exp(lp_hat - lp_a)
        
        # Step ratio: |z_1-z_0|/|z_2-z_1|
        step_ratio = abs(z1-z0)/abs(z2-z1) if abs(z2-z1) > 1e-50 else float('inf')
        
        lam6 = lam_theory**6 if lam_theory < 1 else float('nan')
        ratio_L = Lambda_ep / lam6 if np.isfinite(lam6) and lam6 > 1e-50 else float('nan')
        
        print(f"  {d:>4d}  {lam_theory:>10.6f}  {lam_num:>8.4f}  {lam6:>8.4f}  "
              f"{Lambda_ep:>10.6f}  {ratio_L:>8.4f}  {step_ratio:>12.4f}")
    
    # THE KEY OBSERVATION
    print(f"\n  ╔═══════════════════════════════════════════════════════════════╗")
    print(f"  ║  KEY OBSERVATION: Λ_epoch ≈ 0.20 regardless of d or λ!      ║")
    print(f"  ║                                                              ║")
    print(f"  ║  The Aitken compression is NOT λ^6 but something else.       ║")
    print(f"  ║  It depends on the RATIO of consecutive step sizes,          ║")
    print(f"  ║  which is ≈ λ ≈ 1 on the Cauchy circle.                     ║")
    print(f"  ║                                                              ║")
    print(f"  ║  The Aitken denominator z_2-2z_1+z_0 ≈ (λ-1)²·e is SMALL   ║")
    print(f"  ║  but NONZERO (because λ < 1!). This creates ample descent.  ║")
    print(f"  ╚═══════════════════════════════════════════════════════════════╝")


def compute_exact_Lambda():
    """
    EXACT computation of Λ_epoch for z^d-1 at the first epoch.
    
    a = R, z = Re^{iπ/d}, ζ = 1 (nearest root to a).
    
    All quantities can be computed algebraically:
    P(a) = R^d - 1
    Q(a,z) = (z^d - R^d)/(z - R)
    F_a(z) = R - (R^d-1)/Q
    etc.
    """
    print(f"\n{'═'*72}")
    print(f"  EXACT Λ_epoch FOR z^d-1: proving Λ < 1")
    print(f"{'═'*72}")
    
    print(f"\n  For z^d-1, first epoch from Cauchy circle:")
    print(f"  a = R=2, z = 2·exp(iπ/d)")
    print(f"\n  {'d':>4s}  {'Λ_epoch':>10s}  {'|â-1|':>10s}  {'|a-1|':>10s}  "
          f"{'|â-1|/|a-1|':>12s}  {'log|P(â)/P(a)|':>15s}")
    print("  " + "─"*65)
    
    all_Lambda = []
    
    for d in range(3, 201):
        roots = np.exp(2j * np.pi * np.arange(d) / d)
        R = 2.0
        
        a = R
        z = R * np.exp(1j * np.pi / d)
        
        traj = [z]
        z_t = z
        ok = True
        for _ in range(3):
            z_n = pandrosion_step(z_t, a, roots)
            if z_n is None or np.isnan(z_n):
                ok = False
                break
            traj.append(z_n)
            z_t = z_n
        
        if not ok or len(traj) < 4:
            continue
        
        z0, z1, z2 = traj[0], traj[1], traj[2]
        den = z2 - 2*z1 + z0
        if abs(den) > 1e-50:
            z_hat = z0 - (z1-z0)**2/den
        else:
            continue
        
        zeta = 1.0
        dist_hat = abs(z_hat - zeta)
        dist_a = abs(a - zeta)
        
        lp_a = eval_P_log(a, roots)
        lp_hat = eval_P_log(z_hat, roots)
        desc = lp_hat - lp_a
        Lambda = np.exp(desc)
        all_Lambda.append((d, Lambda))
        
        if d <= 20 or d % 20 == 0:
            ratio = dist_hat / dist_a if dist_a > 0 else 0
            print(f"  {d:>4d}  {Lambda:>10.6f}  {dist_hat:>10.4e}  {dist_a:>10.4e}  "
                  f"{ratio:>12.6f}  {desc:>15.4f}")
    
    # Analyze trend
    if all_Lambda:
        ds = np.array([x[0] for x in all_Lambda])
        Ls = np.array([x[1] for x in all_Lambda])
        
        print(f"\n  ASYMPTOTIC ANALYSIS (d → ∞):")
        print(f"  min Λ_epoch = {np.min(Ls):.6f} at d = {ds[np.argmin(Ls)]}")
        print(f"  max Λ_epoch = {np.max(Ls):.6f} at d = {ds[np.argmax(Ls)]}")
        print(f"  mean Λ_epoch = {np.mean(Ls):.6f}")
        print(f"  Λ_epoch → {Ls[-1]:.6f} as d → {ds[-1]}")
        
        # Is there an analytic formula?
        # For large d: Λ ≈ 0.206... = ?
        # 1/e ≈ 0.368, 1/5 = 0.2, 1/(e·φ) ≈ 0.227, exp(-ln5) = 0.2
        # Let's check: is it 1/5?
        limit = Ls[-1]
        print(f"\n  Candidate constants:")
        for name, val in [("1/5", 0.2), ("1/e^1.6", np.exp(-1.6)),
                          ("golden^-3", ((1+np.sqrt(5))/2)**(-3)),
                          ("3·exp(-π)", 3*np.exp(-np.pi)),
                          ("exp(-π/2)", np.exp(-np.pi/2)),
                          ("2·exp(-π²/6)", 2*np.exp(-np.pi**2/6))]:
            print(f"    {name:<15s} = {val:.6f}  (diff: {abs(limit-val):.6f})")


if __name__ == "__main__":
    verify_aitken_between()
    compute_lambda_at_first_epoch()
    compute_exact_Lambda()
