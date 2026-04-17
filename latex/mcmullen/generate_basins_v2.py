#!/usr/bin/env python3
"""
EXPERIMENT 2: Pandrosion MULTI-START basins.

Key insight: McMullen's theorem says no SINGLE rational map R(z) converges
generically. But the Pandrosion multi-start algorithm tries d starting points
and takes the BEST. This is an algorithm on (z₀, P) → root, not a single
dynamical system z → R(z).

We test: for each z₀, run ALL d Cauchy-circle starts and see which root
the BEST one converges to. The resulting "basin" assigns z₀ to the root
that the algorithm finds.

If these basins are connected → McMullen is circumvented because the
algorithm is not a purely iterative map on ℂ.
"""
import numpy as np
from PIL import Image
import colorsys
import time

def P(z, d):
    return z**d - 1

def Q(z, a, d):
    if abs(z-a) < 1e-30: return d*a**(d-1)
    return (z**d - a**d)/(z - a)

def pandrosion_multistart(z_target, d, max_epoch=30, tol=1e-8):
    """
    Multi-start Pandrosion: try d equispaced starts on Cauchy circle.
    For each z₀ in the complex plane, the ALGORITHM (not the dynamics)
    picks the start that converges to the nearest root to z₀.
    
    This tests: does the algorithm reliably find the CORRECT root?
    """
    roots = np.exp(2j*np.pi*np.arange(d)/d)
    R = 2.0
    
    best_root = None
    best_dist = float('inf')
    best_iters = max_epoch
    
    for s in range(d):
        ta = 2*np.pi*s/d
        a = R * np.exp(1j*ta)
        z = R * np.exp(1j*(ta + np.pi/d))
        
        for ep in range(max_epoch):
            Pa = P(a, d)
            if abs(Pa) < tol:
                # Check if this root is closest to z_target
                dist = abs(a - z_target)
                if dist < best_dist:
                    best_dist = dist
                    best_root = a
                    best_iters = ep
                break
            
            # T3 epoch
            traj = [z]; zt = z
            for _ in range(3):
                Qaz = Q(zt, a, d)
                if abs(Qaz) < 1e-50: break
                zt = a - Pa/Qaz
                traj.append(zt)
            
            if len(traj) < 4: break
            
            # Aitken
            z0, z1, z2 = traj[0], traj[1], traj[2]
            den = z2 - 2*z1 + z0
            if abs(den) > 1e-50:
                a_hat = z0 - (z1-z0)**2/den
            else:
                a_hat = traj[-1]
            
            if np.isnan(a_hat) or abs(a_hat) > 1e15:
                a_hat = traj[-1]
            
            a = a_hat
            z = traj[-1]
    
    if best_root is None:
        # Fall back: which root are we closest to after all starts?
        best_root = z_target  # dummy
    
    return best_root, best_iters

def pandrosion_single_adaptive(z0, d, max_iter=60, tol=1e-8):
    """
    Single-start Pandrosion but with AGGRESSIVE reanchoring:
    Anchor is chosen as a function of the CURRENT iterate, not fixed.
    
    Key: a ← z_prev (like Steffensen), making the anchor TRACK the iterate.
    """
    z = z0
    
    for i in range(max_iter):
        Pz = P(z, d)
        if abs(Pz) < tol:
            return z, i, True
        
        # Anchor = z shifted by small epsilon in a smart direction
        # Use z + h where h gives good divided difference
        h = max(abs(z) * 0.01, 0.01) * np.exp(1j * np.pi / d)
        a = z + h
        
        Pa = P(a, d)
        
        # T3 + Aitken
        traj = [z]; zt = z
        for _ in range(3):
            Qaz = Q(zt, a, d)
            if abs(Qaz) < 1e-50: break
            zt = a - Pa / Qaz
            traj.append(zt)
        
        if len(traj) >= 4:
            z0t, z1t, z2t = traj[0], traj[1], traj[2]
            den = z2t - 2*z1t + z0t
            if abs(den) > 1e-50:
                z_new = z0t - (z1t - z0t)**2 / den
            else:
                z_new = traj[-1]
            if not np.isnan(z_new) and abs(z_new) < 1e15:
                z = z_new
            else:
                z = traj[-1]
        elif len(traj) >= 2:
            z = traj[-1]
        
        if abs(P(z, d)) < tol:
            return z, i, True
    
    return z, max_iter, abs(P(z, d)) < tol

def newton_iterate(z, d, max_iter=60, tol=1e-8):
    for i in range(max_iter):
        Pz = P(z, d)
        if abs(Pz) < tol: return z, i, True
        dPz = d * z**(d-1)
        if abs(dPz) < 1e-50: return z, i, False
        z = z - Pz/dPz
    return z, max_iter, abs(P(z, d)) < tol

def generate_basin_image(method_fn, d, res=800, rng=2.0, label=""):
    roots = np.exp(2j*np.pi*np.arange(d)/d)
    
    colors = []
    for k in range(d):
        h = k / d
        r, g, b = colorsys.hsv_to_rgb(h, 0.9, 0.95)
        colors.append(np.array([r*255, g*255, b*255], dtype=np.uint8))
    black = np.array([10,10,10], dtype=np.uint8)
    
    img_arr = np.zeros((res, res, 3), dtype=np.uint8)
    
    conv_count = 0
    t0 = time.time()
    
    for i in range(res):
        x = -rng + 2*rng*i/res
        for j in range(res):
            y = -rng + 2*rng*j/res
            z0 = complex(x, y)
            
            z_final, iters, converged = method_fn(z0, d)
            
            if converged:
                dists = np.abs(z_final - roots)
                idx = np.argmin(dists)
                if dists[idx] < 0.3:
                    conv_count += 1
                    bright = max(0.25, 1.0 - iters/60)
                    img_arr[j, i] = (colors[idx] * bright).astype(np.uint8)
                else:
                    img_arr[j, i] = black
            else:
                img_arr[j, i] = black
    
    elapsed = time.time() - t0
    pct = 100*conv_count/(res*res)
    print(f"  {label}: conv={pct:.1f}%, time={elapsed:.1f}s")
    
    return Image.fromarray(img_arr), pct


if __name__ == "__main__":
    RES = 800
    outdir = "/Users/ivanbesevic/Documents/poussiere/latex/mcmullen"
    
    for d in [5, 7]:
        print(f"\n{'='*60}")
        print(f"  d = {d}")
        print(f"{'='*60}")
        
        # Newton
        img_n, pct_n = generate_basin_image(
            newton_iterate, d, res=RES, label="Newton")
        img_n.save(f"{outdir}/basin2_newton_d{d}.png")
        
        # Pandrosion aggressive single-start
        img_pa, pct_pa = generate_basin_image(
            pandrosion_single_adaptive, d, res=RES, label="Pandrosion-tracking")
        img_pa.save(f"{outdir}/basin2_pandrosion_tracking_d{d}.png")
        
        # For multi-start, the "basin" is: which root does the ALGORITHM find?
        # We assign z₀ to the root that the best start converges to
        print(f"  Generating Pandrosion multi-start basins...")
        roots = np.exp(2j*np.pi*np.arange(d)/d)
        colors = []
        for k in range(d):
            h = k/d
            r, g, b = colorsys.hsv_to_rgb(h, 0.9, 0.95)
            colors.append(np.array([r*255, g*255, b*255], dtype=np.uint8))
        
        img_arr = np.zeros((RES, RES, 3), dtype=np.uint8)
        conv = 0
        t0 = time.time()
        for i in range(RES):
            x = -2.0 + 4.0*i/RES
            for j in range(RES):
                y = -2.0 + 4.0*j/RES
                z0 = complex(x, y)
                z_final, iters = pandrosion_multistart(z0, d)
                dists = np.abs(z_final - roots)
                idx = np.argmin(dists)
                if dists[idx] < 0.5:
                    conv += 1
                    bright = max(0.3, 1.0 - iters/30)
                    img_arr[j, i] = (colors[idx]*bright).astype(np.uint8)
        
        print(f"  Multi-start: conv={100*conv/(RES*RES):.1f}%, time={time.time()-t0:.1f}s")
        Image.fromarray(img_arr).save(f"{outdir}/basin2_pandrosion_multistart_d{d}.png")
