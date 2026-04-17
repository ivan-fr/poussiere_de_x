#!/usr/bin/env python3
"""
EXPERIMENT: Newton vs Pandrosion basins of attraction for z^d - 1.

Three methods compared:
  1. Newton: z ← z - P(z)/P'(z)
  2. Pandrosion fixed anchor: z ← F_a(z) with a = const (from Paper 2)
  3. Pandrosion adaptive (Steffensen): (a,z) ← (â, z₃)  with reanchoring

Key question: Does adaptive Pandrosion produce CONNECTED basins for d ≥ 4?
If yes → McMullen's impossibility theorem is circumvented.
"""
import numpy as np
from PIL import Image
import colorsys
import time

def roots_of_unity(d):
    return np.exp(2j * np.pi * np.arange(d) / d)

def P(z, d):
    return z**d - 1

def dP(z, d):
    return d * z**(d-1)

def Q_divided_diff(z, a, d):
    """(P(z) - P(a)) / (z - a) for z^d - 1"""
    if abs(z - a) < 1e-30:
        return dP(a, d)
    return (z**d - a**d) / (z - a)

# ═══════════════════════════════════════════════════════
# METHOD 1: Newton
# ═══════════════════════════════════════════════════════
def newton_iterate(z, d, max_iter=50, tol=1e-8):
    for i in range(max_iter):
        Pz = P(z, d)
        dPz = dP(z, d)
        if abs(dPz) < 1e-50:
            return z, i, False
        z = z - Pz / dPz
        if abs(P(z, d)) < tol:
            return z, i, True
    return z, max_iter, False

# ═══════════════════════════════════════════════════════
# METHOD 2: Pandrosion fixed anchor
# ═══════════════════════════════════════════════════════
def pandrosion_fixed_iterate(z, a, d, max_iter=50, tol=1e-8):
    for i in range(max_iter):
        Qaz = Q_divided_diff(z, a, d)
        if abs(Qaz) < 1e-50:
            return z, i, False
        Pa = P(a, d)
        z_new = a - Pa / Qaz
        z = z_new
        if abs(P(z, d)) < tol:
            return z, i, True
    return z, max_iter, False

# ═══════════════════════════════════════════════════════
# METHOD 3: Pandrosion ADAPTIVE (Steffensen reanchoring)
# ═══════════════════════════════════════════════════════
def pandrosion_adaptive_iterate(z0, d, max_iter=50, tol=1e-8):
    """
    The key innovation: after each T3 epoch (3 steps + Aitken),
    we REANCHOR: a ← â, z ← z₃.
    This breaks the holomorphic dynamics on ℂ and operates on ℂ².
    """
    # Initial anchor: use z0 itself shifted slightly
    R = max(2.0, 2*abs(z0))
    # Use the nearest Cauchy circle point as anchor
    if abs(z0) < 1e-10:
        a = R
    else:
        a = R * z0 / abs(z0)
    z = z0
    
    for epoch in range(max_iter):
        Pa = P(a, d)
        if abs(Pa) < tol:
            return a, epoch, True
        
        # T3: three Pandrosion steps
        traj = [z]
        zt = z
        for step in range(3):
            Qaz = Q_divided_diff(zt, a, d)
            if abs(Qaz) < 1e-50:
                break
            zt = a - Pa / Qaz
            traj.append(zt)
        
        if len(traj) < 4:
            return zt, epoch, abs(P(zt, d)) < tol
        
        # Aitken Δ² acceleration
        z0t, z1t, z2t = traj[0], traj[1], traj[2]
        den = z2t - 2*z1t + z0t
        if abs(den) > 1e-50:
            a_hat = z0t - (z1t - z0t)**2 / den
        else:
            a_hat = traj[-1]
        
        if np.isnan(a_hat) or abs(a_hat) > 1e15:
            a_hat = traj[-1]
        
        if abs(P(a_hat, d)) < tol:
            return a_hat, epoch, True
        
        # REANCHOR: this is what breaks McMullen
        a = a_hat
        z = traj[-1]
    
    return a, max_iter, abs(P(a, d)) < tol

# ═══════════════════════════════════════════════════════
# BASIN GENERATION
# ═══════════════════════════════════════════════════════
def classify_root(z_final, roots, tol=0.3):
    """Return index of closest root, or -1 if not converged."""
    dists = np.abs(z_final - roots)
    idx = np.argmin(dists)
    if dists[idx] < tol:
        return idx
    return -1

def generate_basins(method_fn, d, res=800, xrange=(-2,2), yrange=(-2,2), 
                    max_iter=60, label=""):
    roots = roots_of_unity(d)
    
    # Color palette: one hue per root
    colors = []
    for k in range(d):
        h = k / d
        r, g, b = colorsys.hsv_to_rgb(h, 0.85, 0.95)
        colors.append((int(r*255), int(g*255), int(b*255)))
    black = (15, 15, 15)
    
    img = Image.new('RGB', (res, res))
    pixels = img.load()
    
    convergence_count = 0
    boundary_count = 0
    total = res * res
    
    t0 = time.time()
    
    for i in range(res):
        x = xrange[0] + (xrange[1] - xrange[0]) * i / res
        for j in range(res):
            y = yrange[0] + (yrange[1] - yrange[0]) * j / res
            z0 = complex(x, y)
            
            z_final, iters, converged = method_fn(z0, d, max_iter=max_iter)
            
            if converged:
                root_idx = classify_root(z_final, roots)
                if root_idx >= 0:
                    convergence_count += 1
                    # Shade by iteration count  
                    brightness = max(0.3, 1 - iters / max_iter)
                    c = colors[root_idx]
                    pixels[i, j] = (int(c[0]*brightness), int(c[1]*brightness), int(c[2]*brightness))
                else:
                    pixels[i, j] = black
                    boundary_count += 1
            else:
                pixels[i, j] = black
                boundary_count += 1
    
    elapsed = time.time() - t0
    conv_pct = 100 * convergence_count / total
    
    print(f"  {label}: {conv_pct:.1f}% converged, "
          f"boundary pixels: {boundary_count} ({100*boundary_count/total:.1f}%), "
          f"time: {elapsed:.1f}s")
    
    return img, conv_pct, boundary_count

def measure_fractal_dimension(img, threshold=30):
    """Estimate boundary fractal dimension via box-counting."""
    arr = np.array(img)
    gray = arr[:,:,0].astype(float) + arr[:,:,1].astype(float) + arr[:,:,2].astype(float)
    
    # Detect boundaries: pixels where neighbors have different root color
    boundary = np.zeros_like(gray, dtype=bool)
    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:
            if di == 0 and dj == 0:
                continue
            shifted = np.roll(np.roll(gray, di, axis=0), dj, axis=1)
            boundary |= (np.abs(gray - shifted) > threshold)
    
    boundary_count = np.sum(boundary)
    total = gray.shape[0] * gray.shape[1]
    
    # Box-counting at different scales
    sizes = [2, 4, 8, 16, 32, 64, 128]
    counts = []
    for s in sizes:
        count = 0
        for i in range(0, boundary.shape[0], s):
            for j in range(0, boundary.shape[1], s):
                if np.any(boundary[i:i+s, j:j+s]):
                    count += 1
        counts.append(count)
    
    # Linear regression on log-log
    if len(sizes) > 2 and all(c > 0 for c in counts):
        log_s = np.log(sizes)
        log_c = np.log(counts)
        coeffs = np.polyfit(log_s, log_c, 1)
        dim = -coeffs[0]
    else:
        dim = 1.0
    
    return dim, boundary_count


if __name__ == "__main__":
    RES = 600  # Resolution
    
    for d in [3, 4, 5, 7]:
        print(f"\n{'='*70}")
        print(f"  DEGREE d = {d}  (McMullen impossible for d ≥ 4)")
        print(f"{'='*70}")
        
        outdir = "/Users/ivanbesevic/Documents/poussiere/latex/mcmullen"
        
        # 1. Newton
        print(f"\n  Generating Newton basins...")
        img_n, conv_n, bd_n = generate_basins(
            newton_iterate, d, res=RES, label="Newton")
        img_n.save(f"{outdir}/basin_newton_d{d}.png")
        
        # 2. Pandrosion fixed anchor (a = 2)
        print(f"  Generating Pandrosion FIXED basins...")
        def pf(z0, d, max_iter=60):
            return pandrosion_fixed_iterate(z0, 2.0, d, max_iter=max_iter)
        img_pf, conv_pf, bd_pf = generate_basins(
            pf, d, res=RES, label="Pandrosion-fixed(a=2)")
        img_pf.save(f"{outdir}/basin_pandrosion_fixed_d{d}.png")
        
        # 3. Pandrosion ADAPTIVE
        print(f"  Generating Pandrosion ADAPTIVE basins...")
        img_pa, conv_pa, bd_pa = generate_basins(
            pandrosion_adaptive_iterate, d, res=RES, label="Pandrosion-ADAPTIVE")
        img_pa.save(f"{outdir}/basin_pandrosion_adaptive_d{d}.png")
        
        # Fractal dimension
        dim_n, _ = measure_fractal_dimension(img_n)
        dim_pf, _ = measure_fractal_dimension(img_pf)
        dim_pa, _ = measure_fractal_dimension(img_pa)
        
        print(f"\n  Fractal dimension (box-counting):")
        print(f"    Newton:              D = {dim_n:.3f}  {'FRACTAL' if dim_n > 1.1 else 'SMOOTH'}")
        print(f"    Pandrosion-fixed:    D = {dim_pf:.3f}  {'FRACTAL' if dim_pf > 1.1 else 'SMOOTH'}")
        print(f"    Pandrosion-adaptive: D = {dim_pa:.3f}  {'FRACTAL' if dim_pa > 1.1 else 'SMOOTH'}")
        
        print(f"\n  VERDICT for d={d}:")
        if dim_pa < 1.1 and dim_n > 1.1:
            print(f"    >>> PANDROSION ADAPTIVE ELIMINATES FRACTAL BOUNDARIES! <<<")
            print(f"    >>> McMullen's obstruction is CIRCUMVENTED! <<<")
        elif dim_pa < dim_n - 0.1:
            print(f"    >>> Pandrosion adaptive has SMOOTHER boundaries than Newton <<<")
        else:
            print(f"    >>> Both methods have similar boundary structure <<<")
