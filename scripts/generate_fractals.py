#!/usr/bin/env python3
"""
Generate stunning fractal images for Pure Pandrosion T3/T4 with iterated scaling.
Compare with Newton basins. Show the speed of convergence.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os

OUTPUT_DIR = "/Users/ivanbesevic/Documents/poussiere/latex/smale/figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def eval_P_roots(z, roots):
    """Evaluate P(z) = prod(z - roots)"""
    result = np.ones_like(z, dtype=complex)
    for r in roots:
        result *= (z - r)
    return result

def pandrosion_step_grid(z, a, P_a, roots):
    """Vectorized Pandrosion step F_a(z) for a grid.
    a is the anchor (scalar), z is the grid.
    """
    P_z = eval_P_roots(z, roots)
    Q = np.where(np.abs(z - a) > 1e-30,
                 (P_z - P_a) / (z - a),
                 np.zeros_like(z))  # avoid division by zero
    # Where Q is too small, mark as failed
    mask_ok = np.abs(Q) > 1e-50
    result = np.where(mask_ok, a - P_a / Q, z)
    return result, mask_ok

def newton_step_grid(z, roots):
    """Vectorized Newton step for a grid."""
    P_z = eval_P_roots(z, roots)
    P_prime = np.zeros_like(z, dtype=complex)
    for k, r in enumerate(roots):
        prod_others = np.ones_like(z, dtype=complex)
        for j, r2 in enumerate(roots):
            if j != k:
                prod_others *= (z - r2)
        P_prime += prod_others
    mask_ok = np.abs(P_prime) > 1e-50
    result = np.where(mask_ok, z - P_z / P_prime, z)
    return result, mask_ok

def compute_basins(roots, method="pandrosion_T3", 
                   xrange=(-2.5, 2.5), yrange=(-2.5, 2.5), 
                   resolution=800, max_iter=50):
    """Compute convergence basins on a grid."""
    d = len(roots)
    x = np.linspace(xrange[0], xrange[1], resolution)
    y = np.linspace(yrange[0], yrange[1], resolution)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    
    # Which root each point converges to (-1 = not converged)
    root_index = -np.ones(Z.shape, dtype=int)
    # Number of iterations to converge
    iter_count = max_iter * np.ones(Z.shape, dtype=int)
    # Still active mask
    active = np.ones(Z.shape, dtype=bool)
    
    if method.startswith("pandrosion"):
        rho = np.max(np.abs(roots))
        R = 1 + rho
        # Anchor: fixed on Cauchy circle (for grid visualization)
        # For iterated scaling: anchor adapts per pixel
        if "iterated" in method:
            # Each pixel has its own anchor
            A = Z.copy()  # initial anchor = the point itself... no
            # Better: anchor at a nearby Cauchy-circle point
            angles = np.angle(Z)
            A = R * np.exp(1j * angles)
            # But then A could equal Z. Offset slightly:
            A = R * np.exp(1j * (angles + np.pi/(2*d)))
        else:
            A = R * np.exp(1j * 0.3) * np.ones_like(Z)
        
        z_iter = Z.copy()
        a_iter = A.copy()
        
        for iteration in range(max_iter):
            if not np.any(active):
                break
            
            P_a = eval_P_roots(a_iter[active], roots)
            
            # 3 or 4 base-map steps
            n_steps = 4 if "T4" in method else 3
            z_prev = z_iter[active].copy()
            z_cur = z_iter[active].copy()
            z_traj = [z_cur.copy()]
            
            ok = np.ones(z_cur.shape, dtype=bool)
            for step in range(n_steps):
                z_new, step_ok = pandrosion_step_grid(z_cur, a_iter[active], P_a, roots)
                ok &= step_ok
                z_cur = np.where(ok, z_new, z_cur)
                z_traj.append(z_cur.copy())
            
            # Aitken acceleration
            if len(z_traj) >= 3:
                z0, z1, z2 = z_traj[0], z_traj[1], z_traj[2]
                denom = z2 - 2*z1 + z0
                aitken_ok = np.abs(denom) > 1e-50
                z_hat = np.where(aitken_ok,
                                z0 - (z1-z0)**2 / np.where(aitken_ok, denom, 1),
                                z2)
            else:
                z_hat = z_cur
            
            # Update: iterated scaling
            # New anchor = Aitken result, new iterate = last base-map
            new_a = z_hat
            new_z = z_cur
            # Ensure a ≠ z
            too_close = np.abs(new_a - new_z) < 1e-30
            new_z = np.where(too_close, z_traj[-2] if len(z_traj) > 2 else z_prev, new_z)
            
            a_iter[active] = new_a
            z_iter[active] = new_z
            
            # Check convergence
            for k, r in enumerate(roots):
                converged = active & (np.abs(a_iter - r) < 1e-6)
                root_index[converged] = k
                iter_count[converged] = np.minimum(iter_count[converged], iteration + 1)
                active[converged] = False
    
    elif method == "newton":
        z_iter = Z.copy()
        for iteration in range(max_iter):
            if not np.any(active):
                break
            z_new, ok = newton_step_grid(z_iter, roots)
            z_iter = np.where(active[:,:,None] if z_iter.ndim > 2 else active, z_new, z_iter)
            
            for k, r in enumerate(roots):
                converged = active & (np.abs(z_iter - r) < 1e-6)
                root_index[converged] = k
                iter_count[converged] = np.minimum(iter_count[converged], iteration + 1)
                active[converged] = False
    
    return root_index, iter_count, (x, y)

def make_colormap(d):
    """Create a beautiful colormap for d roots."""
    base_colors = [
        '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7',
        '#DDA0DD', '#98D8C8', '#F7DC6F', '#BB8FCE', '#85C1E9',
        '#F8C471', '#82E0AA', '#F1948A', '#AED6F1', '#D5DBDB',
        '#FADBD8', '#D4EFDF', '#FCF3CF', '#D6EAF8', '#F5EEF8',
    ]
    return base_colors[:d]

def plot_basin(root_index, iter_count, xy, roots, title, filename, 
               show_roots=True, speed_mode=False):
    """Plot a beautiful basin of attraction."""
    x, y = xy
    d = len(roots)
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=150)
    
    if speed_mode:
        # Color by number of iterations (speed map)
        cmap = plt.cm.magma_r
        im = ax.imshow(iter_count, extent=[x[0], x[-1], y[0], y[-1]],
                       origin='lower', cmap=cmap, vmin=1, vmax=np.percentile(iter_count, 95))
        plt.colorbar(im, ax=ax, label='Iterations', shrink=0.8)
    else:
        # Color by root (basin map)
        colors = make_colormap(d)
        # Create a colored image
        img = np.zeros((*root_index.shape, 3))
        
        for k in range(d):
            mask = root_index == k
            c = matplotlib.colors.to_rgb(colors[k])
            # Darken by iteration count for depth effect
            brightness = np.clip(1.0 - 0.6 * iter_count / np.max(iter_count), 0.15, 1.0)
            for ch in range(3):
                img[:,:,ch] += mask * c[ch] * brightness
        
        # Not converged: black
        not_conv = root_index == -1
        img[not_conv] = [0.05, 0.05, 0.08]
        
        ax.imshow(img, extent=[x[0], x[-1], y[0], y[-1]], origin='lower')
    
    if show_roots:
        for k, r in enumerate(roots):
            ax.plot(r.real, r.imag, 'w*', markersize=8, markeredgecolor='black', 
                    markeredgewidth=0.5)
    
    ax.set_title(title, fontsize=14, fontweight='bold', pad=10, color='white')
    ax.set_xlabel('Re(z)', fontsize=11, color='white')
    ax.set_ylabel('Im(z)', fontsize=11, color='white')
    ax.tick_params(colors='white')
    fig.patch.set_facecolor('#0a0a0f')
    ax.set_facecolor('#0a0a0f')
    for spine in ax.spines.values():
        spine.set_color('#333')
    
    plt.tight_layout()
    filepath = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(filepath, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Saved: {filepath}")
    return filepath

def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  GENERATING PANDROSION FRACTAL IMAGES                              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    
    RES = 600  # resolution
    
    # ═══════════════════════════════════════════════════════════════════
    # 1. z^5 - 1: Pandrosion T3 basins
    # ═══════════════════════════════════════════════════════════════════
    d = 5
    roots = np.exp(2j * np.pi * np.arange(d) / d)
    
    print(f"\n  Generating z^{d}-1 basins...")
    ri, ic, xy = compute_basins(roots, "pandrosion_T3_iterated", resolution=RES, max_iter=30)
    plot_basin(ri, ic, xy, roots,
               f"Pure Pandrosion T3 — z⁵−1\nIterated Scaling (derivative-free)",
               "pandrosion_T3_z5.png")
    
    ri, ic, xy = compute_basins(roots, "pandrosion_T4_iterated", resolution=RES, max_iter=25)
    plot_basin(ri, ic, xy, roots,
               f"Pure Pandrosion T4 — z⁵−1\nIterated Scaling (derivative-free)",
               "pandrosion_T4_z5.png")
    
    ri, ic, xy = compute_basins(roots, "newton", resolution=RES, max_iter=30)
    plot_basin(ri, ic, xy, roots,
               f"Newton's Method — z⁵−1",
               "newton_z5.png")
    
    # Speed maps
    ri, ic, xy = compute_basins(roots, "pandrosion_T3_iterated", resolution=RES, max_iter=30)
    plot_basin(ri, ic, xy, roots,
               f"Pandrosion T3 Speed — z⁵−1",
               "pandrosion_T3_z5_speed.png", speed_mode=True)
    
    # ═══════════════════════════════════════════════════════════════════
    # 2. z^8 - 1: more complex symmetry
    # ═══════════════════════════════════════════════════════════════════
    d = 8
    roots = np.exp(2j * np.pi * np.arange(d) / d)
    
    print(f"\n  Generating z^{d}-1 basins...")
    ri, ic, xy = compute_basins(roots, "pandrosion_T3_iterated", resolution=RES, max_iter=35)
    plot_basin(ri, ic, xy, roots,
               f"Pure Pandrosion T3 — z⁸−1",
               "pandrosion_T3_z8.png")
    
    ri, ic, xy = compute_basins(roots, "pandrosion_T4_iterated", resolution=RES, max_iter=30)
    plot_basin(ri, ic, xy, roots,
               f"Pure Pandrosion T4 — z⁸−1",
               "pandrosion_T4_z8.png")
    
    ri, ic, xy = compute_basins(roots, "newton", resolution=RES, max_iter=35)
    plot_basin(ri, ic, xy, roots,
               f"Newton's Method — z⁸−1",
               "newton_z8.png")
    
    # ═══════════════════════════════════════════════════════════════════
    # 3. z^12 - 1: high degree
    # ═══════════════════════════════════════════════════════════════════
    d = 12
    roots = np.exp(2j * np.pi * np.arange(d) / d)
    
    print(f"\n  Generating z^{d}-1 basins...")
    ri, ic, xy = compute_basins(roots, "pandrosion_T4_iterated", resolution=RES, max_iter=40)
    plot_basin(ri, ic, xy, roots,
               f"Pure Pandrosion T4 — z¹²−1",
               "pandrosion_T4_z12.png")
    
    ri, ic, xy = compute_basins(roots, "pandrosion_T4_iterated", resolution=RES, max_iter=40)
    plot_basin(ri, ic, xy, roots,
               f"Pandrosion T4 Speed — z¹²−1",
               "pandrosion_T4_z12_speed.png", speed_mode=True)
    
    # ═══════════════════════════════════════════════════════════════════
    # 4. z^3 - 2z + 2: non-symmetric roots
    # ═══════════════════════════════════════════════════════════════════
    roots_cubic = np.roots([1, 0, -2, 2])
    
    print(f"\n  Generating z³-2z+2 basins...")
    ri, ic, xy = compute_basins(roots_cubic, "pandrosion_T3_iterated", 
                                xrange=(-3, 3), yrange=(-3, 3), resolution=RES, max_iter=30)
    plot_basin(ri, ic, xy, roots_cubic,
               f"Pure Pandrosion T3 — z³−2z+2",
               "pandrosion_T3_cubic.png")
    
    ri, ic, xy = compute_basins(roots_cubic, "newton",
                                xrange=(-3, 3), yrange=(-3, 3), resolution=RES, max_iter=30)
    plot_basin(ri, ic, xy, roots_cubic,
               f"Newton — z³−2z+2",
               "newton_cubic.png")
    
    # ═══════════════════════════════════════════════════════════════════
    # 5. z^4 + z^3 + z^2 + z + 1 (5th cyclotomic)
    # ═══════════════════════════════════════════════════════════════════
    roots_cyc = np.exp(2j * np.pi * np.arange(1, 5) / 5)
    
    print(f"\n  Generating 5th cyclotomic basins...")
    ri, ic, xy = compute_basins(roots_cyc, "pandrosion_T4_iterated",
                                xrange=(-2, 2), yrange=(-2, 2), resolution=RES, max_iter=30)
    plot_basin(ri, ic, xy, roots_cyc,
               f"Pandrosion T4 — Φ₅(z)",
               "pandrosion_T4_cyclotomic5.png")
    
    # ═══════════════════════════════════════════════════════════════════
    # 6. Clustered roots: z(z-0.1)(z-0.2)(z-1)(z+1)
    # ═══════════════════════════════════════════════════════════════════
    roots_cluster = np.array([0, 0.1, 0.2, 1, -1], dtype=complex)
    
    print(f"\n  Generating clustered roots basins...")
    ri, ic, xy = compute_basins(roots_cluster, "pandrosion_T4_iterated",
                                xrange=(-2, 2), yrange=(-2, 2), resolution=RES, max_iter=35)
    plot_basin(ri, ic, xy, roots_cluster,
               f"Pandrosion T4 — Clustered roots",
               "pandrosion_T4_clustered.png")
    
    ri, ic, xy = compute_basins(roots_cluster, "newton",
                                xrange=(-2, 2), yrange=(-2, 2), resolution=RES, max_iter=35)
    plot_basin(ri, ic, xy, roots_cluster,
               f"Newton — Clustered roots",
               "newton_clustered.png")
    
    # ═══════════════════════════════════════════════════════════════════
    # 7. High-degree: z^20 - 1 (T4 speed)
    # ═══════════════════════════════════════════════════════════════════
    d = 20
    roots20 = np.exp(2j * np.pi * np.arange(d) / d)
    
    print(f"\n  Generating z^{d}-1 basins...")
    ri, ic, xy = compute_basins(roots20, "pandrosion_T4_iterated",
                                xrange=(-2.5, 2.5), yrange=(-2.5, 2.5), 
                                resolution=RES, max_iter=50)
    plot_basin(ri, ic, xy, roots20,
               f"Pandrosion T4 — z²⁰−1",
               "pandrosion_T4_z20.png")
    
    ri, ic, xy = compute_basins(roots20, "pandrosion_T4_iterated",
                                xrange=(-2.5, 2.5), yrange=(-2.5, 2.5),
                                resolution=RES, max_iter=50)
    plot_basin(ri, ic, xy, roots20,
               f"Pandrosion T4 Speed — z²⁰−1",
               "pandrosion_T4_z20_speed.png", speed_mode=True)
    
    # ═══════════════════════════════════════════════════════════════════
    # 8. Zoomed fractal boundary
    # ═══════════════════════════════════════════════════════════════════
    d = 5
    roots5 = np.exp(2j * np.pi * np.arange(d) / d)
    
    print(f"\n  Generating zoomed boundary...")
    # Zoom near the boundary between two basins
    ri, ic, xy = compute_basins(roots5, "pandrosion_T3_iterated",
                                xrange=(0.0, 1.0), yrange=(0.5, 1.5),
                                resolution=RES, max_iter=50)
    plot_basin(ri, ic, xy, roots5,
               f"Pandrosion T3 — z⁵−1 (zoom)",
               "pandrosion_T3_z5_zoom.png", show_roots=False)
    
    ri, ic, xy = compute_basins(roots5, "newton",
                                xrange=(0.0, 1.0), yrange=(0.5, 1.5),
                                resolution=RES, max_iter=50)
    plot_basin(ri, ic, xy, roots5,
               f"Newton — z⁵−1 (zoom)",
               "newton_z5_zoom.png", show_roots=False)
    
    print(f"\n  All images saved to {OUTPUT_DIR}/")
    print(f"  Total: {len(os.listdir(OUTPUT_DIR))} images")


if __name__ == "__main__":
    main()
