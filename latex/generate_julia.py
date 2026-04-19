"""
Generate stunning Julia Set fractal figures for the Pandrosion paper.
These show the classic Julia chaos that Pandrosion AVOIDS.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

DPI = 300
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 12

# ============================================================
# FIGURE A: Classic Julia Set for Newton z^3 - 1 (ZOOMED boundary)
# ============================================================
def julia_newton_zoom(resolution=1000):
    print("Generating Newton Julia zoom (boundary detail)...")
    # Zoom into the fractal boundary region
    cx, cy = 0.0, 0.0
    span = 0.6  # tight zoom
    x = np.linspace(cx - span, cx + span, resolution)
    y = np.linspace(cy - span, cy + span, resolution)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    
    roots = [1.0, np.exp(2j*np.pi/3), np.exp(-2j*np.pi/3)]
    
    convergence = np.zeros(Z.shape, dtype=int)
    speed = np.zeros(Z.shape, dtype=float)
    
    for iteration in range(80):
        mask = speed == 0
        if not np.any(mask):
            break
        Zm = Z[mask]
        denom = 3 * Zm**2
        safe = np.abs(denom) > 1e-12
        temp = np.copy(Zm)
        temp[safe] = Zm[safe] - (Zm[safe]**3 - 1) / denom[safe]
        Z_new = np.copy(Z)
        Z_new[mask] = temp
        Z = Z_new
        
        for idx, root in enumerate(roots):
            close = np.abs(Z - root) < 1e-8
            newly = close & (speed == 0)
            convergence[newly] = idx + 1
            speed[newly] = iteration + 1
    
    speed[speed == 0] = 80
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Rich colormapping per basin
    img = np.zeros((*Z.shape, 3))
    palettes = [
        # Red basin - warm
        lambda s: [0.1 + 0.8*s, 0.02 + 0.15*s, 0.05 + 0.1*s],
        # Blue basin - cool
        lambda s: [0.02 + 0.1*s, 0.05 + 0.2*s, 0.15 + 0.8*s],
        # Green basin - emerald
        lambda s: [0.02 + 0.12*s, 0.15 + 0.7*s, 0.08 + 0.25*s],
    ]
    
    for idx in range(3):
        mask = convergence == idx + 1
        shade = 1.0 - (speed[mask] / 80.0) * 0.85
        colors = palettes[idx](shade)
        for c in range(3):
            img[:,:,c][mask] = colors[c]
    
    ax.imshow(img, extent=[cx-span, cx+span, cy-span, cy+span], origin='lower')
    ax.set_xlabel(r'$\Re(z)$', fontsize=14)
    ax.set_ylabel(r'$\Im(z)$', fontsize=14)
    ax.set_title(r"Newton's Method: Zoomed Julia Set Boundary" + "\n" +
                 r"$z \mapsto z - \frac{z^3 - 1}{3z^2}$ — Infinite fractal complexity at every scale",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_julia_zoom.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_julia_zoom.pdf")


# ============================================================
# FIGURE B: Classic Julia Set (quadratic z^2 + c)
# ============================================================
def julia_classic(resolution=1000):
    print("Generating classic Julia set z^2 + c...")
    c = -0.7 + 0.27015j  # Famous beautiful Julia parameter
    
    x = np.linspace(-1.5, 1.5, resolution)
    y = np.linspace(-1.5, 1.5, resolution)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    
    escape = np.zeros(Z.shape, dtype=float)
    max_iter = 300
    
    for i in range(max_iter):
        mask = escape == 0
        if not np.any(mask):
            break
        Z_new = np.copy(Z)
        Z_new[mask] = Z[mask]**2 + c
        Z = Z_new
        
        escaped = (np.abs(Z) > 2) & (escape == 0)
        escape[escaped] = i + 1 - np.log2(np.log2(np.abs(Z[escaped]) + 1))
    
    escape[escape == 0] = max_iter
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    cmap = LinearSegmentedColormap.from_list('julia_fire', [
        (0.0, '#000000'),
        (0.05, '#1a0533'),
        (0.1, '#3d0066'),
        (0.15, '#6600aa'),
        (0.2, '#9900cc'),
        (0.3, '#cc3399'),
        (0.4, '#ff3366'),
        (0.5, '#ff6633'),
        (0.6, '#ffaa00'),
        (0.7, '#ffdd33'),
        (0.8, '#ffff66'),
        (0.9, '#ffffcc'),
        (1.0, '#ffffff'),
    ])
    
    im = ax.imshow(escape, extent=[-1.5, 1.5, -1.5, 1.5], origin='lower', 
                   cmap=cmap, vmin=0, vmax=100)
    ax.set_xlabel(r'$\Re(z)$', fontsize=14)
    ax.set_ylabel(r'$\Im(z)$', fontsize=14)
    ax.set_title(r"Classical Julia Set: $z \mapsto z^2 + c$,  $c = -0.7 + 0.27i$" + "\n" +
                 r"The quintessential fractal chaos that Pandrosion eliminates",
                 fontsize=13)
    plt.colorbar(im, ax=ax, shrink=0.75, label='Escape iterations')
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_julia_classic.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_julia_classic.pdf")


# ============================================================
# FIGURE C: Side-by-side Newton vs Pandrosion (HIGH RES)
# ============================================================  
def side_by_side_comparison(resolution=800):
    print("Generating side-by-side Newton vs Pandrosion comparison...")
    
    x = np.linspace(-2, 2, resolution)
    y = np.linspace(-2, 2, resolution)
    X_grid, Y_grid = np.meshgrid(x, y)
    
    target = 1.0
    roots = [1.0, np.exp(2j*np.pi/3), np.exp(-2j*np.pi/3)]
    
    # --- Newton ---
    Z_n = X_grid + 1j * Y_grid
    conv_n = np.zeros(Z_n.shape, dtype=int)
    speed_n = np.zeros(Z_n.shape, dtype=float)
    
    for iteration in range(60):
        mask = speed_n == 0
        if not np.any(mask):
            break
        Zm = Z_n[mask]
        denom = 3 * Zm**2
        safe = np.abs(denom) > 1e-12
        temp = np.copy(Zm)
        temp[safe] = Zm[safe] - (Zm[safe]**3 - target) / denom[safe]
        Z_new = np.copy(Z_n)
        Z_new[mask] = temp
        Z_n = Z_new
        for idx, root in enumerate(roots):
            close = np.abs(Z_n - root) < 1e-6
            newly = close & (speed_n == 0)
            conv_n[newly] = idx + 1
            speed_n[newly] = iteration + 1
    speed_n[speed_n == 0] = 60
    
    # --- Pandrosion ---
    Z_p = X_grid + 1j * Y_grid
    conv_p = np.zeros(Z_p.shape, dtype=int)
    speed_p = np.zeros(Z_p.shape, dtype=float)
    
    for iteration in range(80):
        mask = speed_p == 0
        if not np.any(mask):
            break
        Zm = Z_p[mask]
        denom = 3 * Zm**3 + 2 * target
        safe = np.abs(denom) > 1e-12
        temp = np.copy(Zm)
        temp[safe] = Zm[safe] * (Zm[safe]**3 + 4 * target) / denom[safe]
        Z_new = np.copy(Z_p)
        Z_new[mask] = temp
        Z_p = Z_new
        for idx, root in enumerate(roots):
            close = np.abs(Z_p - root) < 1e-6
            newly = close & (speed_p == 0)
            conv_p[newly] = idx + 1
            speed_p[newly] = iteration + 1
        close_zero = np.abs(Z_p) < 1e-6
        newly_zero = close_zero & (speed_p == 0)
        conv_p[newly_zero] = 4
        speed_p[newly_zero] = iteration + 1
    speed_p[speed_p == 0] = 80
    
    # Build images
    def build_img(conv, speed, max_speed):
        img = np.zeros((*conv.shape, 3))
        palettes = [
            lambda s: [0.1 + 0.85*s, 0.02 + 0.12*s, 0.03 + 0.08*s],
            lambda s: [0.02 + 0.08*s, 0.04 + 0.15*s, 0.12 + 0.85*s],
            lambda s: [0.02 + 0.1*s, 0.12 + 0.75*s, 0.06 + 0.2*s],
            lambda s: [0.05*s, 0.05*s, 0.05*s],
        ]
        for idx in range(4):
            mask = conv == idx + 1
            if not np.any(mask):
                continue
            shade = 1.0 - (speed[mask] / max_speed) * 0.8
            colors = palettes[idx](shade)
            for c in range(3):
                img[:,:,c][mask] = colors[c]
        return img
    
    img_n = build_img(conv_n, speed_n, 60)
    img_p = build_img(conv_p, speed_p, 80)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7.5))
    
    ax1.imshow(img_n, extent=[-2, 2, -2, 2], origin='lower')
    ax1.set_title(r"Newton: $z - \frac{z^3-1}{3z^2}$" + "\nFractal Julia boundaries", fontsize=14)
    ax1.set_xlabel(r'$\Re(z)$', fontsize=13)
    ax1.set_ylabel(r'$\Im(z)$', fontsize=13)
    for root in roots:
        ax1.plot(root.real, root.imag, 'w*', markersize=10, markeredgecolor='k', markeredgewidth=0.5)
    
    ax2.imshow(img_p, extent=[-2, 2, -2, 2], origin='lower')
    ax2.set_title(r"Pandrosion: $z\frac{z^3+4}{3z^3+2}$" + "\nSmooth sterile boundaries", fontsize=14)
    ax2.set_xlabel(r'$\Re(z)$', fontsize=13)
    ax2.set_ylabel(r'$\Im(z)$', fontsize=13)
    for root in roots:
        ax2.plot(root.real, root.imag, 'w*', markersize=10, markeredgecolor='k', markeredgewidth=0.5)
    ax2.plot(0, 0, 'wo', markersize=6, markeredgecolor='k', markeredgewidth=0.5)
    
    plt.suptitle(r'Basin Boundaries: Newton (Chaotic Fractals) vs Pandrosion (Topological Sterility)',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_julia_sidebyside.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_julia_sidebyside.pdf")


# ============================================================
# RUN ALL
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Julia Fractal Generator — Universitas Pandrosion")
    print("=" * 60)
    julia_newton_zoom(resolution=800)
    julia_classic(resolution=800)
    side_by_side_comparison(resolution=700)
    print("=" * 60)
    print("All 3 Julia figures generated!")
    print("=" * 60)
