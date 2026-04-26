"""Hero visualization for paper 101ABC.
Style mimics abc.png: dark navy background, glowing 3D basin, golden particles."""
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap


# Palette matching abc.png aesthetic
BG = '#0c1d2e'
TEAL_OUT = '#1e6b6b'
TEAL_RIM = '#3aa394'
PALE_PLATEAU = '#dfe0ad'
ORANGE_HOT = '#ffa84d'
PURPLE_DEEP = '#3d105a'
GOLD = '#ffd84d'

basin_cmap = LinearSegmentedColormap.from_list(
    'basin',
    [(0.00, PURPLE_DEEP),     # basin floor (lowest z) = deep purple
     (0.10, '#7a266b'),
     (0.25, '#b85e3f'),
     (0.40, ORANGE_HOT),      # basin walls
     (0.55, '#e9c66a'),
     (0.70, PALE_PLATEAU),    # outer plateau
     (0.85, TEAL_RIM),
     (1.00, TEAL_OUT)],       # outer edge of canvas (highest z) = dark teal
    N=512)


def make_hero():
    mpl.rcParams.update({'font.family': 'serif', 'font.size': 11,
                         'mathtext.fontset': 'cm'})
    fig = plt.figure(figsize=(10.5, 9.0), facecolor=BG)
    ax = fig.add_subplot(111, projection='3d', facecolor=BG)

    # ---- Landscape: deep central basin + outer plateau, sloping rim ----
    N = 220
    extent = 2.5
    x = np.linspace(-extent, extent, N)
    y = np.linspace(-extent, extent, N)
    X, Y = np.meshgrid(x, y)
    R2 = X**2 + Y**2
    R  = np.sqrt(R2)

    # Flat plateau at z=0.4 (the pre-basin region)
    plateau = 0.4 + 0.0 * R
    # Deep funnel-shaped basin at origin (steep, narrow cone)
    basin = -2.0 * np.exp(-2.3 * R2)
    # Soft falloff at edges of the canvas (table-like)
    edge = -0.40 * np.maximum(0, R - 2.0)**2
    Z = plateau + basin + edge

    # Three side basins (other roots, smaller)
    for cx, cy, depth in [(1.55, 0.4, 0.30),
                          (-1.35, -0.65, 0.25),
                          (0.85, -1.45, 0.22)]:
        Z -= depth * np.exp(-3.5 * ((X - cx)**2 + (Y - cy)**2))

    # Render surface — use Z range for vivid colour contrast
    Zmin, Zmax = -1.6, 0.5
    ax.plot_surface(X, Y, Z,
                    facecolors=basin_cmap(np.clip((Z - Zmin) / (Zmax - Zmin), 0, 1)),
                    rstride=2, cstride=2,
                    linewidth=0, antialiased=True,
                    alpha=0.95, shade=True)

    # ---- Strategy B' spiral starts ----
    K_max = 14
    R0, q = 2.05, 0.74
    phi = math.pi * (3 - math.sqrt(5))

    def height(sx, sy):
        sr2 = sx*sx + sy*sy
        sR = math.sqrt(sr2)
        h_plat  = 0.4
        h_basin = -2.0 * math.exp(-2.3 * sr2)
        h_edge  = -0.40 * max(0, sR - 2.0)**2
        return h_plat + h_basin + h_edge

    sx_list = []; sy_list = []; sz_list = []
    for k in range(K_max):
        rk = R0 * q**k
        ang = k * phi
        sx = rk * math.cos(ang)
        sy = rk * math.sin(ang)
        sz = height(sx, sy) + 0.03
        sx_list.append(sx); sy_list.append(sy); sz_list.append(sz)

    # Trail
    ax.plot(sx_list, sy_list, sz_list, color=GOLD, lw=2.2, alpha=0.85, zorder=8)

    # Particles
    for k, (sx, sy, sz) in enumerate(zip(sx_list, sy_list, sz_list)):
        size = 280 - 14 * k
        ax.scatter([sx], [sy], [sz], s=size*4, c=GOLD, alpha=0.10, zorder=9)
        ax.scatter([sx], [sy], [sz], s=size, c=GOLD, alpha=0.95,
                   edgecolor='white', linewidth=0.8, zorder=10)

    # Bright convergence point at basin floor
    ax.scatter([0], [0], [-1.55], s=2500, c='white', alpha=0.13, zorder=11)
    ax.scatter([0], [0], [-1.55], s=600, c='white', alpha=0.97,
               edgecolor=GOLD, linewidth=1.6, zorder=12)

    # ---- Decorations ----
    # No axes (cinematic)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.grid(False)
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor((0, 0, 0, 0))

    # Camera — angled view into the basin, like abc.png
    ax.view_init(elev=32, azim=-35)
    ax.set_box_aspect((1, 1, 0.65))

    # ---- Title and side text (abc.png style) ----
    fig.text(0.5, 0.96,
             "Pandrosion synthesis: Strategy B$'$-$T_2$-ELS on Kostlan--Smale",
             ha='center', color='white', fontsize=15, weight='bold')
    fig.text(0.5, 0.925,
             r'Spiral starts descend through the pre-basin plateau into the Smale-$\alpha$ basin at the root',
             ha='center', color='#bcc6c6', fontsize=10.5, style='italic')

    # Legend (top-right, matching abc.png style)
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=GOLD,
                   markersize=11, markeredgecolor='white', markeredgewidth=0.6,
                   label=r'B$\prime$-spiral starts ($z_k = 2 q^k e^{ik\varphi}$)'),
        plt.Line2D([0], [0], marker='o', color='none', markerfacecolor='white',
                   markersize=12, markeredgecolor=GOLD, markeredgewidth=1.5,
                   label=r'Recovered root ($\|F\| < \tau_{\mathrm{tol}}$)'),
        plt.Rectangle((0, 0), 1, 1, fc=ORANGE_HOT, alpha=0.85, ec='white', lw=0.4,
                      label=r'Smale-$\alpha$ basin ($\alpha < \alpha^\star$)'),
        plt.Rectangle((0, 0), 1, 1, fc=PALE_PLATEAU, alpha=0.85, ec='white', lw=0.4,
                      label=r'Pre-basin plateau ($\alpha \geq \alpha^\star$)'),
    ]
    leg = ax.legend(handles=legend_elements, loc='upper right',
                    bbox_to_anchor=(1.18, 0.95),
                    facecolor='#091621', edgecolor='#445a6a', labelcolor='white',
                    framealpha=0.85, fontsize=9.5)
    for text in leg.get_texts():
        text.set_color('white')

    # Bottom-left: the closure result
    fig.text(0.06, 0.075,
             r'$\mathbb{E}_{P\sim\mathrm{KS}}[\mathrm{cost}_{\text{total}}] \;=\; \mathbf{O(d \log d)}$',
             color='white', fontsize=16, weight='bold')
    fig.text(0.06, 0.04,
             r'within $\log d$ of the lower bound $\Omega(d)$  ---  paper 101ABC closure theorem',
             color='#9cc4a8', fontsize=9.5, style='italic')

    # Bottom-right: scale of empirical evidence
    fig.text(0.94, 0.075,
             r'mvKS verified to $D = 4.3 \cdot 10^9$',
             color='white', fontsize=13, weight='bold', ha='right')
    fig.text(0.94, 0.04,
             r'$100\%$ convergence, fit $N_{\mathrm{tot}} = 3.67 + 2.31 \log_2\log_2 D$',
             color='#9cc4a8', fontsize=9.5, style='italic', ha='right')

    out = '/Users/ivanbesevic/Documents/poussiere/latex_legacy/figures/101abc_hero.png'
    fig.savefig(out, dpi=220, bbox_inches='tight', facecolor=BG, edgecolor='none')
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    make_hero()
