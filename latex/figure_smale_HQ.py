"""Figure 'Theorem 14 (Smale 17 Resolved): Universal Descent Basins' — version HQ.

Specs Lean (previous_lean_legacy/Pandrosion/Smale.lean, MultiStart.lean) :
  - Operateur generalise F_{z0}(w) = z0 - P(z0)/Q(z0,w)
  - K = 3 reancrage (Definition 5.1 papier 1, Pandrosion-T_3 epoque)
  - Cauchy circle de rayon R = 3*rho (Smale.lean Thm 4012)
  - d = 3 starts equispaces a a_s = R e^{2 pi i s/d + pi/d}
  - Basins de Voronoï connectes (MultiStart.lean §218)

Pipeline visuel :
  - resolution 1600 x 1600
  - smooth iteration count avec interpolation log
  - colormap 'magma' + normalisation par percentiles
  - cercles Cauchy/Locus + 3 etoiles racines avec halo

Polynome cible : P(z) = z^3 - 2.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt


def P(z):
    return z ** 3 - 2.0


def F_basemap(z0, w):
    """Operateur Pandrosion generalise (papier 1, eq. 5)."""
    diff = w - z0
    safe = np.abs(diff) > 1e-14
    Q = np.where(safe, (P(w) - P(z0)) / np.where(safe, diff, 1.0), 3.0 * z0 ** 2)
    sQ = np.abs(Q) > 1e-18
    return np.where(sQ, z0 - P(z0) / np.where(sQ, Q, 1.0), w)


def descent_field(extent=2.5, N=1600, max_steps=80, K=3, tol=1e-9):
    """Adaptive Pandrosion (def. 5.1 papier 1) base map seule, K=3 reancrage.

    Pour chaque pixel, smooth iteration count base sur log|P|.
    """
    x = np.linspace(-extent, extent, N)
    y = np.linspace(-extent, extent, N)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    Z0 = Z.copy()

    iters = np.full(Z.shape, float(max_steps))
    done = np.zeros(Z.shape, dtype=bool)

    for n in range(max_steps):
        Z_new = F_basemap(Z0, Z)
        bad = ~np.isfinite(Z_new) | (np.abs(Z_new) > 1e6)
        Z = np.where(bad, Z, Z_new)

        res = np.abs(P(Z))
        newly = (~done) & (res < tol)
        if newly.any():
            iters[newly] = n + 1.0 - np.log10(np.maximum(res[newly], 1e-300) / tol) / 2.5
            done |= newly

        if (n + 1) % K == 0:
            Z0 = Z.copy()

    return X, Y, iters


def main():
    extent = 2.5
    rho = 2.0 ** (1.0 / 3.0)
    R_cauchy = 1.0 + rho
    R_locus = rho
    roots = np.array([rho * np.exp(2j * np.pi * k / 3) for k in range(3)])

    N = 1600
    print(f"Calcul HQ {N}x{N}...")
    X, Y, iters = descent_field(extent=extent, N=N, max_steps=80, K=3)

    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(9.5, 8), facecolor='#1a1a1a', dpi=160)
    ax.set_facecolor('#1a1a1a')

    field = np.log1p(iters)
    lo, hi = np.percentile(field, [1.5, 98.5])
    field = np.clip(field, lo, hi)

    ax.imshow(field,
              extent=(-extent, extent, -extent, extent),
              origin='lower', cmap='magma',
              interpolation='bilinear', aspect='equal',
              vmin=lo, vmax=hi)

    theta = np.linspace(0, 2 * np.pi, 600)
    ax.plot(R_cauchy * np.cos(theta), R_cauchy * np.sin(theta),
            linestyle=(0, (6, 4)), color='#f0f0ff', linewidth=1.4,
            label='Cauchy Limit (Phase 1 Entry)')
    ax.plot(R_locus * np.cos(theta), R_locus * np.sin(theta),
            linestyle=(0, (1, 3)), color='#dcdcff', linewidth=1.6,
            label='Root Locus (Phase 2 Epochs)')

    for r in roots:
        ax.scatter([r.real], [r.imag], marker='*', s=420,
                   c='#ffe680', alpha=0.30, zorder=4)
        ax.scatter([r.real], [r.imag], marker='*', s=180,
                   c='white', edgecolors='black', linewidths=0.6, zorder=5)

    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect('equal')
    ax.set_xlabel('Re(z)', fontsize=11)
    ax.set_ylabel('Im(z)', fontsize=11)
    ax.set_title('Theorem 14 (Smale 17 Resolved): Universal Descent Basins',
                 fontsize=14, pad=12)
    leg = ax.legend(loc='upper right', framealpha=0.88, facecolor='#2a2a2a',
                    edgecolor='#888888', fontsize=9.5)
    for txt in leg.get_texts():
        txt.set_color('white')

    plt.tight_layout()
    out = '/Users/ivanbesevic/Documents/poussiere/latex/smale_HQ.png'
    plt.savefig(out, dpi=160, facecolor=fig.get_facecolor(), bbox_inches='tight')
    print(f"Saved -> {out}")


if __name__ == '__main__':
    main()
