"""Reproduit 'Theorem 14 (Smale 17 Resolved): Universal Descent Basins'.

Visuel : meme rendu fractal lisse que smale.jpeg.
Algorithme : adaptive Pandrosion (papier 1, K=3) base sur l'operateur generalise
        F_{z0}(z) = z0 - P(z0)/Q(z0, z),  Q(z0,z) = (P(z) - P(z0))/(z - z0).

Trois pas de la base map par epoque, reancrage z0 <- z apres K=3 pas
(Definition 5.1, Pandrosion-T_3 avec K=3 ; sans acceleration Aitken pour
preserver la structure fractale, cf. Newton = limite K=1).

Coloration : smooth-iteration count (escape time normalise) sur magma,
similaire au standard fractal Newton/Julia.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def P(z):
    return z ** 3 - 2.0


def F_basemap(z0, w):
    """Operateur Pandrosion generalise F_{z0}(w) = z0 - P(z0)/Q(z0, w)."""
    diff = w - z0
    safe = np.abs(diff) > 1e-14
    Q = np.where(safe, (P(w) - P(z0)) / np.where(safe, diff, 1.0), 3.0 * z0 ** 2)
    sQ = np.abs(Q) > 1e-18
    return np.where(sQ, z0 - P(z0) / np.where(sQ, Q, 1.0), w)


def pandrosion_K3(extent=2.5, N=1100, max_steps=60, K=3, tol=1e-8):
    """Adaptive Pandrosion K=3 (def. 5.1) sur une grille NxN.

    Renvoie un champ smooth-iteration : log lisse du temps de descente
    de |P(z_n)|, qui produit le meme genre de bordures fractales que
    le standard Newton-basin plot.
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
        Z_new = np.where(bad, Z, Z_new)
        Z = Z_new

        res = np.abs(P(Z))
        newly = (~done) & (res < tol)
        if newly.any():
            # smooth count : interpolation sur log|P|
            iters[newly] = n + 1.0 - np.log10(np.maximum(res[newly], 1e-300) / tol) / 2.0
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

    print("Calcul des bassins (Pandrosion adaptatif, K=3)...")
    X, Y, iters = pandrosion_K3(extent=extent, N=1100, max_steps=60, K=3)

    # ----- rendu : magma sur fond sombre, contours fractals fins -----
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(8.6, 7.2), facecolor='#1a1a1a')
    ax.set_facecolor('#1a1a1a')

    # log lisse + clip pour saturer les zones de divergence
    field = np.log1p(iters)
    field = np.clip(field, np.percentile(field, 2), np.percentile(field, 98))

    levels = np.linspace(field.min(), field.max(), 80)
    ax.contourf(X, Y, field, levels=levels, cmap='magma')

    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(R_cauchy * np.cos(theta), R_cauchy * np.sin(theta),
            linestyle=(0, (6, 4)), color='#e6e6fa', linewidth=1.2,
            label='Cauchy Limit (Phase 1 Entry)')
    ax.plot(R_locus * np.cos(theta), R_locus * np.sin(theta),
            linestyle=(0, (1, 3)), color='#cfcfff', linewidth=1.4,
            label='Root Locus (Phase 2 Epochs)')

    ax.scatter(roots.real, roots.imag, marker='*', s=170,
               c='white', edgecolors='black', linewidths=0.6, zorder=5)

    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect('equal')
    ax.set_xlabel('Re(z)')
    ax.set_ylabel('Im(z)')
    ax.set_title('Theorem 14 (Smale 17 Resolved): Universal Descent Basins',
                 fontsize=13, pad=10)
    leg = ax.legend(loc='upper right', framealpha=0.85, facecolor='#2a2a2a',
                    edgecolor='#888888', fontsize=9)
    for txt in leg.get_texts():
        txt.set_color('white')

    plt.tight_layout()
    out = Path(__file__).resolve().parents[1] / 'figures' / 'smale_reproduced.png'
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print(f"Saved -> {out}")


if __name__ == '__main__':
    main()
