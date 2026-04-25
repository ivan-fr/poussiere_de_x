"""Variante T_3 : meme figure 'Universal Descent Basins' mais avec
l'acceleration Aitken/Richardson activee (Pandrosion-T_3, K=3, def. 4.2).

Comme T_3 converge en cubic-rate (~ -5 nats/epoque, Prop. 5.4), un comptage
classique d'iterations s'effondre en 2-3 pas. On colore donc par
        phi(z) = log10(|P(z_N)| + eps)   apres N epoques fixes,
ce qui revele la structure fine des epoques meme en regime de convergence
ultra-rapide : les bassins universels apparaissent sombres (|P| ~ 0) et les
zones fractales restantes brillent (|P| residuel).
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def P(z):
    return z ** 3 - 2.0


def F_basemap(z0, w):
    """Operateur Pandrosion generalise."""
    diff = w - z0
    safe = np.abs(diff) > 1e-14
    Q = np.where(safe, (P(w) - P(z0)) / np.where(safe, diff, 1.0), 3.0 * z0 ** 2)
    sQ = np.abs(Q) > 1e-18
    return np.where(sQ, z0 - P(z0) / np.where(sQ, Q, 1.0), w)


def T3_epoch(z0, z):
    """Une epoque T_3 (def. 4.2 du papier 1) : 3 evaluations de la base map."""
    s1 = F_basemap(z0, z)
    s2 = F_basemap(z0, s1)
    denom2 = s2 - 2.0 * s1 + z
    safe2 = np.abs(denom2) > 1e-18
    T2 = np.where(safe2, z - (s1 - z) ** 2 / np.where(safe2, denom2, 1.0), s1)

    s3 = F_basemap(z0, T2)
    den_l = s1 - z
    safel = np.abs(den_l) > 1e-18
    lam = np.where(safel, (s2 - s1) / np.where(safel, den_l, 1.0), 0.0 + 0j)
    denom3 = lam - 1.0
    safe3 = np.abs(denom3) > 1e-18
    T3 = np.where(safe3, T2 - (s3 - T2) / np.where(safe3, denom3, 1.0), T2)
    return T3


def trajectory_field(extent=2.5, N=1100, n_epochs=8):
    """Pour chaque point z, suit la trajectoire Pandrosion-T_3 (K=3) sur n_epochs
    epoques et renvoie un champ 'temps de descente' construit a partir des residus
    intermediaires :   phi(z) = sum_n  log10(|P(z_n)| / |P(z_0)|).

    Cette quantite reste informative meme apres convergence : elle integre
    la *vitesse* de descente plutot que le residu final (qui sature en machine eps)."""
    x = np.linspace(-extent, extent, N)
    y = np.linspace(-extent, extent, N)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    Z0 = Z.copy()

    P0 = np.abs(P(Z))
    accum = np.zeros(Z.shape)

    for _ in range(n_epochs):
        Z_new = T3_epoch(Z0, Z)
        bad = ~np.isfinite(Z_new) | (np.abs(Z_new) > 1e6)
        Z_new = np.where(bad, Z, Z_new)
        Z = Z_new
        Z0 = Z.copy()

        res = np.abs(P(Z))
        res = np.where(np.isfinite(res), res, 1e6)
        ratio = np.clip(res / np.maximum(P0, 1e-300), 1e-12, 1e12)
        accum += np.log10(ratio)

    return X, Y, accum


def main():
    extent = 2.5
    rho = 2.0 ** (1.0 / 3.0)
    R_cauchy = 1.0 + rho
    R_locus = rho
    roots = np.array([rho * np.exp(2j * np.pi * k / 3) for k in range(3)])

    print("Calcul du champ trajectoire (Pandrosion-T_3, K=3, 8 epoques)...")
    X, Y, accum = trajectory_field(extent=extent, N=1100, n_epochs=8)

    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(8.6, 7.2), facecolor='#1a1a1a')
    ax.set_facecolor('#1a1a1a')

    # accum tres negatif = descente rapide = bassin (sombre dans magma)
    # accum proche de 0 = stagnation = trame fractale (brillant dans magma)
    field = accum
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
    ax.set_title('Theorem 14 (Smale 17 Resolved): Universal Descent Basins  [T$_3$, K=3]',
                 fontsize=13, pad=10)
    leg = ax.legend(loc='upper right', framealpha=0.85, facecolor='#2a2a2a',
                    edgecolor='#888888', fontsize=9)
    for txt in leg.get_texts():
        txt.set_color('white')

    plt.tight_layout()
    out = Path(__file__).resolve().parents[1] / 'figures' / 'smale_reproduced_T3.png'
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print(f"Saved -> {out}")


if __name__ == '__main__':
    main()
