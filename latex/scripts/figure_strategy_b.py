"""Publication figure for Part VIII, Strategy B.

Two-panel comparison illustrating Definition 2.1 of Part VIII:

    z_k = rho * q^k * exp(i k phi),  phi = pi(3 - sqrt(5)).

The background of each panel is a fixed-epoch Pandrosion-T2 descent
landscape for P(z) = z^16 - 1.  This polynomial is chosen for visual
clarity: its Cauchy radius is exactly 2 and its roots lie on the unit
circle, so the Cauchy ring, root locus, and start sequences can be
shown on the same scale without artificial rescaling.

Left panel:  Strategy A (equispaced Cauchy ring) -- the Part VII baseline.
Right panel: Strategy B (geometric-decreasing golden-angle spiral).

Each start is shown together with the actual Pandrosion-T2 orbit it
generates, so the reader sees not only the geometry of the start
sequence but the convergence behaviour it induces on the landscape.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-poussiere")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/poussiere-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize


DEGREE = 16
Q = 0.7
RHO = 2.0
EXTENT = 2.42
N_GRID = 960
N_EPOCHS = 7
ORBIT_EPOCHS = 4
EPS = 1e-300

OUT = Path(__file__).resolve().parents[1] / "figures"

BG = "#0c0a14"
PANEL = "#13101d"
INK = "#f6efff"
SUBINK = "#bcb1d4"
ACCENT_A = "#7fdcff"
ACCENT_B = "#ffd166"


def poly(z):
    return z**DEGREE - 1.0


def dpoly(z):
    return DEGREE * z ** (DEGREE - 1)


def pandrosion_base(z0, w):
    """Generalized Pandrosion base map anchored at z0."""
    diff = w - z0
    safe = np.abs(diff) > 1e-13
    q_val = np.where(safe, (poly(w) - poly(z0)) / np.where(safe, diff, 1.0), dpoly(z0))
    ok = np.abs(q_val) > 1e-16
    return np.where(ok, z0 - poly(z0) / np.where(ok, q_val, 1.0), w)


def t2_epoch(z0, z):
    """One Pandrosion-T2 epoch."""
    s1 = pandrosion_base(z0, z)
    s2 = pandrosion_base(z0, s1)
    denom = s2 - 2.0 * s1 + z
    ok = np.abs(denom) > 1e-15
    out = np.where(ok, z - (s1 - z) ** 2 / np.where(ok, denom, 1.0), s1)
    bad = ~np.isfinite(out) | (np.abs(out) > 1e8)
    return np.where(bad, z, out)


def descent_landscape():
    x = np.linspace(-EXTENT, EXTENT, N_GRID)
    y = np.linspace(-EXTENT, EXTENT, N_GRID)
    x_grid, y_grid = np.meshgrid(x, y)
    z = x_grid + 1j * y_grid
    anchor = z.copy()
    initial = np.maximum(np.abs(poly(z)), EPS)

    accum = np.zeros_like(x_grid, dtype=float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for _ in range(N_EPOCHS):
            z = t2_epoch(anchor, z)
            anchor = z.copy()
            residual = np.maximum(np.abs(poly(z)), EPS)
            ratio = np.clip(residual / initial, 1e-14, 1e14)
            accum += np.log10(ratio)

    field = np.clip(accum, np.percentile(accum, 2.0), np.percentile(accum, 98.6))
    return x_grid, y_grid, field


def trace_orbit(z0_complex, n_epochs=ORBIT_EPOCHS):
    """Run a scalar Pandrosion-T2 orbit from z0 and return the trajectory."""
    z = np.array([complex(z0_complex)], dtype=complex)
    anchor = z.copy()
    traj = [complex(z[0])]
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for _ in range(n_epochs):
            z = t2_epoch(anchor, z)
            anchor = z.copy()
            zc = complex(z[0])
            if not np.isfinite(zc.real) or not np.isfinite(zc.imag):
                break
            traj.append(zc)
            if abs(poly(zc)) < 1e-10:
                break
    return np.array(traj, dtype=complex)


def strategy_b_starts(n=DEGREE):
    golden = math.pi * (3.0 - math.sqrt(5.0))
    k = np.arange(n)
    return RHO * (Q**k) * np.exp(1j * k * golden)


def strategy_a_starts(n=DEGREE):
    k = np.arange(n)
    return RHO * np.exp(2j * np.pi * k / n)


def landscape_cmap():
    """Custom violet-magenta-amber colormap tuned for this background."""
    stops = [
        (0.00, "#0a0712"),
        (0.18, "#1a0f2e"),
        (0.38, "#3f1652"),
        (0.58, "#8a2360"),
        (0.78, "#d94d4a"),
        (0.92, "#f6a141"),
        (1.00, "#fde8a4"),
    ]
    positions = [s[0] for s in stops]
    colors = [s[1] for s in stops]
    return LinearSegmentedColormap.from_list("pandrosion_dusk", list(zip(positions, colors)), N=512)


def draw_panel(ax, x_grid, y_grid, field, starts, *, label, accent, cmap, show_spiral_curve):
    ax.set_facecolor(PANEL)
    levels = np.linspace(field.min(), field.max(), 110)
    im = ax.contourf(x_grid, y_grid, field, levels=levels, cmap=cmap, antialiased=True)

    theta = np.linspace(0, 2 * np.pi, 720)

    # Faint geometric radius ladder shared by both panels (shows q-decay scale).
    for j, radius in enumerate(RHO * Q ** np.arange(0, 9)):
        if radius < 0.08:
            continue
        alpha = 0.32 if j == 0 else 0.09
        ax.plot(
            radius * np.cos(theta),
            radius * np.sin(theta),
            color=INK,
            linewidth=0.7,
            alpha=alpha,
            linestyle=(0, (2, 6)),
            zorder=3,
        )

    # Cauchy radius (heavier dashed white).
    ax.plot(
        RHO * np.cos(theta),
        RHO * np.sin(theta),
        color="#efe6ff",
        linewidth=1.55,
        linestyle=(0, (8, 5)),
        alpha=0.92,
        zorder=4,
        label=r"Cauchy ring $\rho$",
    )
    # Root locus (unit circle for z^16 - 1).
    ax.plot(
        np.cos(theta),
        np.sin(theta),
        color="#d6cdff",
        linewidth=1.35,
        linestyle=(0, (1, 3)),
        alpha=0.85,
        zorder=4,
        label="Root locus",
    )

    # Continuous spiral curve underlying Strategy B (only on the right panel).
    if show_spiral_curve:
        golden = math.pi * (3.0 - math.sqrt(5.0))
        kk = np.linspace(0, DEGREE - 1, 800)
        spiral = RHO * (Q**kk) * np.exp(1j * kk * golden)
        # Soft outer halo + crisp inner stroke so the spiral reads as the figure's spine.
        ax.plot(
            spiral.real, spiral.imag,
            color=accent, linewidth=4.5, alpha=0.18,
            zorder=5, solid_capstyle="round",
        )
        ax.plot(
            spiral.real, spiral.imag,
            color=accent, linewidth=1.6, alpha=0.92,
            zorder=6, solid_capstyle="round",
        )

    # Orbit trajectories: faint trails fading along their length, segment by segment.
    norm = Normalize(vmin=0, vmax=max(1, len(starts) - 1))
    trail_cmap = plt.colormaps["plasma"] if show_spiral_curve else plt.colormaps["cool"]

    for idx, z0 in enumerate(starts):
        orbit = trace_orbit(z0)
        if orbit.size < 2:
            continue
        # Clip orbit so wild excursions outside the view do not pollute the figure.
        in_view = (np.abs(orbit.real) <= EXTENT * 1.02) & (np.abs(orbit.imag) <= EXTENT * 1.02)
        if not in_view[0]:
            continue
        last = np.argmax(~in_view) if (~in_view).any() else len(orbit)
        if last == 0:
            last = 1
        orbit = orbit[: max(2, last + 1)] if last < len(orbit) else orbit
        if orbit.size < 2:
            continue
        col = np.array(trail_cmap(0.15 + 0.75 * norm(idx)))
        pts = np.column_stack([orbit.real, orbit.imag])
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        # Fade alpha along the trail: brightest near the start, dimmer near the root.
        n_seg = len(segs)
        alphas = np.linspace(0.85, 0.18, n_seg)
        seg_colors = np.tile(col, (n_seg, 1))
        seg_colors[:, 3] = alphas
        lc = LineCollection(
            segs,
            colors=seg_colors,
            linewidths=np.linspace(1.25, 0.45, n_seg),
            capstyle="round",
            zorder=7,
        )
        ax.add_collection(lc)

    # Start markers: filled disks coloured by index.
    sc = ax.scatter(
        starts.real,
        starts.imag,
        c=np.arange(len(starts)),
        cmap=trail_cmap,
        vmin=0,
        vmax=max(1, len(starts) - 1),
        s=88,
        edgecolors="#0c0813",
        linewidths=0.9,
        zorder=10,
    )

    # Index labels for the first few starts (so the reader sees the order).
    label_indices = (0, 1, 2, 3, 5, len(starts) - 1)
    for idx in label_indices:
        if 0 <= idx < len(starts):
            zk = starts[idx]
            r = abs(zk)
            if r < 1e-9:
                ox, oy = 0.12, 0.12
            else:
                # Push the label radially outward so it never sits on top of the marker.
                ox = (zk.real / r) * 0.16
                oy = (zk.imag / r) * 0.16
            ax.text(
                zk.real + ox,
                zk.imag + oy,
                f"$k={idx}$",
                color=INK,
                fontsize=7.8,
                alpha=0.95,
                ha="center",
                va="center",
                zorder=11,
                path_effects=[pe.withStroke(linewidth=1.8, foreground="#0a0712", alpha=0.85)],
            )

    # Roots of z^16 - 1 as small stars on the unit circle.
    roots = np.exp(2j * np.pi * np.arange(DEGREE) / DEGREE)
    ax.scatter(
        roots.real,
        roots.imag,
        marker="*",
        s=78,
        c="#fff7d6",
        edgecolors="#0c0813",
        linewidths=0.55,
        zorder=12,
        label=r"Roots of $z^{16}-1$",
    )

    ax.set_xlim(-EXTENT, EXTENT)
    ax.set_ylim(-EXTENT, EXTENT)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\mathrm{Re}(z)$", fontsize=10.5, color=SUBINK)
    ax.set_ylabel(r"$\mathrm{Im}(z)$", fontsize=10.5, color=SUBINK)
    ax.tick_params(colors=SUBINK, labelsize=8.5, length=3, width=0.6)
    for spine in ax.spines.values():
        spine.set_color("#3a324a")
        spine.set_linewidth(0.8)

    # Panel label (top-left badge).
    ax.text(
        0.025,
        0.965,
        label,
        transform=ax.transAxes,
        color=accent,
        fontsize=12,
        weight="bold",
        ha="left",
        va="top",
        bbox=dict(
            facecolor="#1a1428",
            edgecolor=accent,
            linewidth=0.9,
            boxstyle="round,pad=0.32",
            alpha=0.92,
        ),
        zorder=15,
    )

    return im, sc


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("Computing Pandrosion-T2 descent landscape for P(z) = z^16 - 1 ...")
    x_grid, y_grid, field = descent_landscape()

    starts_a = strategy_a_starts()
    starts_b = strategy_b_starts()

    cmap = landscape_cmap()

    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
    })

    fig = plt.figure(figsize=(14.6, 8.2), facecolor=BG)
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[1.0, 0.045],
        width_ratios=[1.0, 1.0],
        left=0.055, right=0.985, top=0.905, bottom=0.085,
        wspace=0.10, hspace=0.16,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[1, :])

    im_a, _ = draw_panel(
        ax_a, x_grid, y_grid, field, starts_a,
        label="Strategy A  ·  equispaced Cauchy ring",
        accent=ACCENT_A,
        cmap=cmap,
        show_spiral_curve=False,
    )
    im_b, _ = draw_panel(
        ax_b, x_grid, y_grid, field, starts_b,
        label=r"Strategy B  ·  $\rho q^{k} e^{i k \varphi}$,  $q=0.7$,  $\varphi=\pi(3-\sqrt{5})$",
        accent=ACCENT_B,
        cmap=cmap,
        show_spiral_curve=True,
    )

    # Shared horizontal colorbar at the bottom.
    cbar = fig.colorbar(im_b, cax=cax, orientation="horizontal")
    cbar.set_label(
        r"cumulative $\log_{10}$ residual ratio after %d Pandrosion-$T_2$ epochs   "
        r"(darker = faster contraction toward a root)" % N_EPOCHS,
        color=SUBINK, fontsize=9.5,
    )
    cbar.ax.tick_params(colors=SUBINK, labelsize=8, length=2.5, width=0.5)
    cbar.outline.set_edgecolor("#3a324a")
    cbar.outline.set_linewidth(0.6)

    # Shared legend, harvested from the right panel.
    handles, labels = ax_b.get_legend_handles_labels()
    leg = fig.legend(
        handles, labels,
        loc="lower left",
        bbox_to_anchor=(0.055, 0.005),
        ncol=3,
        frameon=False,
        fontsize=9,
        labelcolor=INK,
    )
    for text in leg.get_texts():
        text.set_color(INK)

    fig.suptitle(
        r"Strategy B vs. Strategy A on a Pandrosion-$T_2$ descent landscape  ($P(z)=z^{16}-1$)",
        color=INK, fontsize=15.5, y=0.965,
    )

    fig.text(
        0.985, 0.012,
        r"Universitas Pandrosion · Part VIII",
        color=SUBINK, fontsize=8.2, ha="right", va="bottom", alpha=0.85,
    )

    png = OUT / "strategy_b_geometric_spiral.png"
    fig.savefig(png, dpi=240, facecolor=fig.get_facecolor())
    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
