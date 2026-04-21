"""
Universitas Pandrosion — paper figure generator.

Produces a family of four paper-quality figures visualizing the
algorithmic geometry and convergence of the Pandrosion-Steffensen
algorithm for z^p = x, matching the 21-module Lean 4 formalization.

Figures emitted (into latex/latex/):
  fig_cyclotomic_basins.png     — cyclotomic anchors + Voronoï basins in C
  fig_loglog_convergence.png    — quadratic loglog tail from inside a basin
  fig_module_graph.png          — 21-module dependency DAG
  fig_efficiency_index.png      — Kung-Traub efficiency frontier, c=2 optimum
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.collections import LineCollection
import numpy as np

OUT = Path(__file__).parent / "latex"
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "mathtext.fontset": "cm",
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})


# ------------------------------------------------------------------
# Figure 1: cyclotomic anchors with Voronoï basins
# ------------------------------------------------------------------
def fig_cyclotomic_basins():
    p = 7
    x_val = 4.0
    alpha_abs = x_val ** (-1.0 / p)
    anchors = np.array([alpha_abs * np.exp(2j * np.pi * k / p) for k in range(p)])

    # Background: nearest-anchor Voronoï field.
    res = 800
    extent = 1.4
    xs = np.linspace(-extent, extent, res)
    ys = np.linspace(-extent, extent, res)
    X, Y = np.meshgrid(xs, ys)
    Z = X + 1j * Y
    dists = np.array([np.abs(Z - a) for a in anchors])
    nearest = np.argmin(dists, axis=0)

    fig, ax = plt.subplots(figsize=(7.5, 7.2))
    cmap = plt.get_cmap("Pastel1", p)
    ax.imshow(
        nearest,
        extent=(-extent, extent, -extent, extent),
        origin="lower",
        cmap=cmap,
        alpha=0.75,
        interpolation="nearest",
    )

    # Voronoï bisectors (explicit pairs).
    for i in range(p):
        for j in range(i + 1, p):
            a, b = anchors[i], anchors[j]
            mid = (a + b) / 2
            direction = b - a
            normal = 1j * direction / abs(direction)
            t = np.linspace(-3, 3, 2)
            seg = np.array([mid + tt * normal for tt in t])
            ax.plot(seg.real, seg.imag, color="white", lw=1.4, alpha=0.9)
            ax.plot(seg.real, seg.imag, color="#333", lw=0.6, alpha=0.75)

    # Anchors γ_k.
    ax.scatter(
        anchors.real,
        anchors.imag,
        s=140,
        c="#111",
        marker="*",
        zorder=5,
        edgecolors="white",
        linewidths=1.0,
    )
    for k, a in enumerate(anchors):
        offset = 0.12 * np.exp(1j * (2 * np.pi * k / p))
        ax.annotate(
            rf"$\gamma_{{{k}}}$",
            xy=(a.real, a.imag),
            xytext=(a.real + offset.real, a.imag + offset.imag),
            fontsize=12,
            ha="center",
            zorder=6,
        )

    # Circle |z| = α = x^{-1/p}.
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(
        alpha_abs * np.cos(theta),
        alpha_abs * np.sin(theta),
        color="#222",
        lw=1.0,
        ls="--",
        alpha=0.55,
        label=rf"$|z| = x^{{-1/p}} \approx {alpha_abs:.3f}$",
    )

    # -- Actual Pandrosion-Steffensen orbit, computed from the real map.
    def Sp(z, p):
        # Sum_{k=0}^{p-1} z^k.
        return sum(z ** k for k in range(p))

    def h(z, x, p):
        return 1 - (x - 1) / (x * Sp(z, p))

    def steffensen(z, x, p):
        hz = h(z, x, p)
        hhz = h(hz, x, p)
        D = hhz - 2 * hz + z
        if abs(D) < 1e-14:
            return z
        return z - (hz - z) ** 2 / D

    # Pick a starting point inside the local Steffensen basin of γ_1.
    # For p=7, x=4, the local super-attracting disk around each γ_k is
    # tight (order ~ 0.04); the orbit below shows 5 real iterates with
    # explicit quadratic collapse:
    #   ‖z_k − γ_1‖ ≈ 3.6e-2, 2.5e-2, 1.3e-2, 2.9e-3, 1.6e-4, 4.7e-7, …
    target = anchors[1]
    z = target + 0.030 + 0.020j
    traj = [z]
    for _ in range(6):
        z = steffensen(z, x_val, p)
        traj.append(z)
        if abs(traj[-1] - traj[-2]) < 1e-12:
            break
    traj = np.array(traj)
    ax.plot(
        traj.real,
        traj.imag,
        color="#b00000",
        lw=1.4,
        marker="o",
        markersize=4.5,
        alpha=0.9,
        zorder=4,
    )
    # Mark the starting point.
    ax.scatter(
        [traj[0].real], [traj[0].imag],
        s=80, marker="o",
        facecolor="white", edgecolor="#b00000", lw=1.6, zorder=5,
    )
    ax.annotate(
        r"$z_0$",
        xy=(traj[0].real, traj[0].imag),
        xytext=(traj[0].real + 0.18, traj[0].imag + 0.20),
        fontsize=12,
        color="#b00000",
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#b00000", lw=1.0),
    )
    # Zoom inset showing the quadratic collapse into γ_1.
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
    axins = inset_axes(
        ax, width="36%", height="36%", loc="lower left",
        bbox_to_anchor=(0.015, 0.015, 1, 1),
        bbox_transform=ax.transAxes, borderpad=1,
    )
    axins.set_facecolor("white")
    axins.plot(
        traj.real, traj.imag,
        color="#b00000", lw=1.5,
        marker="o", markersize=6, alpha=0.95,
        zorder=4,
    )
    axins.scatter(
        [target.real], [target.imag],
        marker="*", s=220,
        color="#111", edgecolor="white", linewidths=1.0, zorder=5,
    )
    axins.annotate(
        r"$\gamma_1$",
        xy=(target.real, target.imag),
        xytext=(target.real - 0.013, target.imag - 0.012),
        fontsize=10,
    )
    # Per-iterate labels z_0,…,z_4 and error.
    for i, zi in enumerate(traj[:5]):
        err = abs(zi - target)
        axins.annotate(
            rf"$z_{{{i}}}$",
            xy=(zi.real, zi.imag),
            xytext=(zi.real + 0.004, zi.imag + 0.005),
            fontsize=8, color="#b00000",
        )
    axins.set_xlim(target.real - 0.025, target.real + 0.050)
    axins.set_ylim(target.imag - 0.020, target.imag + 0.045)
    axins.set_aspect("equal")
    axins.set_xticks([])
    axins.set_yticks([])
    axins.set_title(
        r"zoom: real $\sigma$-orbit, $\|z_k-\gamma_1\| \to 0$ quadratically",
        fontsize=8,
    )
    axins.grid(alpha=0.3)
    for spine in axins.spines.values():
        spine.set_edgecolor("#b00000")
        spine.set_linewidth(1.3)

    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\mathrm{Re}(z)$")
    ax.set_ylabel(r"$\mathrm{Im}(z)$")
    ax.set_title(
        rf"Cyclotomic anchors $\gamma_k = \alpha\,e^{{2\pi i k/p}}$ "
        rf"and Voronoï basins in $\mathbb{{C}}$ "
        rf"($p={p}$, $x={x_val:g}$)"
    )
    ax.legend(loc="lower right", framealpha=0.92)
    ax.grid(alpha=0.15)

    plt.savefig(OUT / "fig_cyclotomic_basins.png")
    plt.close(fig)


# ------------------------------------------------------------------
# Figure 2: loglog convergence tail
# ------------------------------------------------------------------
def fig_loglog_convergence():
    fig, ax = plt.subplots(figsize=(7.5, 4.8))

    K = 1.5
    r = 1.0 / (2 * K)

    # Quadratic recurrence: e_{k+1} = K · e_k^2, starting at e_0 = 0.4·r.
    e_quad = [0.4 * r]
    for _ in range(12):
        e_quad.append(K * e_quad[-1] ** 2)
    e_quad = np.array(e_quad)

    # Linear recurrence for comparison: e_{k+1} = lambda · e_k.
    lam = 0.5
    e_lin = [0.4 * r]
    for _ in range(30):
        e_lin.append(lam * e_lin[-1])
    e_lin = np.array(e_lin)

    ks_quad = np.arange(len(e_quad))
    ks_lin = np.arange(len(e_lin))

    ax.semilogy(
        ks_lin,
        e_lin,
        color="#1f77b4",
        lw=2.0,
        marker="s",
        markersize=4.5,
        label=r"Linear Pandrosion:  $e_{k+1} = \lambda\,e_k$",
    )
    ax.semilogy(
        ks_quad,
        e_quad,
        color="#b00000",
        lw=2.0,
        marker="o",
        markersize=5.5,
        label=r"Quadratic Steffensen:  $e_{k+1} \leq K\,e_k^{\,2}$",
    )

    # Horizontal epsilon targets.
    eps_targets = [1e-4, 1e-8, 1e-16]
    for eps in eps_targets:
        ax.axhline(eps, color="gray", lw=0.6, ls=":", alpha=0.6)
        ax.text(
            0.35,
            eps * 1.6,
            rf"$\varepsilon = 10^{{{int(math.log10(eps))}}}$",
            fontsize=8,
            color="gray",
        )

    ax.set_xlabel(r"iteration count $k$")
    ax.set_ylabel(r"error $\|z_k - \gamma_s\|$")
    ax.set_title(
        r"Quadratic loglog tail inside the basin "
        r"$\mathcal{B}(\gamma_s,\,r_s)$ "
        r"(§34.1 \texttt{quadratic\_loglog\_from\_basin})"
    )
    ax.set_ylim(1e-18, 1.0)
    ax.set_xlim(0, 30)
    ax.grid(which="both", alpha=0.25)
    ax.legend(loc="lower left", framealpha=0.95)

    # Annotation: N_tail = O(log log 1/eps).
    ax.annotate(
        r"$N_{\mathrm{tail}}(\varepsilon) = \mathcal{O}(\log\log(1/\varepsilon))$",
        xy=(6, 1e-12),
        xytext=(12, 1e-6),
        fontsize=11,
        color="#b00000",
        arrowprops=dict(arrowstyle="->", color="#b00000", lw=1.2),
    )

    plt.savefig(OUT / "fig_loglog_convergence.png")
    plt.close(fig)


# ------------------------------------------------------------------
# Figure 3: 21-module dependency graph
# ------------------------------------------------------------------
def fig_module_graph():
    fig, ax = plt.subplots(figsize=(11.5, 8.6))

    # (name, x, y, palette) — laid out by layer.
    nodes = {
        "Foundations":              (0.00,  5.0, "alg"),
        "MultiStartBasins":         (0.00,  4.2, "alg"),
        "VoronoiMeasure":           (3.40,  4.2, "meas"),
        "CyclotomicMcMullen":       (3.40,  3.4, "meas"),
        "QuadraticRate":            (-3.40, 3.4, "cplx"),
        "QuadraticComplexity":      (0.00,  3.4, "cplx"),
        "QuadraticBound":           (-5.50, 2.5, "bnd"),
        "LinearAsymptotics":        (-5.50, 1.5, "bnd"),
        "KungTraub":                (-3.40, 2.5, "opt"),
        "UniformComplexity":        (0.00,  2.5, "cplx"),
        "SuperGrandMaster":         (-3.40, 1.5, "cplx"),
        "UniformContractionRate":   (0.00,  1.5, "cplx"),
        "ConcreteIteration":        (3.40,  1.5, "cplx"),
        "ComplexMultiplier":        (-3.40, 0.5, "dyn"),
        "LocalAttraction":          (0.00,  0.5, "dyn"),
        "DynamicalConvergence":     (3.40,  0.5, "dyn"),
        "SteffensenQuadraticBound": (-3.40, -0.5, "new"),
        "SteffensenMcMullenAE":     (0.00, -0.5, "new"),
        "SteffensenExplicitRate":   (3.40, -0.5, "new"),
        "SteffensenGlobalLoglog":   (0.00, -1.5, "new"),
        "MasterAbsolu":             (3.40, -1.5, "root"),
    }

    palette = {
        "alg":  ("#cfe3ff", "#1e4c8a"),
        "meas": ("#e9d5ff", "#5e2fa4"),
        "cplx": ("#d1f2d1", "#2e7d32"),
        "opt":  ("#fff1b8", "#8a6d0a"),
        "bnd":  ("#ffddee", "#9b2c6f"),
        "dyn":  ("#ffd0c4", "#a6341b"),
        "new":  ("#f5c99b", "#8a4a10"),
        "root": ("#ffb870", "#c34e00"),
    }

    def draw_node(name, x, y, kind):
        face, edge = palette[kind]
        w, h = 2.80, 0.56
        box = FancyBboxPatch(
            (x - w / 2, y - h / 2),
            w, h,
            boxstyle="round,pad=0.02,rounding_size=0.12",
            linewidth=1.2,
            edgecolor=edge,
            facecolor=face,
            zorder=3,
        )
        ax.add_patch(box)
        ax.text(
            x, y, name,
            ha="center", va="center",
            fontsize=9, family="monospace",
            color=edge, zorder=4,
        )

    # Edges (source -> target).
    edges = [
        ("Foundations", "MultiStartBasins"),
        ("Foundations", "ComplexMultiplier"),
        ("Foundations", "QuadraticBound"),
        ("MultiStartBasins", "VoronoiMeasure"),
        ("MultiStartBasins", "QuadraticRate"),
        ("MultiStartBasins", "QuadraticComplexity"),
        ("VoronoiMeasure", "CyclotomicMcMullen"),
        ("QuadraticRate", "KungTraub"),
        ("QuadraticComplexity", "UniformComplexity"),
        ("QuadraticBound", "LinearAsymptotics"),
        ("QuadraticBound", "SteffensenQuadraticBound"),
        ("UniformComplexity", "SuperGrandMaster"),
        ("UniformComplexity", "UniformContractionRate"),
        ("UniformContractionRate", "ConcreteIteration"),
        ("ComplexMultiplier", "LocalAttraction"),
        ("LocalAttraction", "DynamicalConvergence"),
        ("CyclotomicMcMullen", "DynamicalConvergence"),
        ("DynamicalConvergence", "SteffensenQuadraticBound"),
        ("SteffensenQuadraticBound", "SteffensenMcMullenAE"),
        ("SteffensenQuadraticBound", "SteffensenExplicitRate"),
        ("SteffensenExplicitRate", "SteffensenGlobalLoglog"),
        ("SteffensenMcMullenAE", "SteffensenGlobalLoglog"),
        ("KungTraub", "MasterAbsolu"),
        ("SuperGrandMaster", "MasterAbsolu"),
        ("UniformContractionRate", "MasterAbsolu"),
        ("CyclotomicMcMullen", "MasterAbsolu"),
        ("SteffensenGlobalLoglog", "MasterAbsolu"),
    ]

    for src, dst in edges:
        x1, y1, _ = nodes[src]
        x2, y2, _ = nodes[dst]
        ax.annotate(
            "",
            xy=(x2, y2 + 0.28 if y2 > y1 else y2 - 0.28),
            xytext=(x1, y1 - 0.28 if y2 < y1 else y1 + 0.28),
            arrowprops=dict(
                arrowstyle="->",
                color="#555",
                lw=0.8,
                alpha=0.7,
                connectionstyle="arc3,rad=0.08",
            ),
            zorder=1,
        )

    for name, (x, y, kind) in nodes.items():
        draw_node(name, x, y, kind)

    # Legend.
    legend_handles = [
        mpatches.Patch(color=palette[k][0], label=lbl)
        for k, lbl in [
            ("alg",  "Algebraic foundations"),
            ("meas", "Measure / cyclotomic"),
            ("cplx", "Rates & complexity"),
            ("opt",  "Kung-Traub optimality"),
            ("bnd",  "Taylor bounds"),
            ("dyn",  "Complex dynamics"),
            ("new",  "Steffensen unconditional spine (§31-§34)"),
            ("root", "Master Absolu (stitching)"),
        ]
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        framealpha=0.95,
        fontsize=8.5,
        title="Module layer",
        title_fontsize=9,
    )

    ax.set_xlim(-7.5, 8.0)
    ax.set_ylim(-2.4, 5.8)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        "21-module dependency graph of the Pandrosion Lean~4 formalization",
        fontsize=13,
    )

    plt.savefig(OUT / "fig_module_graph.png")
    plt.close(fig)


# ------------------------------------------------------------------
# Figure 4: efficiency index E(q, c) and the Kung-Traub frontier
# ------------------------------------------------------------------
def fig_efficiency_index():
    fig, ax = plt.subplots(figsize=(7.5, 5.2))

    # Kung-Traub frontier: q* = 2^(c-1), E* = 2^((c-1)/c).
    c_grid = np.linspace(1, 8, 400)
    kt_order = 2 ** (c_grid - 1)
    kt_efficiency = 2 ** ((c_grid - 1) / c_grid)

    ax.plot(
        c_grid,
        kt_efficiency,
        color="#b00000",
        lw=2.2,
        label=r"Kung-Traub frontier  $\mathrm{KT}(c) = 2^{(c-1)/c}$",
    )

    # Feasible region (below frontier).
    ax.fill_between(
        c_grid,
        1.0,
        kt_efficiency,
        color="#ffd0c4",
        alpha=0.35,
        label="feasible efficiency region",
    )

    # Discrete Kung-Traub-compliant methods with q = 2^(c-1).
    c_int = np.arange(1, 9)
    kt_int = 2 ** ((c_int - 1) / c_int)
    ax.scatter(
        c_int, kt_int,
        s=60,
        facecolor="white",
        edgecolor="#b00000",
        lw=1.5,
        zorder=5,
        label=r"integer-$c$ optimal methods",
    )

    # --- Pandrosion-Steffensen ---
    pand_E = 2 ** 0.5
    # Horizontal reference line across the plot.
    ax.axhline(
        pand_E,
        color="#1f77b4",
        lw=1.6,
        ls="--",
        alpha=0.85,
        label=rf"Pandrosion-Steffensen  $E=\sqrt{{2}}\approx {pand_E:.3f}$",
    )
    # Vertical drop line from the marker to the x-axis.
    ax.plot([2, 2], [1.0, pand_E], color="#1f77b4", lw=1.1, ls=":", alpha=0.85)
    # Big fat star marker.
    ax.scatter(
        [2], [pand_E],
        marker="*", s=520,
        color="#1f77b4", edgecolor="white",
        linewidths=1.5, zorder=8,
    )
    ax.annotate(
        r"Pandrosion-Steffensen" "\n" r"$(c{=}2,\,q{=}2)$",
        xy=(2, pand_E),
        xytext=(3.5, 1.52),
        fontsize=10, color="#1f77b4", fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=1.4),
    )

    # --- Secant ---
    phi = (1 + 5 ** 0.5) / 2
    sec_E = phi ** 0.5
    ax.axhline(
        sec_E,
        color="#2e7d32",
        lw=1.6,
        ls="--",
        alpha=0.85,
        label=rf"Secant (no-memory model)  $E=\sqrt{{\varphi}}\approx {sec_E:.3f}$",
    )
    ax.plot([2, 2], [1.0, sec_E], color="#2e7d32", lw=1.1, ls=":", alpha=0.85)
    ax.scatter(
        [2], [sec_E],
        marker="D", s=160,
        color="#2e7d32", edgecolor="white",
        linewidths=1.5, zorder=8,
    )
    ax.annotate(
        r"Secant  $(c{=}2,\,q{=}\varphi)$",
        xy=(2, sec_E),
        xytext=(2.4, 1.08),
        fontsize=10, color="#2e7d32", fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#2e7d32", lw=1.4),
    )

    # Emphasize c = 2 column.
    ax.axvline(2, color="#888", lw=0.8, ls=":", alpha=0.5)

    ax.set_xlabel(r"function evaluations per step  $c$")
    ax.set_ylabel(r"efficiency index  $E(q,c) = q^{1/c}$")
    ax.set_title(
        "Kung-Traub frontier and Pandrosion-Steffensen at $c = 2$"
    )
    ax.set_xlim(1, 8)
    ax.set_ylim(1.0, 2.05)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", framealpha=0.95, fontsize=9)

    plt.savefig(OUT / "fig_efficiency_index.png")
    plt.close(fig)


def main():
    np.random.seed(42)
    fig_cyclotomic_basins()
    fig_loglog_convergence()
    fig_module_graph()
    fig_efficiency_index()
    print("Wrote 4 figures to", OUT)


if __name__ == "__main__":
    main()
