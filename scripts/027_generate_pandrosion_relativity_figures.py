from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle


OUT = Path("latex/figures")


def setup() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.facecolor": "#fbfbf7",
            "figure.facecolor": "#fbfbf7",
            "savefig.facecolor": "#fbfbf7",
            "axes.edgecolor": "#253035",
            "axes.labelcolor": "#1d2327",
            "xtick.color": "#293136",
            "ytick.color": "#293136",
            "text.color": "#182126",
            "axes.titleweight": "bold",
            "axes.titlesize": 15,
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=220)
    plt.close(fig)


def fig_causal_probe_lattice() -> None:
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    ax.set_xlim(-6.3, 6.3)
    ax.set_ylim(-0.25, 6.6)
    ax.set_aspect("equal", adjustable="box")

    tmax = 6.2
    cone = Polygon(
        [(0, 0), (-tmax, tmax), (tmax, tmax)],
        closed=True,
        facecolor="#a8dadc",
        edgecolor="none",
        alpha=0.34,
        zorder=0,
    )
    ax.add_patch(cone)
    ax.plot([-tmax, 0, tmax], [tmax, 0, tmax], color="#1d3557", lw=1.7, ls="--")
    ax.text(-4.95, 5.35, "future light cone", fontsize=11, color="#1d3557")

    xs: list[float] = []
    ts: list[float] = []
    for t in np.arange(0.8, 6.1, 0.6):
        for x in np.arange(-5.4, 5.6, 0.6):
            if abs(x) < t - 0.14:
                xs.append(float(x))
                ts.append(float(t))
    ax.scatter(xs, ts, s=13, color="#457b9d", alpha=0.34, edgecolors="none", label="causal candidate events")

    tau = np.linspace(0.0, 6.0, 260)
    x_path = 0.18 * tau + 0.045 * tau**2
    ax.plot(x_path, tau, color="#d1495b", lw=3.0, label="selected Pandrosion worldline")

    anchor_t = 2.4
    anchor_x = 0.18 * anchor_t + 0.045 * anchor_t**2
    ax.scatter([anchor_x], [anchor_t], s=70, color="#d1495b", edgecolor="#111827", zorder=5)
    next_t = anchor_t + 0.8
    candidates = np.array([-0.35, 0.05, 0.45, 0.85, 1.25]) + anchor_x
    for c in candidates:
        ax.add_patch(
            FancyArrowPatch(
                (anchor_x, anchor_t),
                (c, next_t),
                arrowstyle="-|>",
                mutation_scale=10,
                lw=1.1,
                color="#6d597a",
                alpha=0.68,
            )
        )
    ax.scatter(candidates, np.full_like(candidates, next_t), s=45, color="#f4a261", edgecolor="#111827", zorder=4)
    ax.text(anchor_x + 0.25, anchor_t - 0.28, "current event", fontsize=10)
    ax.text(anchor_x + 1.05, next_t + 0.16, "finite-slope probes", fontsize=10, color="#6d597a")

    ax.axhline(0, color="#253035", lw=1.0)
    ax.axvline(0, color="#253035", lw=1.0)
    ax.set_xlabel("space coordinate x")
    ax.set_ylabel("time coordinate ct")
    ax.set_title("Causal Pandrosion probes live inside the future cone")
    ax.grid(True, color="#d8d6ca", alpha=0.6)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", loc="upper left")
    save(fig, "027_causal_probe_lattice")


def fig_dynamic_atlas_horizon() -> None:
    r = np.linspace(1.28, 5.3, 900)
    outside = np.maximum(r - 2.0, 1e-3)
    schwarzschild = 1.0 + 0.32 / outside**2
    schwarzschild = np.clip(schwarzschild, 1.0, 18.0)
    ef_chart = 1.18 + 0.18 * np.exp(-((r - 2.0) / 0.75) ** 2)
    weight = 1.0 / (1.0 + np.exp(-7.5 * (r - 2.55)))
    dynamic = weight * schwarzschild + (1.0 - weight) * ef_chart + 0.12 * weight * (1.0 - weight)

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 7.0), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.15]})

    ax = axes[0]
    ax.axvspan(1.78, 2.36, color="#e9c46a", alpha=0.25, label="atlas transition zone")
    ax.axvline(2.0, color="#1d3557", lw=1.4, ls="--", label="r = 2M")
    ax.plot(r, schwarzschild, color="#d1495b", lw=2.0, label="singular coordinate cost")
    ax.plot(r, ef_chart, color="#2a9d8f", lw=2.0, label="regular infalling chart cost")
    ax.plot(r, dynamic, color="#111827", lw=2.4, label="dynamic Pandrosion atlas cost")
    ax.set_ylim(0.8, 12.0)
    ax.set_ylabel("chart condition cost")
    ax.set_title("Dynamic atlas avoids a coordinate singularity at the horizon")
    ax.grid(True, color="#d8d6ca", alpha=0.65)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", loc="upper right")

    ax2 = axes[1]
    ax2.plot(r, 1.0 - weight, color="#2a9d8f", lw=2.2, label="regular chart weight")
    ax2.plot(r, weight, color="#d1495b", lw=2.2, label="exterior chart weight")
    ax2.axvline(2.0, color="#1d3557", lw=1.2, ls="--")
    ax2.set_ylim(-0.06, 1.06)
    ax2.set_xlabel("areal radius r / M")
    ax2.set_ylabel("chart weight")
    ax2.grid(True, color="#d8d6ca", alpha=0.65)
    ax2.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", loc="center right")
    save(fig, "027_dynamic_atlas_horizon")


def fig_finite_holonomy_curvature() -> None:
    x = np.linspace(-2.8, 2.8, 240)
    y = np.linspace(-2.2, 2.2, 200)
    X, Y = np.meshgrid(x, y)
    K = 0.95 * np.exp(-0.55 * (X**2 + 1.35 * Y**2)) - 0.18 * np.exp(-2.2 * ((X - 1.3) ** 2 + (Y + 0.55) ** 2))

    fig, ax = plt.subplots(figsize=(8.9, 6.2))
    im = ax.contourf(X, Y, K, levels=22, cmap="viridis", alpha=0.82)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("finite curvature density")

    # Draw a small closed causal/spatial probe loop and a transported vector.
    loop = np.array([[-0.95, -0.42], [0.82, -0.25], [1.04, 0.82], [-0.72, 0.65], [-0.95, -0.42]])
    ax.plot(loop[:, 0], loop[:, 1], color="#fefae0", lw=4.2, solid_capstyle="round")
    ax.plot(loop[:, 0], loop[:, 1], color="#111827", lw=1.6, solid_capstyle="round")
    for a, b in zip(loop[:-1], loop[1:]):
        mid = 0.56 * a + 0.44 * b
        ax.add_patch(
            FancyArrowPatch(
                tuple(mid),
                tuple(0.72 * b + 0.28 * a),
                arrowstyle="-|>",
                mutation_scale=13,
                lw=1.5,
                color="#111827",
            )
        )

    origin = np.array([-0.58, -0.12])
    v0 = np.array([0.78, 0.0])
    v1 = np.array([0.68, 0.33])
    ax.add_patch(FancyArrowPatch(tuple(origin), tuple(origin + v0), arrowstyle="-|>", mutation_scale=17, lw=2.2, color="#f4a261"))
    ax.add_patch(FancyArrowPatch(tuple(origin + [0.05, 0.14]), tuple(origin + [0.05, 0.14] + v1), arrowstyle="-|>", mutation_scale=17, lw=2.2, color="#d1495b"))
    ax.text(origin[0] + 0.84, origin[1] - 0.1, "initial vector", fontsize=10, color="#fefae0")
    ax.text(origin[0] + 0.74, origin[1] + 0.55, "after loop", fontsize=10, color="#fefae0")

    rect = Rectangle((-2.58, 1.35), 2.05, 0.62, facecolor="#fbfbf7", edgecolor="#111827", lw=1.0, alpha=0.92)
    ax.add_patch(rect)
    ax.text(-1.555, 1.66, r"$R_h(u,v)w \approx (\Pi_h w-w)/Area_h$", ha="center", va="center", fontsize=11)
    ax.set_xlabel("local spatial probe coordinate 1")
    ax.set_ylabel("local spatial probe coordinate 2")
    ax.set_title("Finite holonomy: curvature measured by a closed Pandrosion probe loop")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-2.8, 2.8)
    ax.set_ylim(-2.2, 2.2)
    save(fig, "027_finite_holonomy_curvature")


def main() -> None:
    setup()
    fig_causal_probe_lattice()
    fig_dynamic_atlas_horizon()
    fig_finite_holonomy_curvature()
    for path in sorted(OUT.glob("027_*")):
        print(path)


if __name__ == "__main__":
    main()
