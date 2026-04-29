#!/usr/bin/env python3
"""Generate a publication-style figure for flow/115.

The script imports the 115 engine, runs a deterministic ks(2,10) extraction,
and renders a compact visual summary of the pipeline, the vectorized slope
construction, the monomial support, the coefficient matrix, and the roots found.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
MPL_DIR = Path("/tmp/poussiere-mplconfig")
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR / "xdg-cache"))
MPL_DIR.mkdir(parents=True, exist_ok=True)
(MPL_DIR / "xdg-cache").mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


def load_engine():
    path = ROOT / "flow" / "115_pandrosion_vectorized_pure_pandrosion.py"
    spec = importlib.util.spec_from_file_location("pandrosion_115_for_figure", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_run(engine):
    args = engine.build_parser().parse_args(
        [
            "--cases",
            "2,10",
            "--count",
            "8",
            "--pool",
            "4096",
            "--outdir",
            str(ROOT / "verification" / "figure_115_tmp"),
        ]
    )
    return engine.run_case(args, "2,10")


def setup_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_node(ax, xy, w, h, title, body, face, edge="#27323f"):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.045",
        linewidth=1.25,
        facecolor=face,
        edgecolor=edge,
        alpha=0.98,
        path_effects=[pe.SimplePatchShadow(offset=(2.2, -2.2), alpha=0.18), pe.Normal()],
    )
    ax.add_patch(patch)
    ax.text(x + 0.035, y + h - 0.065, title, ha="left", va="top", fontsize=10.0, fontweight="bold", color="#17212b")
    ax.text(x + 0.035, y + h - 0.135, body, ha="left", va="top", fontsize=7.2, color="#344150", linespacing=1.1)
    return patch


def draw_arrow(ax, start, end, color="#243447", rad=0.0):
    arr = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=13,
        linewidth=1.55,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
        alpha=0.9,
    )
    ax.add_patch(arr)


def pipeline_panel(ax):
    setup_axis(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.02, 0.965, "115 single geometric flow", fontsize=13.5, fontweight="bold", color="#17212b")
    ax.text(0.02, 0.91, "A direct route from deterministic starts to certified residuals", fontsize=8.6, color="#52606d")

    colors = ["#f6d7bf", "#f0c77e", "#b8d8ba", "#9bd4d7", "#b9c7f5", "#d9c2ee"]
    nodes = [
        ((0.035, 0.63), 0.27, 0.20, "raw direction", "splitmix phases\ncomplex amplitudes"),
        ((0.365, 0.63), 0.27, 0.20, "Mobius chart", "Riemann angle theta\nmoving pole p"),
        ((0.695, 0.63), 0.27, 0.20, "Thales scale", "wide power ladder\nhomothety Lambda"),
        ((0.035, 0.36), 0.27, 0.20, "StartOpt", "radial search\ngeometric descent"),
        ((0.365, 0.36), 0.27, 0.20, "linear chart", "map y to z\nA = identity"),
        ((0.695, 0.36), 0.27, 0.20, "validation", "residual check\ndeduplicate"),
    ]
    centers = []
    for idx, (xy, w, h, title, body) in enumerate(nodes):
        draw_node(ax, xy, w, h, title, body, colors[idx])
        centers.append((xy[0] + w / 2, xy[1] + h / 2))

    draw_arrow(ax, (0.305, 0.73), (0.365, 0.73))
    draw_arrow(ax, (0.635, 0.73), (0.695, 0.73))
    draw_arrow(ax, (0.83, 0.63), (0.83, 0.56), rad=-0.08)
    draw_arrow(ax, (0.695, 0.46), (0.635, 0.46))
    draw_arrow(ax, (0.365, 0.46), (0.305, 0.46))

    q_patch = FancyBboxPatch(
        (0.12, 0.08),
        0.76,
        0.15,
        boxstyle="round,pad=0.02,rounding_size=0.05",
        linewidth=1.3,
        facecolor="#eef0f4",
        edgecolor="#243447",
        path_effects=[pe.SimplePatchShadow(offset=(2.0, -2.0), alpha=0.16), pe.Normal()],
    )
    ax.add_patch(q_patch)
    ax.text(0.5, 0.165, "PURE corrector", ha="center", va="center", fontsize=11.2, fontweight="bold", color="#17212b")
    ax.text(0.5, 0.103, r"build exact slope $Q_G(a,b)$, solve $Q_G\Delta=-G(a)$, damp by Thales line-search", ha="center", va="center", fontsize=8.1, color="#344150")
    draw_arrow(ax, (0.5, 0.36), (0.5, 0.23), color="#b35a2a")


def slope_panel(ax):
    setup_axis(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.02, 0.965, "Vectorized exact slope", fontsize=13.5, fontweight="bold", color="#17212b")
    ax.text(0.02, 0.91, r"No Jacobian in the local corrector: $F(b)-F(a)=Q(a,b)(b-a)$", fontsize=8.6, color="#52606d")

    blocks = [
        (0.07, 0.64, 0.18, 0.18, "#f3b391", "powers", r"$a^m,\ b^m$"),
        (0.31, 0.64, 0.18, 0.18, "#b7d8a7", "prefix", r"$k<j$ factors"),
        (0.55, 0.64, 0.18, 0.18, "#9ccbd1", "suffix", r"$k>j$ factors"),
        (0.79, 0.64, 0.18, 0.18, "#d8c0e9", "slope", r"$S_m$ recurrence"),
    ]
    for x, y, w, h, color, title, body in blocks:
        draw_node(ax, (x, y), w, h, title, body, color)
    for x in [0.25, 0.49, 0.73]:
        draw_arrow(ax, (x, 0.73), (x + 0.055, 0.73), color="#455467")

    ax.text(0.50, 0.49, r"$T_{\cdot j}=P_j^{<}(b)\,P_j^{>}(a)\,S_{j,\alpha_j}(a_j,b_j)$", ha="center", va="center", fontsize=13.0, color="#17212b")
    ax.text(0.50, 0.37, r"$Q_F(a,b)=C\,T$", ha="center", va="center", fontsize=18.5, fontweight="bold", color="#9a4d25")

    matrix_y = 0.12
    ax.imshow(np.outer(np.linspace(0.2, 1.0, 18), np.linspace(1.0, 0.2, 36)), extent=(0.13, 0.43, matrix_y, matrix_y + 0.16), cmap="YlOrBr", aspect="auto", alpha=0.92)
    ax.imshow(np.outer(np.sin(np.linspace(0, 4, 36)) ** 2 + 0.1, np.cos(np.linspace(0, 3, 18)) ** 2 + 0.1).T, extent=(0.57, 0.87, matrix_y, matrix_y + 0.16), cmap="PuBuGn", aspect="auto", alpha=0.92)
    ax.text(0.28, matrix_y - 0.035, r"coefficient matrix $C$", ha="center", va="top", fontsize=8.3, color="#344150")
    ax.text(0.72, matrix_y - 0.035, r"term block $T$", ha="center", va="top", fontsize=8.3, color="#344150")
    draw_arrow(ax, (0.44, matrix_y + 0.08), (0.56, matrix_y + 0.08), color="#9a4d25")


def support_panel(ax, system):
    exps = np.asarray(system.exps)
    totals = exps.sum(axis=1)
    weights = np.asarray(system.weights)
    setup_axis(ax)
    ax.set_facecolor("#fbf8f2")
    ax.set_title("Kostlan support for ks(2,10)", loc="left", fontsize=12.0, fontweight="bold", color="#17212b", pad=10)
    ax.set_xticks(range(0, 11, 2))
    ax.set_yticks(range(0, 11, 2))
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#d7d0c3")
    ax.grid(color="#e5ded2", linewidth=0.7)
    sizes = 22 + 65 * (weights / max(1e-12, weights.max())) ** 0.65
    sc = ax.scatter(exps[:, 0], exps[:, 1], c=totals, s=sizes, cmap="viridis", edgecolor="white", linewidth=0.55, zorder=3)
    ax.set_xlim(-0.7, 10.7)
    ax.set_ylim(-0.7, 10.7)
    ax.set_xlabel(r"$\alpha_1$", fontsize=9)
    ax.set_ylabel(r"$\alpha_2$", fontsize=9)
    cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(r"$|\alpha|$", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.text(0.02, 0.03, f"{system.terms_per_poly} monomials per equation", transform=ax.transAxes, fontsize=8.5, color="#52606d")


def coefficient_panel(ax, system):
    mag = np.log10(1e-14 + np.abs(system.coeff))
    ax.set_title("Vectorized dense coefficient block", loc="left", fontsize=12.0, fontweight="bold", color="#17212b", pad=10)
    im = ax.imshow(mag, aspect="auto", cmap="magma", interpolation="nearest")
    ax.set_yticks([0, 1])
    ax.set_yticklabels([r"$F_1$", r"$F_2$"], fontsize=8)
    ax.set_xlabel("monomial index", fontsize=9)
    ax.set_xticks([0, 15, 30, 45, 60])
    for spine in ax.spines.values():
        spine.set_color("#d7d0c3")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(r"$\log_{10}|C_{i,\alpha}|$", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.text(0.02, 0.86, r"BLAS-backed product: $Q=C\,T$", transform=ax.transAxes, fontsize=8.8, color="white", fontweight="bold")


def roots_panel(ax, result):
    roots = result["roots"]
    z = np.array([[complex(coord[0], coord[1]) for coord in root["z"]] for root in roots], dtype=np.complex128)
    residuals = np.array([float(root["residual"]) for root in roots])
    scores = -np.log10(np.maximum(residuals, 1e-16))
    ax.set_title("Extracted roots, ks(2,10)", loc="left", fontsize=12.0, fontweight="bold", color="#17212b", pad=10)
    ax.axhline(0, color="#d7d0c3", linewidth=0.8)
    ax.axvline(0, color="#d7d0c3", linewidth=0.8)
    ax.grid(color="#eee7db", linewidth=0.7)
    ax.scatter(z[:, 0].real, z[:, 0].imag, s=145, c=scores, cmap="coolwarm", edgecolor="#1f2933", linewidth=0.85, zorder=4)
    for idx, value in enumerate(z[:, 0]):
        ax.text(value.real, value.imag, str(idx + 1), ha="center", va="center", fontsize=7.0, color="white", fontweight="bold", zorder=5)
    ax.set_xlabel(r"Re $z_1$", fontsize=9)
    ax.set_ylabel(r"Im $z_1$", fontsize=9)
    ax.set_aspect("equal", adjustable="datalim")
    margin = 0.2
    ax.set_xlim(float(np.min(z[:, 0].real)) - margin, float(np.max(z[:, 0].real)) + margin)
    ax.set_ylim(float(np.min(z[:, 0].imag)) - margin, float(np.max(z[:, 0].imag)) + margin)
    sm = plt.cm.ScalarMappable(cmap="coolwarm")
    sm.set_array(scores)
    cb = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(r"$-\log_{10}$ residual", fontsize=8)
    cb.ax.tick_params(labelsize=7)


def trial_panel(ax, result):
    trials = result["trials"]
    xs = np.array([int(t["trial"]) for t in trials])
    rs = np.array([float(t.get("r1", math.inf)) for t in trials])
    ys = np.log10(np.maximum(rs, 1e-16))
    status_color = {
        "new-root": "#2c7a7b",
        "duplicate": "#b7791f",
        "converged": "#2c7a7b",
    }
    colors = []
    for t in trials:
        if t.get("status") == "new-root":
            colors.append(status_color["new-root"])
        elif t.get("status") == "duplicate":
            colors.append(status_color["duplicate"])
        elif t.get("accepted"):
            colors.append("#4c9f70")
        else:
            colors.append("#b84a62")
    ax.set_title("Trial residual trace", loc="left", fontsize=12.0, fontweight="bold", color="#17212b", pad=10)
    ax.scatter(xs, ys, s=46, c=colors, edgecolor="white", linewidth=0.55, zorder=3)
    ax.axhline(math.log10(result["parameters"]["accept"]), color="#27323f", linestyle="--", linewidth=1.1)
    ax.text(0.98, 0.13, "acceptance threshold", transform=ax.transAxes, ha="right", va="center", fontsize=8, color="#27323f")
    ax.set_xlabel("trial", fontsize=9)
    ax.set_ylabel("")
    ax.text(0.015, 0.96, r"$\log_{10}\|F(z)\|$", transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color="#27323f")
    ax.grid(color="#eee7db", linewidth=0.7)
    for spine in ax.spines.values():
        spine.set_color("#d7d0c3")


def summary_panel(fig, result):
    s = result["summary"]
    stats = s["eval_stats"]
    best = result["roots"][0]
    text = (
        f"roots {s['unique_roots']}/{s['requested_roots']}   "
        f"trials {s['trials_used']}   "
        f"Q calls {stats['slope_count']}   "
        f"Q time {stats['seconds_slope']:.3f}s   "
        f"best residual {float(best['residual']):.2e}   "
        f"total {s['total_seconds']:.3f}s"
    )
    fig.text(
        0.5,
        0.035,
        text,
        ha="center",
        va="center",
        fontsize=10.3,
        color="#17212b",
        bbox=dict(boxstyle="round,pad=0.42,rounding_size=0.12", facecolor="#fffaf0", edgecolor="#d2b48c", linewidth=1.0),
    )


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    engine = load_engine()
    result = parse_run(engine)
    system = engine.DenseKostlanSystem.make(2, 10, seed_index=0, equation_normalize=False)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "figure.facecolor": "#f8f4ec",
            "axes.facecolor": "#fffdf8",
            "axes.edgecolor": "#d7d0c3",
            "axes.labelcolor": "#27323f",
            "xtick.color": "#52606d",
            "ytick.color": "#52606d",
            "savefig.facecolor": "#f8f4ec",
            "savefig.edgecolor": "#f8f4ec",
        }
    )

    fig = plt.figure(figsize=(18.0, 10.5), dpi=220)
    fig.suptitle("Engine 115: vectorized pure Pandrosion", x=0.5, y=0.975, fontsize=22, fontweight="bold", color="#17212b")
    fig.text(0.5, 0.94, "Exact telescopic slopes, BLAS-level assembly, and geometric Thales starts", ha="center", fontsize=11.2, color="#52606d")

    gs = fig.add_gridspec(2, 3, left=0.045, right=0.985, top=0.895, bottom=0.09, hspace=0.36, wspace=0.30)
    pipeline_panel(fig.add_subplot(gs[0, 0]))
    slope_panel(fig.add_subplot(gs[0, 1]))
    roots_panel(fig.add_subplot(gs[0, 2]), result)
    support_panel(fig.add_subplot(gs[1, 0]), system)
    coefficient_panel(fig.add_subplot(gs[1, 1]), system)
    trial_panel(fig.add_subplot(gs[1, 2]), result)
    summary_panel(fig, result)

    png = FIG_DIR / "115_pandrosion_vectorized_pure_pandrosion.png"
    pdf = FIG_DIR / "115_pandrosion_vectorized_pure_pandrosion.pdf"
    fig.savefig(png, dpi=260)
    fig.savefig(pdf)
    plt.close(fig)
    print(png)
    print(pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
