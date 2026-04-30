#!/usr/bin/env python3
"""Generate a publication-style figure for flow/119 over the 118 engine.

The script imports the 119 multifamily harness, runs a deterministic single
Kostlan family extraction on ks(2,10), and renders the probe-aware PURE Thales
pipeline: Mobius/Thales starts, StartOpt, theorem-guided probe selection,
exact telescopic slope assembly, residual validation, and extracted roots.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
DATA_DIR = ROOT / "verification" / "figure_119_probe_aware"
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


CASE = "2,10"
FAMILY = "ks"
COUNT = 8
POOL = 4096
EPOCHS = 24


def load_engine():
    path = ROOT / "flow" / "119_pandrosion_multifamily_probe_aware_pure_thales_engine.py"
    spec = importlib.util.spec_from_file_location("pandrosion_119_for_figure", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run_case(engine):
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    args = engine.build_parser().parse_args(
        [
            "--cases",
            CASE,
            "--families",
            FAMILY,
            "--count",
            str(COUNT),
            "--pool",
            str(POOL),
            "--epochs",
            str(EPOCHS),
            "--outdir",
            str(DATA_DIR),
        ]
    )
    return engine.run_case(args, CASE, FAMILY)


def setup_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_card(ax, xy, w, h, title, body, face, edge="#263342"):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.04",
        linewidth=1.2,
        facecolor=face,
        edgecolor=edge,
        alpha=0.98,
        path_effects=[pe.SimplePatchShadow(offset=(2.0, -2.0), alpha=0.16), pe.Normal()],
    )
    ax.add_patch(patch)
    ax.text(x + 0.03, y + h - 0.055, title, ha="left", va="top", fontsize=9.8, fontweight="bold", color="#17212b")
    ax.text(x + 0.03, y + h - 0.118, body, ha="left", va="top", fontsize=7.1, color="#344150", linespacing=1.12)
    return patch


def draw_arrow(ax, start, end, color="#263342", rad=0.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.45,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            alpha=0.9,
        )
    )


def pipeline_panel(ax):
    setup_axis(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.02, 0.96, "119 over 118: probe-aware flow", fontsize=13.6, fontweight="bold", color="#17212b")
    ax.text(0.02, 0.905, "Same geometric Thales start; the local endpoint b is selected by residual probes", fontsize=8.4, color="#52606d")

    colors = ["#f1c9ad", "#f1d06f", "#b9d8b6", "#9ed3d8", "#b8c7ef", "#dbc2ea"]
    nodes = [
        ((0.035, 0.63), 0.27, 0.19, "raw direction", "splitmix phases\ncomplex amplitudes"),
        ((0.365, 0.63), 0.27, 0.19, "Mobius chart", "Riemann angle theta\nmoving pole p"),
        ((0.695, 0.63), 0.27, 0.19, "Thales scale", "wide power ladder\nhomothety Lambda"),
        ((0.035, 0.36), 0.27, 0.19, "StartOpt", "radial candidates\nresidual scoring"),
        ((0.365, 0.36), 0.27, 0.19, "probe bank", "self, inertial,\nspiral endpoints"),
        ((0.695, 0.36), 0.27, 0.19, "accept", "residual threshold\ncluster roots"),
    ]
    for idx, (xy, w, h, title, body) in enumerate(nodes):
        draw_card(ax, xy, w, h, title, body, colors[idx])
    draw_arrow(ax, (0.305, 0.725), (0.365, 0.725))
    draw_arrow(ax, (0.635, 0.725), (0.695, 0.725))
    draw_arrow(ax, (0.83, 0.63), (0.83, 0.55), rad=-0.08)
    draw_arrow(ax, (0.305, 0.455), (0.365, 0.455))
    draw_arrow(ax, (0.635, 0.455), (0.695, 0.455))

    patch = FancyBboxPatch(
        (0.105, 0.08),
        0.79,
        0.15,
        boxstyle="round,pad=0.02,rounding_size=0.05",
        linewidth=1.25,
        facecolor="#eef0f4",
        edgecolor="#263342",
        path_effects=[pe.SimplePatchShadow(offset=(2.0, -2.0), alpha=0.16), pe.Normal()],
    )
    ax.add_patch(patch)
    ax.text(0.5, 0.165, "probe-aware PURE corrector", ha="center", va="center", fontsize=11.0, fontweight="bold", color="#17212b")
    ax.text(0.5, 0.103, r"$b_*=\arg\min_b \|G(b)\|$, build exact $Q_G(a,b_*)$, solve $Q_G\Delta=-G(a)$", ha="center", va="center", fontsize=8.1, color="#344150")
    draw_arrow(ax, (0.5, 0.36), (0.5, 0.23), color="#9a4d25")


def slope_probe_panel(ax, result):
    setup_axis(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.02, 0.96, "Finite-slope endpoint selection", fontsize=13.6, fontweight="bold", color="#17212b")
    ax.text(0.02, 0.905, "No Jacobian formula: probes only choose b for the same telescopic slope identity", fontsize=8.4, color="#52606d")
    draw_card(ax, (0.055, 0.60), 0.25, 0.21, "candidate b", "self endpoint\ninertial endpoint\ngeometric spiral", "#f1c9ad")
    draw_card(ax, (0.375, 0.60), 0.25, 0.21, "score", "minimize residual\nchoose best endpoint", "#b9d8b6")
    draw_card(ax, (0.695, 0.60), 0.25, 0.21, "solve", "exact slope matrix\nline search", "#b8c7ef")
    draw_arrow(ax, (0.305, 0.705), (0.375, 0.705))
    draw_arrow(ax, (0.625, 0.705), (0.695, 0.705))
    ax.text(0.50, 0.46, r"$F(b)-F(a)=Q_F(a,b)(b-a)$", ha="center", va="center", fontsize=16.0, fontweight="bold", color="#9a4d25")
    ax.text(0.50, 0.35, r"$Q=C\,T(a,b_*)$ with vectorized prefix/suffix monomial factors", ha="center", va="center", fontsize=9.0, color="#344150")

    trials = result.get("trials", [])
    probe_counts = [int(t.get("probe_total_evals") or 0) for t in trials if t.get("probe_total_evals") is not None]
    slope_counts = result["summary"]["eval_stats"].get("slope_count", 0)
    eval_counts = result["summary"]["eval_stats"].get("eval_count", 0)
    text = (
        f"probe evals in kept trials: {sum(probe_counts)}\n"
        f"slope calls: {slope_counts}\n"
        f"residual evals: {eval_counts}\n"
        f"corrector: {result['roots'][0].get('corrector', 'probe-aware') if result.get('roots') else 'probe-aware'}"
    )
    ax.text(0.50, 0.15, text, ha="center", va="center", fontsize=9.2, color="#17212b",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#fffaf0", edgecolor="#d2b48c"))


def support_panel(ax, system):
    exps = np.asarray(system.exps)
    totals = exps.sum(axis=1)
    weights = np.asarray(system.weights)
    ax.set_title("Kostlan support for ks(2,10)", loc="left", fontsize=12.0, fontweight="bold", color="#17212b", pad=10)
    ax.set_xticks(range(0, 11, 2))
    ax.set_yticks(range(0, 11, 2))
    for spine in ax.spines.values():
        spine.set_color("#d7d0c3")
    ax.grid(color="#e5ded2", linewidth=0.7)
    sizes = 22 + 68 * (weights / max(1e-12, weights.max())) ** 0.65
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
    ax.set_title("Coefficient matrix used by 119", loc="left", fontsize=12.0, fontweight="bold", color="#17212b", pad=10)
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
    ax.text(0.02, 0.86, r"same $C\,T$ slope assembly; b is probe-selected", transform=ax.transAxes, fontsize=8.3, color="white", fontweight="bold")


def roots_panel(ax, result):
    roots = result.get("roots", [])
    ax.set_title("Extracted roots, ks(2,10)", loc="left", fontsize=12.0, fontweight="bold", color="#17212b", pad=10)
    ax.axhline(0, color="#d7d0c3", linewidth=0.8)
    ax.axvline(0, color="#d7d0c3", linewidth=0.8)
    ax.grid(color="#eee7db", linewidth=0.7)
    if not roots:
        ax.text(0.5, 0.5, "no roots", transform=ax.transAxes, ha="center", va="center", fontsize=12, color="#7f1d1d")
        return
    z = np.array([[complex(coord[0], coord[1]) for coord in root["z"]] for root in roots], dtype=np.complex128)
    residuals = np.array([float(root["residual"]) for root in roots])
    scores = -np.log10(np.maximum(residuals, 1e-16))
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
    trials = result.get("trials", [])
    ax.set_title("Trial residuals and probe work", loc="left", fontsize=12.0, fontweight="bold", color="#17212b", pad=10)
    if not trials:
        ax.text(0.5, 0.5, "no trials stored", transform=ax.transAxes, ha="center", va="center", fontsize=12, color="#52606d")
        return
    xs = np.array([int(t["trial"]) for t in trials])
    rs = np.array([float(t.get("r1", math.inf)) for t in trials])
    ys = np.log10(np.maximum(rs, 1e-16))
    probes = np.array([float(t.get("probe_total_evals") or 0.0) for t in trials])
    colors = []
    for t in trials:
        if t.get("status") == "new-root":
            colors.append("#2c7a7b")
        elif t.get("status") == "duplicate":
            colors.append("#b7791f")
        elif t.get("accepted"):
            colors.append("#4c9f70")
        else:
            colors.append("#b84a62")
    sizes = 34 + 80 * probes / max(1.0, float(probes.max()))
    ax.scatter(xs, ys, s=sizes, c=colors, edgecolor="white", linewidth=0.55, zorder=3)
    ax.axhline(math.log10(result["parameters"]["accept"]), color="#27323f", linestyle="--", linewidth=1.1)
    ax.text(0.98, 0.13, "acceptance threshold", transform=ax.transAxes, ha="right", va="center", fontsize=8, color="#27323f")
    ax.set_xlabel("trial", fontsize=9)
    ax.text(0.015, 0.96, r"$\log_{10}\|F(z)\|$", transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color="#27323f")
    ax.grid(color="#eee7db", linewidth=0.7)
    for spine in ax.spines.values():
        spine.set_color("#d7d0c3")


def summary_panel(fig, result):
    s = result["summary"]
    stats = s["eval_stats"]
    best = result["roots"][0] if result.get("roots") else {"residual": float("nan")}
    probe_sum = sum(int(t.get("probe_total_evals") or 0) for t in result.get("trials", []))
    text = (
        f"roots {s['unique_roots']}/{s['requested_roots']}   "
        f"trials {s['trials_used']}   "
        f"Q calls {stats['slope_count']}   "
        f"probe evals {probe_sum}   "
        f"best residual {float(best['residual']):.2e}   "
        f"total {s['total_seconds']:.3f}s"
    )
    fig.text(
        0.5,
        0.035,
        text,
        ha="center",
        va="center",
        fontsize=10.2,
        color="#17212b",
        bbox=dict(boxstyle="round,pad=0.42,rounding_size=0.12", facecolor="#fffaf0", edgecolor="#d2b48c", linewidth=1.0),
    )


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    engine = load_engine()
    result = run_case(engine)
    system = engine.MultiFamilySystem.make(2, 10, seed_index=0, equation_normalize=False, family=FAMILY)

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
    fig.suptitle("Engine 119(118): probe-aware pure Thales Pandrosion", x=0.5, y=0.975, fontsize=21.5, fontweight="bold", color="#17212b")
    fig.text(0.5, 0.94, "Exact telescopic slopes with theorem-guided finite endpoint probes", ha="center", fontsize=11.0, color="#52606d")

    gs = fig.add_gridspec(2, 3, left=0.045, right=0.985, top=0.895, bottom=0.09, hspace=0.36, wspace=0.30)
    pipeline_panel(fig.add_subplot(gs[0, 0]))
    slope_probe_panel(fig.add_subplot(gs[0, 1]), result)
    roots_panel(fig.add_subplot(gs[0, 2]), result)
    support_panel(fig.add_subplot(gs[1, 0]), system)
    coefficient_panel(fig.add_subplot(gs[1, 1]), system)
    trial_panel(fig.add_subplot(gs[1, 2]), result)
    summary_panel(fig, result)

    png = FIG_DIR / "119_probe_aware_pure_thales_pandrosion.png"
    pdf = FIG_DIR / "119_probe_aware_pure_thales_pandrosion.pdf"
    fig.savefig(png, dpi=260)
    fig.savefig(pdf)
    plt.close(fig)
    print(png)
    print(pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
