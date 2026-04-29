#!/usr/bin/env python3
"""Generate a multifamily multivariate figure for flow/116.

The script imports flow/116, runs a deterministic all-family benchmark on
ks(3,4), and renders a publication-style summary.  The case is deliberately
multivariate (n=3) and the family list is exactly the 116 "all" group.
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
DATA_DIR = ROOT / "verification" / "figure_116_multifamily"
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
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


CASE = "3,4"
COUNT = 4
POOL = 512


def load_engine():
    path = ROOT / "flow" / "116_pandrosion_multifamily_vectorized_pure_pandrosion.py"
    spec = importlib.util.spec_from_file_location("pandrosion_116_for_figure", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run_benchmark(engine):
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    args = engine.build_parser().parse_args(
        [
            "--cases",
            CASE,
            "--families",
            "all",
            "--count",
            str(COUNT),
            "--pool",
            str(POOL),
            "--outdir",
            str(DATA_DIR),
        ]
    )
    families = engine.parse_families("all")
    runs = [engine.run_case(args, CASE, family) for family in families]
    data = {
        "script": "116_pandrosion_multifamily_vectorized_pure_pandrosion.py",
        "case": CASE,
        "families": families,
        "runs": runs,
        "summary": {
            "runs": len(runs),
            "successes": int(sum(1 for r in runs if r["summary"]["success"])),
            "total_roots": int(sum(r["summary"]["unique_roots"] for r in runs)),
            "requested_roots_per_family": COUNT,
            "pool": POOL,
        },
    }
    (DATA_DIR / "116_multifamily_3x4_all_for_figure.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    return data


def family_short(name: str) -> str:
    return {
        "degree_shell_ks": "shell",
        "mixed_degree": "mixed",
        "sparse_iid": "sparse iid",
        "dense_iid": "dense iid",
        "ill_scaled": "ill",
    }.get(name, name.replace("_", " "))


def setup_panel(ax):
    for spine in ax.spines.values():
        spine.set_color("#d2cabd")
    ax.tick_params(colors="#52606d", labelsize=8.5)
    ax.set_facecolor("#fffdf8")
    ax.grid(color="#eee6d8", linewidth=0.75, zorder=0)


def card(ax, x, y, w, h, title, body, face, edge="#253140"):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.035",
        linewidth=1.15,
        facecolor=face,
        edgecolor=edge,
        alpha=0.98,
        path_effects=[pe.SimplePatchShadow(offset=(2.0, -2.0), alpha=0.17), pe.Normal()],
    )
    ax.add_patch(patch)
    ax.text(x + 0.025, y + h - 0.045, title, ha="left", va="top", fontsize=9.4, fontweight="bold", color="#17212b")
    ax.text(x + 0.025, y + h - 0.105, body, ha="left", va="top", fontsize=7.1, color="#344150", linespacing=1.12)
    return patch


def arrow(ax, start, end, color="#253140", rad=0.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.35,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            alpha=0.88,
        )
    )


def architecture_panel(ax):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.02, 0.95, "116 multifamily harness", fontsize=13.8, fontweight="bold", color="#17212b")
    ax.text(0.02, 0.895, "Same 115 vectorized pure corrector; generator family changes", fontsize=8.7, color="#52606d")
    card(ax, 0.04, 0.63, 0.22, 0.17, "case", "n = 3, d = 4\nBezout scale 64", "#f4d5bd")
    card(ax, 0.31, 0.63, 0.25, 0.17, "families all", "dense, sparse,\nstructured, ill-scaled", "#f0c96f")
    card(ax, 0.62, 0.63, 0.32, 0.17, "115 engine", "Mobius-Thales starts\nexact slope Q = C T", "#b7d9c0")
    arrow(ax, (0.26, 0.715), (0.31, 0.715))
    arrow(ax, (0.56, 0.715), (0.62, 0.715))

    card(ax, 0.07, 0.32, 0.25, 0.17, "support", "exponents and weights\nper family", "#9ed4d9")
    card(ax, 0.38, 0.32, 0.25, 0.17, "correct", "slope endpoint b\nsolve Q_G Delta = -G", "#b8c7ef")
    card(ax, 0.69, 0.32, 0.25, 0.17, "accept", "residual threshold\ncluster roots", "#d7c1eb")
    arrow(ax, (0.32, 0.405), (0.38, 0.405))
    arrow(ax, (0.63, 0.405), (0.69, 0.405))

    box = FancyBboxPatch(
        (0.16, 0.08),
        0.68,
        0.12,
        boxstyle="round,pad=0.018,rounding_size=0.04",
        linewidth=1.15,
        facecolor="#eef0f4",
        edgecolor="#253140",
    )
    ax.add_patch(box)
    ax.text(0.5, 0.14, "one benchmark, thirteen polynomial landscapes", ha="center", va="center", fontsize=10.2, fontweight="bold", color="#17212b")


def performance_panel(ax, runs):
    setup_panel(ax)
    names = [family_short(r["family"]) for r in runs]
    roots = np.array([r["summary"]["unique_roots"] for r in runs], dtype=float)
    requested = np.array([r["summary"]["requested_roots"] for r in runs], dtype=float)
    failures = np.array([r["summary"]["failures"] for r in runs], dtype=float)
    trials = np.array([max(1, r["summary"]["trials_used"]) for r in runs], dtype=float)
    x = np.arange(len(runs))
    colors = ["#2c7a7b" if r >= q else "#b84a62" for r, q in zip(roots, requested)]
    ax.bar(x, requested, color="#e9dfcf", edgecolor="#d2cabd", linewidth=0.8, label="requested")
    ax.bar(x, roots, color=colors, edgecolor="#253140", linewidth=0.55, label="found", zorder=3)
    ax2 = ax.twinx()
    ax2.plot(x, failures / trials, color="#9a4d25", marker="o", linewidth=1.6, markersize=4.2, label="failure rate")
    ax2.set_ylim(-0.04, 1.04)
    ax2.set_ylabel("failure rate", fontsize=8.5, color="#9a4d25")
    ax2.tick_params(colors="#9a4d25", labelsize=8)
    for spine in ax2.spines.values():
        spine.set_color("#d2cabd")
    ax.set_title("All families on ks(3,4)", loc="left", fontsize=12.2, fontweight="bold", color="#17212b", pad=8)
    ax.set_ylabel("roots found", fontsize=8.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=42, ha="right", fontsize=7.8)
    ax.set_ylim(0, max(requested) + 1.0)
    ax.text(0.02, 0.92, f"count={COUNT}, pool={POOL}", transform=ax.transAxes, fontsize=8.2, color="#52606d")


def root_cloud_panel(ax, runs):
    ax.set_title("3D root cloud: real projection", loc="left", fontsize=12.2, fontweight="bold", color="#17212b", pad=8)
    cmap = plt.get_cmap("tab20")
    xs, ys, zs, cs, ss = [], [], [], [], []
    for i, run in enumerate(runs):
        for root in run.get("roots", []):
            coords = [complex(v[0], v[1]) for v in root["z"]]
            if len(coords) < 3:
                continue
            xs.append(coords[0].real)
            ys.append(coords[1].real)
            zs.append(coords[2].real)
            cs.append(cmap(i % 20))
            ss.append(30 + 15 * max(0.0, -math.log10(max(float(root["residual"]), 1e-16)) - 7.0))
    if xs:
        ax.scatter(xs, ys, zs, s=ss, c=cs, edgecolor="#1f2933", linewidth=0.35, alpha=0.86)
    ax.set_xlabel("Re z1", fontsize=8, labelpad=-1)
    ax.set_ylabel("Re z2", fontsize=8, labelpad=-1)
    ax.set_zlabel("Re z3", fontsize=8, labelpad=-1)
    ax.tick_params(labelsize=7.2, colors="#52606d")
    ax.view_init(elev=23, azim=-47)
    ax.xaxis.pane.set_facecolor((1.0, 0.99, 0.96, 1.0))
    ax.yaxis.pane.set_facecolor((1.0, 0.99, 0.96, 1.0))
    ax.zaxis.pane.set_facecolor((1.0, 0.99, 0.96, 1.0))
    ax.grid(True, color="#eee6d8")


def simplex_panel(ax, engine):
    setup_panel(ax)
    system = engine.MultiFamilySystem.make(3, 4, family="ks", seed_index=0, equation_normalize=False)
    exps = np.asarray(system.exps)
    weights = np.asarray(system.weights)
    # Barycentric coordinates on the total-degree <= 4 simplex.
    x = exps[:, 0] + 0.5 * exps[:, 1]
    y = (math.sqrt(3) / 2.0) * exps[:, 1]
    totals = exps.sum(axis=1)
    sizes = 28 + 105 * (weights / max(1e-12, weights.max())) ** 0.72
    triangle = np.array([[0, 0], [4, 0], [2, 2 * math.sqrt(3)]])
    ax.add_patch(Polygon(triangle, closed=True, facecolor="#f7efe2", edgecolor="#d2cabd", linewidth=1.2, zorder=1))
    sc = ax.scatter(x, y, c=totals, s=sizes, cmap="viridis", edgecolor="white", linewidth=0.55, zorder=3)
    ax.set_title("3-variable exponent simplex", loc="left", fontsize=12.2, fontweight="bold", color="#17212b", pad=8)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.text(0, -0.25, r"$\alpha_1$", ha="center", fontsize=8.5, color="#52606d")
    ax.text(4, -0.25, r"$\alpha_2$", ha="center", fontsize=8.5, color="#52606d")
    ax.text(2, 2 * math.sqrt(3) + 0.15, r"$\alpha_3$", ha="center", fontsize=8.5, color="#52606d")
    ax.text(0.02, 0.04, f"{system.terms_per_poly} monomials", transform=ax.transAxes, fontsize=8.3, color="#52606d")
    cb = plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.015)
    cb.set_label(r"$|\alpha|$", fontsize=8)
    cb.ax.tick_params(labelsize=7.2)


def residual_panel(ax, runs):
    setup_panel(ax)
    names = [family_short(r["family"]) for r in runs]
    max_roots = max([len(r.get("roots", [])) for r in runs] + [1])
    mat = np.full((len(runs), max_roots), np.nan)
    for i, run in enumerate(runs):
        for j, root in enumerate(run.get("roots", [])):
            mat[i, j] = -math.log10(max(float(root["residual"]), 1e-16))
    cmap = LinearSegmentedColormap.from_list("residuals", ["#f3e4cd", "#ef9a62", "#b83251", "#3b1f5f"])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, interpolation="nearest", vmin=7.5, vmax=10.0)
    ax.set_title("Residual quality by family", loc="left", fontsize=12.2, fontweight="bold", color="#17212b", pad=8)
    ax.set_yticks(np.arange(len(runs)))
    ax.set_yticklabels(names, fontsize=7.5)
    ax.set_xticks(np.arange(max_roots))
    ax.set_xticklabels([str(i + 1) for i in range(max_roots)], fontsize=8)
    ax.set_xlabel("root rank", fontsize=8.6)
    ax.set_facecolor("#efe6d8")
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.015)
    cb.set_label(r"$-\log_{10}$ residual", fontsize=8)
    cb.ax.tick_params(labelsize=7.2)
    for i, run in enumerate(runs):
        if len(run.get("roots", [])) == 0:
            ax.text(0, i, "no roots", ha="center", va="center", fontsize=7.5, color="#7f1d1d", fontweight="bold")


def cost_panel(ax, runs):
    setup_panel(ax)
    names = [family_short(r["family"]) for r in runs]
    active = np.array([r["active_terms"] for r in runs], dtype=float)
    q_calls = np.array([r["summary"]["eval_stats"]["slope_count"] for r in runs], dtype=float)
    seconds = np.array([r["summary"]["total_seconds"] for r in runs], dtype=float)
    success = np.array([r["summary"]["success"] for r in runs], dtype=bool)
    colors = np.where(success, "#2c7a7b", "#b84a62")
    sizes = 60 + 650 * seconds / max(float(seconds.max()), 1e-12)
    ax.scatter(active, q_calls, s=sizes, c=colors, edgecolor="#253140", linewidth=0.75, alpha=0.82, zorder=3)
    for idx, (x, y) in enumerate(zip(active, q_calls), start=1):
        ax.text(x, y, str(idx), ha="center", va="center", fontsize=7.1, color="white", fontweight="bold", zorder=4)
    ax.set_title("Support size vs slope workload", loc="left", fontsize=12.2, fontweight="bold", color="#17212b", pad=8)
    ax.set_xlabel("active coefficient terms", fontsize=8.8)
    ax.set_ylabel("Q slope calls", fontsize=8.8)
    ax.text(0.02, 0.92, "marker area = wall time", transform=ax.transAxes, fontsize=8.0, color="#52606d")
    legend_lines = [f"{i + 1}. {name}" for i, name in enumerate(names)]
    left = "\n".join(legend_lines[:7])
    right = "\n".join(legend_lines[7:])
    ax.text(0.02, 0.06, left, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.5, color="#344150",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="#fffaf0", edgecolor="#e0d4c0", alpha=0.92))
    ax.text(0.31, 0.06, right, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.5, color="#344150",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="#fffaf0", edgecolor="#e0d4c0", alpha=0.92))


def summary_footer(fig, data):
    s = data["summary"]
    runs = data["runs"]
    failed = [r["family"] for r in runs if not r["summary"]["success"]]
    footer = (
        f"ks(3,4), families=all   successes {s['successes']}/{s['runs']}   "
        f"roots {s['total_roots']}/{s['runs'] * COUNT}   "
        f"stress failure: {', '.join(failed) if failed else 'none'}   "
        f"flow: 115 exact telescopic slope reused by 116"
    )
    fig.text(
        0.5,
        0.032,
        footer,
        ha="center",
        va="center",
        fontsize=10.2,
        color="#17212b",
        bbox=dict(boxstyle="round,pad=0.42,rounding_size=0.12", facecolor="#fffaf0", edgecolor="#d2b48c", linewidth=1.0),
    )


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    engine = load_engine()
    data = run_benchmark(engine)
    runs = data["runs"]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "figure.facecolor": "#f8f4ec",
            "axes.facecolor": "#fffdf8",
            "axes.edgecolor": "#d2cabd",
            "axes.labelcolor": "#27323f",
            "xtick.color": "#52606d",
            "ytick.color": "#52606d",
            "savefig.facecolor": "#f8f4ec",
            "savefig.edgecolor": "#f8f4ec",
        }
    )
    fig = plt.figure(figsize=(18.0, 11.0), dpi=220)
    fig.suptitle("Engine 116: all-family multivariate pure Pandrosion", x=0.5, y=0.973, fontsize=22, fontweight="bold", color="#17212b")
    fig.text(0.5, 0.942, "flow/116 on ks(3,4): dense, sparse, structured, real, phase, mixed-degree, and ill-scaled systems", ha="center", fontsize=10.9, color="#52606d")

    gs = fig.add_gridspec(2, 3, left=0.045, right=0.985, top=0.895, bottom=0.09, hspace=0.32, wspace=0.25)
    architecture_panel(fig.add_subplot(gs[0, 0]))
    performance_panel(fig.add_subplot(gs[0, 1]), runs)
    root_cloud_panel(fig.add_subplot(gs[0, 2], projection="3d"), runs)
    simplex_panel(fig.add_subplot(gs[1, 0]), engine)
    residual_panel(fig.add_subplot(gs[1, 1]), runs)
    cost_panel(fig.add_subplot(gs[1, 2]), runs)
    summary_footer(fig, data)

    png = FIG_DIR / "116_multifamily_multivariate_all.png"
    pdf = FIG_DIR / "116_multifamily_multivariate_all.pdf"
    fig.savefig(png, dpi=260)
    fig.savefig(pdf)
    plt.close(fig)
    print(png)
    print(pdf)
    print(DATA_DIR / "116_multifamily_3x4_all_for_figure.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
