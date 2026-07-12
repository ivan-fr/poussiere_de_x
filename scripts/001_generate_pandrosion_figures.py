from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon


OUT = Path("latex/figures/001")

INK = "#18212f"
MUTED = "#667085"
GRID = "#d8dee8"
BLUE = "#2463a6"
ORANGE = "#d97706"
GREEN = "#23856d"
RED = "#c43d4d"
PURPLE = "#7257a8"


def setup() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "axes.facecolor": "#ffffff",
            "figure.facecolor": "#ffffff",
            "savefig.facecolor": "#ffffff",
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "axes.linewidth": 0.8,
            "axes.titlesize": 13,
            "axes.titleweight": "regular",
            "legend.framealpha": 1.0,
            "legend.facecolor": "#ffffff",
            "legend.edgecolor": "#c8c8c8",
            "figure.dpi": 140,
            "savefig.dpi": 240,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.035)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", pad_inches=0.035, dpi=240)
    plt.close(fig)


def soften_axes(ax: plt.Axes) -> None:
    ax.grid(True, color=GRID, alpha=0.55, linewidth=0.55)
    for spine in ax.spines.values():
        spine.set_linewidth(0.75)
        spine.set_color(INK)


def sp(alpha: np.ndarray | float, p: int) -> np.ndarray | float:
    total = np.zeros_like(alpha, dtype=float)
    for k in range(p):
        total += alpha**k
    return total


def lambda_py(y: np.ndarray, p: int) -> np.ndarray:
    alpha = np.asarray(y, dtype=float) ** (1.0 / p)
    return 1.0 - p / sp(alpha, p)


def pandrosion_raw_errors(y: float, p: int, n: int = 40) -> tuple[np.ndarray, np.ndarray]:
    alpha = y ** (1.0 / p)

    def h(s: float) -> float:
        ssum = sum(s**k for k in range(p))
        return 1.0 - (y - 1.0) / (y * ssum)

    s = h(1.0)
    xs = []
    errs = []
    for k in range(n + 1):
        u = y * s ** (p - 1)
        xs.append(float(k))
        errs.append(abs(u - alpha) / alpha)
        s = h(s)
    return np.asarray(xs), np.asarray(errs)


def newton_errors(y: float, p: int, n: int = 9) -> tuple[np.ndarray, np.ndarray]:
    alpha = y ** (1.0 / p)
    s0 = 1.0 - (y - 1.0) / (y * sum(1.0**k for k in range(p)))
    u = y * s0 ** (p - 1)
    xs = []
    errs = []
    for k in range(n + 1):
        xs.append(float(k))
        errs.append(abs(u - alpha) / alpha)
        u = u - (u**p - y) / (p * u ** (p - 1))
    return np.asarray(xs), np.asarray(errs)


def select_scale(a: float, p: int, q: float = 1.25, k: int = 4) -> float:
    exponent = np.floor(k * np.log(a) / (p * np.log(q))) / k
    return float(q**exponent)


def irp_layer(u: float, y: float, p: int, q: float = 1.25, k: int = 4) -> float:
    residual = y / u**p
    a = residual if residual >= 1.0 else 1.0 / residual
    b = select_scale(a, p, q, k)
    normalized = a / b**p
    h1 = 1.0 - (normalized - 1.0) / (normalized * p)
    w = b * normalized * h1 ** (p - 1)
    return u * w if residual >= 1.0 else u / w


def irp2_errors(y: float, p: int, n: int = 4) -> tuple[np.ndarray, np.ndarray]:
    alpha = y ** (1.0 / p)
    u = 1.0
    xs: list[float] = []
    errs: list[float] = []
    for step in range(n + 1):
        xs.append(float(step))
        errs.append(abs(u - alpha) / alpha)
        u = irp_layer(irp_layer(u, y, p), y, p)
    return np.asarray(xs), np.asarray(errs)


def local_irp_log_error(ell: np.ndarray, p: int, layers: int) -> np.ndarray:
    current = np.asarray(ell, dtype=float)
    for _ in range(layers):
        a = np.exp(p * np.abs(current))
        h1 = 1.0 - (1.0 - 1.0 / a) / p
        logw = np.log(a) + (p - 1) * np.log(h1)
        current = current - np.sign(current) * logw
    return current


def fig1_pandrosion_construction() -> None:
    fig, ax = plt.subplots(figsize=(7.8, 3.45))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-0.18, 4.1)
    ax.set_ylim(-0.40, 2.55)

    apex = np.array([0.0, 2.15])
    base = np.array([0.0, 1.18, 2.45, 3.82])

    def top_y(x: float) -> float:
        return apex[1] * (1.0 - x / base[-1])

    fills = ["#dceafb", "#fff0d9", "#dcefe9"]
    for left, right, color in zip(base[:-1], base[1:], fills):
        ax.add_patch(
            Polygon(
                [(left, 0.0), (left, top_y(left)), (right, 0.0)],
                closed=True,
                facecolor=color,
                edgecolor="none",
                alpha=0.95,
                zorder=1,
            )
        )

    ax.plot([0, base[-1]], [0, 0], color=INK, lw=1.35, zorder=3)
    ax.plot([0, 0], [0, apex[1]], color=INK, lw=1.35, zorder=3)
    ax.plot([0, base[-1]], [apex[1], 0], color=INK, lw=1.2, zorder=4)
    for x in base[1:-1]:
        ax.plot([0, x], [apex[1], 0], color=INK, lw=0.95, alpha=0.9, zorder=4)

    verticals = [(base[1], ORANGE), (base[2], GREEN), (base[3], RED)]
    for x, color in verticals:
        ax.plot([x, x], [0, top_y(x)], color=color, lw=1.45, ls=(0, (4, 2.5)), zorder=5)
        ax.plot([x - 0.025, x + 0.025], [top_y(x), top_y(x)], color=color, lw=1.0, zorder=5)

    for x in base:
        ax.plot([x, x], [-0.055, 0.055], color=INK, lw=0.8)

    ax.text(-0.05, -0.16, r"$O$", ha="right", va="top", fontsize=9)
    ax.text(0.055, -0.16, r"$1$", ha="left", va="top", fontsize=10.5)
    for x, label in zip(base[1:], [r"$m_1$", r"$m_2$", r"$x$"]):
        ax.text(x, -0.20, label, ha="center", va="top", fontsize=10.5)

    ax.text(1.91, 2.36, r"$1:m_1=m_1:m_2=m_2:x$", ha="center", fontsize=15)
    ax.text(1.91, 2.14, "successive similarities encode the proportional means", ha="center", fontsize=10, color=MUTED)
    save(fig, "fig1_pandrosion_construction")


def fig2_lambda_profile() -> None:
    t = np.linspace(0.0, 5.0, 650)
    y = np.exp(t)
    fig, ax = plt.subplots(figsize=(7.8, 3.75))
    for p, color in [(2, BLUE), (3, ORANGE), (5, GREEN), (8, PURPLE)]:
        ax.plot(t, lambda_py(y, p), lw=1.75, color=color, label=fr"$p={p}$")
    local_t = np.linspace(0.0, 0.9, 120)
    ax.plot(local_t, local_t / 3.0, color=RED, lw=1.15, ls="--", label=r"$\delta/3$ tangent for $p=3$")
    ax.set_xlabel(r"normalized logarithmic distance $\log y$")
    ax.set_ylabel(r"raw multiplier $\lambda_{p,y}$")
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1.02)
    ax.set_title("Residual fingerprint of raw Pandrosion", pad=7)
    soften_axes(ax)
    ax.legend(loc="lower right", fontsize=8.5, ncol=2)
    save(fig, "fig2_lambda_profile")


def lattice_gap(phase: np.ndarray, k: int) -> np.ndarray:
    grid = np.linspace(0.0, 1.0, k + 1)
    distance = np.full_like(phase, np.inf, dtype=float)
    for point in grid:
        distance = np.minimum(distance, np.abs(phase - point))
    return distance


def fig3_scale_lattice() -> None:
    phase = np.linspace(0.0, 1.0, 800)
    fig, axes = plt.subplots(3, 1, figsize=(8.25, 4.85), sharex=True)
    for ax, k in zip(axes, [1, 2, 4]):
        gap = lattice_gap(phase, k)
        bound = 1.0 / (2.0 * k)
        ax.fill_between(phase, 0, gap, color="#d6e6f6", alpha=0.9)
        ax.plot(phase, gap, color=BLUE, lw=1.55)
        ax.axhline(bound, color=RED, lw=0.9, ls="--")
        for point in np.linspace(0.0, 1.0, k + 1):
            ax.axvline(point, color=INK, lw=0.65, alpha=0.65)
        ax.text(0.07, 0.72 * bound + 0.03, fr"$K={k}$", fontsize=11)
        ax.text(0.88, bound + 0.018, r"$h/(2K)$", ha="right", fontsize=9)
        ax.set_ylabel("gap")
        ax.set_ylim(0, 0.55)
        soften_axes(ax)
        ax.grid(True, axis="y", color=GRID, alpha=0.55, linewidth=0.55)
    axes[0].set_title("Scale-lattice geometry: equally spaced phases minimize the worst gap", pad=7)
    axes[-1].set_xlabel(r"logarithmic target phase $(L\ \mathrm{mod}\ h)/h$, with $h=p\log q$")
    save(fig, "fig3_scale_lattice")


def fig4_convergence_comparison() -> None:
    y = 4.0
    p = 3
    xr, er = pandrosion_raw_errors(y, p, 34)
    xn, en = newton_errors(y, p, 8)
    xi, ei = irp2_errors(y, p, 4)
    floor = 1e-15

    fig, ax = plt.subplots(figsize=(7.9, 4.0))
    ax.plot(xr, np.log10(np.maximum(er, floor)), "o-", color=BLUE, lw=1.45, ms=2.8, label="raw Pandrosion")
    ax.plot(xi, np.log10(np.maximum(ei, floor)), "s-", color=ORANGE, lw=1.45, ms=3.2, label="IRP-2")
    ax.plot(xn, np.log10(np.maximum(en, floor)), "^-", color=GREEN, lw=1.45, ms=3.0, label="Newton")
    ax.axhline(np.log10(floor), color=INK, lw=0.85, ls=":")
    ax.text(0.6, -14.55, "display floor", color=MUTED, fontsize=8.5)
    ax.set_xlabel("iteration / macro-update")
    ax.set_ylabel(r"$\log_{10}$ relative error")
    ax.set_xlim(-0.2, 34.5)
    ax.set_ylim(-15.8, 0.4)
    ax.set_title(r"Convergence on the normalized problem $u^3=4$", pad=7)
    soften_axes(ax)
    ax.legend(loc="upper right", fontsize=8.3)
    save(fig, "fig4_convergence_comparison")


def fig5_global_convergence() -> None:
    p = 3
    x = 8.0
    fixed = x ** (-1.0 / p)

    def h(s: np.ndarray | float) -> np.ndarray | float:
        return 1.0 - (x - 1.0) / (x * sp(np.asarray(s), p))

    grid = np.linspace(0.04, 1.12, 500)
    fig, ax = plt.subplots(figsize=(7.4, 4.25))
    ax.plot(grid, grid, color=MUTED, lw=1.15, ls="--", label=r"diagonal $s_{n+1}=s_n$")
    ax.plot(grid, h(grid), color=BLUE, lw=2.0, label=r"$h_{3,8}(s)$")

    for start, color in [(0.12, ORANGE), (1.02, GREEN)]:
        value = start
        for _ in range(5):
            nxt = float(h(value))
            ax.annotate("", xy=(value, nxt), xytext=(value, value),
                        arrowprops=dict(arrowstyle="-", color=color, lw=1.0, alpha=0.9))
            ax.annotate("", xy=(nxt, nxt), xytext=(value, nxt),
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.0, alpha=0.9))
            value = nxt

    ax.scatter([fixed], [fixed], s=46, color=RED, zorder=5)
    ax.annotate(r"$s_*=8^{-1/3}=1/2$", xy=(fixed, fixed), xytext=(0.62, 0.34),
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.9), fontsize=10.5, color=RED)
    ax.set_xlabel(r"current internal coordinate $s_n$")
    ax.set_ylabel(r"next coordinate $s_{n+1}$")
    ax.set_xlim(0.04, 1.12)
    ax.set_ylim(0.04, 1.12)
    ax.set_title("Global monotone convergence on the positive real line", pad=7)
    soften_axes(ax)
    ax.legend(loc="upper left", fontsize=9)
    save(fig, "fig5_global_convergence")


def fig6_riemann_inversion() -> None:
    fig, ax = plt.subplots(figsize=(7.6, 2.65))
    ax.axis("off")
    ax.set_xlim(-3.1, 3.1)
    ax.set_ylim(-1.0, 1.35)

    logx = 0.82
    ticks = np.linspace(-2.4, 2.4, 7)
    for tick in ticks:
        ax.plot([tick, tick], [-0.56, 0.92], color="#d1d1d1", lw=0.75, zorder=1)

    ax.plot([-2.65, 2.65], [0.0, 0.0], color=INK, lw=1.15, zorder=3)
    ax.plot([-logx, -logx], [-0.08, 0.08], color=INK, lw=1.0)
    ax.plot([logx, logx], [-0.08, 0.08], color=INK, lw=1.0)
    ax.plot([0, 0], [-0.14, 0.14], color=INK, lw=0.9)
    ax.scatter([-logx, logx], [0, 0], s=28, color=INK, zorder=4)
    ax.text(-logx, -0.22, r"$-\log x$", ha="center", va="top", fontsize=9)
    ax.text(logx, -0.22, r"$\log x$", ha="center", va="top", fontsize=9)
    ax.text(0, -0.22, r"$0$", ha="center", va="top", fontsize=9)

    for tick in [-2.0, -1.2, -0.35]:
        ax.plot([tick, tick], [0.18, 0.44], color=BLUE, lw=1.55)
        ax.text(tick, 0.50, r"$p\log B$", ha="center", fontsize=7.4, color=MUTED)
    for tick in [0.35, 1.2, 2.0]:
        ax.plot([tick, tick], [0.18, 0.44], color=ORANGE, lw=1.55)
        ax.text(tick, 0.50, r"$p\log B$", ha="center", fontsize=7.4, color=MUTED)

    ax.text(-1.85, 0.72, r"reciprocal residual $\log y_\infty(B)$", color=BLUE, ha="center", fontsize=9.2)
    ax.text(1.85, 0.72, r"direct residual $\log y_+(B)$", color=ORANGE, ha="center", fontsize=9.2)
    box = dict(boxstyle="square,pad=0.18", facecolor="#ffffff", edgecolor="#b8b8b8", linewidth=0.75)
    ax.text(-1.52, -0.58, r"reciprocal chart: $\log y_\infty(B)=-\log x-p\log B$", ha="center", fontsize=8.7, bbox=box)
    ax.text(1.52, -0.58, r"direct chart: $\log y_+(B)=\log x-p\log B$", ha="center", fontsize=8.7, bbox=box)
    ax.text(0, -0.95, "choose the chart and scale whose normalized log-target is closest to 0", ha="center", fontsize=9.4)
    save(fig, "fig6_riemann_inversion")


def fig7_irp_order() -> None:
    p = 3
    # Keep the smallest inputs above the double-precision floor for the
    # eighth-order curve; the slopes remain visible over more than a decade.
    ell = np.logspace(-2.0, -0.65, 260)
    fig, ax = plt.subplots(figsize=(7.4, 4.25))
    for layers, color, marker in [(1, BLUE, "o"), (2, ORANGE, "s"), (3, GREEN, "^")]:
        out = np.abs(local_irp_log_error(ell, p, layers))
        ax.loglog(ell, out, color=color, lw=1.8,
                  marker=marker, markevery=44, ms=3.2,
                  label=fr"{layers} layer{'s' if layers > 1 else ''}: order $2^{layers}$")
    ax.set_xlabel(r"input logarithmic error $|\ell|$")
    ax.set_ylabel(r"output logarithmic error $|\Phi_3^{\circ k}(\ell)|$")
    ax.set_title("Dyadic order amplification of local IRP", pad=7)
    soften_axes(ax)
    ax.legend(loc="lower right", fontsize=9)
    save(fig, "fig7_irp_order")


def main() -> None:
    setup()
    fig1_pandrosion_construction()
    fig2_lambda_profile()
    fig3_scale_lattice()
    fig4_convergence_comparison()
    fig5_global_convergence()
    fig6_riemann_inversion()
    fig7_irp_order()
    for path in sorted(OUT.glob("fig[1-7]_*")):
        print(path)


if __name__ == "__main__":
    main()
