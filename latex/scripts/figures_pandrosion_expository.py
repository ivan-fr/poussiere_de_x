"""
Regenerate the three figures used in pandrosion_en_improved.tex:
    latex/figures/pandrosion_geometry.pdf
    latex/figures/pandrosion_figures.pdf
    latex/figures/pandrosion_steffensen.pdf

Dependencies: numpy, matplotlib.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def h_pandrosion(s: np.ndarray, x: float, p: int) -> np.ndarray:
    s = np.asarray(s, dtype=float)
    Sp = np.where(
        np.isclose(s, 1.0),
        float(p),
        (1.0 - s**p) / (1.0 - s),
    )
    return 1.0 - (x - 1.0) / (x * Sp)


def newton_step(u: float, x: float, p: int) -> float:
    return ((p - 1) * u + x / u ** (p - 1)) / p


def pandrosion_trajectory(s0: float, x: float, p: int, n: int) -> np.ndarray:
    s = np.empty(n + 1)
    s[0] = s0
    for k in range(n):
        s[k + 1] = h_pandrosion(s[k], x, p)
    return s


# ---------------------------------------------------------------------------
# Figure 1 — geometric construction (Thales diagram)
# ---------------------------------------------------------------------------

def make_geometry_figure(path: str) -> None:
    x = 2.0
    p = 3
    alpha = x ** (1.0 / p)
    s_star = 1.0 / alpha
    v_star = x * s_star ** (p - 1)

    # Show s0 clearly below s*. Only label the first three iterates since
    # s_2, s_3 collapse onto s*.
    s_vals = pandrosion_trajectory(0.5, x, p, 3)
    v_vals = x * s_vals ** (p - 1)

    fig, ax = plt.subplots(figsize=(7.4, 7.4))
    ax.set_facecolor("#fafafa")

    # Thales diagonal from (0, 1) to (x, 0)
    ts = np.linspace(0.0, 1.0, 400)
    ax.plot(x * (1.0 - ts), ts, color="black", lw=2.0,
            label="Diagonal (Thales)")

    # Curve v = x * s^{p-1}
    ss = np.linspace(0.0, 1.0, 400)
    ax.plot(x * ss ** (p - 1), ss, color="purple", lw=1.6,
            label=r"Curve $v = x\,s^{p-1}$")

    # Fixed-point markers
    ax.axhline(s_star, color="red", lw=1.2, ls="--",
               label="Fixed point $s^*$")
    ax.axvline(v_star, color="red", lw=1.0, ls=":", alpha=0.6)

    # Horizontal lines for each iterate
    for i, (s, v) in enumerate(zip(s_vals, v_vals)):
        ax.plot([0, v], [s, s], color="#3a78c2", lw=0.9, alpha=0.75)
        ax.plot([v, v], [0, s], color="#3a78c2", lw=0.7, alpha=0.3, ls=":")

    # Legend-visible points (add once)
    ax.plot(v_vals[0], s_vals[0], "o", color="#3a78c2", ms=5,
            label="Curve intersection (v, s)")
    ax.plot(x * (1.0 - s_vals[0]), s_vals[0], "s", color="#3a78c2",
            ms=4, alpha=0.7, label="Diagonal intersection")
    for s, v in zip(s_vals, v_vals):
        ax.plot(v, s, "o", color="#3a78c2", ms=5)
        ax.plot(x * (1.0 - s), s, "s", color="#3a78c2", ms=4, alpha=0.7)

    # Label s_n on the right-hand end of each horizontal line.
    # For the first three iterates the lines are visually separated; s_3
    # merges with s* so we do not label it (shown as merged point only).
    label_texts = [r"$s_0$", r"$s_1$", r"$s_2$"]
    for i in range(3):
        ax.annotate(
            label_texts[i],
            xy=(v_vals[i], s_vals[i]),
            xytext=(v_vals[i] + 0.05, s_vals[i] + 0.018),
            color="#3a78c2",
            fontsize=10,
            fontweight="bold",
            ha="left",
            va="center",
        )

    # Axis styling
    ax.set_xlim(-0.28, x + 0.05)
    ax.set_ylim(-0.22, 1.05)
    ax.set_xticks([0.0, 0.5, 1.0, v_star, 1.5, x])
    ax.set_xticklabels(["0", "0.5", "1.0", r"$\sqrt[3]{2}$", "1.5", "2.0"])
    ax.set_yticks([0.0, 0.5, s_star, 1.0])
    ax.set_yticklabels(["0", "0.5", r"$s^*$", "1"])
    ax.tick_params(axis="both", labelsize=8)

    # Red "s*" bracket on left margin
    ax.add_patch(Rectangle((-0.24, s_star - 0.004), 0.05, 0.008, color="red"))
    ax.text(-0.26, s_star + 0.055, r"$s^*$", color="red", fontsize=10,
            fontweight="bold")

    ax.set_title(
        "Pandrosion's construction  ($p=3$, $x=2$)\n"
        r"Successive parallels --- convergence to $v^* = \sqrt[3]{2}$",
        fontsize=10,
    )

    # Readout strip below the plot
    readout_parts = [
        f"$v_0 = {v_vals[0]:.3f}$",
        f"$v_1 = {v_vals[1]:.3f}$",
        f"$v_2 = {v_vals[2]:.3f}$",
        f"$v_3 = {v_vals[3]:.3f}$",
        f"$s^* = \\sqrt[3]{{2}} \\approx {v_star:.3f}$",
    ]
    ax.text(
        0.05, -0.14,
        "   ".join(readout_parts),
        fontsize=8.8, color="#3a78c2",
    )
    ax.text(
        0.05, -0.20,
        r"Output $v_n = x s_n^{p-1}$  (converges to $\sqrt[3]{2}$)",
        fontsize=8.8, color="black", alpha=0.75,
    )

    # Annotations on the diagonal and curve (positioned to avoid overlap
    # with iterate lines and the legend)
    ax.text(1.35, 0.40, "Diagonal = Thales", rotation=-33, fontsize=8.5,
            color="black", alpha=0.6)
    ax.text(1.78, 0.88, r"Curve $v = x s^{p-1}$",
            rotation=28, fontsize=8.5, color="purple", ha="right")

    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.92,
              bbox_to_anchor=(0.01, 0.99))
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 — cobweb + residual profiles
# ---------------------------------------------------------------------------

def make_cobweb_and_residuals_figure(path: str) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.5))

    # ----- Cobweb diagram (left) -----
    x, p = 2.0, 3
    alpha = x ** (1.0 / p)
    s_star = 1.0 / alpha

    ss = np.linspace(0.4, 1.0, 400)
    ax1.plot(ss, h_pandrosion(ss, x, p), color="#1f4ea1", lw=2.2,
             label=r"$s_{n+1} = h(s_n)$")
    ax1.plot(ss, ss, color="gray", lw=1.0, ls="--", label=r"$s_{n+1} = s_n$")
    ax1.axvline(s_star, color="red", lw=0.8, ls=":")
    ax1.plot(s_star, s_star, "o", color="red", ms=7,
             label=fr"$s^* = 2^{{-1/3}} \approx {s_star:.5f}$")

    # Cobweb steps from s0=0.5 (as in the caption)
    s = 0.5
    ax1.plot([s, s], [0.4, h_pandrosion(s, x, p)], color="#e07a1f", lw=1.0)
    for _ in range(6):
        s_next = h_pandrosion(s, x, p)
        ax1.plot([s, s_next], [s_next, s_next], color="#e07a1f", lw=1.0)
        ax1.plot([s_next, s_next], [s_next, h_pandrosion(s_next, x, p)],
                 color="#e07a1f", lw=1.0)
        s = s_next

    ax1.plot(0.5, 0.5, "s", color="#1b8a5a", ms=9, label=r"$s_0 = 0.5$")

    ax1.set_xlim(0.4, 1.0)
    ax1.set_ylim(0.4, 1.0)
    ax1.set_xlabel("$s_n$")
    ax1.set_ylabel("$s_{n+1}$")
    ax1.set_title("Cobweb diagram --- Pandrosion convergence\n"
                  "($p=3$, $x=2$, $s_0=0.5$)", fontsize=10)
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="upper left", fontsize=8, framealpha=0.9)

    # ----- Residual profiles (right) -----
    max_n = 18
    configs = [
        (3, 2, "#1f4ea1", "o", "Pandrosion $p=3$, $x=2$"),
        (2, 2, "#1b8a5a", "o", "Pandrosion $p=2$, $x=2$"),
        (3, 3, "#e07a1f", "o", "Pandrosion $p=3$, $x=3$"),
        (4, 2, "#7a2fb8", "o", "Pandrosion $p=4$, $x=2$"),
    ]
    for (p_, x_, color, marker, label) in configs:
        alpha_ = x_ ** (1.0 / p_)
        s_star_ = 1.0 / alpha_
        s0 = 0.5
        s = s0
        errs = []
        for n in range(max_n + 1):
            v = x_ * s ** (p_ - 1)
            errs.append(abs(v - alpha_))
            s = h_pandrosion(s, x_, p_)
        errs = np.array(errs)
        errs = np.where(errs < 1e-16, 1e-16, errs)
        ax2.semilogy(range(max_n + 1), errs, color=color, lw=1.4,
                     marker=marker, ms=4, label=label)

    # Newton for p=3, x=2 (doubly-exponential curvature)
    u = 1.5
    newton_errs = [abs(u - 2.0 ** (1 / 3))]
    for _ in range(8):
        u = newton_step(u, 2.0, 3)
        newton_errs.append(max(abs(u - 2.0 ** (1 / 3)), 1e-16))
    ax2.semilogy(range(len(newton_errs)), newton_errs,
                 color="#c72828", lw=1.4, ls="--",
                 marker="s", ms=5, label="Newton $p=3$, $x=2$")

    ax2.set_xlim(0, max_n)
    ax2.set_ylim(1e-13, 1)
    ax2.set_xlabel("Step $n$")
    ax2.set_ylabel(r"$\varepsilon_n = |v_n - \alpha|$")
    ax2.set_title("Residual profiles: Pandrosion vs Newton", fontsize=10)
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3 — Steffensen vs Newton
# ---------------------------------------------------------------------------

def make_steffensen_figure(path: str) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.5))

    x, p = 2.0, 3
    alpha = x ** (1.0 / p)
    s_star = 1.0 / alpha
    floor_tol = 1e-17

    # --- Plain Pandrosion (grey) ---
    s = 0.5
    pandro_err = [abs(x * s ** (p - 1) - alpha)]
    for _ in range(10):
        s = h_pandrosion(s, x, p)
        pandro_err.append(max(abs(x * s ** (p - 1) - alpha), floor_tol))

    # --- Newton ---
    u = 1.5
    newton_err = [abs(u - alpha)]
    for _ in range(6):
        u = newton_step(u, x, p)
        newton_err.append(max(abs(u - alpha), floor_tol))

    # --- Steffensen-Pandrosion ---
    def steffensen_step(s: float) -> float:
        h1 = h_pandrosion(s, x, p)
        h2 = h_pandrosion(h1, x, p)
        denom = h2 - 2.0 * h1 + s
        if abs(denom) < 1e-30:
            return h1
        return s - (h1 - s) ** 2 / denom

    s = 0.875
    steff_err = [abs(x * s ** (p - 1) - alpha)]
    for _ in range(4):
        s = steffensen_step(s)
        steff_err.append(max(abs(x * s ** (p - 1) - alpha), floor_tol))

    ax1.semilogy(range(len(pandro_err)), pandro_err, color="gray", lw=1.0,
                 marker=".", ms=5, alpha=0.6, label="Pandrosion (linear)")
    ax1.semilogy(range(len(newton_err)), newton_err, color="#c72828", lw=1.4,
                 marker="s", ms=6, ls="--", label="Newton")
    ax1.semilogy(range(len(steff_err)), steff_err, color="#1f4ea1", lw=1.8,
                 marker="D", ms=7, label="Steffensen--Pandrosion")

    ax1.axhline(1e-15, color="black", ls=":", lw=0.8, alpha=0.5)
    ax1.text(7.0, 1.3e-15, "Machine precision", fontsize=8, color="black",
             alpha=0.6)

    ax1.annotate(
        "$K_N \\approx 0.794$\n5 steps to $10^{-15}$",
        xy=(5, newton_err[5] if len(newton_err) > 5 else 1e-15),
        xytext=(5.8, 1e-6),
        fontsize=8,
        color="#c72828",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#c72828", alpha=0.9),
    )
    ax1.annotate(
        "$K_S \\approx 0.013$\n3 steps to $10^{-15}$",
        xy=(3, steff_err[3] if len(steff_err) > 3 else 1e-16),
        xytext=(3.3, 1e-10),
        fontsize=8,
        color="#1f4ea1",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#1f4ea1", alpha=0.9),
    )

    ax1.set_xlim(0, 10)
    ax1.set_ylim(1e-17, 1)
    ax1.set_xlabel("Step $n$")
    ax1.set_ylabel(r"$\varepsilon_n = |v_n - \sqrt[3]{2}|$")
    ax1.set_title(r"Convergence comparison ($p=3$, $x=2$)", fontsize=10)
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(loc="lower left", fontsize=9, framealpha=0.9)

    # --- Right panel: quadratic-constant bar chart ---
    K_S = 0.013
    K_N = 0.794
    labels = ["Steffensen--\nPandrosion\n2 evals of $S_p$ / step",
              "Newton\nRequires $f'(x)$\n1 eval / step"]
    values = [K_S, K_N]
    colors = ["#1f4ea1", "#c72828"]

    bars = ax2.bar(labels, values, color=colors, edgecolor="black", width=0.45)
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel("Quadratic constant $K$")
    ax2.set_title(r"Quadratic constant comparison"
                  "\n"
                  r"$\varepsilon_{n+1} \leq K \cdot \varepsilon_n^2$",
                  fontsize=10)
    ax2.grid(True, axis="y", alpha=0.25)

    for bar, val in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 0.02,
                 f"$K \\approx {val}$",
                 ha="center", fontsize=9, fontweight="bold",
                 color=bar.get_facecolor())

    arrow = FancyArrowPatch((0, K_S + 0.03), (0, K_N - 0.03),
                            arrowstyle="<->", mutation_scale=15,
                            color="black", lw=1.2)
    ax2.add_patch(arrow)
    ax2.text(0.08, (K_S + K_N) / 2, r"$\times 61$",
             fontsize=10, fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.3", fc="#fff4b8",
                       ec="black", alpha=0.9))

    ax2.tick_params(axis="x", labelsize=8)

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import os
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
    os.makedirs(outdir, exist_ok=True)
    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "savefig.dpi": 200,
    })
    make_geometry_figure(os.path.join(outdir, "pandrosion_geometry.pdf"))
    make_cobweb_and_residuals_figure(os.path.join(outdir, "pandrosion_figures.pdf"))
    make_steffensen_figure(os.path.join(outdir, "pandrosion_steffensen.pdf"))
    print("Wrote figures/pandrosion_geometry.pdf, figures/pandrosion_figures.pdf, "
          "figures/pandrosion_steffensen.pdf")
