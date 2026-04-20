"""
Generate a ten-panel proof gallery for latex/pandrosion_paper.tex.

The figures are conceptual/numerical illustrations of Lean-certified theorem
families. They are not additional proofs; the captions point to the formal
Lean modules that certify the displayed identities.
"""

from pathlib import Path
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "latex"
DPI = 300

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "figure.titlesize": 14,
    }
)

COLORS = {
    "blue": "#26547C",
    "red": "#C44536",
    "green": "#2A9D8F",
    "gold": "#E9A23B",
    "purple": "#6D597A",
    "ink": "#1F2933",
    "gray": "#6B7280",
}


def save(fig, name):
    path = OUT / name
    fig.tight_layout()
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path.relative_to(ROOT)}")


def pandrosion_step(A, B, X=2):
    return A * (A**3 + 4 * X * B**3), B * (3 * A**3 + 2 * X * B**3)


def phi(A, B, X=2):
    return A**9 - 14 * X * A**6 * B**3 - 20 * X**2 * A**3 * B**6 + 8 * X**3 * B**9


def fig_universal_rate():
    X = np.linspace(0.25, 10, 300)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.plot(X, np.full_like(X, -0.2), color=COLORS["blue"], lw=3, label="Pandrosion: P'(r) = -1/5")
    ax.plot(X, np.zeros_like(X), color=COLORS["red"], lw=2.4, ls="--", label="Newton/Halley local derivative = 0")
    ax.fill_between(X, -1, 1, color=COLORS["green"], alpha=0.08, label="contractive band |lambda| < 1")
    ax.axhline(1, color=COLORS["gray"], lw=1)
    ax.axhline(-1, color=COLORS["gray"], lw=1)
    ax.scatter([2], [-0.2], s=70, color=COLORS["gold"], zorder=4)
    ax.annotate("root-independent\nalternating rate", xy=(2, -0.2), xytext=(3.0, -0.55),
                arrowprops=dict(arrowstyle="->", color=COLORS["ink"]), color=COLORS["ink"])
    ax.set_xlabel("target X")
    ax.set_ylabel("local multiplier at a cube root")
    ax.set_title("Universal derivative certificate")
    ax.set_ylim(-1.05, 0.35)
    ax.grid(alpha=0.2)
    ax.legend(loc="lower right")
    save(fig, "fig_proof_gallery_01_universal_rate.pdf")


def fig_chebyshev_halley_exclusion():
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.2, 1.1)
    ax.hlines(0.35, 0, 1, color=COLORS["ink"], lw=2)
    for x, label, color in [
        (0, "Chebyshev\nalpha=0", COLORS["green"]),
        (0.4, "forced by\ns^3 X: 2/5", COLORS["red"]),
        (0.5, "forced by\ns^6: 1/2", COLORS["blue"]),
        (1, "super-Halley\nalpha=1", COLORS["purple"]),
    ]:
        ax.vlines(x, 0.24, 0.46, color=color, lw=3)
        ax.scatter([x], [0.35], s=80, color=color, zorder=3)
        ax.text(x, 0.55 if x != 0.4 else 0.02, label, ha="center", va="center", color=color)
    ax.annotate(
        "one parameter cannot be\nboth 2/5 and 1/2",
        xy=(0.45, 0.35),
        xytext=(0.62, 0.9),
        arrowprops=dict(arrowstyle="->", color=COLORS["ink"], lw=1.3),
        color=COLORS["ink"],
    )
    ax.text(0.06, 0.82, "6(3 - 2 alpha) = 15 - 6 alpha  ->  alpha = 1/2", color=COLORS["blue"])
    ax.text(0.06, 0.72, "12 alpha = 4 + 2 alpha          ->  alpha = 2/5", color=COLORS["red"])
    ax.text(0.06, -0.08, "Lean conclusion: no alpha satisfies both matching equations.", color=COLORS["ink"])
    ax.set_axis_off()
    ax.set_title("Chebyshev-Halley exclusion by incompatible coefficients")
    save(fig, "fig_proof_gallery_02_ch_exclusion.pdf")


def fig_residual_conservation():
    X = 2.0
    root = X ** (1 / 3)
    s = np.linspace(0.35, 2.4, 600)
    den = (3 * s**3 + 2 * X) ** 3
    factor = s**9 - 14 * s**6 * X - 20 * s**3 * X**2 + 8 * X**3
    multiplier = factor / den
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.plot(s, multiplier, color=COLORS["blue"], lw=2.6)
    ax.axhline(0, color=COLORS["gray"], lw=1)
    ax.axhline(-0.2, color=COLORS["red"], lw=1.8, ls="--", label="limit at root: -1/5")
    ax.axvline(root, color=COLORS["gold"], lw=1.8, ls=":")
    ax.scatter([root], [-0.2], s=80, color=COLORS["gold"], zorder=4)
    ax.annotate("residual multiplier\n(G^3-X)/(s^3-X)", xy=(root, -0.2), xytext=(1.55, 0.15),
                arrowprops=dict(arrowstyle="->", color=COLORS["ink"]), color=COLORS["ink"])
    ax.set_xlabel("state s")
    ax.set_ylabel("exact residual multiplier")
    ax.set_title("Kinematic residual conservation factor for X=2")
    ax.set_ylim(-1.1, 0.55)
    ax.grid(alpha=0.2)
    ax.legend()
    save(fig, "fig_proof_gallery_03_residual_conservation.pdf")


def fig_no_cycles():
    n = np.arange(0, 13)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    for c, color in [(0.2, COLORS["blue"]), (0.5, COLORS["green"]), (0.8, COLORS["gold"])]:
        ax.semilogy(n, c**n, "o-", color=color, lw=2, label=f"c = {c}")
    ax.text(5.9, 0.22, "periodic return would require\n|s-r| <= c^n |s-r| < |s-r|",
            color=COLORS["ink"], bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COLORS["gray"], alpha=0.9))
    ax.set_xlabel("iterate n")
    ax.set_ylabel("upper bound c^n |s-r| / |s-r|")
    ax.set_title("Contraction iterates exclude finite cycles")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    save(fig, "fig_proof_gallery_04_no_cycles.pdf")


def fig_integer_amplification():
    X = 2
    A, B = 1, 1
    d_logs = []
    phi_logs = []
    labels = []
    for k in range(4):
        d = A**3 - X * B**3
        p = phi(A, B, X)
        d_logs.append(math.log10(abs(d)))
        phi_logs.append(math.log10(abs(p)))
        labels.append(str(k))
        A, B = pandrosion_step(A, B, X)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    x = np.arange(len(labels))
    ax.bar(x - 0.18, d_logs, width=0.36, color=COLORS["blue"], label="log10 |d_n|")
    ax.bar(x + 0.18, phi_logs, width=0.36, color=COLORS["red"], label="log10 |Phi_n|")
    ax.set_xticks(x, labels)
    ax.set_xlabel("integer step n")
    ax.set_ylabel("decimal logarithm")
    ax.set_title("Pell-Pandrosion norm amplification for X=2")
    ax.grid(axis="y", alpha=0.2)
    ax.legend()
    save(fig, "fig_proof_gallery_05_integer_amplification.pdf")


def fig_cross_determinant():
    X = 2
    A, B = 1, 1
    det_abs = []
    d_abs = []
    for _ in range(4):
        Ap, Bp = pandrosion_step(A, B, X)
        det_abs.append(abs(A * Bp - B * Ap))
        d_abs.append(abs(A**3 - X * B**3))
        A, B = Ap, Bp
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    y1 = [math.log10(v) for v in det_abs]
    y2 = [math.log10(v) for v in d_abs]
    ax.plot(range(4), y1, "o-", lw=2.4, color=COLORS["purple"], label="log10 |A B' - B A'|")
    ax.plot(range(4), y2, "s--", lw=2.0, color=COLORS["green"], label="log10 |A^3 - X B^3|")
    ax.text(0.35, max(y1) * 0.55, "A B' - B A' = 2 A B (A^3 - X B^3)",
            color=COLORS["ink"], bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COLORS["gray"], alpha=0.9))
    ax.set_xlabel("step n")
    ax.set_ylabel("decimal logarithm")
    ax.set_title("Cross-determinant separation of consecutive approximants")
    ax.grid(alpha=0.2)
    ax.legend(loc="upper left")
    save(fig, "fig_proof_gallery_06_cross_determinant.pdf")


def fig_voronoi_stability():
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    r = np.array([-1.0, 0.0])
    rp = np.array([1.0, 0.0])
    yy = np.linspace(-1.8, 1.8, 200)
    ax.fill_betweenx(yy, -2.2, 0, color=COLORS["blue"], alpha=0.10)
    ax.fill_betweenx(yy, 0, 2.2, color=COLORS["red"], alpha=0.08)
    ax.axvline(0, color=COLORS["gray"], lw=1.8, ls="--", label="Voronoi bisector")
    circle = plt.Circle(tuple(r), 0.55, color=COLORS["blue"], fill=False, lw=2)
    ax.add_patch(circle)
    z = np.array([-1.38, 0.32])
    fz = r + 0.45 * (z - r)
    ax.scatter([r[0], rp[0]], [r[1], rp[1]], s=100, color=[COLORS["blue"], COLORS["red"]], zorder=5)
    ax.text(r[0] - 0.08, r[1] - 0.22, "r", ha="right", color=COLORS["blue"])
    ax.text(rp[0] + 0.08, rp[1] - 0.22, "r'", ha="left", color=COLORS["red"])
    ax.scatter([z[0], fz[0]], [z[1], fz[1]], s=75, color=[COLORS["gold"], COLORS["green"]], zorder=6)
    ax.annotate("", xy=fz, xytext=z, arrowprops=dict(arrowstyle="->", lw=2, color=COLORS["ink"]))
    ax.text(z[0] - 0.05, z[1] + 0.15, "z", color=COLORS["gold"])
    ax.text(fz[0] + 0.08, fz[1] + 0.07, "f(z)", color=COLORS["green"])
    ax.text(-1.85, 1.45, "(1+2c)|z-r| < |z-r'|", color=COLORS["ink"],
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COLORS["gray"], alpha=0.9))
    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(-1.8, 1.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Voronoi basin stability under contraction")
    ax.set_xlabel("real coordinate")
    ax.set_ylabel("imaginary coordinate")
    ax.grid(alpha=0.15)
    save(fig, "fig_proof_gallery_07_voronoi_stability.pdf")


def fig_hermitian_preservation():
    rng = np.random.default_rng(4)
    real_eigs = np.sort(rng.normal(0, 1, 42))
    generic = rng.normal(0, 1, 42) + 1j * rng.normal(0, 0.65, 42)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.4), sharey=True)
    ax1.scatter(generic.real, generic.imag, s=25, color=COLORS["red"], alpha=0.75)
    ax1.axhline(0, color=COLORS["gray"], lw=1)
    ax1.set_title("generic product")
    ax1.set_xlabel("Re(lambda)")
    ax1.set_ylabel("Im(lambda)")
    ax1.grid(alpha=0.2)
    ax2.scatter(real_eigs, np.zeros_like(real_eigs), s=28, color=COLORS["blue"], alpha=0.85)
    ax2.axhline(0, color=COLORS["gray"], lw=1)
    ax2.set_title("Hermitian Pandrosion product")
    ax2.set_xlabel("Re(lambda)")
    ax2.grid(alpha=0.2)
    fig.suptitle("Hermitian preservation confines the spectrum to the real axis")
    save(fig, "fig_proof_gallery_08_hermitian_preservation.pdf")


def fig_dft_spectral():
    p = 12
    theta = np.linspace(0, 2 * np.pi, 500)
    roots = np.exp(2j * np.pi * np.arange(p) / p)
    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    ax.plot(np.cos(theta), np.sin(theta), color=COLORS["gray"], lw=1)
    ax.scatter(roots.real, roots.imag, s=80, color=COLORS["blue"], zorder=4)
    for z in roots:
        ax.plot([0, z.real], [0, z.imag], color=COLORS["blue"], alpha=0.18)
    ax.arrow(0, 0, roots[1].real * 0.8, roots[1].imag * 0.8,
             width=0.01, head_width=0.07, color=COLORS["gold"], length_includes_head=True)
    ax.text(-0.72, -1.18, "sum_j omega^(jm) = 0 unless p divides m", color=COLORS["ink"])
    ax.text(-0.63, 1.08, "character orthogonality", color=COLORS["purple"])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_title("DFT decomposition: roots of unity cancel nonzero modes")
    ax.grid(alpha=0.15)
    save(fig, "fig_proof_gallery_09_dft_spectral.pdf")


def fig_effective_irrationality():
    r = 2 ** (1 / 3)
    B = np.arange(1, 90)
    lower = 1 / (3 * r**2 * B**3)
    thue = 43 / (3 * r**2 * B**3)
    random_like = 1 / B**2.25
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.loglog(B, lower, color=COLORS["blue"], lw=2.5, label="Liouville-type lower envelope")
    ax.loglog(B, thue, color=COLORS["red"], lw=2.0, ls="--", label="first amplified X=2 envelope")
    ax.loglog(B, random_like, color=COLORS["green"], lw=1.8, alpha=0.75, label="reference B^-2.25")
    ax.text(10, thue[9] * 1.8, "|A-rB| >= 1/(A^2 + A rB + (rB)^2)",
            color=COLORS["ink"], bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COLORS["gray"], alpha=0.9))
    ax.set_xlabel("denominator scale B")
    ax.set_ylabel("certified lower distance scale")
    ax.set_title("Effective irrationality bound from cubic factorization")
    ax.grid(alpha=0.22, which="both")
    ax.legend()
    save(fig, "fig_proof_gallery_10_effective_irrationality.pdf")


def main():
    fig_universal_rate()
    fig_chebyshev_halley_exclusion()
    fig_residual_conservation()
    fig_no_cycles()
    fig_integer_amplification()
    fig_cross_determinant()
    fig_voronoi_stability()
    fig_hermitian_preservation()
    fig_dft_spectral()
    fig_effective_irrationality()


if __name__ == "__main__":
    main()
