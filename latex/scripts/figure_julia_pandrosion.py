"""
figure_julia_pandrosion.py

Pandrosion-based analogue of Figure 13 of `p_previous_version.pdf`
(archived as `latex/archive/legacy-pandrosion/pandrosion_paper.pdf`).

Figure 13 of that paper shows the classical Julia set  z -> z^2 + c  with
c = -0.7 + 0.27 i, used there as a *contrast*: a fractal boundary that
Pandrosion avoids (Figures 10-12, 14; Theorems 5.1-5.5).  This script
reproduces the *Pandrosion* counterpart on the same complex parameter
value by using it as the target X in the cubic Pandrosion map

    P_X(z)  =  z * (z^3 + 4 X) / (3 z^3 + 2 X)

whose fixed points are characterised by

    P_X(z) = z  iff  z = 0  or  z^3 = X           (Thm 5.1, Deep.lean)

and whose universal contraction rate at every cube root is

    P'_X(X^{1/3}) = -1/5                          (Thm 2.4, HalleyComparison.lean).

For each pixel z0 in [-1.5, 1.5]^2 we iterate P_X and record:
  * which cube root of X the orbit converges to  (three-basin
    red / green / blue colouring, Figures 10-11 of the paper),
  * how fast it converges  (iteration-count shading, Figure 14 of the
    paper, escape-time aesthetic inherited from classical Julia pictures).

Overlaid: the Voronoi bisectors of the three cube roots (Thm 5.3,
MultiStart.lean).  These affine bisectors are the *Pandrosion* Julia /
basin frontier analogue; per Part I, Section 6.1, for the adaptive
scheme the boundaries empirically vanish altogether.  See also
Paper 0 Section 17 (Lebesgue-null Voronoi frontier on h_{p,x}).

Output: latex/figures/julia_pandrosion_hd.png (HD, 2400 x 2400).
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path


OUTPUT = Path(__file__).resolve().parents[1] / "figures" / "julia_pandrosion_hd.png"


# ------------------------------------------------------------------
# Pandrosion cubic map P_X (eq. 2 of the archived pandrosion_paper.pdf)
# ------------------------------------------------------------------
def P_X(z: np.ndarray, X: complex) -> np.ndarray:
    z3 = z * z * z
    num = z * (z3 + 4.0 * X)
    den = 3.0 * z3 + 2.0 * X
    return num / den


# ------------------------------------------------------------------
# Basin + iteration field
# ------------------------------------------------------------------
def pandrosion_basins(
    X: complex,
    res: int,
    box: float,
    max_iter: int,
    tol: float,
    bailout: float,
):
    axis = np.linspace(-box, box, res)
    Re, Im = np.meshgrid(axis, axis)
    Z = (Re + 1j * Im).astype(np.complex128)

    alpha = X ** (1.0 / 3.0)
    omega = np.exp(2j * np.pi / 3.0)
    roots = np.array([alpha, alpha * omega, alpha * omega * omega])

    basin = np.full(Z.shape, -1, dtype=np.int8)
    iters = np.full(Z.shape, max_iter, dtype=np.int32)

    for k in range(max_iter):
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            Zn = P_X(Z, X)
        for idx, r in enumerate(roots):
            fresh = (np.abs(Zn - r) < tol) & (basin < 0)
            basin[fresh] = idx
            iters[fresh] = k
        Zn = np.where(np.isfinite(Zn) & (np.abs(Zn) < bailout), Zn, 0.0 + 0j)
        Z = Zn

    return basin, iters, roots


# ------------------------------------------------------------------
# Voronoi bisectors (Thm 5.3, MultiStart.lean)
# ------------------------------------------------------------------
def voronoi_bisectors(roots: np.ndarray, box: float):
    segs = []
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            a, b = roots[i], roots[j]
            m = 0.5 * (a + b)
            d = b - a
            perp = 1j * d / abs(d)
            ts = []
            if perp.real != 0:
                for xb in (-box, box):
                    t = (xb - m.real) / perp.real
                    if -box - 1e-9 <= m.imag + t * perp.imag <= box + 1e-9:
                        ts.append(t)
            if perp.imag != 0:
                for yb in (-box, box):
                    t = (yb - m.imag) / perp.imag
                    if -box - 1e-9 <= m.real + t * perp.real <= box + 1e-9:
                        ts.append(t)
            if len(ts) >= 2:
                t1, t2 = min(ts), max(ts)
                p1 = m + t1 * perp
                p2 = m + t2 * perp
                segs.append(([p1.real, p2.real], [p1.imag, p2.imag]))
    return segs


# ------------------------------------------------------------------
# RGB rendering: basin index picks hue, iteration count picks value
# ------------------------------------------------------------------
def basin_rgb(basin: np.ndarray, iters: np.ndarray, max_iter: int) -> np.ndarray:
    palette = np.array(
        [
            [0.95, 0.30, 0.25],  # red   -- root alpha
            [0.25, 0.45, 0.96],  # blue  -- root alpha * omega
            [0.20, 0.80, 0.40],  # green -- root alpha * omega^2
        ],
        dtype=np.float64,
    )
    shade = 1.0 - 0.80 * (iters.astype(np.float64) / max(max_iter - 1, 1))
    shade = np.clip(shade, 0.18, 1.0)
    rgb = np.zeros((*basin.shape, 3), dtype=np.float64)
    for idx in range(3):
        mask = basin == idx
        rgb[mask] = palette[idx] * shade[mask][:, None]
    return np.clip(rgb, 0.0, 1.0)


# ------------------------------------------------------------------
# Render
# ------------------------------------------------------------------
def render(
    X: complex = -0.7 + 0.27j,
    res: int = 2400,
    box: float = 1.5,
    max_iter: int = 60,
    tol: float = 1.0e-8,
    bailout: float = 1.0e6,
    output: str = OUTPUT,
) -> None:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    basin, iters, roots = pandrosion_basins(
        X, res, box, max_iter, tol, bailout
    )
    print(
        "  iterations: min={}, max={}, mean={:.1f} | "
        "basins = ({}, {}, {}), unresolved = {}".format(
            iters.min(),
            iters.max(),
            iters.mean(),
            int((basin == 0).sum()),
            int((basin == 1).sum()),
            int((basin == 2).sum()),
            int((basin < 0).sum()),
        )
    )

    rgb = basin_rgb(basin, iters, max_iter)

    fig, ax = plt.subplots(figsize=(9.3, 8.4), dpi=240)
    ax.imshow(
        rgb,
        origin="lower",
        extent=[-box, box, -box, box],
        interpolation="bilinear",
    )

    for xs, ys in voronoi_bisectors(roots, box):
        ax.plot(xs, ys, "-", color="cyan", linewidth=0.75, alpha=0.85)

    ax.plot(
        roots.real,
        roots.imag,
        marker="*",
        linestyle="",
        color="white",
        markeredgecolor="black",
        markeredgewidth=1.2,
        markersize=15,
    )

    ax.set_xlim(-box, box)
    ax.set_ylim(-box, box)
    ax.set_xlabel(r"$\Re(z)$")
    ax.set_ylabel(r"$\Im(z)$")

    X_txt = f"{X.real:g}{X.imag:+g}i"
    ax.set_title(
        r"Pandrosion Julia set: "
        r"$\mathcal{P}_X(z) = z(z^3 + 4X)/(3 z^3 + 2X)$, "
        rf"$X = {X_txt}$"
        "\n"
        r"Three cube-root basins, shaded by iterations to convergence "
        r"($\mathcal{P}'_X(\!\sqrt[3]{X}) = -1/5$, Thm 2.4)"
    )

    handles = [
        mpatches.Patch(color=(0.95, 0.30, 0.25), label=r"basin of $\sqrt[3]{X}$"),
        mpatches.Patch(color=(0.25, 0.45, 0.96), label=r"basin of $\omega\sqrt[3]{X}$"),
        mpatches.Patch(color=(0.20, 0.80, 0.40), label=r"basin of $\omega^2\sqrt[3]{X}$"),
        mpatches.Patch(
            color="cyan",
            label=r"Voronoi frontier (Thm 5.3)",
        ),
        plt.Line2D(
            [0], [0],
            marker="*",
            color="white",
            markeredgecolor="black",
            linestyle="",
            markersize=10,
            label=r"cube roots of $X$ (Thm 5.1)",
        ),
    ]
    ax.legend(handles=handles, loc="upper right", framealpha=0.88, fontsize=8)

    plt.tight_layout()
    fig.savefig(output, dpi=240)
    plt.close(fig)
    print(f"Saved: {output}  ({res}x{res}, X={X})")


if __name__ == "__main__":
    render()
