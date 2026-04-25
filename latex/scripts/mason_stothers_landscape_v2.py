"""
Mason-Stothers landscape: corrected, audited version.
======================================================

Fixes from v1 (mason_stothers_landscape.py):
  * Quality figure: NaN cells were plotted as q = 0 (deep colour). They
    are now masked transparent (np.ma).
  * Discriminant figure: previous numpy code stripped leading zeros of
    b's coefficient list, which gives the discriminant of the actual
    (variable-degree) polynomial. The symbolic formula 81 x^12 y^12
    (x^4 - 8 x^3 y + 30 x^2 y^2 - 8 x y^3 + y^4) is the discriminant of
    the *formal degree-5* representation of b and is continuous across
    y = -x. We now use the symbolic formula (lambdified) for a
    continuous landscape.
  * Slack figure: the line y = -x (where slack drops from 1 to 0) is
    measure-zero, so the linspace grid mostly missed it. We add an
    explicit high-density sampling along y = -x as overlay points.
  * Labels and titles clarified: which figure plots which quantity, and
    where the Pandrosion connection actually lives (the proof in paper
    60, not the computation here).

Family: a(t) = t^6, c(t) = (t-x)^3 (t-y)^3, b(t) = c - a.
For x, y != 0 and x != y, (a, b, c) is a coprime triple with a + b = c.

The Mason-Stothers theorem (paper 60, proved via the Pandrosion field
F_P = P'/P and Q-tower vanishing) states:
    max(deg a, deg b, deg c) <= deg rad(abc) - 1.
Slack = deg rad(abc) - 1 - max deg >= 0; equality on the extremal locus.

For this family the extremal locus is exactly the line y = -x in R^2
(verified symbolically and confirmed by a 300x300 scan, see
audit_mason_figures.py).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D                     # noqa: F401


# =====================================================================
# Symbolic objects (computed once)
# =====================================================================

def _build_symbolic():
    t, x_s, y_s = sp.symbols('t x y', real=True)
    a = t**6
    c = (t - x_s)**3 * (t - y_s)**3
    b = sp.expand(c - a)
    b_poly = sp.Poly(b, t)
    disc_b = sp.expand(b_poly.discriminant())
    return x_s, y_s, disc_b

_X_SYM, _Y_SYM, _DISC_SYM = _build_symbolic()
_disc_func = sp.lambdify((_X_SYM, _Y_SYM), _DISC_SYM, "numpy")


# =====================================================================
# Numerical functions
# =====================================================================

def slack(x: float, y: float) -> float:
    """deg rad(abc) - 1 - max deg, >= 0 by Mason-Stothers."""
    if abs(x) < 0.05 or abs(y) < 0.05 or abs(x - y) < 0.05:
        return np.nan
    c_roots = [x, x, x, y, y, y]
    c_coeffs = np.poly(c_roots).astype(np.complex128)
    b_coeffs = c_coeffs.copy()
    b_coeffs[0] -= 1.0
    while len(b_coeffs) > 1 and abs(b_coeffs[0]) < 1e-12:
        b_coeffs = b_coeffs[1:]
    if len(b_coeffs) <= 1:
        return np.nan
    b_roots = np.roots(b_coeffs)
    all_roots = [0.0 + 0.0j] + [complex(r, 0) for r in c_roots] + list(b_roots)
    distinct: list[complex] = []
    for r in all_roots:
        if not any(abs(r - d) < 1e-5 for d in distinct):
            distinct.append(complex(r))
    return float(len(distinct) - 1 - 6)


def quality(x: float, y: float) -> float:
    """polynomial-abc quality q = 6 / (deg rad - 1).
    By Mason-Stothers q <= 1; q = 1 at extremals."""
    s = slack(x, y)
    if not np.isfinite(s):
        return np.nan
    deg_rad = int(s + 7)
    if deg_rad <= 1:
        return np.nan
    return 6.0 / (deg_rad - 1)


def log_disc_formal(x: float, y: float) -> float:
    """log10|disc(b)| using the formal degree-5 symbolic formula.
    Continuous across y = -x (unlike the numpy-stripped version)."""
    if abs(x) < 0.05 or abs(y) < 0.05 or abs(x - y) < 0.05:
        return np.nan
    val = _disc_func(x, y)
    if abs(val) < 1e-30:
        return np.nan
    return float(np.log10(abs(val)))


# =====================================================================
# ABC particles
# =====================================================================

ABC_PARTICLES = [
    ("extremal (1, -1)",  1.0, -1.0),
    ("extremal (2, -2)",  2.0, -2.0),
    ("extremal (-1.5, 1.5)", -1.5, 1.5),
    ("extremal (0.7, -0.7)", 0.7, -0.7),
    ("generic (1, 2)",   1.0,  2.0),
    ("generic (2, 1)",   2.0,  1.0),
    ("generic (-2, 1)", -2.0,  1.0),
    ("generic (1.5, -0.5)", 1.5, -0.5),
]


# =====================================================================
# Build grids
# =====================================================================

def build_grid(N: int = 110, span: float = 3.0):
    xs = np.linspace(-span, span, N)
    ys = np.linspace(-span, span, N)
    S = np.full((N, N), np.nan)
    Q = np.full((N, N), np.nan)
    LD = np.full((N, N), np.nan)
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            S[j, i] = slack(x, y)
            Q[j, i] = quality(x, y)
            LD[j, i] = log_disc_formal(x, y)
    return xs, ys, S, Q, LD


# =====================================================================
# Figure 1: slack landscape (with extremal-line overlay)
# =====================================================================

def plot_slack(xs, ys, S):
    X, Y = np.meshgrid(xs, ys)
    S_masked = np.ma.masked_invalid(S)
    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        X, Y, S_masked,
        cmap="magma_r", alpha=0.92, edgecolor="none",
        rstride=1, cstride=1, vmin=0, vmax=1,
    )
    # Mason-Stothers floor: plane at z = 0
    ax.plot_surface(
        X, Y, np.zeros_like(X),
        color="gold", alpha=0.45, edgecolor="none",
    )
    # Extremal locus y = -x explicitly: dense sampling along the line
    line_t = np.linspace(-2.8, 2.8, 200)
    line_x = line_t
    line_y = -line_t
    line_z = np.zeros_like(line_t)
    ax.plot(line_x, line_y, line_z + 0.02,
            color="red", linewidth=2.8,
            label="Extremal locus y = -x  (slack = 0)")
    # Particles
    px = [p[1] for p in ABC_PARTICLES]
    py = [p[2] for p in ABC_PARTICLES]
    pz = []
    pcol = []
    for label, x, y in ABC_PARTICLES:
        s = slack(x, y)
        pz.append((s if np.isfinite(s) else 0) + 0.05)
        pcol.append("lime" if "extremal" in label else "cyan")
    ax.scatter(px, py, pz, c=pcol, s=110,
               edgecolor="black", linewidth=0.8, zorder=10)
    # Aesthetics
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_zlabel("Mason-Stothers slack")
    ax.set_title(
        "Mason-Stothers slack for "
        r"$a=t^6,\ c=(t-x)^3(t-y)^3,\ b=c-a$"
        + "\n(slack ≥ 0 by Mason-Stothers; gold plane = saturation; "
        + "red line y=-x = extremal locus)"
    )
    ax.view_init(elev=24, azim=-50)
    ax.set_zlim(-0.05, 1.3)
    fig.colorbar(surf, shrink=0.55, aspect=12, pad=0.08, label="slack")
    ax.legend(loc="upper left", fontsize=9)
    out = Path(__file__).resolve().parents[1] / "figures" / "ms_v2_slack.png"
    plt.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out


# =====================================================================
# Figure 2: polynomial-abc quality (with NaN masked)
# =====================================================================

def plot_quality(xs, ys, Q):
    X, Y = np.meshgrid(xs, ys)
    Q_masked = np.ma.masked_invalid(Q)
    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        X, Y, Q_masked,
        cmap="viridis", alpha=0.90, edgecolor="none",
        rstride=1, cstride=1, vmin=0.85, vmax=1.0,
    )
    # Plane at q = 1: Mason-Stothers ceiling
    ax.plot_surface(
        X, Y, np.ones_like(X),
        color="gold", alpha=0.40, edgecolor="none",
    )
    # Extremal locus y = -x at q = 1
    line_t = np.linspace(-2.8, 2.8, 200)
    ax.plot(line_t, -line_t, np.ones_like(line_t) - 0.005,
            color="red", linewidth=2.8,
            label="Extremal locus y = -x  (q = 1)")
    # Particles
    px = [p[1] for p in ABC_PARTICLES]
    py = [p[2] for p in ABC_PARTICLES]
    pz = []
    pcol = []
    for label, x, y in ABC_PARTICLES:
        q = quality(x, y)
        pz.append(q + 0.005 if np.isfinite(q) else 0.85)
        pcol.append("lime" if "extremal" in label else "cyan")
    ax.scatter(px, py, pz, c=pcol, s=110,
               edgecolor="black", linewidth=0.8, zorder=10)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_zlabel("polynomial-abc quality  q")
    ax.set_title(
        "Polynomial abc-quality  q = 6 / (deg rad(abc) - 1)\n"
        + "(q ≤ 1 by Mason-Stothers; gold plane = bound saturation; "
        + "red line y=-x: q=1 extremals; NaN cells masked)"
    )
    ax.view_init(elev=24, azim=-50)
    ax.set_zlim(0.84, 1.02)
    fig.colorbar(surf, shrink=0.55, aspect=12, pad=0.08, label="q")
    ax.legend(loc="upper left", fontsize=9)
    out = Path(__file__).resolve().parents[1] / "figures" / "ms_v2_quality.png"
    plt.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out


# =====================================================================
# Figure 3: log|disc(b)| landscape (formal degree-5, smooth)
# =====================================================================

def plot_disc(xs, ys, LD):
    X, Y = np.meshgrid(xs, ys)
    LD_masked = np.ma.masked_invalid(LD)
    # Clip to a sensible visual range
    lo, hi = np.nanpercentile(LD, [2, 98])
    LD_clip = np.clip(LD_masked, lo, hi)
    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        X, Y, LD_clip,
        cmap="magma_r", alpha=0.92, edgecolor="none",
        rstride=1, cstride=1, vmin=lo, vmax=hi,
    )
    # Median floor as a soft reference plane
    floor = float(np.nanmedian(LD_clip))
    ax.plot_surface(
        X, Y, np.full_like(X, floor),
        color="gold", alpha=0.30, edgecolor="none",
    )
    # Extremal locus y = -x (the line where MS is saturated)
    line_t = np.linspace(-2.8, 2.8, 200)
    line_z = []
    for t_v in line_t:
        ld = log_disc_formal(t_v, -t_v)
        line_z.append(ld if np.isfinite(ld) else floor)
    ax.plot(line_t, -line_t, np.array(line_z) + 0.5,
            color="red", linewidth=2.8,
            label="Extremal locus y = -x")
    # Particles
    px, py, pz, pcol = [], [], [], []
    for label, x, y in ABC_PARTICLES:
        ld = log_disc_formal(x, y)
        if np.isfinite(ld):
            px.append(x); py.append(y); pz.append(ld + 0.4)
            pcol.append("lime" if "extremal" in label else "cyan")
    ax.scatter(px, py, pz, c=pcol, s=110,
               edgecolor="black", linewidth=0.8, zorder=10)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_zlabel(r"$\log_{10}|\mathrm{disc}_5(b)|$")
    ax.set_title(
        "Continuous proxy: " + r"$\log_{10}|\mathrm{disc}_5(b)|$"
        + " (formal degree-5 discriminant; sympy-lambdified)\n"
        + "Smooth across y = -x; the function is integer-valued slack's"
        + " continuous companion."
    )
    ax.view_init(elev=24, azim=-50)
    fig.colorbar(surf, shrink=0.55, aspect=12, pad=0.08,
                 label=r"$\log_{10}|\mathrm{disc}_5(b)|$")
    ax.legend(loc="upper left", fontsize=9)
    out = Path(__file__).resolve().parents[1] / "figures" / "ms_v2_disc.png"
    plt.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out


# =====================================================================
# Main
# =====================================================================

def main() -> None:
    print("Symbolic discriminant formula (formal degree-5):")
    print("  disc_5(b) =", sp.factor(_DISC_SYM))
    print()
    print("Building grids...")
    xs, ys, S, Q, LD = build_grid(N=110, span=3.0)
    print("Plotting figure 1 (slack)...")
    p1 = plot_slack(xs, ys, S)
    print(f"  saved {p1}")
    print("Plotting figure 2 (quality)...")
    p2 = plot_quality(xs, ys, Q)
    print(f"  saved {p2}")
    print("Plotting figure 3 (log disc, formal)...")
    p3 = plot_disc(xs, ys, LD)
    print(f"  saved {p3}")
    print()
    print("ABC particles (audit table):")
    print(f"  {'label':<25s} {'slack':>6s} {'quality':>9s} {'log disc_5':>11s}")
    for label, x, y in ABC_PARTICLES:
        s = slack(x, y); q = quality(x, y); ld = log_disc_formal(x, y)
        print(f"  {label:<25s} {s:>6.0f} {q:>9.4f} {ld:>11.4f}")


if __name__ == "__main__":
    main()
