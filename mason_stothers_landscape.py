"""
Topological encoding of the polynomial abc theorem (Mason-Stothers)
====================================================================

This script produces a 3D landscape visualisation of the Mason-Stothers
theorem (= ABC for polynomials over C[t]), in the spirit of the
"Topological Encoding of the ABC Conjecture" image but for the
fully-proved polynomial version.

Family of triples
-----------------
   a(t) = t^6
   c(t) = (t-x)^3 (t-y)^3,   parametrised by (x, y) in R^2
   b(t) = c(t) - a(t)
The three are coprime for x, y != 0 and x != y. They satisfy a + b = c
(the abc relation).

Mason-Stothers theorem
----------------------
    max(deg a, deg b, deg c) <= deg rad(abc) - 1
where rad(f) = product over distinct roots of (t - alpha).
Equivalently the slack
    S(x, y) := deg rad(abc) - 1 - max deg >= 0,
with equality at extremal triples (the polynomial analogue of "high
abc-quality triples").

Pandrosion connection (paper 60: 60pandrosion_mason.pdf)
--------------------------------------------------------
The Wronskian W(a,b) = a'b - a b' decomposes via the Pandrosion field
F_P = P'/P as
    W(a,b) / (ab) = F_a - F_b.
At every multiple root alpha of a (multiplicity m), Q-tower vanishing
Q_a(alpha, alpha) = a'(alpha) = 0 forces (t - alpha)^{m-1} | W. The
proof is finished by counting forced zeros against deg W <= deg a +
deg b - 1.

Hence each ABC particle (specific (a,b,c) triple) is encoded
topologically in this picture by:
  - its position (x, y) in the parameter plane,
  - its slack height S(x, y) above the Mason-Stothers floor,
  - its Wronskian zero count = sum (m_i - 1) over multiple roots of a,b,c.

Visual analogy with the integer ABC figure:
  - Yellow plane at z = 0  ↔  the rad(abc)^{1+eps} bound.
  - Surface above this plane  ↔  the "ABC particle cloud".
  - Sinkholes (where S drops toward 0)  ↔  high-quality triples.
  - Color gradient  ↔  abc-quality of each (x, y) configuration.

Output: PNG figure saved to poussiere/.
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D                      # noqa: F401


# ---------------------------------------------------------------------
# Mason-Stothers slack computation
# ---------------------------------------------------------------------

def slack_3_3(x: float, y: float) -> float | float:
    """
    For a = t^6, c = (t-x)^3 (t-y)^3, b = c - a, return
       slack(x, y) = deg rad(abc) - 1 - max deg.

    Returns NaN at degenerate configurations (x or y near 0, x near y).
    """
    if abs(x) < 0.05 or abs(y) < 0.05 or abs(x - y) < 0.05:
        return np.nan
    # Build c via numpy.poly (coefficients high-degree-first)
    c_roots = [x, x, x, y, y, y]
    c_coeffs = np.poly(c_roots).astype(np.complex128)
    # b = c - t^6: subtract 1 from leading coefficient
    b_coeffs = c_coeffs.copy()
    b_coeffs[0] -= 1.0
    # Strip off any leading zeros (in this family the leading
    # coefficient cancels exactly; the next coefficient is -3(x+y) which
    # is generically nonzero)
    while len(b_coeffs) > 1 and abs(b_coeffs[0]) < 1e-12:
        b_coeffs = b_coeffs[1:]
    if len(b_coeffs) <= 1:
        return np.nan
    b_roots = np.roots(b_coeffs)
    a_roots = [0.0 + 0.0j]
    all_roots = a_roots + [complex(r, 0) for r in c_roots] + list(b_roots)
    # Distinct roots within tolerance
    distinct: list[complex] = []
    tol = 1e-5
    for r in all_roots:
        if not any(abs(r - d) < tol for d in distinct):
            distinct.append(complex(r))
    deg_rad = len(distinct)
    max_deg = 6
    return float(deg_rad - 1 - max_deg)


def quality(x: float, y: float) -> float:
    """polynomial-abc quality = max deg / (deg rad - 1).
    By Mason-Stothers, quality <= 1, with equality at extremals."""
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
    tol = 1e-5
    for r in all_roots:
        if not any(abs(r - d) < tol for d in distinct):
            distinct.append(complex(r))
    deg_rad = len(distinct)
    if deg_rad <= 1:
        return np.nan
    return 6.0 / (deg_rad - 1)


# ---------------------------------------------------------------------
# Build the landscape on a grid
# ---------------------------------------------------------------------

def build_landscape(N: int = 110, span: float = 3.0):
    xs = np.linspace(-span, span, N)
    ys = np.linspace(-span, span, N)
    S = np.full((N, N), np.nan)
    Q = np.full((N, N), np.nan)
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            S[j, i] = slack_3_3(x, y)
            Q[j, i] = quality(x, y)
    return xs, ys, S, Q


# ---------------------------------------------------------------------
# Specific "ABC particles" (concrete polynomial triples)
# ---------------------------------------------------------------------

def abc_particles() -> list[tuple[str, float, float, float]]:
    """Concrete polynomial triples, including extremal Mason-Stothers ones."""
    candidates = [
        # Extremal locus y = -x (slack = 0, polynomial abc bound saturated):
        ("extremal (1, -1)",  1.0, -1.0),
        ("extremal (2, -2)",  2.0, -2.0),
        ("extremal (-1.5, 1.5)", -1.5, 1.5),
        ("extremal (0.7, -0.7)", 0.7, -0.7),
        # Generic triples (slack = 1):
        ("generic (1, 2)",   1.0,  2.0),
        ("generic (2, 1)",   2.0,  1.0),
        ("generic (-2, 1)", -2.0,  1.0),
        ("generic (1.5, -0.5)", 1.5, -0.5),
    ]
    out = []
    for label, x, y in candidates:
        s = slack_3_3(x, y)
        if np.isfinite(s):
            out.append((label, x, y, s))
    return out


def log_disc_b(x: float, y: float) -> float:
    """log |disc(b)| as a continuous proxy for "approach to extremal"."""
    if abs(x) < 0.05 or abs(y) < 0.05 or abs(x - y) < 0.05:
        return np.nan
    c_roots = [x, x, x, y, y, y]
    c_coeffs = np.poly(c_roots).astype(np.complex128)
    b_coeffs = c_coeffs.copy()
    b_coeffs[0] -= 1.0
    while len(b_coeffs) > 1 and abs(b_coeffs[0]) < 1e-12:
        b_coeffs = b_coeffs[1:]
    if len(b_coeffs) <= 2:
        return np.nan
    # discriminant of b(t) via numpy: |disc| = |Res(b, b')| / |lc(b)|
    bp_coeffs = np.polyder(b_coeffs)
    # Resultant via Sylvester determinant of (b, b')
    n, m = len(b_coeffs) - 1, len(bp_coeffs) - 1
    if n < 1 or m < 1:
        return np.nan
    syl = np.zeros((n + m, n + m), dtype=np.complex128)
    for i in range(m):
        syl[i, i:i + n + 1] = b_coeffs
    for i in range(n):
        syl[m + i, i:i + m + 1] = bp_coeffs
    res = np.linalg.det(syl)
    lc = b_coeffs[0]
    disc = (-1) ** (n * (n - 1) // 2) * res / lc
    return float(np.log10(abs(disc) + 1e-15))


# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------

def plot_landscape() -> str:
    xs, ys, S, Q = build_landscape(N=110, span=3.0)
    X, Y = np.meshgrid(xs, ys)

    # Replace NaN by max value (so the singular ridges read as "high"
    # rather than holes)
    S_plot = np.where(np.isnan(S), np.nanmax(S), S)

    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111, projection="3d")

    # Surface: Mason-Stothers slack
    surf = ax.plot_surface(
        X, Y, S_plot,
        cmap="magma_r", alpha=0.92, edgecolor="none",
        rstride=1, cstride=1, antialiased=True,
        vmin=0, vmax=np.nanmax(S),
    )

    # Mason-Stothers floor: plane at z = 0 (the bound saturation)
    plane = np.zeros_like(S_plot)
    ax.plot_surface(
        X, Y, plane,
        color="gold", alpha=0.45, edgecolor="none",
    )

    # Bound line at fixed x, swept across y, at z = 0
    ax.plot(
        xs, [-3.0] * len(xs), [0.0] * len(xs),
        color="goldenrod", linewidth=2.5,
        label="Mason-Stothers bound: max deg = deg rad - 1",
    )

    # Particles: scatter concrete (x, y) triples with their slack height
    parts = abc_particles()
    px = [p[1] for p in parts]
    py = [p[2] for p in parts]
    pz = [p[3] + 0.15 for p in parts]   # tiny lift so dots are visible
    ax.scatter(
        px, py, pz,
        c="lime", s=90, edgecolor="black", linewidth=0.8, zorder=10,
        label="ABC particles (specific triples)",
    )

    # Aesthetics
    ax.set_xlabel("x  (parameter of c = (t-x)^3 (t-y)^3)", fontsize=10)
    ax.set_ylabel("y", fontsize=10)
    ax.set_zlabel("slack  =  deg rad(abc) - 1 - max deg", fontsize=10)
    ax.set_title(
        "Topological encoding of the polynomial abc theorem (Mason-Stothers)\n"
        r"$a(t) = t^6, \ c(t) = (t-x)^3 (t-y)^3, \ b(t) = c(t)-a(t)$"
        + "\nPandrosion-native via field $F_P = P'/P$ and Q-tower vanishing (paper 60)",
        fontsize=11,
    )
    ax.view_init(elev=22, azim=-35)
    ax.set_zlim(-0.5, np.nanmax(S) + 0.5)

    # Colorbar
    cbar = fig.colorbar(surf, shrink=0.55, aspect=12, pad=0.08)
    cbar.set_label("Mason-Stothers slack", fontsize=9)

    # Legend
    ax.legend(loc="upper left", fontsize=9)

    out_path = "/sessions/festive-great-meitner/mnt/poussiere/mason_stothers_landscape.png"
    plt.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_quality_view() -> str:
    """Companion view: polynomial-abc quality q = 6 / (deg rad - 1)."""
    xs, ys, S, Q = build_landscape(N=110, span=3.0)
    X, Y = np.meshgrid(xs, ys)
    Q_plot = np.where(np.isnan(Q), 0.0, Q)

    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111, projection="3d")

    # Quality surface (always <= 1 by Mason-Stothers)
    surf = ax.plot_surface(
        X, Y, Q_plot,
        cmap="viridis", alpha=0.88, edgecolor="none",
        rstride=1, cstride=1, vmin=0, vmax=1,
    )

    # Plane at q = 1: the Mason-Stothers bound (saturation)
    plane = np.ones_like(Q_plot)
    ax.plot_surface(
        X, Y, plane,
        color="gold", alpha=0.40, edgecolor="none",
    )

    # Particles
    parts = abc_particles()
    px = [p[1] for p in parts]
    py = [p[2] for p in parts]
    pz = []
    for _, x, y, _ in parts:
        q = quality(x, y)
        pz.append(q + 0.02 if np.isfinite(q) else 0)
    ax.scatter(
        px, py, pz,
        c="lime", s=90, edgecolor="black", linewidth=0.8, zorder=10,
        label="ABC particles (specific triples)",
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel(r"polynomial-abc quality  $q = 6 / (deg\, rad\, abc - 1)$")
    ax.set_title(
        "Polynomial abc-quality landscape (Mason-Stothers floor at q = 1)\n"
        + r"sinkholes where the bound is saturated; $q \leq 1$ everywhere by Mason-Stothers"
    )
    ax.view_init(elev=22, azim=-40)
    ax.set_zlim(0.4, 1.05)
    fig.colorbar(surf, shrink=0.55, aspect=12, pad=0.08, label="q")
    ax.legend(loc="upper left", fontsize=9)

    out_path = "/sessions/festive-great-meitner/mnt/poussiere/mason_stothers_quality.png"
    plt.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_disc_landscape() -> str:
    """Continuous landscape: log|disc(b)| over (x, y).
    Diverges (--> -infinity) on the extremal locus where b acquires
    a double root, mimicking the funnel of the integer ABC figure.
    """
    N = 110
    span = 3.0
    xs = np.linspace(-span, span, N)
    ys = np.linspace(-span, span, N)
    Z = np.full((N, N), np.nan)
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            Z[j, i] = log_disc_b(x, y)
    X, Y = np.meshgrid(xs, ys)
    Z_plot = np.where(np.isnan(Z), 0.0, Z)
    Z_plot = np.clip(Z_plot, np.nanpercentile(Z, 5), np.nanpercentile(Z, 99))

    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        X, Y, Z_plot,
        cmap="magma_r", alpha=0.92, edgecolor="none",
        rstride=1, cstride=1, antialiased=True,
    )
    # Mason-Stothers floor: a horizontal plane at the median height
    floor_z = float(np.nanmedian(Z_plot))
    ax.plot_surface(
        X, Y, np.full_like(Z_plot, floor_z),
        color="gold", alpha=0.35, edgecolor="none",
    )
    # Particles
    parts = abc_particles()
    px = [p[1] for p in parts]
    py = [p[2] for p in parts]
    pz = [(log_disc_b(x, y) if np.isfinite(log_disc_b(x, y)) else floor_z) + 0.5
          for _, x, y, _ in parts]
    ax.scatter(
        px, py, pz,
        c="lime", s=110, edgecolor="black", linewidth=0.9, zorder=10,
        label="ABC particles (specific triples)",
    )
    # Extremal locus y = -x:
    line_xs = np.linspace(-2.5, 2.5, 100)
    line_xs = line_xs[np.abs(line_xs) > 0.1]
    line_ys = -line_xs
    line_zs = [log_disc_b(x, y) for x, y in zip(line_xs, line_ys)]
    ax.plot(line_xs, line_ys, line_zs,
            color="red", linewidth=3, alpha=0.9,
            label="Extremal Mason-Stothers locus (y = -x)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel(r"$\log_{10}|\mathrm{disc}(b)|$  (continuous height)")
    ax.set_title(
        "Topological encoding of polynomial abc (continuous proxy)\n"
        + r"$a = t^6$, $c = (t-x)^3 (t-y)^3$, $b = c-a$, height = $\log_{10}|\mathrm{disc}(b)|$"
        + "\nValleys (low height) on the extremal locus y = -x where Mason-Stothers is saturated"
    )
    ax.view_init(elev=22, azim=-50)
    fig.colorbar(surf, shrink=0.55, aspect=12, pad=0.08, label="log|disc(b)|")
    ax.legend(loc="upper left", fontsize=9)
    out_path = "/sessions/festive-great-meitner/mnt/poussiere/mason_stothers_disc.png"
    plt.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    print("Building Mason-Stothers landscape (discrete slack)...")
    p1 = plot_landscape()
    print(f"  saved {p1}")
    print("Building polynomial-abc quality view...")
    p2 = plot_quality_view()
    print(f"  saved {p2}")
    print("Building continuous log|disc(b)| view (smoother)...")
    p3 = plot_disc_landscape()
    print(f"  saved {p3}")
    print()
    print("ABC particles (concrete polynomial triples):")
    for label, x, y, s in abc_particles():
        q = quality(x, y)
        ld = log_disc_b(x, y)
        print(f"  {label:<25s}  slack = {s:>3.0f},  quality = {q:.4f},  "
              f"log|disc(b)| = {ld:.3f}" if np.isfinite(ld) else
              f"  {label:<25s}  slack = {s:>3.0f},  quality = {q:.4f},  "
              f"log|disc(b)| = ---")


if __name__ == "__main__":
    main()
