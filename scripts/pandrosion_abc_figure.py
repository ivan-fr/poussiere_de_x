"""
Pandrosion <-> ABC Conjecture : figure HD illustrant les liens formalisés
dans previous_lean_legacy/Pandrosion/Abc*.lean.

Sortie : articles/pandrosion_abc_links.png (3840 x 2160, 300 dpi)
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
import numpy as np


# --------------------------------------------------------------------------
# Palette & style
# --------------------------------------------------------------------------
BG         = "#0B1020"
PANEL      = "#141A2E"
PANEL_EDGE = "#2A3458"
GOLD       = "#E9C46A"
TEAL       = "#2AB7A9"
CORAL      = "#E76F51"
VIOLET     = "#A06CD5"
BLUE       = "#4CC9F0"
GREEN      = "#90E39A"
WHITE      = "#F5F5F7"
DIM        = "#9AA3C4"

plt.rcParams.update({
    "font.family": "DejaVu Serif",
    "mathtext.fontset": "cm",
    "text.color": WHITE,
    "axes.edgecolor": PANEL_EDGE,
    "axes.labelcolor": WHITE,
    "xtick.color": DIM,
    "ytick.color": DIM,
    "savefig.facecolor": BG,
    "figure.facecolor": BG,
    "axes.facecolor": BG,
})


# --------------------------------------------------------------------------
# Figure skeleton
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(19.2, 10.8), dpi=200)  # 3840x2160 at dpi=200

# Main grid : title band + two content rows
gs = fig.add_gridspec(
    nrows=3, ncols=3,
    height_ratios=[0.12, 0.55, 0.45],
    width_ratios=[1.0, 1.15, 1.0],
    hspace=0.30, wspace=0.18,
    left=0.035, right=0.965, top=0.965, bottom=0.045,
)

# ---- Title panel ---------------------------------------------------------
ax_title = fig.add_subplot(gs[0, :])
ax_title.axis("off")
ax_title.text(
    0.5, 0.72,
    "Pandrosion  $\\longleftrightarrow$  Conjecture ABC",
    ha="center", va="center",
    fontsize=34, fontweight="bold", color=GOLD,
    transform=ax_title.transAxes,
)
ax_title.text(
    0.5, 0.22,
    "Formalisation Lean 4  —  le triplet orbital $(X B_n^{3},\\, d_n,\\, A_n^{3})$  comme substrat effectif de  $a+b=c$",
    ha="center", va="center",
    fontsize=16, color=DIM, style="italic",
    transform=ax_title.transAxes,
)


# --------------------------------------------------------------------------
# Helper : rounded panel
# --------------------------------------------------------------------------
def panel(ax, title, color=GOLD):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    # soft rounded background
    rect = FancyBboxPatch(
        (0.01, 0.01), 0.98, 0.98,
        boxstyle="round,pad=0.01,rounding_size=0.025",
        linewidth=1.4, edgecolor=PANEL_EDGE, facecolor=PANEL,
        transform=ax.transAxes, zorder=0,
    )
    ax.add_patch(rect)
    ax.text(
        0.035, 0.935, title,
        transform=ax.transAxes, ha="left", va="center",
        fontsize=15, fontweight="bold", color=color,
    )
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)


# --------------------------------------------------------------------------
# PANEL 1 — Le triplet orbital = triplet abc
# --------------------------------------------------------------------------
ax1 = fig.add_subplot(gs[1, 0])
panel(ax1, "1 · Triplet orbital  =  triplet abc", GOLD)

# schematic triple
cx, cy, R = 0.5, 0.5, 0.26
angles = [np.pi / 2, np.pi / 2 + 2 * np.pi / 3, np.pi / 2 + 4 * np.pi / 3]
labels = [
    ("$c_n = A_n^{3}$", CORAL),
    ("$a_n = X\\, B_n^{3}$", TEAL),
    ("$b_n = d_n$", GOLD),
]
pts = []
for (lbl, col), th in zip(labels, angles):
    x, y = cx + R * np.cos(th), cy + R * np.sin(th)
    pts.append((x, y))
    ax1.add_patch(plt.Circle((x, y), 0.075, color=col, alpha=0.9, zorder=3))
    ax1.text(x, y, "", ha="center", va="center")
    ax1.text(
        x + 0.13 * np.cos(th), y + 0.13 * np.sin(th),
        lbl, ha="center", va="center",
        fontsize=14, color=col, fontweight="bold",
    )

# triangle edges
for i in range(3):
    (x0, y0), (x1, y1) = pts[i], pts[(i + 1) % 3]
    ax1.plot([x0, x1], [y0, y1], color=WHITE, alpha=0.35, lw=1.5, zorder=2)

# central equation
ax1.text(
    0.5, 0.48,
    "$a_n + b_n = c_n$",
    ha="center", va="center",
    fontsize=16, color=WHITE, fontweight="bold",
)

# identity below
ax1.text(
    0.5, 0.18,
    "$A_n^{3} \\;=\\; X\\, B_n^{3} \\;+\\; d_n \\quad(\\text{ring})$",
    ha="center", va="center", fontsize=13, color=DIM,
)
ax1.text(
    0.5, 0.09,
    "Lean : $\\mathtt{orbital\\_triple\\_sum}$  ·  AbcOrbital.lean §1800",
    ha="center", va="center", fontsize=10, color=VIOLET, style="italic",
)


# --------------------------------------------------------------------------
# PANEL 2 — Sandwich abc-orbital + trajectoire X = 2
# --------------------------------------------------------------------------
ax2 = fig.add_subplot(gs[1, 1])
panel(ax2, "2 · Sandwich abc-orbital   $2^{n} \\leq |d_n| \\leq |d_0|\\prod|\\Phi_k|$", GOLD)

# Inset axes inside the panel
inset = ax2.inset_axes([0.09, 0.15, 0.84, 0.68])
inset.set_facecolor(PANEL)

# --- Numerical demo : orbit at X = 2 starting from (A0,B0) = (1,1) -------
def pandrosion_step(A, B, X):
    Ap = A * (A ** 3 + 4 * X * B ** 3)
    Bp = B * (3 * A ** 3 + 2 * X * B ** 3)
    return Ap, Bp

import math

X = 2
A, B = 1, 1
orbit = [(A, B)]
for _ in range(4):
    A, B = pandrosion_step(A, B, X)
    orbit.append((A, B))

d_vals  = [a ** 3 - X * b ** 3 for (a, b) in orbit]
abs_phi = [abs(d_vals[k + 1] // d_vals[k]) for k in range(len(d_vals) - 1)]

n_arr = np.arange(len(d_vals))

# Work in log10 space to avoid overflow: plot log10(|d_n|), log10(2^n), log10(upper)
log_d     = np.array([math.log10(abs(dv)) if abs(dv) > 0 else 0.0 for dv in d_vals])
log_lower = np.array([n * math.log10(2) for n in n_arr])
log_upper = np.array([
    math.log10(abs(d_vals[0])) + sum(math.log10(p) for p in abs_phi[:n])
    for n in n_arr
])

inset.plot(n_arr, log_lower, "-o", color=TEAL, lw=2.2, markersize=7,
           label=r"borne inf.  $\log_{10} 2^{n}$")
inset.plot(n_arr, log_d, "-o", color=GOLD, lw=2.8, markersize=9,
           label=r"$\log_{10}|d_n|$  (orbite $X{=}2$)")
inset.plot(n_arr, log_upper, "-o", color=CORAL, lw=2.2, markersize=7,
           label=r"borne sup.  $\log_{10}(|d_0|\prod|\Phi_k|)$")

# shade the sandwich between lower and upper
inset.fill_between(n_arr, log_lower, log_upper, color=GOLD, alpha=0.08, zorder=0)

inset.set_xlabel("étape  $n$", fontsize=12, color=WHITE)
inset.set_ylabel(r"$\log_{10}$ de l'amplitude", fontsize=12, color=WHITE)
inset.grid(True, which="both", color=PANEL_EDGE, alpha=0.6, lw=0.8)
inset.tick_params(colors=DIM, labelsize=10)
for s in inset.spines.values():
    s.set_color(PANEL_EDGE)
inset.legend(
    facecolor=PANEL, edgecolor=PANEL_EDGE, labelcolor=WHITE,
    fontsize=10, loc="upper left",
)

ax2.text(
    0.5, 0.055,
    "Certificat $X{=}2$, pas 1 :  $\\;686 + 43 = 729\\;$  ($2\\cdot 7^{3} + 43 = 9^{3}$)",
    ha="center", va="center", fontsize=12, color=GREEN, fontweight="bold",
)


# --------------------------------------------------------------------------
# PANEL 3 — Loi de pas radicale : Φ_n concentre rad(abc)
# --------------------------------------------------------------------------
ax3 = fig.add_subplot(gs[1, 2])
panel(ax3, "3 · Loi de pas :  $\\Phi_n$  concentre le radical", GOLD)

# Diagram : (A_n,B_n) -> (A_{n+1},B_{n+1}) with Phi multiplier
def node(ax, x, y, txt, color):
    ax.add_patch(FancyBboxPatch(
        (x - 0.14, y - 0.05), 0.28, 0.10,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        linewidth=1.2, edgecolor=color, facecolor=BG, alpha=0.92,
    ))
    ax.text(x, y, txt, ha="center", va="center",
            fontsize=12, color=color, fontweight="bold")

node(ax3, 0.22, 0.70, "$(A_n,\\,B_n)$", TEAL)
node(ax3, 0.78, 0.70, "$(A_{n+1},\\,B_{n+1})$", CORAL)

arrow = FancyArrowPatch(
    (0.37, 0.70), (0.63, 0.70),
    arrowstyle="-|>", mutation_scale=22,
    color=GOLD, lw=2.2,
)
ax3.add_patch(arrow)
ax3.text(0.5, 0.78, "Pandrosion step", ha="center", fontsize=11, color=DIM, style="italic")

# key identity
ax3.text(
    0.5, 0.52,
    "$A_{n+1}^{3} \\;=\\; X\\, B_{n+1}^{3} \\;+\\; d_n \\cdot \\Phi_n$",
    ha="center", va="center", fontsize=14, color=WHITE, fontweight="bold",
)

ax3.text(
    0.5, 0.40,
    "$d_n \\;=\\; d_0 \\cdot \\prod_{k<n}\\Phi_k$",
    ha="center", va="center", fontsize=13, color=GOLD,
)

# radical concentration text
ax3.text(
    0.5, 0.25,
    "Tout l'accroissement radical vit dans\n$\\Phi_n$  (polynôme de degré 9 en $A,B$)",
    ha="center", va="center", fontsize=11, color=DIM,
)

ax3.text(
    0.5, 0.09,
    "Lean : $\\mathtt{orbital\\_triple\\_step}$, $\\mathtt{abc\\_orbital\\_factorization}$",
    ha="center", va="center", fontsize=10, color=VIOLET, style="italic",
)


# --------------------------------------------------------------------------
# PANEL 4 — Bridge : full abc  <=>  fluid-mechanics escape
# --------------------------------------------------------------------------
ax4 = fig.add_subplot(gs[2, :2])
panel(ax4, "4 · Pont fluido-mécanique :  $\\text{abc complet}\\;\\Leftrightarrow\\;$  borne d'échappement", BLUE)

# Left box : arithmetic side
ax4.add_patch(FancyBboxPatch(
    (0.035, 0.18), 0.36, 0.60,
    boxstyle="round,pad=0.01,rounding_size=0.02",
    linewidth=1.4, edgecolor=TEAL, facecolor=BG, alpha=0.6,
))
ax4.text(0.215, 0.70, "Côté arithmétique", ha="center", fontsize=13, fontweight="bold", color=TEAL)
ax4.text(
    0.215, 0.52,
    "$\\forall\\,\\varepsilon>0\\;\\exists\\,K(\\varepsilon)\\;\\forall\\,(a,b,c)$\ncoprimes, $a+b=c$ :",
    ha="center", va="center", fontsize=11, color=WHITE,
)
ax4.text(
    0.215, 0.33,
    "$c \\;\\leq\\; K(\\varepsilon)\\cdot \\mathrm{rad}(abc)^{1+\\varepsilon}$",
    ha="center", va="center", fontsize=13, color=GOLD, fontweight="bold",
)
ax4.text(0.215, 0.22, "$\\mathtt{full\\_abc\\_conjecture}$",
         ha="center", fontsize=10, color=VIOLET, style="italic")

# Right box : dynamical side
ax4.add_patch(FancyBboxPatch(
    (0.605, 0.18), 0.36, 0.60,
    boxstyle="round,pad=0.01,rounding_size=0.02",
    linewidth=1.4, edgecolor=CORAL, facecolor=BG, alpha=0.6,
))
ax4.text(0.785, 0.70, "Côté dynamique (Pandrosion)", ha="center", fontsize=13, fontweight="bold", color=CORAL)
ax4.text(
    0.785, 0.52,
    "plongement  $P\\mapsto(a,b,c)\\in\\mathbb{R}^{3}$\ndans la métrique Pandrosion",
    ha="center", va="center", fontsize=11, color=WHITE,
)
ax4.text(
    0.785, 0.33,
    "borne d'échappement fluide\nsur toute $\\mathtt{AbcParticle}$",
    ha="center", va="center", fontsize=12, color=GOLD, fontweight="bold",
)
ax4.text(0.785, 0.22, "$\\mathtt{abc\\_fluid\\_mechanics\\_bridge}$",
         ha="center", fontsize=10, color=VIOLET, style="italic")

# Double arrow between them
dbl = FancyArrowPatch(
    (0.405, 0.48), (0.595, 0.48),
    arrowstyle="<|-|>", mutation_scale=28,
    color=BLUE, lw=3.0,
)
ax4.add_patch(dbl)
ax4.text(0.5, 0.57, "$\\Longleftrightarrow$", ha="center", va="center",
         fontsize=22, color=BLUE, fontweight="bold")
ax4.text(0.5, 0.40, "équivalence logique\nprouvée en Lean", ha="center",
         fontsize=10, color=DIM, style="italic")

ax4.text(
    0.5, 0.09,
    "AbcTopologicalEmbedding.lean §4403  ·  AbcGlobalFrontier.lean §4300  (Prop, 0 axiome, 0 sorry)",
    ha="center", va="center", fontsize=10, color=DIM,
)


# --------------------------------------------------------------------------
# PANEL 5 — Graphe de dépendances des modules Lean
# --------------------------------------------------------------------------
ax5 = fig.add_subplot(gs[2, 2])
panel(ax5, "5 · Graphe des modules Lean", VIOLET)

nodes = {
    "ThueBridge":             (0.18, 0.82),
    "HallOrbital":            (0.18, 0.56),
    "AbcOrbital":             (0.52, 0.69),
    "StewartYuOrbital":       (0.85, 0.82),
    "ErdosMahlerOrbital":     (0.85, 0.56),
    "AbcGlobalFrontier":      (0.30, 0.28),
    "AbcTopologicalEmbed.":   (0.72, 0.28),
}
colors = {
    "AbcOrbital":            GOLD,
    "AbcGlobalFrontier":     BLUE,
    "AbcTopologicalEmbed.":  BLUE,
}

edges = [
    ("ThueBridge",         "AbcOrbital"),
    ("HallOrbital",        "AbcOrbital"),
    ("AbcOrbital",         "StewartYuOrbital"),
    ("AbcOrbital",         "ErdosMahlerOrbital"),
    ("AbcGlobalFrontier",  "AbcTopologicalEmbed."),
]

for a, b in edges:
    x0, y0 = nodes[a]; x1, y1 = nodes[b]
    ax5.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="-|>", mutation_scale=14,
        color=DIM, lw=1.2, alpha=0.8,
    ))

for name, (x, y) in nodes.items():
    c = colors.get(name, TEAL)
    ax5.add_patch(FancyBboxPatch(
        (x - 0.145, y - 0.042), 0.29, 0.084,
        boxstyle="round,pad=0.005,rounding_size=0.02",
        linewidth=1.3, edgecolor=c, facecolor=BG, alpha=0.92,
    ))
    ax5.text(x, y, name, ha="center", va="center",
             fontsize=9.5, color=c, fontweight="bold")

# legend
legend_items = [
    Line2D([0], [0], marker="s", color="none", markerfacecolor=GOLD,
           markersize=10, label="triplet orbital (abc effectif)"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=BLUE,
           markersize=10, label="abc complet (frontière)"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=TEAL,
           markersize=10, label="support Thue / Hall"),
]
ax5.legend(
    handles=legend_items, loc="lower center", ncol=1,
    facecolor=PANEL, edgecolor=PANEL_EDGE, labelcolor=WHITE,
    fontsize=9, bbox_to_anchor=(0.5, 0.005),
)


# --------------------------------------------------------------------------
# Footer
# --------------------------------------------------------------------------
fig.text(
    0.5, 0.012,
    "Source : previous_lean_legacy/Pandrosion/{AbcOrbital, AbcGlobalFrontier, AbcTopologicalEmbedding, StewartYuOrbital, ErdosMahlerOrbital, CoprimalityIsolation}.lean",
    ha="center", va="bottom", fontsize=9, color=DIM, style="italic",
)

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
out_dir = "/Users/ivanbesevic/Documents/poussiere/articles"
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "pandrosion_abc_links.png")
fig.savefig(out_path, dpi=200, facecolor=BG, bbox_inches="tight", pad_inches=0.15)
print(f"Saved: {out_path}")
