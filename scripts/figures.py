import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'font.size': 11,
    'text.usetex': False,
})

def S_p(s, p):
    return sum(s**k for k in range(p))

def pandrosion_s(s, p, x):
    return 1 - (x - 1) / (x * S_p(s, p))

def lambda_theoretical(p, x):
    alpha = x**(1/p)
    num = (alpha - 1) * sum(k * alpha**(p-1-k) for k in range(1, p))
    den = sum(alpha**k for k in range(p))
    return num / den


# ═══════════════════════════════════════════════════════════════════
# FIGURE 1  — Pandrosion's geometric construction (p=3, x=2)
# ═══════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(9, 8.5))
fig.patch.set_facecolor('#f5f2eb')
ax.set_facecolor('#f5f2eb')

W, H = 5.0, 5.0
p, x_val = 3, 2

# --- iterations
s = 0.30
s_list = [s]
for _ in range(3):
    s = pandrosion_s(s, p, x_val)
    s_list.append(s)
s_star = x_val**(-1/p)

# --- rectangle with subtle gradient
rect = mpatches.FancyBboxPatch((0, 0), W, H,
    boxstyle='square,pad=0', linewidth=2.0,
    edgecolor='#333333', facecolor='#ede8dc', zorder=1)
ax.add_patch(rect)

# --- main diagonal (Thales) : (0,H) → (W,0)
ax.plot([0, W], [H, 0], color='#444444', linewidth=2.5, zorder=4)

# label diagonale
mid_diag_x, mid_diag_y = W * 0.25, H * 0.75
ax.text(mid_diag_x - 0.55, mid_diag_y + 0.20,
        r'Diagonal — Thales',
        fontsize=9.5, color='#444', ha='left', va='bottom',
        rotation=-45, style='italic')

# --- courbe parabolique v = x · s^(p-1) projetée dans le rectangle
# On trace la courbe y(t) = H(1-t), x(t) = W · x_val · t^(p-1) / x_val = W · t^(p-1)
# Mais v_n = x·s_n^(p-1), donc sur l'axe horizontal la lecture est v_n
# On trace la courbe paramétrique dans le rectangle:
# pour chaque s in [0,1], le point (W · x_val · s^(p-1) / x_val, H·(1-s)) = (W·s^(p-1), H(1-s))
# Mais on veut v_n dans [0, x_val], mappé à [0, W]:  x_plot = v_n * W / x_val
t_curve = np.linspace(0.15, 1.0, 300)
x_curve = x_val * t_curve**(p-1) * (W / x_val)  # = W · t^(p-1)... simplifies to W·t^2
y_curve = H * (1 - t_curve)
# Clip to rectangle
mask = (x_curve >= 0) & (x_curve <= W) & (y_curve >= 0) & (y_curve <= H)
ax.plot(x_curve[mask], y_curve[mask], color='#8855aa', linewidth=2.2, zorder=4,
        linestyle='-', alpha=0.85)
ax.text(W * 0.68, H * 0.08,
        r'Curve $v = x\,s^{p-1}$',
        fontsize=9.5, color='#8855aa', ha='center', va='bottom',
        style='italic')

# --- palette bleue
N = len(s_list)
blue_shades = [plt.cm.Blues(0.30 + 0.60 * i / (N - 1)) for i in range(N)]

for i, s_n in enumerate(s_list):
    u_n = H * (1 - s_n)
    col = blue_shades[i]

    # Intersection with diagonal: diagonal goes from (0,H) to (W,0)
    # y = H - (H/W)*x, donc à y=u_n: x_diag = W*(H-u_n)/H = W*s_n
    x_diag = W * s_n

    # Intersection with parabolic curve:
    # y_curve = H*(1-s_n) = u_n, x_curve = W*s_n^(p-1)
    x_para = W * s_n**(p-1)

    # horizontal line between the two intersection points
    x_left = min(x_diag, x_para)
    x_right = max(x_diag, x_para)
    ax.hlines(u_n, 0, W, colors=col, linewidth=1.2, alpha=0.35, zorder=2,
              linestyle=':')
    # Segment solide entre diagonale et courbe
    ax.plot([x_left, x_right], [u_n, u_n], color=col, linewidth=2.2,
            alpha=0.85, zorder=3)

    # point on diagonal
    ax.plot(x_diag, u_n, 'o', color=col, markersize=8, zorder=6,
            markeredgecolor='white', markeredgewidth=0.8)

    # point on parabolic curve
    ax.plot(x_para, u_n, 's', color=col, markersize=7, zorder=6,
            markeredgecolor='white', markeredgewidth=0.8)

    # label s_n on the left
    ax.text(-0.35, u_n, f'$s_{i}$',
            fontsize=11, va='center', ha='right', color=col)

    # label u_n on the right
    y_off = 0
    if i == 3:
        y_off = -0.22
    ax.text(W + 0.15, u_n + y_off, f'$u_{i}$',
            fontsize=11, va='center', color=col)

    # Projection arrow downward: from parabolic point to x-axis
    v_n = x_val * s_n**(p-1)
    v_x = v_n * (W / x_val)
    ax.annotate('', xy=(v_x, -0.08), xytext=(x_para, u_n),
                arrowprops=dict(arrowstyle='->', color=col,
                                lw=1.0, alpha=0.45,
                                connectionstyle='arc3,rad=0.0',
                                linestyle='dashed'))

    # Label v_n below the rectangle
    ax.text(v_x, -0.30 - 0.28 * i,
            f'$v_{i}={v_n:.3f}$', ha='center', fontsize=9, color=col)

# --- fixed point line u*
u_star = H * (1 - s_star)
ax.hlines(u_star, 0, W, colors='#bb2222', linewidth=2.2,
          linestyle='--', zorder=5)
ax.text(W + 0.15, u_star, r'$u^*$',
        fontsize=13, va='center', color='#bb2222', fontweight='bold')

# Annotations s* and 1-s* on the left edge
ax.annotate('', xy=(-0.10, H), xytext=(-0.10, u_star),
            arrowprops=dict(arrowstyle='<->', color='#bb2222',
                            lw=1.5, mutation_scale=10))
ax.text(-0.60, (H + u_star) / 2, r'$s^*$',
        fontsize=12, va='center', ha='center', color='#bb2222',
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.15', facecolor='#fff8f0',
                  edgecolor='#bb2222', alpha=0.9))

ax.annotate('', xy=(-0.10, u_star), xytext=(-0.10, 0),
            arrowprops=dict(arrowstyle='<->', color='#666',
                            lw=1.2, mutation_scale=10))
ax.text(-0.60, u_star / 2, r'$1\!-\!s^*$',
        fontsize=11, va='center', ha='center', color='#666',
        bbox=dict(boxstyle='round,pad=0.15', facecolor='#f0f0f0',
                  edgecolor='#999', alpha=0.8))

# Fixed point: point on diagonal and on curve
x_diag_star = W * s_star
x_para_star = W * s_star**(p-1)
ax.plot(x_diag_star, u_star, 'o', color='#bb2222', markersize=9, zorder=7,
        markeredgecolor='white', markeredgewidth=1.0)
ax.plot(x_para_star, u_star, 's', color='#bb2222', markersize=8, zorder=7,
        markeredgecolor='white', markeredgewidth=1.0)

# v* below the rectangle
v_star = x_val * s_star**(p-1)
vs_x = v_star * (W / x_val)
ax.annotate('', xy=(vs_x, -0.08), xytext=(x_para_star, u_star),
            arrowprops=dict(arrowstyle='->', color='#bb2222',
                            lw=1.5, alpha=0.6,
                            linestyle='dashed'))
y_vstar = -0.30 - 0.28 * N
ax.text(vs_x, y_vstar,
        r'$v^{\!*}=\sqrt[3]{2}\approx 1.260$', ha='center', fontsize=11,
        color='#bb2222', fontweight='bold')

# --- graduated readout axis (below rectangle)
y_axis = -0.08
ax.plot([0, W], [y_axis, y_axis], 'k-', linewidth=1.5, zorder=2)
for v_tick in [0, 0.5, 1.0, 1.5, 2.0]:
    xt = v_tick * (W / x_val)
    ax.plot(xt, y_axis, 'k|', markersize=7, zorder=3)
    ax.text(xt, y_axis - 0.20, f'{v_tick}', ha='center', fontsize=9.5, color='#444')

# special tick for cbrt(2)
xt_star = v_star * (W / x_val)
ax.plot(xt_star, y_axis, '|', color='#bb2222', markersize=9, zorder=4)
ax.text(xt_star, y_axis - 0.20, r'$\sqrt[3]{2}$', ha='center', fontsize=10,
        color='#bb2222')

# manual legend
legend_items = [
    plt.Line2D([0], [0], color='#444', linewidth=2.2, label='Diagonal (Thales)'),
    plt.Line2D([0], [0], color='#8855aa', linewidth=2.0, label=r'Curve $v=x\,s^{p-1}$'),
    plt.Line2D([0], [0], color='#bb2222', linewidth=2.0, linestyle='--', label=r'Fixed point $u^*$'),
    plt.Line2D([0], [0], marker='o', color='#4477bb', markersize=7,
               linestyle='', label=r'Diagonal intersection'),
    plt.Line2D([0], [0], marker='s', color='#4477bb', markersize=6,
               linestyle='', label=r'Curve intersection'),
]
ax.legend(handles=legend_items, loc='upper right', fontsize=9,
          framealpha=0.85, bbox_to_anchor=(1.0, 1.0))

# --- axes labels
ax.text(W / 2, -0.60 - 0.28 * N,
        r'Output $v_n = x\,s_n^{p-1}$  (converges to $\sqrt[3]{2}$)',
        ha='center', fontsize=11, color='#333', style='italic')

# Dimensions
ax.set_xlim(-0.85, W + 1.05)
ax.set_ylim(-0.60 - 0.28 * (N + 1), H + 0.45)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title(
    r"Pandrosion's construction  ($p=3$, $x=2$)" + '\n'
    r'Successive parallels — convergence to $v^*=\sqrt[3]{2}$',
    fontsize=13, pad=14, fontweight='bold')

plt.tight_layout()
plt.savefig('pandrosion_geometry.pdf', dpi=200, bbox_inches='tight',
            facecolor='#f5f2eb')
plt.savefig('pandrosion_geometry.png', dpi=200, bbox_inches='tight',
            facecolor='#f5f2eb')
print("Geometric figure saved.")
plt.close()


# ═══════════════════════════════════════════════════════════════════
# FIGURE 2  — cobweb + residual profiles
# ═══════════════════════════════════════════════════════════════════

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
fig2.patch.set_facecolor('#fafaf8')

# ── left panel: cobweb ───────────────────────────────────────────
ax1.set_facecolor('#fafaf8')

p_c, x_c = 3, 2
s_star = x_c**(-1/p_c)
s_vals_plot = np.linspace(0.40, 1.00, 800)
h_vals_plot = np.array([pandrosion_s(s, p_c, x_c) for s in s_vals_plot])

ax1.plot(s_vals_plot, h_vals_plot, color='#2255aa', linewidth=2.3,
         label=r'$s_{n+1}=h(s_n)$')
ax1.plot(s_vals_plot, s_vals_plot, 'k--', linewidth=1.2, alpha=0.55,
         label=r'$s_{n+1}=s_n$')
ax1.plot(s_star, s_star, 'o', color='#cc2222', markersize=10, zorder=6,
         label=rf'$s^*=2^{{-1/3}}\approx{s_star:.4f}$')

# cobweb depuis s_0 = 0.50
s0_cob = 0.50
s_cb   = s0_cob
cx, cy = [s_cb], [s_cb]
for _ in range(20):
    hs  = pandrosion_s(s_cb, p_c, x_c)
    cx += [s_cb, hs]
    cy += [hs,   hs]
    s_cb = hs
ax1.plot(cx, cy, '-', color='#e06020', linewidth=1.3, alpha=0.70, zorder=4)
ax1.plot(s0_cob, s0_cob, 's', color='#20a020', markersize=10, zorder=6,
         label=rf'$s_0={s0_cob}$')

lam = lambda_theoretical(p_c, x_c)
ax1.text(0.44, 0.95,
         rf'$\lambda_{{3,2}}\approx{lam:.4f}$',
         fontsize=12, color='#2255aa',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#eef4ff', alpha=0.85))

ax1.set_xlim(0.40, 1.02)
ax1.set_ylim(0.40, 1.02)
ax1.set_xlabel(r'$s_n$', fontsize=13)
ax1.set_ylabel(r'$s_{n+1}$', fontsize=13)
ax1.set_title('Cobweb diagram — Pandrosion convergence\n'
              r'($p=3,\ x=2$,   $s_0=0.5$)', fontsize=11)
ax1.legend(fontsize=9.5, loc='upper left')
ax1.grid(True, alpha=0.25)

# ── right panel: residual profiles ──────────────────────────────
ax2.set_facecolor('#fafaf8')

configs_pand = [
    (3, 2, '#2255aa', r'Pandrosion  $p=3,\ x=2$',  0.50),
    (2, 2, '#22aa55', r'Pandrosion  $p=2,\ x=2$',  0.50),
    (3, 3, '#8833bb', r'Pandrosion  $p=3,\ x=3$',  0.50),
    (4, 2, '#cc8800', r'Pandrosion  $p=4,\ x=2$',  0.50),
]

for p_i, x_i, col, lbl, s0_i in configs_pand:
    alpha_i = x_i**(1/p_i)
    s_i = s0_i
    eps_i = []
    for _ in range(18):
        v_i = x_i * s_i**(p_i - 1)
        eps_i.append(abs(v_i - alpha_i))
        s_i = pandrosion_s(s_i, p_i, x_i)
    eps_arr = np.array(eps_i)
    mask = eps_arr > 1e-14
    ax2.semilogy(np.where(mask)[0], eps_arr[mask],
                 'o-', color=col, label=lbl, linewidth=1.9, markersize=5)

# Newton p=3, x=2
u_n = 1.5
alpha_n = 2**(1/3)
eps_newton = []
for _ in range(8):
    eps_newton.append(abs(u_n - alpha_n))
    u_n = (2*u_n + 2/u_n**2) / 3
ax2.semilogy(range(len(eps_newton)), eps_newton,
             's--', color='#cc2222', label=r'Newton  $p=3,\ x=2$',
             linewidth=2.0, markersize=8)

ax2.set_xlabel(r'Step $n$', fontsize=13)
ax2.set_ylabel(r'$\varepsilon_n = |v_n - \alpha|$', fontsize=13)
ax2.set_title('Residual profiles: Pandrosion vs Newton', fontsize=11)
ax2.legend(fontsize=9.5)
ax2.grid(True, alpha=0.25)

plt.tight_layout(pad=2.2)
plt.savefig('pandrosion_figures.pdf', dpi=200, bbox_inches='tight',
            facecolor='#fafaf8')
plt.savefig('pandrosion_figures.png', dpi=150, bbox_inches='tight',
            facecolor='#fafaf8')
print("Cobweb + residual profiles figure saved.")
plt.close()


# ═══════════════════════════════════════════════════════════════════
# FIGURE 3  — Steffensen–Pandrosion vs Newton (quadratic convergence)
# ═══════════════════════════════════════════════════════════════════

fig3, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 6.5),
                                          gridspec_kw={'width_ratios': [1.2, 1]})
fig3.patch.set_facecolor('#fafaf8')

p_sp, x_sp = 3, 2
alpha_sp = x_sp**(1/p_sp)
s_star_sp = 1 / alpha_sp

# ── left panel: error vs step ────────────────────────────────────
ax_left.set_facecolor('#fafaf8')

# Steffensen–Pandrosion
def steffensen_step(s, p, x):
    s0 = s
    s1 = pandrosion_s(s0, p, x)
    s2 = pandrosion_s(s1, p, x)
    denom = s2 - 2*s1 + s0
    if abs(denom) < 1e-30:
        return s2
    return s0 - (s1 - s0)**2 / denom

s_steff = 0.875
eps_steff = []
s_steff_vals = [s_steff]
for i in range(5):
    v = x_sp * s_steff**(p_sp - 1)
    eps_steff.append(abs(v - alpha_sp))
    if abs(s_steff - s_star_sp) < 1e-16:
        break
    s_steff = steffensen_step(s_steff, p_sp, x_sp)
    s_steff_vals.append(s_steff)
# Final value
v = x_sp * s_steff**(p_sp - 1)
eps_steff.append(abs(v - alpha_sp))

# Newton: u_{n+1} = ((p-1)*u_n + x/u_n^{p-1}) / p
u_newt = 1.5  # same starting approximation as v_0 for Steffensen
eps_newton = []
for i in range(8):
    eps_newton.append(abs(u_newt - alpha_sp))
    u_newt = ((p_sp - 1) * u_newt + x_sp / u_newt**(p_sp - 1)) / p_sp

# Plain Pandrosion (linear)
s_pand = 0.875
eps_pand = []
for i in range(18):
    v = x_sp * s_pand**(p_sp - 1)
    eps_pand.append(abs(v - alpha_sp))
    s_pand = pandrosion_s(s_pand, p_sp, x_sp)

# Plot — filter out zeros for log scale
def plot_filtered(ax, data, **kwargs):
    ns = list(range(len(data)))
    filtered_n = [n for n, e in zip(ns, data) if e > 1e-17]
    filtered_e = [e for e in data if e > 1e-17]
    ax.semilogy(filtered_n, filtered_e, **kwargs)

plot_filtered(ax_left, eps_pand,
              color='#aaaaaa', marker='o', markersize=4, linewidth=1.2,
              alpha=0.5, label=r'Pandrosion (linear)', zorder=3)

plot_filtered(ax_left, eps_newton,
              color='#cc2222', marker='s', markersize=9, linewidth=2.2,
              markeredgecolor='white', markeredgewidth=0.8,
              label=r'Newton', linestyle='--', zorder=4)

plot_filtered(ax_left, eps_steff,
              color='#2255aa', marker='D', markersize=10, linewidth=2.5,
              markeredgecolor='white', markeredgewidth=1.0,
              label=r'Steffensen–Pandrosion', zorder=5)

# Machine precision line
ax_left.axhline(y=2.2e-16, color='#999', linewidth=1.0, linestyle=':',
                alpha=0.7, zorder=1)
ax_left.text(8.2, 3e-16, 'Machine\nprecision', fontsize=8, color='#999',
             style='italic', ha='left', va='bottom')

# Annotations — positioned to avoid overlap
ax_left.annotate(r'$K_S \approx 0.013$' + '\n' + r'3 steps to $10^{-16}$',
                 xy=(3, eps_steff[3] if len(eps_steff) > 3 else 1e-16),
                 xytext=(5.5, 1e-12),
                 fontsize=10, color='#2255aa',
                 arrowprops=dict(arrowstyle='->', color='#2255aa', lw=1.5),
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='#eef4ff',
                          edgecolor='#2255aa', alpha=0.9))

ax_left.annotate(r'$K_N \approx 0.794$' + '\n' + r'5 steps to $10^{-15}$',
                 xy=(4, eps_newton[4] if len(eps_newton) > 4 else 1e-7),
                 xytext=(5.5, 1e-4),
                 fontsize=10, color='#cc2222',
                 arrowprops=dict(arrowstyle='->', color='#cc2222', lw=1.5),
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff0f0',
                          edgecolor='#cc2222', alpha=0.9))

ax_left.set_xlabel(r'Step $n$', fontsize=13)
ax_left.set_ylabel(r'$\varepsilon_n = |v_n - \sqrt[3]{2}|$', fontsize=13)
ax_left.set_title(r'Convergence comparison ($p=3$, $x=2$)', fontsize=13,
                   fontweight='bold')
ax_left.legend(fontsize=10, loc='upper right')
ax_left.set_ylim(1e-17, 5)
ax_left.set_xlim(-0.3, 10)
ax_left.grid(True, alpha=0.25, which='both')

# ── right panel: quadratic constant comparison bar chart ─────────
ax_right.set_facecolor('#fafaf8')

methods = ['Steffensen–\nPandrosion', 'Newton']
K_values = [0.013, 0.794]
colors_bar = ['#2255aa', '#cc2222']
bars = ax_right.bar(methods, K_values, color=colors_bar, width=0.50,
                    edgecolor='white', linewidth=2, zorder=3)

# Value labels on top of bars
for bar, K, col in zip(bars, K_values, colors_bar):
    ax_right.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
                  f'$K \\approx {K}$', ha='center', fontsize=14, fontweight='bold',
                  color=col)

# Ratio annotation arrow between the two bars
ax_right.annotate('',
                  xy=(0.15, 0.025), xycoords='data',
                  xytext=(0.15, 0.78), textcoords='data',
                  arrowprops=dict(arrowstyle='<->', color='#333', lw=2.0,
                                  mutation_scale=15))
ax_right.text(0.15, 0.38, r'$\times\, 61$', fontsize=16, fontweight='bold',
              ha='center', color='#333',
              bbox=dict(boxstyle='round,pad=0.25', facecolor='#ffffdd',
                       edgecolor='#999', alpha=0.95))

# Info boxes below bars
ax_right.text(0, -0.13, 'Derivative-free',
              ha='center', fontsize=10, color='#2255aa', fontweight='bold')
ax_right.text(0, -0.20, r'2 evals of $S_p$ / step',
              ha='center', fontsize=9, color='#2255aa', style='italic')

ax_right.text(1, -0.13, "Requires $f'(x)$",
              ha='center', fontsize=10, color='#cc2222', fontweight='bold')
ax_right.text(1, -0.20, '1 eval / step',
              ha='center', fontsize=9, color='#cc2222', style='italic')

ax_right.set_ylabel(r'Quadratic constant $K$', fontsize=13)
ax_right.set_title('Quadratic constant comparison\n'
                    r'$\varepsilon_{n+1} \leq K \cdot \varepsilon_n^2$',
                    fontsize=13, fontweight='bold')
ax_right.set_ylim(-0.05, 1.0)
ax_right.grid(True, alpha=0.2, axis='y')
ax_right.spines['top'].set_visible(False)
ax_right.spines['right'].set_visible(False)

plt.tight_layout(pad=3.0)
plt.subplots_adjust(bottom=0.15)
plt.savefig('pandrosion_steffensen.pdf', dpi=200, bbox_inches='tight',
            facecolor='#fafaf8')
plt.savefig('pandrosion_steffensen.png', dpi=150, bbox_inches='tight',
            facecolor='#fafaf8')
print("Steffensen vs Newton figure saved.")
plt.close()