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
# FIGURE 1  — construction géométrique de Pandrosion (p=3, x=2)
# ═══════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(9, 8.5))
fig.patch.set_facecolor('#f5f2eb')
ax.set_facecolor('#f5f2eb')

W, H = 5.0, 5.0
p, x_val = 3, 2

# --- itérations
s = 0.30
s_list = [s]
for _ in range(3):
    s = pandrosion_s(s, p, x_val)
    s_list.append(s)
s_star = x_val**(-1/p)

# --- rectangle avec léger gradient
rect = mpatches.FancyBboxPatch((0, 0), W, H,
    boxstyle='square,pad=0', linewidth=2.0,
    edgecolor='#333333', facecolor='#ede8dc', zorder=1)
ax.add_patch(rect)

# --- diagonale principale (Thalès) : (0,H) → (W,0)
ax.plot([0, W], [H, 0], color='#444444', linewidth=2.5, zorder=4)

# label diagonale
mid_diag_x, mid_diag_y = W * 0.25, H * 0.75
ax.text(mid_diag_x - 0.55, mid_diag_y + 0.20,
        r'Diagonale — Thalès',
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
        r'Courbe $v = x\,s^{p-1}$',
        fontsize=9.5, color='#8855aa', ha='center', va='bottom',
        style='italic')

# --- palette bleue
N = len(s_list)
blue_shades = [plt.cm.Blues(0.30 + 0.60 * i / (N - 1)) for i in range(N)]

for i, s_n in enumerate(s_list):
    u_n = H * (1 - s_n)
    col = blue_shades[i]

    # Point d'intersection avec la diagonale: diagonale va de (0,H) à (W,0)
    # y = H - (H/W)*x, donc à y=u_n: x_diag = W*(H-u_n)/H = W*s_n
    x_diag = W * s_n

    # Point d'intersection avec la courbe parabolique:
    # y_curve = H*(1-s_n) = u_n, x_curve = W*s_n^(p-1)
    x_para = W * s_n**(p-1)

    # horizontale entre les deux points d'intersection
    x_left = min(x_diag, x_para)
    x_right = max(x_diag, x_para)
    ax.hlines(u_n, 0, W, colors=col, linewidth=1.2, alpha=0.35, zorder=2,
              linestyle=':')
    # Segment solide entre diagonale et courbe
    ax.plot([x_left, x_right], [u_n, u_n], color=col, linewidth=2.2,
            alpha=0.85, zorder=3)

    # point sur la diagonale
    ax.plot(x_diag, u_n, 'o', color=col, markersize=8, zorder=6,
            markeredgecolor='white', markeredgewidth=0.8)

    # point sur la courbe parabolique
    ax.plot(x_para, u_n, 's', color=col, markersize=7, zorder=6,
            markeredgecolor='white', markeredgewidth=0.8)

    # label s_n sur la gauche : segment annoté
    ax.text(-0.35, u_n, f'$s_{i}$',
            fontsize=11, va='center', ha='right', color=col)

    # label u_n à droite
    y_off = 0
    if i == 3:
        y_off = -0.22
    ax.text(W + 0.15, u_n + y_off, f'$u_{i}$',
            fontsize=11, va='center', color=col)

    # Flèche de projection vers le bas : du point parabolique vers l'axe x
    v_n = x_val * s_n**(p-1)
    v_x = v_n * (W / x_val)
    ax.annotate('', xy=(v_x, -0.08), xytext=(x_para, u_n),
                arrowprops=dict(arrowstyle='->', color=col,
                                lw=1.0, alpha=0.45,
                                connectionstyle='arc3,rad=0.0',
                                linestyle='dashed'))

    # Label v_n sous le rectangle
    ax.text(v_x, -0.30 - 0.28 * i,
            f'$v_{i}={v_n:.3f}$', ha='center', fontsize=9, color=col)

# --- ligne du point fixe u*
u_star = H * (1 - s_star)
ax.hlines(u_star, 0, W, colors='#bb2222', linewidth=2.2,
          linestyle='--', zorder=5)
ax.text(W + 0.15, u_star, r'$u^*$',
        fontsize=13, va='center', color='#bb2222', fontweight='bold')

# Annotations s* et 1-s* sur le bord gauche
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

# Point fixe: point sur la diagonale et sur la courbe
x_diag_star = W * s_star
x_para_star = W * s_star**(p-1)
ax.plot(x_diag_star, u_star, 'o', color='#bb2222', markersize=9, zorder=7,
        markeredgecolor='white', markeredgewidth=1.0)
ax.plot(x_para_star, u_star, 's', color='#bb2222', markersize=8, zorder=7,
        markeredgecolor='white', markeredgewidth=1.0)

# v* sous le rectangle
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

# --- axe de lecture gradué (sous le rectangle)
y_axis = -0.08
ax.plot([0, W], [y_axis, y_axis], 'k-', linewidth=1.5, zorder=2)
for v_tick in [0, 0.5, 1.0, 1.5, 2.0]:
    xt = v_tick * (W / x_val)
    ax.plot(xt, y_axis, 'k|', markersize=7, zorder=3)
    ax.text(xt, y_axis - 0.20, f'{v_tick}', ha='center', fontsize=9.5, color='#444')

# tick spécial pour cbrt(2)
xt_star = v_star * (W / x_val)
ax.plot(xt_star, y_axis, '|', color='#bb2222', markersize=9, zorder=4)
ax.text(xt_star, y_axis - 0.20, r'$\sqrt[3]{2}$', ha='center', fontsize=10,
        color='#bb2222')

# légende manuelle
legend_items = [
    plt.Line2D([0], [0], color='#444', linewidth=2.2, label='Diagonale (Thalès)'),
    plt.Line2D([0], [0], color='#8855aa', linewidth=2.0, label=r'Courbe $v=x\,s^{p-1}$'),
    plt.Line2D([0], [0], color='#bb2222', linewidth=2.0, linestyle='--', label=r'Point fixe $u^*$'),
    plt.Line2D([0], [0], marker='o', color='#4477bb', markersize=7,
               linestyle='', label=r'Intersection diagonale'),
    plt.Line2D([0], [0], marker='s', color='#4477bb', markersize=6,
               linestyle='', label=r'Intersection courbe'),
]
ax.legend(handles=legend_items, loc='upper right', fontsize=9,
          framealpha=0.85, bbox_to_anchor=(1.0, 1.0))

# --- axes labels
ax.text(W / 2, -0.60 - 0.28 * N,
        r'Sortie $v_n = x\,s_n^{p-1}$  (converge vers $\sqrt[3]{2}$)',
        ha='center', fontsize=11, color='#333', style='italic')

# Dimensions
ax.set_xlim(-0.85, W + 1.05)
ax.set_ylim(-0.60 - 0.28 * (N + 1), H + 0.45)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title(
    r'Construction de Pandrosion  ($p=3$, $x=2$)' + '\n'
    r'Parallèles successives — convergence vers $v^*=\sqrt[3]{2}$',
    fontsize=13, pad=14, fontweight='bold')

plt.tight_layout()
plt.savefig('pandrosion_geometry.pdf', dpi=200, bbox_inches='tight',
            facecolor='#f5f2eb')
plt.savefig('pandrosion_geometry.png', dpi=200, bbox_inches='tight',
            facecolor='#f5f2eb')
print("Figure géométrique sauvegardée.")
plt.close()


# ═══════════════════════════════════════════════════════════════════
# FIGURE 2  — cobweb + profils de poussière
# ═══════════════════════════════════════════════════════════════════

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
fig2.patch.set_facecolor('#fafaf8')

# ── gauche : cobweb ──────────────────────────────────────────────
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
ax1.set_title('Diagramme cobweb — convergence de Pandrosion\n'
              r'($p=3,\ x=2$,   $s_0=0{,}5$)', fontsize=11)
ax1.legend(fontsize=9.5, loc='upper left')
ax1.grid(True, alpha=0.25)

# ── droite : profils de poussière ───────────────────────────────
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

ax2.set_xlabel(r'Étape $n$', fontsize=13)
ax2.set_ylabel(r'$\varepsilon_n = |v_n - \alpha|$', fontsize=13)
ax2.set_title('Profils de poussière : Pandrosion vs Newton', fontsize=11)
ax2.legend(fontsize=9.5)
ax2.grid(True, alpha=0.25)

plt.tight_layout(pad=2.2)
plt.savefig('pandrosion_figures.pdf', dpi=200, bbox_inches='tight',
            facecolor='#fafaf8')
plt.savefig('pandrosion_figures.png', dpi=150, bbox_inches='tight',
            facecolor='#fafaf8')
print("Figure cobweb + poussière sauvegardée.")
plt.close()