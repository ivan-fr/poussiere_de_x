"""
Figures for: Pandrosion Iteration in the Complex Plane
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.cm as cm
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'font.size': 11,
    'text.usetex': False,
})

def S_p(s, p):
    if abs(s - 1) < 1e-15:
        return complex(p)
    return (1 - s**p) / (1 - s)

def pandrosion_complex(s, p, x):
    sp = S_p(s, p)
    if abs(sp) < 1e-30:
        return s
    return 1 - (x - 1) / (x * sp)

def lambda_complex(p, x):
    alpha = x**(1/p)
    num = (alpha - 1) * sum(k * alpha**(p-1-k) for k in range(1, p))
    den = sum(alpha**k for k in range(p))
    return num / den


# ═══════════════════════════════════════════════════════════════════
# FIGURE 1 — Basins of attraction (high resolution)
# ═══════════════════════════════════════════════════════════════════

print("Generating Figure 1: Basins of attraction...")

p_basin, x_basin = 3, 2+1j
roots = [x_basin**(1/p_basin) * np.exp(2j*np.pi*k/p_basin) for k in range(p_basin)]
fixed_points = [1/r for r in roots]

resolution = 500
re_range = np.linspace(-1.8, 2.5, resolution)
im_range = np.linspace(-2.2, 2.2, resolution)
basin = np.zeros((resolution, resolution), dtype=int)
speed = np.zeros((resolution, resolution))

for i, re in enumerate(re_range):
    for j, im in enumerate(im_range):
        s = complex(re, im)
        for n in range(150):
            try:
                s_new = pandrosion_complex(s, p_basin, x_basin)
                if abs(s_new) > 1e8:
                    basin[j, i] = -1
                    speed[j, i] = n
                    break
                for k, fp in enumerate(fixed_points):
                    if abs(s_new - fp) < 1e-6:
                        basin[j, i] = k + 1
                        speed[j, i] = n
                        break
                if basin[j, i] != 0:
                    break
                s = s_new
            except:
                basin[j, i] = -1
                speed[j, i] = n
                break
        if basin[j, i] == 0:
            speed[j, i] = 150

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.5))
fig.patch.set_facecolor('#0d1117')

# Basin colors
basin_colors = np.zeros((*basin.shape, 4))
color_map = {
    -1: (0.15, 0.15, 0.15, 1.0),   # diverged: dark grey
    0:  (0.25, 0.25, 0.25, 1.0),    # no convergence: grey
    1:  (0.13, 0.33, 0.67, 1.0),    # blue
    2:  (0.80, 0.13, 0.13, 1.0),    # red
    3:  (0.13, 0.67, 0.33, 1.0),    # green
}
# Add shading by convergence speed
for i in range(resolution):
    for j in range(resolution):
        base = color_map.get(basin[i, j], (0.5, 0.5, 0.5, 1.0))
        # Darken slow-converging pixels
        brightness = max(0.3, 1.0 - speed[i, j] / 80)
        basin_colors[i, j] = (base[0]*brightness, base[1]*brightness, base[2]*brightness, 1.0)

ax1.set_facecolor('#0d1117')
ax1.imshow(basin_colors, extent=[re_range[0], re_range[-1], im_range[0], im_range[-1]],
           origin='lower', aspect='auto')

# Mark fixed points
fp_colors = ['#4488ff', '#ff4444', '#44ff88']
fp_labels = [r'$s_0^*$', r'$s_1^*$', r'$s_2^*$']
for k, (fp, col, lbl) in enumerate(zip(fixed_points, fp_colors, fp_labels)):
    ax1.plot(fp.real, fp.imag, '*', color=col, markersize=18,
             markeredgecolor='white', markeredgewidth=1.5, zorder=10)
    ax1.text(fp.real + 0.12, fp.imag + 0.12, lbl, color='white',
             fontsize=14, fontweight='bold', zorder=11)

ax1.set_xlabel(r'Re($s$)', fontsize=13, color='white')
ax1.set_ylabel(r'Im($s$)', fontsize=13, color='white')
ax1.set_title(f"Basins of attraction\n$p={p_basin}$, $x={x_basin}$",
              fontsize=13, fontweight='bold', color='white')
ax1.tick_params(colors='white')
for spine in ax1.spines.values():
    spine.set_color('#444')

# Convergence speed
ax2.set_facecolor('#0d1117')
im2 = ax2.imshow(speed, extent=[re_range[0], re_range[-1], im_range[0], im_range[-1]],
                  origin='lower', cmap='inferno_r', aspect='auto', vmax=60)
for k, (fp, col) in enumerate(zip(fixed_points, fp_colors)):
    ax2.plot(fp.real, fp.imag, '*', color='cyan', markersize=18,
             markeredgecolor='white', markeredgewidth=1.5, zorder=10)
ax2.set_xlabel(r'Re($s$)', fontsize=13, color='white')
ax2.set_ylabel(r'Im($s$)', fontsize=13, color='white')
ax2.set_title(f"Convergence speed (iterations)\n$p={p_basin}$, $x={x_basin}$",
              fontsize=13, fontweight='bold', color='white')
ax2.tick_params(colors='white')
for spine in ax2.spines.values():
    spine.set_color('#444')
cbar = plt.colorbar(im2, ax=ax2)
cbar.set_label('Iterations', color='white')
cbar.ax.tick_params(colors='white')

plt.tight_layout(pad=2.0)
plt.savefig('pandrosion_complex_basins.png', dpi=200, bbox_inches='tight',
            facecolor='#0d1117')
plt.savefig('pandrosion_complex_basins.pdf', dpi=200, bbox_inches='tight',
            facecolor='#0d1117')
print("  → Saved pandrosion_complex_basins.{png,pdf}")
plt.close()


# ═══════════════════════════════════════════════════════════════════
# FIGURE 2 — |λ| heatmap in complex plane
# ═══════════════════════════════════════════════════════════════════

print("Generating Figure 2: |λ| heatmap...")

p_hm = 3
res_hm = 400
re_hm = np.linspace(-4, 6, res_hm)
im_hm = np.linspace(-4, 4, res_hm)
lambda_grid = np.zeros((res_hm, res_hm))

for i, re in enumerate(re_hm):
    for j, im in enumerate(im_hm):
        x = complex(re, im)
        if abs(x) < 0.05:
            lambda_grid[j, i] = np.nan
            continue
        try:
            lam = lambda_complex(p_hm, x)
            lambda_grid[j, i] = min(abs(lam), 2.0)
        except:
            lambda_grid[j, i] = np.nan

fig, ax = plt.subplots(figsize=(10, 7))
fig.patch.set_facecolor('#0d1117')
ax.set_facecolor('#0d1117')

im = ax.imshow(lambda_grid, extent=[re_hm[0], re_hm[-1], im_hm[0], im_hm[-1]],
               origin='lower', cmap='RdYlBu_r', aspect='auto', vmin=0, vmax=1.5)

# Contour at |λ| = 1 (convergence boundary)
CS = ax.contour(re_hm, im_hm, lambda_grid, levels=[1.0],
                colors=['white'], linewidths=2.5, linestyles='--')
ax.clabel(CS, fmt=r'$|\lambda|=1$', fontsize=11, colors='white')

# Mark some reference points
ref_points = [(2, 0, r'$x=2$'), (0, 1, r'$x=i$'), (-1, 0, r'$x=-1$'), (1, 0, r'$x=1$')]
for rx, iy, lbl in ref_points:
    ax.plot(rx, iy, 'o', color='white', markersize=6, markeredgecolor='black',
            markeredgewidth=0.5, zorder=10)
    ax.text(rx + 0.15, iy + 0.15, lbl, color='white', fontsize=10, zorder=11)

ax.set_xlabel(r'Re($x$)', fontsize=14, color='white')
ax.set_ylabel(r'Im($x$)', fontsize=14, color='white')
ax.set_title(r"Contraction ratio $|\lambda_{3,x}|$ in the complex plane" + "\n"
             r"Convergence region: $|\lambda| < 1$ (blue)",
             fontsize=13, fontweight='bold', color='white')
ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_color('#444')

cbar = plt.colorbar(im, ax=ax, shrink=0.85)
cbar.set_label(r'$|\lambda_{3,x}|$', fontsize=13, color='white')
cbar.ax.tick_params(colors='white')

plt.tight_layout()
plt.savefig('pandrosion_lambda_heatmap.png', dpi=200, bbox_inches='tight',
            facecolor='#0d1117')
plt.savefig('pandrosion_lambda_heatmap.pdf', dpi=200, bbox_inches='tight',
            facecolor='#0d1117')
print("  → Saved pandrosion_lambda_heatmap.{png,pdf}")
plt.close()


# ═══════════════════════════════════════════════════════════════════
# FIGURE 3 — Steffensen convergence in ℂ
# ═══════════════════════════════════════════════════════════════════

print("Generating Figure 3: Steffensen in ℂ...")

def steffensen_step(s, p, x):
    s0 = s
    s1 = pandrosion_complex(s0, p, x)
    s2 = pandrosion_complex(s1, p, x)
    denom = s2 - 2*s1 + s0
    if abs(denom) < 1e-30:
        return s2
    return s0 - (s1 - s0)**2 / denom

fig, ax = plt.subplots(figsize=(10, 6.5))
fig.patch.set_facecolor('#fafaf8')
ax.set_facecolor('#fafaf8')

test_cases = [
    (3, 2+0j, '#2255aa', r'$x=2$ (real)', 'D'),
    (3, 2+1j, '#cc2222', r'$x=2+i$', 's'),
    (3, 1j, '#22aa55', r'$x=i$', '^'),
    (3, 3+2j, '#cc8800', r'$x=3+2i$', 'o'),
    (4, 1+3j, '#8833bb', r'$x=1+3i$, $p=4$', 'v'),
]

for p, x, col, lbl, marker in test_cases:
    alpha = x**(1/p)
    s = 0.5 + 0j
    eps_list = []
    for n in range(7):
        v = x * s**(p-1)
        eps = abs(v - alpha)
        if eps < 1e-17:
            break
        eps_list.append(eps)
        s = steffensen_step(s, p, x)
    # Final point
    v = x * s**(p-1)
    eps = abs(v - alpha)
    if eps > 1e-17:
        eps_list.append(eps)
    
    ns = list(range(len(eps_list)))
    ax.semilogy(ns, eps_list, marker=marker, color=col, markersize=9,
                linewidth=2.2, markeredgecolor='white', markeredgewidth=0.8,
                label=lbl, zorder=4)

ax.axhline(y=2.2e-16, color='#999', linewidth=1.0, linestyle=':', alpha=0.7)
ax.text(5.5, 4e-16, 'Machine precision', fontsize=9, color='#999', style='italic')

ax.set_xlabel(r'Steffensen step $n$', fontsize=13)
ax.set_ylabel(r'$|v_n - \alpha|$', fontsize=13)
ax.set_title(r'Steffensen--Pandrosion convergence in $\mathbb{C}$' + '\n'
             r'Quadratic convergence for various complex targets',
             fontsize=13, fontweight='bold')
ax.legend(fontsize=11, loc='upper right')
ax.set_ylim(1e-17, 5)
ax.set_xlim(-0.3, 6.5)
ax.grid(True, alpha=0.25, which='both')

plt.tight_layout()
plt.savefig('pandrosion_steffensen_complex.png', dpi=200, bbox_inches='tight',
            facecolor='#fafaf8')
plt.savefig('pandrosion_steffensen_complex.pdf', dpi=200, bbox_inches='tight',
            facecolor='#fafaf8')
print("  → Saved pandrosion_steffensen_complex.{png,pdf}")
plt.close()


# ═══════════════════════════════════════════════════════════════════
# FIGURE 4 — Euler product approximation
# ═══════════════════════════════════════════════════════════════════

print("Generating Figure 4: Euler product...")

def sieve(n):
    is_prime = [True] * (n+1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5)+1):
        if is_prime[i]:
            for j in range(i*i, n+1, i):
                is_prime[j] = False
    return [i for i in range(2, n+1) if is_prime[i]]

primes = sieve(1000)

def zeta_direct(s, terms=50000):
    return sum(n**(-s) for n in range(1, terms+1))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#fafaf8')

# Left: partial Euler product convergence to ζ(s) for s=2
ax1.set_facecolor('#fafaf8')
zeta_2 = np.pi**2 / 6  # exact value

partial = 1.0
partials = []
ns_primes = []
for i, p in enumerate(primes[:80]):
    partial *= 1 / (1 - p**(-2))
    partials.append(partial)
    ns_primes.append(i + 1)

errors = [abs(z - zeta_2) for z in partials]
ax1.semilogy(ns_primes, errors, 'o-', color='#2255aa', markersize=4,
             linewidth=1.8, label=r'$|\Pi_N - \zeta(2)|$')
ax1.set_xlabel('Number of primes $N$', fontsize=13)
ax1.set_ylabel(r'$|\Pi_N(2) - \pi^2/6|$', fontsize=13)
ax1.set_title(r"Partial Euler product $\to \zeta(2) = \pi^2/6$" + "\n"
              r"Each factor is a Pandrosion sum $S_\infty(p^{-2})$",
              fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.25)
ax1.legend(fontsize=11)

# Right: Pandrosion ratios on the critical line (unit circle diagram)
ax2.set_facecolor('#fafaf8')
ax2.set_aspect('equal')

# Unit circle
theta = np.linspace(0, 2*np.pi, 200)
ax2.plot(np.cos(theta), np.sin(theta), '--', color='#999', linewidth=1.5, alpha=0.5,
         label=r'$|z|=1$')
ax2.plot(0, 0, '+', color='#333', markersize=12, markeredgewidth=2)

# Plot ratios p^{-s} for s = 1/2 + 14.135i (first zeta zero)
t = 14.134725
colors_cr = plt.cm.viridis(np.linspace(0.2, 0.9, 10))
for idx, p in enumerate(primes[:10]):
    s_val = 0.5 + t*1j
    ratio = p**(-s_val)
    ax2.plot(ratio.real, ratio.imag, 'o', color=colors_cr[idx], markersize=12,
             markeredgecolor='white', markeredgewidth=1.0, zorder=5)
    # Circle of radius p^{-1/2}
    r = p**(-0.5)
    ax2.plot(r*np.cos(theta), r*np.sin(theta), '-', color=colors_cr[idx],
             linewidth=0.8, alpha=0.3)
    ax2.text(ratio.real + 0.03, ratio.imag + 0.04, f'$p={p}$',
             fontsize=8, color=colors_cr[idx])

ax2.set_xlabel(r'Re($p^{-s}$)', fontsize=13)
ax2.set_ylabel(r'Im($p^{-s}$)', fontsize=13)
ax2.set_title(r"Pandrosion ratios $p^{-s}$ on critical line" + "\n"
              r"$s = \frac{1}{2} + 14.135i$ (first $\zeta$ zero)",
              fontsize=12, fontweight='bold')
ax2.set_xlim(-0.85, 0.85)
ax2.set_ylim(-0.85, 0.85)
ax2.grid(True, alpha=0.2)
ax2.legend(fontsize=10, loc='upper right')

plt.tight_layout(pad=2.5)
plt.savefig('pandrosion_euler.png', dpi=200, bbox_inches='tight',
            facecolor='#fafaf8')
plt.savefig('pandrosion_euler.pdf', dpi=200, bbox_inches='tight',
            facecolor='#fafaf8')
print("  → Saved pandrosion_euler.{png,pdf}")
plt.close()

print("\nAll figures generated successfully!")
