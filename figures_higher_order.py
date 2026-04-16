"""Figures for Paper IV: Higher-Order Pandrosion Methods"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 11})

def S_p(s, p):
    if abs(s - 1) < 1e-15: return float(p)
    return (1 - s**p) / (1 - s)

def h(s, p, x):
    sp = S_p(s, p)
    if abs(sp) < 1e-30: return s
    return 1 - (x - 1) / (x * sp)

def T2(s, p, x):
    s0, s1, s2 = s, h(s, p, x), h(h(s, p, x), p, x)
    d = s2 - 2*s1 + s0
    if abs(d) < 1e-30: return s2
    lam_hat = (s2 - s1) / (s1 - s0) if abs(s1 - s0) > 1e-30 else 0
    return s0 - (s1 - s0)**2 / d, lam_hat

def T4(s, p, x):
    t2, lam_hat = T2(s, p, x)
    s3 = h(t2, p, x)
    if abs(lam_hat - 1) < 1e-30: return s3
    return t2 - (s3 - t2) / (lam_hat - 1)

def T8(s, p, x):
    t2, lam_hat = T2(s, p, x)
    s3 = h(t2, p, x)
    if abs(lam_hat - 1) < 1e-30: return s3
    t4 = t2 - (s3 - t2) / (lam_hat - 1)
    s4 = h(t4, p, x)
    if abs(lam_hat - 1) < 1e-30: return s4
    return t4 - (s4 - t4) / (lam_hat - 1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#fafaf8')

# ── LEFT: Convergence traces ──
ax1.set_facecolor('#fafaf8')
p, x = 3, 2.0
alpha = x**(1/p)
s_star = 1/alpha

for name, method, color, marker, ms in [
    ('$T_2$ (order 2)', lambda s: T2(s,p,x)[0], '#2255aa', 's', 8),
    ('$T_4$ (order 4)', lambda s: T4(s,p,x), '#cc2222', 'D', 9),
    ('$T_8$ (order 8)', lambda s: T8(s,p,x), '#22aa55', '^', 9),
]:
    s = 0.5
    eps_list = []
    for n in range(6):
        eps_list.append(abs(s - s_star))
        if abs(s - s_star) < 1e-16: break
        s = method(s)
        if abs(s) > 1e10: break
    eps_plot = [e if e > 1e-17 else 1e-17 for e in eps_list]
    ax1.semilogy(range(len(eps_plot)), eps_plot, f'{marker}-', color=color,
                 markersize=ms, linewidth=2.2, markeredgecolor='white',
                 markeredgewidth=0.8, label=name, zorder=5)

ax1.axhline(y=2.2e-16, color='#999', linewidth=0.8, linestyle=':', alpha=0.6)
ax1.text(3.5, 5e-16, 'Machine $\\varepsilon$', fontsize=9, color='#999', style='italic')
ax1.set_xlabel('Step $n$', fontsize=14)
ax1.set_ylabel('$|s_n - s^*|$', fontsize=14)
ax1.set_title(r'Convergence traces for $\sqrt[3]{2}$' + '\n'
              '$p=3$, $x=2$, $s_0=0.5$', fontsize=12, fontweight='bold')
ax1.legend(fontsize=11)
ax1.set_ylim(1e-18, 2)
ax1.set_xlim(-0.2, 5.2)
ax1.grid(True, alpha=0.2, which='both')

# ── RIGHT: Efficiency index ──
ax2.set_facecolor('#fafaf8')
ns = np.arange(1, 11)
Es = 2**((ns-1)/ns)

ax2.plot(ns, Es, 'o-', color='#cc2222', markersize=10, linewidth=2.5,
         markeredgecolor='white', markeredgewidth=1.0, zorder=5)
ax2.axhline(y=2, color='#999', linewidth=1.0, linestyle='--', alpha=0.5)
ax2.text(8.5, 2.03, '$E = 2$', fontsize=10, color='#999', style='italic')

# Annotate key levels
annotations = [(1,'$h$\n(linear)'), (2,'$T_2$\n(Steffensen)'),
               (3,'$T_4$'), (4,'$T_8$')]
for n_val, label in annotations:
    E_val = 2**((n_val-1)/n_val)
    ax2.annotate(label, xy=(n_val, E_val), xytext=(n_val+0.35, E_val-0.08),
                fontsize=9, color='#333', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#666', lw=0.8))

# Shade optimal region
ax2.fill_between([2.5, 3.5], [1, 1], [2.1, 2.1], alpha=0.08, color='#cc2222')
ax2.text(2.65, 1.05, 'Best\ntradeoff', fontsize=8, color='#cc2222', style='italic')

# Marginal gain arrows
ax2.annotate('', xy=(3, 2**((3-1)/3)), xytext=(2, 2**((2-1)/2)),
             arrowprops=dict(arrowstyle='->', color='#22aa55', lw=2))
ax2.text(2.2, 1.47, '$\\Delta E = 0.17$', fontsize=9, color='#22aa55', fontweight='bold')

ax2.set_xlabel('Evaluations per step $n$', fontsize=14)
ax2.set_ylabel('Efficiency index $E = (2^{n-1})^{1/n}$', fontsize=14)
ax2.set_title('Kung--Traub efficiency hierarchy\n'
              'Pandrosion saturates every level', fontsize=12, fontweight='bold')
ax2.set_xlim(0.5, 10.5)
ax2.set_ylim(0.9, 2.15)
ax2.set_xticks(range(1, 11))
ax2.grid(True, alpha=0.2)

plt.tight_layout(pad=2.0)
plt.savefig('pandrosion_higher_order.png', dpi=200, bbox_inches='tight', facecolor='#fafaf8')
plt.savefig('pandrosion_higher_order.pdf', dpi=200, bbox_inches='tight', facecolor='#fafaf8')
print("Saved pandrosion_higher_order.{png,pdf}")
plt.close()
