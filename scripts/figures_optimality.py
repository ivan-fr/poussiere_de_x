"""
Figures for: Kung-Traub Optimality of Steffensen-Pandrosion
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 11})

def S_p(s, p):
    if abs(s - 1) < 1e-15: return float(p)
    return (1 - s**p) / (1 - s)

def pandrosion(s, p, x):
    sp = S_p(s, p)
    if abs(sp) < 1e-30: return s
    return 1 - (x - 1) / (x * sp)

def steffensen_pandrosion(s, p, x):
    s0 = s
    s1 = pandrosion(s0, p, x)
    s2 = pandrosion(s1, p, x)
    d = s2 - 2*s1 + s0
    if abs(d) < 1e-30: return s2
    return s0 - (s1 - s0)**2 / d

def steffensen_generic(u, p, x):
    fu = u**p - x
    if abs(fu) < 1e-30: return u
    w = u + fu
    fw = w**p - x
    d = fw - fu
    if abs(d) < 1e-30: return u
    return u - fu**2 / d

def newton_step(u, p, x):
    return ((p-1)*u + x/u**(p-1)) / p

def secant_step(u0, u1, p, x):
    f0, f1 = u0**p - x, u1**p - x
    if abs(f1 - f0) < 1e-30: return u1
    return u1 - f1*(u1-u0)/(f1-f0)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#fafaf8')

# ── LEFT: K_SP/K_N ratio vs x ──
ax1.set_facecolor('#fafaf8')

for p, col, marker in [(2,'#2255aa','o'), (3,'#cc2222','s'), (5,'#22aa55','^'), (10,'#cc8800','D')]:
    xs = np.linspace(1.1, 12, 50)
    ratios = []
    for x in xs:
        alpha = x**(1/p)
        K_N = (p-1) / (2*alpha)
        s = 0.5; s_star = 1/alpha
        eps_list = []
        for n in range(10):
            eps_list.append(abs(s - s_star))
            if abs(s - s_star) < 1e-15: break
            s = steffensen_pandrosion(s, p, x)
        Ks = [eps_list[i+1]/eps_list[i]**2 for i in range(len(eps_list)-1) if eps_list[i] > 1e-6]
        if Ks:
            K_SP = np.mean(Ks[-2:]) if len(Ks) >= 2 else Ks[0]
            ratios.append(K_SP / K_N)
        else:
            ratios.append(np.nan)
    ax1.semilogy(xs, ratios, marker=marker, color=col, markersize=4, linewidth=1.8,
                 markevery=5, label=f'$p={p}$')

ax1.axhline(y=1, color='#999', linewidth=1.0, linestyle='--', alpha=0.6)
ax1.text(11, 1.15, '$K_{SP}=K_N$', fontsize=9, color='#999', style='italic')
ax1.fill_between([1, 13], [0.001, 0.001], [1, 1], alpha=0.06, color='green')
ax1.text(2, 0.005, 'Pandrosion advantage', fontsize=9, color='#228833', style='italic')

ax1.set_xlabel('$x$', fontsize=14)
ax1.set_ylabel('$K_{SP} / K_N$', fontsize=14)
ax1.set_title('Quadratic constant ratio\n'
              r'$K_{SP}/K_N < 1$ means Pandrosion wins',
              fontsize=12, fontweight='bold')
ax1.legend(fontsize=11, title='Order $p$')
ax1.set_xlim(1, 12.5)
ax1.set_ylim(0.005, 50)
ax1.grid(True, alpha=0.2, which='both')

# ── RIGHT: Convergence traces ──
ax2.set_facecolor('#fafaf8')
p, x = 3, 2.0
alpha = x**(1/p)

# Newton
u = 0.5
eps_newton = []
for n in range(8):
    eps_newton.append(abs(u - alpha))
    if abs(u - alpha) < 1e-16: break
    u = newton_step(u, p, x)

# Steffensen-Pandrosion
s = 0.5
eps_sp = []
for n in range(6):
    v = x * s**(p-1)
    eps_sp.append(abs(v - alpha))
    if abs(v - alpha) < 1e-16: break
    s = steffensen_pandrosion(s, p, x)

# Generic Steffensen
u = 1.5
eps_gs = []
for n in range(12):
    eps_gs.append(abs(u - alpha))
    if abs(u - alpha) < 1e-16: break
    try:
        u_new = steffensen_generic(u, p, x)
        if abs(u_new) > 1e6: break
        u = u_new
    except: break

# Secant
u0, u1 = 0.5, 1.5
eps_sec = [abs(u1 - alpha)]
for n in range(15):
    u_new = secant_step(u0, u1, p, x)
    eps_sec.append(abs(u_new - alpha))
    if abs(u_new - alpha) < 1e-16: break
    u0, u1 = u1, u_new

ax2.semilogy(range(len(eps_sp)), eps_sp, 'D-', color='#cc2222', markersize=9,
             linewidth=2.5, markeredgecolor='white', markeredgewidth=0.8,
             label='Steffensen-Pandrosion', zorder=5)
ax2.semilogy(range(len(eps_newton)), eps_newton, 's-', color='#2255aa', markersize=8,
             linewidth=2.0, markeredgecolor='white', markeredgewidth=0.8,
             label='Newton')
ax2.semilogy(range(len(eps_gs)), eps_gs, '^-', color='#888888', markersize=8,
             linewidth=1.5, markeredgecolor='white', markeredgewidth=0.8,
             label='Generic Steffensen')
ax2.semilogy(range(len(eps_sec)), eps_sec, 'o-', color='#22aa55', markersize=6,
             linewidth=1.5, markeredgecolor='white', markeredgewidth=0.8,
             label='Secant', alpha=0.7)

ax2.axhline(y=2.2e-16, color='#999', linewidth=0.8, linestyle=':', alpha=0.6)
ax2.text(6, 5e-16, 'Machine $\\varepsilon$', fontsize=9, color='#999', style='italic')

ax2.set_xlabel('Step $n$', fontsize=14)
ax2.set_ylabel('$|u_n - \\alpha|$', fontsize=14)
ax2.set_title(r'Convergence comparison for $\sqrt[3]{2}$' + '\n'
              r'$p=3$, $x=2$', fontsize=12, fontweight='bold')
ax2.legend(fontsize=10, loc='upper right')
ax2.set_ylim(1e-17, 5)
ax2.set_xlim(-0.3, 10)
ax2.grid(True, alpha=0.2, which='both')

plt.tight_layout(pad=2.0)
plt.savefig('pandrosion_optimality_figures.png', dpi=200, bbox_inches='tight', facecolor='#fafaf8')
plt.savefig('pandrosion_optimality_figures.pdf', dpi=200, bbox_inches='tight', facecolor='#fafaf8')
print("Saved pandrosion_optimality_figures.{png,pdf}")
plt.close()
