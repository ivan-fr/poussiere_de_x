"""
Figures for Paper V: Analog Pandrosion (final clean)
Adaptive MC samples + theoretical bias lines to avoid artifacts.
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 11})
np.random.seed(42)

def S_p_noisy(s, p, sigma):
    result = 1.0; power = s
    for k in range(1, p):
        if k > 1: power = power * s * (1 + np.random.normal(0, sigma))
        result += power * (1 + np.random.normal(0, sigma))
    return result * (1 + np.random.normal(0, sigma))

def h_noisy(s, p, x, sigma):
    sp = S_p_noisy(s, p, sigma)
    if abs(sp) < 1e-30: return s
    num = (x - 1) * (1 + np.random.normal(0, sigma))
    den = x * sp * (1 + np.random.normal(0, sigma))
    return (1 - num / den) * (1 + np.random.normal(0, sigma))

def h_iter_noisy(s, p, x, sigma, n):
    for _ in range(n):
        s = h_noisy(s, p, x, sigma)
    return s

def newton_1step(u, p, x, sigma):
    up1 = u**(p-1) * (1 + np.random.normal(0, sigma))
    fp = p * up1 * (1 + np.random.normal(0, sigma))
    fu = (u * up1 - x) * (1 + np.random.normal(0, sigma))
    if abs(fp) < 1e-30: return u
    return (u - fu / fp) * (1 + np.random.normal(0, sigma))

# Noiseless for theoretical lines
def S_p(s,p): return (1-s**p)/(1-s) if abs(s-1)>1e-15 else float(p)
def h(s,p,x): return 1-(x-1)/(x*S_p(s,p))
def h_iter(s,p,x,n):
    for _ in range(n): s = h(s,p,x)
    return s

p, x = 3, 2.0
alpha = x**(1/p)
s0_opt = 1 - (x-1)/(x*p)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#fafaf8')

sigmas = np.logspace(-4.5, -0.7, 20)

configs = [
    ('$h^1$  (10 ops)', 1, '#ff8800', '^', 1.6, 5),
    ('$h^2$  (15 ops)', 2, '#22aa55', 'o', 2.4, 6),
    ('$h^3$  (20 ops)', 3, '#cc2222', 'D', 2.0, 5),
    ("Newton (6 ops, needs $f'$)", 0, '#2255aa', 's', 2.0, 6),
]

# Theoretical (noise-free) biases
theo_biases = {}
for _, n_h, _, _, _, _ in configs:
    if n_h > 0:
        s_n = h_iter(s0_opt, p, x, n_h)
        theo_biases[n_h] = abs(x * s_n**(p-1) - alpha)
    else:
        u1 = ((p-1)*1.0 + x/1.0**(p-1))/p
        theo_biases[0] = abs(u1 - alpha)

# Collect MC data with adaptive N
data = {}
for label, n_h, color, marker, lw, ms in configs:
    biases, bits_list = [], []
    for sig in sigmas:
        # More samples at high noise to stabilize bias estimation
        N = max(5000, int(15000 * (sig / 1e-3)**0.5))
        N = min(N, 50000)
        if n_h > 0:
            out = np.array([x * h_iter_noisy(s0_opt, p, x, sig, n_h)**(p-1)
                           for _ in range(N)])
        else:
            out = np.array([newton_1step(1.0, p, x, sig) for _ in range(N)])
        std = np.std(out)
        biases.append(abs(np.mean(out) - alpha))
        bits_list.append(-np.log2(std/alpha) if std > 0 else 0)
    data[n_h] = (biases, bits_list)

# ── LEFT: Bias ──
ax1.set_facecolor('#fafaf8')

for label, n_h, color, marker, lw, ms in configs:
    biases, _ = data[n_h]
    zorder = 7 if n_h == 2 else (6 if n_h == 3 else (5 if n_h == 1 else 4))
    ax1.loglog(sigmas, biases, f'{marker}-', color=color, markersize=ms,
               linewidth=lw, markeredgecolor='white', markeredgewidth=0.8,
               label=label, zorder=zorder)
    # Theoretical bias as dashed line
    ax1.axhline(y=theo_biases[n_h], color=color, linewidth=0.8,
                linestyle='--', alpha=0.4, zorder=2)

# Shade commercial zone
ax1.fill_between([1e-4, 1e-3], [3e-7, 3e-7], [3e-1, 3e-1], alpha=0.06, color='#22aa55')
ax1.text(1.3e-4, 5e-7, 'Commercial\nanalog ICs', fontsize=8, color='#228833', style='italic')

# Bias advantage annotation
ax1.annotate('', xy=(2e-4, 6e-3), xytext=(2e-4, 7e-2),
            arrowprops=dict(arrowstyle='<->', color='#22aa55', lw=2.2))
ax1.text(2.6e-4, 2.3e-2, '$12\\times$', fontsize=15, color='#22aa55', fontweight='bold')

ax1.set_xlabel('Per-component noise $\\sigma$', fontsize=13)
ax1.set_ylabel('Systematic bias $|E[\\hat{\\alpha}] - \\alpha|$', fontsize=13)
ax1.set_title('Bias: Pandrosion $h^n$ beats Newton\n'
              'dashed lines = noise-free math residual',
              fontsize=12, fontweight='bold')
ax1.legend(fontsize=9.5, loc='lower right',
           framealpha=0.9, edgecolor='#ddd')
ax1.set_ylim(3e-7, 3e-1)
ax1.set_xlim(2e-5, 3e-1)
ax1.grid(True, alpha=0.15, which='both')

# ── RIGHT: Bits ──
ax2.set_facecolor('#fafaf8')

for label, n_h, color, marker, lw, ms in configs:
    _, bits_list = data[n_h]
    zorder = 7 if n_h == 2 else (6 if n_h == 3 else (5 if n_h == 1 else 4))
    ax2.semilogx(sigmas, bits_list, f'{marker}-', color=color, markersize=ms,
                 linewidth=lw, markeredgecolor='white', markeredgewidth=0.8,
                 label=label, zorder=zorder)

for b, lab in [(8, '8-bit'), (12, '12-bit')]:
    ax2.axhline(y=b, color='#aaa', linewidth=0.8, linestyle=':', alpha=0.5)
    ax2.text(3e-5, b + 0.3, lab, fontsize=9, color='#888', style='italic')

ax2.fill_between([1e-4, 1e-3], [0, 0], [20, 20], alpha=0.06, color='#22aa55')

# Annotation
ax2.annotate('Newton: $\\sim$1 extra bit\n(fewer ops, but 7% bias)',
             xy=(3e-4, 12.5), fontsize=8, color='#2255aa',
             style='italic', ha='center')

ax2.set_xlabel('Per-component noise $\\sigma$', fontsize=13)
ax2.set_ylabel('Output precision (bits)', fontsize=13)
ax2.set_title('Precision: graceful degradation\n'
              'all methods scale as $-\\log_2(\\sqrt{N}\\,\\sigma)$',
              fontsize=12, fontweight='bold')
ax2.legend(fontsize=9.5, loc='upper right',
           framealpha=0.9, edgecolor='#ddd')
ax2.set_ylim(0, 16)
ax2.set_xlim(2e-5, 3e-1)
ax2.grid(True, alpha=0.15)

plt.tight_layout(pad=2.0)
plt.savefig('pandrosion_analog.png', dpi=200, bbox_inches='tight', facecolor='#fafaf8')
plt.savefig('pandrosion_analog.pdf', dpi=200, bbox_inches='tight', facecolor='#fafaf8')
print("Saved pandrosion_analog.{png,pdf}")
plt.close()
