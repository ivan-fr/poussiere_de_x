"""
Generate figures for the new theorems: contraction, orbit, oscillation.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

DPI = 300
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 12

# ============================================================
# FIGURE: Contraction Factor λ(s) = (x-1)/[x(1+s)(1+s*)]
# ============================================================
def contraction_factor():
    print("Generating contraction factor plot...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    
    # LEFT: λ as a function of s for various x
    s = np.linspace(0, 3, 300)
    for x_val, color, label in [(1.5, '#2060cc', 'x=1.5'), 
                                  (2.0, '#cc3030', 'x=2'), 
                                  (5.0, '#30a030', 'x=5'),
                                  (10.0, '#9030cc', 'x=10')]:
        sstar = (1/x_val)**0.5  # p=2
        lam = (x_val - 1) / (x_val * (1 + s) * (1 + sstar))
        ax1.plot(s, lam, color=color, linewidth=2, label=label)
        ax1.axvline(sstar, color=color, linestyle=':', alpha=0.4, linewidth=1)
    
    ax1.axhline(1, color='black', linestyle='--', linewidth=1, alpha=0.5, label=r'$\lambda = 1$')
    ax1.set_xlabel(r'$s$', fontsize=14)
    ax1.set_ylabel(r'$\lambda(s)$', fontsize=14)
    ax1.set_title(r'Contraction factor $\lambda = \frac{x-1}{x(1+s)(1+s^*)}$', fontsize=13)
    ax1.set_ylim(-0.05, 1.1)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.2)
    ax1.fill_between(s, 0, 1, alpha=0.05, color='green')
    ax1.text(1.5, 0.5, r'$\lambda < 1$: contraction', fontsize=11, 
             ha='center', color='green', alpha=0.7)
    
    # RIGHT: Distance |h(s) - s*| vs |s - s*| for p=2, x=2
    x_val = 2.0
    sstar = (1/x_val)**0.5
    s_range = np.linspace(0.01, 2.5, 500)
    dist_before = np.abs(s_range - sstar)
    h_s = 1 - (x_val - 1) / (x_val * (1 + s_range))
    dist_after = np.abs(h_s - sstar)
    
    ax2.plot(dist_before, dist_after, 'b-', linewidth=2.5, label=r'$|h(s)-s^*|$ vs $|s-s^*|$')
    max_d = 2.0
    ax2.plot([0, max_d], [0, max_d], 'k--', linewidth=1, alpha=0.5, label=r'$y = x$ (no contraction)')
    ax2.set_xlabel(r'$|s - s^*|$', fontsize=14)
    ax2.set_ylabel(r'$|h(s) - s^*|$', fontsize=14)
    ax2.set_title(r'Distance decrease (Theorem 2.5), $x=2$', fontsize=13)
    ax2.set_xlim(0, max_d)
    ax2.set_ylim(0, max_d)
    ax2.set_aspect('equal')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.2)
    ax2.fill_between([0, max_d], [0, max_d], [0, 0], alpha=0.05, color='green')
    
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_contraction.pdf', 
                dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_contraction.pdf")


# ============================================================
# FIGURE: Orbit in (0,1) + oscillation signature
# ============================================================
def orbit_and_oscillation():
    print("Generating orbit + oscillation plot...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    
    # LEFT: Orbit staying in (0,1) for p=2, various x
    for x_val, color, label in [(2, '#2060cc', 'x=2'), 
                                  (5, '#cc3030', 'x=5'), 
                                  (10, '#30a030', 'x=10')]:
        s = 0.9  # start near 1
        orbit = [s]
        for _ in range(15):
            s = 1 - (x_val - 1) / (x_val * (1 + s))
            orbit.append(s)
        ax1.plot(range(len(orbit)), orbit, 'o-', color=color, linewidth=1.5, 
                markersize=4, label=label)
        sstar = (1/x_val)**0.5
        ax1.axhline(sstar, color=color, linestyle=':', alpha=0.4)
    
    ax1.axhline(0, color='gray', linewidth=1)
    ax1.axhline(1, color='gray', linewidth=1)
    ax1.fill_between(range(16), 0, 1, alpha=0.05, color='blue')
    ax1.set_xlabel('Iteration $n$', fontsize=13)
    ax1.set_ylabel(r'$h^n(s_0)$', fontsize=13)
    ax1.set_title(r'Orbit invariance (Theorem 2.6): $h^n(s_0) \in (0,1)$', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.2)
    ax1.set_ylim(-0.05, 1.05)
    
    # RIGHT: Oscillation signature for p=3
    X_val = 2.0
    root = X_val**(1/3)
    s = 2.5  # start above
    errors = []
    iterates = []
    for _ in range(12):
        iterates.append(s)
        errors.append(s - root)
        denom = 3*s**3 + 2*X_val
        s = s * (s**3 + 4*X_val) / denom
    iterates.append(s)
    errors.append(s - root)
    
    n_range = range(len(errors))
    colors_bar = ['#cc3030' if e > 0 else '#2060cc' for e in errors]
    
    ax2.bar(n_range, errors, color=colors_bar, alpha=0.7, width=0.6)
    ax2.axhline(0, color='black', linewidth=0.8)
    ax2.set_xlabel('Iteration $n$', fontsize=13)
    ax2.set_ylabel(r'$s_n - \sqrt[3]{2}$', fontsize=13)
    ax2.set_title(r'Oscillation signature (Theorem 2.7): alternating error', fontsize=13)
    ax2.grid(True, alpha=0.2)
    
    # Annotate ratio
    if len(errors) >= 3 and abs(errors[0]) > 0:
        ratio = errors[1] / errors[0]
        ax2.annotate(f'$\\lambda_3 \\approx {ratio:.3f}$', 
                    xy=(1, errors[1]), xytext=(3, errors[0]*0.5),
                    fontsize=12, color='navy',
                    arrowprops=dict(arrowstyle='->', color='navy', lw=1.5))
    
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_orbit_oscillation.pdf', 
                dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_orbit_oscillation.pdf")


# ============================================================
# FIGURE: Pell-Pandrosion multiplicative error growth
# ============================================================
def pell_multiplicative():
    print("Generating Pell-Pandrosion multiplicative error...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    
    X = 2
    A, B = 1, 1
    
    M_vals = []
    Phi_vals = []
    A_vals = [A]
    B_vals = [B]
    
    for i in range(5):
        M = A**3 - X * B**3
        M_vals.append(M)
        Phi = A**9 - 14*X*A**6*B**3 - 20*X**2*A**3*B**6 + 8*X**3*B**9
        Phi_vals.append(Phi)
        A_new = A**4 + 4*X*A*B**3
        B_new = 3*A**3*B + 2*X*B**4
        A, B = A_new, B_new
        A_vals.append(A)
        B_vals.append(B)
    M_vals.append(A**3 - X * B**3)
    
    # LEFT: |M_n| on log scale
    abs_M = [abs(m) if m != 0 else 1 for m in M_vals]
    import math
    log_M = [math.log10(m) if m > 0 else 0 for m in abs_M]
    
    ax1.bar(range(len(log_M)), log_M, color='#2060cc', alpha=0.8, width=0.6)
    ax1.set_xlabel('Iteration $n$', fontsize=13)
    ax1.set_ylabel(r'$\log_{10}|M_n|$', fontsize=13)
    ax1.set_title(r'$|M_n| = |A_n^3 - 2B_n^3|$: multiplicative growth', fontsize=13)
    ax1.grid(True, alpha=0.2)
    
    # Annotate ratios
    for i in range(min(3, len(M_vals)-1)):
        if M_vals[i] != 0:
            ratio_digits = len(str(abs(M_vals[i+1]))) / max(len(str(abs(M_vals[i]))), 1)
            ax1.annotate(f'×{ratio_digits:.0f} digits', 
                        xy=(i+0.5, (log_M[i]+log_M[i+1])/2),
                        fontsize=9, color='navy', ha='center')
    
    # RIGHT: |Φ_n| showing the polynomial factor
    log_Phi = [math.log10(abs(p)) if p != 0 else 0 for p in Phi_vals]
    ax2.bar(range(len(log_Phi)), log_Phi, color='#cc5030', alpha=0.8, width=0.6)
    ax2.set_xlabel('Iteration $n$', fontsize=13)
    ax2.set_ylabel(r'$\log_{10}|\Phi_n|$', fontsize=13)
    ax2.set_title(r'Polynomial factor $|\Phi_n|$ (degree 9)', fontsize=13)
    ax2.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_pell_identity.pdf', 
                dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_pell_identity.pdf")


if __name__ == '__main__':
    print("=" * 60)
    print("New Theorem Figures Generator")
    print("=" * 60)
    contraction_factor()
    orbit_and_oscillation()
    pell_multiplicative()
    print("=" * 60)
    print("All 3 new figures generated!")
    print("=" * 60)
