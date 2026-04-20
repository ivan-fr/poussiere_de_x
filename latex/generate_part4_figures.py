import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Wedge
import matplotlib.gridspec as gridspec

plt.style.use('dark_background')
GOLD = '#FFD700'
CYAN = '#00FFFF'
MAGENTA = '#FF00FF'
NEON_GREEN = '#39FF14'
RED_DUST = '#FF0055'
VIOLET = '#8B5CF6'
ORANGE = '#FF6B35'

def create_kronecker():
    print("Generating Kronecker Roots of Unity...")
    fig, ax = plt.subplots(figsize=(10, 10), facecolor='#0a0a0a')
    ax.set_facecolor('#050505')
    
    theta = np.linspace(0, 2*np.pi, 300)
    # Unit circle
    ax.plot(np.cos(theta), np.sin(theta), color=GOLD, lw=3, label="Unit Circle $|z|=1$")
    
    # Lehmer void (from Module 116)
    c_lehmer = 1.176
    ax.plot(c_lehmer*np.cos(theta), c_lehmer*np.sin(theta), color=RED_DUST, lw=1.5, linestyle='--', alpha=0.5)
    void_wedge = Wedge((0,0), c_lehmer, 0, 360, width=(c_lehmer-1), color=RED_DUST, alpha=0.08)
    ax.add_patch(void_wedge)
    
    # Roots of unity (the ONLY algebraic integers on the circle)
    for n in [3, 4, 5, 6, 7, 8, 12]:
        angles = np.linspace(0, 2*np.pi, n, endpoint=False)
        xs = np.cos(angles)
        ys = np.sin(angles)
        ax.scatter(xs, ys, s=120, color=NEON_GREEN, zorder=5, edgecolors='white', linewidths=0.5)
        # Connect them
        for i in range(n):
            ax.plot([xs[i], xs[(i+1)%n]], [ys[i], ys[(i+1)%n]], color=NEON_GREEN, alpha=0.15, lw=1)
    
    ax.scatter([], [], color=NEON_GREEN, s=120, label="Roots of Unity (Cyclotomic)")
    
    # Arrows showing collapse: any algebraic integer ON circle → must be root of unity
    for angle in [0.7, 2.1, 3.8, 5.2]:
        ax.annotate('', xy=(np.cos(angle), np.sin(angle)),
                    xytext=(1.35*np.cos(angle), 1.35*np.sin(angle)),
                    arrowprops=dict(arrowstyle='->', color=CYAN, lw=1.5))
    
    ax.text(0, 0, "Mahler\n$M=1$", ha='center', va='center', color=GOLD, fontsize=11, fontweight='bold')
    
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.8, 1.8)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title("117. Kronecker: All Conjugates on $|z|=1$ $\\Rightarrow$ Root of Unity", color='white', fontsize=14)
    ax.legend(loc='lower left', facecolor='#0a0a0a', edgecolor='none', labelcolor='white')
    
    plt.tight_layout()
    plt.savefig('part4_kronecker.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.savefig('part4_kronecker.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

def create_weyl():
    print("Generating Weyl Equidistribution...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor='#0a0a0a')
    
    alpha = (1 + np.sqrt(5)) / 2  # golden ratio
    N = 500
    orbit = np.array([(n * alpha) % 1 for n in range(1, N+1)])
    
    # Left: orbit scatter on [0,1)
    ax = axes[0]
    ax.set_facecolor('#050505')
    ax.scatter(range(1, N+1), orbit, s=1.5, color=CYAN, alpha=0.7)
    ax.set_xlabel("Iteration $n$", color='white', fontsize=11)
    ax.set_ylabel("$n\\alpha$ mod 1", color='white', fontsize=11)
    ax.set_title("Orbital Trace $(n\\varphi$ mod $1)$", color=GOLD, fontsize=13)
    for spine in ax.spines.values():
        spine.set_color('#333')
    ax.tick_params(colors='#888')
    
    # Right: histogram showing uniform distribution
    ax2 = axes[1]
    ax2.set_facecolor('#050505')
    counts, bins, patches = ax2.hist(orbit, bins=30, color=CYAN, alpha=0.7, edgecolor='#050505')
    expected = N / 30
    ax2.axhline(expected, color=GOLD, lw=2, linestyle='--', label=f"Expected Uniform: {expected:.1f}")
    ax2.set_xlabel("Position in $[0,1)$", color='white', fontsize=11)
    ax2.set_ylabel("Count", color='white', fontsize=11)
    ax2.set_title("Tensor Conservation $\\Rightarrow$ Uniform Density", color=GOLD, fontsize=13)
    ax2.legend(facecolor='#0a0a0a', edgecolor='none', labelcolor='white')
    for spine in ax2.spines.values():
        spine.set_color('#333')
    ax2.tick_params(colors='#888')
    
    fig.suptitle("118. Weyl Equidistribution: Phase Space Volume Conservation", color='white', fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig('part4_weyl.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.savefig('part4_weyl.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

def create_gelfond_schneider():
    print("Generating Gelfond-Schneider Transcendence...")
    fig, ax = plt.subplots(figsize=(10, 8), facecolor='#0a0a0a')
    ax.set_facecolor('#050505')
    
    q_vals = np.linspace(1, 200, 500)
    
    # Algebraic bound (Roth limit)
    roth = 1.0 / (q_vals ** 2)
    ax.plot(q_vals, roth, color=MAGENTA, lw=2.5, label="Algebraic Roth Limit: $1/q^2$")
    ax.fill_between(q_vals, 0, roth, color=MAGENTA, alpha=0.1)
    
    # Transcendental compression rate (exceeds Roth)
    trans_rate = 1.0 / (q_vals ** 1.3)  # slower decay = better approximation
    ax.plot(q_vals, trans_rate, color=NEON_GREEN, lw=2.5, linestyle='--',
            label="$2^{\\sqrt{2}}$ Compression Rate")
    
    # The gap zone
    ax.fill_between(q_vals, roth, trans_rate, color=NEON_GREEN, alpha=0.08,
                    label="Transcendence Gap (Violation Zone)")
    
    # Scattered approximation points for 2^sqrt(2)
    np.random.seed(119)
    for _ in range(40):
        q = np.random.uniform(5, 180)
        approx = 1.0 / (q ** (1.3 + np.random.uniform(-0.15, 0.15)))
        ax.scatter(q, approx, color=ORANGE, s=35, alpha=0.8, zorder=4)
    ax.scatter([], [], color=ORANGE, s=35, label="Approximations of $\\alpha^\\beta$")
    
    ax.set_yscale('log')
    ax.set_xlim(1, 200)
    ax.set_ylim(1e-6, 1)
    ax.set_xlabel("Denominator $q$", color='white', fontsize=12)
    ax.set_ylabel("Approximation Error", color='white', fontsize=12)
    ax.set_title("119. Gelfond-Schneider: $\\alpha^\\beta$ Breaks the Algebraic Barrier", color='white', fontsize=14)
    ax.legend(loc='upper right', facecolor='#0a0a0a', edgecolor='none', labelcolor='white')
    for spine in ax.spines.values():
        spine.set_color('#333')
    ax.tick_params(colors='#888')
    
    plt.tight_layout()
    plt.savefig('part4_gelfond.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.savefig('part4_gelfond.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

def create_baker():
    print("Generating Baker Linear Forms...")
    fig, ax = plt.subplots(figsize=(10, 8), facecolor='#0a0a0a')
    ax.set_facecolor('#050505')
    
    B_vals = np.linspace(1, 100, 300)
    
    # Baker lower bound for different dimensions n
    for n, color, alpha_v in [(2, CYAN, 0.9), (3, NEON_GREEN, 0.7), (5, MAGENTA, 0.5)]:
        C = 0.5 / (n**2)
        kappa = n + 1
        baker_bound = C / (B_vals ** kappa)
        ax.plot(B_vals, baker_bound, color=color, lw=2.5, 
                label=f"$n={n}$: $C/B^{{{kappa}}}$")
        ax.fill_between(B_vals, 0, baker_bound, color=color, alpha=0.05)
    
    # Forbidden zone
    ax.fill_between(B_vals, 0, 1e-15, color=RED_DUST, alpha=0.3)
    ax.text(50, 3e-14, "VANISHING FORBIDDEN", ha='center', color=RED_DUST, fontsize=11, fontweight='bold')
    
    ax.set_yscale('log')
    ax.set_xlim(1, 100)
    ax.set_ylim(1e-15, 1)
    ax.set_xlabel("Coefficient Height $B = \\max|b_i|$", color='white', fontsize=12)
    ax.set_ylabel("$|b_1 \\log \\alpha_1 + \\cdots + b_n \\log \\alpha_n|$", color='white', fontsize=12)
    ax.set_title("120. Baker: Multi-Dimensional Logarithmic Tensor Bounds", color='white', fontsize=14)
    ax.legend(loc='upper right', facecolor='#0a0a0a', edgecolor='none', labelcolor='white')
    for spine in ax.spines.values():
        spine.set_color('#333')
    ax.tick_params(colors='#888')
    
    plt.tight_layout()
    plt.savefig('part4_baker.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.savefig('part4_baker.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

def create_schmidt():
    print("Generating Schmidt Subspace...")
    fig, ax = plt.subplots(figsize=(10, 10), facecolor='#0a0a0a')
    ax.set_facecolor('#050505')
    
    # Generate random integer-like points in 2D
    np.random.seed(121)
    
    # The finite union of subspaces (lines through the origin)
    angles = [np.pi/6, np.pi/3, 2*np.pi/3, 5*np.pi/6]
    t_vals = np.linspace(-8, 8, 200)
    
    for i, angle in enumerate(angles):
        xs = t_vals * np.cos(angle)
        ys = t_vals * np.sin(angle)
        ax.plot(xs, ys, color=GOLD, lw=2, alpha=0.6)
        # Points clustered on subspaces
        for _ in range(15):
            t = np.random.uniform(-6, 6)
            jitter = np.random.normal(0, 0.05, 2)
            ax.scatter(t*np.cos(angle) + jitter[0], t*np.sin(angle) + jitter[1],
                      color=CYAN, s=50, zorder=5, alpha=0.9)
    
    ax.plot([], [], color=GOLD, lw=2, label="Finite Subspaces (Invariant Manifolds)")
    ax.scatter([], [], color=CYAN, s=50, label="Integer Solutions (Collapsed)")
    
    # Show the "forbidden" full-dimensional zone
    circle = plt.Circle((0, 0), 7, color=RED_DUST, fill=False, lw=1.5, linestyle=':', alpha=0.5)
    ax.add_patch(circle)
    
    # Scattered "impossible" full-dimensional points (very few, fading)
    for _ in range(8):
        x = np.random.uniform(-5, 5)
        y = np.random.uniform(-5, 5)
        ax.scatter(x, y, color=RED_DUST, s=30, alpha=0.3, marker='x', zorder=3)
    ax.scatter([], [], color=RED_DUST, s=30, marker='x', label="Full-Rank Exceptions (Finite)")
    
    ax.set_xlim(-8, 8)
    ax.set_ylim(-8, 8)
    ax.set_aspect('equal')
    ax.set_title("121. Schmidt Subspace: Solutions Collapse onto Finite Invariant Manifolds", color='white', fontsize=14)
    ax.legend(loc='lower left', facecolor='#0a0a0a', edgecolor='none', labelcolor='white')
    for spine in ax.spines.values():
        spine.set_color('#333')
    ax.tick_params(colors='#888')
    
    plt.tight_layout()
    plt.savefig('part4_schmidt.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.savefig('part4_schmidt.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    create_kronecker()
    create_weyl()
    create_gelfond_schneider()
    create_baker()
    create_schmidt()
    print("All Part IV figures generated!")
