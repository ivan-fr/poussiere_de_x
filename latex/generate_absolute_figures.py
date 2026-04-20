import numpy as np
import matplotlib.pyplot as plt

# Global aesthetic settings
plt.style.use('dark_background')
GOLD = '#FFD700'
CYAN = '#00FFFF'
MAGENTA = '#FF00FF'
NEON_GREEN = '#39FF14'
RED_DUST = '#FF0055'

def create_sml_finiteness():
    print("Generating Skolem-Mahler-Lech Finiteness...")
    fig, ax = plt.subplots(figsize=(10, 8), facecolor='#0a0a0a')
    ax.set_facecolor('#050505')
    
    n_iters = np.arange(0, 50, 0.1)
    # The linear recurrence trace oscillating around 0
    # It has a periodic mode (sin(n)) and a decaying divergence mode
    trace = np.sin(2 * n_iters) * np.exp(-0.15 * n_iters) * 10
    
    # Faltings Extinction Boundary
    faltings_env = 10 * np.exp(-0.15 * n_iters)
    
    # Granular resolution threshold (beneath this, integers cannot cross zero)
    horizon_threshold = 0.5
    N_horizon = 20 # where envelope crosses the threshold

    # Plot trace and envelope
    ax.plot(n_iters, trace, color=CYAN, lw=2, alpha=0.8, label="Recurrence Trace $a_n$")
    ax.plot(n_iters, faltings_env, color=RED_DUST, linestyle='--', lw=2, label="Faltings Phase 2 Extinction")
    ax.plot(n_iters, -faltings_env, color=RED_DUST, linestyle='--', lw=2)
    
    # The discrete exact roots
    discrete_n = np.arange(0, 50, 1)
    discrete_trace = np.sin(2 * discrete_n) * np.exp(-0.15 * discrete_n) * 10
    
    # Identify zeros
    zero_nodes = np.arange(0, 50, np.pi/2)
    valid_zeros = zero_nodes[zero_nodes < N_horizon]
    periodic_zeros = np.arange(0, 50, np.pi/2) * 0 # Flatline after N_horizon
    
    ax.scatter(valid_zeros, [0]*len(valid_zeros), color=CYAN, s=80, marker='x', lw=2, label="Finite Divergent Zeroes")
    
    # Periodic Szemeredi Nodes (standing waves)
    szemeredi = np.arange(22, 50, 3)
    ax.scatter(szemeredi, [0]*len(szemeredi), color=NEON_GREEN, s=100, zorder=5, label="Periodic Roots (Szemerédi Prog.)")
    
    # Faltings Horizon line
    ax.axvline(N_horizon, color=GOLD, lw=2, linestyle=':', label="$N_{horizon}$ (Turing Depth Limit)")
    ax.axhspan(-horizon_threshold, horizon_threshold, color='#ffffff', alpha=0.05, label="Granular Integer Threshold")
    
    ax.set_xlim(0, 40)
    ax.set_ylim(-12, 12)
    for spine in ax.spines.values():
        spine.set_color('#444')
    ax.tick_params(colors='#888')
    ax.set_title("114. Skolem-Mahler-Lech: Phase 2 Finiteness & Szemerédi Periodicity", color='white', fontsize=15)
    ax.legend(loc='upper right', facecolor='#0a0a0a', edgecolor='none', labelcolor='white')
    ax.axhline(0, color='white', alpha=0.3, lw=1)
    
    plt.tight_layout()
    plt.savefig('absolute_sml_extinction.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.savefig('absolute_sml_extinction.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

def create_roth_irrationality():
    print("Generating Liouville-Roth Bound...")
    fig, ax = plt.subplots(figsize=(10, 8), facecolor='#0a0a0a')
    ax.set_facecolor('#050505')
    
    q_vals = np.linspace(1, 100, 400)
    # The absolute Roth limit: error = 1/q^2
    roth_limit = 100 / (q_vals ** 2)
    
    # Shade the forbidden topological zone below Roth
    ax.fill_between(q_vals, 0, roth_limit, color=RED_DUST, alpha=0.2, hatch='///', label="Forbidden Jacobian Compression")
    ax.plot(q_vals, roth_limit, color=MAGENTA, lw=3, label="Absolute Limit: $1/q^{2+\epsilon}$")
    
    # Plot generated approximations sitting safely above the bound
    np.random.seed(42)
    for _ in range(50):
        q = np.random.uniform(2, 90)
        # ensure error is strictly greater than 100/q^2
        min_err = 100 / (q ** 2)
        err = min_err + np.random.uniform(0.1, 5) * (1/np.sqrt(q))
        ax.scatter(q, err, color=CYAN, s=40, alpha=0.8)
        
    ax.scatter([], [], color=CYAN, s=40, label="Diophantine Orbital Nodes $p/q$")
    
    ax.set_xlim(1, 90)
    ax.set_ylim(0, 10)
    ax.set_yscale('symlog', linthresh=0.1)
    ax.set_xlabel("Denominator Scale ($q$)", color='white', fontsize=12)
    ax.set_ylabel("Compression Error $|x - p/q|$", color='white', fontsize=12)
    
    ax.set_title("115. Liouville-Roth: Jacobian Impossibility of Infinite Compression", color='white', fontsize=15)
    ax.legend(loc='upper right', facecolor='#0a0a0a', edgecolor='none', labelcolor='white')
    
    for spine in ax.spines.values():
        spine.set_color('#444')
    ax.tick_params(colors='#888')
    
    plt.tight_layout()
    plt.savefig('absolute_roth_bound.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.savefig('absolute_roth_bound.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

def create_lehmer_void():
    print("Generating Lehmer's Spectral Void...")
    fig, ax = plt.subplots(figsize=(10, 10), facecolor='#0a0a0a')
    ax.set_facecolor('#050505')
    
    # Unit Circle (Basis measure 1)
    theta = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(theta), np.sin(theta), color=GOLD, lw=3, label="Cyclotomic Base ($\lambda=1$)")
    
    # Lehmer Limit (c = 1.17...)
    c_lehmer = 1.176
    ax.plot(c_lehmer * np.cos(theta), c_lehmer * np.sin(theta), color=CYAN, lw=2, linestyle='--', label=r"Lehmer Limit $c > 1$")
    
    # Shade the Trace Gap (The Void)
    ax.fill_between(theta, 1, c_lehmer, color='#000000', alpha=0.9)
    # create patches to fake fill_between for polar in cartesian using 2 circles
    circle_in = plt.Circle((0,0), 1, color='black', fill=False)
    circle_out = plt.Circle((0,0), c_lehmer, color='black', fill=False)
    from matplotlib.patches import Wedge
    void_wedge = Wedge((0,0), c_lehmer, 0, 360, width=(c_lehmer-1), color=RED_DUST, alpha=0.15, label="Absolute Spectral Void (Trace Gap)")
    ax.add_patch(void_wedge)
    
    # Plot Roots (Dust) strictly outside the Lehmer boundary
    np.random.seed(116)
    for _ in range(80):
        t = np.random.uniform(0, 2*np.pi)
        r = np.random.uniform(c_lehmer + 0.02, 1.8)
        ax.scatter(r*np.cos(t), r*np.sin(t), color=CYAN, s=30, alpha=0.8)
        
    ax.scatter([], [], color=CYAN, s=30, label="Pisot-Salem Asymmetric Dust")
    
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.axis('off')
    
    ax.set_title("116. Lehmer's Limit: Non-Cyclotomic Trace Gap Conservation", color='white', fontsize=15)
    ax.legend(loc='lower left', facecolor='#0a0a0a', edgecolor='none', labelcolor='white')
    
    plt.tight_layout()
    plt.savefig('absolute_lehmer_void.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.savefig('absolute_lehmer_void.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    create_sml_finiteness()
    create_roth_irrationality()
    create_lehmer_void()
    print("All Part III Absolute Vision Theorems instantiated successfully!")
