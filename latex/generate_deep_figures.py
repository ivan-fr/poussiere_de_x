import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import logging

logging.basicConfig(level=logging.INFO)

# ==========================================
# Figure 1: Deep Theorem 14 (Smale 17 Basin)
# ==========================================
def pandrosion_map(z, p=3, X=2):
    # Newton step: N(z) = z - (z^p - X)/(p z^{p-1})
    # Pandrosion corrects this. For p=3: 
    # A' = A(A^3 + 4X), B' = B(3A^3 + 2X) for A, B where z = A/B
    # In continuous complex variables:
    z3 = z**3
    numerator = z * (z3 + 4*X)
    denominator = 3*z3 + 2*X
    
    # Avoid division by zero
    mask = np.abs(denominator) > 1e-12
    out = np.copy(z)
    out[mask] = numerator[mask] / denominator[mask]
    return out

def generate_smale_basin():
    logging.info("Generating Smale-17 Basin (Deep Theorem 14)...")
    res = 800
    x = np.linspace(-2.5, 2.5, res)
    y = np.linspace(-2.5, 2.5, res)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    
    iterations = np.zeros_like(Z, dtype=int)
    max_iter = 15
    
    roots = [
        2**(1/3),
        2**(1/3) * np.exp(2j * np.pi / 3),
        2**(1/3) * np.exp(4j * np.pi / 3)
    ]
    
    for i in range(max_iter):
        Z = pandrosion_map(Z)
        for idx, r in enumerate(roots):
            converged = np.abs(Z - r) < 1e-2
            # Update iterations only if not already converged
            update_mask = converged & (iterations == 0)
            iterations[update_mask] = i + 1

    plt.figure(figsize=(10, 8), dpi=300)
    
    # Custom colormap for spectacular deep space look
    cmap = plt.cm.magma
    plt.imshow(iterations, extent=[-2.5, 2.5, -2.5, 2.5], cmap=cmap, origin='lower')
    
    # Overlay theoretical Phase 1 / Phase 2 markers
    circle1 = plt.Circle((0, 0), 2**(1/3) * 1.5, color='white', fill=False, linestyle='--', alpha=0.5, label='Cauchy Limit (Phase 1 Entry)')
    circle2 = plt.Circle((0, 0), 2**(1/3), color='cyan', fill=False, linestyle=':', alpha=0.8, label='Root Locus (Phase 2 Epochs)')
    plt.gca().add_artist(circle1)
    plt.gca().add_artist(circle2)
    
    for r in roots:
        plt.plot(r.real, r.imag, 'w*', markersize=12, markeredgecolor='black')
        
    plt.title("Theorem 14 (Smale 17 Resolved): Universal Descent Basins", fontsize=16, color='white')
    plt.xlabel("Re(z)", fontsize=12)
    plt.ylabel("Im(z)", fontsize=12)
    
    # Dark theme formatting
    plt.gca().set_facecolor('black')
    plt.gcf().patch.set_facecolor('#1c1c1c')
    plt.gca().tick_params(colors='white')
    plt.gca().xaxis.label.set_color('white')
    plt.gca().yaxis.label.set_color('white')
    legend = plt.legend(loc='upper right', facecolor='black', edgecolor='white', labelcolor='white')
    
    plt.tight_layout()
    plt.savefig('latex/deep_theorem_14_smale.png', facecolor='#1c1c1c')
    plt.close()
    logging.info("Figure 1 saved to latex/deep_theorem_14_smale.png")

# ==========================================
# Figure 2: Deep Theorems 11, 13, 15 (Heights)
# ==========================================
def generate_vojta_heights():
    logging.info("Generating Vojta Heights (Deep Theorems 11, 13, 15)...")
    
    # Compute d_n = A_n^3 - 2 B_n^3
    n_steps = 10
    
    # Start 1: (1, 1) -> target 2^(1/3)
    A1, B1 = 1, 1
    d_vals_1 = []
    
    # Start 2: (5, 4) -> closer start
    A2, B2 = 5, 4
    d_vals_2 = []
    
    for _ in range(n_steps):
        d1 = abs(A1**3 - 2*B1**3)
        d_vals_1.append(float(d1) if d1 < 1e300 else 1e300)
        A1_new = A1 * (A1**3 + 8*B1**3)
        B1_new = B1 * (3*A1**3 + 4*B1**3)
        A1, B1 = A1_new, B1_new
        
        d2 = abs(A2**3 - 2*B2**3)
        d_vals_2.append(float(d2) if d2 < 1e300 else 1e300)
        A2_new = A2 * (A2**3 + 8*B2**3)
        B2_new = B2 * (3*A2**3 + 4*B2**3)
        A2, B2 = A2_new, B2_new

    # Vojta theoretical bound: log |d_n| >= n log 2 + log |d_0|
    n_arr = np.arange(n_steps)
    vojta_bound_1 = d_vals_1[0] * (2 ** n_arr)
    vojta_bound_2 = d_vals_2[0] * (2 ** n_arr)
    
    plt.figure(figsize=(10, 6), dpi=300)
    plt.gcf().patch.set_facecolor('#f4f4f4')
    
    plt.plot(n_arr, np.log10(d_vals_1), 'o-', color='#e3e', linewidth=2.5, markersize=8, label=r'Orbital Height: Start (1,1)')
    plt.plot(n_arr, np.log10(vojta_bound_1), '--', color='#e3e', alpha=0.5, label=r'Vojta Bound (Thm 15): $n \log_{10} 2 + \log_{10} |d_0|$')
    
    plt.plot(n_arr, np.log10(d_vals_2), 's-', color='#3ae', linewidth=2.5, markersize=8, label=r'Orbital Height: Start (5,4)')
    plt.plot(n_arr, np.log10(vojta_bound_2), '--', color='#3ae', alpha=0.5, label=r'Vojta Bound (Start 5,4)')
    
    plt.title("Theorems 11, 13, 15: Thue Amplification & Vojta Height Growth", fontsize=16, fontweight='bold', color='#222')
    plt.xlabel("Iteration $n$", fontsize=14)
    plt.ylabel(r"$\log_{10} |d_n|$ (Orbital Complexity)", fontsize=14)
    plt.grid(True, linestyle=':', alpha=0.7, color='gray')
    plt.legend(loc='upper left', fontsize=11, framealpha=0.9)
    
    plt.xticks(n_arr)
    plt.tight_layout()
    plt.savefig('latex/deep_theorem_15_vojta.png')
    plt.close()
    logging.info("Figure 2 saved to latex/deep_theorem_15_vojta.png")

if __name__ == "__main__":
    generate_smale_basin()
    generate_vojta_heights()
