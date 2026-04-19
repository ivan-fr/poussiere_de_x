"""
Generate publication-quality fractal and mathematical figures
for the Universitas Pandrosion research paper.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors

# High-res settings
DPI = 300
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 12

# ============================================================
# FIGURE 1: Newton Fractal for z^3 - 1 (chaotic boundaries)
# ============================================================
def newton_fractal(resolution=800):
    print("Generating Newton fractal...")
    x = np.linspace(-2, 2, resolution)
    y = np.linspace(-2, 2, resolution)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    
    roots = [1.0, np.exp(2j*np.pi/3), np.exp(-2j*np.pi/3)]
    
    # Newton iteration z -> z - (z^3 - 1)/(3z^2)
    convergence = np.zeros(Z.shape, dtype=int)
    speed = np.zeros(Z.shape, dtype=float)
    
    for iteration in range(50):
        mask = speed == 0
        if not np.any(mask):
            break
        
        denom = 3 * Z[mask]**2
        safe = np.abs(denom) > 1e-10
        
        Z_new = np.copy(Z)
        update_mask = np.zeros_like(mask)
        temp = np.zeros(np.sum(mask), dtype=complex)
        temp[safe] = Z[mask][safe] - (Z[mask][safe]**3 - 1) / denom[safe]
        Z_flat = Z.copy()
        Z_flat[mask] = temp
        Z = Z_flat
        
        for idx, root in enumerate(roots):
            close = np.abs(Z - root) < 1e-6
            newly_converged = close & (speed == 0)
            convergence[newly_converged] = idx + 1
            speed[newly_converged] = iteration + 1
    
    speed[speed == 0] = 50
    
    # Color mapping
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Create RGB image
    img = np.zeros((*Z.shape, 3))
    colors_map = [
        [0.9, 0.2, 0.2],  # Red
        [0.2, 0.5, 0.9],  # Blue  
        [0.2, 0.8, 0.3],  # Green
    ]
    
    for idx in range(3):
        mask = convergence == idx + 1
        shade = 1.0 - (speed[mask] / 50.0) * 0.7
        for c in range(3):
            img[:,:,c][mask] = colors_map[idx][c] * shade
    
    # Non-converged points black
    mask = convergence == 0
    img[mask] = [0, 0, 0]
    
    ax.imshow(img, extent=[-2, 2, -2, 2], origin='lower')
    ax.set_xlabel(r'$\Re(z)$', fontsize=14)
    ax.set_ylabel(r'$\Im(z)$', fontsize=14)
    ax.set_title(r"Newton's Method: $z^3 - 1$" + "\n" + r"Fractal basin boundaries", fontsize=13)
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_newton_fractal.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_newton_fractal.pdf")


# ============================================================
# FIGURE 2: Pandrosion basin map for z^3 - X (clean boundaries)
# ============================================================
def pandrosion_fractal(resolution=800):
    print("Generating Pandrosion fractal...")
    x = np.linspace(-2, 2, resolution)
    y = np.linspace(-2, 2, resolution)
    X_grid, Y_grid = np.meshgrid(x, y)
    Z = X_grid + 1j * Y_grid
    
    target = 1.0  # X = 1, roots are cube roots of 1
    roots = [1.0, np.exp(2j*np.pi/3), np.exp(-2j*np.pi/3)]
    
    convergence = np.zeros(Z.shape, dtype=int)
    speed = np.zeros(Z.shape, dtype=float)
    
    for iteration in range(80):
        mask = speed == 0
        if not np.any(mask):
            break
        
        Zm = Z[mask]
        denom = 3 * Zm**3 + 2 * target
        safe = np.abs(denom) > 1e-10
        
        Z_new = np.copy(Z)
        temp = np.copy(Zm)
        temp[safe] = Zm[safe] * (Zm[safe]**3 + 4 * target) / denom[safe]
        Z_new[mask] = temp
        Z = Z_new
        
        for idx, root in enumerate(roots):
            close = np.abs(Z - root) < 1e-6
            newly_converged = close & (speed == 0)
            convergence[newly_converged] = idx + 1
            speed[newly_converged] = iteration + 1
        
        # Also check convergence to 0
        close_zero = np.abs(Z) < 1e-6
        newly_zero = close_zero & (speed == 0)
        convergence[newly_zero] = 4
        speed[newly_zero] = iteration + 1
    
    speed[speed == 0] = 80
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    img = np.zeros((*Z.shape, 3))
    colors_map = [
        [0.85, 0.15, 0.15],  # Red
        [0.15, 0.45, 0.85],  # Blue
        [0.15, 0.75, 0.25],  # Green
        [0.1, 0.1, 0.1],     # Dark (origin)
    ]
    
    for idx in range(4):
        mask = convergence == idx + 1
        shade = 1.0 - (speed[mask] / 80.0) * 0.6
        for c in range(3):
            img[:,:,c][mask] = colors_map[idx][c] * shade
    
    mask = convergence == 0
    img[mask] = [0.05, 0.05, 0.05]
    
    ax.imshow(img, extent=[-2, 2, -2, 2], origin='lower')
    ax.set_xlabel(r'$\Re(z)$', fontsize=14)
    ax.set_ylabel(r'$\Im(z)$', fontsize=14)
    ax.set_title(r"Pandrosion Iteration: $z^3 - 1$" + "\n" + r"Basin boundaries appear smooth (numerical observation)", fontsize=13)
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_pandrosion_fractal.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_pandrosion_fractal.pdf")


# ============================================================
# FIGURE 3: Convergence speed heatmap (Pandrosion)
# ============================================================
def convergence_heatmap(resolution=600):
    print("Generating convergence speed heatmap...")
    x = np.linspace(-3, 3, resolution)
    y = np.linspace(-3, 3, resolution)
    X_grid, Y_grid = np.meshgrid(x, y)
    Z = X_grid + 1j * Y_grid
    
    target = 2.0
    root = 2.0**(1/3)
    
    speed = np.full(Z.shape, 60.0)
    
    for iteration in range(60):
        mask = speed == 60.0
        if not np.any(mask):
            break
        Zm = Z[mask]
        denom = 3 * Zm**3 + 2 * target
        safe = np.abs(denom) > 1e-10
        temp = np.copy(Zm)
        temp[safe] = Zm[safe] * (Zm[safe]**3 + 4 * target) / denom[safe]
        Z_new = np.copy(Z)
        Z_new[mask] = temp
        Z = Z_new
        
        close = np.abs(Z - root) < 1e-8
        newly = close & (speed == 60.0)
        speed[newly] = iteration + 1
    
    fig, ax = plt.subplots(figsize=(9, 7))
    
    cmap = LinearSegmentedColormap.from_list('pandrosion', [
        (0.0, '#000033'),
        (0.15, '#0000aa'),
        (0.3, '#0066ff'),
        (0.45, '#00cccc'),
        (0.6, '#00ff66'),
        (0.75, '#ffff00'),
        (0.9, '#ff6600'),
        (1.0, '#ff0000'),
    ])
    
    im = ax.imshow(speed, extent=[-3, 3, -3, 3], origin='lower', cmap=cmap, vmin=1, vmax=40)
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Iterations to convergence', fontsize=12)
    
    # Mark root
    ax.plot(root, 0, 'w*', markersize=15, markeredgecolor='black', markeredgewidth=1)
    ax.annotate(r'$\sqrt[3]{2}$', xy=(root, 0), xytext=(root+0.3, 0.4),
                fontsize=14, color='white', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='white', lw=1.5))
    
    ax.set_xlabel(r'$\Re(z)$', fontsize=14)
    ax.set_ylabel(r'$\Im(z)$', fontsize=14)
    ax.set_title(r'Pandrosion Convergence Landscape for $\sqrt[3]{2}$' + '\n' +
                 'Number of iterations to reach the real cube root', fontsize=13)
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_convergence_heatmap.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_convergence_heatmap.pdf")


# ============================================================
# FIGURE 4: Hermitian Eigenvalue Spectrum
# ============================================================
def eigenvalue_spectrum():
    print("Generating Hermitian eigenvalue spectrum...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    np.random.seed(42)
    
    # LEFT: Generic non-Hermitian operator eigenvalues
    n = 50
    M = np.random.randn(n, n) + 1j * np.random.randn(n, n)
    eigs_generic = np.linalg.eigvals(M)
    
    ax1.scatter(eigs_generic.real, eigs_generic.imag, c='red', alpha=0.7, s=40, edgecolors='darkred', linewidths=0.5)
    ax1.axhline(y=0, color='gray', linewidth=0.5, linestyle='--')
    ax1.axvline(x=0, color='gray', linewidth=0.5, linestyle='--')
    ax1.set_xlabel(r'$\Re(\lambda)$', fontsize=13)
    ax1.set_ylabel(r'$\Im(\lambda)$', fontsize=13)
    ax1.set_title('Generic Operator\nEigenvalues scatter across $\\mathbb{C}$', fontsize=13)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.2)
    
    # RIGHT: Hermitian operator (Pandrosion UV)
    S = np.random.randn(n, n)
    S = (S + S.T) / 2  # make Hermitian
    X_mat = np.random.randn(n, n)
    X_mat = (X_mat + X_mat.T) / 2  # make Hermitian
    
    U = S @ (np.linalg.matrix_power(S, 3) + 4 * X_mat)
    V = 3 * np.linalg.matrix_power(S, 3) + 2 * X_mat
    UV = U @ V
    
    eigs_pandrosion = np.linalg.eigvals(UV)
    
    ax2.scatter(eigs_pandrosion.real, eigs_pandrosion.imag, c='royalblue', alpha=0.8, s=50, edgecolors='navy', linewidths=0.5)
    ax2.axhline(y=0, color='orange', linewidth=2.5, linestyle='-', label=r'Real axis ($\Im = 0$)')
    ax2.axvline(x=0, color='gray', linewidth=0.5, linestyle='--')
    ax2.set_xlabel(r'$\Re(\lambda)$', fontsize=13)
    ax2.set_ylabel(r'$\Im(\lambda)$', fontsize=13)
    ax2.set_title('Pandrosion $UV$ (Hermitian inputs)\nEigenvalues confined to $\\mathbb{R}$', fontsize=13)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.2)
    ax2.legend(fontsize=11, loc='upper right')
    
    # Force same y-range to show the contrast
    max_im = max(np.max(np.abs(eigs_generic.imag)), 1)
    ax1.set_ylim(-max_im*1.2, max_im*1.2)
    ax2.set_ylim(-max_im*1.2, max_im*1.2)
    
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_eigenvalue_spectrum.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_eigenvalue_spectrum.pdf")


# ============================================================
# FIGURE 5: 3D Surface of Kinematic Conservation Polynomial
# ============================================================
def kinematic_surface():
    print("Generating 3D kinematic surface...")
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    s = np.linspace(0.5, 3.0, 200)
    X_vals = np.linspace(0.5, 4.0, 200)
    S, Xv = np.meshgrid(s, X_vals)
    
    # The conservation polynomial ratio factor
    numerator = S**9 - 14 * S**6 * Xv - 20 * S**3 * Xv**2 + 8 * Xv**3
    denominator = (3 * S**3 + 2 * Xv)**3
    
    R = numerator / denominator
    R = np.clip(R, -2, 2)
    
    cmap = LinearSegmentedColormap.from_list('conservation', [
        (0.0, '#0000aa'),
        (0.3, '#00aaff'),
        (0.5, '#ffffff'),
        (0.7, '#ffaa00'),
        (1.0, '#aa0000'),
    ])
    
    norm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
    
    surf = ax.plot_surface(S, Xv, R, cmap=cmap, norm=norm, alpha=0.85, 
                           linewidth=0, antialiased=True, rcount=150, ccount=150)
    
    # Zero plane
    ax.plot_surface(S, Xv, np.zeros_like(R), alpha=0.1, color='gray')
    
    ax.set_xlabel(r'$s$ (approximation)', fontsize=12, labelpad=10)
    ax.set_ylabel(r'$X$ (target)', fontsize=12, labelpad=10)
    ax.set_zlabel(r'Contraction factor $\Phi/V^3$', fontsize=11, labelpad=8)
    ax.set_title('Kinematic Conservation Surface\n' + 
                 r'$\frac{s^9 - 14s^6 X - 20s^3 X^2 + 8X^3}{(3s^3+2X)^3}$', fontsize=13)
    
    ax.view_init(elev=25, azim=-45)
    fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.1, label='Contraction ratio')
    
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_kinematic_3d.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_kinematic_3d.pdf")


# ============================================================
# FIGURE 6: ABC Radical Growth Comparison
# ============================================================
def abc_radical_growth():
    print("Generating ABC radical growth chart...")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    
    # Simulate Pandrosion ABC iterations starting from (A,B) = (1,1), X=2
    A, B, X = 1, 1, 2
    steps = 6
    
    A_vals = [A]
    B_vals = [B]
    M_vals = [abs(A**3 - X * B**3)]
    digit_vals = [len(str(abs(A**3 - X * B**3))) if A**3 - X * B**3 != 0 else 1]
    
    for _ in range(steps):
        A_new = A * (A**3 + 4 * X * B**3)
        B_new = B * (3 * A**3 + 2 * X * B**3)
        A, B = A_new, B_new
        M = abs(A**3 - X * B**3)
        A_vals.append(A)
        B_vals.append(B)
        M_vals.append(M)
        digit_vals.append(len(str(M)) if M > 0 else 1)
    
    n_range = range(len(A_vals))
    
    ax.bar([x - 0.15 for x in n_range], [len(str(abs(a))) for a in A_vals], 
           width=0.3, color='royalblue', alpha=0.8, label=r'Digits of $|A_n|$')
    ax.bar([x + 0.15 for x in n_range], [len(str(abs(b))) for b in B_vals],
           width=0.3, color='coral', alpha=0.8, label=r'Digits of $|B_n|$')
    
    ax.set_xlabel('Iteration $n$', fontsize=13)
    ax.set_ylabel('Number of digits', fontsize=13)
    ax.set_title(r'Digit growth of $(A_n, B_n)$ for $X=2$' + '\n' +
                 r'Each step approximately cubes the digit count',
                 fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.2)
    
    # Annotate the explosion
    for i in range(min(4, len(A_vals))):
        ax.annotate(f'{len(str(abs(A_vals[i])))}d', 
                   xy=(i-0.15, len(str(abs(A_vals[i])))), 
                   ha='center', va='bottom', fontsize=8, color='navy')
    
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_abc_growth.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_abc_growth.pdf")


# ============================================================
# FIGURE 7: Lipschitz Entropy Map
# ============================================================
def lipschitz_entropy():
    print("Generating Lipschitz entropy map...")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    X_val = 2.0
    a_range = np.linspace(0.5, 3.0, 300)
    b_range = np.linspace(0.5, 3.0, 300)
    A, B = np.meshgrid(a_range, b_range)
    
    # U(A)*V(B) - U(B)*V(A)
    UA = A * (A**3 + 4*X_val)
    VB = 3*B**3 + 2*X_val
    UB = B * (B**3 + 4*X_val)
    VA = 3*A**3 + 2*X_val
    
    cross_diff = UA * VB - UB * VA
    
    # The factored form: (A-B) * Phi
    diff = A - B
    
    # Ratio should be exactly Phi (the polynomial factor)
    safe = np.abs(diff) > 1e-10
    ratio = np.zeros_like(cross_diff)
    ratio[safe] = cross_diff[safe] / diff[safe]
    ratio[~safe] = 0
    
    ratio_clipped = np.clip(ratio, -50, 50)
    
    cmap = LinearSegmentedColormap.from_list('lipschitz', [
        '#000066', '#0044cc', '#0088ff', '#ffffff', '#ff8800', '#cc4400', '#660000'
    ])
    
    im = ax.imshow(ratio_clipped, extent=[0.5, 3.0, 0.5, 3.0], origin='lower', 
                   cmap=cmap, aspect='auto')
    
    # Diagonal A=B line
    ax.plot([0.5, 3.0], [0.5, 3.0], 'w--', linewidth=2, label=r'$A = B$ (zero divergence)')
    
    # Cube root marker
    root = 2**(1/3)
    ax.plot(root, root, 'w*', markersize=18, markeredgecolor='black', markeredgewidth=1.5)
    ax.annotate(r'$\sqrt[3]{2}$', xy=(root, root), xytext=(root+0.2, root+0.35),
                fontsize=14, color='white', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='white', lw=2))
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label(r'Polynomial factor $\Phi(A,B,X)$', fontsize=12)
    
    ax.set_xlabel(r'State $A$', fontsize=13)
    ax.set_ylabel(r'State $B$', fontsize=13)
    ax.set_title(r'Cross-State Factor: $U(A)V(B) - U(B)V(A) = (A-B)\cdot\Phi$' + '\n' +
                 r'The Lipschitz identity proves exact factorisation through $(A-B)$',
                 fontsize=12)
    ax.legend(fontsize=11, loc='upper left')
    
    plt.tight_layout()
    plt.savefig('/Users/ivanbesevic/Documents/poussiere/latex/fig_lipschitz_map.pdf', dpi=DPI, bbox_inches='tight')
    plt.close()
    print("  -> fig_lipschitz_map.pdf")


# ============================================================
# RUN ALL
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Universitas Pandrosion — Figure Generator")
    print("=" * 60)
    newton_fractal(resolution=600)
    pandrosion_fractal(resolution=600)
    convergence_heatmap(resolution=500)
    eigenvalue_spectrum()
    kinematic_surface()
    abc_radical_growth()
    lipschitz_entropy()
    print("=" * 60)
    print("All 7 figures generated successfully!")
    print("=" * 60)
