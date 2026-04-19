import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Use a professional style suitable for a thesis/paper
plt.style.use('seaborn-v0_8-paper')
mpl.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 300,
    'font.family': 'serif'
})

def compute_pandrosion_scan(poly, max_r, d, N_samples=200):
    """
    Simulates the continuous Pandrosion scanner on the circle R.
    For the continuous profile we use a fine grid.
    For the DFT we use exactly d samples.
    """
    # Continuous profile for visualization
    theta = np.linspace(0, 2*np.pi, N_samples)
    a_cont = max_r * np.exp(1j * theta)
    z_cont = max_r * np.exp(1j * (theta + np.pi/d))
    ratio_cont = np.abs(poly(z_cont) / poly(a_cont))
    
    # Discrete Pandrosion vector for the DFT
    s = np.arange(d)
    a_disc = max_r * np.exp(1j * 2 * np.pi * s / d)
    z_disc = max_r * np.exp(1j * (2 * np.pi * s / d + np.pi/d))
    r_vec = poly(z_disc) / poly(a_disc)
    
    # FFT
    r_hat = np.fft.fft(r_vec) / d
    return theta, ratio_cont, np.abs(r_hat)

fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# Define two Belyi polynomials (or just nice drawings)
# 1. Symmetric Star: P(z) = z^4 (passport: complete symmetry)
p_sym = lambda z: z**4
deg_sym = 4

# 2. Asymmetric Tree: P(z) = (256/27) z^3 (1-z) 
p_asym = lambda z: (256/27) * z**3 * (1 - z)
deg_asym = 4
R_cauchy = 2.0

# --- PANEL 1: Symmetric Root Geometry ---
ax1 = axs[0, 0]
# Roots are just 0 (quadruple)
ax1.plot(0, 0, 'ko', markersize=12, label='Degenerate Root (x4)')
# Draw Cauchy Circle
circle = plt.Circle((0, 0), R_cauchy, color='blue', fill=False, linestyle='--', linewidth=1.5, label=f'Cauchy Scanner R={R_cauchy}')
ax1.add_patch(circle)
ax1.set_xlim(-2.5, 2.5)
ax1.set_ylim(-2.5, 2.5)
ax1.set_aspect('equal')
ax1.set_title("Symmetric Star Dessin ($\\beta(z) = z^4$)")
ax1.legend(loc='upper right')
ax1.grid(True, linestyle=':', alpha=0.6)

# --- PANEL 2: Asymmetric Root Geometry ---
ax2 = axs[0, 1]
# Roots: 0 (triple), 1 (simple)
ax2.plot(0, 0, 'ko', markersize=10, label='Root (x3)')
ax2.plot(1, 0, 'ko', markersize=6, label='Root (x1)')
circle2 = plt.Circle((0, 0), R_cauchy, color='red', fill=False, linestyle='--', linewidth=1.5, label=f'Cauchy Scanner R={R_cauchy}')
ax2.add_patch(circle2)
ax2.set_xlim(-2.5, 2.5)
ax2.set_ylim(-2.5, 2.5)
ax2.set_aspect('equal')
ax2.set_title(r"Asymmetric Tree ($\beta(z) = \frac{256}{27} z^3(1-z)$)")
ax2.legend(loc='upper right')
ax2.grid(True, linestyle=':', alpha=0.6)

# --- PANEL 3 & 4: Pandrosion Scans and Spectral Energy ---
theta_s, ratio_s, r_hat_s = compute_pandrosion_scan(p_sym, R_cauchy, deg_sym)
theta_a, ratio_a, r_hat_a = compute_pandrosion_scan(p_asym, R_cauchy, deg_asym)

ax3 = axs[1, 0]
# Plot the magnitude of the continuous ratio along the circle
ax3.plot(theta_s, ratio_s, 'b-', linewidth=2, label='|r(θ)| on Cauchy Circle')
ax3.set_title("Symmetric Scan: Flat Ratio Envelope")
ax3.set_xlabel("Angle θ")
ax3.set_ylabel(r"Amplitude $|P(z)/P(a)|$")
# Set y-limits to show it's flat
mean_val = np.mean(ratio_s)
ax3.set_ylim(mean_val - 0.5, mean_val + 0.5)
ax3.grid(True, linestyle=':', alpha=0.6)
ax3.legend()

ax4 = axs[1, 1]
# For the asymmetric case, plot the spectral energy modes (bar chart)
k_modes = np.arange(1, deg_asym) # Skip DC component
energy_modes = r_hat_a[1:deg_asym]**2 # Spectral Energy per mode
ax4.bar(k_modes, energy_modes, color='red', alpha=0.7)
ax4.set_title(r"Asymmetric Spectrum: $\mathcal{E}(R)$ Signatures")
ax4.set_xlabel("Fourier Frequency Mode $k$")
ax4.set_ylabel(r"Spectral Energy $|\hat{r}_k|^2$")
ax4.set_xticks(k_modes)
ax4.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig("../../figures/spectral_dessins_plot.png", bbox_inches='tight', dpi=300)
print("Saved figure to ../../figures/spectral_dessins_plot.png")
