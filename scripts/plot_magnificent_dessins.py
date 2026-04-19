import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Use a professional style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 300,
    'font.family': 'serif',
    'axes.facecolor': '#111111'
})

def domain_coloring(w):
    """
    Computes an RGB image representing the complex field w.
    Phase maps to Hue. Magnitude maps to Lightness (dark at 0, bright at infinity).
    """
    # Hue: Phase of w
    H = (np.angle(w) + np.pi) / (2 * np.pi)
    
    # Value (brightness): Dark at zeros (roots), bright normally.
    # We want roots to be visibly black.
    mag = np.abs(w)
    V = 1.0 - 1.0 / (1.0 + mag**1.5)
    
    # Saturation: High generally, drop off slightly as it gets very large
    S = 1.0 * np.ones_like(H)
    
    HSV = np.dstack((H, S, V))
    RGB = mcolors.hsv_to_rgb(HSV)
    return RGB

def compute_pandrosion_scan(poly, max_r, d, N_samples=200):
    theta = np.linspace(0, 2*np.pi, N_samples)
    a_cont = max_r * np.exp(1j * theta)
    z_cont = max_r * np.exp(1j * (theta + np.pi/d))
    ratio_cont = np.abs(poly(z_cont) / poly(a_cont))
    
    s = np.arange(d)
    a_disc = max_r * np.exp(1j * 2 * np.pi * s / d)
    z_disc = max_r * np.exp(1j * (2 * np.pi * s / d + np.pi/d))
    r_vec = poly(z_disc) / poly(a_disc)
    
    r_hat = np.fft.fft(r_vec) / d
    return theta, ratio_cont, np.abs(r_hat)

fig, axs = plt.subplots(2, 2, figsize=(12, 9))

# --- Polynomial Definitions ---
poly_sym = lambda z: z**4
deg_sym = 4
roots_sym_0 = [0]
roots_sym_1 = np.roots([1, 0, 0, 0, -1]) # z^4 - 1 = 0

def poly_asym(z):
    return (256/27) * (z**3 - z**4)
deg_asym = 4
roots_asym_0 = [0, 1]
roots_asym_1 = np.roots([-256/27, 256/27, 0, 0, -1]) 

R_cauchy = 2.0
LIMIT = 2.5
x = np.linspace(-LIMIT, LIMIT, 800)
y = np.linspace(-LIMIT, LIMIT, 800)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y

# --- PANEL 1: Symmetric Phase Portrait ---
W_sym = poly_sym(Z)
RGB_sym = domain_coloring(W_sym)

ax1 = axs[0, 0]
ax1.imshow(RGB_sym, origin='lower', extent=[-LIMIT, LIMIT, -LIMIT, LIMIT])
ax1.set_title("Classical Star Dessin ($\\beta(z) = z^4$)")

# Plot roots of beta (black)
for r in roots_sym_0:
    ax1.plot(np.real(r), np.imag(r), 'ko', markersize=14, markeredgecolor='white', markeredgewidth=2)
# Plot roots of beta-1 (white)
for r in roots_sym_1:
    ax1.plot(np.real(r), np.imag(r), 'wo', markersize=10, markeredgecolor='black', markeredgewidth=1.5)

circle1 = plt.Circle((0, 0), R_cauchy, color='cyan', fill=False, linestyle='--', linewidth=2.5, label=f'Pandrosion Scanner R={R_cauchy}')
ax1.add_patch(circle1)
ax1.set_xlim(-LIMIT, LIMIT)
ax1.set_ylim(-LIMIT, LIMIT)

# --- PANEL 2: Asymmetric Phase Portrait ---
W_asym = poly_asym(Z)
RGB_asym = domain_coloring(W_asym)

ax2 = axs[0, 1]
ax2.imshow(RGB_asym, origin='lower', extent=[-LIMIT, LIMIT, -LIMIT, LIMIT])
ax2.set_title(r"Asymmetric Tree ($\beta(z) = \frac{256}{27} z^3(1-z)$)")

# Plot roots of beta (black)
for r in roots_asym_0:
    ax2.plot(np.real(r), np.imag(r), 'ko', markersize=14, markeredgecolor='white', markeredgewidth=2)
# Plot roots of beta-1 (white)
for r in roots_asym_1:
    ax2.plot(np.real(r), np.imag(r), 'wo', markersize=10, markeredgecolor='black', markeredgewidth=1.5)

circle2 = plt.Circle((0, 0), R_cauchy, color='cyan', fill=False, linestyle='--', linewidth=2.5)
ax2.add_patch(circle2)
ax2.set_xlim(-LIMIT, LIMIT)
ax2.set_ylim(-LIMIT, LIMIT)

# --- PANEL 3 & 4: Pandrosion Spectra (as before) ---
theta_s, ratio_s, r_hat_s = compute_pandrosion_scan(poly_sym, R_cauchy, deg_sym)
theta_a, ratio_a, r_hat_a = compute_pandrosion_scan(poly_asym, R_cauchy, deg_asym)

ax3 = axs[1, 0]
ax3.set_facecolor('#ffffff') # Explicitly white for graph readability
ax3.plot(theta_s, ratio_s, 'b-', linewidth=3)
ax3.set_title("Macroscopic Scan $|r(\\theta)|$ (Perfect Symmetry)")
ax3.set_xlabel("Angle θ on Cauchy Circle")
ax3.set_ylabel(r"Amplitude $|P(z)/P(a)|$")
mean_val = np.mean(ratio_s)
ax3.set_ylim(mean_val - 0.5, mean_val + 0.5)
ax3.grid(True, linestyle=':', alpha=0.6, color='black')

ax4 = axs[1, 1]
ax4.set_facecolor('#ffffff')
k_modes = np.arange(1, deg_asym) 
energy_modes = r_hat_a[1:deg_asym]**2 
ax4.bar(k_modes, energy_modes, color='#d62728', alpha=0.8, edgecolor='black', zorder=3)
ax4.set_title(r"Spectral Histogram $\mathcal{E}(R)$: Asymmetric Signature")
ax4.set_xlabel("Fourier Frequency Mode $k$")
ax4.set_ylabel(r"Spectral Energy $|\hat{r}_k|^2$")
ax4.set_xticks(k_modes)
ax4.grid(True, linestyle='-', alpha=0.3, zorder=0, color='black')

plt.tight_layout()
plt.savefig("../figures/magnificent_dessins.png", bbox_inches='tight', dpi=300)
print("Saved figure to ../figures/magnificent_dessins.png")
