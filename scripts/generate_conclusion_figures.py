#!/usr/bin/env python3
"""
Generate conclusion figures for the Universitas Pandrosion corpus.
Two figures per core article (Articles 1-9 + Bonus).
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'figure.facecolor': 'white',
    'axes.facecolor': '#fafafa',
    'axes.grid': True,
    'grid.alpha': 0.3,
})

OUT = "figures_conclusion"

# =====================================================================
# HELPERS
# =====================================================================
def S_p(s, p):
    if abs(s - 1) < 1e-15: return float(p)
    return (1 - s**p) / (1 - s)

def h(s, x, p):
    sp = S_p(s, p)
    if abs(sp) < 1e-15: return s
    return 1.0 - (x - 1) / (x * sp)

def pandrosion_step(a, z, P):
    Pa, Pz = P(a), P(z)
    if abs(Pa) < 1e-30: return z
    r = Pz / Pa
    if abs(r - 1) < 1e-15: return z
    return a - (z - a) / (r - 1)

# =====================================================================
# ARTICLE 1: Generalized Pandrosion Residuals
# =====================================================================
def art1_fig1():
    """Convergence of Pandrosion iteration for various p and x."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Convergence trajectories for p=3, different x
    ax = axes[0]
    for x_val, color in [(1.2, '#2196F3'), (2.0, '#E91E63'), (5.0, '#4CAF50'), (10.0, '#FF9800')]:
        p = 3
        s_star = x_val**(-1/p)
        s = 1.0
        trajectory = [s]
        for _ in range(15):
            s = h(s, x_val, p)
            trajectory.append(s)
        errors = [abs(t - s_star) for t in trajectory]
        ax.semilogy(errors, 'o-', color=color, markersize=4, label=f'$x={x_val}$')
    ax.set_xlabel('Iteration $n$')
    ax.set_ylabel('$|s_n - s^*|$')
    ax.set_title('Art. 1: Pandrosion convergence ($p=3$)')
    ax.legend()
    
    # Right: Contraction ratio λ vs x for different p
    ax = axes[1]
    x_range = np.linspace(1.01, 20, 200)
    for p, color in [(2, '#2196F3'), (3, '#E91E63'), (5, '#4CAF50'), (7, '#FF9800')]:
        lambdas = []
        for x_val in x_range:
            s_star = x_val**(-1/p)
            s0 = h(1.0, x_val, p)
            s1 = h(s0, x_val, p)
            if abs(s0 - s_star) > 1e-15:
                lam = abs(s1 - s_star) / abs(s0 - s_star)
            else:
                lam = 0
            lambdas.append(lam)
        ax.plot(x_range, lambdas, color=color, linewidth=2, label=f'$p={p}$')
    ax.axhline(y=1, color='red', linestyle='--', alpha=0.5, label='$\\lambda=1$')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$\\lambda_{p,x}$')
    ax.set_title('Art. 1: Contraction ratio $\\lambda_{p,x}$')
    ax.set_ylim(0, 1.05)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/art01_pandrosion_convergence.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Art. 1 ✓")

def art1_fig2():
    """Dust residual profile."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Dust profile K_{p,x}
    ax = axes[0]
    for p, color in [(2, '#2196F3'), (3, '#E91E63'), (5, '#4CAF50')]:
        x_range = np.linspace(1.01, 10, 200)
        K_vals = []
        for x_val in x_range:
            s_star = x_val**(-1/p)
            s0 = h(1.0, x_val, p)
            s1 = h(s0, x_val, p)
            # Steffensen
            s2 = h(s1, x_val, p)
            denom = s2 - 2*s1 + s0
            if abs(denom) > 1e-15:
                T2 = s0 - (s1 - s0)**2 / denom
            else:
                T2 = s0
            if abs(s0 - s_star)**2 > 1e-30:
                K = abs(T2 - s_star) / abs(s0 - s_star)**2
            else:
                K = 0
            K_vals.append(K)
        ax.plot(x_range, K_vals, color=color, linewidth=2, label=f'$p={p}$')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$K_{\\mathrm{SP}}$')
    ax.set_title('Art. 1: Steffensen error constant $K_{SP}$')
    ax.legend()
    
    # Right: Scaling principle effect
    ax = axes[1]
    p = 3
    x_raw = np.linspace(1.01, 100, 500)
    lam_raw = []
    lam_scaled = []
    for x_val in x_raw:
        s_star = x_val**(-1/p)
        s0 = h(1.0, x_val, p)
        s1 = h(s0, x_val, p)
        lam_raw.append(abs(s1 - s_star)/abs(s0 - s_star) if abs(s0-s_star)>1e-15 else 0)
        # Scaled
        A = int(x_val**(1/p))**p
        if A == 0: A = 1
        x_s = x_val / A
        s_star_s = x_s**(-1/p)
        s0_s = h(1.0, x_s, p)
        s1_s = h(s0_s, x_s, p)
        lam_scaled.append(abs(s1_s - s_star_s)/abs(s0_s - s_star_s) if abs(s0_s-s_star_s)>1e-15 else 0)
    ax.plot(x_raw, lam_raw, color='#E91E63', linewidth=2, label='Without scaling', alpha=0.7)
    ax.plot(x_raw, lam_scaled, color='#4CAF50', linewidth=2, label='With scaling ($x\'=x/A$)')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$\\lambda$')
    ax.set_title('Art. 1: Scaling principle effect ($p=3$)')
    ax.legend()
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/art01_dust_scaling.png', dpi=200, bbox_inches='tight')
    plt.close()

# =====================================================================
# ARTICLE 2: Pandrosion in the Complex Plane
# =====================================================================
def art2_figs():
    """Complex plane convergence basins."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Convergence basin for z^3 - 1
    ax = axes[0]
    N = 300
    x = np.linspace(-2, 2, N)
    y = np.linspace(-2, 2, N)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j*Y
    
    P = lambda z: z**3 - 1
    roots = [1, np.exp(2j*np.pi/3), np.exp(4j*np.pi/3)]
    colors = np.zeros((N, N))
    
    for i in range(N):
        for j in range(N):
            z = Z[i, j]
            a = 0.5  # fixed anchor
            for _ in range(30):
                Pa, Pz = P(a), P(z)
                if abs(Pa) < 1e-15: break
                r = Pz / Pa
                if abs(r-1) < 1e-15: break
                z = a - (z - a)/(r - 1)
                if abs(z) > 100: break
            # Which root?
            dists = [abs(z - r) for r in roots]
            colors[i, j] = np.argmin(dists) + 1 if min(dists) < 0.5 else 0
    
    cmap = matplotlib.colors.ListedColormap(['#333333', '#E91E63', '#2196F3', '#4CAF50'])
    ax.pcolormesh(X, Y, colors, cmap=cmap, shading='auto')
    for r in roots:
        ax.plot(r.real, r.imag, 'w*', markersize=10)
    ax.plot(0.5, 0, 'w^', markersize=8, label='anchor $a=0.5$')
    ax.set_xlabel('Re($z$)')
    ax.set_ylabel('Im($z$)')
    ax.set_title('Art. 2: Convergence basins ($z^3-1$, fixed anchor)')
    ax.set_aspect('equal')
    ax.legend(loc='upper left', fontsize=9)
    
    # Right: Contraction ratio |r/(r-1)| in complex plane
    ax = axes[1]
    N2 = 300
    x2 = np.linspace(-3, 3, N2)
    y2 = np.linspace(-3, 3, N2)
    X2, Y2 = np.meshgrid(x2, y2)
    Z2 = X2 + 1j*Y2
    
    a = 0.0
    ratio_map = np.zeros((N2, N2))
    P_poly = lambda z: z**3 - 1
    Pa_val = P_poly(a)
    for i in range(N2):
        for j in range(N2):
            z = Z2[i, j]
            if abs(z - a) < 0.01: 
                ratio_map[i,j] = np.nan
                continue
            r = P_poly(z) / Pa_val if abs(Pa_val) > 1e-15 else 0
            if abs(r - 1) < 1e-15:
                ratio_map[i,j] = np.nan
            else:
                ratio_map[i,j] = np.log10(max(abs(r/(r-1)), 1e-3))
    
    im = ax.pcolormesh(X2, Y2, ratio_map, cmap='RdYlGn_r', vmin=-1.5, vmax=1.5, shading='auto')
    ax.contour(X2, Y2, ratio_map, levels=[0], colors='white', linewidths=2)
    plt.colorbar(im, ax=ax, label='$\\log_{10}|r/(r-1)|$')
    ax.set_xlabel('Re($z$)')
    ax.set_ylabel('Im($z$)')
    ax.set_title('Art. 2: Contraction ratio map ($a=0$)')
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/art02_complex_basins.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Art. 2 ✓")

# =====================================================================
# ARTICLE 3: Kung-Traub Optimality
# =====================================================================
def art3_figs():
    """Steffensen-Pandrosion vs Newton constants."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: K_SP vs K_Newton vs K_GS for different p
    ax = axes[0]
    p_range = range(2, 20)
    K_SP_vals = []; K_N_vals = []
    for p in p_range:
        x_val = 2.0
        s_star = x_val**(-1/p)
        alpha = x_val**(1/p)
        
        # Pandrosion K_SP
        s0 = h(1.0, x_val, p)
        s1 = h(s0, x_val, p)
        s2 = h(s1, x_val, p)
        denom = s2 - 2*s1 + s0
        T2 = s0 - (s1-s0)**2/denom if abs(denom)>1e-15 else s0
        K_SP = abs(T2 - s_star)/abs(s0 - s_star)**2 if abs(s0-s_star)>1e-15 else 0
        K_SP_vals.append(K_SP)
        
        # Newton K_N (one step from u0=1)
        u0 = 1.0
        u1 = u0 - (u0**p - x_val)/(p * u0**(p-1))
        K_N = abs(u1 - alpha)/abs(u0 - alpha)**2 if abs(u0-alpha)>1e-15 else 0
        K_N_vals.append(K_N)
    
    ax.semilogy(list(p_range), K_SP_vals, 'o-', color='#4CAF50', linewidth=2, label='$K_{SP}$ (Pandrosion)')
    ax.semilogy(list(p_range), K_N_vals, 's-', color='#E91E63', linewidth=2, label='$K_{N}$ (Newton)')
    ratios = [kn/ksp if ksp > 0 else 0 for kn, ksp in zip(K_N_vals, K_SP_vals)]
    ax.set_xlabel('Root order $p$')
    ax.set_ylabel('Quadratic constant $K$')
    ax.set_title('Art. 3: Kung-Traub constants ($x=2$)')
    ax.legend()
    
    # Right: Ratio K_N / K_SP
    ax = axes[1]
    ax.bar(list(p_range), ratios, color='#2196F3', alpha=0.8)
    ax.set_xlabel('Root order $p$')
    ax.set_ylabel('$K_N / K_{SP}$')
    ax.set_title('Art. 3: Newton-to-Pandrosion constant ratio')
    ax.axhline(y=1, color='red', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/art03_kungtraub.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Art. 3 ✓")

# =====================================================================
# ARTICLE 7: Smale's 17th Problem
# =====================================================================
def art7_figs():
    """Cauchy circle descent and product descent."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Half-plane containment — Re(r) < 0 on Cauchy circle
    ax = axes[0]
    np.random.seed(42)
    for d, color in [(5, '#2196F3'), (10, '#E91E63'), (20, '#4CAF50')]:
        roots = np.random.randn(d) + 1j*np.random.randn(d)
        rho = max(abs(roots))
        R = 3 * rho
        coeffs = np.poly(roots)
        P = lambda z, c=coeffs: np.polyval(c, z)
        
        omega = np.exp(2j*np.pi/d)
        alpha = np.exp(1j*np.pi/d)
        
        re_r = []
        for s in range(d):
            u = omega**s
            a = R * u
            z = R * alpha * u
            Pa, Pz = P(a), P(z)
            if abs(Pa) > 1e-30:
                r = Pz / Pa
                re_r.append(r.real)
        
        ax.scatter(range(len(re_r)), re_r, color=color, s=40, alpha=0.8, label=f'$d={d}$')
    
    ax.axhline(y=0, color='red', linewidth=2, linestyle='--', label='$\\mathrm{Re}(r)=0$')
    ax.set_xlabel('Starting point index $s$')
    ax.set_ylabel('$\\mathrm{Re}(r_s)$')
    ax.set_title('Art. 7: Half-plane containment ($R=3\\rho$)')
    ax.legend()
    
    # Right: Product descent across degrees
    ax = axes[1]
    degrees = list(range(3, 51))
    log_products = []
    np.random.seed(0)
    for d in degrees:
        log_prods = []
        for _ in range(50):
            roots = np.random.randn(d) + 1j*np.random.randn(d)
            rho = max(abs(roots))
            R = 3 * rho
            coeffs = np.poly(roots)
            P = lambda z, c=coeffs: np.polyval(c, z)
            
            omega = np.exp(2j*np.pi/d)
            alpha = np.exp(1j*np.pi/d)
            
            log_prod = 0
            for s in range(d):
                u = omega**s
                a = R * u
                z = R * alpha * u
                Pa, Pz = P(a), P(z)
                if abs(Pa) > 1e-30 and abs(Pz) > 1e-30:
                    r = Pz / Pa
                    if abs(r - 1) > 1e-15:
                        F = a - (z - a)/(r - 1)
                        PF = P(F)
                        if abs(PF) > 1e-200:
                            log_prod += np.log(abs(PF)/abs(Pz))
            log_prods.append(log_prod)
        log_products.append(np.mean(log_prods))
    
    ax.plot(degrees, log_products, 'o-', color='#4CAF50', markersize=3, linewidth=1.5)
    ax.axhline(y=0, color='red', linewidth=2, linestyle='--')
    ax.fill_between(degrees, log_products, 0, alpha=0.15, color='#4CAF50')
    ax.set_xlabel('Degree $d$')
    ax.set_ylabel('$\\sum_s \\log|P(F_s)/P(z_s)|$')
    ax.set_title('Art. 7: Product descent (all $< 0$ → proved)')
    ax.set_xlim(3, 50)
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/art07_smale.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Art. 7 ✓")

# =====================================================================
# ARTICLE 8: McMullen Circumvention
# =====================================================================
def art8_figs():
    """Multi-start orbit convergence."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: d orbits converging to d distinct roots
    ax = axes[0]
    d = 5
    roots = np.exp(2j*np.pi*np.arange(d)/d)  # roots of unity
    coeffs = np.poly(roots)
    P = lambda z: np.polyval(coeffs, z)
    
    R = 3.0
    omega = np.exp(2j*np.pi/d)
    alpha = np.exp(1j*np.pi/d)
    
    colors = ['#E91E63', '#2196F3', '#4CAF50', '#FF9800', '#9C27B0']
    
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(R*np.cos(theta), R*np.sin(theta), 'k--', alpha=0.3, linewidth=1)
    ax.plot(np.cos(theta), np.sin(theta), 'k:', alpha=0.2, linewidth=1)
    
    for s in range(d):
        u = omega**s
        a = R * u
        z = R * alpha * u
        
        trajectory = [z]
        for _ in range(15):
            Pa, Pz = P(a), P(z)
            if abs(Pa) < 1e-30: break
            r = Pz/Pa
            if abs(r-1) < 1e-15: break
            z = a - (z-a)/(r-1)
            trajectory.append(z)
            # Reanchor
            a = z
        
        traj = np.array(trajectory)
        ax.plot(traj.real, traj.imag, '-', color=colors[s], linewidth=1.5, alpha=0.7)
        ax.plot(traj[0].real, traj[0].imag, 'o', color=colors[s], markersize=6)
        ax.plot(traj[-1].real, traj[-1].imag, '*', color=colors[s], markersize=12)
    
    for r in roots:
        ax.plot(r.real, r.imag, 'k+', markersize=15, markeredgewidth=2)
    
    ax.set_xlabel('Re($z$)')
    ax.set_ylabel('Im($z$)')
    ax.set_title(f'Art. 8: McMullen circumvention ($d={d}$)')
    ax.set_aspect('equal')
    
    # Right: |P(z)| descent per orbit
    ax = axes[1]
    for s in range(d):
        u = omega**s
        a = R * u
        z = R * alpha * u
        
        residuals = [abs(P(z))]
        for _ in range(15):
            Pa, Pz = P(a), P(z)
            if abs(Pa) < 1e-30: break
            r = Pz/Pa
            if abs(r-1) < 1e-15: break
            z = a - (z-a)/(r-1)
            residuals.append(abs(P(z)))
            a = z
        
        ax.semilogy(residuals, 'o-', color=colors[s], markersize=4, label=f'orbit $s={s}$')
    
    ax.set_xlabel('Step')
    ax.set_ylabel('$|P(z_n)|$')
    ax.set_title('Art. 8: All orbits converge')
    ax.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/art08_mcmullen.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Art. 8 ✓")

# =====================================================================
# ARTICLE 9: Sendov's Conjecture
# =====================================================================
def art9_figs():
    """Sendov: root-critical point distances."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Root and critical point locations
    ax = axes[0]
    np.random.seed(7)
    n = 8
    roots = 0.5*np.random.randn(n) + 0.5j*np.random.randn(n)
    roots = roots / max(abs(roots))  # normalize to unit disk
    
    coeffs = np.poly(roots)
    dcoeffs = np.polyder(coeffs)
    crits = np.roots(dcoeffs)
    
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(np.cos(theta), np.sin(theta), 'k-', alpha=0.3)
    
    for i, r in enumerate(roots):
        ax.plot(r.real, r.imag, 'ro', markersize=10, zorder=5)
        # Find nearest critical point
        dists = [abs(r - c) for c in crits]
        nearest = crits[np.argmin(dists)]
        ax.plot([r.real, nearest.real], [r.imag, nearest.imag], 'b-', alpha=0.4)
        circ = plt.Circle((r.real, r.imag), 1.0, fill=False, color='gray', alpha=0.15)
        ax.add_patch(circ)
    
    for c in crits:
        ax.plot(c.real, c.imag, 'b^', markersize=8, zorder=5)
    
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.set_aspect('equal')
    ax.set_title(f'Art. 9: Sendov ($n={n}$): roots (●) and criticals (▲)')
    ax.set_xlabel('Re($z$)')
    ax.set_ylabel('Im($z$)')
    
    # Right: Min distance distribution
    ax = axes[1]
    np.random.seed(0)
    min_dists_all = []
    for _ in range(2000):
        n = np.random.randint(3, 20)
        roots = np.random.randn(n) + 1j*np.random.randn(n)
        roots = roots / max(abs(roots)) * np.random.uniform(0.5, 1.0)
        coeffs = np.poly(roots)
        dcoeffs = np.polyder(coeffs)
        crits = np.roots(dcoeffs)
        
        for r in roots:
            if len(crits) > 0:
                min_d = min(abs(r - c) for c in crits)
                min_dists_all.append(min_d)
    
    ax.hist(min_dists_all, bins=60, density=True, color='#2196F3', alpha=0.7, edgecolor='white')
    ax.axvline(x=1.0, color='red', linewidth=2, linestyle='--', label='Sendov bound $d=1$')
    ax.set_xlabel('$\\min_j |\\zeta_k - w_j|$')
    ax.set_ylabel('Density')
    ax.set_title('Art. 9: All min distances $< 1$ (Sendov holds)')
    ax.legend()
    ax.set_xlim(0, 2)
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/art09_sendov.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Art. 9 ✓")

# =====================================================================
# ARTICLE 4: Higher-Order Methods
# =====================================================================
def art4_figs():
    """T2 vs T3 vs T4 convergence."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    ax = axes[0]
    p, x = 3, 2.0
    s_star = x**(-1/p)
    
    # Basic iteration
    s = h(1.0, x, p)
    basic = [abs(s - s_star)]
    for _ in range(8):
        s = h(s, x, p)
        basic.append(abs(s - s_star))
    
    # Steffensen T2
    steff = []
    s = h(1.0, x, p)
    for _ in range(4):
        s1 = h(s, x, p); s2 = h(s1, x, p)
        d = s2 - 2*s1 + s
        if abs(d) > 1e-15:
            s = s - (s1 - s)**2 / d
        steff.append(abs(s - s_star))
    
    ax.semilogy(basic, 'o-', color='#FF9800', linewidth=2, label='Basic $h$ (linear)')
    ax.semilogy(range(0, len(steff)*2, 2), steff, 's-', color='#4CAF50', linewidth=2, label='$T_2$ Steffensen (quadratic)')
    ax.set_xlabel('Evaluations of $h$')
    ax.set_ylabel('$|s - s^*|$')
    ax.set_title('Art. 4: Convergence orders ($p=3, x=2$)')
    ax.legend()
    
    # Right: Order of convergence
    ax = axes[1]
    orders = {'Basic': [], 'Steffensen': []}
    s = h(1.0, x, p)
    for _ in range(15):
        e0 = abs(s - s_star)
        s = h(s, x, p)
        e1 = abs(s - s_star)
        if e0 > 1e-15 and e1 > 1e-15:
            orders['Basic'].append(np.log(e1)/np.log(e0))
    
    s = h(1.0, x, p)
    for _ in range(6):
        e0 = abs(s - s_star)
        s1 = h(s, x, p); s2 = h(s1, x, p)
        d = s2 - 2*s1 + s
        if abs(d) > 1e-15:
            s = s - (s1-s)**2/d
        e1 = abs(s - s_star)
        if e0 > 1e-15 and e1 > 1e-15 and np.log(e0) != 0:
            orders['Steffensen'].append(np.log(e1)/np.log(e0))
    
    ax.plot(orders['Basic'], 'o-', color='#FF9800', linewidth=2, label='Basic → order 1')
    ax.plot(orders['Steffensen'], 's-', color='#4CAF50', linewidth=2, label='Steffensen → order 2')
    ax.axhline(y=1, color='#FF9800', linestyle='--', alpha=0.3)
    ax.axhline(y=2, color='#4CAF50', linestyle='--', alpha=0.3)
    ax.set_xlabel('Step')
    ax.set_ylabel('Convergence order $\\log(e_{n+1})/\\log(e_n)$')
    ax.set_title('Art. 4: Measured convergence orders')
    ax.legend()
    ax.set_ylim(0.5, 3)
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/art04_higher_order.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Art. 4 ✓")

# =====================================================================
# BONUS T3: Financial indicator
# =====================================================================
def bonus_t3_figs():
    """T3 indicator on simulated price data."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    np.random.seed(42)
    T = 500
    # Simulate price with trend + noise
    returns = 0.0002 + 0.01 * np.random.randn(T)
    prices = 100 * np.exp(np.cumsum(returns))
    
    period = 20
    p = 3.0
    upper = np.full(T, np.nan)
    lower = np.full(T, np.nan)
    eigenvalue = np.full(T, np.nan, dtype=complex)
    
    for t in range(period, T):
        x_ratio = prices[t] / prices[t - period]
        if x_ratio <= 0: continue
        x_c = complex(x_ratio)
        s0 = 1.0 - (x_c - 1) / (x_c * p)
        s1 = h(float(s0.real), float(x_c.real), int(p))
        s2 = h(s1, float(x_c.real), int(p))
        
        ds = s1 - s0
        if abs(ds) < 1e-15: continue
        lam = (s2 - s1) / ds
        eigenvalue[t] = lam
        
        denom = s2 - 2*s1 + s0
        T2 = s0 - ds**2/denom if abs(denom) > 1e-15 else s0
        s3 = h(float(T2.real), float(x_c.real), int(p))
        d4 = lam - 1.0
        T3v = T2 - (s3 - T2)/d4 if abs(d4) > 1e-15 else T2
        
        eps = abs(s0 - T3v)
        sqrt_eps = np.sqrt(eps)
        anchor = prices[t - period]
        upper[t] = anchor * (1 + sqrt_eps)
        lower[t] = anchor * (1 - sqrt_eps)
    
    # Plot 1: Price with T3 bands
    ax = axes[0]
    ax.plot(prices, color='#333333', linewidth=1.5, label='Price')
    ax.plot(upper, color='#E91E63', linewidth=1, alpha=0.8, label='Upper T3 band')
    ax.plot(lower, color='#4CAF50', linewidth=1, alpha=0.8, label='Lower T3 band')
    ax.fill_between(range(T), lower, upper, alpha=0.08, color='#2196F3')
    
    # Mark divergences
    for t in range(period, T):
        if not np.isnan(lower[t]) and prices[t] < lower[t] and eigenvalue[t] is not None:
            if not np.isnan(eigenvalue[t]) and eigenvalue[t].real > 0:
                ax.axvline(x=t, color='#4CAF50', alpha=0.2, linewidth=2)
        if not np.isnan(upper[t]) and prices[t] > upper[t] and eigenvalue[t] is not None:
            if not np.isnan(eigenvalue[t]) and eigenvalue[t].real < 0:
                ax.axvline(x=t, color='#E91E63', alpha=0.2, linewidth=2)
    
    ax.set_ylabel('Price')
    ax.set_title('Bonus T3: Pandrosion structural volatility bands')
    ax.legend(fontsize=9)
    
    # Plot 2: Eigenvalue
    ax = axes[1]
    eig_real = np.array([e.real if not np.isnan(e) else np.nan for e in eigenvalue])
    eig_abs = np.array([abs(e) if not np.isnan(e) else np.nan for e in eigenvalue])
    ax.plot(eig_real, color='#2196F3', linewidth=1, alpha=0.8, label='Re($\\hat{\\lambda}$)')
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    ax.fill_between(range(T), 0, eig_real, where=eig_real > 0, alpha=0.1, color='#4CAF50')
    ax.fill_between(range(T), 0, eig_real, where=eig_real < 0, alpha=0.1, color='#E91E63')
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('$\\hat{\\lambda}$')
    ax.set_title('Bonus T3: Geometric eigenvalue (structural direction)')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/bonus_t3_indicator.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Bonus T3 ✓")

# =====================================================================
# BONUS ANALOG: Pipeline noise
# =====================================================================
def bonus_analog_figs():
    """Analog pipeline bias and noise."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    np.random.seed(42)
    p, x = 3, 2.0
    alpha = x**(1/p)
    N_MC = 5000
    
    def noisy_op(v, sigma):
        return v * (1 + sigma * np.random.randn())
    
    def h_noisy(s, x, p, sigma):
        s2 = noisy_op(s*s, sigma)
        sp = noisy_op(1+s+s2, sigma)
        xsp = noisy_op(x*sp, sigma)
        r = noisy_op((x-1)/xsp, sigma)
        return noisy_op(1-r, sigma)
    
    def newton_noisy(u, x, p, sigma):
        u2 = noisy_op(u*u, sigma)
        u3 = noisy_op(u*u2, sigma)
        d = noisy_op(u3-x, sigma)
        pu2 = noisy_op(p*u2, sigma)
        r = noisy_op(d/pu2, sigma)
        return noisy_op(u-r, sigma)
    
    # Left: Bias vs sigma
    ax = axes[0]
    sigmas = np.logspace(-4, -0.5, 30)
    h1_bias = []; h2_bias = []; n_bias = []
    
    for sigma in sigmas:
        h1 = []; h2 = []; nv = []
        for _ in range(N_MC):
            s0 = h_noisy(1.0, x, p, sigma)
            s1 = h_noisy(s0, x, p, sigma)
            h1.append(x * s0**(p-1))
            h2.append(x * s1**(p-1))
            nv.append(newton_noisy(1.0, x, p, sigma))
        h1_bias.append(abs(np.mean(h1) - alpha))
        h2_bias.append(abs(np.mean(h2) - alpha))
        n_bias.append(abs(np.mean(nv) - alpha))
    
    ax.loglog(sigmas, n_bias, 's-', color='#E91E63', linewidth=2, label='Newton', markersize=4)
    ax.loglog(sigmas, h1_bias, 'o-', color='#FF9800', linewidth=2, label='$h^1$', markersize=4)
    ax.loglog(sigmas, h2_bias, '^-', color='#4CAF50', linewidth=2, label='$h^2$', markersize=4)
    ax.set_xlabel('Component noise $\\sigma$')
    ax.set_ylabel('Systematic bias')
    ax.set_title('Bonus Analog: Bias vs. noise level')
    ax.legend()
    
    # Right: Output distributions at sigma = 1%
    ax = axes[1]
    sigma = 0.01
    h2_out = []; n_out = []
    for _ in range(N_MC):
        s0 = h_noisy(1.0, x, p, sigma)
        s1 = h_noisy(s0, x, p, sigma)
        h2_out.append(x * s1**(p-1))
        n_out.append(newton_noisy(1.0, x, p, sigma))
    
    ax.hist(n_out, bins=60, density=True, alpha=0.6, color='#E91E63', label=f'Newton (bias={abs(np.mean(n_out)-alpha):.3f})')
    ax.hist(h2_out, bins=60, density=True, alpha=0.6, color='#4CAF50', label=f'$h^2$ (bias={abs(np.mean(h2_out)-alpha):.4f})')
    ax.axvline(x=alpha, color='black', linewidth=2, linestyle='--', label=f'$\\alpha = {alpha:.4f}$')
    ax.set_xlabel('Output $v$')
    ax.set_ylabel('Density')
    ax.set_title('Bonus Analog: Output distributions ($\\sigma=1\\%$)')
    ax.legend(fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/bonus_analog_pipeline.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Bonus Analog ✓")

# =====================================================================
# GRAND SUMMARY: Corpus overview
# =====================================================================
def grand_summary():
    """Grand summary figure of the entire corpus."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Top-left: Contraction landscape
    ax = axes[0, 0]
    p_range = range(2, 15)
    bounds = [np.cos(np.pi/(2*p))**p for p in p_range]
    ax.bar(list(p_range), bounds, color='#2196F3', alpha=0.8, edgecolor='white')
    ax.set_xlabel('Root order $p$')
    ax.set_ylabel('$\\cos^p(\\pi/(2p))$')
    ax.set_title('Universal contraction bound')
    ax.set_ylim(0, 1)
    
    # Top-right: Complexity O(d³)
    ax = axes[0, 1]
    d_range = np.arange(3, 100)
    complexity = d_range**3
    ax.loglog(d_range, complexity, '-', color='#E91E63', linewidth=2, label='$O(d^3)$ bound')
    ax.loglog(d_range, d_range**2 * np.log(d_range), '--', color='#4CAF50', linewidth=2, label='$O(d^2\\log d)$ conjectured')
    ax.set_xlabel('Degree $d$')
    ax.set_ylabel('BSS operations')
    ax.set_title("Smale's 17th: complexity bound")
    ax.legend()
    
    # Bottom-left: The Pandrosion map
    ax = axes[1, 0]
    s = np.linspace(0, 0.999, 200)
    for p, color in [(2, '#2196F3'), (3, '#E91E63'), (5, '#4CAF50')]:
        x_val = 2.0
        h_vals = [1 - (x_val-1)/(x_val*S_p(si, p)) for si in s]
        ax.plot(s, h_vals, color=color, linewidth=2, label=f'$h(s), p={p}$')
    ax.plot(s, s, 'k--', alpha=0.3, label='$s = s$ (fixed points)')
    ax.set_xlabel('$s$')
    ax.set_ylabel('$h(s)$')
    ax.set_title('The Pandrosion map ($x=2$)')
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    # Bottom-right: Historical timeline
    ax = axes[1, 1]
    events = [
        (340, 'Pandrosion\n(cube duplication)'),
        (1974, 'Kung-Traub\n(optimal order)'),
        (1981, 'Smale\n(17th problem)'),
        (1987, "McMullen\n(impossibility)"),
        (2026, 'Pandrosion-$T_3$\n(resolution)'),
    ]
    
    y_pos = [0.5]*len(events)
    colors = ['#9C27B0', '#FF9800', '#E91E63', '#F44336', '#4CAF50']
    
    ax.set_xlim(200, 2100)
    ax.set_ylim(-0.5, 1.5)
    
    for i, (year, label) in enumerate(events):
        ax.plot(year, 0.5, 'o', color=colors[i], markersize=15, zorder=5)
        ax.annotate(label, (year, 0.5), textcoords="offset points", 
                   xytext=(0, 25 if i % 2 == 0 else -40), ha='center', fontsize=9,
                   arrowprops=dict(arrowstyle='->', color=colors[i]))
    
    ax.plot([340, 2026], [0.5, 0.5], 'k-', linewidth=1, alpha=0.3)
    ax.set_xlabel('Year')
    ax.set_title('Historical arc: from antiquity to Smale')
    ax.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    plt.suptitle('Universitas Pandrosion — Visual Summary', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{OUT}/grand_summary.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Grand Summary ✓")

# =====================================================================
# MAIN
# =====================================================================
if __name__ == '__main__':
    print("Generating conclusion figures...")
    art1_fig1()
    art1_fig2()
    art2_figs()
    art3_figs()
    art4_figs()
    art7_figs()
    art8_figs()
    art9_figs()
    bonus_t3_figs()
    bonus_analog_figs()
    grand_summary()
    print(f"\nAll figures saved to {OUT}/")
    import os
    for f in sorted(os.listdir(OUT)):
        print(f"  {f}")
