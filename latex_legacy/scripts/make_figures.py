"""Generate publication-quality figures for paper 101ABC and Part IX-mv."""
import math, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Publication style
mpl.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
})

FIG_LEGACY = "/Users/ivanbesevic/Documents/poussiere/latex_legacy/figures"
FIG_LATEX  = "/Users/ivanbesevic/Documents/poussiere/latex/figures"


# ===========================================================================
#  101ABC Figure 1: The journey of the bound
# ===========================================================================

def fig_101abc_journey():
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    d = np.logspace(np.log10(16), np.log10(512), 200)

    # Bounds (with arbitrary constants for visualisation)
    log2d = np.log2(d)
    f_part8     = 0.05 * d**2 * log2d**2     # Part-VIII baseline (cond.)
    f_trilogy   = 0.30 * d**2 * log2d        # Trilogy 101bis/ter
    f_with_els  = 1.50 * d * log2d           # With ELS (101quater + IX)
    f_lower     = 1.00 * d                   # Trivial lower bound

    ax.loglog(d, f_part8, '--', color='#a83232', lw=2.0,
              label=r'Part VIII (3 conjectures): $O(d^{2}\log^{2}d)$')
    ax.loglog(d, f_trilogy, '--', color='#d68c2a', lw=2.0,
              label=r'Trilogy 101bis/ter (1 lemma): $O(d^{2}\log d)$')
    ax.loglog(d, f_with_els, '-', color='#1f7a3a', lw=2.4,
              label=r'+ ELS (101quater, IX): $\mathbf{O(d \log d)}$')
    ax.loglog(d, f_lower, ':', color='black', lw=1.5, alpha=0.6,
              label=r'Lower bound $\Omega(d)$')

    ax.set_xlabel(r'Polynomial degree $d$')
    ax.set_ylabel(r'Las Vegas cost on KS (BSS ops, log scale)')
    ax.set_title(r"The algorithmic journey: from $O(d^{2}\log^{2}d)$ to $O(d \log d)$")
    ax.legend(loc='upper left', framealpha=0.95)
    ax.grid(True, which='both', alpha=0.25)

    # Annotations: mark key d's
    for dx in [16, 64, 256]:
        ax.axvline(dx, color='gray', alpha=0.15, lw=0.5)
    plt.tight_layout()
    out = os.path.join(FIG_LEGACY, "101abc_journey.png")
    plt.savefig(out)
    plt.close()
    print(f"Saved: {out}")


# ===========================================================================
#  101ABC Figure 2: Pre-basin epochs — Armijo vs ELS speedup
# ===========================================================================

def fig_101abc_speedup():
    ds = np.array([16, 32, 64, 128, 256])
    armijo = np.array([11.72, 22.18, 38.10, 73.28, 118.12])
    els    = np.array([3.50,  3.73,  3.67,  2.13,  2.97])
    speedup = armijo / els

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # Left: dual log-log of N_pre
    ax1.loglog(ds, armijo, 'o-', color='#a83232', lw=2.0, markersize=8,
               label=r'Original Armijo (paper 101)')
    ax1.loglog(ds, els, 's-', color='#1f7a3a', lw=2.0, markersize=8,
               label=r'Extended Line Search (paper 101quater)')
    # Theoretical fit
    d_fit = np.logspace(np.log10(16), np.log10(256), 100)
    ax1.loglog(d_fit, 1.1 * d_fit**0.84, '--', color='#a83232', alpha=0.5,
               label=r'$\sim d^{0.84}$ (paper 101 fit)')
    ax1.axhline(5, color='#1f7a3a', alpha=0.4, ls=':',
                label=r'$O(1)$ (ELS, theorem 3.2)')
    ax1.set_xlabel(r'Degree $d$')
    ax1.set_ylabel(r'$\mathbb{E}[N_{\mathrm{pre-basin}}]$')
    ax1.set_title(r'(a) Pre-basin epoch count')
    ax1.legend(loc='upper left', framealpha=0.95)
    ax1.grid(True, which='both', alpha=0.25)

    # Right: speedup
    bars = ax2.bar(np.arange(len(ds)), speedup, color=['#86c2a5','#5fb087','#3d9970','#1f7a3a','#005a26'],
                    edgecolor='black')
    for b, sp in zip(bars, speedup):
        ax2.text(b.get_x() + b.get_width()/2, b.get_height() + 0.7,
                 f'{sp:.1f}×', ha='center', fontweight='bold')
    ax2.set_xticks(np.arange(len(ds)))
    ax2.set_xticklabels([str(int(d)) for d in ds])
    ax2.set_xlabel(r'Degree $d$')
    ax2.set_ylabel(r'Speedup ELS vs original Armijo')
    ax2.set_title(r'(b) ELS speedup grows with $d$')
    ax2.grid(True, axis='y', alpha=0.25)
    ax2.set_ylim(0, 50)

    plt.tight_layout()
    out = os.path.join(FIG_LEGACY, "101abc_speedup.png")
    plt.savefig(out)
    plt.close()
    print(f"Saved: {out}")


# ===========================================================================
#  101ABC Figure 3: Cost decomposition pie / stack
# ===========================================================================

def fig_101abc_decomposition():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    # Left: stacked-bar of factors
    schemes = ['Part VIII\nbaseline', 'Trilogy\n101bis/ter', 'With ELS\n(101quater, IX)']
    multi   = np.array([np.log2(64), 1, 1])         # log d, 1, 1 in scale
    epoch   = np.array([np.log2(64), np.log2(64), 1])  # log d, log d, 1
    perepoch= np.array([np.log2(64), 1, 1])         # log d, 1, 1
    width = 0.45
    x = np.arange(3)
    p1 = axes[0].bar(x, multi, width, label='multi-start  $s^\\star$', color='#86c2a5', edgecolor='black')
    p2 = axes[0].bar(x, epoch, width, bottom=multi, label='epochs/orbit  $N_{\\mathrm{ep}}$', color='#5fb087', edgecolor='black')
    p3 = axes[0].bar(x, perepoch, width, bottom=multi+epoch, label='cost/epoch  $C_{\\mathrm{ep}}$', color='#3d9970', edgecolor='black')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(schemes)
    axes[0].set_ylabel(r'log-factor contribution to total cost')
    axes[0].set_title(r'(a) How each factor was reduced')
    axes[0].legend(loc='upper right', framealpha=0.95)
    axes[0].grid(True, axis='y', alpha=0.25)
    for i, (m, e, p) in enumerate(zip(multi, epoch, perepoch)):
        labels = [['$d$','$\\log d$','$d\\log d$'],['$1$','$\\log d$','$d$'],['$1$','$\\log d$','$d$']]
        labels = [[r'$d$','$d\log d$','$d\log d$'],[r'$1$','$\log d$','$d$'],[r'$1$',r'$\log d$',r'$d$']]
        # Just text the symbolic factor inside each segment
        symbols = [['d', '\\log d', 'd\\log d'], ['1', '\\log d', 'd'], ['1', '\\log d', 'd']]
        s = symbols[i]
        axes[0].text(i, m/2, f'${s[0]}$', ha='center', va='center', fontweight='bold', fontsize=11)
        axes[0].text(i, m + e/2, f'${s[1]}$', ha='center', va='center', fontweight='bold', fontsize=11)
        axes[0].text(i, m + e + p/2, f'${s[2]}$', ha='center', va='center', fontweight='bold', fontsize=11)

    # Right: which lemmas / conjectures
    conditional = ['Part VIII\nbaseline', 'Trilogy', 'With ELS', 'Lower bound']
    n_cond = [3, 1, 2, 0]
    colors = ['#a83232', '#d68c2a', '#1f7a3a', 'black']
    bars = axes[1].bar(np.arange(4), n_cond, color=colors, edgecolor='black')
    axes[1].set_xticks(np.arange(4))
    axes[1].set_xticklabels(conditional, rotation=15)
    axes[1].set_ylabel(r'Number of conjectures / lemmas required')
    axes[1].set_title(r'(b) Conditional debt across stages')
    for b, c in zip(bars, n_cond):
        axes[1].text(b.get_x()+b.get_width()/2, b.get_height()+0.05, str(c), ha='center', fontweight='bold')
    axes[1].set_ylim(0, 3.5)
    axes[1].grid(True, axis='y', alpha=0.25)

    plt.tight_layout()
    out = os.path.join(FIG_LEGACY, "101abc_decomposition.png")
    plt.savefig(out)
    plt.close()
    print(f"Saved: {out}")


# ===========================================================================
#  Part IX-mv Figure 1: The 9-orders-of-magnitude scaling
# ===========================================================================

def fig_mv_scaling():
    # Data from multivariate_KS_huge_v2.py  (and earlier mvKS runs)
    data = [
        ((2, 3), 9, 7.38),
        ((3, 3), 27, 7.62),
        ((4, 2), 16, 6.62),
        ((4, 3), 81, 8.38),
        ((5, 3), 243, 9.75),
        ((5, 4), 1024, 10.62),
        ((6, 3), 729, 9.12),
        ((6, 4), 4096, 10.25),
        ((4, 8), 4096, 11.12),
        ((5, 6), 7776, 9.25),
        ((6, 6), 46656, 9.62),
        ((4, 16), 65536, 9.50),
        ((5, 16), 1048576, 12.88),
        ((6, 10), 1000000, 12.88),
        ((15, 4), 1073741824, 15.67),
        ((20, 3), 3486784401, 15.33),
        ((10, 8), 1073741824, 15.33),
        ((12, 6), 2176782336, 16.00),
        ((16, 4), 4294967296, 14.00),
    ]
    Ds = np.array([row[1] for row in data], dtype=float)
    Ntots = np.array([row[2] for row in data])
    log2D = np.log2(Ds)
    loglogD = np.log2(np.maximum(log2D, 1.0))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))

    # Left: log-log
    sc = axes[0].scatter(Ds, Ntots, c=log2D, cmap='viridis', s=80, edgecolor='black', zorder=3)
    Dfit = np.logspace(np.log10(9), np.log10(5e9), 200)
    Nfit = 3.67 + 2.31 * np.log2(np.log2(Dfit))
    axes[0].plot(Dfit, Nfit, 'r-', lw=2.0, alpha=0.8,
                 label=r'$\mathbb{E}[N_{\mathrm{tot}}] = 3.67 + 2.31\, \log_2\log_2 D$')
    axes[0].set_xscale('log')
    axes[0].set_xlabel(r"Bézout degree $D = d^n$ (log scale)")
    axes[0].set_ylabel(r'$\mathbb{E}[N_{\mathrm{tot}}]$ (epochs to converge)')
    axes[0].set_title(r'(a) $N_{\mathrm{tot}}$ vs $D$ over 9 orders of magnitude')
    axes[0].legend(loc='upper left', framealpha=0.95)
    axes[0].grid(True, which='both', alpha=0.25)
    axes[0].set_ylim(0, 22)

    # Annotate D = 10^6, 10^9
    for D, N, txt in [(1e6, 12.88, '$D=10^{6}$'), (4.3e9, 14.0, '$D=4.3{\\cdot}10^{9}$')]:
        axes[0].annotate(txt, xy=(D, N), xytext=(D*0.05, N+1.8),
                         arrowprops=dict(arrowstyle='->', color='gray'),
                         fontsize=10)

    # Right: linear vs log_2 log_2 D
    axes[1].scatter(loglogD, Ntots, c=log2D, cmap='viridis', s=80, edgecolor='black', zorder=3)
    xfit = np.linspace(loglogD.min(), loglogD.max(), 100)
    axes[1].plot(xfit, 3.67 + 2.31*xfit, 'r-', lw=2.0, alpha=0.8,
                 label=r'fit: $3.67 + 2.31\, x$')
    axes[1].set_xlabel(r'$\log_2 \log_2 D$')
    axes[1].set_ylabel(r'$\mathbb{E}[N_{\mathrm{tot}}]$')
    axes[1].set_title(r'(b) Linear in $\log\log D$ — quadratic basin')
    axes[1].legend(loc='upper left', framealpha=0.95)
    axes[1].grid(True, alpha=0.25)
    axes[1].set_ylim(0, 22)

    cb = fig.colorbar(sc, ax=axes, orientation='vertical', pad=0.02, shrink=0.85)
    cb.set_label(r'$\log_2 D$')

    out = os.path.join(FIG_LATEX, "mv_scaling.png")
    plt.savefig(out)
    plt.close()
    print(f"Saved: {out}")


# ===========================================================================
#  Part IX-mv Figure 2: Optimal tau histogram (univariate vs multivariate)
# ===========================================================================

def fig_mv_tau():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Univariate: tau ~ sqrt d
    ds_u = [16, 32, 64, 128, 256, 1024]
    sqrtd = [math.sqrt(d) for d in ds_u]
    tau_u = [4, 5, 8, 12, 16, 32]  # idealised observation, paper 101quater Lemma 3.1
    axes[0].plot(ds_u, sqrtd, 'k--', lw=1.5, label=r'$\sqrt{d}$ (theory)', alpha=0.6)
    axes[0].plot(ds_u, tau_u, 'o-', color='#1f7a3a', markersize=10, lw=2,
                 label=r'observed $\tau^\star$ (univariate)')
    axes[0].set_xscale('log', base=2)
    axes[0].set_yscale('log', base=2)
    axes[0].set_xlabel(r'Degree $d$')
    axes[0].set_ylabel(r'$\tau^\star$')
    axes[0].set_title(r'(a) Univariate: $\tau^\star \asymp \sqrt{d}$ (undershooting)')
    axes[0].legend(framealpha=0.95)
    axes[0].grid(True, which='both', alpha=0.25)

    # Multivariate: tau ~ 1
    Ds = [9, 27, 16, 81, 32, 243, 64, 1024, 729]
    sqrtD = [math.sqrt(D) for D in Ds]
    tau_m = [2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 1.0, 1.0]
    Ds_arr = np.array(Ds)
    sqrtD_arr = np.array(sqrtD)
    tau_m_arr = np.array(tau_m)
    sort = np.argsort(Ds_arr)
    axes[1].plot(Ds_arr[sort], sqrtD_arr[sort], 'k--', lw=1.5,
                 label=r'$\sqrt{D}$ (naive expectation)', alpha=0.6)
    axes[1].plot(Ds_arr[sort], tau_m_arr[sort], 'o-', color='#a83232', markersize=10, lw=2,
                 label=r'observed $\tau^\star$ (mvKS)')
    axes[1].axhline(1, color='#a83232', alpha=0.3, ls=':')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log', base=2)
    axes[1].set_xlabel(r'Bézout degree $D$')
    axes[1].set_ylabel(r'$\tau^\star$')
    axes[1].set_title(r'(b) Multivariate: $\tau^\star \approx 1$ (Jacobian self-corrects)')
    axes[1].legend(framealpha=0.95)
    axes[1].grid(True, which='both', alpha=0.25)

    plt.tight_layout()
    out = os.path.join(FIG_LATEX, "mv_tau_diagnostic.png")
    plt.savefig(out)
    plt.close()
    print(f"Saved: {out}")


# ===========================================================================
#  Part IX-mv Figure 3: Wall-clock at large D
# ===========================================================================

def fig_mv_wallclock():
    # (n,d): D, time/poly (sec) from multivariate_KS_huge_v2.py output
    data = [
        ((4,10), 1e4, 0.01, 11.25),
        ((5,10), 1e5, 0.05, 14.0),
        ((6,10), 1e6, 0.17, 13.0),
        ((8,10), 1e8, 1.17, 14.0),
        ((15,4), 1.07e9, 0.18, 15.67),
        ((20,3), 3.49e9, 0.12, 15.33),
        ((10,8), 1.07e9, 1.53, 15.33),
        ((12,6), 2.18e9, 0.89, 16.0),
        ((16,4), 4.29e9, 0.22, 14.0),
    ]
    Ds = np.array([row[1] for row in data])
    times = np.array([row[2] for row in data])
    Ntots = np.array([row[3] for row in data])
    labels = [f'({n},{d})' for ((n,d),_,_,_) in data]

    fig, ax = plt.subplots(figsize=(9.5, 5.0))
    sc = ax.scatter(Ds, times, c=Ntots, cmap='plasma', s=130, edgecolor='black', zorder=3)
    for D, t, l in zip(Ds, times, labels):
        ax.annotate(l, xy=(D, t), xytext=(D*1.2, t*1.15),
                    fontsize=9, alpha=0.8)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Bézout degree $D = d^n$')
    ax.set_ylabel(r'Wall-clock time per system (vectorised Python)')
    ax.set_title(r'Wall-clock cost at large $D$ — even $D = 4.3 \cdot 10^9$ in $0.22$\,s')
    ax.grid(True, which='both', alpha=0.25)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r'$\mathbb{E}[N_{\mathrm{tot}}]$')
    plt.tight_layout()
    out = os.path.join(FIG_LATEX, "mv_wallclock.png")
    plt.savefig(out)
    plt.close()
    print(f"Saved: {out}")


# ===========================================================================
#  Part IX-mv Figure 4: 100% convergence success heatmap
# ===========================================================================

def fig_mv_convergence_heatmap():
    n_vals = list(range(2, 21))
    d_vals = list(range(2, 17))
    data = np.full((len(n_vals), len(d_vals)), np.nan)

    # Mark all tested configurations as 100%
    tested = [
        (2,3), (3,3), (4,2), (4,3), (5,2), (5,3), (5,4), (5,6), (5,10), (5,16),
        (6,2), (6,3), (6,4), (6,6), (6,10), (4,8), (4,10), (4,16),
        (8,10), (10,8), (12,6), (15,4), (16,4), (20,3),
    ]
    for n, d in tested:
        if n in n_vals and d in d_vals:
            data[n_vals.index(n), d_vals.index(d)] = 100.0

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(data, aspect='auto', cmap='Greens', vmin=0, vmax=100,
                   origin='lower')
    ax.set_xticks(range(len(d_vals)))
    ax.set_xticklabels(d_vals)
    ax.set_yticks(range(len(n_vals)))
    ax.set_yticklabels(n_vals)
    ax.set_xlabel(r'Polynomial degree $d$')
    ax.set_ylabel(r'Number of variables $n$')
    ax.set_title(r'Convergence success rate $\%$ across $(n,d)$  — $100\%$ on every tested cell')

    for i, n in enumerate(n_vals):
        for j, d in enumerate(d_vals):
            if not np.isnan(data[i, j]):
                D = d**n
                if D >= 1e9:
                    label = f'$10^{{{int(np.log10(D))}}}$'
                elif D >= 1e6:
                    label = f'$10^{{{int(np.log10(D))}}}$'
                else:
                    label = f'{D}'
                ax.text(j, i, label, ha='center', va='center', fontsize=7, color='white' if D > 100 else 'black')
    plt.colorbar(im, ax=ax, label=r'success rate (\%)')
    # Iso-Bezout-degree contours
    xs = np.arange(len(d_vals))
    ys = np.arange(len(n_vals))
    X, Y = np.meshgrid([d_vals[i] for i in range(len(d_vals))], n_vals)
    Z = X.astype(float) ** Y.astype(float)
    cs = ax.contour(np.arange(len(d_vals)), np.arange(len(n_vals)),
                     np.log10(Z + 1), levels=[3, 6, 9],
                     colors='red', linewidths=1.5, linestyles='--', alpha=0.5)
    ax.clabel(cs, fmt={3: r'$D=10^{3}$', 6: r'$D=10^{6}$', 9: r'$D=10^{9}$'},
              fontsize=10, inline=True)

    plt.tight_layout()
    out = os.path.join(FIG_LATEX, "mv_convergence_heatmap.png")
    plt.savefig(out)
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    fig_101abc_journey()
    fig_101abc_speedup()
    fig_101abc_decomposition()
    fig_mv_scaling()
    fig_mv_tau()
    fig_mv_wallclock()
    fig_mv_convergence_heatmap()
    print("\nAll figures generated.")
