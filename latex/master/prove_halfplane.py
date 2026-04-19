"""
PROOF of the Existential Half-Plane Theorem.

STRATEGY: We prove the UNIVERSAL version for R ≥ 3ρ, which is stronger.

The key inequality: for |w| ≤ ρ/R ≤ 1/3 and θ = π/d,

  |Φ| ≤ Σ_{n≥1} (nπ/d · d · (ρ/R)^n) / n = π · Σ (ρ/R)^n = πρ/R / (1-ρ/R)

For ρ/R ≤ 1/3:  |Φ| ≤ (π/3)/(2/3) = π/2.

The inequality sin(x) < x (strict for x > 0) makes it STRICT: |Φ| < π/2.

This gives Re(r_s) < 0 for ALL s.
"""

import numpy as np

def compute_rs_and_phase(roots, R):
    """Compute r_s and correction phases for all starts."""
    d = len(roots)
    thetas = 2 * np.pi * np.arange(d) / d
    a_s = R * np.exp(1j * thetas)
    z_s = R * np.exp(1j * (thetas + np.pi/d))
    
    results = []
    for s in range(d):
        # Direct computation
        Pa = np.prod(a_s[s] - roots)
        Pz = np.prod(z_s[s] - roots)
        r_s = Pz / Pa
        
        # Correction phase via factorization
        w = roots / (R * np.exp(1j * thetas[s]))
        factors = (1 - w * np.exp(-1j * np.pi/d)) / (1 - w)
        Pi = np.prod(factors)
        Phi = np.angle(Pi)
        
        # Verify: r_s = -Pi
        r_check = -Pi
        
        results.append({
            'r_s': r_s, 'Phi': Phi, 'max_w': np.max(np.abs(w)),
            'r_check': r_check
        })
    return results

# ============================================================
print("=" * 70)
print("THEOREM VERIFICATION: |Φ| < π·(ρ/R)/(1-ρ/R)")
print("=" * 70)

print("\n--- Analytical bound vs actual |Φ| ---")
print(f"{'d':>4} {'ρ/R':>8} {'bound':>10} {'π/2':>8} {'actual':>10} {'safe':>6}")
print("-" * 55)

for d in [3, 5, 7, 10, 20, 50, 100, 200, 500]:
    for rho_ratio in [1/3, 0.3, 0.25, 0.2, 0.49]:
        R = 1/rho_ratio  # so ρ/R = rho_ratio (with ρ=1)
        roots = np.exp(2j * np.pi * np.arange(d) / d) * 0.99
        rho = max(abs(roots))
        R_actual = rho / rho_ratio
        
        # Only test the worst case: all roots at same point
        roots_worst = np.ones(d) * rho
        
        analytical_bound = np.pi * rho_ratio / (1 - rho_ratio)
        
        # Compute actual phase for worst case
        results = compute_rs_and_phase(roots_worst, R_actual)
        max_Phi = max(abs(r['Phi']) for r in results)
        max_w = max(r['max_w'] for r in results)
        
        safe = max_Phi < np.pi/2
        
        if d in [3, 10, 50, 200] or (rho_ratio == 1/3):
            print(f"{d:>4} {rho_ratio:>8.4f} {analytical_bound:>10.4f} {np.pi/2:>8.4f} "
                  f"{max_Phi:>10.4f} {'✓' if safe else '✗'}")

# ============================================================
print("\n" + "=" * 70)
print("KEY VERIFICATION: convergence to π/2 as d → ∞ with ρ/R = 1/3")
print("  (this confirms the bound is TIGHT)")
print("=" * 70)

print(f"\n{'d':>5} {'actual Φ':>12} {'π/2':>8} {'gap':>10} {'gap → 0?':>10}")
print("-" * 50)

for d in [3, 5, 10, 20, 50, 100, 200, 500, 1000]:
    rho = 1.0
    R = 3.0 * rho
    
    # Worst case: all d roots at same point ζ = ρ = 1
    # Then S_n = d · (1/3)^n for all n
    # Φ = d · Σ (1/3)^n sin(nπ/d) / n
    
    Phi_exact = 0
    for n in range(1, 5000):
        Phi_exact += (1/3)**n * np.sin(n * np.pi / d) / n
    Phi_exact *= d
    
    gap = np.pi/2 - Phi_exact
    print(f"{d:>5} {Phi_exact:>12.6f} {np.pi/2:>8.6f} {gap:>10.6f} {'→ 0' if gap < 0.01 else ''}")

# ============================================================
print("\n" + "=" * 70)
print("STRESS TEST: R = 3ρ on ALL adversarial families")
print("=" * 70)

np.random.seed(42)
families = []

for d in [3, 5, 7, 10, 15, 20, 30, 50, 100]:
    rho = 1.0
    R = 3.0 * rho
    
    # 1. All roots at ρ (worst case for the bound)
    families.append((d, np.ones(d) * rho, R, f"AllAt1 d={d}"))
    
    # 2. Roots of unity
    families.append((d, rho * np.exp(2j*np.pi*np.arange(d)/d), R, f"Unity d={d}"))
    
    # 3. Clustered near +1
    roots = 0.99 + 0.01j * np.linspace(-1, 1, d)
    families.append((d, roots, R, f"NearDeg d={d}"))
    
    # 4. Clustered near -1
    roots = -0.99 + 0.01j * np.linspace(-1, 1, d)
    families.append((d, roots, R, f"Near-1 d={d}"))
    
    # 5. Random in unit disk (10 trials)
    for trial in range(10):
        roots = np.random.rand(d) * np.exp(2j*np.pi*np.random.rand(d))
        families.append((d, roots, R, f"Rand d={d}"))

violations = 0
worst_margin = np.inf
for d, roots, R, label in families:
    rho = max(abs(roots))
    if R < 3 * rho:
        R = 3 * rho  # Ensure theorem conditions
    
    results = compute_rs_and_phase(roots, R)
    max_Phi = max(abs(r['Phi']) for r in results)
    max_Re = max(np.real(r['r_s']) for r in results)
    margin = np.pi/2 - max_Phi
    
    if max_Re >= 0:
        violations += 1
        print(f"  VIOLATION: {label} max Re = {max_Re:.6f}")
    
    worst_margin = min(worst_margin, margin)

print(f"\nTotal: {len(families)} families tested")
print(f"Violations: {violations}")
print(f"Worst margin (π/2 - max|Φ|): {worst_margin:.6f} rad = {np.degrees(worst_margin):.3f}°")

# ============================================================
print("\n" + "=" * 70)
print("THE PROOF (LaTeX-ready)")
print("=" * 70)
print(r"""
\begin{theorem}[Universal half-plane containment]
\label{thm:halfplane_proved}
Let $P(z) = \prod_{k=1}^d(z-\zeta_k)$ be monic of degree $d \geq 3$ with
$|\zeta_k| \leq \rho$.  For $R \geq 3\rho$, define equispaced starts
$a_s = Re^{2\pi is/d}$, $z_s = Re^{i(2\pi s/d + \pi/d)}$, $r_s = P(z_s)/P(a_s)$.
Then $\operatorname{Re}(r_s) < 0$ for every $s \in \{0,\ldots,d-1\}$.
\end{theorem}

\begin{proof}
\textbf{Step 1: Factorisation.}
Each ratio factorises as
\[
  r_s = -\prod_{k=1}^d \frac{1 - w_{k,s}\,e^{-i\pi/d}}{1 - w_{k,s}},
  \qquad w_{k,s} = \frac{\zeta_k}{R\,e^{i\theta_s}},
  \quad |w_{k,s}| \leq \frac{\rho}{R} \leq \frac{1}{3}.
\]
It suffices to show the correction phase
$\Phi_s := \arg\bigl(\prod_k(1-w_k e^{-i\pi/d})/(1-w_k)\bigr)$
satisfies $|\Phi_s| < \pi/2$.

\medskip
\textbf{Step 2: Power series.}
Since $|w_{k,s}| \leq 1/3 < 1$, the logarithm converges:
\[
  \log\frac{1 - w e^{-i\theta}}{1-w}
  = -\sum_{n=1}^{\infty} \frac{w^n(e^{-in\theta}-1)}{n},
  \qquad \theta = \pi/d.
\]
Summing over~$k$ and taking the imaginary part:
\[
  \Phi_s = -\operatorname{Im}\sum_{n=1}^{\infty}
  \frac{(e^{-in\pi/d}-1)\,S_n}{n},
  \qquad S_n = \sum_{k=1}^d w_{k,s}^n.
\]

\medskip
\textbf{Step 3: Triangle inequality.}
$|S_n| \leq d\,(\rho/R)^n$ and
$|e^{-in\pi/d}-1| = 2\bigl|\sin(n\pi/(2d))\bigr| < n\pi/d$
(strict, since $\sin x < x$ for $x > 0$).  Therefore:
\begin{align*}
  |\Phi_s|
  &< \sum_{n=1}^{\infty}
     \frac{(n\pi/d)\,d\,(\rho/R)^n}{n}
  = \pi\sum_{n=1}^{\infty}\Bigl(\frac{\rho}{R}\Bigr)^{\!n}
  = \frac{\pi\,\rho/R}{1-\rho/R}.
\end{align*}
For $R \geq 3\rho$: $\rho/R \leq 1/3$, so
$|\Phi_s| < \dfrac{\pi/3}{2/3} = \dfrac{\pi}{2}$.

\medskip
\textbf{Step 4: Conclusion.}
Since $r_s = -e^{i\Phi_s}|\Pi_s|$ with $|\Phi_s| < \pi/2$,
we have $\arg(r_s) = \pi + \Phi_s \in (\pi/2, 3\pi/2)$,
giving $\operatorname{Re}(r_s) < 0$.
\end{proof}

\begin{remark}[Tightness]
The bound is tight: for $P(z) = (z-\rho)^d$ with $R = 3\rho$,
the correction phase $\Phi_s \to \pi/2^-$ as $d \to \infty$.
For $R \geq 4\rho$, the sharper bound
$|\Phi_s| < \pi/3$ holds.
\end{remark}
""")
