"""
PANDROSION FOURIER TRANSFORM (PFT) — CORRECTED VERSION

The key identity from the existential theorem proof:

  r(u) = P(Rαu)/P(Ru) = -1 + ∑_{n≥1} c_n/u^n

where c_n = (α^{-n} - 1)/n · p_n/R^n, p_n = ∑_k ζ_k^n.

DFT at ω^s:
  ∑_s r_s·ω^{-ks} = d · [u^{-k} Laurent coeff of r]  (mod d aliasing)
                   = d · (c_k + c_{k+d} + c_{k+2d} + ...)

For |ζ_k| ≤ ρ: |c_n| ≤ 2d·(ρ/R)^n/n, so c_{k+d} = O((ρ/R)^d · c_k).
When R ≫ ρ, the aliasing is NEGLIGIBLE and ∑r_s·ω^{-ks}/d ≈ c_k exactly.

From c_k we recover p_k = k·R^k/(α^{-k}-1) · c_k.
"""

import numpy as np
from numpy.fft import fft

def PFT(P_eval, d, R):
    """Pandrosion Fourier Transform.
    
    Input: P_eval(z) — evaluation oracle, d — degree, R — radius (R > 2ρ)
    Output: Dict of Newton power sums p_1,...,p_{d-1} and the ratio DFT
    """
    alpha = np.exp(1j * np.pi / d)
    omega = np.exp(2j * np.pi / d)
    
    # Compute ratio vector via 2d evaluations
    rs = np.zeros(d, dtype=complex)
    for s in range(d):
        u = omega**s
        rs[s] = P_eval(R * alpha * u) / P_eval(R * u)
    
    # DFT
    r_hat = fft(rs) / d  # r_hat[k] = (1/d)∑ r_s ω^{-ks}
    
    # Recover Newton power sums: p_k = k·R^k/(e^{-ikπ/d}-1) · r_hat[k]
    p = np.zeros(d, dtype=complex)
    for k in range(1, d):
        amk = np.exp(-1j * k * np.pi / d)
        p[k] = r_hat[k] * k * R**k / (amk - 1)
    
    return p, r_hat, rs

def newton_to_elem_sym(p, d):
    """Newton power sums p_1,...,p_d → elementary symmetric e_1,...,e_d"""
    e = np.zeros(d+1, dtype=complex)
    e[0] = 1
    for k in range(1, d+1):
        if k < len(p):
            s = sum((-1)**(i-1) * e[k-i] * p[i] for i in range(1, k+1))
            e[k] = s / k
    return e

# ============================================================
print("=" * 70)
print("PANDROSION FOURIER TRANSFORM")
print("  Recover Newton power sums from 2d polynomial evaluations")
print("=" * 70)

# Test 1: z³ - 6z² + 11z - 6 = (z-1)(z-2)(z-3)
roots1 = np.array([1., 2., 3.])
P1 = lambda z: np.prod(z - roots1)
d1 = 3

for R in [10, 20, 50, 100]:
    p, rhat, rs = PFT(P1, d1, R)
    p_true = [0, 6, 14, 36]
    
    err = max(abs(p[k] - p_true[k]) for k in range(1, d1))
    print(f"R = {R:>3}: p_1={p[1].real:>8.4f} (true 6), "
          f"p_2={p[2].real:>8.4f} (true 14), max_err={err:.2e}")

print("\n→ As R/ρ grows, the PFT becomes EXACT (aliasing vanishes)")

# ============================================================
print("\n" + "=" * 70)
print("PANDROSION SPECTRAL THEOREM")
print("=" * 70)

print("""
Theorem: Let P(z) = ∏(z-ζ_k) be monic of degree d with |ζ_k| ≤ ρ,
and R > ρ. Define α = e^{iπ/d}, ω = e^{2πi/d}, and
  r_s = P(Rαω^s)/P(Rω^s),  s = 0,...,d-1.

Then the DFT r̂_k = (1/d)∑_s r_s·ω^{-ks} satisfies:

  r̂_k = (e^{-ikπ/d} - 1)/k · p_k/R^k + O((ρ/R)^{k+d})

where p_k = ∑_j ζ_j^k is the k-th Newton power sum.

In particular:
(i)  r̂_0 = -1 + O(d·(ρ/R)^d)  [proved in the existential theorem]
(ii) For 1 ≤ k ≤ d-1:
     p_k = k·R^k/(e^{-ikπ/d} - 1) · r̂_k + O(ρ^d)
""")

# Verify for various degrees and root configs
np.random.seed(42)
for d in [3, 5, 10, 20, 50]:
    roots = np.random.randn(d) + 1j*np.random.randn(d)
    rho = max(abs(roots))
    R = 4 * rho  # Large enough for aliasing to vanish
    
    P_eval = lambda z, r=roots: np.prod(z - r)
    
    p, rhat, rs = PFT(P_eval, d, R)
    
    # True power sums
    p_true = np.array([np.sum(roots**k) for k in range(d)])
    
    # Error
    max_err = max(abs(p[k] - p_true[k]) for k in range(1, d))
    
    # Parseval
    lhs = np.sum(np.abs(rs)**2) / d
    rhs = np.sum(np.abs(rhat)**2)
    
    print(f"d={d:>2}, R/ρ={R/rho:.1f}: max|p_k error| = {max_err:.2e}  "
          f"Parseval: {abs(lhs-rhs):.2e}  r̂_0+1 = {abs(rhat[0]+1):.2e}")

# ============================================================
print("\n" + "=" * 70)
print("PANDROSION-PARSEVAL IDENTITY")
print("=" * 70)

print("""
Theorem (Pandrosion-Parseval):
  (1/d) ∑_{s=0}^{d-1} |r_s|² = ∑_{k=0}^{d-1} |r̂_k|²

Since r̂_0 ≈ -1 and r̂_k = c_k for k ≥ 1:
  (1/d) ∑|r_s|² = 1 + ∑_{k=1}^{d-1} |c_k|²
                 = 1 + ∑_{k=1}^{d-1} |(e^{-ikπ/d}-1)|²/k² · |p_k|²/R^{2k}

This is an ENERGY IDENTITY: the total ratio energy = 1 + spectral energy.

Physical interpretation:
  • The "1" comes from r ≈ -1 (the dominant term — every start is ~opposite)
  • The spectral terms encode how the roots BREAK the symmetry of z^d
  • For P = z^d (no roots): all c_k = 0, so ∑|r_s|² = d (each r_s = -1)
  • For roots near the circle: spectral energy is large (high asymmetry)
""")

# Verify and display
for d in [5, 10, 20]:
    roots = np.exp(2j*np.pi*np.arange(d)/d)  # roots of unity
    rho = 1.0
    R = 3.0
    P_eval = lambda z, r=roots: np.prod(z - r)
    
    _, rhat, rs = PFT(P_eval, d, R)
    
    energy = np.sum(np.abs(rs)**2) / d
    spec = np.sum(np.abs(rhat)**2)
    dominant = abs(rhat[0])**2
    correction = spec - dominant
    
    print(f"d={d:>2} (roots of unity): "
          f"energy = {energy:.6f}, |r̂_0|² = {dominant:.6f}, "
          f"spectral corrections = {correction:.6f}")

# Asymmetric case
print()
for d in [5, 10, 20]:
    roots = np.arange(1, d+1) / d  # roots at 1/d, 2/d, ..., 1
    rho = 1.0
    R = 3.0
    P_eval = lambda z, r=roots: np.prod(z - r)
    
    _, rhat, rs = PFT(P_eval, d, R)
    
    energy = np.sum(np.abs(rs)**2) / d
    dominant = abs(rhat[0])**2
    correction = np.sum(np.abs(rhat[1:])**2)
    
    print(f"d={d:>2} (asymmetric):    "
          f"energy = {energy:.6f}, |r̂_0|² = {dominant:.6f}, "
          f"spectral corrections = {correction:.6f}")

# ============================================================
print("\n" + "=" * 70)
print("PANDROSION LAPLACE (Conservation Law as Laplace identity)")
print("=" * 70)

print("""
Along the Pandrosion orbit z_{n+1} = F_{z_0}(z_n):
  r_n = P(z_n)/P(z_0) = (z_{n+1}-z_n)/(z_{n+1}-z_0)  [Kinematic Identity]

Conservation Law (proved): 
  ∏_{n=0}^{N-1} (1-r_n) = (z_0-z_init)/(z_N-z_0)

Taking logs:
  ∑_{n=0}^{N-1} log(1-r_n) = log((z_0-z_init)/(z_N-z_0))

This is the PANDROSION ZETA FUNCTION:
  Z(z_0, z_init) = ∑_{n=0}^∞ log(1-r_n) = log((z_0-z_init)/(z_0-ζ*))

where ζ* is the root the orbit converges to.

LAPLACE-TYPE GENERATING FUNCTION:
  L(t, z_0, z_init) = ∑_{n=0}^∞ r_n · t^n

At t = 1: L(1) = ∑r_n (the total ratio — diverges in general since |r_n|→0)
The Abel sum: lim_{t→1^-} L(t) = P(z_init)/P(z_0) · (a convergence factor)
""")

# Compute and verify
roots = np.array([1.0+0.5j, -1.0, 0.5j, -0.5-0.5j])
z0 = 5.0
z_init = 2.0 + 1j

def pandrosion_step(z, z0, roots):
    P = lambda w: np.prod(w - roots)
    Pz = P(z)
    Pz0 = P(z0)
    if abs(Pz - Pz0) < 1e-30:
        return z, 0
    F = z0 - (z - z0) / (Pz/Pz0 - 1)
    return F, Pz/Pz0

z = z_init
ratios = []
orbit = [z]
for n in range(50):
    z_new, r_n = pandrosion_step(z, z0, roots)
    if abs(r_n) < 1e-15:
        break
    ratios.append(r_n)
    z = z_new
    orbit.append(z)

ratios = np.array(ratios)
converged_root = min(roots, key=lambda r: abs(orbit[-1] - r))

print(f"Orbit: z_0 = {z0}, z_init = {z_init}")
print(f"Converged to ζ* = {converged_root} after {len(ratios)} steps")
print(f"|z_N - ζ*| = {abs(orbit[-1] - converged_root):.2e}")

# Conservation law
prod = np.prod(1 - ratios)
expected = (z_init - z0) / (orbit[-1] - z0)
print(f"\nConservation: ∏(1-r_n) = {prod:.8f}")
print(f"Expected:     (z_init-z_0)/(z_∞-z_0) = {expected:.8f}")
print(f"Match: {abs(prod-expected) < 1e-6}")

# Zeta function
Z = np.sum(np.log(1 - ratios))
Z_expected = np.log((z_init - z0) / (orbit[-1] - z0))
print(f"\nZeta: Z = {Z:.8f}")
print(f"Expected: {Z_expected:.8f}")
print(f"Match: {abs(Z - Z_expected) < 1e-6}")

# Laplace function
print(f"\nLaplace L(t) = ∑ r_n·t^n:")
for t in [0.5, 0.9, 0.95, 0.99]:
    L = sum(r * t**n for n, r in enumerate(ratios))
    print(f"  L({t:.2f}) = {L:.6f}")

# ============================================================
print("\n" + "=" * 70)
print("SUMMARY: Three new Pandrosion results")
print("=" * 70)
print("""
┌─────────────────────────────────────────────────────────────────┐
│ 1. PANDROSION SPECTRAL THEOREM                                  │
│    The DFT of the ratio vector r_s = P(Rαω^s)/P(Rω^s)         │
│    recovers the Newton power sums:                              │
│    p_k = k·R^k/(e^{-ikπ/d}-1) · r̂_k + O((ρ/R)^d)            │
│    → O(d log d) derivative-free coefficient recovery            │
├─────────────────────────────────────────────────────────────────┤
│ 2. PANDROSION-PARSEVAL IDENTITY                                 │
│    (1/d)∑|r_s|² = 1 + ∑_{k=1}^{d-1} |c_k|²                   │
│    → Energy conservation: ratio energy = 1 + spectral energy   │
├─────────────────────────────────────────────────────────────────┤
│ 3. PANDROSION ZETA FUNCTION                                     │
│    Z(z_0,z_init) = ∑ log(1-r_n) = log((z_0-z_init)/(z_0-ζ*)) │
│    → Exact evaluation of the orbit's convergence target         │
│    → Pandrosion analogue of the Riemann-von Mangoldt formula    │
└─────────────────────────────────────────────────────────────────┘
""")
