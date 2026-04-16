"""
Research: Can we get order 3 or 4 convergence from Pandrosion?

Kung-Traub bounds:
  2 evals → max order 2 (Steffensen achieves this) ✅
  3 evals → max order 4 (can we achieve this?)
  4 evals → max order 8

Approach: Apply higher-order Aitken/Neville extrapolation to the Pandrosion
fixed-point iteration h(s).

Methods to test:
  1. Double Steffensen (Steffensen on Steffensen)
  2. 3-point Neville-Aitken extrapolation
  3. Traub-type composition: h³ with Richardson extrapolation
  4. Householder-like derivative-free methods
"""

import numpy as np

def S_p(s, p):
    if abs(s - 1) < 1e-15: return float(p)
    return (1 - s**p) / (1 - s)

def h(s, p, x):
    """Pandrosion iteration."""
    sp = S_p(s, p)
    if abs(sp) < 1e-30: return s
    return 1 - (x - 1) / (x * sp)

def steffensen(s, p, x):
    """Order 2: Aitken Δ² on h. Uses 2 evals of h."""
    s0 = s
    s1 = h(s0, p, x)
    s2 = h(s1, p, x)
    d = s2 - 2*s1 + s0
    if abs(d) < 1e-30: return s2
    return s0 - (s1 - s0)**2 / d


# ═══════════════════════════════════════════════════════════════════
# METHOD A: Iterated Aitken (3 evals of h → order ??)
# ═══════════════════════════════════════════════════════════════════

def iterated_aitken(s, p, x):
    """Apply h three times, then use 4-point Aitken-Neville."""
    s0 = s
    s1 = h(s0, p, x)  # eval 1
    s2 = h(s1, p, x)  # eval 2
    s3 = h(s2, p, x)  # eval 3
    
    # Level 1: two Aitken Δ² extrapolations
    d1 = s2 - 2*s1 + s0
    if abs(d1) < 1e-30: return s3
    t01 = s0 - (s1 - s0)**2 / d1
    
    d2 = s3 - 2*s2 + s1
    if abs(d2) < 1e-30: return s3
    t12 = s1 - (s2 - s1)**2 / d2
    
    # Level 2: Aitken again on t01, t12
    # But we need a third t value... use direct approach instead
    # Actually: if t01 and t12 are both O(ε²) approximations,
    # another Aitken won't help without knowing their "linear rate"
    
    # Alternative: Richardson-like combination
    # If s_n - s* ~ A·λ^n, then Aitken gives t ~ s* + B·λ^{2n}
    # So t01 - s* ~ C·λ^0 = C (order 2 from s0)
    # And t12 - s* ~ C·λ^2 (shifted)
    # Ratio: (t12 - s*)/(t01 - s*) ~ λ² → can extrapolate!
    
    if abs(t12 - t01) < 1e-30: return t12
    # This is "iterated Aitken": 
    return t01 - (t12 - t01)**2 / (t12 - 2*t12 + t01) if abs(t12 - 2*t12 + t01) > 1e-30 else t12


# ═══════════════════════════════════════════════════════════════════
# METHOD B: Steffensen + one Newton-like correction (3 evals)
# ═══════════════════════════════════════════════════════════════════

def steffensen_plus(s, p, x):
    """Steffensen (2 evals) + one more h eval for higher order."""
    # Step 1: Steffensen gives s_new with error O(ε²)
    s0 = s
    s1 = h(s0, p, x)    # eval 1
    s2 = h(s1, p, x)    # eval 2
    d = s2 - 2*s1 + s0
    if abs(d) < 1e-30: return s2
    s_steff = s0 - (s1 - s0)**2 / d
    
    # Step 2: One more h evaluation at the Steffensen result
    s3 = h(s_steff, p, x)  # eval 3
    
    # Now s_steff has error ε², and h(s_steff) has error λ·ε²
    # Aitken-like correction: s* ≈ s_steff - (s3 - s_steff)/(1 - λ_approx)
    # where λ_approx = (s3 - s_steff) / (s_steff - s0)... but this is tricky
    
    # Simpler: use the "secant correction"
    # s_steff - s* ~ K·(s - s*)²
    # s3 - s* ~ λ·K·(s - s*)²
    # So s3 - s_steff ~ (λ-1)·K·ε²
    # And s* ~ s_steff - (s3 - s_steff)/(λ-1)
    # But we don't know λ... approximate from the Steffensen data
    
    # λ ≈ (s2-s1)/(s1-s0) (classic Aitken estimate)
    lam_est = (s2 - s1) / (s1 - s0) if abs(s1 - s0) > 1e-30 else 0
    
    if abs(lam_est - 1) < 1e-10:
        return s3
    
    return s_steff - (s3 - s_steff) / (lam_est - 1)


# ═══════════════════════════════════════════════════════════════════
# METHOD C: Traub-Steffensen (optimal 3-eval, order 4)
# ═══════════════════════════════════════════════════════════════════

def traub_steffensen(s, p, x):
    """
    Traub's optimal 3-eval method adapted to fixed-point form.
    Uses h(s), h(h(s)), h(h(h(s))) to build order 4.
    
    Key idea: define g(s) = s - (h(s)-s)/(h'_approx(s)-1)
    then apply Steffensen to g.
    """
    s0 = s
    s1 = h(s0, p, x)     # eval 1
    s2 = h(s1, p, x)     # eval 2
    
    # Steffensen correction (order 2)
    d = s2 - 2*s1 + s0
    if abs(d) < 1e-30: return s2
    s_steff = s0 - (s1 - s0)**2 / d
    
    # Now evaluate h at the Steffensen point
    s3 = h(s_steff, p, x)  # eval 3
    
    # Use (s_steff, s3) as a new pair for extrapolation
    # s_steff - s* ~ K·ε²
    # s3 = h(s_steff) ~ s* + λ·(s_steff - s*) + h''(s*)/2·(s_steff-s*)² + ...
    # So s3 - s* ~ λ·K·ε² + h''/2·K²·ε⁴
    # The ratio (s3 - s*)/(s_steff - s*) ~ λ (same λ as before!)
    
    # So we can do another Aitken on {s_steff, s3, h(s3)} but we only have 3 evals
    # Instead: use the linear extrapolation
    # s* = s_steff - (s3 - s_steff) · s_steff_err / (s3_err - s_steff_err)
    # where errors have ratio λ... but λ is known approximately
    
    lam_est = (s1 - s0)
    if abs(s1 - s0) > 1e-30:
        lam_est = (s2 - s1) / (s1 - s0)
    
    # Newton-like correction using estimated slope
    if abs(1 - lam_est) > 1e-10:
        return s_steff - (s3 - s_steff) / (lam_est - 1)
    return s3


# ═══════════════════════════════════════════════════════════════════
# METHOD D: Direct Householder-Pandrosion (uses S_p structure)
# ═══════════════════════════════════════════════════════════════════

def householder_pandrosion(s, p, x):
    """
    Householder-like correction using 3 evals of S_p.
    Key: use the Pandrosion structure directly.
    
    Define F(s) = h(s) - s = -(x-1)/(x·S_p(s)) - s + 1
    F'(s) ≈ [F(s+δ) - F(s-δ)]/(2δ) using h evals
    F''(s) ≈ [F(s+δ) - 2F(s) + F(s-δ)]/δ²
    
    Halley: s - F/(F' - F·F''/(2F'))
    """
    F0 = h(s, p, x) - s       # eval 1
    delta = max(abs(F0), 1e-4) * 0.01  # adaptive step
    
    Fp = h(s + delta, p, x) - (s + delta)  # eval 2
    Fm = h(s - delta, p, x) - (s - delta)  # eval 3
    
    Fprime = (Fp - Fm) / (2 * delta)
    Fdoubleprime = (Fp - 2*F0 + Fm) / (delta**2)
    
    if abs(Fprime) < 1e-30: return s + F0
    
    # Halley formula: cubic convergence
    denom = Fprime - F0 * Fdoubleprime / (2 * Fprime)
    if abs(denom) < 1e-30: return s - F0 / Fprime  # fallback to Newton-like
    
    return s - F0 / denom


# ═══════════════════════════════════════════════════════════════════
# TESTING
# ═══════════════════════════════════════════════════════════════════

p, x = 3, 2.0
alpha = x**(1/p)
s_star = 1/alpha

print("=" * 70)
print(f"Testing higher-order methods for p={p}, x={x} (cube root of 2)")
print(f"α = {alpha:.15f}, s* = {s_star:.15f}")
print("=" * 70)

methods = [
    ("Steffensen (order 2, 2 evals)", steffensen, 2),
    ("Iterated Aitken (3 evals)", iterated_aitken, 3),
    ("Steffensen+ (3 evals)", steffensen_plus, 3),
    ("Traub-Steffensen (3 evals)", traub_steffensen, 3),
    ("Householder-Pandrosion (3 evals)", householder_pandrosion, 3),
]

for name, method, n_evals in methods:
    print(f"\n  {name}:")
    s = 0.5
    eps_list = []
    for n in range(8):
        eps = abs(s - s_star)
        eps_list.append(eps)
        if eps < 1e-16:
            break
        try:
            s = method(s, p, x)
            if abs(s) > 1e10:
                print(f"    DIVERGED at step {n}")
                break
        except Exception as e:
            print(f"    ERROR at step {n}: {e}")
            break
    
    # Print trace
    for i, e in enumerate(eps_list):
        order_str = ""
        if i > 0 and eps_list[i-1] > 1e-14:
            ratio = np.log(e) / np.log(eps_list[i-1]) if eps_list[i-1] > 1e-15 and e > 1e-300 else 0
            order_str = f"  (apparent order ≈ {ratio:.2f})"
        print(f"    Step {i}: |s-s*| = {e:.2e}{order_str}")


print("\n" + "=" * 70)
print("CONVERGENCE ORDER ANALYSIS")
print("=" * 70)

for name, method, n_evals in methods:
    s = 0.8  # closer starting point for cleaner order measurement
    eps_list = []
    for n in range(10):
        eps = abs(s - s_star)
        eps_list.append(eps)
        if eps < 1e-16: break
        try:
            s = method(s, p, x)
            if abs(s) > 1e10: break
        except: break
    
    # Measure order: log(ε_{n+1})/log(ε_n)
    orders = []
    for i in range(1, len(eps_list)):
        if eps_list[i-1] > 1e-14 and eps_list[i] > 1e-300 and eps_list[i-1] < 0.5:
            orders.append(np.log(eps_list[i]) / np.log(eps_list[i-1]))
    
    if orders:
        avg_order = np.mean(orders[-2:]) if len(orders) >= 2 else orders[0]
        eff = avg_order**(1/n_evals)
        steps = len([e for e in eps_list if e > 1e-14])
        print(f"  {name}")
        print(f"    Measured order: {avg_order:.2f} | Evals/step: {n_evals} | E = {eff:.3f} | Steps to 10^-14: {steps}")


print("\n" + "=" * 70)
print("BEST METHOD?")
print("=" * 70)
print("""
Kung-Traub bounds:
  2 evals → order 2 max  (E = 1.414)
  3 evals → order 4 max  (E = 1.587)
  
If we achieve order 4 with 3 Pandrosion evals: E = 4^{1/3} = 1.587 > 1.414
This would BEAT Steffensen-Pandrosion in efficiency index!
""")
