"""
Research script for new sections of the Pandrosion paper.
#1: Aitken acceleration ∘ Pandrosion → Newton?
#2: Extension to x < 1
#3: Non-asymptotic bounds
#4: Functional properties of λ_{p,x}
#6: Optimal starting point
"""
import numpy as np
from fractions import Fraction

def S_p(s, p):
    return sum(s**k for k in range(p))

def pandrosion_s(s, p, x):
    return 1 - (x - 1) / (x * S_p(s, p))

def lambda_theoretical(p, x):
    alpha = x**(1/p)
    num = (alpha - 1) * sum(k * alpha**(p-1-k) for k in range(1, p))
    den = sum(alpha**k for k in range(p))
    return num / den

def newton_step(u, p, x):
    """Newton iteration for u^p = x"""
    return ((p-1)*u + x/u**(p-1)) / p

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  #1: AITKEN Δ² ACCELERATION OF PANDROSION                              ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("=" * 78)
print("#1: AITKEN Δ² ACCELERATION OF PANDROSION")
print("=" * 78)

def aitken_acceleration(s0, s1, s2):
    """Aitken Δ² process"""
    denom = s2 - 2*s1 + s0
    if abs(denom) < 1e-15:
        return s2
    return s0 - (s1 - s0)**2 / denom

print("\n--- Testing Aitken on Pandrosion for p=3, x=2 ---")
p, x = 3, 2
alpha = x**(1/p)
s = 0.875
print(f"Target: α = {alpha:.15f}, s* = {1/alpha:.15f}")
print(f"\n{'n':>3} {'s_n':>18} {'v_n=x·s^(p-1)':>18} {'ε_n(Pandrosion)':>18} {'s_aitken':>18} {'v_aitken':>18} {'ε_aitken':>18}")

s_history = [s]
for i in range(15):
    s = pandrosion_s(s, p, x)
    s_history.append(s)

# Apply Aitken Δ² to consecutive triples
print("\nAitken accelerated sequence:")
aitken_eps = []
pandrosion_eps = []
for i in range(len(s_history) - 2):
    s_a = aitken_acceleration(s_history[i], s_history[i+1], s_history[i+2])
    v_a = x * s_a**(p-1)
    v_p = x * s_history[i+2]**(p-1)
    eps_a = abs(v_a - alpha)
    eps_p = abs(v_p - alpha)
    aitken_eps.append(eps_a)
    pandrosion_eps.append(eps_p)
    if i < 10:
        print(f"  n={i}: Pandrosion ε={eps_p:.6e}  |  Aitken ε={eps_a:.6e}  |  speedup={eps_p/max(eps_a,1e-16):.1f}x")

# Check if Aitken convergence is quadratic
print("\nAitken convergence order analysis:")
for i in range(1, min(8, len(aitken_eps)-1)):
    if aitken_eps[i-1] > 1e-14 and aitken_eps[i] > 1e-14:
        # Check ε_{n+1}/ε_n (linear) vs ε_{n+1}/ε_n² (quadratic)
        ratio_lin = aitken_eps[i] / aitken_eps[i-1]
        ratio_quad = aitken_eps[i] / aitken_eps[i-1]**2 if aitken_eps[i-1] > 1e-10 else float('inf')
        print(f"  n={i}: ε_{i+1}/ε_{i} = {ratio_lin:.6f}  |  ε_{i+1}/ε_{i}² = {ratio_quad:.4f}")

# Compare Aitken-Pandrosion output with Newton output
print("\n--- Compare Aitken(Pandrosion) vs Newton step-by-step ---")
print("Starting from same v_0:")
v_pand = x * 0.875**(p-1)  # Pandrosion v_0
u_newt = v_pand  # Newton starting from same point

s = 0.875
s_hist = [s]
for i in range(8):
    s = pandrosion_s(s, p, x)
    s_hist.append(s)

print(f"\n{'step':>5} {'Newton v_n':>18} {'Aitken(Pand) v_n':>18} {'Diff':>15}")
for i in range(min(6, len(s_hist)-2)):
    u_newt = newton_step(u_newt, p, x)
    s_a = aitken_acceleration(s_hist[i], s_hist[i+1], s_hist[i+2])
    v_a = x * s_a**(p-1)
    print(f"  {i+1:>3}   {u_newt:.15f}   {v_a:.15f}   {abs(u_newt - v_a):.2e}")

# Test for p=2
print("\n--- Testing Aitken on Pandrosion for p=2, x=2 ---")
p, x = 2, 2
alpha = x**(1/p)
s = 0.75
s_hist = [s]
for i in range(12):
    s = pandrosion_s(s, p, x)
    s_hist.append(s)

print("Aitken convergence:")
aitken_eps_p2 = []
for i in range(len(s_hist) - 2):
    s_a = aitken_acceleration(s_hist[i], s_hist[i+1], s_hist[i+2])
    v_a = x * s_a**(p-1)
    eps_a = abs(v_a - alpha)
    aitken_eps_p2.append(eps_a)
    if i < 8:
        print(f"  n={i}: ε_aitken = {eps_a:.6e}")

print("\nAitken order for p=2:")
for i in range(1, min(6, len(aitken_eps_p2)-1)):
    if aitken_eps_p2[i-1] > 1e-14 and aitken_eps_p2[i] > 1e-14:
        ratio_quad = aitken_eps_p2[i] / aitken_eps_p2[i-1]**2
        print(f"  n={i}: ε_{i+1}/ε_{i}² = {ratio_quad:.6f}")

# Steffensen's method (iterated Aitken)
print("\n--- Steffensen's method (iterated Aitken) for p=3, x=2 ---")
p, x = 3, 2
alpha = x**(1/p)
s = 0.875
print(f"  Target s* = {1/alpha:.15f}")
for i in range(8):
    s0 = s
    s1 = pandrosion_s(s0, p, x)
    s2 = pandrosion_s(s1, p, x)
    s_new = aitken_acceleration(s0, s1, s2)
    v = x * s_new**(p-1)
    eps = abs(v - alpha)
    print(f"  Steffensen step {i}: s = {s_new:.15f}, v = {v:.15f}, ε = {eps:.2e}")
    if eps < 1e-14:
        print("  → Machine precision reached!")
        break
    s = s_new

# Steffensen convergence order
print("\nSteffensen convergence order check:")
p, x = 3, 2
alpha = x**(1/p)
s = 0.875
steff_eps = []
for i in range(8):
    s0 = s
    s1 = pandrosion_s(s0, p, x)
    s2 = pandrosion_s(s1, p, x)
    s_new = aitken_acceleration(s0, s1, s2)
    v = x * s_new**(p-1)
    eps = abs(v - alpha)
    steff_eps.append(eps)
    if eps < 1e-15:
        break
    s = s_new

for i in range(1, len(steff_eps)):
    if steff_eps[i-1] > 1e-14:
        ratio = steff_eps[i] / steff_eps[i-1]**2
        print(f"  ε_{i}/ε_{i-1}² = {ratio:.6f}  (constant ⇒ quadratic)")

# What is K_steffensen? Compare with K_Newton = (p-1)/(2α)
K_newton = (p-1) / (2*alpha)
print(f"\n  K_Newton = (p-1)/(2α) = {K_newton:.6f}")
if len(steff_eps) >= 2 and steff_eps[0] > 1e-10:
    K_steff = steff_eps[1] / steff_eps[0]**2
    print(f"  K_Steffensen ≈ {K_steff:.6f}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  #2: EXTENSION TO x < 1                                                 ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n\n" + "=" * 78)
print("#2: EXTENSION TO x < 1")
print("=" * 78)

# For x < 1: s* = x^{-1/p} > 1, so s* is outside (0,1)!
# Need to rethink the iteration domain
print("\n--- Fixed point for x < 1 ---")
for p, x in [(2, 0.5), (3, 0.5), (2, 0.25), (3, 0.125)]:
    alpha = x**(1/p)
    s_star = x**(-1/p)
    print(f"  p={p}, x={x}: α={alpha:.6f}, s*={s_star:.6f} {'(in (0,1))' if 0 < s_star < 1 else '(OUTSIDE (0,1)!)'}")

# For x < 1: we need s* = x^{-1/p} = α^{-1} > 1
# The iteration h(s) = 1 - (x-1)/(x·S_p(s)) with x<1 means x-1 < 0
# So h(s) = 1 + (1-x)/(x·S_p(s)) > 1 for s > 0
# The iteration goes to s > 1

print("\n--- Testing iteration directly for x < 1 ---")
for p, x in [(2, 0.5), (3, 0.5), (3, 0.125)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    print(f"\n  p={p}, x={x}: α={alpha:.6f}, s*=1/α={s_star:.6f}")
    s = 1.5  # start > 1
    for i in range(30):
        s = pandrosion_s(s, p, x)
    v = x * s**(p-1)
    print(f"  After 30 iter: s={s:.10f}, v={v:.10f}, target α={alpha:.10f}")
    print(f"  Converged: {abs(v - alpha) < 1e-8}")
    
    # Check λ formula still works
    lam = lambda_theoretical(p, x)
    print(f"  λ_{p},{x} = {lam:.6f}")
    if 0 < lam < 1:
        print(f"  λ ∈ (0,1) ✓")
    else:
        print(f"  λ NOT in (0,1)!")

# Numerical convergence rate for x < 1
print("\n--- Convergence rate verification for x < 1 ---")
for p, x in [(2, 0.5), (3, 0.5), (2, 0.25), (3, 0.125)]:
    alpha = x**(1/p)
    s = 1.5
    lam_th = lambda_theoretical(p, x)
    ratios = []
    for i in range(60):
        v = x * s**(p-1)
        eps = abs(v - alpha)
        s_new = pandrosion_s(s, p, x)
        v_new = x * s_new**(p-1)
        eps_new = abs(v_new - alpha)
        if eps > 1e-13 and eps_new > 1e-13:
            ratios.append(eps_new / eps)
        s = s_new
    if len(ratios) >= 5:
        ratio_final = np.mean(ratios[-5:])
        match = abs(ratio_final - abs(lam_th)) < 0.01
        print(f"  p={p}, x={x}: λ_th={lam_th:.6f}, ratio_final={ratio_final:.6f} {'✓' if match else '✗'}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  #3: NON-ASYMPTOTIC BOUNDS                                              ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n\n" + "=" * 78)
print("#3: NON-ASYMPTOTIC BOUNDS")
print("=" * 78)

# For a contraction mapping h on [a,b] with |h'(s)| ≤ λ for all s ∈ [a,b],
# we have |s_n - s*| ≤ λ^n |s_0 - s*|
# Need to show h' is bounded by some Λ on the entire basin

print("\n--- Check if h' is monotone and bounded on (0,1) for x > 1 ---")
for p, x in [(3, 2), (2, 2), (4, 2)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    lam = lambda_theoretical(p, x)
    
    # Sample h'(s) on (0,1)
    ss = np.linspace(0.01, 0.99, 1000)
    h_primes = []
    for s in ss:
        eps = 1e-8
        hp = (pandrosion_s(s+eps, p, x) - pandrosion_s(s-eps, p, x)) / (2*eps)
        h_primes.append(hp)
    
    max_hp = max(h_primes)
    min_hp = min(h_primes)
    hp_at_star = h_primes[np.argmin(np.abs(ss - s_star))]
    
    print(f"  p={p}, x={x}: h'(s*) = {lam:.6f}")
    print(f"    max h'(s) on (0,1) = {max_hp:.6f} (at s={ss[np.argmax(h_primes)]:.4f})")
    print(f"    min h'(s) on (0,1) = {min_hp:.6f}")
    print(f"    h' monotone increasing: {all(h_primes[i] <= h_primes[i+1]+1e-6 for i in range(len(h_primes)-1))}")

# For non-asymptotic bound: need sup |h'(s)| on [s_0, s*] or [s*, s_0]
print("\n--- Non-asymptotic bound: Λ = sup h' on basin ---")
for p, x in [(3, 2), (2, 2)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    s0 = 0.5  # typical starting point
    
    # Basin: [s0, s*] since s0 < s*
    ss = np.linspace(s0, 0.99, 500)
    h_primes = []
    for s in ss:
        eps = 1e-8
        hp = (pandrosion_s(s+eps, p, x) - pandrosion_s(s-eps, p, x)) / (2*eps)
        h_primes.append(hp)
    
    Lambda = max(h_primes)
    lam = lambda_theoretical(p, x)
    print(f"  p={p}, x={x}: λ(s*) = {lam:.6f}, Λ = sup h' on [{s0}, 1) = {Lambda:.6f}")
    print(f"    Non-asymptotic: |s_n - s*| ≤ Λ^n · |s_0 - s*|")
    
    # Verify numerically
    s = s0
    for n in range(10):
        err = abs(s - s_star)
        bound = Lambda**n * abs(s0 - s_star)
        ok = err <= bound * 1.001
        if n < 6:
            print(f"    n={n}: err={err:.6e}, bound={bound:.6e} {'✓' if ok else '✗'}")
        s = pandrosion_s(s, p, x)

# Check if h' is maximized at s=1 (boundary)
print("\n--- h'(s) behavior near s=1 ---")
for p, x in [(3, 2), (2, 2), (4, 2), (3, 8)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    # h'(1^-) analytically
    # h(s) = 1 - (x-1)/(x·S_p(s))
    # S_p(1) = p
    # S'_p(1) = sum(k for k=1..p-1) = p(p-1)/2
    # h'(1) = (x-1)/x · S'_p(1)/S_p(1)^2 = (x-1)/x · p(p-1)/2 / p^2 = (x-1)(p-1)/(2xp)
    hp_at_1 = (x-1)*(p-1) / (2*x*p)
    lam = lambda_theoretical(p, x)
    print(f"  p={p}, x={x}: h'(1) = {hp_at_1:.6f}, h'(s*) = λ = {lam:.6f}, h'(1) < λ: {hp_at_1 < lam}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  #4: FUNCTIONAL PROPERTIES OF λ_{p,x}                                   ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n\n" + "=" * 78)
print("#4: FUNCTIONAL PROPERTIES OF λ_{p,x}")
print("=" * 78)

# Monotonicity in x (fixed p)
print("\n--- Monotonicity: λ_{p,x} increasing in x for fixed p ---")
for p in [2, 3, 4, 5]:
    xs = [1.01, 1.1, 1.5, 2, 3, 5, 10, 50, 100]
    lams = [lambda_theoretical(p, x) for x in xs]
    increasing = all(lams[i] < lams[i+1] for i in range(len(lams)-1))
    print(f"  p={p}: λ strictly increasing in x: {increasing}")
    print(f"    values: {[f'{l:.4f}' for l in lams]}")

# Monotonicity in p (fixed x)
print("\n--- Monotonicity: λ_{p,x} in p for fixed x ---")
for x in [2, 3, 10]:
    ps = list(range(2, 10))
    lams = [lambda_theoretical(p, x) for p in ps]
    increasing = all(lams[i] < lams[i+1] for i in range(len(lams)-1))
    print(f"  x={x}: λ values for p=2..9: {[f'{l:.4f}' for l in lams]}")
    print(f"    Increasing in p: {increasing}")

# Limit as x → 1+
print("\n--- Limit λ_{p,x} as x → 1+ ---")
for p in [2, 3, 4]:
    for x in [1.001, 1.0001, 1.00001]:
        lam = lambda_theoretical(p, x)
        alpha = x**(1/p)
        # Expected: λ ≈ (p-1)/2 · (α-1) for α near 1?
        approx = (p-1) * (alpha - 1) / 2  # ??? Let's check
        print(f"  p={p}, x={x}: λ={lam:.8f}, (p-1)(α-1)/2 = {approx:.8f}, ratio = {lam/approx:.6f}")

# Better: λ_{p,x} ~ (p-1)(α-1)/(p) as α → 1
print("\n--- Better asymptotic: λ ~ (p-1)(α-1)/... as α → 1 ---")
for p in [2, 3, 4, 5]:
    for x in [1.0001]:
        lam = lambda_theoretical(p, x)
        alpha = x**(1/p)
        # At α=1: num = (α-1)·Σk·1^{p-1-k} = (α-1)·Σk = (α-1)·p(p-1)/2
        # den = Σ1^k = p
        # So λ ~ (α-1)·p(p-1)/(2p) = (α-1)(p-1)/2
        approx = (alpha - 1) * (p-1) / 2
        print(f"  p={p}: λ={lam:.10f}, (α-1)(p-1)/2 = {approx:.10f}, ratio = {lam/approx:.6f}")

# Limit as x → ∞ (α → ∞) for fixed p
print("\n--- Limit λ_{p,x} as α → ∞ ---")
for p in [2, 3, 4]:
    for x in [10**6, 10**9, 10**12]:
        lam = lambda_theoretical(p, x)
        print(f"  p={p}, x={x:.0e}: λ = {lam:.10f}, 1-λ = {1-lam:.6e}")
    # Expect λ → 1, and 1-λ → ?
    # At large α: num ≈ α · (1·α^{p-2}) = α^{p-1} (dominant term)
    # den ≈ α^{p-1}
    # So num/den ≈ (α-1)·1·α^{p-2}/α^{p-1} = (α-1)/α → 1
    # More precisely: 1-λ ~ ?

# Limit as p → ∞ for fixed x
print("\n--- Limit λ_{p,x} as p → ∞ for fixed x ---")
for x in [2, 3, 10]:
    ps = [2, 5, 10, 20, 50, 100, 200]
    lams = [lambda_theoretical(p, x) for p in ps]
    print(f"  x={x}: p → λ:")
    for p, l in zip(ps, lams):
        alpha = x**(1/p)
        print(f"    p={p:>3}: λ={l:.6f}, α={alpha:.6f}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  #6: OPTIMAL STARTING POINT                                             ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n\n" + "=" * 78)
print("#6: OPTIMAL STARTING POINT")
print("=" * 78)

# The number of iterations is n ~ ln(ε_0/δ) / ln(1/λ)
# where ε_0 = |v_0 - α| = |x·s_0^{p-1} - α|
# To minimize n, minimize ε_0 = |x·s_0^{p-1} - α|
# This is minimized when s_0 = (α/x)^{1/(p-1)} = (x^{1/p}/x)^{1/(p-1)} = x^{(1/p-1)/(p-1)}
# = x^{(1-p)/(p(p-1))} = x^{-1/p} = s*
# But that's the fixed point itself! So any s_0 ≠ s* has ε_0 > 0.
# The question is: for a "natural" starting point, what is best?

# Actually the question is more nuanced: the convergence is not exactly
# ε_{n+1} = λ·ε_n from the start. The non-asymptotic bound uses Λ = sup h'.
# So the effective bound is Λ^n·ε_0 where Λ depends on [s_0, s*].
# Optimal s_0 minimizes Λ^n·ε_0 subject to Λ = sup_{[s_0,s*]} |h'|.

# Simpler approach: for what s_0 is the convergence fastest?
# Since h' is increasing on (0,1), starting closer to s* gives smaller Λ.
# So the optimal starting point is as close to s* as possible.
# But we don't know s* (that's what we're computing!)

# Natural starting points:
print("\n--- Natural starting points comparison ---")
for p, x in [(3, 2), (2, 2), (4, 2), (3, 8)]:
    alpha = x**(1/p)
    s_star = 1 / alpha
    
    # Option 1: s_0 = 1 - 1/x (from h(0) = 1/x, so first iterate)
    s0_opt1 = 1/x
    # Option 2: s_0 = 1 - (x-1)/(x·p) = h(1^-) asymptotic
    s0_opt2 = 1 - (x-1)/(x*p)
    # Option 3: s_0 = (p-1)/p (heuristic for large p)
    s0_opt3 = (p-1)/p
    # Option 4: s_0 = 1/x^{1/p} ≈ 1 - ln(x)/p for x near 1 (linear approx)
    s0_opt4 = 1 - np.log(x)/p  # first-order Taylor approx of x^{-1/p}
    # Option 5: s_0 = (1 + 1/x)/2 (average of boundary values h(0) and fixed point)
    s0_opt5 = 1 - (x-1)/(x*p)

    starts = {
        "s₀ = 1/2": 0.5,
        "s₀ = 1/x": 1/x,
        "s₀ = 1-(x-1)/(xp)": s0_opt2,
        "s₀ = 1-ln(x)/p": s0_opt4,
        "s₀ = (p-1)/p": s0_opt3,
    }
    
    print(f"\n  p={p}, x={x} (s*={s_star:.6f}):")
    for name, s0 in starts.items():
        if s0 <= 0 or s0 >= 2:
            continue
        eps0 = abs(x * s0**(p-1) - alpha)
        # Count iterations to reach ε < 1e-10
        s = s0
        n_iter = 0
        for _ in range(200):
            v = x * s**(p-1)
            if abs(v - alpha) < 1e-10:
                break
            s = pandrosion_s(s, p, x)
            n_iter += 1
        print(f"    {name:>25}: |s₀-s*|={abs(s0-s_star):.4f}, ε₀={eps0:.4e}, iters={n_iter}")

# Closed-form optimal s_0
# Since we want v_0 closest to α, and v_0 = x·s_0^{p-1},
# the ideal would be s_0 = (α/x)^{1/(p-1)} = s*. 
# But without knowing α, a natural approximation is:
# Use s_0 = 1 - (x-1)/(xp) which is h(1), the image of s=1
# This is actually the "best naive guess" since it uses the iteration itself
print("\n--- Theoretical optimal: s_0 = h(1) = 1 - (x-1)/(xp) ---")
for p, x in [(3, 2), (2, 2), (3, 8), (4, 16)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    s_h1 = 1 - (x-1)/(x*p)
    print(f"  p={p}, x={x}: h(1)={s_h1:.6f}, s*={s_star:.6f}, gap={abs(s_h1-s_star):.4f}")


print("\n\n=== RESEARCH COMPLETE ===")
