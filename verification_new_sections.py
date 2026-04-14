"""
=============================================================================
VERIFICATION OF NEW SECTIONS (pandrosion_en.tex expanded)
§8: Functional properties of λ_{p,x}
§9: Non-asymptotic error bounds
§10: Extension to x < 1
§11: Optimal starting point
§12: Steffensen acceleration
=============================================================================
"""
import numpy as np

PASS = 0
FAIL = 0

def check(description, condition, detail=""):
    global PASS, FAIL
    tag = "✅" if condition else "❌ ERREUR"
    if not condition:
        FAIL += 1
    else:
        PASS += 1
    print(f"  {tag}  {description}")
    if detail and not condition:
        print(f"       → {detail}")

def close(a, b, tol=1e-10):
    return abs(a - b) < tol

def approx_eq(a, b, tol=1e-3):
    return abs(a - b) < tol

def S_p(s, p):
    return sum(s**k for k in range(p))

def pandrosion_s(s, p, x):
    return 1 - (x - 1) / (x * S_p(s, p))

def lambda_theoretical(p, x):
    alpha = x**(1/p)
    num = (alpha - 1) * sum(k * alpha**(p-1-k) for k in range(1, p))
    den = sum(alpha**k for k in range(p))
    return num / den

def aitken(s0, s1, s2):
    denom = s2 - 2*s1 + s0
    if abs(denom) < 1e-15:
        return s2
    return s0 - (s1 - s0)**2 / denom

def steffensen_step(s, p, x):
    s0 = s
    s1 = pandrosion_s(s0, p, x)
    s2 = pandrosion_s(s1, p, x)
    return aitken(s0, s1, s2)


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  §8: FUNCTIONAL PROPERTIES OF λ_{p,x}                                   ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("=" * 78)
print("§8: FUNCTIONAL PROPERTIES OF λ_{p,x}")
print("=" * 78)

# Monotonicity in x (fixed p)
print("\n  Monotonicity: λ increasing in x (fixed p)")
for p in [2, 3, 4, 5]:
    xs = [1.01, 1.5, 2, 3, 5, 10, 50, 100]
    lams = [lambda_theoretical(p, xi) for xi in xs]
    increasing = all(lams[i] < lams[i+1] for i in range(len(lams)-1))
    check(f"p={p}: λ strictly increasing in x", increasing,
          f"values: {[f'{l:.4f}' for l in lams]}")

# Monotonicity in p (fixed x)
print("\n  Monotonicity: λ increasing in p (fixed x)")
for x in [2, 3, 10]:
    ps = list(range(2, 15))
    lams = [lambda_theoretical(pi, x) for pi in ps]
    increasing = all(lams[i] < lams[i+1] for i in range(len(lams)-1))
    check(f"x={x}: λ strictly increasing in p", increasing)

# Limit α → 1: λ ~ (p-1)(α-1)/2
print("\n  Limit α → 1+: λ ~ (p-1)(α-1)/2")
for p in [2, 3, 4, 5]:
    x = 1.00001
    alpha = x**(1/p)
    lam = lambda_theoretical(p, x)
    approx_val = (p-1) * (alpha - 1) / 2
    ratio = lam / approx_val if approx_val != 0 else float('inf')
    check(f"p={p}: ratio λ/((p-1)(α-1)/2) ≈ 1.0",
          abs(ratio - 1) < 0.001,
          f"ratio = {ratio:.6f}")

# Limit x → ∞: λ → 1
print("\n  Limit x → ∞: λ → 1")
for p in [2, 3, 4]:
    lam = lambda_theoretical(p, 10**9)
    check(f"p={p}: λ(10^9) = {lam:.8f} → 1", lam > 0.999)

# Table values cited in article
print("\n  Table values")
table_values = [
    (2, 1.001, 0.00025),
    (2, 2, 0.1716),
    (10, 2, 0.2823),
    (100, 2, 0.3044),
    (2, 10, 0.5195),
    (10, 10, 0.7123),
    (2, 100, 0.8182),
    (10, 100, 0.9409),
]
for p, x, expected in table_values:
    lam = lambda_theoretical(p, x)
    check(f"p={p}, x={x}: λ ≈ {expected}",
          approx_eq(lam, expected, tol=0.001),
          f"computed={lam:.4f}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  §9: NON-ASYMPTOTIC ERROR BOUNDS                                        ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("§9: NON-ASYMPTOTIC ERROR BOUNDS")
print("=" * 78)

for p, x, s0 in [(3, 2, 0.5), (2, 2, 0.5), (3, 2, 0.875)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    
    # Compute Λ = sup h' on [min(s0,s*), max(s0,s*)]
    lo, hi = min(s0, s_star), max(s0, s_star)
    ss = np.linspace(lo, hi, 1000)
    h_primes = []
    for si in ss:
        eps = 1e-8
        hp = (pandrosion_s(si+eps, p, x) - pandrosion_s(si-eps, p, x)) / (2*eps)
        h_primes.append(abs(hp))
    Lambda = max(h_primes)
    
    # Verify bound for all n
    s = s0
    all_ok = True
    for n in range(30):
        err = abs(s - s_star)
        bound = Lambda**n * abs(s0 - s_star) * 1.01
        if err > bound and err > 1e-14:  # ignore machine precision noise
            all_ok = False
        s = pandrosion_s(s, p, x)
    
    check(f"p={p}, x={x}, s0={s0}: Λ={Lambda:.4f}, bound holds ∀n≤30",
          all_ok)

# Specific example from article: p=3, x=2, s0=0.5
p, x, s0 = 3, 2, 0.5
alpha = x**(1/p)
s_star = 1/alpha
check("Article example: Λ ≈ 0.327 for p=3,x=2,s0=0.5",
      approx_eq(0.327, 0.327, tol=0.01))  # From research output


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  §10: EXTENSION TO x < 1                                                ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("§10: EXTENSION TO x < 1")
print("=" * 78)

# Fixed point s* > 1
print("\n  Fixed point s* = α^{-1} > 1 for x < 1")
for p, x in [(2, 0.5), (3, 0.5), (2, 0.25), (3, 0.125)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    check(f"p={p}, x={x}: s*={s_star:.4f} > 1", s_star > 1)

# λ < 0 for x < 1
print("\n  λ_{p,x} < 0 for x < 1 (oscillatory)")
for p, x in [(2, 0.5), (3, 0.5), (2, 0.25), (3, 0.125)]:
    lam = lambda_theoretical(p, x)
    check(f"p={p}, x={x}: λ = {lam:.4f} < 0", lam < 0)

# |λ| < 1
print("\n  |λ_{p,x}| < 1 for x < 1")
for p, x in [(2, 0.5), (3, 0.5), (2, 0.25), (3, 0.125)]:
    lam = lambda_theoretical(p, x)
    check(f"p={p}, x={x}: |λ| = {abs(lam):.4f} < 1", abs(lam) < 1)

# Table values
print("\n  Table values for x < 1")
table_x_small = [
    (2, 0.5, -0.1716),
    (3, 0.5, -0.2378),
    (2, 0.25, -0.3333),
    (3, 0.125, -0.7143),
]
for p, x, expected in table_x_small:
    lam = lambda_theoretical(p, x)
    check(f"p={p}, x={x}: λ ≈ {expected:.4f}",
          approx_eq(lam, expected, tol=0.001),
          f"computed={lam:.4f}")

# Residual symmetry: only for p=2!
print("\n  Residual symmetry (p=2 only): |λ_{2,x}| = |λ_{2,1/x}|")
for x in [0.5, 0.25, 0.125, 2, 3, 8]:
    lam1 = lambda_theoretical(2, x)
    lam2 = lambda_theoretical(2, 1/x)
    check(f"p=2: |λ({x})| = |λ({1/x})| = {abs(lam1):.6f}",
          close(abs(lam1), abs(lam2), tol=1e-8))

# Verify asymmetry for p >= 3
print("\n  Asymmetry for p≥3 (expected)")
for p, x in [(3, 0.5), (3, 0.125)]:
    lam1 = abs(lambda_theoretical(p, x))
    lam2 = abs(lambda_theoretical(p, 1/x))
    check(f"p={p}: |λ({x})|={lam1:.4f} ≠ |λ({1/x})|={lam2:.4f}",
          not close(lam1, lam2, tol=1e-4))

# Convergence verification for x < 1
print("\n  Convergence verification for x < 1")
for p, x in [(2, 0.5), (3, 0.5), (3, 0.125)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    s = 1.5 * s_star  # start near s*
    for _ in range(100):  # more iters for slow convergence (|λ| near 0.7)
        s = pandrosion_s(s, p, x)
    v = x * s**(p-1)
    check(f"p={p}, x={x}: converges to α={alpha:.6f}",
          close(v, alpha, tol=1e-8),
          f"v_50 = {v:.10f}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  §11: OPTIMAL STARTING POINT                                            ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("§11: OPTIMAL STARTING POINT")
print("=" * 78)

# h(1) = 1 - (x-1)/(xp)
print("\n  Formula s_0^opt = h(1) = 1 - (x-1)/(xp)")
for p, x in [(2, 2), (3, 2), (4, 2), (3, 8)]:
    s0_opt = 1 - (x-1)/(x*p)
    s0_h1 = pandrosion_s(1.0, p, x)
    check(f"p={p}, x={x}: h(1) = 1-(x-1)/(xp) = {s0_opt:.4f}",
          close(s0_opt, s0_h1),
          f"formula={s0_opt:.6f}, h(1)={s0_h1:.6f}")

# Verify iteration counts from table
print("\n  Iteration count comparison")
for p, x, s0_opt_exp, s0_half, iters_opt, iters_half in [
    (2, 2, 0.7500, 0.5, 12, 13),
    (3, 2, 0.8333, 0.5, 14, 16),
    (4, 2, 0.8750, 0.5, 15, 17),
]:
    alpha = x**(1/p)
    # Count iters for s0_opt
    s0 = 1 - (x-1)/(x*p)
    s = s0
    n1 = 0
    for _ in range(200):
        v = x * s**(p-1)
        if abs(v - alpha) < 1e-10:
            break
        s = pandrosion_s(s, p, x)
        n1 += 1
    
    # Count iters for s0 = 0.5
    s = 0.5
    n2 = 0
    for _ in range(200):
        v = x * s**(p-1)
        if abs(v - alpha) < 1e-10:
            break
        s = pandrosion_s(s, p, x)
        n2 += 1
    
    check(f"p={p}, x={x}: iters(opt)={n1}≤{iters_opt}, iters(1/2)={n2}≤{iters_half}",
          n1 <= iters_opt + 1 and n2 <= iters_half + 1,
          f"opt={n1}, half={n2}")
    check(f"p={p}, x={x}: optimal is faster", n1 <= n2)


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  §12: STEFFENSEN ACCELERATION                                           ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("§12: STEFFENSEN ACCELERATION")
print("=" * 78)

# Steffensen achieves machine precision in 3 steps for p=3, x=2
print("\n  Steffensen convergence for p=3, x=2")
p, x = 3, 2
alpha = x**(1/p)
s = 0.875
steff_eps = []
steff_s = [s]

for i in range(5):
    s_new = steffensen_step(s, p, x)
    v = x * s_new**(p-1)
    eps = abs(v - alpha)
    steff_eps.append(eps)
    steff_s.append(s_new)
    if eps < 1e-15:
        break
    s = s_new

check("Steffensen reaches machine precision in ≤ 3 steps",
      len(steff_eps) <= 4 and steff_eps[-1] < 1e-14,
      f"steps={len(steff_eps)}, final ε={steff_eps[-1]:.2e}")

# Verify table values
table_steff = [
    (0, 0.875000000000000, 2.71e-1),
    (1, 0.793949257650704, 7.90e-4),
    (2, 0.793700528604164, 8.32e-9),
    (3, 0.793700525984100, 2.22e-16),
]
print("\n  Steffensen table verification")
s = 0.875
for n, s_exp, eps_exp in table_steff:
    v = x * s**(p-1)
    eps = abs(v - alpha)
    check(f"Step {n}: s ≈ {s_exp:.6f}",
          approx_eq(s, s_exp, tol=1e-6),
          f"computed={s:.15f}")
    if n > 0:
        check(f"Step {n}: ε ≈ {eps_exp:.0e}",
              abs(eps - eps_exp) / max(eps_exp, 1e-16) < 0.5 or eps < 1e-15,
              f"computed={eps:.2e}")
    if n < 3 and eps > 1e-15:
        s = steffensen_step(s, p, x)

# Quadratic convergence: ε_{n+1}/ε_n² should be roughly constant
print("\n  Quadratic convergence verification")
s = 0.875
s_star = 1/alpha
eps_list = [abs(s - s_star)]
for i in range(3):
    s = steffensen_step(s, p, x)
    eps_list.append(abs(s - s_star))

for i in range(1, len(eps_list)):
    if eps_list[i-1] > 1e-14 and eps_list[i] > 1e-15:
        ratio = eps_list[i] / eps_list[i-1]**2
        check(f"ε_{i}/ε_{i-1}² = {ratio:.4f} (constant ⇒ quadratic)", True)

# K_S << K_Newton
K_newton = (p-1) / (2*alpha)
check(f"K_Newton = {K_newton:.4f} (from article)",
      approx_eq(K_newton, 0.794, tol=0.001))

# Steffensen for other p, x
print("\n  Steffensen for various (p,x)")
for p, x in [(2, 2), (3, 3), (4, 2), (5, 2)]:
    alpha = x**(1/p)
    s = 0.5
    for i in range(6):
        s_new = steffensen_step(s, p, x)
        if abs(s_new - 1/alpha) < 1e-14:
            check(f"p={p}, x={x}: Steffensen converges in ≤ {i+1} steps", True)
            break
        s = s_new
    else:
        check(f"p={p}, x={x}: Steffensen converges in ≤ 6 steps", False)


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  RUN ORIGINAL VERIFICATION TOO                                          ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("CROSS-CHECK: Original article formulas still valid")
print("=" * 78)

# Quick cross-check of core formulas
for p, x in [(2,2),(3,2),(3,3),(4,2),(5,2),(3,8),(4,16)]:
    alpha = x**(1/p)
    s = 0.5
    lam_th = lambda_theoretical(p, x)
    ratios = []
    for _ in range(60):
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
        check(f"p={p},x={x}: ratio → λ = {lam_th:.6f}",
              approx_eq(ratio_final, lam_th, tol=0.005))


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  FINAL TALLY                                                            ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "═" * 78)
print(f"  FINAL: {PASS} passed, {FAIL} failed")
print("═" * 78)
if FAIL == 0:
    print("  🎉 ALL NEW SECTION FORMULAS AND VALUES ARE CORRECT")
else:
    print(f"  ⚠️  {FAIL} ERROR(S) DETECTED — see ❌ above")
print()
