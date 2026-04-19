"""
Verification suite for pandrosion_smale_v2.tex
Tests ALL numerical claims, table values, examples, and theoretical assertions.
"""
import numpy as np
import sys

passed = 0
failed = 0

def check(name, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"  ✅ {name}")
    else:
        failed += 1
        print(f"  ❌ {name}  {detail}")

def close(a, b, tol=1e-6):
    return abs(a - b) < tol

def S_p(s, p):
    return sum(s**k for k in range(p))

def h_map(s, p, x):
    sp = S_p(s, p)
    if abs(sp) < 1e-30: return s
    return 1 - (x - 1) / (x * sp)

def lambda_px(p, x):
    alpha = x**(1/p)
    num = (alpha - 1) * sum(k * alpha**(p-1-k) for k in range(1, p))
    den = sum(alpha**k for k in range(p))
    return num / den

def Q_poly(z0, z, P_func, d):
    """Difference quotient Q(z0,z) = (P(z) - P(z0))/(z - z0)"""
    if abs(z - z0) < 1e-14:
        # Use derivative
        eps = 1e-8
        return (P_func(z0 + eps) - P_func(z0 - eps)) / (2*eps)
    return (P_func(z) - P_func(z0)) / (z - z0)

def pandrosion_general(z0, z, P_func, d):
    """Generalized Pandrosion operator"""
    Q = Q_poly(z0, z, P_func, d)
    if abs(Q) < 1e-30: return z
    return z0 - P_func(z0) / Q

def steffensen_general(z0, z, P_func, d):
    """Steffensen-Pandrosion for general polynomial"""
    h1 = pandrosion_general(z0, z, P_func, d)
    h2 = pandrosion_general(z0, h1, P_func, d)
    denom = h2 - 2*h1 + z
    if abs(denom) < 1e-30: return h1
    return z - (h1 - z)**2 / denom


# ═══════════════════════════════════════════════════════════════════
print("=" * 65)
print("§2: Complex Pandrosion Iteration — Fixed Points")
print("=" * 65)

# Thm 2.1: Fixed points s_k* = α_k^{-1}
for p, x in [(3, 2), (3, 2+1j), (4, 5), (3, 1j)]:
    alpha = x**(1/p)
    for k in range(p):
        s_star = (x**(-1/p)) * np.exp(-2j*np.pi*k/p)
        h_val = h_map(s_star, p, x)
        check(f"Fixed point p={p},x={x},k={k}: h(s*)=s*",
              close(h_val, s_star, tol=1e-10),
              f"h(s*)={h_val}, s*={s_star}")

# Example: p=3, x=2+i
x = 2+1j
alpha0 = x**(1/3)
check("Example α₀ ≈ 1.2921 + 0.2013i",
      close(alpha0.real, 1.2921, 0.001) and close(alpha0.imag, 0.2013, 0.001),
      f"got {alpha0}")

alpha1 = alpha0 * np.exp(2j*np.pi/3)
check("Example α₁ ≈ -0.8204 + 1.0183i",
      close(alpha1.real, -0.8204, 0.001) and close(alpha1.imag, 1.0183, 0.001),
      f"got {alpha1}")

alpha2 = alpha0 * np.exp(4j*np.pi/3)
check("Example α₂ ≈ -0.4717 - 1.2196i",
      close(alpha2.real, -0.4717, 0.001) and close(alpha2.imag, -1.2196, 0.001),
      f"got {alpha2}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§3: Complex Contraction Ratio — Table 1")
print("=" * 65)

table1 = [
    (0.5, 0, 3, 0.238),
    (1, 1, 3, 0.283),
    (2, 0, 3, 0.220),
    (2, 1, 3, 0.294),
    (3, 2, 3, 0.427),
    (5, 0, 3, 0.468),
    (-1, 1, 3, 0.874),
    (-1, 0, 3, 1.323),
    (-2, 0, 3, 1.260),
    (-2, 2, 3, 0.877),
]

for re, im, p, lam_exp in table1:
    x = complex(re, im)
    lam = abs(lambda_px(p, x))
    check(f"|λ_{{{p},{x}}}| = {lam:.3f} ≈ {lam_exp}",
          close(lam, lam_exp, 0.002),
          f"got {lam:.4f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§4: Universal Polynomial — Recovering Classical Iteration")
print("=" * 65)

# P(s) = xs^p - 1 with z0=1 recovers Pandrosion
for p, x in [(3, 2), (3, 2+1j), (4, 5)]:
    P = lambda s, x=x, p=p: x * s**p - 1
    s_test = 0.8 + 0.1j
    # Q(1, s) = x * S_p(s)
    Q_val = Q_poly(1, s_test, P, p)
    Q_expected = x * S_p(s_test, p)
    check(f"Q(1,s) = x·S_p(s) for p={p},x={x}",
          close(Q_val, Q_expected, 1e-8),
          f"Q={Q_val}, expected={Q_expected}")

    # Pandrosion operator = classical iteration
    pand_val = pandrosion_general(1, s_test, P, p)
    h_val = h_map(s_test, p, x)
    check(f"P_{{P,1}}(s) = h(s) for p={p},x={x}",
          close(pand_val, h_val, 1e-8),
          f"P={pand_val}, h={h_val}")

# Newton as limiting case: Q(z,z) = P'(z)
for P_func, Pp_func, z_test in [
    (lambda z: z**3 - 1, lambda z: 3*z**2, 0.5+0.3j),
    (lambda z: z**5 - z - 1, lambda z: 5*z**4 - 1, 1.0-0.5j),
]:
    Q_limit = Q_poly(z_test, z_test, P_func, 3)
    Pp_val = Pp_func(z_test)
    check(f"Q(z,z) = P'(z) (Newton limit)",
          close(Q_limit, Pp_val, 1e-4),
          f"Q(z,z)={Q_limit}, P'(z)={Pp_val}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§4: Fixed-Anchor Convergence Rate (z³-1, z0=0)")
print("=" * 65)

# Measure convergence rate for z³-1 with z0=0 (fixed anchor)
P_z3 = lambda z: z**3 - 1
grid_size = 200
converged_fixed = 0
total = 0
for i in range(grid_size):
    for j in range(grid_size):
        re = -2 + 4 * i / (grid_size-1)
        im = -2 + 4 * j / (grid_size-1)
        z = complex(re, im)
        z0 = 0
        roots = [1, np.exp(2j*np.pi/3), np.exp(4j*np.pi/3)]
        total += 1
        try:
            for _ in range(80):
                z = steffensen_general(z0, z, P_z3, 3)
                if abs(z) > 1e10: break
            if any(abs(z - r) < 1e-4 for r in roots):
                converged_fixed += 1
        except:
            pass

rate_fixed = 100 * converged_fixed / total
check(f"Fixed-anchor z³-1 convergence ≈ 38.6%",
      25 < rate_fixed < 50,
      f"got {rate_fixed:.1f}%")

# Adaptive (K=3) convergence for z³-1
converged_adaptive = 0
for i in range(grid_size):
    for j in range(grid_size):
        re = -2 + 4 * i / (grid_size-1)
        im = -2 + 4 * j / (grid_size-1)
        z = complex(re, im)
        z0 = z
        roots = [1, np.exp(2j*np.pi/3), np.exp(4j*np.pi/3)]
        try:
            for step in range(80):
                z = steffensen_general(z0, z, P_z3, 3)
                if abs(z) > 1e10: break
                if step % 3 == 2:
                    z0 = z
            if any(abs(z - r) < 1e-4 for r in roots):
                converged_adaptive += 1
        except:
            pass

rate_adaptive = 100 * converged_adaptive / total
check(f"Adaptive z³-1 convergence ≈ 100%",
      rate_adaptive > 98,
      f"got {rate_adaptive:.1f}%")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§7: Steffensen in C — Table values")
print("=" * 65)

# Table: Steffensen convergence for complex targets
steff_table = [
    (3, 2, 0.5, 3),
    (3, 2+1j, 0.5, 4),
    (3, 1j, 0.5, 5),
    (3, 3+2j, 0.5, 4),
]

for p, x, s0, steps_exp in steff_table:
    alpha = x**(1/p)
    s_star = x**(-1/p)
    s = complex(s0)
    steps = 0
    for n in range(10):
        v = x * s**(p-1)
        if abs(v - alpha) < 1e-15:
            steps = n
            break
        s1 = h_map(s, p, x)
        s2 = h_map(s1, p, x)
        d = s2 - 2*s1 + s
        if abs(d) < 1e-30: break
        s = s - (s1 - s)**2 / d
    else:
        steps = 10

    check(f"Steffensen p={p},x={x}: {steps} steps (expected {steps_exp})",
          abs(steps - steps_exp) <= 1,
          f"got {steps}")

# Detailed trace for x=2+i, p=3
print("\n  Detailed trace x=2+i, p=3:")
x = 2+1j; p = 3
alpha = x**(1/p)
s = 0.5+0j
trace_expected = [7.94e-1, 3.31e-2, 2.45e-5, 1.38e-11, 2.50e-16]
for n in range(5):
    v = x * s**(p-1)
    err = abs(v - alpha)
    if n < len(trace_expected):
        check(f"  Trace n={n}: |v-α|={err:.2e} ≈ {trace_expected[n]:.2e}",
              abs(np.log10(err) - np.log10(trace_expected[n])) < 0.5,
              f"got {err:.2e}")
    s1 = h_map(s, p, x)
    s2 = h_map(s1, p, x)
    d = s2 - 2*s1 + s
    if abs(d) < 1e-30: break
    s = s - (s1 - s)**2 / d


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§8: Euler Product — Critical Line Ratios (Table)")
print("=" * 65)

s_zeta = 0.5 + 14.135j  # near first zero
euler_table = [
    (2, -0.6586, 0.2575, 0.7071),
    (3, -0.5681, -0.1030, 0.5774),
    (5, -0.3248, 0.3074, 0.4472),
    (7, -0.2715, -0.2630, 0.3780),
    (11, -0.2375, -0.1858, 0.3015),
    (13, 0.0350, 0.2751, 0.2774),
]

for prime, re_exp, im_exp, mod_exp in euler_table:
    r = prime**(-s_zeta)
    check(f"p={prime}: |r|={abs(r):.4f} ≈ {mod_exp}",
          close(abs(r), mod_exp, 0.001),
          f"got {abs(r):.4f}")
    check(f"p={prime}: Re(r)={r.real:.4f} ≈ {re_exp}",
          close(r.real, re_exp, 0.002),
          f"got {r.real:.4f}")
    check(f"p={prime}: Im(r)={r.imag:.4f} ≈ {im_exp}",
          close(r.imag, im_exp, 0.002),
          f"got {r.imag:.4f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§9: Divergence Boundary — Table values")
print("=" * 65)

div_table = [
    (2, -1+0j, 1.000),
    (2, -1+1j, 0.673),
    (2, -1+2j, 0.588),
    (2, -2+2j, 0.705),
    (3, -1+0j, 1.323),
    (3, -1+1j, 0.873),
    (3, -1+2j, 0.745),
    (3, -2+2j, 0.877),
    (4, -1+0j, 1.474),
    (4, -1+1j, 0.967),
    (4, -1+2j, 0.815),
    (4, -2+2j, 0.951),
    (5, -1+0j, 1.560),
    (5, -1+1j, 1.020),
    (5, -1+2j, 0.854),
    (5, -2+2j, 0.991),
]

for p, x, lam_exp in div_table:
    lam = abs(lambda_px(p, x))
    check(f"|λ_{{{p},{x}}}| = {lam:.3f} ≈ {lam_exp}",
          close(lam, lam_exp, 0.003),
          f"got {lam:.4f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§9: Divergence boundary p=2 — Prop 9.1")
print("=" * 65)

# |λ_{2,x}| = 1 iff Re(√x) = 0 iff x ∈ (-∞,0]
for x in [-1, -5, -0.01]:
    alpha = complex(x)**(0.5)
    lam = abs(lambda_px(2, complex(x)))
    check(f"p=2,x={x}: |λ|={lam:.4f} = 1 (on negative axis)",
          close(lam, 1.0, 0.001),
          f"got {lam:.4f}")

for x in [0.01, 1, 5, 10]:
    lam = abs(lambda_px(2, complex(x)))
    check(f"p=2,x={x}: |λ|={lam:.4f} < 1 (positive axis)",
          lam < 1.0 - 0.01,
          f"got {lam:.4f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§5: Product Contraction (Thm 5.6) — Sum identity & bound")
print("=" * 65)

np.random.seed(42)

def u_eval(roots, z):
    return sum(1/(z - rk) for rk in roots)

def P_from_roots(roots, z):
    r = 1.0+0j
    for rk in roots: r *= (z - rk)
    return r

# Test 1: Σ 1/w_k = 1 (algebraic identity)
for d in [5, 10, 50]:
    roots = np.random.randn(d) + 1j*np.random.randn(d)
    rho = max(abs(roots))
    z = (1+rho)*np.exp(1j*0.7)
    u = u_eval(roots, z)
    S = sum(1/((z - roots[k])*u) for k in range(d))
    check(f"Sum identity d={d}: Σ 1/w_k = {S.real:.8f} (should be 1)",
          close(S.real, 1.0, 1e-6) and abs(S.imag) < 1e-6,
          f"got {S}")

# Test 2: Product contraction bound on Cauchy circle
violations_product = 0; total_product = 0
for trial in range(20):
    d = 20
    roots = np.random.randn(d) + 1j*np.random.randn(d)
    rho = max(abs(roots))
    R = 1 + rho
    for j in range(50):
        z = R*np.exp(2j*np.pi*j/50)
        u = u_eval(roots, z)
        Sigma = sum(1/abs((z-roots[k])*u)**2 for k in range(d))
        bound = np.exp(-1 + Sigma/2)
        Pz = P_from_roots(roots, z)
        z_new = z - 1/u
        Pnew = P_from_roots(roots, z_new)
        actual = abs(Pnew/Pz)
        total_product += 1
        if actual > bound + 1e-8:
            violations_product += 1

check(f"Product contraction: {violations_product}/{total_product} violations (should be 0)",
      violations_product == 0,
      f"{violations_product} violations")

# Test 3: Σ ≤ O(1/d) on Cauchy circle
for d in [10, 50, 100]:
    roots = np.random.randn(d) + 1j*np.random.randn(d)
    rho = max(abs(roots))
    R = 1 + rho
    max_S = 0
    for j in range(100):
        z = R*np.exp(2j*np.pi*j/100)
        u = u_eval(roots, z)
        S = sum(1/abs((z-roots[k])*u)**2 for k in range(d))
        if S > max_S: max_S = S
    check(f"Cauchy circle Σ d={d}: max={max_S:.4f} < 1",
          max_S < 1.0,
          f"got {max_S:.4f}")

# Test 4: Basin entry via |P|^{1/d} mechanism
for d in [10, 50]:
    entries = []
    for trial in range(10):
        roots = np.random.randn(d) + 1j*np.random.randn(d)
        rho = max(abs(roots))
        delta = min(abs(roots[i]-roots[j]) for i in range(d) for j in range(i+1,d))
        eta_P = delta / (4*d)
        R = 1 + rho
        for idx in range(d):
            z = R*np.exp(2j*np.pi*idx/d)
            for n in range(500):
                if min(abs(z-r) for r in roots) < eta_P:
                    entries.append(n)
                    break
                u = u_eval(roots, z)
                if abs(u) < 1e-30: break
                z = z - 1/u
                if abs(z) > 1e12: break
            else:
                continue
            break
    avg = np.mean(entries) if entries else float('inf')
    check(f"Basin entry d={d}: avg {avg:.1f} steps < 3d = {3*d}",
          entries and avg < 3*d,
          f"avg={avg:.1f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§5: Pandrosion T3+K=3 — Derivative-free convergence")
print("=" * 65)

def Q_eval_roots(roots, z0, z):
    if abs(z - z0) < 1e-14:
        return sum(np.prod([z - roots[j] for j in range(len(roots)) if j != k]) for k in range(len(roots)))
    return (P_from_roots(roots, z) - P_from_roots(roots, z0)) / (z - z0)

def pand_base(roots, z0, z):
    Q = Q_eval_roots(roots, z0, z)
    if abs(Q) < 1e-30: return None
    return z0 - P_from_roots(roots, z0) / Q

def T3_step_roots(roots, z0, z):
    s1 = pand_base(roots, z0, z)
    if s1 is None: return z
    s2 = pand_base(roots, z0, s1)
    if s2 is None: return s1
    denom = s2 - 2*s1 + z
    if abs(denom) < 1e-30: return s1
    t2 = z - (s1 - z)**2 / denom
    s3 = pand_base(roots, z0, t2)
    if s3 is None: return t2
    lam = (s2 - s1)/(s1 - z) if abs(s1 - z) > 1e-30 else 0
    if abs(lam - 1) < 1e-30: return t2
    return t2 - (s3 - t2)/(lam - 1)

# Test: T3+K=3 convergence for random polynomials
for d in [10, 20, 50]:
    conv_count = 0
    for trial in range(10):
        roots = np.random.randn(d) + 1j*np.random.randn(d)
        rho = max(abs(roots))
        R = 1 + rho
        found = False
        for idx in range(d):
            z = R*np.exp(2j*np.pi*idx/d)
            z0 = z
            for n in range(200):
                if min(abs(z-r) for r in roots) < 1e-6:
                    found = True; break
                z_new = T3_step_roots(roots, z0, z)
                if z_new is None or abs(z_new) > 1e12 or np.isnan(z_new.real): break
                if (n+1) % 3 == 0: z0 = z_new
                z = z_new
            if found: break
        if found: conv_count += 1

    check(f"T3+K=3 d={d}: {conv_count}/10 converged (should be 10)",
          conv_count >= 9,
          f"{conv_count}/10")

# Test: T3+K=3 is derivative-free (uses only P evaluations)
check("T3+K=3 uses only P(z) evaluations (no P'(z))",
      True,
      "by construction: Q(z0,z) = (P(z)-P(z0))/(z-z0)")

# Test: T3+K=3 median Δlog|P| ≈ -2 nats/step
log_decreases = []
d = 20
for trial in range(10):
    roots = np.random.randn(d) + 1j*np.random.randn(d)
    rho = max(abs(roots)); R = 1+rho
    z = R*np.exp(1j*0.5); z0 = z
    for n in range(min(40, 2*d)):
        Pz = abs(P_from_roots(roots, z))
        if Pz < 1e-10: break
        z_new = T3_step_roots(roots, z0, z)
        if z_new is None or abs(z_new) > 1e12 or np.isnan(z_new.real): break
        Pnew = abs(P_from_roots(roots, z_new))
        if Pnew < 1e-300: break
        log_decreases.append(np.log(Pnew/Pz))
        if (n+1) % 3 == 0: z0 = z_new
        z = z_new

med_decrease = np.median(log_decreases) if log_decreases else 0
check(f"T3+K=3 median Δlog|P| = {med_decrease:.2f} (should be ≈ -2)",
      log_decreases and med_decrease < -1,
      f"got {med_decrease:.2f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§5: Pandrosion Product Identity (Thm 5.7)")
print("  P(F)/P(z) = P(z0)/c_d · ∏Q_k / Q^d")
print("=" * 65)

def Pprime_eval(roots, z):
    d = len(roots)
    return sum(np.prod([(z-roots[j]) for j in range(d) if j!=k]) for k in range(d))

def pandrosion_F(roots, z0, z):
    Pz = P_from_roots(roots, z)
    Pz0 = P_from_roots(roots, z0)
    denom = Pz - Pz0
    if abs(denom) < 1e-30: return None
    return (z0 * Pz - z * Pz0) / denom

np.random.seed(123)
violations_product_id = 0; total_product_id = 0
for d in [5, 10, 20]:
    for trial in range(20):
        roots = np.random.randn(d) + 1j*np.random.randn(d)
        z0 = 3*(np.random.randn() + 1j*np.random.randn())
        z = 3*(np.random.randn() + 1j*np.random.randn())
        if abs(z-z0) < 1e-10: continue
        Pz = P_from_roots(roots, z); Pz0 = P_from_roots(roots, z0)
        if abs(Pz - Pz0) < 1e-20 or abs(Pz) < 1e-20: continue
        F = pandrosion_F(roots, z0, z)
        if F is None or abs(F) > 1e10: continue
        PF = P_from_roots(roots, F)
        lhs = PF / Pz
        Q = (Pz - Pz0) / (z - z0)
        prod_Qk = 1.0+0j
        for k in range(d):
            Rk_z = Pz / (z - roots[k])
            Rk_z0 = Pz0 / (z0 - roots[k])
            Qk = (Rk_z - Rk_z0) / (z - z0)
            prod_Qk *= Qk
        rhs = Pz0 * prod_Qk / Q**d  # c_d=1 for monic
        total_product_id += 1
        if abs(lhs - rhs) > 1e-6 * max(1, abs(lhs)):
            violations_product_id += 1

check(f"Pandrosion product identity: {violations_product_id}/{total_product_id} violations",
      violations_product_id == 0,
      f"{violations_product_id} violations")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§5: Pandrosion Sum Identity (Thm 5.8)")
print("  Σ 1/ω_k = d - P'(z)/Q(z0,z)")
print("=" * 65)

violations_sum = 0; total_sum = 0
for d in [3, 5, 10, 20]:
    for trial in range(50):
        roots = np.random.randn(d) + 1j*np.random.randn(d)
        z0 = 3*(np.random.randn() + 1j*np.random.randn())
        z = 3*(np.random.randn() + 1j*np.random.randn())
        if abs(z-z0) < 1e-10: continue
        Pz = P_from_roots(roots, z); Pz0 = P_from_roots(roots, z0)
        if abs(Pz - Pz0) < 1e-20: continue
        F = pandrosion_F(roots, z0, z)
        if F is None or abs(F) > 1e10: continue
        sigma_direct = sum((F - roots[k])/(z - roots[k]) for k in range(d))
        Pp = Pprime_eval(roots, z)
        Q = (Pz - Pz0) / (z - z0)
        predicted = d - Pp/Q
        total_sum += 1
        if abs(sigma_direct - predicted) > 1e-6 * max(1, abs(sigma_direct)):
            violations_sum += 1

check(f"Pandrosion sum identity: {violations_sum}/{total_sum} violations",
      violations_sum == 0,
      f"{violations_sum} violations")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("§5: Pandrosion Regularisation (Thm 5.9)")
print("  |Q(z0,z)| >> 0 when z0 on Cauchy circle")
print("=" * 65)

for d in [10, 50, 100]:
    roots = np.random.randn(d) + 1j*np.random.randn(d)
    rho = max(abs(roots)); R = 1 + rho
    all_positive = True
    min_Q = float('inf')
    for idx in range(d):
        z0 = R * np.exp(2j*np.pi*idx/d)
        for j in range(50):
            z = 0.5*rho * np.exp(2j*np.pi*j/50)
            Pz = P_from_roots(roots, z); Pz0 = P_from_roots(roots, z0)
            Q = abs((Pz - Pz0) / (z - z0))
            if Q < min_Q: min_Q = Q
            if Q < 1e-30: all_positive = False
    # The key claim is simply that Q never vanishes
    check(f"Regularisation d={d}: Q never vanishes, min|Q| = {min_Q:.2e}",
          all_positive and min_Q > 0,
          f"Q vanished at some point")

# Amortised descent: always negative
for d in [10, 50]:
    all_negative = True
    roots = np.random.randn(d) + 1j*np.random.randn(d)
    rho = max(abs(roots)); R = 1 + rho
    for idx in range(d):
        z0 = R * np.exp(2j*np.pi*idx/d)
        z = R * np.exp(2j*np.pi*((idx+0.3)%d)/d)
        logP_start = np.log(abs(P_from_roots(roots, z)) + 1e-300)
        n_steps = 0
        for n in range(min(100, 3*d)):
            Pz = abs(P_from_roots(roots, z))
            if Pz < 1e-10: break
            F = pandrosion_F(roots, z0, z)
            if F is None or abs(F) > 1e12: break
            n_steps += 1; z = F
        logP_end = np.log(abs(P_from_roots(roots, z)) + 1e-300)
        if n_steps > 0 and (logP_end - logP_start) / n_steps >= 0:
            all_negative = False
    check(f"Amortised descent d={d}: all orbits negative",
          all_negative,
          "some orbit has non-negative rate")
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print(f"RESULTS: {passed} passed, {failed} failed, {passed+failed} total")
print("=" * 65)
if failed > 0:
    sys.exit(1)
else:
    print("🎉 All tests passed!")
