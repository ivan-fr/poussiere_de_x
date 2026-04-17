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
print(f"RESULTS: {passed} passed, {failed} failed, {passed+failed} total")
print("=" * 65)
if failed > 0:
    sys.exit(1)
else:
    print("🎉 All tests passed!")
