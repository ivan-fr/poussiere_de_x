"""
Verification suite for: Pandrosion Iteration in the Complex Plane
Tests all theorems, tables, and numerical claims in the article.
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

def S_p(s, p):
    if abs(s - 1) < 1e-15:
        return complex(p)
    return (1 - s**p) / (1 - s)

def pandrosion_complex(s, p, x):
    sp = S_p(s, p)
    if abs(sp) < 1e-30:
        return s
    return 1 - (x - 1) / (x * sp)

def lambda_complex(p, x):
    alpha = x**(1/p)
    num = (alpha - 1) * sum(k * alpha**(p-1-k) for k in range(1, p))
    den = sum(alpha**k for k in range(p))
    return num / den

def steffensen_step(s, p, x):
    s0 = s
    s1 = pandrosion_complex(s0, p, x)
    s2 = pandrosion_complex(s1, p, x)
    denom = s2 - 2*s1 + s0
    if abs(denom) < 1e-30:
        return s2
    return s0 - (s1 - s0)**2 / denom


# ═══════════════════════════════════════════════════════════════════
print("=" * 60)
print("§2: Complex Pandrosion iteration — Fixed points")
print("=" * 60)

# Theorem 2.1: p fixed points are s*_k = α_k^{-1}
for p in [2, 3, 4, 5]:
    for x in [2+1j, 3+2j, 1j, 0.5+0j, 5+0j]:
        for k in range(p):
            alpha_k = x**(1/p) * np.exp(2j*np.pi*k/p)
            s_star_k = 1 / alpha_k
            # Verify it's a fixed point
            h_s = pandrosion_complex(s_star_k, p, x)
            check(f"Fixed point p={p}, x={x}, k={k}",
                  abs(h_s - s_star_k) < 1e-10,
                  f"|h(s*) - s*| = {abs(h_s - s_star_k):.2e}")

# Fixed points lie on circle of radius |x|^{-1/p}
for p, x in [(3, 2+1j), (4, 3+2j), (2, 1j)]:
    radius = abs(x)**(-1/p)
    for k in range(p):
        s_star_k = (1/x)**(1/p) * np.exp(-2j*np.pi*k/p)
        # Use x^{-1/p} principal root and rotate
        alpha_k = x**(1/p) * np.exp(2j*np.pi*k/p)
        s_star_k = 1 / alpha_k
        check(f"Circle radius p={p}, x={x}, k={k}",
              abs(abs(s_star_k) - radius) < 1e-10,
              f"|s*| = {abs(s_star_k):.6f}, expected {radius:.6f}")

# Output v* = alpha_k
for p, x in [(3, 2+1j), (3, 1j), (2, 3+2j)]:
    for k in range(p):
        alpha_k = x**(1/p) * np.exp(2j*np.pi*k/p)
        s_star_k = 1 / alpha_k
        v_star = x * s_star_k**(p-1)
        check(f"Output v*=α_k p={p}, x={x}, k={k}",
              abs(v_star - alpha_k) < 1e-10,
              f"|v*-α| = {abs(v_star - alpha_k):.2e}")

# Example values: p=3, x=2+i
alpha_0 = (2+1j)**(1/3)
check("Example: α_0 ≈ 1.2921+0.2013i",
      abs(alpha_0 - (1.2921 + 0.2013j)) < 0.001)
alpha_1 = alpha_0 * np.exp(2j*np.pi/3)
check("Example: α_1 ≈ -0.8204+1.0183i",
      abs(alpha_1 - (-0.8204 + 1.0183j)) < 0.001)
alpha_2 = alpha_0 * np.exp(4j*np.pi/3)
check("Example: α_2 ≈ -0.4717-1.2196i",
      abs(alpha_2 - (-0.4717 - 1.2196j)) < 0.001)


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§3: Complex contraction ratio")
print("=" * 60)

# Theorem 3.1: λ formula matches h'(s*) by finite difference
for p, x in [(3, 2+1j), (3, 1j), (2, 3+2j), (4, 1+3j), (5, 2+2j), (3, 0.5+0j)]:
    lam = lambda_complex(p, x)
    s_star = 1 / x**(1/p)
    epsilon = 1e-8
    h_plus = pandrosion_complex(s_star + epsilon, p, x)
    h_minus = pandrosion_complex(s_star - epsilon, p, x)
    h_prime_numerical = (h_plus - h_minus) / (2*epsilon)
    check(f"λ = h'(s*) p={p}, x={x}",
          abs(lam - h_prime_numerical) < 1e-5,
          f"formula={lam:.6f}, numerical={h_prime_numerical:.6f}")

# Table 1: |λ_{3,x}| values
table_1 = [
    (0.5, 0, 0.238), (1, 1, 0.283), (2, 0, 0.220), (2, 1, 0.294),
    (3, 2, 0.427), (5, 0, 0.468), (-1, 1, 0.874), (-1, 0, 1.323),
    (-2, 0, 1.260), (-2, 2, 0.877),
]
for re, im, expected in table_1:
    x = complex(re, im)
    lam = lambda_complex(3, x)
    check(f"Table 1: |λ_{{3,{x}}}| ≈ {expected}",
          abs(abs(lam) - expected) < 0.002,
          f"got {abs(lam):.3f}")

# Conjugate symmetry: λ(x̄) = λ̄(x)
for p, x in [(3, 2+1j), (3, 3+2j), (4, 1+3j)]:
    lam_x = lambda_complex(p, x)
    lam_xbar = lambda_complex(p, np.conj(x))
    check(f"Conjugate symmetry p={p}, x={x}",
          abs(lam_xbar - np.conj(lam_x)) < 1e-10)

# Positive real axis ⊂ C_p
for p in [2, 3, 4, 5]:
    for x_real in [0.5, 1.5, 2, 5, 10, 100]:
        lam = lambda_complex(p, complex(x_real))
        check(f"x={x_real} > 0 → |λ|<1 (p={p})",
              abs(lam) < 1,
              f"|λ| = {abs(lam):.4f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§4: Basins of attraction")
print("=" * 60)

# Theorem 4.1: (0,1) ⊂ B_0 for x with Re(x) > 0
for x in [2+1j, 3+2j, 5+0j, 1+0.5j]:
    alpha_0 = x**(1/3)
    s_star_0 = 1 / alpha_0
    for s0_real in [0.1, 0.3, 0.5, 0.7, 0.9]:
        s = complex(s0_real)
        converged = False
        for n in range(200):
            s = pandrosion_complex(s, 3, x)
            if abs(s - s_star_0) < 1e-10:
                converged = True
                break
        check(f"(0,1)⊂B_0: s0={s0_real}, x={x}",
              converged,
              "did not converge to s*_0")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§5: Steffensen acceleration in ℂ")
print("=" * 60)

# Theorem 5.1: quadratic convergence
steffensen_cases = [
    (3, 2+0j, 3), (3, 2+1j, 4), (3, 1j, 5), (3, 3+2j, 4), (4, 1+3j, 5),
]
for p, x, expected_steps in steffensen_cases:
    alpha = x**(1/p)
    s = 0.5 + 0j
    steps_to_precision = None
    for n in range(8):
        v = x * s**(p-1)
        eps = abs(v - alpha)
        if eps < 1e-15:
            steps_to_precision = n
            break
        s = steffensen_step(s, p, x)
    check(f"Steffensen p={p}, x={x}: {expected_steps} steps",
          steps_to_precision is not None and steps_to_precision <= expected_steps + 1,
          f"got {steps_to_precision} steps")

# Verify quadratic rate: eps_{n+1} ≈ K * eps_n^2
for p, x in [(3, 2+1j), (3, 3+2j)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    s = 0.5 + 0j
    epsilons = []
    for n in range(5):
        epsilons.append(abs(s - s_star))
        if abs(s - s_star) < 1e-15:
            break
        s = steffensen_step(s, p, x)
    # Check quadratic: eps_{n+1}/eps_n^2 should be roughly constant
    if len(epsilons) >= 3:
        ratios = [epsilons[i+1] / epsilons[i]**2 for i in range(len(epsilons)-1) if epsilons[i] > 1e-6]
        if len(ratios) >= 2:
            # Ratios should stabilize (ignoring near-machine-precision values)
            check(f"Quadratic rate p={p}, x={x}",
                  all(r < 100 for r in ratios),
                  f"ratios = {[f'{r:.2f}' for r in ratios]}")

# Example trace: x=2+i
x_trace = 2+1j
alpha_trace = x_trace**(1/3)
s = 0.5 + 0j
expected_eps = [7.94e-1, 3.31e-2, 2.45e-5, 1.38e-11, 2.50e-16]
for n in range(5):
    v = x_trace * s**(3-1)  # p=3
    eps = abs(v - alpha_trace)
    if n < len(expected_eps):
        # Order of magnitude match
        if expected_eps[n] > 1e-15:
            check(f"Trace x=2+i, step {n}: |v-α| ≈ {expected_eps[n]:.2e}",
                  abs(np.log10(eps) - np.log10(expected_eps[n])) < 1.0,
                  f"got {eps:.2e}")
    if eps < 1e-16:
        break
    s = steffensen_step(s, 3, x_trace)


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§6: Euler product connection")
print("=" * 60)

def sieve(n):
    is_prime = [True] * (n+1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5)+1):
        if is_prime[i]:
            for j in range(i*i, n+1, i):
                is_prime[j] = False
    return [i for i in range(2, n+1) if is_prime[i]]

primes = sieve(1000)

# Proposition 6.1: Euler product identity
# Verify partial Euler product ≈ ζ(s) for Re(s) > 1
def zeta_direct(s, terms=50000):
    return sum(n**(-s) for n in range(1, terms+1))

for s_val in [2, 3, 4, 2+1j, 3+2j]:
    z_direct = zeta_direct(s_val, 50000)
    z_euler = 1.0 + 0j
    for p in primes[:168]:  # 1000th prime = 997
        z_euler *= 1 / (1 - p**(-s_val))
    check(f"Euler product ≈ ζ({s_val})",
          abs(z_direct - z_euler) / abs(z_direct) < 0.01,
          f"relative error = {abs(z_direct - z_euler)/abs(z_direct):.2e}")

# Each factor = lim S_m(p^{-s})
for prime in [2, 3, 5, 7]:
    s_val = 2
    ratio = prime**(-s_val)
    factor_geometric = sum(ratio**k for k in range(100))  # S_100
    factor_exact = 1 / (1 - ratio)
    check(f"S_100({prime}^{{-2}}) ≈ 1/(1-{prime}^{{-2}})",
          abs(factor_geometric - factor_exact) < 1e-10)

# Proposition 6.2: |p^{-s}| = p^{-1/2} on critical line
t = 14.134725
s_crit = 0.5 + t*1j
for prime in primes[:20]:
    ratio = prime**(-s_crit)
    expected_mod = prime**(-0.5)
    check(f"|{prime}^{{-s}}| = {prime}^{{-1/2}} on critical line",
          abs(abs(ratio) - expected_mod) < 1e-10)

# Table 5: specific ratio values
table_5 = [
    (2, -0.6586, 0.2575, 0.7071),
    (3, -0.5681, -0.1030, 0.5774),
    (5, -0.3248, 0.3074, 0.4472),
    (7, -0.2715, -0.2630, 0.3780),
    (11, -0.2375, -0.1858, 0.3015),
    (13, 0.0350, 0.2751, 0.2774),
]
for prime, re_exp, im_exp, mod_exp in table_5:
    ratio = prime**(-s_crit)
    check(f"Table 5: {prime}^{{-s}} ≈ {re_exp}+{im_exp}i",
          abs(ratio.real - re_exp) < 0.001 and abs(ratio.imag - im_exp) < 0.001,
          f"got {ratio.real:.4f}+{ratio.imag:.4f}i")
    check(f"Table 5: |{prime}^{{-s}}| ≈ {mod_exp}",
          abs(abs(ratio) - mod_exp) < 0.001)


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§7: Divergence boundary")
print("=" * 60)

# Proposition 7.1: |λ_{2,x}| = 1 iff Re(√x) = 0 iff x ∈ (-∞,0]
# On negative real axis: |λ| = 1
for x_neg in [-0.5, -1, -2, -5]:
    alpha = complex(x_neg)**(0.5)
    check(f"|λ_{{2,{x_neg}}}| = 1 (negative real)",
          abs(abs(lambda_complex(2, complex(x_neg))) - 1.0) < 0.01,
          f"|λ| = {abs(lambda_complex(2, complex(x_neg))):.4f}")

# Off negative real: |λ| < 1
for x in [2+0j, 1+1j, -1+1j, -2+2j, 0.5+0j]:
    if x.imag != 0 or x.real > 0:
        lam = lambda_complex(2, x)
        check(f"|λ_{{2,{x}}}| < 1",
              abs(lam) < 1,
              f"|λ| = {abs(lam):.4f}")

# Table 6: divergence boundary values
table_6 = [
    (2, -1+0j, 1.000), (2, -1+1j, 0.673), (2, -1+2j, 0.588), (2, -2+2j, 0.705),
    (3, -1+0j, 1.323), (3, -1+1j, 0.873), (3, -1+2j, 0.745), (3, -2+2j, 0.877),
    (4, -1+0j, 1.474), (4, -1+1j, 0.967), (4, -1+2j, 0.815), (4, -2+2j, 0.951),
    (5, -1+0j, 1.560), (5, -1+1j, 1.020), (5, -1+2j, 0.854), (5, -2+2j, 0.991),
]
for p, x, expected in table_6:
    lam = lambda_complex(p, x)
    check(f"Table 6: |λ_{{{p},{x}}}| ≈ {expected}",
          abs(abs(lam) - expected) < 0.002,
          f"got {abs(lam):.3f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("CONSISTENCY CHECKS")
print("=" * 60)

# Real case recovery: complex λ matches real λ for real x > 0
for p in [2, 3, 4, 5]:
    for x_real in [2, 3, 5, 10]:
        alpha_real = x_real**(1/p)
        lam_real_formula = ((alpha_real - 1) * sum(k * alpha_real**(p-1-k) for k in range(1, p))
                           / sum(alpha_real**k for k in range(p)))
        lam_complex_val = lambda_complex(p, complex(x_real))
        check(f"Real recovery: p={p}, x={x_real}",
              abs(lam_complex_val.imag) < 1e-10 and abs(lam_complex_val.real - lam_real_formula) < 1e-10)

# x=1 → λ=0 (trivial case)
for p in [2, 3, 4, 5]:
    lam = lambda_complex(p, 1.0 + 0j)
    check(f"x=1 → λ=0 (p={p})",
          abs(lam) < 1e-10,
          f"|λ| = {abs(lam):.2e}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print(f"RESULTS: {passed} passed, {failed} failed, {passed+failed} total")
print("=" * 60)

if failed > 0:
    sys.exit(1)
else:
    print("🎉 All tests passed!")
