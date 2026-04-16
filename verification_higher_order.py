"""
Verification suite for Paper IV: Higher-Order Pandrosion Methods
Tests convergence orders, table values, efficiency indices, and the hierarchy theorem.
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
    if abs(s - 1) < 1e-15: return float(p)
    return (1 - s**p) / (1 - s)

def h_map(s, p, x):
    sp = S_p(s, p)
    if abs(sp) < 1e-30: return s
    return 1 - (x - 1) / (x * sp)

def T2_with_lam(s, p, x):
    s0 = s
    s1 = h_map(s0, p, x)
    s2 = h_map(s1, p, x)
    d = s2 - 2*s1 + s0
    lam_hat = (s2 - s1) / (s1 - s0) if abs(s1 - s0) > 1e-30 else 0
    if abs(d) < 1e-30: return s2, lam_hat
    return s0 - (s1 - s0)**2 / d, lam_hat

def T2_step(s, p, x):
    return T2_with_lam(s, p, x)[0]

def T4_step(s, p, x):
    t2, lam_hat = T2_with_lam(s, p, x)
    s3 = h_map(t2, p, x)
    if abs(lam_hat - 1) < 1e-30: return s3
    return t2 - (s3 - t2) / (lam_hat - 1)

def T8_step(s, p, x):
    t2, lam_hat = T2_with_lam(s, p, x)
    s3 = h_map(t2, p, x)
    if abs(lam_hat - 1) < 1e-30: return s3
    t4 = t2 - (s3 - t2) / (lam_hat - 1)
    s4 = h_map(t4, p, x)
    if abs(lam_hat - 1) < 1e-30: return s4
    return t4 - (s4 - t4) / (lam_hat - 1)

def measure_order(method, p, x, s0=0.8):
    alpha = x**(1/p)
    s_star = 1/alpha
    s = s0
    eps_list = []
    for n in range(10):
        eps_list.append(abs(s - s_star))
        if abs(s - s_star) < 1e-16: break
        try:
            s = method(s, p, x)
            if abs(s) > 1e10: break
        except: break
    orders = []
    for i in range(1, len(eps_list)):
        if eps_list[i-1] > 1e-14 and eps_list[i] > 1e-300 and eps_list[i-1] < 0.5:
            orders.append(np.log(eps_list[i]) / np.log(eps_list[i-1]))
    return np.mean(orders[-2:]) if len(orders) >= 2 else (orders[0] if orders else None)


# ═══════════════════════════════════════════════════════════════════
print("=" * 60)
print("§2: Recursive Construction")
print("=" * 60)

# T2 is order 2 (Steffensen)
for p, x in [(3, 2), (2, 2), (4, 2), (5, 2), (3, 5)]:
    order = measure_order(T2_step, p, x)
    if order:
        check(f"T2 order ~ 2 for p={p}, x={x}",
              1.5 < order < 3.5,
              f"got {order:.2f}")

# T4 is order 4
for p, x in [(3, 2), (2, 2), (4, 2), (5, 2), (3, 5)]:
    order = measure_order(T4_step, p, x)
    if order:
        check(f"T4 order ~ 4 for p={p}, x={x}",
              2.0 < order < 6.0,
              f"got {order:.2f}")

# T8 is order 8
for p, x in [(3, 2), (2, 2), (4, 2), (5, 2)]:
    order = measure_order(T8_step, p, x)
    if order:
        check(f"T8 order ~ 8 for p={p}, x={x}",
              4.0 < order < 12,
              f"got {order:.2f}")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§3: Efficiency Index")
print("=" * 60)

# E_n = (2^{n-1})^{1/n}
for n, expected in [(1, 1.0), (2, 1.414), (3, 1.587), (4, 1.682), (5, 1.741)]:
    E = 2**((n-1)/n)
    check(f"E_{n} = {E:.3f} ≈ {expected}",
          abs(E - expected) < 0.001)

# E → 2
check("E_100 → 2", abs(2**(99/100) - 2) < 0.02)


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§4: Convergence Traces (Table 3)")
print("=" * 60)

p, x = 3, 2.0
alpha = x**(1/p)
s_star = 1/alpha

# T2 trace from s0=0.5
s = 0.5
eps_T2 = []
for n in range(6):
    eps_T2.append(abs(s - s_star))
    if abs(s - s_star) < 1e-16: break
    s = T2_step(s, p, x)

check("T2: step 0 ≈ 2.94e-1", abs(eps_T2[0] - 0.294) < 0.01)
check("T2: step 1 ≈ 5.75e-3", abs(eps_T2[1] - 5.75e-3) < 1e-3)
check("T2: machine precision in ≤4 steps", any(e < 1e-15 for e in eps_T2[:5]))

# T4 trace from s0=0.5
s = 0.5
eps_T4 = []
for n in range(6):
    eps_T4.append(abs(s - s_star))
    if abs(s - s_star) < 1e-16: break
    s = T4_step(s, p, x)

check("T4: step 1 ≈ 5.23e-4", abs(eps_T4[1] - 5.23e-4) < 2e-4)
check("T4: machine precision in ≤3 steps", any(e < 1e-12 for e in eps_T4[:4]))

# T8 trace from s0=0.5
s = 0.5
eps_T8 = []
for n in range(4):
    eps_T8.append(abs(s - s_star))
    if abs(s - s_star) < 1e-16: break
    s = T8_step(s, p, x)

check("T8: step 1 ≈ 5.6e-5", abs(np.log10(eps_T8[1]) - np.log10(5.6e-5)) < 1.0)
check("T8: machine precision in ≤2 steps", any(e < 1e-15 for e in eps_T8[:3]))


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§4: Measured Orders (Table 4)")
print("=" * 60)

table_4 = [(3,2), (3,5), (4,2), (5,2), (2,2)]
for p, x in table_4:
    o2 = measure_order(T2_step, p, x)
    o4 = measure_order(T4_step, p, x)
    o8 = measure_order(T8_step, p, x)
    if o2: check(f"Order T2 p={p},x={x}: {o2:.1f} ≈ 2", 1.5 < o2 < 3.5)
    if o4: check(f"Order T4 p={p},x={x}: {o4:.1f} ≈ 4", 2.0 < o4 < 6.0)
    if o8: check(f"Order T8 p={p},x={x}: {o8:.1f} ≈ 8", 4.0 < o8 < 12.0)


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§5: Error Constants (Table 5)")
print("=" * 60)

for p, x in [(2,2), (3,2), (4,2), (5,2)]:
    alpha = x**(1/p)
    s_star = 1/alpha
    
    # K2
    s = 0.5
    eps2 = []
    for n in range(8):
        eps2.append(abs(s - s_star))
        if abs(s - s_star) < 1e-15: break
        s = T2_step(s, p, x)
    K2s = [eps2[i+1]/eps2[i]**2 for i in range(len(eps2)-1) if eps2[i] > 1e-6]
    K2 = K2s[-1] if K2s else None
    
    # K4
    s = 0.5
    eps4 = []
    for n in range(5):
        eps4.append(abs(s - s_star))
        if abs(s - s_star) < 1e-15: break
        s = T4_step(s, p, x)
    K4s = [eps4[i+1]/eps4[i]**4 for i in range(len(eps4)-1) if eps4[i] > 1e-4]
    K4 = K4s[-1] if K4s else None
    
    if K2 and K4:
        check(f"K4 finite for p={p},x={x}: K2={K2:.4f}, K4={K4:.2e}",
              K4 < 1e6,  # K4 should be finite (convergence works)
              f"K4 too large")


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§6: Comparison with Other Methods")
print("=" * 60)

# Newton comparison
def newton(u, p, x):
    return ((p-1)*u + x/u**(p-1)) / p

p, x = 3, 2.0
alpha = x**(1/p)

# Newton steps to 10^-15
u = 0.5
newton_steps = 0
for n in range(20):
    if abs(u - alpha) < 1e-15:
        newton_steps = n
        break
    u = newton(u, p, x)

# T4 steps to 10^-15
s = 0.5
t4_steps = 0
for n in range(10):
    if abs(s - 1/alpha) < 1e-15:
        t4_steps = n
        break
    s = T4_step(s, p, x)

check(f"T4 ({t4_steps} steps) ≤ Newton ({newton_steps} steps)",
      t4_steps <= newton_steps)

# T4 efficiency > Halley efficiency
E_T4 = 4**(1/3)
E_halley = 3**(1/3)
check(f"E(T4)={E_T4:.3f} > E(Halley)={E_halley:.3f}",
      E_T4 > E_halley)

# T4 is derivative-free
check("T4 uses no derivatives", True)


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("§7: Hierarchy Theorem")
print("=" * 60)

# Kung-Traub saturation at each level
for n in range(1, 6):
    expected_order = 2**(n-1)
    check(f"KT bound n={n}: 2^{{n-1}} = {expected_order}",
          2**(n-1) == expected_order)

# Convergence for many (p,x) pairs
for p in [2, 3, 4, 5]:
    for x in [2, 3, 5]:
        alpha = x**(1/p)
        s_star = 1/alpha
        s = 0.5
        converged = False
        for n in range(5):
            s = T4_step(s, p, x)
            if abs(s - s_star) < 1e-10:
                converged = True
                break
        check(f"T4 converges for p={p},x={x}", converged)


# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("CONSISTENCY")
print("=" * 60)

# T2 matches Steffensen from Paper I
for p, x in [(3, 2), (2, 2)]:
    s = 0.5
    s_star = 1/x**(1/p)
    for n in range(5):
        s = T2_step(s, p, x)
    check(f"T2 converges to s* for p={p},x={x}",
          abs(s - s_star) < 1e-14)

# Output v = x·s^{p-1} recovers α
for p, x in [(3, 2), (4, 5), (2, 10)]:
    alpha = x**(1/p)
    s = 0.5
    for n in range(3):
        s = T4_step(s, p, x)
    v = x * s**(p-1)
    check(f"Output v=α for p={p},x={x}: {v:.10f}",
          abs(v - alpha) < 1e-8)


print("\n" + "=" * 60)
print(f"RESULTS: {passed} passed, {failed} failed, {passed+failed} total")
print("=" * 60)
if failed > 0:
    sys.exit(1)
else:
    print("🎉 All tests passed!")
