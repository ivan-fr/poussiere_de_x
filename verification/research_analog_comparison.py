"""Fair comparison: Pandrosion vs Newton in analog (same eval budget)"""
import numpy as np
np.random.seed(42)

def S_p_noisy(s, p, sigma):
    result = 1.0; power = s
    for k in range(1, p):
        if k > 1: power = power * s * (1 + np.random.normal(0, sigma))
        result += power * (1 + np.random.normal(0, sigma))
    return result * (1 + np.random.normal(0, sigma))

def h_noisy(s, p, x, sigma):
    sp = S_p_noisy(s, p, sigma)
    if abs(sp) < 1e-30: return s
    num = (x - 1) * (1 + np.random.normal(0, sigma))
    den = x * sp * (1 + np.random.normal(0, sigma))
    return (1 - num / den) * (1 + np.random.normal(0, sigma))

def T2_noisy(s, p, x, sigma):
    s0, s1 = s, h_noisy(s, p, x, sigma)
    s2 = h_noisy(s1, p, x, sigma)
    d = s2 - 2*s1 + s0
    if abs(d) < 1e-30: return s2
    return s0 - (s1 - s0)**2 / d

def T4_noisy(s, p, x, sigma):
    s0, s1 = s, h_noisy(s, p, x, sigma)
    s2 = h_noisy(s1, p, x, sigma)
    d1, d2 = s1 - s0, s2 - s1
    lam = d2 / d1 if abs(d1) > 1e-30 else 0
    d = (s2 - 2*s1 + s0) * (1 + np.random.normal(0, sigma))
    if abs(d) < 1e-30: return s2
    t2 = s0 - d1**2 / d
    s3 = h_noisy(t2, p, x, sigma)
    if abs(lam - 1) < 1e-30: return s3
    return t2 - (s3 - t2) / (lam - 1) * (1 + np.random.normal(0, sigma))

def newton_noisy(u, p, x, sigma, n_steps=1):
    for _ in range(n_steps):
        up1 = u**(p-1) * (1 + np.random.normal(0, sigma))
        fp = p * up1 * (1 + np.random.normal(0, sigma))
        fu = (u * up1 - x) * (1 + np.random.normal(0, sigma))
        if abs(fp) < 1e-30: return u
        u = (u - fu / fp) * (1 + np.random.normal(0, sigma))
    return u

p, x = 3, 2.0
alpha = x**(1/p)
N = 5000

print("FAIR COMPARISON: same evaluation budget")
print("=" * 65)

for sigma in [1e-3, 1e-4]:
    print(f"\nsigma = {sigma:.0e}:")
    
    # Newton 1-step from u=1.0, 2 evals (f, f')
    out = [newton_noisy(1.0, p, x, sigma, 1) for _ in range(N)]
    std, bias = np.std(out), abs(np.mean(out) - alpha)
    bits = -np.log2(std/alpha) if std > 0 else 0
    print(f"  Newton 1-step u=1.0 (2 evals): {bits:.1f} bits, bias={bias:.2e}")
    
    # Newton 1-step from u=0.5, 2 evals (same bad start)
    out = [newton_noisy(0.5, p, x, sigma, 1) for _ in range(N)]
    std, bias = np.std(out), abs(np.mean(out) - alpha)
    bits = -np.log2(std/alpha) if std > 0 else 0
    print(f"  Newton 1-step u=0.5 (2 evals): {bits:.1f} bits, bias={bias:.2e}")
    
    # T2 from s=0.5, 2 evals of h
    out = [x * T2_noisy(0.5, p, x, sigma)**(p-1) for _ in range(N)]
    std, bias = np.std(out), abs(np.mean(out) - alpha)
    bits = -np.log2(std/alpha) if std > 0 else 0
    print(f"  Pandrosion T2 s=0.5  (2 evals): {bits:.1f} bits, bias={bias:.2e}")
    
    # T4 single pass from s=0.5, 3 evals of h
    out = [x * T4_noisy(0.5, p, x, sigma)**(p-1) for _ in range(N)]
    std, bias = np.std(out), abs(np.mean(out) - alpha)
    bits = -np.log2(std/alpha) if std > 0 else 0
    print(f"  Pandrosion T4 s=0.5  (3 evals): {bits:.1f} bits, bias={bias:.2e}")
    
    # Newton 3-step from u=1.0, 6 evals (unfair but show it)
    out = [newton_noisy(1.0, p, x, sigma, 3) for _ in range(N)]
    std, bias = np.std(out), abs(np.mean(out) - alpha)
    bits = -np.log2(std/alpha) if std > 0 else 0
    print(f"  Newton 3-step u=1.0 (6 evals): {bits:.1f} bits, bias={bias:.2e}")

print("\n" + "=" * 65)
print("MATH ERROR (no noise):")
print("=" * 65)

def S_p(s, p):
    return (1-s**p)/(1-s) if abs(s-1)>1e-15 else float(p)
def h(s, p, x):
    return 1-(x-1)/(x*S_p(s,p))

s_star = 1/alpha

# T2
s1 = h(0.5, p, x); s2 = h(s1, p, x)
t2 = 0.5 - (s1-0.5)**2 / (s2 - 2*s1 + 0.5)
print(f"  T2(0.5):  |s - s*| = {abs(t2 - s_star):.2e}")

# T4
lam = (s2-s1)/(s1-0.5)
s3 = h(t2, p, x)
t4 = t2 - (s3-t2)/(lam-1)
print(f"  T4(0.5):  |s - s*| = {abs(t4 - s_star):.2e}")

# Newton from 1.0
un = ((p-1)*1.0 + x/1.0**(p-1))/p
print(f"  Newton(1.0): |u - a| = {abs(un - alpha):.2e}")

# Newton from 0.5
un = ((p-1)*0.5 + x/0.5**(p-1))/p
print(f"  Newton(0.5): |u - a| = {abs(un - alpha):.2e}")
