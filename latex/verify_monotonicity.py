"""
Numerically verify monotonicity of lambda_{p,x} in (p, alpha).
Check several proof strategies.
"""
import math


def lam(p, alpha):
    num = (alpha - 1.0) * sum(k * alpha ** (p - 1 - k) for k in range(1, p))
    den = sum(alpha ** k for k in range(p))
    return num / den


# (1) Monotonicity in alpha for fixed p
print("Monotonicity in alpha (fixed p):")
for p in [2, 3, 5, 10, 20]:
    vals = [(alpha, lam(p, alpha))
            for alpha in [1.01, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 100.0]]
    inc = all(vals[i+1][1] > vals[i][1] for i in range(len(vals) - 1))
    print(f"  p={p}: strictly increasing = {inc}")

# (2) Monotonicity in p for fixed alpha
print("\nMonotonicity in p (fixed alpha):")
for alpha in [1.1, 1.5, 2.0, 3.0, 10.0, 100.0]:
    vals = [(p, lam(p, alpha)) for p in range(2, 51)]
    inc = all(vals[i+1][1] > vals[i][1] for i in range(len(vals) - 1))
    print(f"  alpha={alpha}: strictly increasing in p = {inc}")

# (3) Check the log-derivative sign in alpha
print("\nLog-derivative d/d alpha ln(lambda) > 0 ?  (p=3 closed form)")
print("  For p=3: d/d alpha ln lambda = 1/(alpha-1) + 1/(alpha+2) - (2 alpha+1)/(alpha^2+alpha+1)")
for a in [1.1, 1.5, 2.0, 5.0, 20.0, 100.0]:
    v = 1/(a-1) + 1/(a+2) - (2*a+1)/(a**2+a+1)
    print(f"  alpha={a}: value = {v:+.6f}")

# (4) alpha^p - 1 = (alpha-1) D_p   and   p D_p vs alpha^p N_p
print("\nSign of p D_p - alpha^p N_p  (for step p -> p+1 argument):")
for p in [2, 3, 4, 5, 7, 10]:
    for alpha in [1.5, 2.0, 5.0]:
        Dp = sum(alpha ** k for k in range(p))
        Np = sum(k * alpha ** (p - 1 - k) for k in range(1, p))
        diff = p * Dp - alpha ** p * Np
        print(f"  p={p}, alpha={alpha}: p D_p - alpha^p N_p = {diff:+.4f}")
