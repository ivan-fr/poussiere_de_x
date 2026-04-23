"""
Numerically verify monotonicity of lambda_{p,x} in (p, alpha).
Also verifies the closed-form identities:
  B_p(alpha) = (alpha-1) N_p D_p + p D_p - alpha^p N_p
            = D_p(alpha)^2 - alpha^p N_p(alpha)
            = sum_{j=1}^{p} j alpha^{j-1}
  P_p(alpha) = numerator of d/dalpha log(lambda_{p, alpha^p}) after
              common-denom reduction by (alpha-1) N_p D_p
            = p * D_p'(alpha)
            = p * sum_{j=1}^{p-1} j alpha^{j-1}
symbolically for p up to 20 using SymPy, if available.
"""
import math
try:
    import sympy as sp
    _HAS_SYMPY = True
except ImportError:
    _HAS_SYMPY = False


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


# ---------------------------------------------------------------------------
# (5) Symbolic verification of B_p and P_p closed forms
# ---------------------------------------------------------------------------
if _HAS_SYMPY:
    print()
    print("=" * 72)
    print("(5) Symbolic verification of closed forms (requires SymPy):")
    alpha = sp.Symbol("alpha")
    def N(p): return sum(k * alpha ** (p - 1 - k) for k in range(1, p))
    def D(p): return sum(alpha ** k for k in range(p))
    def Bp(p):
        return sp.expand((alpha - 1) * N(p) * D(p) + p * D(p) - alpha ** p * N(p))
    def Pp(p):
        logderiv = sp.diff(sp.log((alpha - 1) * N(p) / D(p)), alpha)
        return sp.expand(sp.simplify(logderiv * (alpha - 1) * N(p) * D(p)))
    def closed_B(p): return sp.expand(sum(j * alpha ** (j - 1) for j in range(1, p + 1)))
    def closed_P(p): return sp.expand(p * sum(j * alpha ** (j - 1) for j in range(1, p)))

    all_ok = True
    for p in range(2, 21):
        b_ok = sp.simplify(Bp(p) - closed_B(p)) == 0
        p_ok = sp.simplify(Pp(p) - closed_P(p)) == 0
        all_ok = all_ok and b_ok and p_ok
    print(f"  B_p = sum_{{j=1}}^{{p}} j alpha^{{j-1}}  : verified for 2 <= p <= 20 -> {all_ok}")
    print(f"  P_p = p * sum_{{j=1}}^{{p-1}} j alpha^{{j-1}}  : verified for 2 <= p <= 20 -> {all_ok}")
else:
    print()
    print("SymPy not available; skipping symbolic closed-form verification.")
