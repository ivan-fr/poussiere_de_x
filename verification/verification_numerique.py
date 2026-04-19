"""
Script de verification numerique des tables de l'article
Poussiere de Pandrosion generalisee
"""
import numpy as np

def S_p(s, p):
    return sum(s**k for k in range(p))

def pandrosion_s(s, p, x):
    return 1 - (x - 1) / (x * S_p(s, p))

def lambda_theoretical(p, x):
    alpha = x**(1/p)
    num = (alpha - 1) * sum(k * alpha**(p-1-k) for k in range(1, p))
    den = sum(alpha**k for k in range(p))
    return num / den

def lambda_numerical(p, x, s0=0.9):
    """Compute lambda as the limiting ratio eps_{n+1}/eps_n."""
    alpha = x**(1/p)
    s = s0
    ratios = []
    for _ in range(80):
        v = x * s**(p-1)
        eps = abs(v - alpha)
        s_new = pandrosion_s(s, p, x)
        v_new = x * s_new**(p-1)
        eps_new = abs(v_new - alpha)
        if eps > 1e-13 and eps_new > 1e-13:
            ratios.append(eps_new / eps)
        s = s_new
    # Take the last 5 stable ratios
    return np.mean(ratios[-5:]) if len(ratios) >= 5 else None

def h_prime(p, x):
    """Derivee de h en s* — verifie le coefficient theorique."""
    alpha = x**(1/p)
    s_star = 1.0 / alpha
    eps = 1e-7
    return (pandrosion_s(s_star + eps, p, x) - pandrosion_s(s_star - eps, p, x)) / (2*eps)

# ── Table 1 : coefficients lambda ───────────────────────────────────────────
print("=" * 75)
print("TABLE 1 : Coefficients de poussiere lambda_{p,x}")
print("=" * 75)
print(f"{'p':>3} {'x':>4} {'alpha':>10} {'lambda_theory':>16} {'lambda_numeric':>16} {'h_prime':>10} {'OK?':>6}")
print("-" * 75)

configs = [
    (2, 2), (2, 3), (2, 5),
    (3, 2), (3, 3), (3, 8),
    (4, 2), (4, 16),
    (5, 2),
]

all_ok = True
for p, x in configs:
    alpha = x**(1/p)
    lt  = lambda_theoretical(p, x)
    ln  = lambda_numerical(p, x, s0=0.9 if x <= 3 else 0.7)
    hp  = h_prime(p, x)
    ok  = ln is not None and abs(lt - ln) / lt < 0.01 and abs(lt - hp) / lt < 1e-4
    if not ok:
        all_ok = False
    ln_str = f"{ln:.8f}" if ln is not None else "N/A"
    tag = "OK" if ok else "ERREUR"
    print(f"{p:>3} {x:>4} {alpha:>10.6f} {lt:>16.8f} {ln_str:>16} {hp:>10.6f} {tag:>6}")

print("-" * 75)
print(f"Resultat global : {'TOUTES LES VALEURS SONT CORRECTES' if all_ok else 'CERTAINES VALEURS ECHOUENT'}")

# ── Formules closes : p=2 et p=3 ────────────────────────────────────────────
print()
print("=" * 75)
print("VERIFICATION DES FORMULES CLOSES")
print("=" * 75)

print("\np=2 : lambda_{2,x} = (sqrt(x)-1)/(sqrt(x)+1)")
for x in [2, 3, 5, 10]:
    alpha = x**0.5
    closed = (alpha - 1) / (alpha + 1)
    general = lambda_theoretical(2, x)
    print(f"  x={x:>2}: closed={closed:.8f}, general={general:.8f}, diff={abs(closed-general):.2e}  {'OK' if abs(closed-general)<1e-12 else 'ECART'}")

print("\np=3 : lambda_{3,x} = (cbrt(x)-1)*(cbrt(x)+2) / (cbrt(x^2)+cbrt(x)+1)")
for x in [2, 3, 8, 27]:
    alpha = x**(1/3)
    closed = (alpha - 1) * (alpha + 2) / (alpha**2 + alpha + 1)
    general = lambda_theoretical(3, x)
    print(f"  x={x:>2}: closed={closed:.8f}, general={general:.8f}, diff={abs(closed-general):.2e}  {'OK' if abs(closed-general)<1e-12 else 'ECART'}")

# ── Convergence lineaire : verification du taux ──────────────────────────────
print()
print("=" * 75)
print("VERIFICATION DE LA CONVERGENCE LINEAIRE (vs quadratique Newton)")
print("=" * 75)

p, x = 3, 2
alpha = x**(1/p)
s = 0.875
lam = lambda_theoretical(p, x)

print(f"\nPandrosion p={p}, x={x}, lambda_theorique={lam:.6f}")
print(f"{'n':>4} {'v_n':>14} {'eps_n':>16} {'eps_{n+1}/eps_n':>18} {'lin_predict':>15}")
print("-" * 72)

s = 0.875
eps_prev = None
for n in range(12):
    v = x * s**(p-1)
    eps = abs(v - alpha)
    if eps_prev is not None and eps_prev > 1e-14:
        ratio = eps / eps_prev
        pred  = lam * eps_prev
        print(f"{n:>4} {v:>14.10f} {eps:>16.6e} {ratio:>18.6f} {pred:>15.6e}")
    else:
        print(f"{n:>4} {v:>14.10f} {eps:>16.6e}")
    eps_prev = eps
    s = pandrosion_s(s, p, x)

print()
print(f"\nNewton p={p}, x={x}")
print(f"{'n':>4} {'u_n':>14} {'eps_n':>16} {'eps_{n+1}/eps_n^2':>20}")
print("-" * 55)
u = 1.5
eps_prev = None
for n in range(7):
    eps = abs(u - alpha)
    if eps_prev is not None and eps_prev > 1e-14:
        ratio2 = eps / eps_prev**2
        print(f"{n:>4} {u:>14.10f} {eps:>16.6e} {ratio2:>20.6f}")
    else:
        print(f"{n:>4} {u:>14.10f} {eps:>16.6e}")
    eps_prev = eps
    u = (2*u + 2/u**2) / 3

# ── Test du point fixe ────────────────────────────────────────────────────────
print()
print("=" * 75)
print("VERIFICATION DES POINTS FIXES")
print("=" * 75)
for p, x in [(2,2),(3,2),(3,3),(4,2),(5,2)]:
    s_star = x**(-1/p)
    s_fp = pandrosion_s(s_star, p, x)
    v_star = x * s_star**(p-1)
    print(f"p={p}, x={x}: s*={s_star:.8f}, h(s*)={s_fp:.8f}, diff={abs(s_fp-s_star):.2e}, v*={v_star:.8f} (alpha={x**(1/p):.8f})")

print("\nToutes les verifications sont completes.")
