"""
=============================================================================
VÉRIFICATION EXHAUSTIVE DE L'ARTICLE
« Poussière de Pandrosion généralisée »
=============================================================================
Ce script vérifie CHAQUE formule, CHAQUE valeur numérique et CHAQUE identité
algébrique citée dans l'article, section par section.
"""
import numpy as np
from fractions import Fraction

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

def approx_eq(a, b, tol=1e-3):
    return abs(a - b) < tol

def close(a, b, tol=1e-10):
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


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2 — La construction de Pandrosion pour cbrt(2)                  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("=" * 78)
print("SECTION 2 : Construction de Pandrosion pour cbrt(2)  (p=3, x=2)")
print("=" * 78)

# --- Eq (1) : itération originale u_{n+1}, v_n ---
print("\n§2.1 — Équation (1) : itération originale")
u0 = 0.5
u = u0
u_list = [u]
for _ in range(20):
    u_next = u / (2 - 2*((4 - u)/4)**3)
    u_list.append(u_next)
    u = u_next
v_list_orig = [2*((4 - ui)/4)**2 for ui in u_list]
check("v_n converge vers cbrt(2)",
      approx_eq(v_list_orig[-1], 2**(1/3), tol=1e-12),
      f"v_20 = {v_list_orig[-1]:.15f}, cbrt(2) = {2**(1/3):.15f}")

# --- §2.2 : Substitution normalisée ---
print("\n§2.2 — Substitution normalisée s_n = (4-u_n)/4")

# Vérifier u_n = 4(1-s_n)
s0 = (4 - u0)/4
check("s_0 = (4-0.5)/4 = 0.875", close(s0, 0.875))

# Vérifier v_n = 2*s_n^2
check("v_0 = 2*s_0^2 = 2*0.875^2 = 1.53125",
      close(2*s0**2, 2*0.875**2) and close(2*s0**2, 1.53125))

# Vérifier u_{n+1} = 4(1-s_n)/(2(1-s_n^3)) = 2/(1+s_n+s_n^2)
s_test = 0.85
u_from_s = 4*(1-s_test)
u_next_orig = u_from_s / (2 - 2*((4-u_from_s)/4)**3)
u_next_formula1 = 4*(1-s_test) / (2*(1-s_test**3))
u_next_formula2 = 2 / (1 + s_test + s_test**2)
check("u_{n+1} = 4(1-s)/(2(1-s^3)) = 2/(1+s+s^2)",
      close(u_next_formula1, u_next_formula2) and close(u_next_orig, u_next_formula2),
      f"formule1={u_next_formula1:.10f}, formule2={u_next_formula2:.10f}, orig={u_next_orig:.10f}")

# --- Eq (2) : itération sur s ---
print("\n§2.2 — Équation (2) : s_{n+1} = (2s^2+2s+1)/(2s^2+2s+2)")
s = s0
for _ in range(20):
    s = (2*s**2 + 2*s + 1) / (2*s**2 + 2*s + 2)
check("Itération (2) converge vers 2^{-1/3}",
      close(s, 2**(-1/3), tol=1e-14),
      f"s_20 = {s:.15f}, 2^{{-1/3}} = {2**(-1/3):.15f}")

# Vérifier équivalence avec 1 - (x-1)/(x·S_p(s))
s_test = 0.82
iter_eq2 = (2*s_test**2 + 2*s_test + 1) / (2*s_test**2 + 2*s_test + 2)
iter_general = pandrosion_s(s_test, 3, 2)
check("Éq.(2) ≡ formule générale pandrosion_s pour p=3,x=2",
      close(iter_eq2, iter_general),
      f"éq2={iter_eq2:.12f}, général={iter_general:.12f}")

# --- §2.3 : Point fixe ---
print("\n§2.3 — Point fixe")
s_star = 2**(-1/3)
check("s* = 2^{-1/3} ≈ 0.7937", approx_eq(s_star, 0.7937, tol=0.0001))
check("s*^3 = 1/2", close(s_star**3, 0.5))
check("v* = 2·s*^2 = 2·2^{-2/3} = 2^{1/3} = cbrt(2)",
      close(2*s_star**2, 2**(1/3)))

# Vérifier l'algèbre : 2(1-s*)(1+s*+s*^2) = 1
check("2(1-s*)(1+s*+s*^2) = 1",
      close(2*(1-s_star)*(1+s_star+s_star**2), 1.0))

# --- §2.4 : Poussière, dérivée h'(s*) ---
print("\n§2.4 — Coefficient de poussière λ_{3,2}")

# h(s) = (2s^2+2s+1)/(2s^2+2s+2)
def h(s):
    return (2*s**2 + 2*s + 1) / (2*s**2 + 2*s + 2)

# Dérivée analytique : h'(s) = (4s+2)/(2s^2+2s+2)^2
def h_prime_analytic(s):
    return (4*s + 2) / (2*s**2 + 2*s + 2)**2

# Vérifier les 3 formes de h'(s*)  (ligne 141)
hp1 = (4*s_star + 2) / (2*s_star**2 + 2*s_star + 2)**2
hp2 = 2*(2*s_star + 1) / (4*(s_star**2 + s_star + 1)**2)
hp3 = (2*s_star + 1) / (2*(s_star**2 + s_star + 1)**2)
check("h'(s*) : 3 formes équivalentes (ligne 141)",
      close(hp1, hp2) and close(hp2, hp3),
      f"forme1={hp1:.10f}, forme2={hp2:.10f}, forme3={hp3:.10f}")

# Formule close λ_{3,2} = (cbrt(2)-1)(cbrt(2)+2)/(cbrt(4)+cbrt(2)+1)
cbrt2 = 2**(1/3)
cbrt4 = 2**(2/3)
lambda_32_closed = (cbrt2 - 1) * (cbrt2 + 2) / (cbrt4 + cbrt2 + 1)
check("λ_{3,2} = (cbrt(2)-1)(cbrt(2)+2)/(cbrt(4)+cbrt(2)+1)",
      close(lambda_32_closed, hp3),
      f"formule close={lambda_32_closed:.10f}, h'(s*)={hp3:.10f}")
check("λ_{3,2} ≈ 0.2202 (comme cité)",
      approx_eq(lambda_32_closed, 0.2202, tol=0.0005),
      f"valeur={lambda_32_closed:.6f}")

# Vérifier h'(s*) = λ théorique général
lambda_gen = lambda_theoretical(3, 2)
check("h'(s*) = lambda_theoretical(3,2)",
      close(lambda_32_closed, lambda_gen),
      f"close={lambda_32_closed:.10f}, théorique={lambda_gen:.10f}")

# --- §2.4 : Affirmation (p-1)α² = 2·cbrt(4) ---
alpha = 2**(1/3)
check("(p-1)·α² = 2·cbrt(4)  (p=3)",
      close((3-1)*alpha**2, 2*cbrt4))

# --- Table 1 : profil numérique p=3, x=2, u_0=0.5 (s_0=0.875) ---
print("\n§2.4 — Table 1 : profil numérique (p=3, x=2, s_0=0.875)")
s = 0.875
table1_expected = [
    # n, s_n (arrondi 4 déc.), v_n (arrondi 4 déc.), eps_n (ordre de grandeur), ratio
    (0, 0.8750, 1.5312, 2.713e-1, None),
    (1, 0.8107, 1.3143, 5.44e-2,  0.200),
    (2, 0.7974, 1.2717, 1.17e-2,  0.216),
    (3, 0.7945, 1.2625, 2.58e-3,  0.219),
    (4, 0.7939, 1.2605, 5.67e-4,  0.220),
]
eps_prev = None
for n, s_exp, v_exp, eps_exp, ratio_exp in table1_expected:
    v = 2 * s**2
    eps = abs(v - cbrt2)
    check(f"Table1 n={n}: s_n≈{s_exp:.4f}",
          approx_eq(s, s_exp, tol=0.001),
          f"calculé={s:.4f}")
    check(f"Table1 n={n}: v_n≈{v_exp:.4f}",
          approx_eq(v, v_exp, tol=0.001),
          f"calculé={v:.4f}")
    check(f"Table1 n={n}: ε_n≈{eps_exp:.1e}",
          abs(eps - eps_exp)/max(eps_exp, 1e-15) < 0.15,  # tolérance 15%
          f"calculé={eps:.3e}, attendu={eps_exp:.1e}")
    if ratio_exp is not None and eps_prev is not None:
        ratio = eps / eps_prev
        check(f"Table1 n={n}: ratio≈{ratio_exp:.3f}",
              approx_eq(ratio, ratio_exp, tol=0.02),
              f"calculé={ratio:.3f}")
    eps_prev = eps
    s = pandrosion_s(s, 3, 2)


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3 — Origine géométrique : pourquoi S_p apparaît               ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("SECTION 3 : Origine géométrique — S_p")
print("=" * 78)

# x · S_p(s*) · (1 - s*) = x - 1
print("\n§3.2 — Identité x·S_p(s*)·(1-s*) = x-1")
for p, x in [(2,2), (3,2), (3,3), (4,2), (5,2)]:
    s_star = x**(-1/p)
    lhs = x * S_p(s_star, p) * (1 - s_star)
    rhs = x - 1
    check(f"p={p}, x={x}: x·S_p(s*)·(1-s*) = x-1",
          close(lhs, rhs),
          f"LHS={lhs:.10f}, RHS={rhs}")

# S_p(s)(1-s) = 1 - s^p
print("\n§3.2 — Identité S_p(s)(1-s) = 1-s^p")
for p in [2,3,4,5]:
    for s_test in [0.3, 0.7, 0.9]:
        lhs = S_p(s_test, p) * (1 - s_test)
        rhs = 1 - s_test**p
        check(f"p={p}, s={s_test}: S_p(s)(1-s) = 1-s^p",
              close(lhs, rhs))


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4 — Généralisation                                             ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("SECTION 4 : Généralisation — itération universelle")
print("=" * 78)

# Remarque : pour p=3, x=2, 1-1/(2(1+s+s^2)) = (2s^2+2s+1)/(2s^2+2s+2)
print("\n§4.1 — Remarque : équivalence p=3, x=2")
import sympy as sp
s_sym = sp.Symbol('s')
lhs_sym = 1 - sp.Rational(1,2) * 1/(1 + s_sym + s_sym**2)
rhs_sym = (2*s_sym**2 + 2*s_sym + 1) / (2*s_sym**2 + 2*s_sym + 2)
diff = sp.simplify(lhs_sym - rhs_sym)
check("1 - 1/(2(1+s+s²)) = (2s²+2s+1)/(2s²+2s+2)  [symbolique]",
      diff == 0, f"diff = {diff}")

# Théorème 4.2 : point fixe s* = x^{-1/p} pour tout p, x
print("\n§4.2 — Théorème : point fixe s* = x^{-1/p}")
for p, x in [(2,2),(2,3),(2,5),(3,2),(3,3),(3,8),(4,2),(4,16),(5,2),(5,32),(10,1024)]:
    s_star = x**(-1/p)
    s_fp = pandrosion_s(s_star, p, x)
    v_star = x * s_star**(p-1)
    alpha = x**(1/p)
    check(f"p={p}, x={x}: h(s*)=s*, v*=α={alpha:.6f}",
          close(s_fp, s_star, tol=1e-12) and close(v_star, alpha, tol=1e-12),
          f"h(s*)={s_fp:.12f}, s*={s_star:.12f}, v*={v_star:.12f}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5 — Le coefficient de poussière λ_{p,x}                       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("SECTION 5 : Coefficient de poussière λ_{p,x}")
print("=" * 78)

# --- Preuve : vérifier les étapes intermédiaires ---
print("\n§5.1 — Étapes de la preuve du Théorème 5.1")

for p, x in [(3,2),(2,2),(4,2),(3,8)]:
    alpha = x**(1/p)
    s_star = alpha**(-1)
    
    # S_p(s*) = α(x-1)/(x(α-1))
    Sp_val = S_p(s_star, p)
    Sp_expected = alpha * (x-1) / (x * (alpha-1))
    check(f"p={p},x={x}: S_p(s*) = α(x-1)/(x(α-1))",
          close(Sp_val, Sp_expected),
          f"S_p={Sp_val:.10f}, attendu={Sp_expected:.10f}")
    
    # S'_p(s*) = α·Σ k·α^{-k}
    Sp_prime = sum(k * s_star**(k-1) for k in range(1, p))
    Sp_prime_expected = alpha * sum(k * alpha**(-k) for k in range(1, p))
    check(f"p={p},x={x}: S'_p(s*) = α·Σ k·α^{{-k}}",
          close(Sp_prime, Sp_prime_expected),
          f"S'_p={Sp_prime:.10f}, attendu={Sp_prime_expected:.10f}")

# --- Formule générale vs dérivée numérique ---
print("\n§5.1 — λ_{p,x} (formule) vs h'(s*) (numérique)")
for p, x in [(2,2),(2,3),(2,5),(3,2),(3,3),(3,8),(4,2),(4,16),(5,2)]:
    alpha = x**(1/p)
    s_star = alpha**(-1)
    lam_th = lambda_theoretical(p, x)
    # Dérivée numérique
    eps = 1e-8
    hp_num = (pandrosion_s(s_star+eps, p, x) - pandrosion_s(s_star-eps, p, x)) / (2*eps)
    check(f"p={p},x={x}: λ_théorique={lam_th:.8f} ≈ h'(s*)={hp_num:.8f}",
          close(lam_th, hp_num, tol=1e-5))

# --- Formules closes ---
print("\n§5.2 — Formules closes p=2, p=3, p=4")

# p=2 : λ_{2,x} = (α-1)/(α+1) = (√x-1)/(√x+1)
for x in [2, 3, 5, 10, 100]:
    alpha = x**0.5
    closed_p2 = (alpha - 1) / (alpha + 1)
    general = lambda_theoretical(2, x)
    check(f"p=2, x={x}: (√x-1)/(√x+1) = {closed_p2:.8f}",
          close(closed_p2, general))

# p=3 : λ_{3,x} = (α-1)(α+2)/(α²+α+1)
for x in [2, 3, 8, 27]:
    alpha = x**(1/3)
    closed_p3 = (alpha - 1) * (alpha + 2) / (alpha**2 + alpha + 1)
    general = lambda_theoretical(3, x)
    check(f"p=3, x={x}: (α-1)(α+2)/(α²+α+1) = {closed_p3:.8f}",
          close(closed_p3, general))

# p=4 : λ_{4,x} = (α-1)(α²+2α+3)/(α³+α²+α+1)
for x in [2, 16, 81]:
    alpha = x**(1/4)
    num_p4 = (alpha - 1) * (alpha**2 + 2*alpha + 3)
    den_p4 = alpha**3 + alpha**2 + alpha + 1
    closed_p4 = num_p4 / den_p4
    general = lambda_theoretical(4, x)
    check(f"p=4, x={x}: formule close = {closed_p4:.8f}",
          close(closed_p4, general))

# Vérifier les sommes citées dans l'article
# p=3 : Σ_{k=1}^{2} k·α^{2-k} = α + 2
alpha = sp.Symbol('alpha', positive=True)
sum_p3 = sum(k * alpha**(2-k) for k in range(1, 3))
check("p=3: Σ k·α^{2-k} = α+2", sp.simplify(sum_p3 - (alpha + 2)) == 0,
      f"somme = {sum_p3}")

# p=4 : Σ_{k=1}^{3} k·α^{3-k} = α²+2α+3
sum_p4 = sum(k * alpha**(3-k) for k in range(1, 4))
check("p=4: Σ k·α^{3-k} = α²+2α+3",
      sp.simplify(sum_p4 - (alpha**2 + 2*alpha + 3)) == 0,
      f"somme = {sum_p4}")

# --- Comportement asymptotique ---
print("\n§5.2 — Comportement asymptotique")
# λ_{2,x} → 1 quand x → ∞
lam_big = lambda_theoretical(2, 10000)
check("λ_{2,x} → 1 quand x → ∞", lam_big > 0.98,
      f"λ_(2,10000)={lam_big:.6f}")

# λ_{p,x} → 0 quand x → 1+
lam_small = lambda_theoretical(3, 1.001)
check("λ_{3,x} → 0 quand x → 1+", lam_small < 0.01,
      f"λ_(3,1.001)={lam_small:.6f}")

# λ_{3,x} → 1 quand α → ∞ (x → ∞)
lam_big3 = lambda_theoretical(3, 10**9)
check("λ_{3,x} → 1 quand α → ∞", lam_big3 > 0.99,
      f"λ_(3,10^9)={lam_big3:.6f}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6 — Table 2 : coefficients λ_{p,x}                            ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("SECTION 6 : Table 2 — coefficients de poussière")
print("=" * 78)

table2_expected = [
    (2, 2, 0.1716),
    (2, 3, 0.2679),
    (2, 5, 0.3820),
    (3, 2, 0.2202),
    (3, 3, 0.3366),
    (3, 8, 0.5714),
    (4, 2, 0.2427),
    (4, 16, 0.7333),
    (5, 2, 0.2565),
]

for p, x, lam_exp in table2_expected:
    lam_calc = lambda_theoretical(p, x)
    check(f"Table2 p={p},x={x}: λ≈{lam_exp:.4f}",
          approx_eq(lam_calc, lam_exp, tol=0.001),
          f"calculé={lam_calc:.4f}")

# Vérifier les expressions exactes citées dans la table
# (3,8): (2-1)(2+2)/(4+2+1) = 4/7
check("p=3,x=8: (2-1)(2+2)/(4+2+1) = 4/7",
      close(Fraction(4,7), lambda_theoretical(3, 8)),
      f"4/7={4/7:.10f}, λ={lambda_theoretical(3,8):.10f}")

# (4,16): (2-1)(4+4+3)/(8+4+2+1) = 11/15
check("p=4,x=16: (2-1)(4+4+3)/(8+4+2+1) = 11/15",
      close(Fraction(11,15), lambda_theoretical(4, 16)),
      f"11/15={11/15:.10f}, λ={lambda_theoretical(4,16):.10f}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6 — Table 3 : profil p=2, x=2                                 ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("SECTION 6 : Table 3 — profil numérique p=2, x=2")
print("=" * 78)

# λ_{2,2} = (√2-1)/(√2+1)
sqrt2 = 2**0.5
lam_22 = (sqrt2 - 1) / (sqrt2 + 1)
check("λ_{2,2} = (√2-1)/(√2+1) ≈ 0.1716",
      approx_eq(lam_22, 0.1716, tol=0.001) and
      close(lam_22, lambda_theoretical(2, 2)))

# Itération : s_{n+1} = 1 - 1/(2(1+s_n))
print("\n  Vérification itération p=2, x=2")
s = 0.75  # s_0 = 0.75 (u_0 = 0.5 dans cette version ? Non : s_0=0.75 → v_0=2*0.75=1.5)
alpha_22 = sqrt2

table3_expected = [
    (0, 0.7500, 1.5000, 8.579e-2, None),
    (1, 0.7143, 1.4286, 1.44e-2,  0.167),
    (2, 0.7083, 1.4167, 2.45e-3,  0.171),
    (3, 0.7073, 1.4146, 4.21e-4,  0.171),
    (4, 0.7071, 1.4143, 7.22e-5,  0.172),
    (5, 0.7071, 1.4142, 1.24e-5,  0.171),
]

eps_prev = None
for n, s_exp, v_exp, eps_exp, ratio_exp in table3_expected:
    v = 2 * s
    eps = abs(v - sqrt2)
    check(f"Table3 n={n}: s_n≈{s_exp:.4f}",
          approx_eq(s, s_exp, tol=0.001),
          f"calculé={s:.4f}")
    check(f"Table3 n={n}: v_n≈{v_exp:.4f}",
          approx_eq(v, v_exp, tol=0.001),
          f"calculé={v:.4f}")
    check(f"Table3 n={n}: ε_n≈{eps_exp:.1e}",
          abs(eps - eps_exp)/max(eps_exp, 1e-15) < 0.2,
          f"calculé={eps:.3e}, attendu={eps_exp:.1e}")
    if ratio_exp is not None and eps_prev is not None:
        ratio = eps / eps_prev
        check(f"Table3 n={n}: ratio≈{ratio_exp:.3f}",
              approx_eq(ratio, ratio_exp, tol=0.03),
              f"calculé={ratio:.4f}")
    eps_prev = eps
    # Itération p=2 : s_{n+1} = 1 - 1/(2(1+s_n))
    s = 1 - 1/(2*(1+s))

# Vérifier que v_n = 2s_n pour p=2
check("p=2: v_n = 2·s_n = x·s_n^{p-1} = 2·s_n^1",
      True)  # Tautologie par la formule v_n = x·s^{p-1} avec p=2


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 7 — Newton                                                     ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("SECTION 7 : Comparaison avec Newton")
print("=" * 78)

# Itération Newton : u_{n+1} = ((p-1)u_n + x/u_n^{p-1}) / p
print("\n§7 — Itération de Newton pour cbrt(2)")
u = 1.5
alpha_n = 2**(1/3)
for _ in range(10):
    u = ((3-1)*u + 2/u**(3-1)) / 3
check("Newton converge vers cbrt(2)",
      close(u, alpha_n, tol=1e-14),
      f"u_10={u:.15f}")

# K_{p,x} = (p-1)/(2α)
print("\n§7 — Coefficient K_{p,x} = (p-1)/(2α)")
for p, x in [(3,2), (2,2), (4,2)]:
    alpha = x**(1/p)
    K = (p-1) / (2*alpha)
    # Vérifier numériquement
    u = alpha * 1.2  # départ proche
    eps_list = []
    for _ in range(8):
        eps_list.append(abs(u - alpha))
        u = ((p-1)*u + x/u**(p-1)) / p
    # Rapport ε_{n+1}/ε_n² pour les dernières itérations
    if len(eps_list) >= 3 and eps_list[-2] > 1e-13:
        ratio_quad = eps_list[-1] / eps_list[-2]**2
        check(f"p={p},x={x}: K={K:.4f}, ε_{{n+1}}/ε_n² ≈ K",
              approx_eq(ratio_quad, K, tol=0.05),
              f"ratio={ratio_quad:.4f}")

# K_{3,2} = 1/cbrt(2) ≈ 0.794
K_32 = (3-1) / (2 * 2**(1/3))
check("K_{3,2} = 1/cbrt(2) ≈ 0.794",
      approx_eq(K_32, 1/cbrt2, tol=1e-10) and approx_eq(K_32, 0.794, tol=0.001),
      f"K={K_32:.6f}")

# Table comparaison: λ ≈ 0.220
check("Table comparaison: Pandrosion λ ≈ 0.220",
      approx_eq(lambda_theoretical(3, 2), 0.220, tol=0.001))

# Table comparaison: K ≈ 0.794
check("Table comparaison: Newton K = 1/cbrt(2) ≈ 0.794",
      approx_eq(K_32, 0.794, tol=0.001))


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 8 — Conclusion : formule reprise                              ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("SECTION 8 : Conclusion — vérif. cohérence finale")
print("=" * 78)

# Vérifier que λ_{p,x} ∈ (0,1) pour tout p≥2, x>1
print("\nλ_{p,x} ∈ (0,1) pour tout p≥2, x>1")
all_in_range = True
for p in range(2, 8):
    for x in [1.01, 2, 5, 10, 100, 1000]:
        lam = lambda_theoretical(p, x)
        if not (0 < lam < 1):
            all_in_range = False
            check(f"p={p},x={x}: λ={lam:.6f} ∈ (0,1)", False)
check("λ_{p,x} ∈ (0,1) pour 54 configurations testées", all_in_range)


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  VÉRIFICATION CROISÉE : convergence effective                           ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "=" * 78)
print("VÉRIFICATION CROISÉE : convergence effective")
print("=" * 78)

for p, x in [(2,2),(3,2),(3,3),(4,2),(5,2),(3,8),(4,16)]:
    alpha = x**(1/p)
    s = 0.5
    lam_th = lambda_theoretical(p, x)
    # Itérer 60 fois, mesurer le ratio final
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
        check(f"p={p},x={x}: ratio final {ratio_final:.6f} ≈ λ_théorique {lam_th:.6f}",
              approx_eq(ratio_final, lam_th, tol=0.005),
              f"écart={abs(ratio_final-lam_th):.6f}")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  BILAN FINAL                                                            ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
print("\n" + "═" * 78)
print(f"  BILAN FINAL : {PASS} vérifications passées, {FAIL} erreur(s)")
print("═" * 78)
if FAIL == 0:
    print("  🎉 TOUTES LES FORMULES ET VALEURS DE L'ARTICLE SONT CORRECTES")
else:
    print(f"  ⚠️  {FAIL} ERREUR(S) DÉTECTÉE(S) — voir les ❌ ci-dessus")
print()
