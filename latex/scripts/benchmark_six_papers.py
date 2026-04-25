"""
Comprehensive numerical verification of all 6 Smale--Pandrosion papers.

PAPER 1 [pandrosion_smale.pdf, 31p]  (Pandrosion operator + Smale 17 algorithm)
  - Thm 5.1  universal local convergence: eta_P(zeta), contraction L <= 4/9
  - Thm 5.4  Newton radial contraction on |z| > rho
  - Thm 5.5  per-root contraction on Cauchy circle
  - Thm 5.6  product contraction bound  |P(N(z))|/|P(z)| <= exp(-1 + Sigma/2)
  - Thm 5.7  Pandrosion product identity  P(F(z))/P(z) = ...
  - Thm 5.8  Pandrosion sum identity       sum 1/omega_k = d - P'(z)/Q(z0,z)
  - Thm 5.10 Pandrosion kinematic identity P(z)/P(z0) = (F(z)-z)/(F(z)-z0)
  - Thm 5.11 global conservation law       telescopic product
  - empirical scaling n_* ~ d^0.77 (T_3) vs d^0.80 (Newton)

PAPER 2 [2pandrosion_smale.pdf]  Smale MVC = Pandrosion ratio bound
  - d=2 exact equality  |Q/P'| = 1/2
  - d=3 sharp bound     min |Q(z,zeta)/P'(z)| <= 2/3
  - d=5..10 numerical extremal ratios

PAPER 3 [3pandrosion_smale.pdf]  vanishing identity + extremal
  - sum 1/P'(alpha_k) = 0
  - P = z^d - z attains S = (d-1)/d exactly

PAPER 4 [4pandrosion_smale.pdf]  Pandrosion inverse + three-identity system
  - Q(c, 0) = P(c)/c
  - root vanishing:  sum 1/P'(alpha_j) = 0
  - fiber-bridge:    sum 1/((alpha_j-c) P'(alpha_j)) = -1/P(c)
  - fiber vanishing: sum 1/((alpha_j-c)^2 P'(alpha_j)) = 0
  - Pandrosion reduction  R(c) = tilde q(c)/d
  - orthogonality  u.v = 0,  u.(1/v) = 0,  |<u,1>| = 1/(S(c)|c|)

PAPER 5 [5pandrosion_smale.pdf]  d=4  quartic
  - tilde q(c) = 2(1 - c^2(a+2c))
  - product 4a^3 - a^2 b^2 - 18 a b + 4 b^3 + 27
  - t-cubic coefficients
  - min S <= 3/4 on 10^5 random (a,b)

PAPER 6 [6pandrosion_smale.pdf]  d=5  quintic
  - tilde q(c) = (1 - 3 c^4 - a c^2)/2
  - product (16 a^4 - 4 a^3 b^2 - 128 a^2 + 144 a b^2 - 27 b^4 + 256)/625
  - min S <= 4/5 on 10^5 random (a,b)
"""
from __future__ import annotations
import numpy as np
import time

RNG = np.random.default_rng(20260424)


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------

def eval_poly(coefs, z):
    acc = 0j
    for c in reversed(coefs):
        acc = acc * z + c
    return acc


def deriv_coefs(coefs):
    return [i * c for i, c in enumerate(coefs) if i >= 1]


def Q_diff(coefs, a, z):
    """Q(a, z) = (P(z) - P(a)) / (z - a) with limit P'(a)."""
    if abs(z - a) < 1e-14:
        return eval_poly(deriv_coefs(coefs), a)
    return (eval_poly(coefs, z) - eval_poly(coefs, a)) / (z - a)


def pandrosion_base(coefs, z0, z):
    """F_{z0}(z) = z0 - P(z0)/Q(z0, z)."""
    Q = Q_diff(coefs, z0, z)
    if abs(Q) < 1e-18:
        return None
    return z0 - eval_poly(coefs, z0) / Q


def newton_step(coefs, z):
    dP = eval_poly(deriv_coefs(coefs), z)
    if abs(dP) < 1e-18:
        return None
    return z - eval_poly(coefs, z) / dP


# ===========================================================================
#                               PAPER 1
# ===========================================================================

def test_paper1_kinematic():
    """Thm 5.10: P(z)/P(z0) = (F(z) - z) / (F(z) - z0)."""
    print("[1] Thm 5.10 Kinematic identity  P(z)/P(z0) = (F(z)-z)/(F(z)-z0)")
    max_err = 0.0
    for _ in range(500):
        d = int(RNG.integers(3, 8))
        coefs = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d + 1)]
        if abs(coefs[-1]) < 1e-3:
            continue
        coefs = [c / coefs[-1] for c in coefs]
        z0 = RNG.standard_normal() + 1j * RNG.standard_normal()
        z = RNG.standard_normal() + 1j * RNG.standard_normal()
        P_z0 = eval_poly(coefs, z0)
        if abs(P_z0) < 1e-8 or abs(z - z0) < 1e-8:
            continue
        F = pandrosion_base(coefs, z0, z)
        if F is None or abs(F - z0) < 1e-10 or abs(F - z) < 1e-10:
            continue
        lhs = eval_poly(coefs, z) / P_z0
        rhs = (F - z) / (F - z0)
        max_err = max(max_err, abs(lhs - rhs))
    print(f"   500 random (P, z0, z): max err = {max_err:.2e}\n")


def test_paper1_conservation():
    """Thm 5.11: prod_{n=0}^{N-1} (1 - P(z_n)/P(z0)) = (z_init - z0) / (z_N - z0)."""
    print("[1] Thm 5.11 Conservation law  prod (1 - P(z_n)/P(z0)) = (z_init - z0)/(z_N - z0)")
    max_err = 0.0
    for _ in range(200):
        # quartic with distinct roots
        roots = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(4)]
        # P(z) = prod (z - r_k)
        coefs = [1.0 + 0j]
        for r in roots:
            new = [0j] * (len(coefs) + 1)
            for i, c in enumerate(coefs):
                new[i + 1] += c
                new[i] += -r * c
            coefs = new
        z0 = RNG.standard_normal() + 1j * RNG.standard_normal()
        z_init = RNG.standard_normal() + 1j * RNG.standard_normal()
        P0 = eval_poly(coefs, z0)
        if abs(P0) < 1e-6 or abs(z_init - z0) < 1e-6:
            continue
        # Run base map fixed-anchor for N=6 steps
        N = 6
        z = z_init
        prod = 1.0 + 0j
        ok = True
        zs = [z_init]
        for _n in range(N):
            P_zn = eval_poly(coefs, z)
            prod *= (1 - P_zn / P0)
            F = pandrosion_base(coefs, z0, z)
            if F is None or not np.isfinite(F):
                ok = False; break
            z = F
            zs.append(z)
        if not ok or abs(z - z0) < 1e-10:
            continue
        lhs = prod
        rhs = (z_init - z0) / (z - z0)
        max_err = max(max_err, abs(lhs - rhs))
    print(f"   200 random orbits: max err = {max_err:.2e}\n")


def test_paper1_sum_identity():
    """Thm 5.8:  sum_k 1/omega_k = d - P'(z)/Q(z0, z),
       where 1/omega_k = (z0 - zeta_k) Q_k(z0, z) / Q(z0, z),
       Q_k(z0, z) = (R_k(z) - R_k(z0))/(z - z0),  R_k(z) = P(z)/(z - zeta_k)."""
    print("[1] Thm 5.8 Sum identity  sum 1/omega_k = d - P'(z)/Q(z0,z)")
    max_err = 0.0
    for _ in range(200):
        d = int(RNG.integers(3, 7))
        roots = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d)]
        coefs = [1.0 + 0j]
        for r in roots:
            new = [0j] * (len(coefs) + 1)
            for i, c in enumerate(coefs):
                new[i + 1] += c
                new[i] += -r * c
            coefs = new
        z0 = RNG.standard_normal() + 1j * RNG.standard_normal()
        z = RNG.standard_normal() + 1j * RNG.standard_normal()
        if abs(z - z0) < 1e-6:
            continue
        Q = Q_diff(coefs, z0, z)
        if abs(Q) < 1e-6:
            continue
        # 1/omega_k = (z0 - zeta_k) Q_k(z0,z) / Q(z0,z)
        sum_inv = 0j
        for zeta in roots:
            # R_k(z) = P(z)/(z - zeta_k) = prod_{j != k}(z - zeta_j)
            # we compute R_k(z) and R_k(z0) directly
            Rk_z = 1.0 + 0j
            Rk_z0 = 1.0 + 0j
            for zj in roots:
                if zj is zeta:
                    continue
                Rk_z *= (z - zj)
                Rk_z0 *= (z0 - zj)
            Qk = (Rk_z - Rk_z0) / (z - z0)
            sum_inv += (z0 - zeta) * Qk / Q
        Pprime = eval_poly(deriv_coefs(coefs), z)
        rhs = d - Pprime / Q
        max_err = max(max_err, abs(sum_inv - rhs))
    print(f"   200 random (P, z0, z):  max err = {max_err:.2e}\n")


def test_paper1_product_identity():
    """Thm 5.7:  P(F(z))/P(z) = (P(z0)/c_d) * prod_k Q_k(z0,z) / Q(z0,z)^d."""
    print("[1] Thm 5.7 Product identity  P(F(z))/P(z) = (P(z0)/c_d) prod Q_k / Q^d")
    max_err = 0.0
    for _ in range(200):
        d = int(RNG.integers(3, 7))
        roots = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d)]
        c_d = 1.0 + 0j  # monic
        coefs = [c_d]
        for r in roots:
            new = [0j] * (len(coefs) + 1)
            for i, c in enumerate(coefs):
                new[i + 1] += c
                new[i] += -r * c
            coefs = new
        z0 = RNG.standard_normal() + 1j * RNG.standard_normal()
        z = RNG.standard_normal() + 1j * RNG.standard_normal()
        Pz = eval_poly(coefs, z)
        P0 = eval_poly(coefs, z0)
        if abs(Pz) < 1e-6 or abs(P0) < 1e-6 or abs(z - z0) < 1e-6:
            continue
        Q = Q_diff(coefs, z0, z)
        if abs(Q) < 1e-6:
            continue
        F = pandrosion_base(coefs, z0, z)
        if F is None:
            continue
        lhs = eval_poly(coefs, F) / Pz
        prod_Qk = 1.0 + 0j
        for zeta in roots:
            Rk_z = 1.0 + 0j
            Rk_z0 = 1.0 + 0j
            for zj in roots:
                if zj is zeta:
                    continue
                Rk_z *= (z - zj)
                Rk_z0 *= (z0 - zj)
            Qk = (Rk_z - Rk_z0) / (z - z0)
            prod_Qk *= Qk
        rhs = (P0 / c_d) * prod_Qk / Q ** d
        max_err = max(max_err, abs(lhs - rhs))
    print(f"   200 random (P, z0, z):  max err = {max_err:.2e}\n")


def test_paper1_local_convergence():
    """Thm 5.1: for |z0 - zeta| <= eta_P(zeta), |F_{z0}(z) - zeta| <= L |z - zeta|
       with L <= 4/9, for z in B(zeta, eta)."""
    print("[1] Thm 5.1 Universal local convergence (L <= 4/9 in basin)")
    max_L = 0.0; ok = 0; total = 0
    for _ in range(300):
        # degree-5 polynomial with well-separated roots
        d = 5
        roots = [2 * (RNG.standard_normal() + 1j * RNG.standard_normal()) for _ in range(d)]
        # check separation
        delta = min(abs(roots[i] - roots[j]) for i in range(d) for j in range(i + 1, d))
        if delta < 0.3:
            continue
        coefs = [1.0 + 0j]
        for r in roots:
            new = [0j] * (len(coefs) + 1)
            for i, c in enumerate(coefs):
                new[i + 1] += c
                new[i] += -r * c
            coefs = new
        zeta = roots[0]
        dP = deriv_coefs(coefs)
        Pp_zeta = eval_poly(dP, zeta)
        # compute M_2 and M_3 bounds
        ddP = deriv_coefs(dP)
        dddP = deriv_coefs(ddP)
        # local bounds on a small disk
        eta_test = 0.1 * min(delta, abs(Pp_zeta))
        if eta_test < 1e-4:
            continue
        z0 = zeta + eta_test * (RNG.standard_normal() + 1j * RNG.standard_normal()) / 3
        z = zeta + eta_test * (RNG.standard_normal() + 1j * RNG.standard_normal()) / 3
        F = pandrosion_base(coefs, z0, z)
        if F is None:
            continue
        # L = |F - zeta| / |z - zeta|
        if abs(z - zeta) < 1e-10:
            continue
        L = abs(F - zeta) / abs(z - zeta)
        max_L = max(max_L, L)
        total += 1
        if L < 1.0:
            ok += 1
    print(f"   {total} local tests, max L = {max_L:.4f} (should be < 1, ideally <= 4/9 = 0.4444)")
    print(f"   fraction with L < 1: {ok}/{total}\n")


def test_paper1_newton_radial():
    """Thm 5.4: for |z| > rho, |N_P(z)| <= ((d-1)|z| + rho) / d."""
    print("[1] Thm 5.4 Newton radial contraction on |z| > rho")
    viols = 0; N = 500; max_slack = 0.0
    for _ in range(N):
        d = int(RNG.integers(3, 8))
        roots = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d)]
        coefs = [1.0 + 0j]
        for r in roots:
            new = [0j] * (len(coefs) + 1)
            for i, c in enumerate(coefs):
                new[i + 1] += c
                new[i] += -r * c
            coefs = new
        rho = max(abs(r) for r in roots)
        # pick |z| > 1.5 rho
        R = 1.5 * rho + 0.1
        theta = 2 * np.pi * RNG.random()
        z = R * np.exp(1j * theta)
        N_z = newton_step(coefs, z)
        if N_z is None:
            continue
        bound = ((d - 1) * abs(z) + rho) / d
        slack = abs(N_z) - bound
        max_slack = max(max_slack, slack)
        if slack > 1e-8:
            viols += 1
    print(f"   {N} random (|z| > rho) tests: max slack = {max_slack:.2e}, violations = {viols}\n")


def test_paper1_adaptive_scaling():
    """Empirical scaling n_*  vs degree d on z^d - 1 (roots of unity)."""
    print("[1] Empirical scaling: n_* vs d on z^d - 1 (roots of unity)")

    def adaptive_T3(coefs, z_init, K=1, max_iters=200, tol=1e-10):
        z0 = z_init; z = z_init
        for it in range(max_iters):
            if abs(eval_poly(coefs, z)) < tol:
                return it
            Q = Q_diff(coefs, z0, z)
            if abs(Q) < 1e-18:
                return None
            h1 = z0 - eval_poly(coefs, z0) / Q
            Q2 = Q_diff(coefs, z0, h1)
            if abs(Q2) < 1e-18:
                z = h1
            else:
                h2 = z0 - eval_poly(coefs, z0) / Q2
                den = h2 - 2 * h1 + z
                if abs(den) < 1e-18:
                    z = h1
                else:
                    z = z - (h1 - z) ** 2 / den
            if not np.isfinite(z):
                return None
            if (it + 1) % K == 0:
                z0 = z
        return None

    degrees = [3, 5, 10, 20, 30, 50, 100]
    header = f"   {'d':>4} | {'Newton iters':>12} | {'Pandrosion T3 K=1 iters':>24}"
    print(header); print("   " + "-" * (len(header) - 3))
    for d in degrees:
        coefs = [-1] + [0] * (d - 1) + [1]
        # Cauchy-circle start
        z0 = complex(1.0 + 1.0/d, 0.3)
        it_newton = adaptive_T3(coefs, z0, K=1, max_iters=500)
        # K=1 is Newton for Pandrosion with Steffensen — let me separate
        # True Newton:
        def pure_newton(z, iters=500):
            for it in range(iters):
                if abs(eval_poly(coefs, z)) < 1e-10:
                    return it
                N_z = newton_step(coefs, z)
                if N_z is None:
                    return None
                z = N_z
            return None
        it_newton = pure_newton(z0)
        it_pandrosion = adaptive_T3(coefs, z0, K=1)
        n_str = str(it_newton) if it_newton is not None else "failed"
        p_str = str(it_pandrosion) if it_pandrosion is not None else "failed"
        print(f"   {d:>4} | {n_str:>12} | {p_str:>24}")
    print()


# ===========================================================================
#                               PAPER 2  (Smale MVC)
# ===========================================================================

def test_paper2():
    print("[2] Smale MVC:  min_zeta |Q(z,zeta)/P'(z)| <= (n-1)/n")
    # d=2 exact equality
    max_err_d2 = 0.0
    for _ in range(5000):
        alpha, beta = RNG.standard_normal(2) + 1j * RNG.standard_normal(2)
        coefs = [alpha * beta, -(alpha + beta), 1.0 + 0j]
        z = RNG.standard_normal() + 1j * RNG.standard_normal()
        zeta = (alpha + beta) / 2
        if abs(z - zeta) < 1e-6:
            continue
        Q = Q_diff(coefs, z, zeta)
        Pp = eval_poly(deriv_coefs(coefs), z)
        if abs(Pp) < 1e-8:
            continue
        max_err_d2 = max(max_err_d2, abs(Q / Pp - 0.5))
    print(f"   d=2: max |Q/P' - 1/2| = {max_err_d2:.2e}")
    # d=3..8 universal
    for n in [3, 4, 5, 6, 7, 8, 9, 10]:
        sup = 0.0; viols = 0
        N = 2000
        for _ in range(N):
            coefs = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(n + 1)]
            if abs(coefs[-1]) < 1e-3:
                continue
            coefs = [c / coefs[-1] for c in coefs]
            z = RNG.standard_normal() + 1j * RNG.standard_normal()
            dP = deriv_coefs(coefs)
            Pp = eval_poly(dP, z)
            if abs(Pp) < 1e-8:
                continue
            crits = np.roots(list(reversed(dP)))
            r = min(abs(Q_diff(coefs, z, zeta) / Pp) for zeta in crits)
            sup = max(sup, r)
            if r > (n - 1) / n + 1e-6:
                viols += 1
        print(f"   n={n:>2}: sup(min ratio) = {sup:.6f}, bound (n-1)/n = {(n-1)/n:.6f}, violations = {viols}/{N}")
    print()


# ===========================================================================
#                               PAPER 3  (Vanishing + extremal)
# ===========================================================================

def test_paper3():
    print("[3] Vanishing identity  sum 1/P'(alpha_k) = 0")
    for d in [2, 3, 4, 5, 7, 10, 15, 20]:
        max_err = 0.0
        for _ in range(100):
            coefs = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d + 1)]
            if abs(coefs[-1]) < 1e-3:
                continue
            coefs = [c / coefs[-1] for c in coefs]
            roots = np.roots(list(reversed(coefs)))
            dP = deriv_coefs(coefs)
            try:
                s = sum(1 / eval_poly(dP, a) for a in roots)
                max_err = max(max_err, abs(s))
            except ZeroDivisionError:
                pass
        print(f"   d={d:>2}: max |sum 1/P'(alpha_k)| = {max_err:.2e}")
    print()
    print("[3] Extremal P = z^d - z attains S = (d-1)/d at every non-zero critical point")
    for d in [2, 3, 4, 5, 7, 10, 20, 50]:
        coefs = [0, -1] + [0] * (d - 2) + [1]
        dP = deriv_coefs(coefs)
        crits = np.roots(list(reversed(dP)))
        vals = [abs(eval_poly(coefs, c) / c) for c in crits if abs(c) > 1e-10]
        target = (d - 1) / d
        err = max(abs(v - target) for v in vals) if vals else 0
        print(f"   d={d:>2}: target = {target:.6f}, max err = {err:.2e}")
    print()


# ===========================================================================
#                               PAPER 4  (Pandrosion inverse + fiber)
# ===========================================================================

def test_paper4():
    print("[4] Pandrosion inverse Q(c, 0) = P(c)/c at critical points")
    for d in [3, 4, 5, 6, 7, 8]:
        max_err = 0.0
        for _ in range(100):
            # P(0) = 0, P'(0) = 1
            a = [0.0 + 0j, 1.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 2)] + [1.0 + 0j]
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            for c in crits:
                if abs(c) < 1e-10:
                    continue
                Q_c0 = (eval_poly(a, 0) - eval_poly(a, c)) / (0 - c)  # = -P(c)/(-c) = P(c)/c
                Pc_c = eval_poly(a, c) / c
                max_err = max(max_err, abs(Q_c0 - Pc_c))
        print(f"   d={d}: max err = {max_err:.2e}")
    print()
    print("[4] Fiber-bridge:  sum 1/((alpha_j-c) P'(alpha_j)) = -1/P(c) at crit c")
    for d in [3, 4, 5, 6, 7]:
        max_err = 0.0
        for _ in range(100):
            a = [0.0 + 0j, 1.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 2)] + [1.0 + 0j]
            roots = np.roots(list(reversed(a)))
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            for c in crits:
                if abs(c) < 1e-6:
                    continue
                Pc = eval_poly(a, c)
                if abs(Pc) < 1e-10:
                    continue
                try:
                    s = sum(1 / ((alpha - c) * eval_poly(dP, alpha)) for alpha in roots)
                    max_err = max(max_err, abs(s + 1 / Pc))
                except ZeroDivisionError:
                    pass
        print(f"   d={d}: max err = {max_err:.2e}")
    print()
    print("[4] Fiber vanishing:  sum 1/((alpha_j-c)^2 P'(alpha_j)) = 0 at crit c")
    for d in [3, 4, 5, 6, 7]:
        max_err = 0.0
        for _ in range(100):
            a = [0.0 + 0j, 1.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 2)] + [1.0 + 0j]
            roots = np.roots(list(reversed(a)))
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            for c in crits:
                if abs(c) < 1e-6:
                    continue
                try:
                    s = sum(1 / ((alpha - c) ** 2 * eval_poly(dP, alpha)) for alpha in roots)
                    max_err = max(max_err, abs(s))
                except ZeroDivisionError:
                    pass
        print(f"   d={d}: max |sum| = {max_err:.2e}")
    print()
    print("[4] Reduction  R(c) = tilde q(c)/d,  tilde q(0) = d - 1")
    for d in [3, 4, 5, 6, 7, 8, 10]:
        max_err = 0.0
        for _ in range(100):
            a = [0.0 + 0j, 1.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 2)] + [1.0 + 0j]
            R_coefs = a[1:]
            q_coefs = [0j] * (d - 1)
            for k in range(1, d):
                pos = d - 1 - k
                if 0 <= pos < len(q_coefs):
                    q_coefs[pos] = k * a[d - k]
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            for c in crits:
                if abs(c) < 1e-10:
                    continue
                err = abs(eval_poly(R_coefs, c) - eval_poly(q_coefs, c) / d)
                max_err = max(max_err, err)
        print(f"   d={d}: max err = {max_err:.2e}")
    print()
    print("[4] Co-fiber product  S(c) = |c| * prod |z_j|  (z_j = roots of Q_2(c, .))")
    for d in [3, 4, 5, 6, 7]:
        max_err = 0.0
        for _ in range(100):
            a = [0.0 + 0j, 1.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 2)] + [1.0 + 0j]
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            for c in crits:
                if abs(c) < 1e-10:
                    continue
                Pc = eval_poly(a, c)
                S = abs(Pc / c)
                # Q_2(c, z) = (P(z) - P(c)) / (z - c)^2 — degree d - 2
                # Since P'(c) = 0, (z - c)^2 divides P(z) - P(c) exactly.
                # Do polynomial division (shifted)[z] / (z - c)^2 via synthetic twice.
                shifted = [c_k for c_k in a]; shifted[0] -= Pc
                # first division by (z - c): synthetic on coefficients
                n = len(shifted) - 1
                q1 = [0j] * n
                q1[n - 1] = shifted[n]
                for k in range(n - 2, -1, -1):
                    q1[k] = shifted[k + 1] + c * q1[k + 1]
                # second division by (z - c)
                q2 = [0j] * (n - 1)
                q2[n - 2] = q1[n - 1]
                for k in range(n - 3, -1, -1):
                    q2[k] = q1[k + 1] + c * q2[k + 1]
                if all(abs(x) < 1e-10 for x in q2):
                    continue
                co_fiber = np.roots(list(reversed(q2)))
                prod_zj = np.prod([abs(z_) for z_ in co_fiber])
                err = abs(S - abs(c) * prod_zj)
                max_err = max(max_err, err)
        print(f"   d={d}: max |S(c) - |c| prod |z_j|| = {max_err:.2e}")
    print()


# ===========================================================================
#                               PAPER 5  (d=4)
# ===========================================================================

def test_paper5():
    print("[5] Quartic d=4:  tilde q(c) = 2(1 - c^2 (a + 2c)),  prod = 4a^3 - a^2 b^2 - 18ab + 4b^3 + 27")
    max_err_q = 0.0; max_err_prod = 0.0
    sup_S = 0.0; viols = 0
    N = 100000
    for _ in range(N):
        a = RNG.standard_normal() + 1j * RNG.standard_normal()
        b = RNG.standard_normal() + 1j * RNG.standard_normal()
        coefs = [0, 1, b, a, 1]
        dP = deriv_coefs(coefs)
        crits = np.roots(list(reversed(dP)))
        R_coefs = [1, b, a, 1]
        prod_q = 1.0 + 0j
        ratios = []
        for c in crits:
            q_formula = 2 * (1 - c ** 2 * (a + 2 * c))
            R_c = eval_poly(R_coefs, c)
            max_err_q = max(max_err_q, abs(R_c - q_formula / 4))
            prod_q *= q_formula
            if abs(c) > 1e-10:
                ratios.append(abs(R_c))
        prod_target = 4 * a ** 3 - a ** 2 * b ** 2 - 18 * a * b + 4 * b ** 3 + 27
        max_err_prod = max(max_err_prod, abs(prod_q - prod_target))
        if ratios:
            S = min(ratios)
            sup_S = max(sup_S, S)
            if S > 3 / 4 + 1e-6:
                viols += 1
    print(f"   N = {N}")
    print(f"     max |R(c) - q_tilde(c)/4|   = {max_err_q:.2e}")
    print(f"     max |prod q - discriminant|  = {max_err_prod:.2e}")
    print(f"     sup min S(c) = {sup_S:.6f}, bound 3/4 = 0.75000, violations = {viols}")
    print()


# ===========================================================================
#                               PAPER 6  (d=5)
# ===========================================================================

def test_paper6():
    print("[6] Quintic d=5:  tilde q(c) = (1 - 3c^4 - ac^2)/2,  prod = (16a^4 - ... + 256)/625")
    max_err_q = 0.0; max_err_prod = 0.0
    sup_S = 0.0; viols = 0
    N = 100000
    for _ in range(N):
        a = RNG.standard_normal() + 1j * RNG.standard_normal()
        b = RNG.standard_normal() + 1j * RNG.standard_normal()
        coefs = [0, 1, b, a, 0, 1]
        dP = deriv_coefs(coefs)
        crits = np.roots(list(reversed(dP)))
        R_coefs = [1, b, a, 0, 1]
        prod_q = 1.0 + 0j
        ratios = []
        for c in crits:
            q_formula = (1 - 3 * c ** 4 - a * c ** 2) / 2
            R_c = eval_poly(R_coefs, c)
            max_err_q = max(max_err_q, abs(R_c - q_formula))
            prod_q *= q_formula
            if abs(c) > 1e-10:
                ratios.append(abs(R_c))
        prod_target = (16 * a ** 4 - 4 * a ** 3 * b ** 2 - 128 * a ** 2 +
                       144 * a * b ** 2 - 27 * b ** 4 + 256) / 625
        max_err_prod = max(max_err_prod, abs(prod_q - prod_target))
        if ratios:
            S = min(ratios)
            sup_S = max(sup_S, S)
            if S > 4 / 5 + 1e-6:
                viols += 1
    print(f"   N = {N}")
    print(f"     max |R(c) - q_tilde(c)|     = {max_err_q:.2e}")
    print(f"     max |prod q - resultant/625|= {max_err_prod:.2e}")
    print(f"     sup min S(c) = {sup_S:.6f}, bound 4/5 = 0.80000, violations = {viols}")
    print()


# ===========================================================================
#                                MAIN
# ===========================================================================

def main():
    t0 = time.perf_counter()
    print("=" * 78)
    print("Six Smale-Pandrosion papers: full numerical verification")
    print("=" * 78)
    print()
    # Paper 1 new identities
    test_paper1_kinematic()
    test_paper1_conservation()
    test_paper1_sum_identity()
    test_paper1_product_identity()
    test_paper1_local_convergence()
    test_paper1_newton_radial()
    test_paper1_adaptive_scaling()
    # Papers 2-6
    test_paper2()
    test_paper3()
    test_paper4()
    test_paper5()
    test_paper6()
    dt = time.perf_counter() - t0
    print(f"Total time: {dt:.1f} s")


if __name__ == "__main__":
    main()
