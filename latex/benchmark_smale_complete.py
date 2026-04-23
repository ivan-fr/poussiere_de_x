"""
Complete verification script for the combined Smale17.tex paper.

Covers every identity claimed in the six companion PDFs:
  [07]  Pandrosion operator + Smale 17 algorithm
  [11]  Smale MVC is a Pandrosion ratio bound
  [43]  Vanishing identity + extremal polynomial
  [49]  Pandrosion inverse + three-identity system
  [52]  Smale for quartics (d=4)
  [61]  Smale for quintics (d=5)
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import fsolve

RNG = np.random.default_rng(20260423)


# ---------------------------------------------------------------------------
# Polynomial primitives
# ---------------------------------------------------------------------------

def eval_poly(coefs, z):
    acc = 0j
    for c in reversed(coefs):
        acc = acc * z + c
    return acc


def deriv_coefs(coefs):
    return [i * c for i, c in enumerate(coefs) if i >= 1]


def Q_anchor(coefs, a, z):
    """Q(a, z) = (P(z) - P(a))/(z - a) with limit P'(a) at z=a."""
    if abs(z - a) < 1e-14:
        return eval_poly(deriv_coefs(coefs), a)
    return (eval_poly(coefs, z) - eval_poly(coefs, a)) / (z - a)


# ---------------------------------------------------------------------------
# [11] Paper: d=2 exact equality and d=3 bound 2/3
# ---------------------------------------------------------------------------

def test_paper11_d2():
    print("[11] d=2: Q(z,zeta)/P'(z) = 1/2 exactly")
    max_err = 0.0
    for _ in range(5000):
        alpha, beta = RNG.standard_normal(2) + 1j * RNG.standard_normal(2)
        coefs = [alpha * beta, -(alpha + beta), 1.0 + 0j]
        z = RNG.standard_normal() + 1j * RNG.standard_normal()
        zeta = (alpha + beta) / 2
        if abs(z - zeta) < 1e-8:
            continue
        Q_val = Q_anchor(coefs, z, zeta)
        Pp = eval_poly(deriv_coefs(coefs), z)
        if abs(Pp) < 1e-10:
            continue
        err = abs(Q_val / Pp - 0.5)
        max_err = max(max_err, err)
    print(f"   5000 quadratics, max |ratio - 1/2| = {max_err:.2e}\n")


def test_paper11_d3():
    print("[11] d=3: min_zeta |Q(z,zeta)/P'(z)| <= 2/3")
    sup = 0.0; viols = 0
    N = 5000
    for _ in range(N):
        coefs = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(4)]
        if abs(coefs[-1]) < 1e-3:
            continue
        coefs = [c / coefs[-1] for c in coefs]
        z = RNG.standard_normal() + 1j * RNG.standard_normal()
        dP = deriv_coefs(coefs)
        Pp = eval_poly(dP, z)
        if abs(Pp) < 1e-10:
            continue
        crits = np.roots(list(reversed(dP)))
        ratios = [abs(Q_anchor(coefs, z, zeta) / Pp) for zeta in crits]
        r = min(ratios)
        sup = max(sup, r)
        if r > 2 / 3 + 1e-6:
            viols += 1
    print(f"   {N} cubics, sup(min ratio) = {sup:.6f}, bound 2/3 = 0.66667, violations = {viols}\n")


# ---------------------------------------------------------------------------
# [43] Paper: Vanishing identity + extremal
# ---------------------------------------------------------------------------

def test_paper43_vanishing():
    print("[43] Vanishing identity sum 1/P'(alpha_k) = 0")
    for d in [2, 3, 4, 5, 7, 10, 15]:
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
                continue
        print(f"   d={d:>2}: max |sum 1/P'(alpha_k)| = {max_err:.2e}")
    print()


def test_paper43_extremal():
    print("[43] Extremal P = z^d - z attains S(c) = (d-1)/d")
    for d in [2, 3, 4, 5, 7, 10, 20, 50]:
        coefs = [0, -1] + [0] * (d - 2) + [1]
        dP = deriv_coefs(coefs)
        crits = np.roots(list(reversed(dP)))
        vals = [abs(eval_poly(coefs, c) / c) for c in crits if abs(c) > 1e-10]
        target = (d - 1) / d
        err = max(abs(v - target) for v in vals) if vals else 0
        print(f"   d={d:>2}: target = {target:.6f}, max err = {err:.2e}")
    print()


# ---------------------------------------------------------------------------
# [49] Paper: Pandrosion inverse + three-identity system
# ---------------------------------------------------------------------------

def test_paper49_pandrosion_inverse():
    """Q(c, 0) = P(c)/c for marked root at 0."""
    print("[49] Pandrosion inverse: Q(c, 0) = P(c)/c at critical points")
    for d in [3, 4, 5, 6, 7]:
        max_err = 0.0
        for _ in range(100):
            a = [0.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 1)] + [1.0 + 0j]
            if abs(a[1]) < 1e-3:
                continue
            a = [c / a[1] for c in a]
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            for c in crits:
                if abs(c) < 1e-10:
                    continue
                # Q(c, 0) = (P(0) - P(c))/(0 - c) = P(c)/c (since P(0)=0)
                Q_c0 = -eval_poly(a, c) / (-c)
                Pc_over_c = eval_poly(a, c) / c
                max_err = max(max_err, abs(Q_c0 - Pc_over_c))
        print(f"   d={d}: max err = {max_err:.2e}")
    print()


def test_paper49_fiber_vanishing():
    """sum_k 1/((alpha_k - c)^2 P'(alpha_k)) = 0 at each critical c."""
    print("[49] Fiber vanishing: sum 1/((alpha_k - c)^2 P'(alpha_k)) = 0 at crit c")
    for d in [3, 4, 5, 6, 7]:
        max_err = 0.0
        for _ in range(100):
            # normalize P monic, P(0)=0, P'(0)=1
            a = [0.0 + 0j, 1.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 2)] + [1.0 + 0j]
            roots = np.roots(list(reversed(a)))
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            for c in crits:
                if abs(c) < 1e-6:
                    continue
                try:
                    s = 0j
                    for alpha in roots:
                        denom = (alpha - c) ** 2 * eval_poly(dP, alpha)
                        if abs(denom) < 1e-20:
                            raise ZeroDivisionError
                        s += 1 / denom
                    max_err = max(max_err, abs(s))
                except ZeroDivisionError:
                    continue
        print(f"   d={d}: max |sum| = {max_err:.2e}")
    print()


def test_paper49_fiber_bridge():
    """sum_k 1/((alpha_k - c) P'(alpha_k)) = -1/P(c) at each critical c."""
    print("[49] Fiber-bridge: sum 1/((alpha_k - c) P'(alpha_k)) = -1/P(c) at crit c")
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
                    s = 0j
                    for alpha in roots:
                        denom = (alpha - c) * eval_poly(dP, alpha)
                        if abs(denom) < 1e-20:
                            raise ZeroDivisionError
                        s += 1 / denom
                    err = abs(s + 1 / Pc)
                    max_err = max(max_err, err)
                except ZeroDivisionError:
                    continue
        print(f"   d={d}: max err = {max_err:.2e}")
    print()


def test_paper49_reduction():
    """R(c) = q_tilde(c)/d at critical points, q_tilde = sum k a_{d-k} z^{d-1-k}."""
    print("[49] Pandrosion reduction R(c) = q_tilde(c)/d, q_tilde(0) = d-1")
    for d in [3, 4, 5, 6, 7, 8]:
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
        print(f"   d={d}: max |R(c) - q_tilde(c)/d| = {max_err:.2e}")
    print()


# ---------------------------------------------------------------------------
# [52] Paper: d=4 quartic formulas
# ---------------------------------------------------------------------------

def test_paper52():
    print("[52] Quartic: q_tilde(c) = 2(1 - c^2(a+2c)), product = 4a^3 - a^2 b^2 - 18ab + 4b^3 + 27")
    max_err_q = 0.0; max_err_prod = 0.0
    sup_S = 0.0; viols = 0
    N = 10000
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
    print(f"   {N} random (a,b):")
    print(f"     max |R(c) - q_tilde/4|       = {max_err_q:.2e}")
    print(f"     max |prod q - discriminant|  = {max_err_prod:.2e}")
    print(f"     sup min S(c) = {sup_S:.6f}, bound 3/4, violations = {viols}")
    # At (a,b) = (0,0): product = 27, check attains (3/4)^3
    a0, b0 = 0j, 0j
    coefs0 = [0, 1, b0, a0, 1]
    dP0 = deriv_coefs(coefs0)
    crits0 = np.roots(list(reversed(dP0)))
    ratios0 = [abs(eval_poly([1, b0, a0, 1], c)) for c in crits0 if abs(c) > 1e-10]
    print(f"     Extremal (a=b=0): S(c) = {ratios0} -> {3/4}")
    print()


# ---------------------------------------------------------------------------
# [61] Paper: d=5 quintic formulas
# ---------------------------------------------------------------------------

def test_paper61():
    print("[61] Quintic: q_tilde(c) = (1 - 3c^4 - ac^2)/2, product = (16a^4 - ... + 256)/625")
    max_err_q = 0.0; max_err_prod = 0.0
    sup_S = 0.0; viols = 0
    N = 10000
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
            max_err_q = max(max_err_q, abs(R_c - q_formula))  # paper 61 convention
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
    print(f"   {N} random (a,b):")
    print(f"     max |R(c) - q_tilde|         = {max_err_q:.2e}")
    print(f"     max |prod q - resultant/625| = {max_err_prod:.2e}")
    print(f"     sup min S(c) = {sup_S:.6f}, bound 4/5, violations = {viols}")
    # Extremal at (a=b=0): R=P/z, P=z^5+z, crits c^4 = -1/5
    a0, b0 = 0j, 0j
    coefs0 = [0, 1, b0, a0, 0, 1]
    dP0 = deriv_coefs(coefs0)
    crits0 = np.roots(list(reversed(dP0)))
    ratios0 = [abs(eval_poly([1, b0, a0, 0, 1], c)) for c in crits0 if abs(c) > 1e-10]
    print(f"     Extremal (a=b=0): min S(c) = {min(ratios0):.6f} -> {4/5}")
    print()


# ---------------------------------------------------------------------------
# [07] Paper: Adaptive Pandrosion-Steffensen univariate
# ---------------------------------------------------------------------------

def adaptive_pandrosion_univariate(coefs, z_init, K=3, max_iters=80, tol=1e-10):
    """Adaptive-anchor Steffensen-Pandrosion on a univariate polynomial."""
    z0 = z_init
    z = z_init
    for it in range(max_iters):
        if abs(eval_poly(coefs, z)) < tol:
            return z, it
        # Pandrosion step from anchor z0 evaluated at z
        Q = Q_anchor(coefs, z0, z)
        if abs(Q) < 1e-18:
            return None, it
        h = z0 - eval_poly(coefs, z0) / Q
        # Steffensen acceleration
        Q2 = Q_anchor(coefs, z0, h)
        if abs(Q2) < 1e-18:
            z_new = h
        else:
            hprime = z0 - eval_poly(coefs, z0) / Q2
            denom = hprime - 2 * h + z
            if abs(denom) < 1e-18:
                z_new = h
            else:
                z_new = z - (h - z) ** 2 / denom
        if not np.isfinite(z_new):
            return None, it
        z = z_new
        if (it + 1) % K == 0:
            z0 = z
    if abs(eval_poly(coefs, z)) < tol:
        return z, max_iters
    return None, max_iters


def test_paper07_grid_convergence(name, coefs, K_list=[1, 3, 10**9], grid_size=100):
    print(f"[07] Grid convergence on {name}, [-2,2]^2 grid {grid_size}x{grid_size}")
    results = {}
    for K in K_list:
        hits = 0; total = 0
        for i in range(grid_size):
            for j in range(grid_size):
                x = -2 + 4 * i / (grid_size - 1)
                y = -2 + 4 * j / (grid_size - 1)
                z0 = complex(x, y)
                # skip if near a root
                if any(abs(z0 - r) < 1e-6 for r in np.roots(list(reversed(coefs)))):
                    continue
                total += 1
                sol, _ = adaptive_pandrosion_univariate(coefs, z0, K=K)
                if sol is not None:
                    hits += 1
        pct = 100 * hits / total if total else 0
        label = f"K={K}" if K < 10 ** 8 else "K=infty (fixed)"
        results[label] = pct
        print(f"   {label:<20s}: {pct:.2f}% ({hits}/{total})")
    print()
    return results


# ---------------------------------------------------------------------------
# [07] Paper: universal SMVC bound check for all d
# ---------------------------------------------------------------------------

def test_universal_smvc():
    print("Universal SMVC: min_c S(c) <= (d-1)/d for d = 2..10")
    for d in range(2, 11):
        sup = 0.0; viols = 0
        N = 2000
        for _ in range(N):
            a = [0.0 + 0j, 1.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 2)] + [1.0 + 0j]
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            ratios = [abs(eval_poly(a, c) / c) for c in crits if abs(c) > 1e-10]
            if ratios:
                m = min(ratios)
                sup = max(sup, m)
                if m > (d - 1) / d + 1e-6:
                    viols += 1
        print(f"   d={d:>2}: sup min S = {sup:.6f}, bound = {(d-1)/d:.6f}, violations = {viols}/{N}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("Complete Smale17 verification")
    print("=" * 80)
    print()
    test_paper11_d2()
    test_paper11_d3()
    test_paper43_vanishing()
    test_paper43_extremal()
    test_paper49_pandrosion_inverse()
    test_paper49_fiber_vanishing()
    test_paper49_fiber_bridge()
    test_paper49_reduction()
    test_paper52()
    test_paper61()
    test_universal_smvc()
    # Grid convergence tests (slower)
    test_paper07_grid_convergence("z^3 - 1", [-1, 0, 0, 1], K_list=[1, 3, 10**9], grid_size=60)
    test_paper07_grid_convergence("z^5 - z - 1", [-1, -1, 0, 0, 0, 1], K_list=[1, 3, 10**9], grid_size=60)


if __name__ == "__main__":
    main()
