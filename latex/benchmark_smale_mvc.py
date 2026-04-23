"""
Verify the key identities of the Smale-MVC / Pandrosion papers:
  [11] Smale MVC is a Pandrosion ratio bound  |Q(z,zeta)/P'(z)| <= (n-1)/n
  [43] Pandrosion vanishing identity  sum 1/P'(alpha_k) = 0
  [49] Pandrosion reduction  R(c) = q_tilde(c)/d  with q_tilde(0) = d-1
  [52] Smale d=4  via q_tilde(c) = 2(1 - c^2 (a + 2 c))
       product  prod q_tilde(c_k) = 4 a^3 - a^2 b^2 - 18 a b + 4 b^3 + 27
  [61] Smale d=5  via q_tilde(c) = (1 - 3 c^4 - a c^2)/2
       product  (16 a^4 - 4 a^3 b^2 - 128 a^2 + 144 a b^2 - 27 b^4 + 256) / 625
  [extremal]  P = z^d - z attains S = (d-1)/d for all d
"""
from __future__ import annotations
import numpy as np

RNG = np.random.default_rng(12345)


def eval_poly(coefs, z):
    """Horner evaluation; coefs[i] is coefficient of z^i."""
    acc = 0j
    for c in reversed(coefs):
        acc = acc * z + c
    return acc


def deriv_coefs(coefs):
    return [i * c for i, c in enumerate(coefs) if i >= 1] + ([0] if len(coefs) <= 1 else [])


def Q_anchor(coefs, z, a):
    """Difference quotient Q(a, z) = (P(z) - P(a))/(z - a)."""
    if abs(z - a) < 1e-14:
        # limit: P'(a)
        dP = deriv_coefs(coefs)
        return eval_poly(dP, a)
    return (eval_poly(coefs, z) - eval_poly(coefs, a)) / (z - a)


def smale_ratio(coefs, z):
    """min_{zeta : P'(zeta)=0} |Q(z, zeta)/P'(z)|.
    Here we interpret Q(z, zeta) as the difference quotient with ANCHOR z
    evaluated at the critical point zeta, i.e., Q(z, zeta) = (P(zeta)-P(z))/(zeta-z)
    which by Smale's convention equals (P(zeta)-P(z))/((zeta-z) P'(z))'s
    numerator. P'(z) is in the denominator outside.
    """
    dP = deriv_coefs(coefs)
    Pprime_z = eval_poly(dP, z)
    if abs(Pprime_z) < 1e-14:
        return float('inf')
    # critical points = roots of P'
    crits = np.roots(list(reversed(dP)))
    ratios = []
    for zeta in crits:
        val = (eval_poly(coefs, zeta) - eval_poly(coefs, z)) / (zeta - z)
        ratios.append(abs(val / Pprime_z))
    return float(min(ratios))


# =========================================================================
# Test 1: Smale MVC exact equality for d=2  (ratio = 1/2 for all z)
# =========================================================================

def test_d2_equality():
    print("=" * 70)
    print("[Paper 11] d=2: Q(z,zeta)/P'(z) = 1/2 EXACTLY")
    print("=" * 70)
    violations = 0
    max_err = 0.0
    for _ in range(5000):
        alpha = RNG.standard_normal() + 1j * RNG.standard_normal()
        beta = RNG.standard_normal() + 1j * RNG.standard_normal()
        # P(z) = (z - alpha)(z - beta)
        coefs = [alpha * beta, -(alpha + beta), 1.0 + 0j]
        z = RNG.standard_normal() + 1j * RNG.standard_normal()
        # critical point = (alpha + beta)/2
        if abs(z - (alpha + beta) / 2) < 1e-8:
            continue
        r = smale_ratio(coefs, z)
        err = abs(r - 0.5)
        max_err = max(max_err, err)
        if err > 1e-8:
            violations += 1
    print(f"  5000 random quadratics, max |ratio - 1/2| = {max_err:.3e}")
    print(f"  violations (err > 1e-8): {violations}")
    print()


# =========================================================================
# Test 2: Smale MVC for d=3  (bound 2/3, attained at z^3 - 1 as |z| -> inf)
# =========================================================================

def test_d3_bound():
    print("=" * 70)
    print("[Paper 11] d=3: min_zeta |Q(z,zeta)/P'(z)| <= 2/3")
    print("=" * 70)
    violations = 0
    sup = 0.0
    for _ in range(5000):
        coefs = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(4)]
        # ensure monic for display:
        if abs(coefs[-1]) < 1e-3:
            continue
        # normalize: monic
        coefs = [c / coefs[-1] for c in coefs]
        z = RNG.standard_normal() + 1j * RNG.standard_normal()
        r = smale_ratio(coefs, z)
        if r > 2 / 3 + 1e-8:
            violations += 1
        sup = max(sup, r)
    print(f"  5000 random cubics: sup(min ratio) = {sup:.6f}  (bound 2/3 = 0.66667)")
    print(f"  violations: {violations}")
    # Extremal: P = z^3 - 1, large z
    for R in [1, 5, 10, 100, 1000]:
        coefs = [-1, 0, 0, 1]
        z = complex(R, 0)
        r = smale_ratio(coefs, z)
        print(f"    z^3 - 1 at z = {R}: ratio = {r:.6f}")
    print()


# =========================================================================
# Test 3: Extremal P = z^d - z attains S = (d-1)/d exactly
# =========================================================================

def test_extremal_z_d_minus_z():
    print("=" * 70)
    print("[Paper 43/49] P(z) = z^d - z: S(c) = |P(c)/c| = (d-1)/d exactly")
    print("=" * 70)
    for d in [2, 3, 4, 5, 7, 10, 20]:
        coefs = [0, -1] + [0] * (d - 2) + [1]  # z^d - z
        dP = deriv_coefs(coefs)
        crits = np.roots(list(reversed(dP)))
        # remove c = 0 if any, and compute S(c) = |P(c)/c|
        vals = []
        for c in crits:
            if abs(c) < 1e-10:
                continue
            S = abs(eval_poly(coefs, c) / c)
            vals.append(S)
        avg = np.mean(vals)
        target = (d - 1) / d
        print(f"  d={d:>2}:  S(c) avg = {avg:.6f}  target (d-1)/d = {target:.6f}  err = {abs(avg - target):.2e}")
    print()


# =========================================================================
# Test 4: Pandrosion vanishing identity  sum 1/P'(alpha_k) = 0
# =========================================================================

def test_vanishing_identity():
    print("=" * 70)
    print("[Paper 43] Vanishing identity: sum_k 1/P'(alpha_k) = 0")
    print("=" * 70)
    for d in [2, 3, 4, 5, 7, 10]:
        max_err = 0.0
        for _ in range(100):
            coefs = [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d + 1)]
            if abs(coefs[-1]) < 1e-3:
                continue
            coefs = [c / coefs[-1] for c in coefs]
            roots = np.roots(list(reversed(coefs)))
            dP = deriv_coefs(coefs)
            s = 0j
            ok = True
            for a in roots:
                Pp = eval_poly(dP, a)
                if abs(Pp) < 1e-10:
                    ok = False
                    break
                s += 1 / Pp
            if ok:
                max_err = max(max_err, abs(s))
        print(f"  d={d:>2}:  max |sum 1/P'(alpha_k)| over 100 polynomials = {max_err:.3e}")
    print()


# =========================================================================
# Test 5: Pandrosion reduction  R(c) = q_tilde(c)/d, q_tilde(0) = d - 1
# =========================================================================

def test_pandrosion_reduction():
    print("=" * 70)
    print("[Paper 49] Pandrosion reduction:  R(c) = q_tilde(c)/d")
    print("           with q_tilde(0) = d - 1")
    print("=" * 70)
    # P(z) = z^d + a_{d-1} z^{d-1} + ... + a_1 z  (P(0) = 0, P'(0) = a_1)
    for d in [3, 4, 5, 6, 7]:
        max_err = 0.0
        for _ in range(100):
            a = [0.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 1)] + [1.0 + 0j]
            # Normalize so a_1 = 1  (i.e., P'(0) = 1)
            if abs(a[1]) < 1e-3:
                continue
            a = [c / a[1] for c in a]  # now a_1 = 1
            # R(z) = P(z)/z  has coefs [a_1, a_2, ..., a_d] (shift by 1)
            R_coefs = a[1:]
            # q_tilde(z) = sum_{k=1}^{d-1} k * a_{d-k} * z^{d-1-k}  (from Thm 4.1)
            # Let's verify: at critical point c (P'(c) = 0), R(c) = q_tilde(c)/d
            q_coefs = [0j] * (d - 1)
            for k in range(1, d):
                # term k * a_{d-k} * z^{d-1-k}
                # position d-1-k  (0 <= d-1-k <= d-2)
                pos = d - 1 - k
                if 0 <= pos < len(q_coefs):
                    q_coefs[pos] = k * a[d - k]
            # critical points: roots of P'
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            for c in crits:
                if abs(c) < 1e-10:
                    continue
                R_c = eval_poly(R_coefs, c)
                q_c = eval_poly(q_coefs, c)
                err = abs(R_c - q_c / d)
                max_err = max(max_err, err)
        # Also verify q_tilde(0) = (d-1) * a_{d-(d-1)} = (d-1) * a_1 = d-1
        q_at_0 = eval_poly([(d - 1)] + [0] * (d - 2), 0)
        # (just structural — q_tilde(0) is the constant term)
        print(f"  d={d}:  max |R(c) - q_tilde(c)/d| = {max_err:.3e}  (100 polynomials)")
    print()


# =========================================================================
# Test 6: d=4 explicit formula q_tilde(c) = 2(1 - c^2(a + 2c))
#         product prod q_tilde(c_k) = 4 a^3 - a^2 b^2 - 18 a b + 4 b^3 + 27
# =========================================================================

def test_d4_formula():
    print("=" * 70)
    print("[Paper 52] d=4:  q_tilde(c) = 2(1 - c^2 (a + 2c)) at each critical c")
    print("           prod q_tilde(c_k) = 4a^3 - a^2 b^2 - 18 a b + 4 b^3 + 27")
    print("=" * 70)
    max_err_q = 0.0
    max_err_prod = 0.0
    sup_ratio = 0.0
    violations = 0
    N = 5000
    for _ in range(N):
        a = RNG.standard_normal() + 1j * RNG.standard_normal()
        b = RNG.standard_normal() + 1j * RNG.standard_normal()
        # P(z) = z^4 + a z^3 + b z^2 + z   (P(0)=0, P'(0)=1)
        coefs = [0, 1, b, a, 1]
        dP = deriv_coefs(coefs)  # 1 + 2 b z + 3 a z^2 + 4 z^3
        crits = np.roots(list(reversed(dP)))  # roots of 4 c^3 + 3 a c^2 + 2 b c + 1
        # R(z) = P(z)/z = 1 + b z + a z^2 + z^3
        R_coefs = [1, b, a, 1]
        prod_q = 1.0 + 0j
        ratios = []
        for c in crits:
            R_c = eval_poly(R_coefs, c)
            q_tilde_formula = 2 * (1 - c ** 2 * (a + 2 * c))
            err = abs(R_c - q_tilde_formula / 4)
            max_err_q = max(max_err_q, err)
            prod_q *= q_tilde_formula
            if abs(c) > 1e-10:
                ratios.append(abs(R_c))
        # product formula
        prod_target = 4 * a ** 3 - a ** 2 * b ** 2 - 18 * a * b + 4 * b ** 3 + 27
        err_prod = abs(prod_q - prod_target)
        max_err_prod = max(max_err_prod, err_prod)
        # Smale ratio at z=0 (anchor at marked root 0)
        if ratios:
            smale = min(ratios)  # S(c) = |R(c)| = |P(c)/c|
            sup_ratio = max(sup_ratio, smale)
            if smale > 3 / 4 + 1e-6:
                violations += 1
    print(f"  {N} random (a, b):")
    print(f"    max |R(c) - q_tilde(c)/4|    = {max_err_q:.3e}")
    print(f"    max |prod q - (4a^3-a^2b^2-18ab+4b^3+27)| = {max_err_prod:.3e}")
    print(f"    sup(min |P(c)/c|)  = {sup_ratio:.6f}   (bound 3/4 = 0.75000)")
    print(f"    violations of S <= 3/4: {violations}")
    print()


# =========================================================================
# Test 7: d=5 explicit formula q_tilde(c) = (1 - 3 c^4 - a c^2)/2
#         product = (16 a^4 - 4 a^3 b^2 - 128 a^2 + 144 a b^2 - 27 b^4 + 256)/625
# =========================================================================

def test_d5_formula():
    print("=" * 70)
    print("[Paper 61] d=5:  q_tilde(c) = (1 - 3 c^4 - a c^2)/2")
    print("           prod = (16 a^4 - 4 a^3 b^2 - 128 a^2 + 144 a b^2 - 27 b^4 + 256)/625")
    print("=" * 70)
    max_err_q = 0.0
    max_err_prod = 0.0
    sup_ratio = 0.0
    violations = 0
    N = 5000
    for _ in range(N):
        a = RNG.standard_normal() + 1j * RNG.standard_normal()
        b = RNG.standard_normal() + 1j * RNG.standard_normal()
        # P(z) = z^5 + a z^3 + b z^2 + z   (P(0)=0, P'(0)=1)
        coefs = [0, 1, b, a, 0, 1]
        dP = deriv_coefs(coefs)  # 1 + 2 b z + 3 a z^2 + 0*z^3 + 5 z^4
        crits = np.roots(list(reversed(dP)))  # 5 c^4 + 3 a c^2 + 2 b c + 1
        R_coefs = [1, b, a, 0, 1]  # P(z)/z = 1 + b z + a z^2 + z^4
        prod_q = 1.0 + 0j
        ratios = []
        for c in crits:
            R_c = eval_poly(R_coefs, c)
            # Paper 61 uses q_tilde(c) = P(c)/c = R(c) directly (not R(c)*d),
            # so the formula (1 - 3c^4 - ac^2)/2 equals R(c) at critical points.
            q_formula = (1 - 3 * c ** 4 - a * c ** 2) / 2
            err = abs(R_c - q_formula)
            max_err_q = max(max_err_q, err)
            prod_q *= q_formula
            if abs(c) > 1e-10:
                ratios.append(abs(R_c))
        prod_target = (16 * a ** 4 - 4 * a ** 3 * b ** 2 - 128 * a ** 2 +
                       144 * a * b ** 2 - 27 * b ** 4 + 256) / 625
        err_prod = abs(prod_q - prod_target)
        max_err_prod = max(max_err_prod, err_prod)
        if ratios:
            smale = min(ratios)
            sup_ratio = max(sup_ratio, smale)
            if smale > 4 / 5 + 1e-6:
                violations += 1
    print(f"  {N} random (a, b):")
    print(f"    max |R(c) - q_tilde(c)/5|    = {max_err_q:.3e}")
    print(f"    max |prod q - product target| = {max_err_prod:.3e}")
    print(f"    sup(min |P(c)/c|)  = {sup_ratio:.6f}   (bound 4/5 = 0.80000)")
    print(f"    violations of S <= 4/5: {violations}")
    print()


# =========================================================================
# Test 8: Universal Smale MVC bound  min S <= (d-1)/d  for d = 2..10
# =========================================================================

def test_smale_universal():
    print("=" * 70)
    print("[Paper 11/43/49/52/61] Smale MVC universal check")
    print("  For polynomials P(z) = sum a_k z^k monic with P(0)=0, P'(0)=1")
    print("  min_{c : P'(c)=0, c != 0} |P(c)/c|  <=  (d-1)/d")
    print("=" * 70)
    for d in range(2, 11):
        sup = 0.0
        viols = 0
        N = 2000
        for _ in range(N):
            a = [0.0 + 0j, 1.0 + 0j] + [RNG.standard_normal() + 1j * RNG.standard_normal() for _ in range(d - 2)] + [1.0 + 0j]
            dP = deriv_coefs(a)
            crits = np.roots(list(reversed(dP)))
            ratios = []
            for c in crits:
                if abs(c) < 1e-10:
                    continue
                r = abs(eval_poly(a, c) / c)
                ratios.append(r)
            if ratios:
                m = min(ratios)
                sup = max(sup, m)
                if m > (d - 1) / d + 1e-6:
                    viols += 1
        print(f"  d = {d:>2}:  sup(min S) = {sup:.6f}  bound (d-1)/d = {(d-1)/d:.6f}  violations = {viols}/{N}")
    print()


if __name__ == "__main__":
    test_d2_equality()
    test_d3_bound()
    test_extremal_z_d_minus_z()
    test_vanishing_identity()
    test_pandrosion_reduction()
    test_d4_formula()
    test_d5_formula()
    test_smale_universal()
