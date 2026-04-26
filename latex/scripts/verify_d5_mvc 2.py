"""
Independent verification of Smale's MVC at d=5 on the centered slice.
=====================================================================

We re-derive everything from scratch, without copying the existing scripts in
latex/scripts/. The aim is to verify (or refute) the proof candidate of
6pandrosion_smale.tex in three layers:

   Layer A (pure algebra, sympy): prove the Pandrosion reduction,
       qtil(c) = (1 - 3c^4 - a c^2) / 2  at every critical point of P',
       and the resultant identity prod_j qtil(c_j) = N(a,b)/625.
   Layer B (symbolic Hessian, sympy): compute the full first- and
       second-order Taylor of  Phi_k(d) = |qtil(c_k(d))|^2  in the
       four real coordinates  d = (da_re, da_im, db_re, db_im)
       around the extremal (a, b) = (0, 0). The paper's first-order
       certificate covers three of four axes; we extend to second
       order to attack the cone obstruction around the fourth axis.
   Layer C (numerical stress test, numpy + scipy.optimize):
       - dense grid;
       - random sampling at multiple scales in C^2;
       - aggressive Nelder-Mead and SLSQP maximisation of
         min_j |qtil(c_j(a,b))|  to find the supremum on the slice.
       Goal: any (a, b) where the ratio exceeds 4/5 = 0.8 would
       refute the centered-slice conjecture.

The slice studied is  P(z) = z^5 + a z^3 + b z^2 + z,  with marked point
0 and normalisation P'(0) = 1.

This file is independent of the existing scripts. The structure was
inspired by them but every formula is rederived here.
"""

from __future__ import annotations
import time
from typing import List, Tuple
import numpy as np
import sympy as sp
import mpmath as mp


# =====================================================================
# Layer A — pure algebra in sympy, from scratch
# =====================================================================

def layerA_pandrosion_reduction() -> None:
    print("=" * 78)
    print("LAYER A.1 — Pandrosion reduction (qtil at critical points)")
    print("=" * 78)
    z, a, b = sp.symbols('z a b')
    P = z ** 5 + a * z ** 3 + b * z ** 2 + z
    Pp = sp.diff(P, z)
    print(f"  P(z)   = {P}")
    print(f"  P'(z)  = {Pp}")
    # At a critical point c we have 5c^4 + 3a c^2 + 2b c + 1 = 0.
    # We want to show that  P(c)/c  equals (1 - 3 c^4 - a c^2)/2.
    c = sp.symbols('c')
    Pc_over_c = sp.expand(P.subs(z, c) / c)             # c^4 + a c^2 + b c + 1
    qtil_claim = (1 - 3 * c ** 4 - a * c ** 2) / 2
    diff = sp.expand(Pc_over_c - qtil_claim)
    print(f"  P(c)/c = {Pc_over_c}")
    print(f"  qtil   = {qtil_claim}  [claim of paper 6, Theorem 1.1]")
    print(f"  difference (without using P'(c)=0) = {diff}")
    # Now reduce modulo the critical equation
    crit = 5 * c ** 4 + 3 * a * c ** 2 + 2 * b * c + 1
    # diff should be a multiple of crit
    quot, rem = sp.div(diff, crit, c)
    print(f"  diff = ({quot}) * (5c^4+3ac^2+2bc+1)  +  ({rem})")
    print(f"  remainder modulo critical eq.: {sp.simplify(rem)}")
    if sp.simplify(rem) == 0:
        print("  => Layer A.1 PASS: qtil(c) = (1 - 3c^4 - ac^2)/2 at critical points.")
    else:
        print("  => FAIL")


def layerA_resultant_product() -> None:
    print("\n" + "=" * 78)
    print("LAYER A.2 — Resultant product formula  prod qtil(c_j) = N(a,b)/625")
    print("=" * 78)
    z, a, b = sp.symbols('z a b')
    Pp = 5 * z ** 4 + 3 * a * z ** 2 + 2 * b * z + 1     # critical polynomial
    qtil = (1 - 3 * z ** 4 - a * z ** 2) / 2
    # prod over critical points c_j of qtil(c_j) = Res(qtil, Pp) / lc(Pp)^deg(qtil)
    # (with appropriate sign / normalisation). Compute via Sylvester resultant.
    R = sp.resultant(qtil, Pp, z)                        # symbolic resultant
    R_simplified = sp.expand(R)
    # qtil has leading coeff -3/2, degree 4.  Pp has leading coeff 5, degree 4.
    # Res(f, g) = lc(f)^{deg g} * prod_{f(c)=0} g(c)
    # We want  prod_{Pp(c)=0} qtil(c)
    #        = Res(Pp, qtil) / lc(Pp)^{deg qtil}
    #        = Res(Pp, qtil) / 5^4
    # Note: Res(f, g) = (-1)^{deg f * deg g} Res(g, f).
    R2 = sp.resultant(Pp, qtil, z)
    prod_qtil = sp.expand(R2 / sp.Integer(5) ** 4)       # 5^4 = 625
    print(f"  Sylvester resultant / 625 (symbolic):")
    print(f"    {prod_qtil}")
    # Compare with paper 6's claim: N(a,b)/625 with
    # N(a,b) = 16 a^4 - 4 a^3 b^2 - 128 a^2 + 144 a b^2 - 27 b^4 + 256
    N = 16 * a ** 4 - 4 * a ** 3 * b ** 2 - 128 * a ** 2 + 144 * a * b ** 2 \
        - 27 * b ** 4 + 256
    claim = N / 625
    discrepancy = sp.expand(prod_qtil - claim)
    print(f"  Paper 6 claim N(a,b)/625 = {claim}")
    print(f"  computed - claim         = {discrepancy}")
    if sp.simplify(discrepancy) == 0:
        print("  => Layer A.2 PASS: resultant product matches paper 6 exactly.")
    else:
        print("  => FAIL or scaling issue")
    # Numerical sanity at (a,b) = (0,0): N(0,0)/625 = 256/625
    print(f"  At (a,b)=(0,0): N/625 = 256/625 = {sp.Rational(256, 625)} = (4/5)^4")


# =====================================================================
# Layer B — symbolic 1st + 2nd order Taylor of Phi_k = |qtil|^2 at (0,0)
# =====================================================================

def layerB_full_hessian():
    """
    Compute the full second-order Taylor expansion of
        Phi_k(d) = |qtil(c_k(d))|^2
    around d = 0, where d = (da_re, da_im, db_re, db_im).

    The paper's existing analysis stops at first order (3 of 4 axes covered)
    plus an axis-only quadratic computation. Here we go further: we get the
    full real Hessian H_k for each k, and then verify whether the Hessian
    on the kernel of the linear forms (i.e. the da_re axis and a small cone
    around it) is uniformly negative-definite.
    """
    print("\n" + "=" * 78)
    print("LAYER B — Full 1st + 2nd order Taylor of Phi_k = |qtil|^2 at (0,0)")
    print("=" * 78)

    da_re, da_im, db_re, db_im = sp.symbols('da_re da_im db_re db_im', real=True)
    eps = sp.symbols('eps', positive=True)
    # Real perturbation at order eps
    da = eps * (da_re + sp.I * da_im)
    db = eps * (db_re + sp.I * db_im)
    # Critical equation: 5 c^4 + 3 a c^2 + 2 b c + 1 = 0
    # At eps = 0:  5 c0^4 = -1, so c0^4 = -1/5.
    # Four critical points  c0_k = rho * w^(2k+1)  with rho = 5^(-1/4),
    # w = exp(i pi/4).  Equivalent: c0_k = (-1/5)^(1/4) * i^k.
    rho = sp.Rational(1, 5) ** sp.Rational(1, 4)
    omega = sp.exp(sp.I * sp.pi * sp.Rational(1, 4))
    c0_list = [rho * omega * sp.exp(sp.I * sp.pi * sp.Rational(k, 2))
               for k in range(4)]
    # Series expansion of c_k(eps) to order eps^2
    # c = c0 + eps * c1 + eps^2 * c2 + ...
    # Substitute into 5c^4 + 3 a c^2 + 2 b c + 1 = 0 and expand.
    print("  (computing series for each c_k to order eps^2 ...)")
    Hessian_results = []
    rates = []     # symbolic linear forms in (da_re, da_im, db_re, db_im)
    quads = []     # symbolic quadratic forms in same vars
    for k, c0 in enumerate(c0_list):
        c1 = sp.symbols(f'c1_{k}')
        c2 = sp.symbols(f'c2_{k}')
        c_series = c0 + eps * c1 + eps ** 2 * c2
        eq = 5 * c_series ** 4 + 3 * da * c_series ** 2 + 2 * db * c_series + 1
        eq = sp.expand(eq)
        # Coefficient of eps^1: linear equation in c1
        coef1 = sp.expand(eq.coeff(eps, 1))
        # Solve for c1 in terms of (da_re, da_im, db_re, db_im)
        c1_sol = sp.solve(coef1, c1)[0]
        # Coefficient of eps^2: quadratic + linear-in-c2 equation
        coef2 = sp.expand(eq.coeff(eps, 2)).subs(c1, c1_sol)
        c2_sol = sp.solve(coef2, c2)[0]
        # Now compute qtil(c, a) = (1 - 3 c^4 - a c^2)/2 and its expansion
        a_full = da
        qtil_full = (1 - 3 * c_series ** 4 - a_full * c_series ** 2) / 2
        qtil_full = qtil_full.subs([(c1, c1_sol), (c2, c2_sol)])
        qtil_full = sp.expand(sp.series(qtil_full, eps, 0, 3).removeO())
        # |qtil|^2 = qtil * conj(qtil); at eps=0 this is (4/5)^2 = 16/25
        # but qtil is *complex*, so we use sympy's conjugate carefully.
        # Since (da_re, da_im, db_re, db_im) are declared REAL, conj(da)=da_re-i*da_im.
        qtil_bar = sp.conjugate(qtil_full)
        Phi = sp.expand(qtil_full * qtil_bar)
        # Series in eps to order 2
        Phi = sp.expand(sp.series(Phi, eps, 0, 3).removeO())
        # Extract coefficients
        Phi0 = Phi.subs(eps, 0)                                      # = 16/25
        Phi1 = sp.expand(sp.diff(Phi, eps).subs(eps, 0))             # linear form
        Phi2 = sp.expand(sp.diff(Phi, eps, 2).subs(eps, 0)) / 2     # quadratic form
        Phi0 = sp.simplify(Phi0)
        Phi1 = sp.simplify(Phi1)
        Phi2 = sp.simplify(Phi2)
        Hessian_results.append((Phi0, Phi1, Phi2))
        rates.append(Phi1)
        quads.append(Phi2)
        if k == 0:
            print(f"  c_0 = {sp.nsimplify(c0)}")
            print(f"     Phi(0)   = {Phi0}                (must be 16/25)")
            print(f"     Phi'(0)  = {Phi1}")
            print(f"     Phi''(0)/2 = {Phi2}")
    # Sanity: Phi0 = 16/25 for all k
    for k, (P0, _, _) in enumerate(Hessian_results):
        assert sp.simplify(P0 - sp.Rational(16, 25)) == 0, f"Phi_0({k}) wrong"
    print("  All 4 Phi_k(0) = 16/25 confirmed.")
    # Sum of 4 linear rates: should be zero
    sum_rates = sp.simplify(sum(rates))
    print(f"\n  sum_k (rate_k) = {sum_rates}     (must be 0 by C4 symmetry)")
    # Sum of 4 quadratic rates: a real quadratic form Q(d)
    sum_quads = sp.simplify(sum(quads))
    print(f"  sum_k (quadratic_k) (raw) ="
          f"\n    {sp.expand(sum_quads)}")
    # Now compute min_k Phi_k(d) at second order.
    # We use the structure: at second order,
    #     Phi_k(d) = 16/25 + L_k(d) + Q_k(d) + O(d^3),
    # so min_k Phi_k(d) <= min_k (16/25 + L_k + Q_k).
    return rates, quads


def layerB_axis_quadratic(rates: List[sp.Expr], quads: List[sp.Expr]) -> None:
    """On the da_re axis, all L_k vanish; what is min_k Q_k(da_re)?"""
    print("\n  -- Restriction to the pure da_re axis --")
    da_re, da_im, db_re, db_im = sp.symbols('da_re da_im db_re db_im', real=True)
    sub_axis = {da_im: 0, db_re: 0, db_im: 0}
    print(f"    rate_k on (da_re, 0, 0, 0):")
    for k, r in enumerate(rates):
        print(f"      k={k}: {sp.simplify(r.subs(sub_axis))}")
    print(f"    quadratic_k on (da_re, 0, 0, 0):")
    quad_axis = []
    for k, q in enumerate(quads):
        v = sp.simplify(q.subs(sub_axis))
        quad_axis.append(v)
        print(f"      k={k}: {v}")
    # The minimum of the four quadratic coefficients, divided by da_re^2:
    print(f"\n    min over k of [quad_k / da_re^2]: {sp.simplify(sp.Min(*[v / da_re ** 2 for v in quad_axis]))}")


def layerB_cone_certificate(rates: List[sp.Expr], quads: List[sp.Expr]) -> None:
    """
    Numerical cone certificate around the da_re axis.
    Strategy: parameterise directions on the unit sphere in R^4 and check that

         min_k  [ L_k(d) + Q_k(d) ]  <  - C ||d||^2

    for some absolute constant C > 0, uniformly in directions.

    If this holds, the second-order Taylor descent is UNIFORM — closing the
    "second-order cone" gap of paper 6 numerically (still not a symbolic SOS
    proof, but a strong quantitative certificate).
    """
    print("\n" + "=" * 78)
    print("LAYER B+ — Numerical cone certificate (full 2nd-order test)")
    print("=" * 78)
    da_re, da_im, db_re, db_im = sp.symbols('da_re da_im db_re db_im', real=True)
    # Lambdify the four functions  16/25 + L_k + Q_k
    expr_k = [sp.Rational(16, 25) + r + q for r, q in zip(rates, quads)]
    fns = [sp.lambdify((da_re, da_im, db_re, db_im), e, 'numpy')
           for e in expr_k]
    # Sample a dense grid on S^3
    rng = np.random.default_rng(42)
    N = 200_000
    # Uniform on S^3
    d = rng.standard_normal((N, 4))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    # Evaluate min_k Phi_k(d) at unit norm; record max
    # Phi_k = 16/25 + 1*L_k(d) + 1^2*Q_k(d) since Taylor is in eps=1
    # but we sample at unit norm, so eps = 1.
    vals = np.stack([f(d[:, 0], d[:, 1], d[:, 2], d[:, 3]) for f in fns],
                    axis=1)
    min_per_dir = vals.min(axis=1)
    sup_min = min_per_dir.max()
    arg = d[min_per_dir.argmax()]
    print(f"  Sampled {N} unit directions in R^4 (S^3).")
    print(f"  sup over directions of [16/25 + L_min + Q_min] (truncated 2nd order):")
    print(f"     value = {sup_min:.10f}    (16/25 = {16/25:.10f})")
    print(f"     gap from 16/25 = {sup_min - 16/25:+.4e}")
    print(f"     attained at d ≈ ({arg[0]:+.4f}, {arg[1]:+.4f}, "
          f"{arg[2]:+.4f}, {arg[3]:+.4f})")
    # The cone certificate would say: this sup is < 16/25, with margin C * ||d||^2.
    # Since we sample at ||d||=1, sup < 16/25 means second-order descent uniform.
    if sup_min < 16 / 25 - 1e-6:
        margin = 16 / 25 - sup_min
        print(f"  => 2nd-order CONE CERTIFICATE *appears* to hold:")
        print(f"     uniform descent with constant C >= {margin:.4f} per ||d||^2.")
        print(f"  (Caveat: sampling-based, not a symbolic SOS proof.)")
    else:
        print(f"  WARNING: at second order alone we do NOT see uniform descent.")
        print(f"  Need to check higher-order terms or use SOS decomposition.")
    # Also: refine search via local optimisation around the worst direction
    return sup_min, arg


# =====================================================================
# Layer C — independent numerical stress test
# =====================================================================

def critical_pts(a: complex, b: complex) -> np.ndarray:
    """Roots of 5z^4 + 3 a z^2 + 2 b z + 1 = 0."""
    return np.roots([5.0, 0.0, 3.0 * a, 2.0 * b, 1.0])


def qtil(c: complex, a: complex) -> complex:
    return (1.0 - 3.0 * c ** 4 - a * c ** 2) / 2.0


def smale_score(a: complex, b: complex) -> float:
    """min_j |qtil(c_j(a,b))|"""
    cs = critical_pts(a, b)
    return float(min(abs(qtil(c, a)) for c in cs))


def layerC_grid_and_random() -> Tuple[float, complex, complex]:
    print("\n" + "=" * 78)
    print("LAYER C.1 — Independent grid + random sweep (numpy)")
    print("=" * 78)
    rng = np.random.default_rng(2026)
    sup_global = 0.0
    arg_global = (0+0j, 0+0j)

    # (a) Real grid [-2,2]^2
    n = 250
    grid = np.linspace(-2.0, 2.0, n)
    sup_local = 0.0
    arg_local = (0.0, 0.0)
    for ar in grid:
        for br in grid:
            s = smale_score(ar, br)
            if s > sup_local:
                sup_local = s
                arg_local = (ar, br)
    print(f"  Real grid [-2,2]^2 ({n}x{n}):  sup = {sup_local:.10f}, "
          f"at (a,b)=({arg_local[0]:+.4f}, {arg_local[1]:+.4f})")
    if sup_local > sup_global:
        sup_global = sup_local
        arg_global = (complex(arg_local[0]), complex(arg_local[1]))

    # (b) Complex random at multiple scales
    for scale in [0.05, 0.5, 5.0, 50.0, 500.0]:
        sup_s = 0.0
        arg_s = (0+0j, 0+0j)
        for _ in range(8000):
            a = (rng.standard_normal() + 1j * rng.standard_normal()) * scale
            b = (rng.standard_normal() + 1j * rng.standard_normal()) * scale
            s = smale_score(a, b)
            if s > sup_s:
                sup_s = s
                arg_s = (a, b)
        print(f"  Complex random scale={scale:>6.2f}: sup = {sup_s:.10f}")
        if sup_s > sup_global:
            sup_global = sup_s
            arg_global = arg_s

    return sup_global, arg_global[0], arg_global[1]


def layerC_local_optimisation(a0: complex, b0: complex) -> Tuple[float, complex, complex]:
    """Aggressive local maximisation of smale_score using SLSQP / Nelder-Mead."""
    from scipy.optimize import minimize
    print("\n" + "=" * 78)
    print("LAYER C.2 — Local maximisation (SLSQP + Nelder-Mead, multistart)")
    print("=" * 78)

    def neg_score(params):
        a = params[0] + 1j * params[1]
        b = params[2] + 1j * params[3]
        try:
            return -smale_score(a, b)
        except Exception:
            return 0.0

    rng = np.random.default_rng(7)
    sup = 0.0
    best_arg = (0+0j, 0+0j)
    starts = [
        (0.0, 0.0, 0.0, 0.0),
        (0.01, 0.0, 0.0, 0.0),
        (a0.real, a0.imag, b0.real, b0.imag),
    ]
    # plus 50 random starts at multiple scales
    for scale in [0.01, 0.1, 1.0, 10.0]:
        for _ in range(15):
            starts.append(tuple(rng.standard_normal(4) * scale))
    for x0 in starts:
        try:
            res = minimize(neg_score, x0, method='Nelder-Mead',
                           options={'xatol': 1e-10, 'fatol': 1e-10,
                                    'maxiter': 5000})
            score = -res.fun
            if score > sup:
                sup = score
                a = res.x[0] + 1j * res.x[1]
                b = res.x[2] + 1j * res.x[3]
                best_arg = (a, b)
        except Exception:
            continue
    a, b = best_arg
    print(f"  Multi-start Nelder-Mead sup over slice: {sup:.12f}")
    print(f"  Compared with 4/5 = 0.800000000000")
    print(f"  Margin to 4/5: {4/5 - sup:+.4e}")
    print(f"  Best (a, b) ≈ ({a.real:+.5f} + {a.imag:+.5f}j, "
          f"{b.real:+.5f} + {b.imag:+.5f}j)")
    if sup <= 4.0 / 5.0 + 1e-9:
        print("  => Layer C: NO violation of the 4/5 bound discovered.")
    else:
        print(f"  => Layer C: ALERT, value above 4/5 by {sup - 4/5:.4e}!")
    return sup, a, b


def layerC_high_precision_extremal() -> None:
    """At (a,b)=(0,0), confirm the four |qtil(c_k)| are exactly 4/5 to 50 digits."""
    print("\n" + "=" * 78)
    print("LAYER C.3 — High-precision check at the extremal (mpmath)")
    print("=" * 78)
    mp.mp.dps = 60
    a_mp, b_mp = mp.mpf(0), mp.mpf(0)
    coeffs = [mp.mpf(5), mp.mpf(0), 3 * a_mp, 2 * b_mp, mp.mpf(1)]
    # Solve quartic 5c^4 + 1 = 0 by hand: c^4 = -1/5
    c4 = -mp.mpf(1) / 5
    cs = [mp.power(-mp.mpf(1) / 5, mp.mpf(1) / 4) * mp.exp(1j * mp.pi * mp.mpf(k) / 2)
          for k in range(4)]
    print(f"  At (0,0): the 4 critical points are c_k = (-1/5)^(1/4) * i^k.")
    qs = []
    for k, c in enumerate(cs):
        q = (1 - 3 * c ** 4 - a_mp * c ** 2) / 2
        qs.append(q)
        print(f"    k={k}: |qtil(c_k)| = {mp.nstr(abs(q), 30)}")
    target = mp.mpf(4) / 5
    max_err = max(abs(abs(q) - target) for q in qs)
    print(f"  max |.|qtil| - 4/5| = {mp.nstr(max_err, 6)}")
    print(f"  => Extremal achieves 4/5 to within {mp.nstr(max_err, 6)}.")


# =====================================================================
# Main
# =====================================================================

def main() -> None:
    t0 = time.time()
    layerA_pandrosion_reduction()
    layerA_resultant_product()
    rates, quads = layerB_full_hessian()
    layerB_axis_quadratic(rates, quads)
    sup2, arg2 = layerB_cone_certificate(rates, quads)
    sup, a0, b0 = layerC_grid_and_random()
    sup_loc, _, _ = layerC_local_optimisation(a0, b0)
    layerC_high_precision_extremal()
    print("\n" + "=" * 78)
    print(f"OVERALL ELAPSED: {time.time()-t0:.1f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()
