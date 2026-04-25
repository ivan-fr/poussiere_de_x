"""
Symbolic first-order certificate for the d=5 centered-slice MVC program.

Goal: compute analytically, using sympy, the first-order variation of
Phi(a, b) = min_j |q̃(c_j(a,b))|^2 at the extremal (a, b) = (0, 0).

Important status:
    This is not a complete strict-local-maximum proof. The first-order certificate
    proves descent in the da_im, db_re, and db_im directions. The da_re direction
    has zero first-order rate; this script closes that one-axis obstruction by
    an exact b=0, a real computation:
        min |q_j|^2 = 16/25 - 4 a^2/25.
    A full uniform local proof still requires the mixed second-order estimate
    in a cone around the da_re axis.
"""
from __future__ import annotations
import sympy as sp


def main():
    # Real coordinates: a = a_re + i a_im, b = b_re + i b_im
    a_re, a_im, b_re, b_im = sp.symbols('a_re a_im b_re b_im', real=True)
    a = a_re + sp.I * a_im
    b = b_re + sp.I * b_im

    print("=" * 78, flush=True)
    print("Symbolic first-order rates for Phi(a,b)=min_j |q̃(c_j)|^2 at (0,0)", flush=True)
    print("=" * 78, flush=True)

    # Critical points at (0,0): c^4 = -1/5
    # The 4 critical points are c_k = rho * exp(i pi (2k+1)/4) where rho = 5^(-1/4)
    rho_val = sp.Rational(1, 5) ** sp.Rational(1, 4)
    print(f"\nAt (a,b)=(0,0): c^4 = -1/5, so |c| = 5^(-1/4) ≈ {float(rho_val):.4f}",
          flush=True)

    # First-order expansion of c around c_0 in (delta_a, delta_b):
    # P'(c) = 5c^4 + 3 a c^2 + 2 b c + 1 = 0
    # delta c = -(3 c_0^2 delta_a + 2 c_0 delta_b) / (20 c_0^3 + 0 + 0)
    #         = -(3 c_0^2 delta_a + 2 c_0 delta_b) / (20 c_0^3)
    #
    # q̃(c) = (1 - 3 c^4 - a c^2)/2
    # delta q̃ = (-12 c_0^3 delta_c - delta_a * c_0^2 - 0)/2
    #         = (-12 c_0^3 * delta_c - c_0^2 * delta_a)/2

    # We compute the linearised |q̃(c_k)|^2 for each critical point, then min.

    print("\n" + "-" * 78, flush=True)
    print("Step 1 (perturbation rate at each c_k for direction (delta_a, delta_b)):",
          flush=True)
    print("-" * 78, flush=True)

    # The 4 critical points (symbolic)
    c0_vals = [rho_val * sp.exp(sp.I * sp.pi * (2 * k + 1) / 4) for k in range(4)]

    da, db = sp.symbols('da db', complex=True)
    rates = []
    for k, c0 in enumerate(c0_vals):
        c0 = sp.simplify(c0)
        # delta c (first order)
        denom = 20 * c0 ** 3
        delta_c = -(3 * c0 ** 2 * da + 2 * c0 * db) / denom
        # q_tilde(c_0) = 4/5
        q0 = sp.Rational(4, 5)
        # delta q̃ (first order)
        delta_q = (-12 * c0 ** 3 * delta_c - c0 ** 2 * da) / 2
        # |q̃|^2 = q̃ * conj(q̃);  at q̃_0 = 4/5 (real positive),
        # delta |q̃|^2 = 2 * (4/5) * Re(delta q̃)
        rate_lin = 2 * q0 * sp.re(delta_q)
        rate_lin = sp.simplify(rate_lin)
        rates.append(rate_lin)
        print(f"  c_{k} = {sp.nsimplify(c0)}  ->  d|q̃|²/d(da,db) = {rate_lin}",
              flush=True)

    # The minimum over k: we need to show that for ANY (da, db) != 0, at least one
    # rate_k is strictly negative.
    print("\n" + "-" * 78, flush=True)
    print("Step 2 (min_k of linear forms = strictly negative for any direction)",
          flush=True)
    print("-" * 78, flush=True)

    # Substitute da = da_re + I da_im, db = db_re + I db_im
    da_re, da_im, db_re, db_im = sp.symbols('da_re da_im db_re db_im', real=True)
    da_sub = da_re + sp.I * da_im
    db_sub = db_re + sp.I * db_im

    rates_real = []
    for k, r in enumerate(rates):
        r_real = r.subs({da: da_sub, db: db_sub})
        r_real = sp.simplify(sp.expand(r_real))
        rates_real.append(r_real)
        print(f"  rate_{k}(da_re, da_im, db_re, db_im) = {r_real}", flush=True)

    # Sum of all 4 rates: should be negative for any direction (Hessian sum)
    sum_rates = sum(rates_real)
    sum_rates = sp.simplify(sum_rates)
    print(f"\n  sum_k rate_k = {sum_rates}", flush=True)

    # Sum of squares: positive contribution (definite-ness witness)
    sum_sq = sum(r ** 2 for r in rates_real)
    sum_sq = sp.simplify(sum_sq)
    print(f"\n  Sum of squares (witness of non-degeneracy):", flush=True)
    print(f"  sum rate_k^2 = ", sp.expand(sum_sq), flush=True)

    print("\n" + "-" * 78, flush=True)
    print("Step 3 (compactness: min over directions of unit norm)", flush=True)
    print("-" * 78, flush=True)

    # We want to show min_k rate_k(theta) < 0 for any theta on S^3.
    # Direction = (cos a1 cos a2, sin a1 cos a2, cos a3 sin a2, sin a3 sin a2)
    # in (da_re, da_im, db_re, db_im).
    # Equivalent: show max over directions of min_k rate_k < 0.

    # We compute the directional rate for several test directions
    # to verify negativity numerically.
    test_dirs = [
        ('(da_re=1, others=0)', {da_re: 1, da_im: 0, db_re: 0, db_im: 0}),
        ('(da_im=1)',           {da_re: 0, da_im: 1, db_re: 0, db_im: 0}),
        ('(db_re=1)',           {da_re: 0, da_im: 0, db_re: 1, db_im: 0}),
        ('(db_im=1)',           {da_re: 0, da_im: 0, db_re: 0, db_im: 1}),
        ('(1,1,1,1)/2',         {da_re: sp.Rational(1, 2),
                                 da_im: sp.Rational(1, 2),
                                 db_re: sp.Rational(1, 2),
                                 db_im: sp.Rational(1, 2)}),
        ('(1,-1,1,-1)/2',       {da_re: sp.Rational(1, 2),
                                 da_im: -sp.Rational(1, 2),
                                 db_re: sp.Rational(1, 2),
                                 db_im: -sp.Rational(1, 2)}),
    ]
    print(f"\n  {'direction':<25} {'rates per c_k':<60}  {'min':>10}", flush=True)
    print(f"  {'-' * 25} {'-' * 60}  {'-' * 10}", flush=True)
    all_min_negative = True
    for name, subs in test_dirs:
        vals = [float(r.subs(subs)) for r in rates_real]
        m = min(vals)
        rates_str = ', '.join(f'{v:+.4f}' for v in vals)
        flag = ""
        if m >= 0:
            flag = " <-- NON-NEGATIVE!"
            all_min_negative = False
        print(f"  {name:<25} ({rates_str})  {m:>+8.4f}{flag}", flush=True)

    print(f"\n  Conclusion: first-order descent in every tested direction: "
          f"{'PASS' if all_min_negative else 'FAIL'}", flush=True)

    print("\n" + "=" * 78, flush=True)
    print("CERTIFICATE STATUS: partial first-order certificate", flush=True)
    print("=" * 78, flush=True)
    print("""
The four first-order rates sum to zero, and their sum of squares is

  1024/3125 * da_im^2 + 1152*sqrt(5)/3125 * (db_re^2 + db_im^2).

Therefore, whenever (da_im, db_re, db_im) is nonzero, at least one rate is
strictly negative and Phi decreases at first order in that direction.

The pure da_re direction has all first-order rates equal to zero. On that axis
we can compute exactly: with b=0 and a real, y=c^2 satisfies

  5 y^2 + 3 a y + 1 = 0,        q = 4/5 + (2/5) a y.

For the two y-values, y_1+y_2=-3a/5 and y_1 y_2=1/5, hence

  q_1 q_2 = 16/25 - 4 a^2/25.

For |a| < 2/sqrt(3), the two y-values are conjugate up to the real scaling and
the two q-values have equal modulus, so

  min_j |q_j|^2 = 16/25 - 4 a^2/25 < 16/25  for a != 0.

Thus the missing coordinate axis descends quadratically. This still does NOT
prove the full local maximum uniformly in mixed directions; that requires a
second-order cone estimate around the da_re axis.
""", flush=True)


if __name__ == "__main__":
    main()
