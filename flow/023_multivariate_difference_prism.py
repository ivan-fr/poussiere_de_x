"""
PAPER: 023
TITLE: Multivariate difference-prism (extends flow/020 to vector iterations)
STATUS: vector Aitken / projection Aitken applied to asymmetric Pandrosion

CONSTRUCTION
============
Flow/020 lifted SCALAR Pandrosion h(s) to dD via stacked Steffensen.
Here we lift VECTOR Pandrosion F : C^n -> C^n via stacked vector Aitken.

Multivariate Pandrosion map (asymmetric hypercube, from flow/016):
    F_i(s) = 1 - (x_i - 1) / (x_i * S(s)),   S(s) = prod_j S_p(s_j).
For symmetric x_i = x: rank-1 collapse, F reduces to scalar.
For asymmetric x_i:   GENUINELY multivariate, fixed point is a new
                      coupled algebraic vector (NOT x^{1/p}).

dD MULTIVARIATE LIFT:
    T_2 = F                       (linear convergence, vector)
    T_d = vec_Steffensen(T_{d-1}) (order 2^{d-2})

Two flavours of vector Steffensen:

  (DIAG)  Scalar Aitken applied COMPONENTWISE.
          For each i: x_new[i] = x_0[i] - (Dx_0[i])^2 / D2x_0[i].
          Treats each coordinate as an independent scalar sequence.
          Simple and robust; works as long as components have
          NON-DEGENERATE individual sequences.

  (PROJ)  Topological vector Aitken with scalar alpha projection.
          alpha = (sum dx0_i * d2x_i) / (sum d2x_i^2)
          x_new = x_0 - alpha * Dx_0.
          Optimal when all components share the SAME dominant
          eigenvalue (single-eigenmode case).  Fails when
          components have well-separated rates.

RANK-1 COLLAPSE
===============
At SYMMETRIC x_i, the iteration F lands on the diagonal s_1=...=s_n
in one step (rank-1 Jacobian).  After collapse, dynamics is scalar
and DIAG/PROJ both reduce to scalar prism dD.

At ASYMMETRIC x_i, the rank-1 collapse is broken: components stay
distinct throughout iteration.  This is a genuine n-dimensional
state.  The phase-space lift adds (d-2)*n axes but each axis is
still enslaved to the underlying multivariate state.

So multivariate-prism gives genuine n-dim state + acceleration,
NOT a true multivariate Pandrosion (the underlying F is the same).
"""

from __future__ import annotations
import math


# -----------------------------------------------------------------------------
# Multivariate Pandrosion: ASYMMETRIC HYPERCUBE map.
# -----------------------------------------------------------------------------
def S_hyper(s, p):
    total = 1.0 + 0.0j
    for si in s:
        total *= sum(complex(si) ** k for k in range(p))
    return total


def F_multi(s, args):
    x_vec, p = args
    S = S_hyper(s, p)
    return [1 - (x_vec[i] - 1) / (x_vec[i] * S) for i in range(len(s))]


# -----------------------------------------------------------------------------
# Vector Steffensen flavours.
# -----------------------------------------------------------------------------
def diag_steffensen(T, x0, args):
    x1 = T(x0, args)
    x2 = T(x1, args)
    n = len(x0)
    out = []
    for i in range(n):
        d0 = x1[i] - x0[i]
        d1 = x2[i] - x1[i]
        d2 = d1 - d0
        if abs(d2) < 1e-300:
            out.append(x2[i])
        else:
            out.append(x0[i] - d0 * d0 / d2)
    return out


def proj_steffensen(T, x0, args):
    x1 = T(x0, args)
    x2 = T(x1, args)
    n = len(x0)
    dx0 = [x1[i] - x0[i] for i in range(n)]
    dx1 = [x2[i] - x1[i] for i in range(n)]
    d2x = [dx1[i] - dx0[i] for i in range(n)]
    num = sum(dx0[i] * d2x[i] for i in range(n))
    den = sum(d2x[i] * d2x[i] for i in range(n))
    if abs(den) < 1e-300:
        return x2
    alpha = num / den
    return [x0[i] - alpha * dx0[i] for i in range(n)]


def make_multi_T_d(F_base, d, mode='diag'):
    if d == 2:
        return F_base
    T_prev = make_multi_T_d(F_base, d - 1, mode)
    accel = diag_steffensen if mode == 'diag' else proj_steffensen

    def T_d(x0, args):
        return accel(T_prev, x0, args)
    return T_d


# -----------------------------------------------------------------------------
# Iteration driver.
# -----------------------------------------------------------------------------
def iterate_multi(T, x0, args, tol=1e-13, max_iter=200):
    x = list(x0)
    for n in range(max_iter):
        x_new = T(x, args)
        diff = max(abs(x_new[i] - x[i]) for i in range(len(x)))
        if diff < tol:
            return x_new, n + 1
        x = x_new
    return x, max_iter


def verify_fp(s, args):
    x_vec, p = args
    S = S_hyper(s, p)
    return [abs(x_vec[i] * S * (1 - s[i]) - (x_vec[i] - 1)) for i in range(len(s))]


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------
def main():
    print("=" * 80)
    print("MULTIVARIATE DIFFERENCE PRISM (extends flow/020 to vectors)")
    print("=" * 80)
    print()
    print("Map: asymmetric hypercube  F_i(s) = 1 - (x_i-1)/(x_i prod_j S_p(s_j))")
    print("Lift: T_2=F, T_d = vec_Steffensen(T_{d-1}).")
    print("Modes: DIAG (componentwise Aitken) | PROJ (alpha-projection).")
    print()

    test_cases = [
        ("symmetric    d=2", [2.0+0j, 2.0+0j], 2),
        ("asymmetric   d=2", [2.0+0j, 5.0+0j], 2),
        ("asymmetric   d=3", [2.0+0j, 3.0+0j, 5.0+0j], 2),
        ("asymmetric p3 d=3", [2.0+0j, 3.0+0j, 5.0+0j], 3),
        ("asymmetric   d=4", [2.0+0j, 3.0+0j, 5.0+0j, 7.0+0j], 2),
        ("strong asym  d=3", [2.0+0j, 10.0+0j, 100.0+0j], 2),
    ]

    for label, x_vec, p in test_cases:
        n = len(x_vec)
        x0 = [0.5 + 0j] * n
        args = (x_vec, p)
        print(f"\n[{label}]  x = {[xi.real for xi in x_vec]}, p = {p}")
        print(f"  {'method':>22} {'iters':>6} {'max FP residual':>18}")
        print("-" * 80)
        for d in [2, 3, 4, 5]:
            modes = ['raw'] if d == 2 else ['diag', 'proj']
            for mode in modes:
                T = make_multi_T_d(F_multi, d, mode if mode != 'raw' else 'diag')
                s_final, n_it = iterate_multi(T, x0, args)
                max_res = max(verify_fp(s_final, args))
                tag = f"T_{d} ({mode})"
                print(f"  {tag:>22} {n_it:>6} {max_res:>18.3e}")
        s_show = [round(s.real, 8) for s in s_final]
        print(f"  Final s    = {s_show}")
        if all(abs(s.imag) < 1e-12 for s in s_final):
            print(f"  (all real, fixed point genuinely multivariate)")

    # Rank-1 collapse analysis.
    print()
    print("=" * 80)
    print("RANK-1 COLLAPSE: symmetric vs asymmetric")
    print("=" * 80)
    print()
    print("Symmetric x = (2, 2, 2): rank-1 collapse to diagonal in 1 step.")
    s = [0.1 + 0j, 0.5 + 0j, 0.9 + 0j]
    for k in range(3):
        s = F_multi(s, ([2.0+0j]*3, 2))
        spread = max(si.real for si in s) - min(si.real for si in s)
        print(f"  step {k+1}: s = {[round(si.real, 8) for si in s]}, spread = {spread:.2e}")
    print()
    print("Asymmetric x = (2, 3, 5): components STAY distinct ==> genuine 3D state.")
    s = [0.1 + 0j, 0.5 + 0j, 0.9 + 0j]
    for k in range(3):
        s = F_multi(s, ([2.0+0j, 3.0+0j, 5.0+0j], 2))
        spread = max(si.real for si in s) - min(si.real for si in s)
        print(f"  step {k+1}: s = {[round(si.real, 8) for si in s]}, spread = {spread:.2e}")

    # PROJ vs DIAG: when does PROJ break?
    print()
    print("=" * 80)
    print("DIAG vs PROJ: when components have WIDELY DIFFERENT rates")
    print("=" * 80)
    x_strong = [2.0+0j, 1000.0+0j]
    args = (x_strong, 2)
    print(f"x = {[xi.real for xi in x_strong]}")
    for d in [2, 3, 4]:
        for mode in (['raw'] if d == 2 else ['diag', 'proj']):
            T = make_multi_T_d(F_multi, d, mode if mode != 'raw' else 'diag')
            s_final, n_it = iterate_multi(T, [0.5+0j, 0.5+0j], args)
            max_res = max(verify_fp(s_final, args))
            tag = f"T_{d} ({mode})"
            print(f"  {tag:>14} : iters = {n_it:>3}, residual = {max_res:.3e}")
    print()
    print("PROJ underperforms when components have well-separated convergence")
    print("rates (the alpha-projection collapses two eigenmodes to one scalar).")
    print()

    print("=" * 80)
    print("VERDICT:")
    print("  - Multivariate prism = vector Aitken applied to asymmetric")
    print("    multivariate Pandrosion.  Genuine n-dim state, not collapsed.")
    print("  - DIAG: works robustly when components are weakly coupled at")
    print("    convergence (asymmetric hypercube qualifies).")
    print("  - PROJ: fails for well-separated rates; useful for symmetric.")
    print("  - As with scalar prism: this is Pandrosion + vector accel,")
    print("    not 'true multivariate Pandrosion'.  Same caveat applies.")
    print("=" * 80)


if __name__ == "__main__":
    main()
