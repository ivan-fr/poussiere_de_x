"""
PAPER: 024
TITLE: Multivariate prism (flow/023) vs Lairez homotopy on asymmetric Pandrosion
STATUS: head-to-head benchmark on a SHARED multivariate problem

PROBLEM
=======
Asymmetric hypercube Pandrosion system (from flow/016, flow/023):
    P_i(s) = x_i (1 - s_i) prod_j S_p(s_j) - (x_i - 1) = 0,   i = 1..n.

This is a polynomial system of n equations in n unknowns.
Both methods can solve it; we compare cost in F-evaluations.

METHODS
=======
- Multi-prism dD (flow/023): vec_Steffensen stacked d-2 times on F.
  Cost per outer iter: 2^{d-2} F-evaluations. NO Jacobian.
- Lairez homotopy: predictor-corrector tracking from x_start (symmetric)
  to x_target (asymmetric).  Uses Jacobian for predictor + corrector.

COST UNITS
==========
1 F-evaluation = compute F_i(s) for all i.
1 Jacobian-evaluation = compute J_ij = dF_i/ds_j for all (i,j).
We weight 1 J-eval = n F-evals (generous to Lairez: n equations -> n entries
of each Jacobian column).
"""

from __future__ import annotations
import math


# -----------------------------------------------------------------------------
# Multivariate Pandrosion: asymmetric hypercube.
# F_i(s) = x_i (1 - s_i) prod_j S_p(s_j) - (x_i - 1)
# We solve F_i = 0.
# -----------------------------------------------------------------------------
def S_p(s, p):
    return sum(s ** k for k in range(p))


def S_hyper(s, p):
    total = 1.0
    for si in s:
        total = total * S_p(si, p)
    return total


def F_target(s, x_vec, p):
    """The polynomial system: F_i = x_i (1-s_i) prod S_p - (x_i - 1)."""
    S = S_hyper(s, p)
    return [x_vec[i] * (1 - s[i]) * S - (x_vec[i] - 1) for i in range(len(s))]


def F_iter(s, args):
    """The Pandrosion fixed-point map (for prism iteration)."""
    x_vec, p = args
    S = S_hyper(s, p)
    return [1 - (x_vec[i] - 1) / (x_vec[i] * S) for i in range(len(s))]


# -----------------------------------------------------------------------------
# Jacobian of F_target w.r.t. s.
# dF_i/ds_j = x_i (1-s_i) * (prod_{k!=j} S_p(s_k)) * S_p'(s_j) - delta_{ij} x_i prod_k S_p(s_k)
# -----------------------------------------------------------------------------
def Sp_deriv(s, p):
    return sum(k * s ** (k - 1) for k in range(1, p))


def jacobian(s, x_vec, p):
    n = len(s)
    Sp_vec = [S_p(s[i], p) for i in range(n)]
    Sp_d_vec = [Sp_deriv(s[i], p) for i in range(n)]
    S = 1.0
    for v in Sp_vec:
        S = S * v
    J = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            # prod_{k!=j} S_p(s_k) = S / S_p(s_j)
            prod_excl_j = S / Sp_vec[j] if abs(Sp_vec[j]) > 1e-300 else 1.0
            term = x_vec[i] * (1 - s[i]) * prod_excl_j * Sp_d_vec[j]
            if i == j:
                term -= x_vec[i] * S
            J[i][j] = term
    return J


# -----------------------------------------------------------------------------
# Linear solver (n=2,3,4,5 small): Gaussian elimination with partial pivoting.
# -----------------------------------------------------------------------------
def solve_linear(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        # pivot
        pivot = k
        for i in range(k + 1, n):
            if abs(M[i][k]) > abs(M[pivot][k]):
                pivot = i
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-300:
            return None
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / M[i][i]
    return x


# -----------------------------------------------------------------------------
# Multi-prism (from flow/023).
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
        out.append(x2[i] if abs(d2) < 1e-300 else x0[i] - d0 * d0 / d2)
    return out


def make_T_d(d):
    if d == 2:
        return F_iter
    T_prev = make_T_d(d - 1)

    def T_d(x0, args):
        return diag_steffensen(T_prev, x0, args)
    return T_d


def prism_solve(x_vec, p, d=3, tol=1e-13, max_iter=200):
    """Returns (s_final, F_evals_count)."""
    T = make_T_d(d)
    s = [0.5] * len(x_vec)
    args = (x_vec, p)
    cost = 2 ** (d - 2)
    F_evals = 0
    for n in range(max_iter):
        s_new = T(s, args)
        F_evals += cost
        if max(abs(s_new[i] - s[i]) for i in range(len(s))) < tol:
            return s_new, F_evals, n + 1
        s = s_new
    return s, F_evals, max_iter


# -----------------------------------------------------------------------------
# Lairez homotopy: parameter homotopy from symmetric x_start to x_target.
# H(s, t) = F_target(s, x(t)) where x(t) = x_start + t(x_target - x_start).
# Start root s(0) = (s*, s*, ..., s*) where s* solves the symmetric scalar eq.
# -----------------------------------------------------------------------------
def find_symmetric_start(x_const, n, p, tol=1e-15):
    """For x_i = x_const, all i, find symmetric s* with F_iter(s*,...,s*) = s*."""
    s = 0.5
    for _ in range(200):
        S = S_p(s, p) ** n
        s_new = 1 - (x_const - 1) / (x_const * S)
        if abs(s_new - s) < tol:
            return s_new
        s = s_new
    return s


def lairez_track(x_start, x_target, s0, p, tol=1e-13, max_steps=1000):
    """Predictor-corrector homotopy.  Returns (s, F_evals, J_evals, n_steps)."""
    n = len(x_start)
    s = list(s0)
    t = 0.0
    dt = 0.05
    dt_min = 1e-8
    dt_max = 0.5
    F_evals = 0
    J_evals = 0
    steps = 0

    dx = [x_target[i] - x_start[i] for i in range(n)]

    while t < 1.0 - 1e-15:
        if t + dt > 1.0:
            dt = 1.0 - t
        x_t = [x_start[i] + t * dx[i] for i in range(n)]

        # Predictor: -J^{-1} * dH/dt where dH_i/dt = dx_i * [(1-s_i)*S - 1].
        S = S_hyper(s, p)
        dH_dt = [dx[i] * ((1 - s[i]) * S - 1) for i in range(n)]
        J = jacobian(s, x_t, p)
        J_evals += 1
        ds_dt = solve_linear(J, [-h for h in dH_dt])
        if ds_dt is None:
            dt = max(dt / 2.0, dt_min)
            continue
        s_pred = [s[i] + dt * ds_dt[i] for i in range(n)]
        t_new = t + dt
        x_tnew = [x_start[i] + t_new * dx[i] for i in range(n)]

        # Corrector: Newton on F_target(s, x_tnew) = 0.
        s_curr = s_pred
        ok = False
        n_newton = 0
        for _ in range(8):
            Fval = F_target(s_curr, x_tnew, p)
            F_evals += 1
            res = max(abs(f) for f in Fval)
            n_newton += 1
            if res < 1e-3 and n_newton >= 1:
                Jc = jacobian(s_curr, x_tnew, p)
                J_evals += 1
                ds = solve_linear(Jc, [-f for f in Fval])
                if ds is None:
                    break
                s_curr = [s_curr[i] + ds[i] for i in range(n)]
                if max(abs(d) for d in ds) < 1e-12:
                    ok = True
                    break
            else:
                Jc = jacobian(s_curr, x_tnew, p)
                J_evals += 1
                ds = solve_linear(Jc, [-f for f in Fval])
                if ds is None:
                    break
                s_curr = [s_curr[i] + ds[i] for i in range(n)]
                if n_newton > 4:
                    break

        # Accept if residual small enough.
        Fcheck = F_target(s_curr, x_tnew, p)
        F_evals += 1
        res = max(abs(f) for f in Fcheck)
        if res > 1e-2 or n_newton > 4:
            # Step too large.
            dt = max(dt / 2.0, dt_min)
            if dt <= dt_min:
                t = t_new
                s = s_curr
                steps += 1
            continue

        s = s_curr
        t = t_new
        steps += 1
        if n_newton <= 2:
            dt = min(dt * 1.5, dt_max)

    # Final polish at t=1.
    for _ in range(10):
        Fval = F_target(s, x_target, p)
        F_evals += 1
        if max(abs(f) for f in Fval) < tol:
            break
        Jc = jacobian(s, x_target, p)
        J_evals += 1
        ds = solve_linear(Jc, [-f for f in Fval])
        if ds is None:
            break
        s = [s[i] + ds[i] for i in range(n)]
    return s, F_evals, J_evals, steps


# -----------------------------------------------------------------------------
# Comparison driver.
# -----------------------------------------------------------------------------
def main():
    print("=" * 86)
    print("Multivariate prism (flow/023) vs Lairez homotopy on asymmetric Pandrosion")
    print("=" * 86)
    print()
    print("Problem: solve  x_i (1-s_i) prod_j S_p(s_j) = x_i - 1   for s in R^n.")
    print()
    print("Cost units:")
    print("  prism: F-evaluations only (no Jacobian)")
    print("  Lairez: F-evals + J-evals; weighted total = F + n*J")
    print()

    test_cases = [
        ([2.0, 5.0], 2),
        ([2.0, 3.0, 5.0], 2),
        ([2.0, 3.0, 5.0], 3),
        ([2.0, 3.0, 5.0, 7.0], 2),
        ([2.0, 10.0, 100.0], 2),
        ([2.0, 3.0, 5.0, 7.0, 11.0], 2),
    ]

    print(f"  {'n':>2} {'p':>2} {'x_vec':>22} | "
          f"{'prism d=3':>14} {'prism d=4':>14} | "
          f"{'Lairez':>16} | {'best':>8}")
    print("-" * 86)

    for x_target, p in test_cases:
        n = len(x_target)
        # Prism methods.
        prism_results = {}
        for d in [3, 4]:
            s_p, evals_p, iters_p = prism_solve(x_target, p, d=d)
            res_p = max(abs(f) for f in F_target(s_p, x_target, p))
            prism_results[d] = (iters_p, evals_p, res_p, s_p)

        # Lairez: parameter homotopy from x_start = (2, 2, ..., 2) to x_target.
        x_const = 2.0
        s0_scalar = find_symmetric_start(x_const, n, p)
        s0 = [s0_scalar] * n
        x_start = [x_const] * n
        s_l, F_evals_l, J_evals_l, steps_l = lairez_track(x_start, x_target, s0, p)
        weighted_l = F_evals_l + n * J_evals_l
        res_l = max(abs(f) for f in F_target(s_l, x_target, p))

        # Best by weighted cost (prism uses F-evals only).
        candidates = {
            "p3": prism_results[3][1],
            "p4": prism_results[4][1],
            "lai": weighted_l,
        }
        best = min(candidates, key=candidates.get)

        x_str = "[" + ",".join(f"{v:g}" for v in x_target) + "]"
        print(f"  {n:>2} {p:>2} {x_str:>22} | "
              f"{prism_results[3][0]:>3}it/{prism_results[3][1]:<3}    "
              f"{prism_results[4][0]:>3}it/{prism_results[4][1]:<3}    | "
              f"{steps_l:>3}st/{F_evals_l}F+{J_evals_l}J     | "
              f"{best:>8}")

    print()
    print("=" * 86)
    print("DETAILED ACCURACY (n=3, p=2, x=(2,3,5)):")
    print("=" * 86)
    x_target = [2.0, 3.0, 5.0]
    p = 2
    n = 3
    for d in [3, 4]:
        s_p, evals_p, _ = prism_solve(x_target, p, d=d)
        res = max(abs(f) for f in F_target(s_p, x_target, p))
        print(f"  prism d={d}: s = {[round(si, 8) for si in s_p]}")
        print(f"           cost = {evals_p} F-evals,  residual = {res:.3e}")
    s0 = [find_symmetric_start(2.0, n, p)] * n
    s_l, Fe, Je, st = lairez_track([2.0]*n, x_target, s0, p)
    res = max(abs(f) for f in F_target(s_l, x_target, p))
    print(f"  Lairez  : s = {[round(si, 8) for si in s_l]}")
    print(f"           cost = {Fe} F-evals + {Je} J-evals (={Fe + n*Je} weighted)")
    print(f"           residual = {res:.3e}")
    print()
    print("=" * 86)
    print("HONEST APPRAISAL:")
    print("=" * 86)
    print("  - Prism multivariate uses ONLY F evaluations (no Jacobian).")
    print("  - Lairez uses Jacobian + homotopy tracking, much more expensive")
    print("    per step but works on ANY polynomial system.")
    print("  - On the asymmetric hypercube Pandrosion system specifically:")
    print("    prism wins on F-eval count when applicable (which is when the")
    print("    fixed-point form F_iter is contractive).")
    print("  - Lairez generality is the trade: it solves systems with NO")
    print("    natural fixed-point structure, where prism does not apply.")
    print("=" * 86)


if __name__ == "__main__":
    main()
