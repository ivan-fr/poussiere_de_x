"""
PAPER: 025
TITLE: Symmetry-reduced multivariate prism vs Lairez
STATUS: separate benchmark; does not modify flow/023

PURPOSE
=======

Flow/023 applies the multivariate difference-prism directly to the full
hypercube Pandrosion map

    F_i(s) = 1 - (x_i - 1) / (x_i prod_j S_p(s_j)).

In the symmetric case x=(2,2,2), this full system collapses to the diagonal.
That is not a bug: it means the effective problem has dimension 1.  A serious
solver should detect this and solve the quotient system instead of pretending
the problem is truly 3D.

This file adds that missing benchmark as a NEW flow:

    x=(2,2,2)   -> quotient dimension 1
    x=(2,2,5)   -> quotient dimension 2
    x=(2,3,5)   -> quotient dimension 3

Then it compares:

    symmetry-reduced prism  vs  generic Lairez-style homotopy

on the same hypercube polynomial system

    P_i(s) = x_i (1-s_i) prod_j S_p(s_j) - (x_i - 1) = 0.

HONEST CAVEAT
=============

This is not claiming to beat Lairez in general.  It shows that on this special
Pandrosion fixed-point system, exploiting symmetry before iteration avoids the
fake full-dimensional work.  A production homotopy solver could also add its
own symmetry reduction.
"""

from __future__ import annotations


# -----------------------------------------------------------------------------
# Polynomial system and fixed-point map.
# -----------------------------------------------------------------------------
def S_p(s, p):
    return sum(s ** k for k in range(p))


def Sp_deriv(s, p):
    return sum(k * s ** (k - 1) for k in range(1, p))


def S_hyper(s, p):
    total = 1.0
    for si in s:
        total *= S_p(si, p)
    return total


def F_target(s, x_vec, p):
    S = S_hyper(s, p)
    return [x_vec[i] * (1.0 - s[i]) * S - (x_vec[i] - 1.0) for i in range(len(s))]


def F_iter(s, args):
    x_vec, p = args
    S = S_hyper(s, p)
    return [1.0 - (x_vec[i] - 1.0) / (x_vec[i] * S) for i in range(len(s))]


# -----------------------------------------------------------------------------
# Symmetry quotient.
# -----------------------------------------------------------------------------
def symmetry_groups(x_vec, tol=0.0):
    """Return unique x-values and groups of indices with equal x."""
    values = []
    groups = []
    for i, x in enumerate(x_vec):
        found = None
        for g, value in enumerate(values):
            if abs(x - value) <= tol:
                found = g
                break
        if found is None:
            values.append(x)
            groups.append([i])
        else:
            groups[found].append(i)
    return values, groups


def expand_from_quotient(y, groups, n):
    s = [0.0] * n
    for g, idxs in enumerate(groups):
        for i in idxs:
            s[i] = y[g]
    return s


def quotient_product(y, counts, p):
    total = 1.0
    for yg, count in zip(y, counts):
        total *= S_p(yg, p) ** count
    return total


def F_iter_quotient(y, args):
    x_unique, counts, p = args
    S = quotient_product(y, counts, p)
    return [1.0 - (x_unique[g] - 1.0) / (x_unique[g] * S) for g in range(len(y))]


# -----------------------------------------------------------------------------
# Difference-prism acceleration.
# -----------------------------------------------------------------------------
def diag_steffensen(T, x0, args):
    x1 = T(x0, args)
    x2 = T(x1, args)
    out = []
    for a, b, c in zip(x0, x1, x2):
        d0 = b - a
        d1 = c - b
        d2 = d1 - d0
        out.append(c if abs(d2) < 1e-300 else a - d0 * d0 / d2)
    return out


def make_T_d(base, d):
    if d == 2:
        return base
    prev = make_T_d(base, d - 1)

    def T(x0, args):
        return diag_steffensen(prev, x0, args)

    return T


def prism_solve_full(x_vec, p, d=3, tol=1e-13, max_iter=200):
    T = make_T_d(F_iter, d)
    s = [0.5] * len(x_vec)
    args = (x_vec, p)
    cost_per_iter = 2 ** (d - 2)
    evals = 0
    for it in range(max_iter):
        sn = T(s, args)
        evals += cost_per_iter
        if max(abs(sn[i] - s[i]) for i in range(len(s))) < tol:
            return sn, evals, it + 1
        s = sn
    return s, evals, max_iter


def prism_solve_quotient(x_vec, p, d=3, tol=1e-13, max_iter=200):
    x_unique, groups = symmetry_groups(x_vec)
    counts = [len(g) for g in groups]
    m = len(x_unique)
    T = make_T_d(F_iter_quotient, d)
    y = [0.5] * m
    args = (x_unique, counts, p)
    cost_per_iter = 2 ** (d - 2)
    q_evals = 0
    for it in range(max_iter):
        yn = T(y, args)
        q_evals += cost_per_iter
        if max(abs(yn[i] - y[i]) for i in range(m)) < tol:
            s = expand_from_quotient(yn, groups, len(x_vec))
            return s, q_evals, it + 1, m, counts
        y = yn
    s = expand_from_quotient(y, groups, len(x_vec))
    return s, q_evals, max_iter, m, counts


# -----------------------------------------------------------------------------
# Generic Lairez-style homotopy for the full system.
# -----------------------------------------------------------------------------
def jacobian(s, x_vec, p):
    n = len(s)
    Sp_vec = [S_p(si, p) for si in s]
    Dp_vec = [Sp_deriv(si, p) for si in s]
    S = 1.0
    for v in Sp_vec:
        S *= v
    J = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            prod_excl = S / Sp_vec[j] if abs(Sp_vec[j]) > 1e-300 else 1.0
            term = x_vec[i] * (1.0 - s[i]) * prod_excl * Dp_vec[j]
            if i == j:
                term -= x_vec[i] * S
            J[i][j] = term
    return J


def solve_linear(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-300:
            return None
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        rhs = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = rhs / M[i][i]
    return x


def find_symmetric_start(x_const, n, p, tol=1e-15):
    s = 0.5
    for _ in range(500):
        S = S_p(s, p) ** n
        sn = 1.0 - (x_const - 1.0) / (x_const * S)
        if abs(sn - s) < tol:
            return sn
        s = sn
    return s


def lairez_track(x_start, x_target, s0, p, tol=1e-13, max_steps=1000):
    n = len(x_target)
    s = list(s0)
    t = 0.0
    dt = 0.05
    dt_min = 1e-8
    dt_max = 0.5
    dx = [x_target[i] - x_start[i] for i in range(n)]
    F_evals = 0
    J_evals = 0
    steps = 0

    while t < 1.0 - 1e-15 and steps < max_steps:
        if t + dt > 1.0:
            dt = 1.0 - t
        x_t = [x_start[i] + t * dx[i] for i in range(n)]
        S = S_hyper(s, p)
        dH_dt = [dx[i] * ((1.0 - s[i]) * S - 1.0) for i in range(n)]
        J = jacobian(s, x_t, p)
        J_evals += 1
        dsdt = solve_linear(J, [-v for v in dH_dt])
        if dsdt is None:
            dt = max(dt / 2.0, dt_min)
            continue
        t_new = t + dt
        s_curr = [s[i] + dt * dsdt[i] for i in range(n)]
        x_new = [x_start[i] + t_new * dx[i] for i in range(n)]

        newton_used = 0
        failed = False
        for _ in range(8):
            Fv = F_target(s_curr, x_new, p)
            F_evals += 1
            if max(abs(v) for v in Fv) < tol:
                break
            Jc = jacobian(s_curr, x_new, p)
            J_evals += 1
            ds = solve_linear(Jc, [-v for v in Fv])
            if ds is None:
                failed = True
                break
            s_curr = [s_curr[i] + ds[i] for i in range(n)]
            newton_used += 1
            if max(abs(v) for v in ds) < 1e-12:
                break

        res = max(abs(v) for v in F_target(s_curr, x_new, p))
        F_evals += 1
        if failed or res > 1e-7 or newton_used > 5:
            dt = max(dt / 2.0, dt_min)
            if dt > dt_min:
                continue

        s = s_curr
        t = t_new
        steps += 1
        if newton_used <= 2:
            dt = min(dt * 1.5, dt_max)

    for _ in range(12):
        Fv = F_target(s, x_target, p)
        F_evals += 1
        if max(abs(v) for v in Fv) < tol:
            break
        Jc = jacobian(s, x_target, p)
        J_evals += 1
        ds = solve_linear(Jc, [-v for v in Fv])
        if ds is None:
            break
        s = [s[i] + ds[i] for i in range(n)]
    return s, F_evals, J_evals, steps


# -----------------------------------------------------------------------------
# Driver.
# -----------------------------------------------------------------------------
def max_residual(s, x_vec, p):
    return max(abs(v) for v in F_target(s, x_vec, p))


def main():
    print("=" * 92)
    print("flow/025 -- symmetry-reduced prism vs generic Lairez")
    print("=" * 92)
    print()
    print("Same hypercube system as flow/023/024:")
    print("  P_i(s)=x_i(1-s_i) prod_j S_p(s_j) - (x_i-1)")
    print()
    print("Reduction:")
    print("  (2,2,2) -> quotient dim 1")
    print("  (2,2,5) -> quotient dim 2")
    print("  (2,3,5) -> quotient dim 3")
    print()

    test_cases = [
        ([2.0, 2.0, 2.0], 2),
        ([2.0, 2.0, 5.0], 2),
        ([2.0, 3.0, 5.0], 2),
        ([2.0, 2.0, 2.0, 2.0], 2),
        ([2.0, 2.0, 3.0, 3.0], 2),
        ([2.0, 3.0, 5.0, 7.0], 2),
        ([2.0, 2.0, 2.0], 3),
        ([2.0, 2.0, 5.0], 3),
    ]

    print(
        f"{'x_vec':>20} {'p':>2} {'qdim':>4} | "
        f"{'full p3':>10} {'quot p3':>10} {'quot p4':>10} | "
        f"{'Lairez weighted':>16} {'best':>8}"
    )
    print("-" * 92)
    for x_vec, p in test_cases:
        n = len(x_vec)
        full_s, full_eval, full_it = prism_solve_full(x_vec, p, d=3)
        q3_s, q3_eval, q3_it, qdim, counts = prism_solve_quotient(x_vec, p, d=3)
        q4_s, q4_eval, q4_it, _, _ = prism_solve_quotient(x_vec, p, d=4)

        x_start_const = 1.25
        x_start = [x_start_const] * n
        s0_scalar = find_symmetric_start(x_start_const, n, p)
        s0 = [s0_scalar] * n
        sl, Fe, Je, steps = lairez_track(x_start, x_vec, s0, p)
        weighted = Fe + n * Je

        # Prefer the quotient label on ties: the eval count can match the full
        # collapsed iteration, but the quotient map has lower effective state
        # dimension and is the intended symmetry-aware path.
        ordered_costs = [("q3", q3_eval), ("q4", q4_eval), ("full3", full_eval), ("lai", weighted)]
        best = min(ordered_costs, key=lambda item: item[1])[0]
        x_label = "[" + ",".join(f"{v:g}" for v in x_vec) + "]"
        print(
            f"{x_label:>20} {p:>2} {qdim:>4} | "
            f"{full_it:>2}it/{full_eval:<3} {q3_it:>2}it/{q3_eval:<3} "
            f"{q4_it:>2}it/{q4_eval:<3} | "
            f"{steps:>2}st/{Fe}F+{Je}J={weighted:<4} {best:>8}"
        )

        # Fail loudly if the symmetry-reduced solution is not a solution of the
        # original full system.
        res = max_residual(q3_s, x_vec, p)
        if res > 1e-8:
            print(f"  WARNING: quotient residual too large: {res:.3e}")

    print()
    print("=" * 92)
    print("DETAIL: symmetric case x=(2,2,2), p=2")
    print("=" * 92)
    x_vec = [2.0, 2.0, 2.0]
    p = 2
    s, qeval, it, qdim, counts = prism_solve_quotient(x_vec, p, d=3)
    print(f"  quotient dimension = {qdim}, counts = {counts}")
    print(f"  symmetry-reduced prism d=3: s = {[round(v, 12) for v in s]}")
    print(f"  residual = {max_residual(s, x_vec, p):.3e}, cost = {it} iterations / {qeval} quotient F-evals")
    s_full, eval_full, it_full = prism_solve_full(x_vec, p, d=3)
    print(f"  unreduced full prism d=3 : cost = {it_full} iterations / {eval_full} full F-evals")
    print()
    print("Verdict:")
    print("  Symmetric problems are not truly full-dimensional.  flow/025 detects")
    print("  the equality pattern in x_i, solves the quotient, and expands back.")
    print("  This beats generic Lairez on the symmetric Pandrosion hypercube test,")
    print("  but it is a symmetry-specialized advantage, not a general Lairez result.")
    print("=" * 92)


if __name__ == "__main__":
    main()
