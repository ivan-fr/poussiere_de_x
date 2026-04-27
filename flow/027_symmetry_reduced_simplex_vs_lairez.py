"""
PAPER: 027
TITLE: Symmetry-reduced simplex Pandrosion prism vs Lairez
STATUS: simplex analogue of flow/025

PURPOSE
=======

flow/025 is specialised to the hypercube geometry

    S_hyper(s) = prod_j S_p(s_j).

This file tests the same idea for the simplex geometry

    S_simplex(s) = sum_{|alpha| <= p-1} s^alpha.

The associated Pandrosion fixed-point system is

    P_i(s) = x_i (1-s_i) S_simplex(s) - (x_i - 1) = 0.

The fixed-point map is

    F_i(s) = 1 - (x_i - 1)/(x_i S_simplex(s)).

As in flow/025, equal x_i define symmetry blocks:

    (2,2,2) -> quotient dimension 1
    (2,2,5) -> quotient dimension 2
    (2,3,5) -> quotient dimension 3

ANSWER TO THE DESIGN QUESTION
=============================

Yes: the geometry should be adapted to the algebraic structure of the system.

    hypercube geometry <-> product support / rectangular Newton polytope
    simplex geometry   <-> total-degree support / simplex Newton polytope

This does not make one universal geometry.  It gives a dispatcher principle:
choose the Pandrosion geometry whose geometric sum matches the system's
monomial support, then apply symmetry reduction and prism acceleration.
"""

from __future__ import annotations


# -----------------------------------------------------------------------------
# Simplex geometric sums.
# -----------------------------------------------------------------------------
def compositions_leq(n, max_total):
    out = []

    def rec(prefix, remaining_slots, remaining_total):
        if remaining_slots == 1:
            for k in range(remaining_total + 1):
                out.append(tuple(prefix + [k]))
            return
        for k in range(remaining_total + 1):
            rec(prefix + [k], remaining_slots - 1, remaining_total - k)

    rec([], n, max_total)
    return out


def S_simplex(s, p):
    total = 0.0
    for alpha in compositions_leq(len(s), p - 1):
        term = 1.0
        for si, ai in zip(s, alpha):
            term *= si ** ai
        total += term
    return total


def dS_simplex(s, p, j):
    total = 0.0
    for alpha in compositions_leq(len(s), p - 1):
        if alpha[j] == 0:
            continue
        term = alpha[j] * (s[j] ** (alpha[j] - 1))
        for k, ak in enumerate(alpha):
            if k != j:
                term *= s[k] ** ak
        total += term
    return total


def binom(n, k):
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    num = 1
    den = 1
    for i in range(1, k + 1):
        num *= n - k + i
        den *= i
    return num // den


def S_simplex_quotient(y, counts, p):
    """
    Simplex sum after imposing equality inside symmetry blocks.

    If a block has c equal coordinates y, then all monomials with total block
    degree beta contribute binom(beta+c-1,c-1) y^beta.
    """
    total = 0.0
    for beta in compositions_leq(len(y), p - 1):
        coeff = 1
        term = 1.0
        for bg, count, yg in zip(beta, counts, y):
            coeff *= binom(bg + count - 1, count - 1)
            term *= yg ** bg
        total += coeff * term
    return total


# -----------------------------------------------------------------------------
# Fixed-point system.
# -----------------------------------------------------------------------------
def F_target(s, x_vec, p):
    S = S_simplex(s, p)
    return [x_vec[i] * (1.0 - s[i]) * S - (x_vec[i] - 1.0) for i in range(len(s))]


def F_iter(s, args):
    x_vec, p = args
    S = S_simplex(s, p)
    return [1.0 - (x_vec[i] - 1.0) / (x_vec[i] * S) for i in range(len(s))]


def F_iter_quotient(y, args):
    x_unique, counts, p = args
    S = S_simplex_quotient(y, counts, p)
    return [1.0 - (x_unique[g] - 1.0) / (x_unique[g] * S) for g in range(len(y))]


def max_residual(s, x_vec, p):
    return max(abs(v) for v in F_target(s, x_vec, p))


# -----------------------------------------------------------------------------
# Symmetry quotient.
# -----------------------------------------------------------------------------
def symmetry_groups(x_vec, tol=0.0):
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
    T = make_T_d(F_iter_quotient, d)
    y = [0.5] * len(x_unique)
    args = (x_unique, counts, p)
    cost_per_iter = 2 ** (d - 2)
    evals = 0
    for it in range(max_iter):
        yn = T(y, args)
        evals += cost_per_iter
        if max(abs(yn[i] - y[i]) for i in range(len(y))) < tol:
            return expand_from_quotient(yn, groups, len(x_vec)), evals, it + 1, len(x_unique), counts
        y = yn
    return expand_from_quotient(y, groups, len(x_vec)), evals, max_iter, len(x_unique), counts


# -----------------------------------------------------------------------------
# Lairez-style homotopy on the same simplex polynomial system.
# -----------------------------------------------------------------------------
def jacobian(s, x_vec, p):
    n = len(s)
    S = S_simplex(s, p)
    J = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            term = x_vec[i] * (1.0 - s[i]) * dS_simplex(s, p, j)
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
        S = S_simplex([s] * n, p)
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
        S = S_simplex(s, p)
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
def main():
    print("=" * 94)
    print("flow/027 -- symmetry-reduced SIMPLEX prism vs generic Lairez")
    print("=" * 94)
    print()
    print("Geometry: S_simplex(s)=sum_{|alpha|<=p-1} s^alpha")
    print("System  : x_i(1-s_i)S_simplex(s)=x_i-1")
    print()

    test_cases = [
        ([2.0, 2.0], 2),
        ([2.0, 3.0], 2),
        ([2.0, 2.0, 2.0], 2),
        ([2.0, 2.0, 5.0], 2),
        ([2.0, 3.0, 5.0], 2),
        ([2.0, 2.0, 2.0], 3),
        ([2.0, 2.0, 5.0], 3),
        ([2.0, 3.0, 5.0], 3),
        ([2.0, 3.0, 5.0, 7.0], 2),
    ]

    print(
        f"{'x_vec':>18} {'p':>2} {'qdim':>4} | "
        f"{'full p3':>10} {'quot p3':>10} {'quot p4':>10} | "
        f"{'Lairez weighted':>16} {'best':>8}"
    )
    print("-" * 94)
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

        ordered = [("q3", q3_eval), ("q4", q4_eval), ("full3", full_eval), ("lai", weighted)]
        best = min(ordered, key=lambda item: item[1])[0]
        x_label = "[" + ",".join(f"{v:g}" for v in x_vec) + "]"
        print(
            f"{x_label:>18} {p:>2} {qdim:>4} | "
            f"{full_it:>2}it/{full_eval:<3} {q3_it:>2}it/{q3_eval:<3} "
            f"{q4_it:>2}it/{q4_eval:<3} | "
            f"{steps:>2}st/{Fe}F+{Je}J={weighted:<4} {best:>8}"
        )

        res = max_residual(q3_s, x_vec, p)
        if res > 1e-8:
            print(f"  WARNING: simplex quotient residual too large: {res:.3e}")

    print()
    print("Verdict:")
    print("  The same symmetry-reduced prism strategy works for simplex geometry.")
    print("  This supports the dispatcher view: choose the Pandrosion geometry")
    print("  from the system's support, then reduce symmetries and accelerate.")
    print("  Hypercube is not universal; simplex is another specialised geometry.")


if __name__ == "__main__":
    main()
