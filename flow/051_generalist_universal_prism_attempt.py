"""
PAPER: 051
TITLE: Generalist universal-prism attempt (hypercube/simplex/sparse Q_F)
STATUS: honest attempt to go beyond exact Pandrosion S_A systems

MISSION
=======

The user wants the fully general version:

    universal-prism T_2/T_3 + reanchor + multistart
    with latex/1 hypercube/simplex/sparse universal geometries
    on arbitrary polynomial systems.

This flow implements a geometry-portfolio solver for a general sparse
polynomial map F: R^n -> R^n:

    H_a(z) = a - Q_F(a,z)^(-1) F(a)

Candidate universal geometries:
    path      : one Schmidt hypercube path
    cube      : averaged Schmidt hypercube paths
    simplex   : simplex edge frame + exact radial correction
    sparse    : simplex edge frame + exact support-weighted correction

All Q geometries satisfy the exact telescope:

    Q_F(a,z)(z-a) = F(z)-F(a)

so regular zeros remain fixed points.  We then apply T_2/T_3 prism,
K=3 reanchor, and multistart.  The solver picks the best converged attempt.

HONEST EXPECTATION
==================

This is the right generalist Pandrosion-based object, but it cannot be assumed
to beat Lairez everywhere.  On arbitrary polynomial systems the special S_A
closed-form contraction is gone; the method becomes a derivative-free secant
family competing with Newton/homotopy, not a guaranteed replacement.
"""

from __future__ import annotations

import itertools
import math
import random


# -----------------------------------------------------------------------------
# Sparse polynomial utilities.
# -----------------------------------------------------------------------------
def eval_poly(poly, z):
    total = 0.0
    try:
        for exp, coeff in poly.items():
            term = coeff
            for zi, ai in zip(z, exp):
                term *= zi ** ai
            total += term
    except OverflowError:
        return float("inf")
    return total if math.isfinite(total) else float("inf")


def F_eval(polys, z):
    return [eval_poly(p, z) for p in polys]


def residual_norm(polys, z):
    if any(not math.isfinite(v) for v in z):
        return float("inf")
    vals = F_eval(polys, z)
    if any(not math.isfinite(v) for v in vals):
        return float("inf")
    return max(abs(v) for v in vals)


def add_term(poly, exp, coeff):
    poly[exp] = poly.get(exp, 0.0) + coeff


def shift(exp, i):
    out = list(exp)
    out[i] += 1
    return tuple(out)


def clean(poly, eps=1e-12):
    return {exp: coeff for exp, coeff in poly.items() if abs(coeff) > eps}


def monomial_support_activity(polys, n):
    activity = [0.0] * n
    for poly in polys:
        for exp, coeff in poly.items():
            mag = abs(coeff)
            for j, e in enumerate(exp):
                if e:
                    activity[j] += mag * e
    if max(activity) == 0.0:
        return [1.0] * n
    return [a if a > 0.0 else 0.1 * max(activity) for a in activity]


# -----------------------------------------------------------------------------
# Linear algebra.
# -----------------------------------------------------------------------------
def solve_linear(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-14:
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


def matvec(A, x):
    return [sum(row[j] * x[j] for j in range(len(x))) for row in A]


def add_matrix(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def scale_matrix(A, c):
    return [[c * v for v in row] for row in A]


# -----------------------------------------------------------------------------
# Universal Q geometries.
# -----------------------------------------------------------------------------
def path_slope(polys, anchor, z, order):
    n = len(z)
    Q = [[0.0] * n for _ in range(n)]
    cur = list(z)
    F_cur = F_eval(polys, cur)
    cost = 1
    for j in order:
        nxt = list(cur)
        nxt[j] = anchor[j]
        F_next = F_eval(polys, nxt)
        cost += 1
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300:
            zp = list(cur)
            zp[j] += 1e-7
            fp = F_eval(polys, zp)
            cost += 1
            for i in range(n):
                Q[i][j] = (fp[i] - F_cur[i]) / 1e-7
        else:
            for i in range(n):
                Q[i][j] = (F_cur[i] - F_next[i]) / denom
        cur = nxt
        F_cur = F_next
    return Q, cost


def cube_slope(polys, anchor, z, max_orders=24):
    n = len(z)
    perms = list(itertools.permutations(range(n)))
    if len(perms) > max_orders:
        rng = random.Random(5100 + n)
        orders = [tuple(range(n)), tuple(reversed(range(n)))]
        while len(orders) < max_orders:
            p = tuple(rng.sample(range(n), n))
            if p not in orders:
                orders.append(p)
    else:
        orders = perms
    Q_sum = [[0.0] * n for _ in range(n)]
    cost = 0
    for order in orders:
        Q, c = path_slope(polys, anchor, z, order)
        Q_sum = add_matrix(Q_sum, Q)
        cost += c
    return scale_matrix(Q_sum, 1.0 / len(orders)), cost


def edge_frame(polys, anchor, z, F_anchor):
    n = len(z)
    B = [[0.0] * n for _ in range(n)]
    cost = 0
    for j in range(n):
        edge = list(anchor)
        edge[j] = z[j]
        F_edge = F_eval(polys, edge)
        cost += 1
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300:
            edge[j] = anchor[j] + 1e-7
            F_fd = F_eval(polys, edge)
            cost += 1
            for i in range(n):
                B[i][j] = (F_fd[i] - F_anchor[i]) / 1e-7
        else:
            for i in range(n):
                B[i][j] = (F_edge[i] - F_anchor[i]) / denom
    return B, cost


def corrected_edge_slope(polys, anchor, z, F_anchor, mode):
    n = len(z)
    delta = [z[i] - anchor[i] for i in range(n)]
    norm2 = sum(v * v for v in delta)
    if norm2 < 1e-28:
        return fd_jacobian(polys, anchor)
    B, cost = edge_frame(polys, anchor, z, F_anchor)
    F_z = F_eval(polys, z)
    cost += 1
    defect = [F_z[i] - F_anchor[i] - bd for i, bd in enumerate(matvec(B, delta))]

    if mode == "sparse":
        activity = monomial_support_activity(polys, n)
        w = [activity[j] * delta[j] for j in range(n)]
    else:
        w = list(delta)
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        w = list(delta)
        denom = norm2

    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q, cost


def fd_jacobian(polys, z, h=1e-7):
    n = len(z)
    f0 = F_eval(polys, z)
    J = [[0.0] * n for _ in range(n)]
    cost = 1
    for j in range(n):
        zp = list(z)
        zp[j] += h
        fp = F_eval(polys, zp)
        cost += 1
        for i in range(n):
            J[i][j] = (fp[i] - f0[i]) / h
    return J, cost


def Q_geometry(polys, anchor, z, F_anchor, mode):
    if mode == "path":
        return path_slope(polys, anchor, z, tuple(range(len(z))))
    if mode == "cube":
        return cube_slope(polys, anchor, z)
    if mode == "simplex":
        return corrected_edge_slope(polys, anchor, z, F_anchor, "simplex")
    if mode == "sparse":
        return corrected_edge_slope(polys, anchor, z, F_anchor, "sparse")
    raise ValueError(mode)


# -----------------------------------------------------------------------------
# Universal-prism with reanchor.
# -----------------------------------------------------------------------------
def universal_h(polys, anchor, z, F_anchor, mode):
    Q, cost = Q_geometry(polys, anchor, z, F_anchor, mode)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), cost, False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if any(not math.isfinite(v) for v in out):
        return list(z), cost, False
    return out, cost, True


def t2_step(polys, anchor, z, F_anchor, mode):
    s1, c1, ok1 = universal_h(polys, anchor, z, F_anchor, mode)
    if not ok1:
        return list(z), c1, False, list(z)
    s2, c2, ok2 = universal_h(polys, anchor, s1, F_anchor, mode)
    if not ok2:
        return s1, c1 + c2, True, s1
    out = []
    for a, b, c in zip(z, s1, s2):
        d0 = b - a
        d2 = c - 2.0 * b + a
        out.append(c if abs(d2) < 1e-300 else a - d0 * d0 / d2)
    return out, c1 + c2, True, s1


def t3_step(polys, anchor, z, F_anchor, mode):
    s1, c1, ok1 = universal_h(polys, anchor, z, F_anchor, mode)
    if not ok1:
        return list(z), c1, False, list(z)
    s2, c2, ok2 = universal_h(polys, anchor, s1, F_anchor, mode)
    if not ok2:
        return s1, c1 + c2, True, s1
    t2 = []
    for a, b, c in zip(z, s1, s2):
        d0 = b - a
        d2 = c - 2.0 * b + a
        t2.append(c if abs(d2) < 1e-300 else a - d0 * d0 / d2)
    s3, c3, ok3 = universal_h(polys, anchor, t2, F_anchor, mode)
    if not ok3:
        return t2, c1 + c2 + c3, True, s2
    lam_num = sum(s2[k] - s1[k] for k in range(len(z)))
    lam_den = sum(s1[k] - z[k] for k in range(len(z)))
    if abs(lam_den) < 1e-300:
        return t2, c1 + c2 + c3, True, s2
    lam = lam_num / lam_den
    if abs(lam - 1.0) < 1e-12:
        return s3, c1 + c2 + c3, True, t2
    out = [t2[k] - (s3[k] - t2[k]) / (lam - 1.0) for k in range(len(z))]
    return out, c1 + c2 + c3, True, t2


def accept_descent(polys, z, proposal, fallback, r0):
    candidates = []
    for cand in [proposal, fallback]:
        r = residual_norm(polys, cand)
        if math.isfinite(r):
            candidates.append((r, cand))
    for lam in [0.5, 0.25, 0.125, 0.0625]:
        cand = [z[i] + lam * (proposal[i] - z[i]) for i in range(len(z))]
        r = residual_norm(polys, cand)
        if math.isfinite(r):
            candidates.append((r, cand))
    if not candidates:
        return list(z), r0
    best_r, best_z = min(candidates, key=lambda item: item[0])
    return (list(best_z), best_r) if best_r <= r0 or math.isfinite(best_r) else (list(z), r0)


def universal_reanchor(polys, start, mode, accel="T2", K=3, tol=1e-10, max_cycles=50):
    n = len(start)
    z = list(start)
    anchor = [1.0] * n
    # If anchor is too close to start, perturb it deterministically.
    if max(abs(anchor[i] - z[i]) for i in range(n)) < 1e-8:
        anchor = [1.0 + 0.13 * (i + 1) for i in range(n)]
    evals = 0
    steps = 0
    for cycle in range(max_cycles):
        F_anchor = F_eval(polys, anchor)
        evals += 1
        if max(abs(v) for v in F_anchor) < tol:
                return {"z": anchor, "evals": evals, "steps": steps, "cycles": cycle, "ok": True, "res": 0.0, "mode": f"{mode}/{accel}"}
        for _ in range(K):
            r0 = residual_norm(polys, z)
            evals += 1
            if r0 < tol:
                return {"z": z, "evals": evals, "steps": steps, "cycles": cycle, "ok": True, "res": r0, "mode": f"{mode}/{accel}"}
            if accel == "T3":
                proposal, c, ok, fallback = t3_step(polys, anchor, z, F_anchor, mode)
            else:
                proposal, c, ok, fallback = t2_step(polys, anchor, z, F_anchor, mode)
            evals += c
            if not ok:
                return {"z": z, "evals": evals, "steps": steps, "cycles": cycle, "ok": False, "res": r0, "mode": f"{mode}/{accel}"}
            z, _ = accept_descent(polys, z, proposal, fallback, r0)
            steps += 1
        anchor = list(z)
    r = residual_norm(polys, z)
    return {"z": z, "evals": evals, "steps": steps, "cycles": max_cycles, "ok": r < tol, "res": r, "mode": f"{mode}/{accel}"}


def make_starts(n, count=16):
    starts = [[0.5] * n, [0.2] * n, [0.8] * n, [-0.5] * n]
    rng = random.Random(51051 + n)
    while len(starts) < count:
        starts.append([rng.uniform(-1.2, 1.2) for _ in range(n)])
    return starts[:count]


def generalist_solve(polys, n_starts=12):
    modes = ["path", "simplex", "sparse"]
    if len(next(iter(polys[0]))) <= 4:
        modes.append("cube")
    results = []
    for start in make_starts(len(next(iter(polys[0]))), n_starts):
        for mode in modes:
            results.append(universal_reanchor(polys, start, mode, accel="T2"))
            results.append(universal_reanchor(polys, start, mode, accel="T3"))
    ok = [r for r in results if r["ok"]]
    best = min(ok, key=lambda r: r["evals"]) if ok else min(results, key=lambda r: r["res"])
    roots = []
    for r in ok:
        z = r["z"]
        if not any(max(abs(z[i] - q[i]) for i in range(len(z))) < 1e-7 for q in roots):
            roots.append(z)
    best = dict(best)
    best["roots"] = roots
    best["successes"] = len(ok)
    best["attempts"] = len(results)
    return best


# -----------------------------------------------------------------------------
# Newton baseline from same starts.
# -----------------------------------------------------------------------------
def newton_solve(polys, start, tol=1e-10, max_iter=60):
    z = list(start)
    n = len(z)
    F_evals = 0
    J_evals = 0
    for it in range(max_iter):
        r = F_eval(polys, z)
        F_evals += 1
        rn = max(abs(v) for v in r)
        if rn < tol:
            return {"z": z, "F": F_evals, "J": J_evals, "weighted": F_evals + n * J_evals, "iters": it, "ok": True, "res": rn}
        J, c = fd_jacobian(polys, z)
        F_evals += c - 1
        J_evals += 1
        step = solve_linear(J, [-v for v in r])
        if step is None:
            return {"z": z, "F": F_evals, "J": J_evals, "weighted": F_evals + n * J_evals, "iters": it, "ok": False, "res": rn}
        lam = 1.0
        accepted = False
        for _ in range(15):
            cand = [z[i] + lam * step[i] for i in range(n)]
            rc = residual_norm(polys, cand)
            F_evals += 1
            if math.isfinite(rc) and rc < rn:
                z = cand
                accepted = True
                break
            lam *= 0.5
        if not accepted:
            return {"z": z, "F": F_evals, "J": J_evals, "weighted": F_evals + n * J_evals, "iters": it, "ok": False, "res": rn}
    rn = residual_norm(polys, z)
    return {"z": z, "F": F_evals, "J": J_evals, "weighted": F_evals + n * J_evals, "iters": max_iter, "ok": rn < tol, "res": rn}


def newton_multistart(polys, n_starts=12):
    results = [newton_solve(polys, s) for s in make_starts(len(next(iter(polys[0]))), n_starts)]
    ok = [r for r in results if r["ok"]]
    best = min(ok, key=lambda r: r["weighted"]) if ok else min(results, key=lambda r: r["res"])
    best = dict(best)
    best["successes"] = len(ok)
    best["attempts"] = len(results)
    return best


# -----------------------------------------------------------------------------
# Systems.
# -----------------------------------------------------------------------------
def box_support(n, p):
    out = [tuple()]
    for _ in range(n):
        out = [t + (k,) for t in out for k in range(p)]
    return out


def build_pandrosion_system(name, x_vec, A):
    n = len(x_vec)
    zero = tuple([0] * n)
    polys = []
    for i, x in enumerate(x_vec):
        poly = {}
        for alpha in A:
            add_term(poly, alpha, x)
            add_term(poly, shift(alpha, i), -x)
        add_term(poly, zero, -(x - 1.0))
        polys.append(clean(poly))
    return {"name": name, "n": n, "polys": polys, "kind": "Pandrosion-form"}


def random_known_root_system(name, n, degree, root, density=0.55, seed=0):
    rng = random.Random(seed)
    exps = []

    def rec(prefix, slots, remaining):
        if slots == 1:
            for k in range(remaining + 1):
                exp = tuple(prefix + [k])
                if sum(exp) > 0:
                    exps.append(exp)
            return
        for k in range(remaining + 1):
            rec(prefix + [k], slots - 1, remaining - k)

    rec([], n, degree)
    polys = []
    zero = tuple([0] * n)
    for i in range(n):
        poly = {}
        for exp in exps:
            if rng.random() < density or exp == tuple(degree if j == i else 0 for j in range(n)):
                poly[exp] = rng.uniform(-1.0, 1.0)
        # Dominant pure term to avoid too many singular systems.
        poly[tuple(degree if j == i else 0 for j in range(n))] = 1.0
        value_at_root = eval_poly(poly, root)
        poly[zero] = poly.get(zero, 0.0) - value_at_root
        polys.append(clean(poly))
    return {"name": name, "n": n, "polys": polys, "kind": "generic-known-root", "root": root}


def systems():
    return [
        build_pandrosion_system("S_A hyper n3", [2.0, 3.0, 5.0], box_support(3, 2)),
        random_known_root_system("generic n2 d3", 2, 3, [0.31, -0.42], density=0.70, seed=3),
        random_known_root_system("generic n3 d3", 3, 3, [0.25, -0.35, 0.45], density=0.60, seed=7),
        random_known_root_system("dense n2 d5", 2, 5, [0.2, -0.3], density=0.85, seed=11),
        random_known_root_system("sparse n4 d2", 4, 2, [0.2, -0.1, 0.3, -0.25], density=0.30, seed=17),
    ]


def main():
    print("=" * 124)
    print("flow/051 -- generalist universal-prism attempt")
    print("=" * 124)
    print("Portfolio: path/cube/simplex/sparse universal Q_F + T2/T3 prism + K=3 reanchor + multistart.")
    print("Baseline: finite-difference Newton from the same starts.  This is a proxy, not real Lairez.")
    print()
    print(f"{'system':>16} | {'kind':>18} | {'universal portfolio':>34} | {'Newton FD':>30} | {'winner':>10}")
    print("-" * 124)
    wins = 0
    total = 0
    for system in systems():
        total += 1
        polys = system["polys"]
        U = generalist_solve(polys, n_starts=10)
        N = newton_multistart(polys, n_starts=10)
        winner = "Universal" if U["ok"] and (not N["ok"] or U["evals"] < N["weighted"]) else "Newton"
        if winner == "Universal":
            wins += 1
        utag = f"{U['mode']} {U['steps']}st/{U['evals']} r={U['res']:.0e} {U['successes']}/{U['attempts']} roots={len(U['roots'])}"
        ntag = f"{N['iters']}it/w={N['weighted']} r={N['res']:.0e} {N['successes']}/{N['attempts']}"
        print(f"{system['name']:>16} | {system['kind']:>18} | {utag:>34} | {ntag:>30} | {winner:>10}")

    print()
    print("=" * 124)
    print("VERDICT")
    print("=" * 124)
    print(f"  Universal portfolio wins by this cost rule: {wins}/{total}")
    print()
    print("  Honest reading:")
    print("  - This is the generalist Pandrosion-based solver: no S_A closed form is")
    print("    required; only polynomial evaluations are used to build Q_F.")
    print("  - It still does NOT imply 'beats Lairez everywhere'.  Once the system is")
    print("    generic, the exploitable Pandrosion contraction disappears and Newton-like")
    print("    local methods can be cheaper.")
    print("  - The right claim is now precise: Pandrosion has a universal derivative-free")
    print("    extension, but global dominance over Lairez would require a theorem about")
    print("    basins/condition numbers that we do not have.")
    print("=" * 124)


if __name__ == "__main__":
    main()
