"""
PAPER: 059
TITLE: Total-degree homotopy with 051 universal-Pandrosion corrector
STATUS: homotopy completeness scaffold, Pandrosion-based local correction

MISSION
=======

057 improved root coverage by using complex Bezout-aware multistart, but it
still relied on attraction basins.  This flow changes the game:

    start system G_i(z) = z_i^{d_i} - 1
    all Bezout roots of G are known exactly
    H_t(z) = (1 - t) G(z) + t gamma F(z)
    each path is tracked from t=0 to t=1

The corrector is NOT Newton.  It is the universal 051 geometry portfolio,
ported to complex arithmetic:

    path / cube / simplex / sparse Q_F
    T2 / T3 prism acceleration
    local reanchor

This is the honest "homotopy + 051" version: homotopy supplies the all-root
coverage mechanism; Pandrosion supplies the path corrector.
"""

from __future__ import annotations

import cmath
import itertools
import math
import random
import time
from itertools import product as iprod


# -----------------------------------------------------------------------------
# Complex sparse polynomial primitives.
# -----------------------------------------------------------------------------
def eval_poly(poly, z):
    total = 0.0 + 0.0j
    try:
        for exp, coeff in poly.items():
            term = complex(coeff)
            for zi, ai in zip(z, exp):
                term *= zi ** ai
            total += term
    except (OverflowError, ZeroDivisionError):
        return complex("inf")
    if not (math.isfinite(total.real) and math.isfinite(total.imag)):
        return complex("inf")
    return total


def F_eval(polys, z):
    return [eval_poly(poly, z) for poly in polys]


def is_finite_vec(z):
    return all(math.isfinite(complex(v).real) and math.isfinite(complex(v).imag) for v in z)


def residual_norm(polys, z):
    if not is_finite_vec(z):
        return float("inf")
    vals = F_eval(polys, z)
    if any(not (math.isfinite(v.real) and math.isfinite(v.imag)) for v in vals):
        return float("inf")
    return max(abs(v) for v in vals)


def clean(poly, eps=1e-14):
    return {exp: coeff for exp, coeff in poly.items() if abs(coeff) > eps}


def degree(poly):
    return max((sum(exp) for exp in poly), default=0)


def bezout_estimate(polys):
    out = 1
    for poly in polys:
        out *= max(1, degree(poly))
    return out


def add_scaled(out, poly, scale):
    for exp, coeff in poly.items():
        out[exp] = out.get(exp, 0.0 + 0.0j) + scale * coeff


def start_system(degrees, n):
    polys = []
    zero = tuple([0] * n)
    for i, d in enumerate(degrees):
        exp = tuple(d if j == i else 0 for j in range(n))
        polys.append({exp: 1.0 + 0.0j, zero: -1.0 + 0.0j})
    return polys


def homotopy_polys(target, start, gamma, t):
    polys = []
    for f, g in zip(target, start):
        poly = {}
        add_scaled(poly, g, 1.0 - t)
        add_scaled(poly, f, gamma * t)
        polys.append(clean(poly))
    return polys


# -----------------------------------------------------------------------------
# Linear algebra.
# -----------------------------------------------------------------------------
def solve_linear(A, b):
    n = len(A)
    M = [list(row) + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-13:
            return None
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]
    x = [0.0 + 0.0j] * n
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
# 051 universal Q geometries, complex port.
# -----------------------------------------------------------------------------
def path_slope(polys, anchor, z, order):
    n = len(z)
    Q = [[0.0 + 0.0j] * n for _ in range(n)]
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


def cube_slope(polys, anchor, z, max_orders=12):
    n = len(z)
    perms = list(itertools.permutations(range(n)))
    orders = perms if len(perms) <= max_orders else perms[:max_orders]
    Q_sum = [[0.0 + 0.0j] * n for _ in range(n)]
    cost = 0
    for order in orders:
        Q, c = path_slope(polys, anchor, z, order)
        Q_sum = add_matrix(Q_sum, Q)
        cost += c
    return scale_matrix(Q_sum, 1.0 / len(orders)), cost


def support_activity(polys, n):
    activity = [0.0] * n
    for poly in polys:
        for exp, coeff in poly.items():
            mag = abs(coeff)
            for j, e in enumerate(exp):
                activity[j] += mag * e
    m = max(activity) if activity else 0.0
    return [a if a > 0.0 else max(1.0, 0.1 * m) for a in activity]


def edge_frame(polys, anchor, z, F_anchor):
    n = len(z)
    B = [[0.0 + 0.0j] * n for _ in range(n)]
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
    norm2 = sum(abs(v) ** 2 for v in delta)
    if norm2 < 1e-28:
        return path_slope(polys, anchor, z, tuple(range(n)))
    B, cost = edge_frame(polys, anchor, z, F_anchor)
    F_z = F_eval(polys, z)
    cost += 1
    defect = [F_z[i] - F_anchor[i] - bd for i, bd in enumerate(matvec(B, delta))]
    if mode == "sparse":
        activity = support_activity(polys, n)
        w = [activity[j] * delta[j].conjugate() for j in range(n)]
    else:
        w = [v.conjugate() for v in delta]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        denom = sum(delta[j].conjugate() * delta[j] for j in range(n))
        w = [v.conjugate() for v in delta]
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q, cost


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


def universal_h(polys, anchor, z, F_anchor, mode):
    Q, cost = Q_geometry(polys, anchor, z, F_anchor, mode)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), cost, False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if not is_finite_vec(out):
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
    lam_den = sum(s1[k] - z[k] for k in range(len(z)))
    if abs(lam_den) < 1e-300:
        return t2, c1 + c2 + c3, True, s2
    lam = sum(s2[k] - s1[k] for k in range(len(z))) / lam_den
    if abs(lam - 1.0) < 1e-12:
        return s3, c1 + c2 + c3, True, t2
    out = [t2[k] - (s3[k] - t2[k]) / (lam - 1.0) for k in range(len(z))]
    return out, c1 + c2 + c3, True, t2


def accept_best(polys, z, proposal, fallback, r0):
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
    return (list(best_z), best_r) if best_r <= r0 else (list(z), r0)


def correct_051(polys, start, mode="path", accel="T2", K=2, tol=1e-10, max_cycles=10):
    n = len(start)
    z = list(start)
    anchor = [z[i] + complex(0.17 * (i + 1), -0.11) for i in range(n)]
    evals = 0
    steps = 0
    for cycle in range(max_cycles):
        F_anchor = F_eval(polys, anchor)
        evals += 1
        if max(abs(v) for v in F_anchor) < tol:
            return {"z": anchor, "ok": True, "res": 0.0, "evals": evals, "steps": steps}
        for _ in range(K):
            r0 = residual_norm(polys, z)
            evals += 1
            if r0 < tol:
                return {"z": z, "ok": True, "res": r0, "evals": evals, "steps": steps}
            if accel == "T3":
                proposal, c, ok, fallback = t3_step(polys, anchor, z, F_anchor, mode)
            else:
                proposal, c, ok, fallback = t2_step(polys, anchor, z, F_anchor, mode)
            evals += c
            if not ok:
                return {"z": z, "ok": False, "res": r0, "evals": evals, "steps": steps}
            z, r1 = accept_best(polys, z, proposal, fallback, r0)
            steps += 1
            if r1 < tol:
                return {"z": z, "ok": True, "res": r1, "evals": evals, "steps": steps}
        anchor = list(z)
    r = residual_norm(polys, z)
    return {"z": z, "ok": r < tol, "res": r, "evals": evals, "steps": steps}


def correct_portfolio(polys, start, tol=1e-10):
    modes = ["path", "simplex", "sparse"]
    if len(start) <= 3:
        modes.append("cube")
    attempts = []
    for mode in modes:
        for accel in ["T2", "T3"]:
            result = correct_051(polys, start, mode=mode, accel=accel, tol=tol)
            result["mode"] = f"{mode}/{accel}"
            attempts.append(result)
            if result["ok"]:
                return result
    return min(attempts, key=lambda r: r["res"])


# -----------------------------------------------------------------------------
# Newton corrector baseline inside the same homotopy tracker.
# -----------------------------------------------------------------------------
def fd_jacobian(polys, z, h=1e-7):
    n = len(z)
    f0 = F_eval(polys, z)
    J = [[0.0 + 0.0j] * n for _ in range(n)]
    cost = 1
    for j in range(n):
        zp = list(z)
        zp[j] += h
        fp = F_eval(polys, zp)
        cost += 1
        for i in range(n):
            J[i][j] = (fp[i] - f0[i]) / h
    return J, cost


def correct_newton(polys, start, tol=1e-10, max_iter=10):
    z = list(start)
    evals = 0
    weighted = 0
    n = len(z)
    for it in range(max_iter):
        vals = F_eval(polys, z)
        evals += 1
        r0 = max(abs(v) for v in vals)
        if r0 < tol:
            return {"z": z, "ok": True, "res": r0, "evals": evals, "weighted": weighted + evals, "steps": it}
        J, c = fd_jacobian(polys, z)
        evals += c - 1
        weighted += n
        step = solve_linear(J, [-v for v in vals])
        if step is None:
            return {"z": z, "ok": False, "res": r0, "evals": evals, "weighted": weighted + evals, "steps": it}
        accepted = False
        lam = 1.0
        for _ in range(10):
            cand = [z[i] + lam * step[i] for i in range(n)]
            rc = residual_norm(polys, cand)
            evals += 1
            if math.isfinite(rc) and rc < r0:
                z = cand
                accepted = True
                break
            lam *= 0.5
        if not accepted:
            return {"z": z, "ok": False, "res": r0, "evals": evals, "weighted": weighted + evals, "steps": it}
    r = residual_norm(polys, z)
    return {"z": z, "ok": r < tol, "res": r, "evals": evals, "weighted": weighted + evals, "steps": max_iter}


# -----------------------------------------------------------------------------
# Homotopy tracking.
# -----------------------------------------------------------------------------
def roots_of_unity(d):
    return [cmath.exp(2j * math.pi * k / d) for k in range(d)]


def start_roots(degrees):
    return [list(root) for root in iprod(*[roots_of_unity(d) for d in degrees])]


def track_one_path(target, start, gamma, z0, corrector, tol=1e-9):
    t = 0.0
    dt = 0.003
    z = list(z0)
    prev_z = None
    prev_t = None
    evals = 0
    steps = 0
    failures = 0
    last_mode = ""
    while t < 1.0 - 1e-15:
        trial_dt = min(dt, 1.0 - t)
        t_next = t + trial_dt
        if prev_z is None:
            pred = list(z)
        else:
            slope = [(z[i] - prev_z[i]) / (t - prev_t) for i in range(len(z))]
            pred = [z[i] + trial_dt * slope[i] for i in range(len(z))]
        H = homotopy_polys(target, start, gamma, t_next)
        result = corrector(H, pred, tol=tol)
        evals += result.get("evals", 0)
        if result["ok"] and residual_norm(H, result["z"]) < 10.0 * tol:
            prev_z, prev_t = list(z), t
            z = result["z"]
            t = t_next
            steps += 1
            dt = min(0.01, dt * 1.08)
            last_mode = result.get("mode", last_mode)
            continue
        failures += 1
        dt *= 0.5
        if dt < 1e-5 or failures > 80:
            return {"z": z, "ok": False, "res": residual_norm(target, z), "evals": evals, "steps": steps, "failures": failures, "mode": last_mode}
    r = residual_norm(target, z)
    return {"z": z, "ok": r < 1e-7, "res": r, "evals": evals, "steps": steps, "failures": failures, "mode": last_mode}


def is_new_root(z, roots, tol=1e-5):
    for root in roots:
        if max(abs(z[i] - root[i]) for i in range(len(z))) < tol:
            return False
    return True


def track_all_roots(target, use="051", tol=1e-9):
    n = len(target)
    degrees = [max(1, degree(poly)) for poly in target]
    start = start_system(degrees, n)
    gamma = cmath.exp(0.73j)
    roots0 = start_roots(degrees)
    corrector = correct_portfolio if use == "051" else correct_newton
    roots = []
    total_evals = 0
    ok_paths = 0
    failed_paths = 0
    mode_counts = {}
    for z0 in roots0:
        result = track_one_path(target, start, gamma, z0, corrector, tol=tol)
        total_evals += result["evals"] if use == "051" else result.get("weighted", result["evals"])
        if result["ok"]:
            ok_paths += 1
            if is_new_root(result["z"], roots):
                roots.append(result["z"])
            mode = result.get("mode", "")
            if mode:
                mode_counts[mode] = mode_counts.get(mode, 0) + 1
        else:
            failed_paths += 1
    bez = math.prod(degrees)
    return {
        "roots": roots,
        "bezout": bez,
        "paths": len(roots0),
        "ok_paths": ok_paths,
        "failed_paths": failed_paths,
        "evals": total_evals,
        "coverage": len(roots) / max(bez, 1),
        "modes": mode_counts,
    }


# -----------------------------------------------------------------------------
# Test systems.
# -----------------------------------------------------------------------------
def gen_random_poly_system(n, d, seed):
    rng = random.Random(seed)
    polys = []
    for i in range(n):
        poly = {}
        for alpha in iprod(*[range(d + 1) for _ in range(n)]):
            if sum(alpha) > d:
                continue
            if rng.random() < 0.75:
                poly[tuple(alpha)] = complex(rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0) * 0.15)
        diag = tuple(d if k == i else 0 for k in range(n))
        poly[diag] = poly.get(diag, 0.0) + 1.0
        m = max(abs(v) for v in poly.values())
        polys.append({exp: coeff / m for exp, coeff in poly.items()})
    return polys


def main():
    print("=" * 132)
    print("flow/059 -- Homotopy + 051 universal-Pandrosion corrector (all roots)")
    print("=" * 132)
    print("Start system: G_i=z_i^d-1.  Paths: Bezout-many.  Corrector: 051 path/cube/simplex/sparse T2/T3.")
    print()
    print(f"{'case':>8} {'Bezout':>7} | {'051 roots':>9} {'cov%':>7} {'paths ok/fail':>14} {'051 evals':>10} {'time':>7} | {'Newton roots':>12} {'cov%':>7} {'weighted':>9}")
    print("-" * 132)
    cases = [
        (2, 2),
        (2, 3),
        (2, 4),
        (2, 5),
        (2, 6),
        (3, 2),
        (3, 3),
    ]
    for n, d in cases:
        target = gen_random_poly_system(n, d, seed=59000 + 100 * n + d)
        t0 = time.time()
        p051 = track_all_roots(target, use="051", tol=1e-9)
        elapsed = time.time() - t0
        newt = track_all_roots(target, use="newton", tol=1e-9)
        print(
            f"({n},{d}) {p051['bezout']:>9} | "
            f"{len(p051['roots']):>4}/{p051['bezout']:<4} {100*p051['coverage']:>6.1f}% "
            f"{p051['ok_paths']:>4}/{p051['failed_paths']:<4} {p051['evals']:>10} {elapsed:>6.1f}s | "
            f"{len(newt['roots']):>4}/{newt['bezout']:<4} {100*newt['coverage']:>6.1f}% {newt['evals']:>9}"
        )
        if p051["modes"]:
            modes = " ".join(f"{k}:{v}" for k, v in sorted(p051["modes"].items()))
            print(" " * 20 + "051 modes " + modes)
    print()
    print("=" * 132)
    print("VERDICT")
    print("=" * 132)
    print("  - This is not random multistart coverage: every total-degree start root is tracked.")
    print("  - Completeness now comes from homotopy; the local correction remains Pandrosion/051.")
    print("  - If a path fails, the missing piece is certified adaptive continuation, not more random starts.")
    print("=" * 132)


if __name__ == "__main__":
    main()
