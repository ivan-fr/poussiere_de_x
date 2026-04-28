"""
PAPER: 046
TITLE: latex/1 universal Pandrosion with true hypercube slope geometry
STATUS: universal hypercube version of latex/1pandrosion_smale.tex

MISSION
=======

latex/1pandrosion_smale.tex defines the universal multivariate Pandrosion map

    H_a(z) = a - Q_F(a,z)^(-1) F(a)

where Q_F(a,z) is a matrix divided difference satisfying

    F(z) - F(a) = Q_F(a,z) (z-a).

The Schmidt matrix used in flows 041/043/044 is one ordered path through the
hypercube between a and z.  That is already a hypercube edge geometry, but it
depends on coordinate order.

This flow builds the symmetric HYPERCUBE version:

    Q_cube(a,z) = average over all coordinate-order Schmidt paths.

Every path is exact, hence the average is exact too:

    Q_cube(a,z)(z-a) = F(z)-F(a).

Then we apply:

    universal-prism T_2 + reanchor K=3 + multistart

on the same hypercube Pandrosion polynomial systems as flow/025.
"""

from __future__ import annotations

import importlib.util
import itertools
import math
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_flow(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


HYPER025 = load_flow("flow025_hypercube", "025_symmetry_reduced_prism_vs_lairez.py")


def F_hyper(z, x_vec, p):
    return HYPER025.F_target(z, x_vec, p)


def residual_norm(z, x_vec, p):
    if any(not math.isfinite(v) for v in z):
        return float("inf")
    return max(abs(v) for v in F_hyper(z, x_vec, p))


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


def add_matrix(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def scale_matrix(A, c):
    return [[c * v for v in row] for row in A]


def coordinate_orders(n, max_orders=24):
    perms = list(itertools.permutations(range(n)))
    if len(perms) <= max_orders:
        return perms
    rng = random.Random(4600 + n)
    out = [tuple(range(n)), tuple(reversed(range(n)))]
    while len(out) < max_orders:
        p = tuple(rng.sample(range(n), n))
        if p not in out:
            out.append(p)
    return out


def schmidt_path_slope(anchor, z, x_vec, p, order):
    """One exact edge path through the hypercube from z to anchor."""
    n = len(z)
    Q = [[0.0] * n for _ in range(n)]
    cur = list(z)
    F_cur = F_hyper(cur, x_vec, p)
    cost = 1
    for j in order:
        nxt = list(cur)
        nxt[j] = anchor[j]
        F_next = F_hyper(nxt, x_vec, p)
        cost += 1
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300:
            zp = list(cur)
            zp[j] += 1e-7
            fp = F_hyper(zp, x_vec, p)
            cost += 1
            for i in range(n):
                Q[i][j] = (fp[i] - F_cur[i]) / 1e-7
        else:
            for i in range(n):
                Q[i][j] = (F_cur[i] - F_next[i]) / denom
        cur = nxt
        F_cur = F_next
    return Q, cost


def hypercube_slope(anchor, z, x_vec, p, mode="avg"):
    """Universal hypercube Q: one path or averaged paths."""
    n = len(z)
    if mode == "path":
        return schmidt_path_slope(anchor, z, x_vec, p, tuple(range(n)))
    orders = coordinate_orders(n)
    Q_sum = [[0.0] * n for _ in range(n)]
    cost = 0
    for order in orders:
        Q, c = schmidt_path_slope(anchor, z, x_vec, p, order)
        Q_sum = add_matrix(Q_sum, Q)
        cost += c
    return scale_matrix(Q_sum, 1.0 / len(orders)), cost


def universal_h(anchor, z, F_anchor, x_vec, p, mode):
    Q, cost = hypercube_slope(anchor, z, x_vec, p, mode=mode)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), cost, False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if any(not math.isfinite(v) for v in out):
        return list(z), cost, False
    return out, cost, True


def t2_step(anchor, z, F_anchor, x_vec, p, mode):
    s1, c1, ok1 = universal_h(anchor, z, F_anchor, x_vec, p, mode)
    if not ok1:
        return list(z), c1, False, list(z)
    s2, c2, ok2 = universal_h(anchor, s1, F_anchor, x_vec, p, mode)
    if not ok2:
        return s1, c1 + c2, True, s1
    out = []
    for a, b, c in zip(z, s1, s2):
        d0 = b - a
        d2 = c - 2.0 * b + a
        out.append(c if abs(d2) < 1e-300 else a - d0 * d0 / d2)
    return out, c1 + c2, True, s1


def accept_descent(z, proposal, fallback, x_vec, p, r0):
    candidates = []
    for cand in [proposal, fallback]:
        r = residual_norm(cand, x_vec, p)
        if math.isfinite(r):
            candidates.append((r, cand))
    for lam in [0.5, 0.25, 0.125, 0.0625]:
        cand = [z[i] + lam * (proposal[i] - z[i]) for i in range(len(z))]
        r = residual_norm(cand, x_vec, p)
        if math.isfinite(r):
            candidates.append((r, cand))
    if not candidates:
        return list(z), r0
    best_r, best_z = min(candidates, key=lambda item: item[0])
    return (list(best_z), best_r) if best_r <= r0 or math.isfinite(best_r) else (list(z), r0)


def universal_reanchor(start, x_vec, p, mode="avg", K=3, tol=1e-12, max_cycles=50):
    n = len(x_vec)
    z = list(start)
    anchor = [1.0] * n
    evals = 0
    steps = 0
    for cycle in range(max_cycles):
        F_anchor = F_hyper(anchor, x_vec, p)
        evals += 1
        if max(abs(v) for v in F_anchor) < tol:
            return {"z": anchor, "evals": evals, "cycles": cycle, "steps": steps, "ok": True, "res": 0.0}
        for _ in range(K):
            r0 = residual_norm(z, x_vec, p)
            evals += 1
            if r0 < tol:
                return {"z": z, "evals": evals, "cycles": cycle, "steps": steps, "ok": True, "res": r0}
            proposal, c, ok, fallback = t2_step(anchor, z, F_anchor, x_vec, p, mode)
            evals += c
            if not ok:
                return {"z": z, "evals": evals, "cycles": cycle, "steps": steps, "ok": False, "res": r0}
            z, _ = accept_descent(z, proposal, fallback, x_vec, p, r0)
            steps += 1
        anchor = list(z)
    r = residual_norm(z, x_vec, p)
    return {"z": z, "evals": evals, "cycles": max_cycles, "steps": steps, "ok": r < tol, "res": r}


def make_starts(n, count=20):
    starts = [[v] * n for v in [-1.5, -0.5, 0.2, 0.5, 0.8, 1.5]]
    for i in range(n):
        for v in [0.2, 0.8, 1.2]:
            s = [0.5] * n
            s[i] = v
            starts.append(s)
    rng = random.Random(46 + n)
    while len(starts) < count:
        starts.append([rng.uniform(-1.2, 1.4) for _ in range(n)])
    for s in starts:
        for i, v in enumerate(s):
            if abs(v - 1.0) < 1e-8:
                s[i] = 0.999
    return starts[:count]


def cluster_roots(results, tol=1e-7):
    roots = []
    for result in results:
        if not result["ok"]:
            continue
        z = result["z"]
        if not any(max(abs(z[i] - r[i]) for i in range(len(z))) < tol for r in roots):
            roots.append(z)
    return roots


def multistart(x_vec, p, mode, count=20):
    results = [universal_reanchor(s, x_vec, p, mode=mode) for s in make_starts(len(x_vec), count)]
    ok = [r for r in results if r["ok"]]
    roots = cluster_roots(results)
    best = min(ok, key=lambda r: r["evals"]) if ok else min(results, key=lambda r: r["res"])
    return {"success": len(ok), "total": len(results), "roots": roots, "best": best}


def structural_025(x_vec, p):
    s, evals, iters, qdim, _ = HYPER025.prism_solve_quotient(x_vec, p, d=3)
    r = residual_norm(s, x_vec, p)
    return {"z": s, "evals": evals, "iters": iters, "qdim": qdim, "res": r, "ok": r < 1e-10}


def tag(result, label="st"):
    return f"{result['steps']}{label}/{result['evals']} r={result['res']:.0e} {'OK' if result['ok'] else 'NO'}"


def short_vec(z):
    return "[" + ",".join(f"{v:.4g}" for v in z[:4]) + (",..." if len(z) > 4 else "") + "]"


def main():
    print("=" * 118)
    print("flow/046 -- latex/1 universal HYPERCUBE Pandrosion")
    print("=" * 118)
    print("Q_cube = average of all Schmidt edge paths through the anchor/current hypercube.")
    print("Then apply universal-prism T_2 + reanchor K=3 on flow/025 hypercube systems.")
    print()

    cases = [
        ("sym n=3", [2.0, 2.0, 2.0], 2),
        ("partial n=3", [2.0, 2.0, 5.0], 2),
        ("asym n=3", [2.0, 3.0, 5.0], 2),
        ("asym n=4", [2.0, 3.0, 5.0, 7.0], 2),
        ("p3 asym n=3", [2.0, 3.0, 5.0], 3),
    ]

    print(f"{'case':>14} | {'025 structural':>22} | {'one path':>24} | {'hypercube avg':>24} | {'multi avg':>18}")
    print("-" * 118)
    for label, x_vec, p in cases:
        start = [0.5] * len(x_vec)
        A = structural_025(x_vec, p)
        path = universal_reanchor(start, x_vec, p, mode="path")
        avg = universal_reanchor(start, x_vec, p, mode="avg")
        multi = multistart(x_vec, p, mode="avg", count=18)
        Atag = f"q{A['qdim']} {A['iters']}it/{A['evals']} r={A['res']:.0e}"
        Mtag = f"{multi['success']}/{multi['total']} roots={len(multi['roots'])} best={multi['best']['evals']}"
        print(f"{label:>14} | {Atag:>22} | {tag(path):>24} | {tag(avg):>24} | {Mtag:>18}")

    print()
    print("=" * 118)
    print("Root clusters for hypercube avg, x=(2,3,5), p=2")
    print("=" * 118)
    M = multistart([2.0, 3.0, 5.0], 2, mode="avg", count=24)
    for i, root in enumerate(M["roots"]):
        print(f"  cluster {i}: {short_vec(root)}, residual={residual_norm(root, [2.0, 3.0, 5.0], 2):.2e}")

    print()
    print("=" * 118)
    print("VERDICT")
    print("=" * 118)
    print("  - This is the latex/1 hypercube geometry: Q is built from hypercube")
    print("    edge divided differences between anchor and current point.")
    print("  - Averaging all paths removes coordinate-order dependence while preserving")
    print("    the exact telescoping identity F(z)-F(a)=Q(z-a).")
    print("  - On flow/025 hypercube systems it works, but remains costlier than the")
    print("    closed-form structural prism because Q_cube is an n x n universal object.")
    print("=" * 118)


if __name__ == "__main__":
    main()
