"""
PAPER: 047
TITLE: latex/1 universal Pandrosion with simplex slope geometry
STATUS: simplex analogue of flow/046

MISSION
=======

flow/046 built the hypercube version of latex/1:

    Q_cube(a,z) = average of exact Schmidt edge paths through the cube.

This flow builds the SIMPLEX version.  The naive simplex secant matrix from
the n edges

    a -> a + (z_j-a_j)e_j

is not exact for mixed polynomial terms: B(z-a) != F(z)-F(a).  If we used it
directly, regular zeros would not be guaranteed fixed points of

    H_a(z) = a - Q(a,z)^(-1)F(a).

So we use the exact simplex-corrected matrix:

    B_j = [F(a + (z_j-a_j)e_j) - F(a)] / (z_j-a_j)
    r   = F(z) - F(a) - B(z-a)
    Q_simplex = B + r (z-a)^T / ||z-a||^2.

Then

    Q_simplex(z-a) = F(z)-F(a)

exactly.  This is the simplex geometry from latex/1: n simplex edges from the
anchor, plus the unique radial correction needed to hit the terminal vertex.

We apply universal-prism T_2 + reanchor K=3 + multistart on the simplex
Pandrosion systems from flow/027.
"""

from __future__ import annotations

import importlib.util
import math
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_flow(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SIMPLEX027 = load_flow("flow027_simplex", "027_symmetry_reduced_simplex_vs_lairez.py")


def F_simplex(z, x_vec, p):
    return SIMPLEX027.F_target(z, x_vec, p)


def residual_norm(z, x_vec, p):
    if any(not math.isfinite(v) for v in z):
        return float("inf")
    return max(abs(v) for v in F_simplex(z, x_vec, p))


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


def fd_jacobian(anchor, x_vec, p, h=1e-7):
    n = len(anchor)
    f0 = F_simplex(anchor, x_vec, p)
    J = [[0.0] * n for _ in range(n)]
    cost = 1
    for j in range(n):
        zp = list(anchor)
        zp[j] += h
        fp = F_simplex(zp, x_vec, p)
        cost += 1
        for i in range(n):
            J[i][j] = (fp[i] - f0[i]) / h
    return J, cost


def simplex_slope(anchor, z, F_anchor, x_vec, p):
    """Exact simplex-corrected divided-difference matrix."""
    n = len(z)
    delta = [z[i] - anchor[i] for i in range(n)]
    norm2 = sum(v * v for v in delta)
    if norm2 < 1e-28:
        return fd_jacobian(anchor, x_vec, p)

    B = [[0.0] * n for _ in range(n)]
    cost = 0
    for j in range(n):
        edge = list(anchor)
        edge[j] = z[j]
        F_edge = F_simplex(edge, x_vec, p)
        cost += 1
        denom = delta[j]
        if abs(denom) < 1e-300:
            edge[j] = anchor[j] + 1e-7
            F_fd = F_simplex(edge, x_vec, p)
            cost += 1
            for i in range(n):
                B[i][j] = (F_fd[i] - F_anchor[i]) / 1e-7
        else:
            for i in range(n):
                B[i][j] = (F_edge[i] - F_anchor[i]) / denom

    F_z = F_simplex(z, x_vec, p)
    cost += 1
    B_delta = matvec(B, delta)
    defect = [F_z[i] - F_anchor[i] - B_delta[i] for i in range(n)]
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * delta[j] / norm2
    return Q, cost


def universal_h(anchor, z, F_anchor, x_vec, p):
    Q, cost = simplex_slope(anchor, z, F_anchor, x_vec, p)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), cost, False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if any(not math.isfinite(v) for v in out):
        return list(z), cost, False
    return out, cost, True


def t2_step(anchor, z, F_anchor, x_vec, p):
    s1, c1, ok1 = universal_h(anchor, z, F_anchor, x_vec, p)
    if not ok1:
        return list(z), c1, False, list(z)
    s2, c2, ok2 = universal_h(anchor, s1, F_anchor, x_vec, p)
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


def universal_reanchor(start, x_vec, p, K=3, tol=1e-12, max_cycles=50):
    n = len(x_vec)
    z = list(start)
    anchor = [1.0] * n
    evals = 0
    steps = 0
    for cycle in range(max_cycles):
        F_anchor = F_simplex(anchor, x_vec, p)
        evals += 1
        if max(abs(v) for v in F_anchor) < tol:
            return {"z": anchor, "evals": evals, "cycles": cycle, "steps": steps, "ok": True, "res": 0.0}
        for _ in range(K):
            r0 = residual_norm(z, x_vec, p)
            evals += 1
            if r0 < tol:
                return {"z": z, "evals": evals, "cycles": cycle, "steps": steps, "ok": True, "res": r0}
            proposal, c, ok, fallback = t2_step(anchor, z, F_anchor, x_vec, p)
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
    rng = random.Random(47 + n)
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


def multistart(x_vec, p, count=20):
    results = [universal_reanchor(s, x_vec, p) for s in make_starts(len(x_vec), count)]
    ok = [r for r in results if r["ok"]]
    roots = cluster_roots(results)
    best = min(ok, key=lambda r: r["evals"]) if ok else min(results, key=lambda r: r["res"])
    return {"success": len(ok), "total": len(results), "roots": roots, "best": best}


def structural_027(x_vec, p):
    s, evals, iters, qdim, _ = SIMPLEX027.prism_solve_quotient(x_vec, p, d=3)
    r = residual_norm(s, x_vec, p)
    return {"z": s, "evals": evals, "iters": iters, "qdim": qdim, "res": r, "ok": r < 1e-10}


def tag(result):
    return f"{result['steps']}st/{result['evals']} r={result['res']:.0e} {'OK' if result['ok'] else 'NO'}"


def short_vec(z):
    return "[" + ",".join(f"{v:.4g}" for v in z[:4]) + (",..." if len(z) > 4 else "") + "]"


def main():
    print("=" * 118)
    print("flow/047 -- latex/1 universal SIMPLEX Pandrosion")
    print("=" * 118)
    print("Q_simplex = simplex edge secants from anchor + exact radial correction.")
    print("Then apply universal-prism T_2 + reanchor K=3 on flow/027 simplex systems.")
    print()

    cases = [
        ("sym n=3", [2.0, 2.0, 2.0], 2),
        ("partial n=3", [2.0, 2.0, 5.0], 2),
        ("asym n=3", [2.0, 3.0, 5.0], 2),
        ("asym n=4", [2.0, 3.0, 5.0, 7.0], 2),
        ("p3 asym n=3", [2.0, 3.0, 5.0], 3),
    ]

    print(f"{'case':>14} | {'027 structural':>22} | {'simplex universal':>24} | {'multistart':>18}")
    print("-" * 118)
    for label, x_vec, p in cases:
        start = [0.5] * len(x_vec)
        A = structural_027(x_vec, p)
        U = universal_reanchor(start, x_vec, p)
        M = multistart(x_vec, p, count=18)
        Atag = f"q{A['qdim']} {A['iters']}it/{A['evals']} r={A['res']:.0e}"
        Mtag = f"{M['success']}/{M['total']} roots={len(M['roots'])} best={M['best']['evals']}"
        print(f"{label:>14} | {Atag:>22} | {tag(U):>24} | {Mtag:>18}")

    print()
    print("=" * 118)
    print("Root clusters for simplex universal, x=(2,3,5), p=2")
    print("=" * 118)
    M = multistart([2.0, 3.0, 5.0], 2, count=24)
    for i, root in enumerate(M["roots"]):
        print(f"  cluster {i}: {short_vec(root)}, residual={residual_norm(root, [2.0, 3.0, 5.0], 2):.2e}")

    print()
    print("=" * 118)
    print("VERDICT")
    print("=" * 118)
    print("  - This is the latex/1 simplex geometry: Q starts from the simplex")
    print("    secant frame at the anchor and adds the minimal radial correction")
    print("    needed for exact telescoping.")
    print("  - Exact telescoping is the key difference from a naive simplex secant:")
    print("    without it, roots would not be guaranteed fixed points.")
    print("  - On flow/027 simplex systems it works, but the structural prism remains")
    print("    cheaper when the closed-form simplex F_iter is available.")
    print("=" * 118)


if __name__ == "__main__":
    main()
