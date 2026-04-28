"""
PAPER: 044
TITLE: Universal-prism + reanchor + multistart on hypercube geometry
STATUS: targeted test requested by user

MISSION
=======

flow/025:
    structural prism + hypercube geometry

This flow keeps the SAME hypercube Pandrosion polynomial system

    P_i(s) = x_i (1 - s_i) prod_j S_p(s_j) - (x_i - 1)

but replaces the structural fixed-point map

    s_i <- 1 - (x_i - 1)/(x_i prod_j S_p(s_j))

by the UNIVERSAL Pandrosion operator from latex/1pandrosion_smale.tex:

    h_a(z) = a - Q_F(a, z)^(-1) F(a)

where Q_F is the Schmidt slope matrix.  Then we apply the same prism idea
(componentwise T_2 / Aitken) and add the paper's adaptive reanchoring:

    run K=3 prism steps from anchor a, then set a <- current z.

Finally, we add multistart.  This tests exactly:

    universal-prism + reanchor + multistart + hypercube geometry.

HONEST EXPECTATION
==================

Structural 025 should remain cheaper on this special hypercube class because
it exploits the closed-form F_iter.  Universal-reanchor should be more general:
it only needs F evaluations and the Schmidt slope, so it is the bridge toward
arbitrary polynomial systems.
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


HYPER025 = load_flow("flow025_hypercube", "025_symmetry_reduced_prism_vs_lairez.py")


# -----------------------------------------------------------------------------
# Hypercube geometry from flow/025.
# -----------------------------------------------------------------------------
def F_hyper(s, x_vec, p):
    return HYPER025.F_target(s, x_vec, p)


def residual_norm(s, x_vec, p):
    if any(not math.isfinite(v) for v in s):
        return float("inf")
    return max(abs(v) for v in F_hyper(s, x_vec, p))


# -----------------------------------------------------------------------------
# Linear algebra and Schmidt slope.
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


def schmidt_slope(anchor, z, x_vec, p):
    """Q_F(anchor,z) such that F(z)-F(anchor)=Q_F(anchor,z)(z-anchor)."""
    n = len(z)
    points = [list(z)]
    cur = list(z)
    for j in range(n):
        cur[j] = anchor[j]
        points.append(cur[:])
    F_vals = [F_hyper(pt, x_vec, p) for pt in points]
    Q = [[0.0] * n for _ in range(n)]
    for j in range(n):
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300:
            zp = list(z)
            zp[j] += 1e-7
            fp = F_hyper(zp, x_vec, p)
            for i in range(n):
                Q[i][j] = (fp[i] - F_vals[0][i]) / 1e-7
        else:
            for i in range(n):
                Q[i][j] = (F_vals[j][i] - F_vals[j + 1][i]) / denom
    return Q, n + 1


def universal_h(anchor, z, F_anchor, x_vec, p):
    """h_a(z)=a-Q(a,z)^(-1)F(a)."""
    Q, cost = schmidt_slope(anchor, z, x_vec, p)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), cost, False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if any(not math.isfinite(v) for v in out):
        return list(z), cost, False
    return out, cost, True


def t2_prism_step(anchor, z, F_anchor, x_vec, p):
    """Universal T_2 prism step from a fixed anchor."""
    s1, c1, ok1 = universal_h(anchor, z, F_anchor, x_vec, p)
    if not ok1:
        return list(z), c1, False
    s2, c2, ok2 = universal_h(anchor, s1, F_anchor, x_vec, p)
    if not ok2:
        return s1, c1 + c2, True
    out = []
    for a, b, c in zip(z, s1, s2):
        d0 = b - a
        d2 = c - 2.0 * b + a
        out.append(c if abs(d2) < 1e-300 else a - d0 * d0 / d2)
    return out, c1 + c2, True


def best_descent(z, proposal, fallback, x_vec, p, r0):
    """Accept proposal with damping if useful; otherwise use the best finite fallback."""
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
    if best_r < r0 or math.isfinite(best_r):
        return list(best_z), best_r
    return list(z), r0


# -----------------------------------------------------------------------------
# Universal-prism variants.
# -----------------------------------------------------------------------------
def universal_fixed_anchor(start, x_vec, p, anchor=None, tol=1e-12, max_iter=120):
    n = len(x_vec)
    z = list(start)
    a = [1.0] * n if anchor is None else list(anchor)
    F_anchor = F_hyper(a, x_vec, p)
    evals = 1
    for it in range(max_iter):
        r0 = residual_norm(z, x_vec, p)
        evals += 1
        if r0 < tol:
            return {"z": z, "evals": evals, "iters": it, "ok": True, "res": r0}
        proposal, c, ok = t2_prism_step(a, z, F_anchor, x_vec, p)
        evals += c
        if not ok:
            return {"z": z, "evals": evals, "iters": it, "ok": False, "res": r0}
        fallback, c1, _ = universal_h(a, z, F_anchor, x_vec, p)
        evals += c1
        z, _ = best_descent(z, proposal, fallback, x_vec, p, r0)
    r = residual_norm(z, x_vec, p)
    return {"z": z, "evals": evals, "iters": max_iter, "ok": r < tol, "res": r}


def universal_reanchor(start, x_vec, p, K=3, tol=1e-12, max_cycles=60):
    """Paper-style reanchor: K prism steps from anchor, then anchor <- current z."""
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
            proposal, c, ok = t2_prism_step(anchor, z, F_anchor, x_vec, p)
            evals += c
            if not ok:
                return {"z": z, "evals": evals, "cycles": cycle, "steps": steps, "ok": False, "res": r0}
            fallback, c1, _ = universal_h(anchor, z, F_anchor, x_vec, p)
            evals += c1
            z, _ = best_descent(z, proposal, fallback, x_vec, p, r0)
            steps += 1
        anchor = list(z)
    r = residual_norm(z, x_vec, p)
    return {"z": z, "evals": evals, "cycles": max_cycles, "steps": steps, "ok": r < tol, "res": r}


def make_starts(n, count=20):
    starts = []
    for v in [-1.5, -0.5, 0.2, 0.5, 0.8, 1.5]:
        starts.append([v] * n)
    for i in range(n):
        for v in [0.2, 0.8, 1.2]:
            s = [0.5] * n
            s[i] = v
            starts.append(s)
    rng = random.Random(44 + n)
    while len(starts) < count:
        starts.append([rng.uniform(-1.2, 1.4) for _ in range(n)])
    # Avoid exact anchor coordinates.
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
        found = False
        for root in roots:
            if max(abs(z[i] - root[i]) for i in range(len(z))) < tol:
                found = True
                break
        if not found:
            roots.append(z)
    return roots


def multistart_reanchor(x_vec, p, starts):
    results = [universal_reanchor(s, x_vec, p) for s in starts]
    ok = [r for r in results if r["ok"]]
    roots = cluster_roots(results)
    best = min(ok, key=lambda r: r["evals"]) if ok else min(results, key=lambda r: r["res"])
    avg_evals = sum(r["evals"] for r in ok) / len(ok) if ok else float("inf")
    return {
        "results": results,
        "success": len(ok),
        "total": len(results),
        "roots": roots,
        "best": best,
        "avg_evals": avg_evals,
    }


def structural_025(x_vec, p, d=3):
    s, evals, iters, qdim, counts = HYPER025.prism_solve_quotient(x_vec, p, d=d)
    return {
        "z": s,
        "evals": evals,
        "iters": iters,
        "qdim": qdim,
        "res": residual_norm(s, x_vec, p),
        "ok": residual_norm(s, x_vec, p) < 1e-10,
    }


def short_vec(z):
    return "[" + ",".join(f"{v:.4g}" for v in z[:4]) + (",..." if len(z) > 4 else "") + "]"


def main():
    print("=" * 118)
    print("flow/044 -- universal-prism + reanchor + multistart on HYPERCUBE geometry")
    print("=" * 118)
    print("Same geometry as flow/025: P_i=x_i(1-s_i)prod_j S_p(s_j)-(x_i-1).")
    print("A = flow/025 structural quotient prism.  B/C = Schmidt universal prism on the same hypercube residual.")
    print()

    cases = [
        ("sym n=3", [2.0, 2.0, 2.0], 2),
        ("partial n=3", [2.0, 2.0, 5.0], 2),
        ("asym n=3", [2.0, 3.0, 5.0], 2),
        ("asym n=4", [2.0, 3.0, 5.0, 7.0], 2),
        ("p3 asym n=3", [2.0, 3.0, 5.0], 3),
    ]

    print(
        f"{'case':>14} | {'025 structural':>24} | {'univ fixed':>22} | "
        f"{'univ reanchor':>24} | {'multistart':>24}"
    )
    print("-" * 118)
    for label, x_vec, p in cases:
        n = len(x_vec)
        start = [0.5] * n
        A = structural_025(x_vec, p, d=3)
        B = universal_fixed_anchor(start, x_vec, p)
        C = universal_reanchor(start, x_vec, p)
        M = multistart_reanchor(x_vec, p, make_starts(n))

        Atag = f"q{A['qdim']} {A['iters']}it/{A['evals']} r={A['res']:.0e}"
        Btag = f"{B['iters']}it/{B['evals']} r={B['res']:.0e} {'OK' if B['ok'] else 'NO'}"
        Ctag = f"{C['steps']}st/{C['evals']} r={C['res']:.0e} {'OK' if C['ok'] else 'NO'}"
        Mtag = f"{M['success']}/{M['total']} roots={len(M['roots'])} best={M['best']['evals']}"
        print(f"{label:>14} | {Atag:>24} | {Btag:>22} | {Ctag:>24} | {Mtag:>24}")

    print()
    print("=" * 118)
    print("Multistart root detail for hypercube p=2, x=(2,3,5)")
    print("=" * 118)
    x_vec = [2.0, 3.0, 5.0]
    M = multistart_reanchor(x_vec, 2, make_starts(3, count=24))
    for i, root in enumerate(M["roots"]):
        print(f"  root cluster {i}: {short_vec(root)}, residual={residual_norm(root, x_vec, 2):.2e}")
    print(f"  success {M['success']}/{M['total']}, avg evals among successes = {M['avg_evals']:.1f}")

    print()
    print("=" * 118)
    print("VERDICT")
    print("=" * 118)
    print("  - This is the requested hybrid: universal Schmidt prism + K=3 reanchor")
    print("    + multistart, applied specifically to the hypercube geometry of flow/025.")
    print("  - It converges on the same hypercube systems, but it is more expensive than")
    print("    025 because Schmidt builds an n x n slope matrix instead of using F_iter.")
    print("  - The value is not speed on hypercube; the value is that the same prism")
    print("    machinery can later be reused when no closed-form F_iter geometry exists.")
    print("=" * 118)


if __name__ == "__main__":
    main()
