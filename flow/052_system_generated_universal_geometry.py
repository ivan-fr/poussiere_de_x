"""
PAPER: 052
TITLE: System-generated universal Pandrosion geometry
STATUS: flow/040 idea applied to latex/1 universal Q_F geometries

MISSION
=======

flow/040 generated a closed-form Pandrosion geometry S_A from the input
system.  flow/051 used a fixed portfolio of universal Q_F geometries.

This flow combines both ideas:

    input polynomial system
      -> analyse support/coupling/aniso
      -> generate the most adequate universal geometry Q_F
      -> universal-prism T_2/T_3 + reanchor + multistart

The generated geometry is NOT a fixed catalog pick only.  It synthesises:
    - coordinate order for path/hypercube geometry
    - sparse correction weights from monomial activity
    - coupling blocks from the support interaction graph
    - whether to use path, simplex, sparse, or block geometry

All generated Q geometries still satisfy exact telescoping:

    Q_F(a,z)(z-a) = F(z)-F(a)

so regular zeros remain fixed points.
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


G051 = load_flow("flow051_generalist", "051_generalist_universal_prism_attempt.py")


# Reuse primitives from 051.
eval_poly = G051.eval_poly
F_eval = G051.F_eval
residual_norm = G051.residual_norm
solve_linear = G051.solve_linear
matvec = G051.matvec
fd_jacobian = G051.fd_jacobian
t2_step = None  # local generated version below


def system_dimension(polys):
    return len(next(iter(polys[0])))


def support_stats(polys):
    n = system_dimension(polys)
    activity = [0.0] * n
    interaction = [[0.0] * n for _ in range(n)]
    monomials = 0
    mixed = 0
    max_total = 0

    for poly in polys:
        for exp, coeff in poly.items():
            if abs(coeff) < 1e-12:
                continue
            monomials += 1
            active = [j for j, e in enumerate(exp) if e]
            max_total = max(max_total, sum(exp))
            if len(active) >= 2:
                mixed += 1
            mag = abs(coeff)
            for j in active:
                activity[j] += mag * exp[j]
            for a in active:
                for b in active:
                    if a != b:
                        interaction[a][b] += mag

    if max(activity) == 0:
        activity = [1.0] * n
    else:
        floor = 0.05 * max(activity)
        activity = [v if v > 0 else floor for v in activity]

    density = monomials / max(1, len(polys) * max(1, (max_total + 1) ** n))
    mixed_ratio = mixed / max(1, monomials)
    return {
        "n": n,
        "activity": activity,
        "interaction": interaction,
        "monomials": monomials,
        "max_total": max_total,
        "density": density,
        "mixed_ratio": mixed_ratio,
    }


def connected_blocks(interaction, threshold=0.05):
    n = len(interaction)
    max_w = max([interaction[i][j] for i in range(n) for j in range(n) if i != j] + [0.0])
    if max_w <= 0:
        return [[i] for i in range(n)]
    cutoff = threshold * max_w
    seen = [False] * n
    blocks = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        block = []
        while stack:
            i = stack.pop()
            block.append(i)
            for j in range(n):
                if not seen[j] and (interaction[i][j] + interaction[j][i]) >= cutoff:
                    seen[j] = True
                    stack.append(j)
        blocks.append(sorted(block))
    return sorted(blocks, key=lambda b: (len(b), b[0]), reverse=True)


def synthesize_universal_geometry(polys):
    stats = support_stats(polys)
    n = stats["n"]
    activity = stats["activity"]
    blocks = connected_blocks(stats["interaction"])
    order = sorted(range(n), key=lambda j: (-activity[j], j))
    anisotropy = max(activity) / max(min(activity), 1e-12)

    if len(blocks) > 1 and max(len(b) for b in blocks) < n:
        family = "block"
    elif stats["mixed_ratio"] < 0.20:
        family = "sparse"
    elif anisotropy > 8.0:
        family = "path"
    elif stats["density"] > 0.35 and n <= 2:
        # Full cube symmetrisation is geometrically clean but expensive.
        # Keep it only for very small dense systems.
        family = "cube"
    else:
        family = "simplex"

    return {
        "family": family,
        "order": order,
        "weights": activity,
        "blocks": blocks,
        "stats": stats,
    }


def path_slope_ordered(polys, anchor, z, order):
    return G051.path_slope(polys, anchor, z, tuple(order))


def edge_frame(polys, anchor, z, F_anchor):
    return G051.edge_frame(polys, anchor, z, F_anchor)


def generated_corrected_slope(polys, anchor, z, F_anchor, geom):
    n = len(z)
    delta = [z[i] - anchor[i] for i in range(n)]
    norm2 = sum(v * v for v in delta)
    if norm2 < 1e-28:
        return fd_jacobian(polys, anchor)
    B, cost = edge_frame(polys, anchor, z, F_anchor)
    F_z = F_eval(polys, z)
    cost += 1
    defect = [F_z[i] - F_anchor[i] - bd for i, bd in enumerate(matvec(B, delta))]
    weights = geom["weights"]
    w = [weights[j] * delta[j] for j in range(n)]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        w = list(delta)
        denom = norm2
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q, cost


def block_slope(polys, anchor, z, F_anchor, geom):
    """Generate Q by averaging exact ordered paths inside coupling blocks."""
    n = len(z)
    orders = []
    base_order = []
    for block in geom["blocks"]:
        block_order = sorted(block, key=lambda j: (-geom["weights"][j], j))
        base_order.extend(block_order)
    orders.append(tuple(base_order))
    orders.append(tuple(reversed(base_order)))
    # Add one order with blocks reversed, preserving inside-block order.
    rev_blocks = []
    for block in reversed(geom["blocks"]):
        rev_blocks.extend(sorted(block, key=lambda j: (-geom["weights"][j], j)))
    if tuple(rev_blocks) not in orders:
        orders.append(tuple(rev_blocks))

    Q_sum = [[0.0] * n for _ in range(n)]
    cost = 0
    for order in orders:
        Q, c = path_slope_ordered(polys, anchor, z, order)
        Q_sum = [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(Q_sum, Q)]
        cost += c
    return [[v / len(orders) for v in row] for row in Q_sum], cost


def Q_generated(polys, anchor, z, F_anchor, geom):
    family = geom["family"]
    if family == "path":
        return path_slope_ordered(polys, anchor, z, geom["order"])
    if family == "cube":
        return G051.cube_slope(polys, anchor, z, max_orders=12)
    if family in ("simplex", "sparse"):
        return generated_corrected_slope(polys, anchor, z, F_anchor, geom)
    if family == "block":
        return block_slope(polys, anchor, z, F_anchor, geom)
    raise ValueError(family)


def universal_h_generated(polys, anchor, z, F_anchor, geom):
    Q, cost = Q_generated(polys, anchor, z, F_anchor, geom)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), cost, False
    out = [anchor[i] + step[i] for i in range(len(z))]
    if any(not math.isfinite(v) for v in out):
        return list(z), cost, False
    return out, cost, True


def t2_generated(polys, anchor, z, F_anchor, geom):
    s1, c1, ok1 = universal_h_generated(polys, anchor, z, F_anchor, geom)
    if not ok1:
        return list(z), c1, False, list(z)
    s2, c2, ok2 = universal_h_generated(polys, anchor, s1, F_anchor, geom)
    if not ok2:
        return s1, c1 + c2, True, s1
    out = []
    for a, b, c in zip(z, s1, s2):
        d0 = b - a
        d2 = c - 2.0 * b + a
        out.append(c if abs(d2) < 1e-300 else a - d0 * d0 / d2)
    return out, c1 + c2, True, s1


def t3_generated(polys, anchor, z, F_anchor, geom):
    s1, c1, ok1 = universal_h_generated(polys, anchor, z, F_anchor, geom)
    if not ok1:
        return list(z), c1, False, list(z)
    s2, c2, ok2 = universal_h_generated(polys, anchor, s1, F_anchor, geom)
    if not ok2:
        return s1, c1 + c2, True, s1
    t2, _, _, _ = t2_generated(polys, anchor, z, F_anchor, geom)
    s3, c3, ok3 = universal_h_generated(polys, anchor, t2, F_anchor, geom)
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
    return G051.accept_descent(polys, z, proposal, fallback, r0)


def universal_reanchor_generated(polys, start, geom, accel="T2", K=3, tol=1e-10, max_cycles=50):
    n = len(start)
    z = list(start)
    anchor = [1.0] * n
    if max(abs(anchor[i] - z[i]) for i in range(n)) < 1e-8:
        anchor = [1.0 + 0.13 * (i + 1) for i in range(n)]
    evals = 0
    steps = 0
    for cycle in range(max_cycles):
        F_anchor = F_eval(polys, anchor)
        evals += 1
        if max(abs(v) for v in F_anchor) < tol:
            return {"z": anchor, "evals": evals, "steps": steps, "cycles": cycle, "ok": True, "res": 0.0}
        for _ in range(K):
            r0 = residual_norm(polys, z)
            evals += 1
            if r0 < tol:
                return {"z": z, "evals": evals, "steps": steps, "cycles": cycle, "ok": True, "res": r0}
            if accel == "T3":
                proposal, c, ok, fallback = t3_generated(polys, anchor, z, F_anchor, geom)
            else:
                proposal, c, ok, fallback = t2_generated(polys, anchor, z, F_anchor, geom)
            evals += c
            if not ok:
                return {"z": z, "evals": evals, "steps": steps, "cycles": cycle, "ok": False, "res": r0}
            z, _ = accept_descent(polys, z, proposal, fallback, r0)
            steps += 1
        anchor = list(z)
    r = residual_norm(polys, z)
    return {"z": z, "evals": evals, "steps": steps, "cycles": max_cycles, "ok": r < tol, "res": r}


def make_starts(n, count=10):
    return G051.make_starts(n, count)


def solve_generated(polys, n_starts=10):
    geom = synthesize_universal_geometry(polys)
    results = []
    n = system_dimension(polys)
    for start in make_starts(n, n_starts):
        results.append(universal_reanchor_generated(polys, start, geom, accel="T2"))
        results.append(universal_reanchor_generated(polys, start, geom, accel="T3"))
    ok = [r for r in results if r["ok"]]
    best = min(ok, key=lambda r: r["evals"]) if ok else min(results, key=lambda r: r["res"])
    roots = []
    for r in ok:
        z = r["z"]
        if not any(max(abs(z[i] - q[i]) for i in range(len(z))) < 1e-7 for q in roots):
            roots.append(z)
    best = dict(best)
    best["geom"] = geom
    best["successes"] = len(ok)
    best["attempts"] = len(results)
    best["roots"] = roots
    return best


def main():
    print("=" * 126)
    print("flow/052 -- system-generated universal Pandrosion geometry")
    print("=" * 126)
    print("Input system -> synthesize universal Q geometry -> T2/T3 reanchor multistart.")
    print("Compared against fixed portfolio from flow/051 and Newton FD.")
    print()
    print(f"{'system':>16} | {'generated geom':>20} | {'generated':>30} | {'051 portfolio':>30} | {'Newton':>22}")
    print("-" * 126)

    for system in G051.systems():
        polys = system["polys"]
        gen = solve_generated(polys, n_starts=10)
        port = G051.generalist_solve(polys, n_starts=10)
        newt = G051.newton_multistart(polys, n_starts=10)
        geom = gen["geom"]
        gdesc = f"{geom['family']} ord={geom['order']} b={list(map(len, geom['blocks']))}"
        gtag = f"{gen['steps']}st/{gen['evals']} r={gen['res']:.0e} {gen['successes']}/{gen['attempts']}"
        ptag = f"{port['mode']} {port['steps']}st/{port['evals']} r={port['res']:.0e}"
        ntag = f"{newt['iters']}it/w={newt['weighted']} r={newt['res']:.0e}"
        print(f"{system['name']:>16} | {gdesc:>20} | {gtag:>30} | {ptag:>30} | {ntag:>22}")

    print()
    print("=" * 126)
    print("VERDICT")
    print("=" * 126)
    print("  - The system can generate a universal Pandrosion geometry: order, weights,")
    print("    blocks, and family are inferred from support/coupling.")
    print("  - This reduces brute-force portfolio search, but does not magically create")
    print("    the S_A contraction on generic systems.")
    print("  - If generated geometry beats the portfolio on some cases, it is a useful")
    print("    heuristic.  If not, the honest next step is learning/tuning the synthesis")
    print("    rule, not claiming universal dominance.")
    print("=" * 126)


if __name__ == "__main__":
    main()
