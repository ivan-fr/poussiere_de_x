"""
PAPER: 056
TITLE: Hard KS benchmark for flow/051 universal geometry portfolio
STATUS: targeted test of hypercube/simplex/sparse/path universal geometries

MISSION
=======

The user asked for a hard benchmark of flow/051 specifically:

    path      = one hypercube/Schmidt path
    cube      = averaged hypercube paths
    simplex   = simplex edge frame + exact radial correction
    sparse    = support-weighted universal correction

This flow does NOT call flow/052.  It directly tests the flow/051 portfolio on
harder multivariate Kostlan-Smale style systems, including degree 128.

Important calibration:
    Systems have a planted root r, but each equation is rescaled using probe
    starts so that residuals away from r are O(1).  This prevents the false
    easy case where every start already has residual ~1e-18.

Compared:
    A. flow/051 universal geometries, shown mode-by-mode
    B. best flow/051 portfolio result
    C. finite-difference Newton baseline from same starts
"""

from __future__ import annotations

import importlib.util
import math
import random
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_flow(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


G051 = load_flow("flow051_generalist", "051_generalist_universal_prism_attempt.py")


eval_poly = G051.eval_poly
F_eval = G051.F_eval
residual_norm = G051.residual_norm


def compositions_leq(n, degree):
    out = []

    def rec(prefix, slots, remaining):
        if slots == 1:
            for k in range(remaining + 1):
                out.append(tuple(prefix + [k]))
            return
        for k in range(remaining + 1):
            rec(prefix + [k], slots - 1, remaining - k)

    rec([], n, degree)
    return out


def log_ks_weight(degree, exp):
    total = sum(exp)
    out = math.lgamma(degree + 1) - math.lgamma(degree - total + 1)
    for e in exp:
        out -= math.lgamma(e + 1)
    return 0.5 * out


def ks_exponents(n, degree, mode, sample_size, seed):
    exps = [e for e in compositions_leq(n, degree) if sum(e) > 0]
    if mode == "dense" or len(exps) <= sample_size:
        return exps
    rng = random.Random(seed)
    required = {tuple(degree if i == j else 0 for i in range(n)) for j in range(n)}
    pool = [e for e in exps if e not in required]
    sampled = rng.sample(pool, max(0, sample_size - len(required)))
    return sorted(required.union(sampled))


def probe_starts(root, n):
    starts = [
        [0.5] * n,
        [0.2] * n,
        [0.8] * n,
        [-0.5] * n,
        [root[i] + 0.35 for i in range(n)],
        [root[i] - 0.35 for i in range(n)],
    ]
    rng = random.Random(56000 + n)
    while len(starts) < 10:
        starts.append([rng.uniform(-0.9, 0.9) for _ in range(n)])
    return starts


def calibrate_polys(polys, root):
    probes = probe_starts(root, len(root))
    out = []
    for poly in polys:
        scale = max(abs(eval_poly(poly, p)) for p in probes)
        if scale < 1e-12:
            scale = 1.0
        out.append({e: c / scale for e, c in poly.items()})
    return out


def gen_ks_planted_system(name, n, degree, root, mode="dense", sample_size=1600, seed=0):
    rng = random.Random(seed)
    exps = ks_exponents(n, degree, mode, sample_size, seed)
    max_log = max(log_ks_weight(degree, e) for e in exps)
    zero = tuple([0] * n)
    polys = []
    for i in range(n):
        poly = {}
        for exp in exps:
            # Relative KS weights to avoid overflow, then probe-calibration later.
            weight = math.exp(log_ks_weight(degree, exp) - max_log)
            poly[exp] = rng.gauss(0.0, 1.0) * weight
        pure = tuple(degree if j == i else 0 for j in range(n))
        poly[pure] = poly.get(pure, 0.0) + 1.0
        value_at_root = eval_poly(poly, root)
        poly[zero] = poly.get(zero, 0.0) - value_at_root
        polys.append({e: c for e, c in poly.items() if abs(c) > 1e-14})
    polys = calibrate_polys(polys, root)
    return {
        "name": name,
        "n": n,
        "degree": degree,
        "mode": mode,
        "root": root,
        "polys": polys,
        "terms": sum(len(p) for p in polys),
    }


def hard_starts(root, n, count=8):
    starts = [
        [0.5] * n,
        [0.8] * n,
        [-0.5] * n,
        [root[i] + 0.35 for i in range(n)],
        [root[i] - 0.35 for i in range(n)],
    ]
    rng = random.Random(56100 + n + count)
    while len(starts) < count:
        starts.append([rng.uniform(-0.9, 0.9) for _ in range(n)])
    return starts[:count]


def filtered_hard_starts(polys, root, n, count=3, min_res=1e-3, max_res=1e8):
    starts = []

    def consider(candidate):
        if len(starts) >= count:
            return
        if max(abs(candidate[i] - root[i]) for i in range(n)) < 0.05:
            return
        r = residual_norm(polys, candidate)
        if min_res <= r <= max_res:
            starts.append(candidate)

    for candidate in hard_starts(root, n, count=12):
        consider(candidate)

    rng = random.Random(56200 + n + count)
    attempts = 0
    while len(starts) < count and attempts < 1000:
        candidate = [rng.uniform(-0.95, 0.95) for _ in range(n)]
        consider(candidate)
        attempts += 1

    if len(starts) < count:
        fallback = []
        for candidate in hard_starts(root, n, count=20):
            r = residual_norm(polys, candidate)
            if math.isfinite(r) and r > 0.0:
                fallback.append((r, candidate))
        fallback.sort(reverse=True, key=lambda item: item[0])
        for _, candidate in fallback:
            if len(starts) >= count:
                break
            if candidate not in starts:
                starts.append(candidate)

    return starts[:count]


def run_mode(polys, starts, mode, accel="T2", tol=1e-8, max_cycles=8):
    results = []
    t0 = time.time()
    for start in starts:
        results.append(G051.universal_reanchor(polys, start, mode, accel=accel, tol=tol, max_cycles=max_cycles))
    ok = [r for r in results if r["ok"]]
    best = min(ok, key=lambda r: r["evals"]) if ok else min(results, key=lambda r: r["res"])
    roots = []
    for r in ok:
        z = r["z"]
        if not any(max(abs(z[i] - q[i]) for i in range(len(z))) < 1e-7 for q in roots):
            roots.append(z)
    best = dict(best)
    best["successes"] = len(ok)
    best["attempts"] = len(results)
    best["roots"] = roots
    best["time"] = time.time() - t0
    return best


def run_flow051_portfolio(polys, starts, tol=1e-8):
    modes = ["path", "simplex", "sparse", "cube"]
    rows = {}
    for mode in modes:
        # Cube can be expensive; keep it, because this flow is explicitly about
        # testing the full 051 portfolio.
        t2 = run_mode(polys, starts, mode, accel="T2", tol=tol)
        t3 = run_mode(polys, starts, mode, accel="T3", tol=tol)
        rows[f"{mode}/T2"] = t2
        rows[f"{mode}/T3"] = t3
    ok = [(k, v) for k, v in rows.items() if v["ok"]]
    if ok:
        best_name, best = min(ok, key=lambda kv: kv[1]["evals"])
    else:
        best_name, best = min(rows.items(), key=lambda kv: kv[1]["res"])
    return rows, best_name, best


def run_newton(polys, starts, tol=1e-8, max_iter=25):
    results = []
    t0 = time.time()
    for start in starts:
        results.append(G051.newton_solve(polys, start, tol=tol, max_iter=max_iter))
    ok = [r for r in results if r["ok"]]
    best = min(ok, key=lambda r: r["weighted"]) if ok else min(results, key=lambda r: r["res"])
    best = dict(best)
    best["successes"] = len(ok)
    best["attempts"] = len(results)
    best["time"] = time.time() - t0
    return best


def fmt_mode(name, result):
    return f"{name}:{result['steps']}st/{result['evals']} r={result['res']:.1e} {result['successes']}/{result['attempts']} t={result['time']:.1f}s"


def fmt_newton(result):
    return f"{result['iters']}it/w={result['weighted']} r={result['res']:.1e} {result['successes']}/{result['attempts']} t={result['time']:.1f}s"


def systems():
    return [
        gen_ks_planted_system("KS dense n2 d32", 2, 32, [0.21, -0.27], "dense", seed=5632),
        gen_ks_planted_system("KS dense n2 d128", 2, 128, [0.21, -0.27], "dense", seed=56128),
        gen_ks_planted_system("KS dense n3 d32", 3, 32, [0.19, -0.23, 0.31], "dense", seed=56332),
        gen_ks_planted_system("KS sparse n3 d128", 3, 128, [0.17, -0.22, 0.28], "sparse", sample_size=1400, seed=561283),
        gen_ks_planted_system("KS sparse n4 d64", 4, 64, [0.12, -0.16, 0.21, -0.25], "sparse", sample_size=1200, seed=56464),
    ]


def main():
    print("=" * 132, flush=True)
    print("flow/056 -- HARD KS test of flow/051 universal geometry portfolio", flush=True)
    print("=" * 132, flush=True)
    print("Tests flow/051 only: path / cube / simplex / sparse, each with T2/T3, reanchor, multistart.", flush=True)
    print("KS systems are probe-calibrated so non-root starts have O(1) residual.", flush=True)
    print(flush=True)
    print(
        f"{'system':>20} | {'n':>2} {'d':>4} {'terms':>7} {'start min':>9} {'start max':>9} {'root r':>9} | "
        f"{'best 051':>42} | {'Newton FD':>30}",
        flush=True,
    )
    print("-" * 132, flush=True)

    for system in systems():
        starts = filtered_hard_starts(system["polys"], system["root"], system["n"], count=3)
        start_res = [residual_norm(system["polys"], s) for s in starts]
        start_min = min(start_res)
        start_max = max(start_res)
        root_r = residual_norm(system["polys"], system["root"])
        rows, best_name, best = run_flow051_portfolio(system["polys"], starts)
        newt = run_newton(system["polys"], starts)
        print(
            f"{system['name']:>20} | {system['n']:>2} {system['degree']:>4} {system['terms']:>7} "
            f"{start_min:>9.1e} {start_max:>9.1e} {root_r:>9.1e} | "
            f"{fmt_mode(best_name, best):>42} | {fmt_newton(newt):>30}",
            flush=True,
        )

        # Compact per-geometry detail.
        detail = []
        for key in ["path/T2", "cube/T2", "simplex/T2", "sparse/T2"]:
            r = rows[key]
            detail.append(f"{key}:{'OK' if r['ok'] else 'NO'}:{r['evals']}:{r['res']:.0e}")
        print(" " * 31 + "modes " + " | ".join(detail), flush=True)

    print()
    print("=" * 132)
    print("VERDICT")
    print("=" * 132)
    print("  - This is the requested 051 benchmark: all universal geometries are tested.")
    print("  - The start residual column verifies the systems are no longer trivially solved.")
    print("  - Degree 128 is included for dense n=2 and sparse n=3.")
    print("  - Read the winner honestly: if Newton wins, universal Q_F is working but not")
    print("    automatically better than derivative/Jacobian methods on generic KS.")
    print("=" * 132)


if __name__ == "__main__":
    main()
