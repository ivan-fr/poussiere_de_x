"""
PAPER: 054
TITLE: Hard multivariate Kostlan-Smale benchmark for universal Pandrosion geometry
STATUS: stress test after 052; much harder than toy generic systems

MISSION
=======

The previous "generic" systems were too easy: degrees 2,3,5 and small support.
This flow raises the bar with multivariate Kostlan-Smale style systems:

    n = 2, degree 128, dense total-degree support
    n = 3, degree 32, dense total-degree support
    n = 4, degree 16, dense total-degree support
    n = 3, degree 128, sparse sampled KS support
    n = 5, degree 32, sparse sampled KS support

Each system has a planted regular root r by subtracting F(r) from the constant
term.  This keeps the systems high-degree/high-support while giving an
objective success criterion.

Compared methods:
    A. generated universal Pandrosion geometry from flow/052
    B. fixed portfolio universal Pandrosion from flow/051
    C. finite-difference Newton from same starts

This is intentionally not a friendly benchmark.  If universal Pandrosion loses
on these, that is useful information.
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
G052 = load_flow("flow052_generated", "052_system_generated_universal_geometry.py")


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
    all_exps = compositions_leq(n, degree)
    all_exps = [e for e in all_exps if sum(e) > 0]
    if mode == "dense" or len(all_exps) <= sample_size:
        return all_exps
    rng = random.Random(seed)
    required = {tuple(degree if i == j else 0 for i in range(n)) for j in range(n)}
    pool = [e for e in all_exps if e not in required]
    sampled = rng.sample(pool, max(0, sample_size - len(required)))
    return sorted(required.union(sampled))


def gen_ks_planted_system(name, n, degree, root, mode="dense", sample_size=2000, seed=0):
    rng = random.Random(seed)
    exps = ks_exponents(n, degree, mode, sample_size, seed)
    max_log = max(log_ks_weight(degree, e) for e in exps)
    zero = tuple([0] * n)
    polys = []
    for i in range(n):
        poly = {}
        for exp in exps:
            weight = math.exp(log_ks_weight(degree, exp) - max_log)
            coeff = rng.gauss(0.0, 1.0) * weight
            poly[exp] = coeff
        # Ensure each equation has its pure leading direction.
        pure = tuple(degree if j == i else 0 for j in range(n))
        poly[pure] = poly.get(pure, 0.0) + 1.0
        value = eval_poly(poly, root)
        poly[zero] = poly.get(zero, 0.0) - value
        polys.append({e: c for e, c in poly.items() if abs(c) > 1e-14})
    return {
        "name": name,
        "n": n,
        "degree": degree,
        "mode": mode,
        "root": root,
        "polys": polys,
        "monomials": sum(len(p) for p in polys),
    }


def hard_starts(root, n, count):
    starts = [
        [0.5] * n,
        [0.2] * n,
        [0.8] * n,
        [-0.5] * n,
        [root[i] + 0.35 for i in range(n)],
        [root[i] - 0.35 for i in range(n)],
    ]
    rng = random.Random(54000 + n + count)
    while len(starts) < count:
        starts.append([rng.uniform(-0.9, 0.9) for _ in range(n)])
    return starts[:count]


def generated_budget_solve(system, starts, tol=1e-8, max_cycles=8):
    polys = system["polys"]
    geom = G052.synthesize_universal_geometry(polys)
    results = []
    t0 = time.time()
    for start in starts:
        results.append(G052.universal_reanchor_generated(polys, start, geom, accel="T2", tol=tol, max_cycles=max_cycles))
        results.append(G052.universal_reanchor_generated(polys, start, geom, accel="T3", tol=tol, max_cycles=max_cycles))
    ok = [r for r in results if r["ok"]]
    best = min(ok, key=lambda r: r["evals"]) if ok else min(results, key=lambda r: r["res"])
    best = dict(best)
    best["geom"] = geom
    best["successes"] = len(ok)
    best["attempts"] = len(results)
    best["time"] = time.time() - t0
    return best


def portfolio_budget_solve(system, starts, tol=1e-8, max_cycles=6):
    polys = system["polys"]
    # Cube averaging is intentionally skipped here: on high-degree KS systems
    # it is too expensive for a stress benchmark and was already isolated in 046.
    modes = ["path", "simplex", "sparse"]
    results = []
    t0 = time.time()
    for start in starts:
        for mode in modes:
            results.append(G051.universal_reanchor(polys, start, mode, accel="T2", tol=tol, max_cycles=max_cycles))
            results.append(G051.universal_reanchor(polys, start, mode, accel="T3", tol=tol, max_cycles=max_cycles))
    ok = [r for r in results if r["ok"]]
    best = min(ok, key=lambda r: r["evals"]) if ok else min(results, key=lambda r: r["res"])
    best = dict(best)
    best["successes"] = len(ok)
    best["attempts"] = len(results)
    best["time"] = time.time() - t0
    return best


def newton_budget_solve(system, starts, tol=1e-8, max_iter=20):
    polys = system["polys"]
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


def known_root_residual(system):
    return residual_norm(system["polys"], system["root"])


def systems():
    return [
        gen_ks_planted_system("KS dense n2 d128", 2, 128, [0.21, -0.27], "dense", seed=1282),
        gen_ks_planted_system("KS dense n3 d32", 3, 32, [0.19, -0.23, 0.31], "dense", seed=3203),
        gen_ks_planted_system("KS dense n4 d16", 4, 16, [0.15, -0.18, 0.24, -0.29], "dense", seed=1604),
        gen_ks_planted_system("KS sparse n3 d128", 3, 128, [0.17, -0.22, 0.28], "sparse", sample_size=2200, seed=1283),
        gen_ks_planted_system("KS sparse n5 d32", 5, 32, [0.12, -0.16, 0.21, -0.25, 0.29], "sparse", sample_size=1800, seed=3205),
    ]


def fmt_pand(result):
    geom = result.get("geom")
    if geom:
        mode = geom["family"]
    else:
        mode = result.get("mode", "?")
    return f"{mode} {result['steps']}st/{result['evals']} r={result['res']:.1e} {result['successes']}/{result['attempts']} t={result['time']:.1f}s"


def fmt_newton(result):
    return f"{result['iters']}it/w={result['weighted']} r={result['res']:.1e} {result['successes']}/{result['attempts']} t={result['time']:.1f}s"


def main():
    print("=" * 132, flush=True)
    print("flow/054 -- HARD multivariate KS benchmark (degrees up to 128)", flush=True)
    print("=" * 132, flush=True)
    print("KS-style affine systems with planted root.  tol=1e-8, deliberately tight budget, starts are not exact roots.", flush=True)
    print("Budget: generated max_cycles=8, portfolio max_cycles=6, Newton max_iter=20, starts=3.", flush=True)
    print(flush=True)
    print(
        f"{'system':>20} | {'n':>2} {'d':>4} {'terms':>7} {'root res':>9} | "
        f"{'generated geom':>34} | {'portfolio':>34} | {'Newton FD':>28}"
    , flush=True)
    print("-" * 132, flush=True)

    for system in systems():
        starts = hard_starts(system["root"], system["n"], count=3)
        gen = generated_budget_solve(system, starts)
        port = portfolio_budget_solve(system, starts)
        newt = newton_budget_solve(system, starts)
        print(
            f"{system['name']:>20} | {system['n']:>2} {system['degree']:>4} {system['monomials']:>7} "
            f"{known_root_residual(system):>9.1e} | {fmt_pand(gen):>34} | {fmt_pand(port):>34} | {fmt_newton(newt):>28}"
        , flush=True)

    print(flush=True)
    print("=" * 132, flush=True)
    print("VERDICT", flush=True)
    print("=" * 132, flush=True)
    print("  - These are much harder than flow/051: high degree, thousands of terms,", flush=True)
    print("    strict residual, and starts not placed exactly on the planted root.", flush=True)
    print("  - This benchmark is meant to expose scaling and failure modes, not to be", flush=True)
    print("    friendly to Pandrosion.", flush=True)
    print("  - If Newton wins here, that is evidence that universal Q_F alone does not", flush=True)
    print("    provide the special contraction we get on exact S_A systems.", flush=True)
    print("=" * 132, flush=True)


if __name__ == "__main__":
    main()
