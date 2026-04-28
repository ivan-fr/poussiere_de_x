"""
PAPER: 068
TITLE: dD-simplex universal Pandrosion bridge — A (multi-anchor collinear,
       exact telescoping) and B (simplex-perturbed multi-anchor average,
       approximate telescoping) compared to first-order paper-1 path.

THEORY
======
Companion to /Users/ivanbesevic/Documents/poussiere/latex/dD_simplex_universal_bridge.tex .

Paper 0 dD: S_p^{dD}(s_1,...,s_{d-1}) = sum_{|alpha|<=p-1} prod s_i^{alpha_i}
            = [1, s_1, ..., s_{d-1}](z^p)  (Newton divided difference).

Paper 1: Q_F(z_0, z) Schmidt slope matrix; F(z) - F(z_0) = Q_F (z - z_0).

CONSTRUCTION A (collinear multi-anchor, exact telescoping)
-----------------------------------------------------------
For order d >= 2, anchors w_k = z_0 + tau_k * (z - z_0), tau_k = k/(d-1):
    Q_A^{(d)}(z_0, z) = sum_{k=0}^{d-2} (tau_{k+1} - tau_k) * Q_F(w_k, w_{k+1})
Telescoping holds exactly (Theorem 4.1 of bridge.tex).

CONSTRUCTION B (simplex-perturbed average, approximate telescoping)
-------------------------------------------------------------------
For lattice Lambda_{p-1}^{(n)} = {alpha in N^n : |alpha| <= p-1}:
    delta_alpha = epsilon * (alpha - bar_alpha) / (p-1)
    w_alpha = z_0 + delta_alpha
    Q_B(z_0, z) = (1/|Lambda|) sum_{alpha} Q_F(w_alpha, z)
Telescoping is approximate (off by O(epsilon)).

OUTER SCHEME
============
T_2 K=1 (Aitken Delta^2 of P_{F,z_0}) with anchor reset every step.
Same scaffold as flow/067.

GEOMETRIES TESTED
=================
  - "path"          : paper-1 first-order canonical (control)
  - "A_2"           : collinear order 2 = paper-1 path (sanity check)
  - "A_3", "A_4"    : collinear order 3, 4 (multi-anchor)
  - "B_p2"          : simplex-perturbed lattice with p=2
  - "B_p3"          : simplex-perturbed lattice with p=3

BENCH
=====
Configs: (2,3), (3,2), (4,2), (3,3) — Bezout 9, 8, 16, 27.
3 seeds per config.
Persistence: /tmp/bench_068_results.json.
"""
from __future__ import annotations
import cmath, itertools, json, math, os, random, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bench_064_vs_059 import (
    F_eval, eval_poly, residual_norm, is_finite_vec, degree, solve_linear,
    schmidt_path, edge_frame, gen_random_poly_system,
)


# ============================================================================
# Construction A — collinear multi-anchor, exact telescoping
# ============================================================================
def Q_A(polys, z0, z, F_z0, d, counter=None):
    """Construction A of order d (Theorem 4.1, bridge.tex).

    Anchors w_k = z_0 + (k/(d-1)) * (z - z_0), k = 0, ..., d-1.
    Q_A^{(d)} = sum_{k=0}^{d-2} (1/(d-1)) * Q_F(w_k, w_{k+1}).

    For d = 2: collapses to standard paper-1 path Schmidt slope.
    """
    n = len(z)
    if d <= 2:
        return schmidt_path(polys, z0, z, tuple(range(n)), counter=counter)
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    weight = 1.0 / (d - 1)  # uniform tau spacing
    for k in range(d - 1):
        tau_k = k / (d - 1)
        tau_k1 = (k + 1) / (d - 1)
        w_k = [z0[i] + tau_k * (z[i] - z0[i]) for i in range(n)]
        w_k1 = [z0[i] + tau_k1 * (z[i] - z0[i]) for i in range(n)]
        Qk = schmidt_path(polys, w_k, w_k1, tuple(range(n)), counter=counter)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += weight * Qk[i][j]
    return Q_avg


# ============================================================================
# Construction B — simplex-perturbed lattice, approximate telescoping
# ============================================================================
def _simplex_lattice(n, p_minus_1):
    """Yield all alpha in N^n with sum(alpha) <= p_minus_1."""
    if n == 0:
        yield ()
        return
    for k in range(p_minus_1 + 1):
        for rest in _simplex_lattice(n - 1, p_minus_1 - k):
            yield (k,) + rest


def Q_B(polys, z0, z, F_z0, p, eps, counter=None):
    """Construction B of order p (Definition 5.1, bridge.tex).

    Lattice Lambda_{p-1}^{(n)} = {alpha : |alpha| <= p-1}.
    Each anchor w_alpha = z_0 + eps * (alpha - bar_alpha) / (p-1).
    Q_B = (1/|Lambda|) sum_alpha Q_F(w_alpha, z).
    """
    n = len(z)
    lattice = list(_simplex_lattice(n, p - 1))
    L = len(lattice)
    if L == 0:
        return schmidt_path(polys, z0, z, tuple(range(n)), counter=counter)
    bar_alpha = [0.0] * n
    for alpha in lattice:
        for i in range(n):
            bar_alpha[i] += alpha[i] / L
    denom = max(p - 1, 1)
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for alpha in lattice:
        delta = [eps * (alpha[i] - bar_alpha[i]) / denom for i in range(n)]
        w_alpha = [z0[i] + complex(delta[i], 0) for i in range(n)]
        Qa = schmidt_path(polys, w_alpha, z, tuple(range(n)), counter=counter)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Qa[i][j] / L
    return Q_avg


def Q_path(polys, z0, z, F_z0, counter=None):
    return schmidt_path(polys, z0, z, tuple(range(len(z))), counter=counter)


# ============================================================================
# Pandrosion operator + T_2 corrector (paper 1 §3.5 + paper 8)
# ============================================================================
def pandrosion_h(polys, z0, z, F_z0, Q_func, counter=None):
    Q = Q_func(polys, z0, z, F_z0, counter=counter)
    step = solve_linear(Q, [-v for v in F_z0])
    if step is None:
        return list(z), False
    out = [z0[i] + step[i] for i in range(len(z))]
    return (out, True) if is_finite_vec(out) else (list(z), False)


def T2_step(polys, z0, z, F_z0, Q_func, counter=None):
    s1, ok1 = pandrosion_h(polys, z0, z, F_z0, Q_func, counter)
    if not ok1:
        return list(z), False
    s2, ok2 = pandrosion_h(polys, z0, s1, F_z0, Q_func, counter)
    if not ok2:
        return s1, True
    n = len(z); out = []
    for k in range(n):
        d0 = s1[k] - z[k]; d2 = s2[k] - 2 * s1[k] + z[k]
        out.append(s2[k] if abs(d2) < 1e-300 else z[k] - d0 * d0 / d2)
    return out, True


# ============================================================================
# Geometry registry
# ============================================================================
def geom_factory(name):
    if name == "path":
        return lambda polys, z0, z, F_z0, counter=None: Q_path(polys, z0, z, F_z0, counter)
    if name == "A_2":
        return lambda polys, z0, z, F_z0, counter=None: Q_A(polys, z0, z, F_z0, 2, counter)
    if name == "A_3":
        return lambda polys, z0, z, F_z0, counter=None: Q_A(polys, z0, z, F_z0, 3, counter)
    if name == "A_4":
        return lambda polys, z0, z, F_z0, counter=None: Q_A(polys, z0, z, F_z0, 4, counter)
    if name == "A_5":
        return lambda polys, z0, z, F_z0, counter=None: Q_A(polys, z0, z, F_z0, 5, counter)
    if name == "B_p2":
        return lambda polys, z0, z, F_z0, counter=None: Q_B(polys, z0, z, F_z0, 2, 0.05, counter)
    if name == "B_p3":
        return lambda polys, z0, z, F_z0, counter=None: Q_B(polys, z0, z, F_z0, 3, 0.05, counter)
    raise ValueError(name)


GEOM_LIST = ["path", "A_2", "A_3", "A_4", "A_5", "B_p2", "B_p3"]


# ============================================================================
# Single-geometry orbit
# ============================================================================
def orbit(polys, z_init, geom_name, tol=1e-9, max_epochs=60, counter=None):
    n = len(z_init); z = list(z_init)
    anchor = [zi + complex(0.13 * (i + 1), 0.07) for i, zi in enumerate(z_init)]
    rng = random.Random(68000 + n + hash(geom_name) % 1000)
    Q_func = geom_factory(geom_name)
    for _ in range(max_epochs):
        F_anchor = F_eval(polys, anchor, counter)
        if is_finite_vec(F_anchor) and max(abs(v) for v in F_anchor) < tol:
            return anchor, True
        r_z = residual_norm(polys, z, counter)
        if r_z < tol:
            return z, True
        z_new, ok = T2_step(polys, anchor, z, F_anchor, Q_func, counter)
        r_new = residual_norm(polys, z_new, counter) if ok else float("inf")
        if not ok or not math.isfinite(r_new):
            anchor = [a + complex(0.05 * rng.gauss(0, 1), 0.05 * rng.gauss(0, 1))
                      for a in anchor]
            continue
        if r_new < r_z:
            anchor = list(z); z = z_new
        else:
            z_half = [z[k] + 0.5 * (z_new[k] - z[k]) for k in range(n)]
            r_half = residual_norm(polys, z_half, counter)
            if r_half < r_z:
                anchor = list(z); z = z_half
            else:
                anchor = [a + complex(0.05 * rng.gauss(0, 1), 0.05 * rng.gauss(0, 1))
                          for a in anchor]
    return z, residual_norm(polys, z, counter) < tol


def gen_starts(n, count, seed=20260428):
    rng = random.Random(seed); starts = []
    for k in range(count):
        u = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
        norm = sum(abs(v) ** 2 for v in u) ** 0.5 or 1.0
        u = [v / norm for v in u]
        rho = 2.0 * (0.7 ** k)
        if rho < 0.05:
            rho = 0.05 + 0.5 * rng.random()
        starts.append([rho * v for v in u])
    return starts


def is_new(z, found, tol=1e-4):
    return all(max(abs(z[i] - r[i]) for i in range(len(z))) > tol for r in found)


def coverage_run(polys, geom_name, n_factor=3, n_min=12, tol=1e-9):
    n = len(next(iter(polys[0])))
    bez = 1
    for p in polys:
        bez *= max(degree(p), 1)
    n_starts = max(n_min, n_factor * bez)
    starts = gen_starts(n, n_starts)
    counter = [0]; found = []; ok = fail = 0
    t0 = time.time()
    for z0 in starts:
        z, conv = orbit(polys, z0, geom_name, tol=tol, counter=counter)
        if conv:
            ok += 1
            if is_new(z, found):
                found.append(z)
        else:
            fail += 1
    elapsed = time.time() - t0
    return {"bezout": bez, "found": len(found),
            "coverage": len(found) / max(bez, 1),
            "time": elapsed, "F_evals": counter[0],
            "paths_ok": ok, "paths_fail": fail,
            "collisions": ok - len(found)}


# ============================================================================
# Bench harness
# ============================================================================
RESULTS_FILE = "/tmp/bench_068_results.json"


def load_results():
    if os.path.exists(RESULTS_FILE):
        with open(RESULTS_FILE) as f:
            return json.load(f)
    return []


def save_results(rs):
    with open(RESULTS_FILE, "w") as f:
        json.dump(rs, f, indent=1)


def already_done(results, n, d, seed, geom):
    return any(r["n"] == n and r["d"] == d and r["seed"] == seed
               and r["geom"] == geom for r in results)


def run_chunk(configs, seeds_per_config=3, geom_list=None, timeout=40):
    if geom_list is None:
        geom_list = GEOM_LIST
    results = load_results()
    t_start = time.time()
    print(f"  {'(n,d,s)':>9} {'geom':>8} {'Bez':>4} | "
          f"{'cov%':>6} {'time':>7} {'Fevals':>9} {'ok/fail':>10} {'coll':>4}",
          flush=True)
    for (n, d) in configs:
        for seed_off in range(seeds_per_config):
            seed = 68000 + 1000 * n + 100 * d + seed_off
            polys = gen_random_poly_system(n, d, seed=seed)
            for geom in geom_list:
                if already_done(results, n, d, seed_off, geom):
                    continue
                if time.time() - t_start > timeout - 5:
                    print("  -- chunk budget exhausted --", flush=True)
                    save_results(results); return False
                r = coverage_run(polys, geom)
                row = {"n": n, "d": d, "seed": seed_off, "geom": geom, **r}
                results.append(row); save_results(results)
                print(f"  ({n:>2},{d:>2},{seed_off}) {geom:>8} "
                      f"{r['bezout']:>4} | "
                      f"{100*r['coverage']:>5.1f}% {r['time']:>6.2f}s "
                      f"{r['F_evals']:>9} {r['paths_ok']:>4}/{r['paths_fail']:<3} "
                      f"{r['collisions']:>4}", flush=True)
    save_results(results)
    return True


def report():
    results = load_results()
    if not results:
        print("No results yet."); return
    by_key = {}
    for r in results:
        key = (r["geom"], r["n"], r["d"])
        by_key.setdefault(key, []).append(r)
    print()
    print("AGGREGATED 068 RESULTS  (T_2 K=1 + collinear A or simplex-perturbed B)")
    print("-" * 90)
    print(f"  {'geom':>8} {'(n,d)':>7} {'Bez':>4} {'#s':>3} | "
          f"{'cov%':>6} {'time':>7} {'F-evals':>10} {'eff(c/F)':>9}")
    print("-" * 90)
    for key in sorted(by_key.keys(), key=lambda k: (k[1], k[2], k[0])):
        geom, n, d = key
        rs = by_key[key]
        cov = 100 * sum(r["coverage"] for r in rs) / len(rs)
        tt = sum(r["time"] for r in rs) / len(rs)
        fe = sum(r["F_evals"] for r in rs) / len(rs)
        eff = (cov / 100.0) / max(fe, 1) * 1e6
        print(f"  {geom:>8} ({n:>2},{d:>2}) {rs[0]['bezout']:>4} "
              f"{len(rs):>3} | {cov:>5.1f}% {tt:>6.2f}s {fe:>10.0f} {eff:>8.2f}")


CONFIG_CHUNKS = {
    "1": [(2, 3)],
    "2": [(3, 2)],
    "3": [(4, 2)],
    "4": [(3, 3)],
    "all": [(2, 3), (3, 2), (4, 2), (3, 3)],
}


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "1"
    seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    if cmd == "report":
        report()
    elif cmd == "reset":
        if os.path.exists(RESULTS_FILE):
            os.remove(RESULTS_FILE); print("results reset")
    else:
        configs = CONFIG_CHUNKS.get(cmd, [(2, 3)])
        print(f"Chunk {cmd}: configs {configs}, {seeds} seeds each", flush=True)
        run_chunk(configs, seeds_per_config=seeds)
