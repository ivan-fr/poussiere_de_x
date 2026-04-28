"""
Vary the *internal dimension* of each universal geometry and measure
coverage / time / F-evaluations.  Goal: find the sweet spot per geometry
and identify if there is a universally dominant choice.

Definitions of "dimension" per geometry:
  - path     : T_k accelerator order (paper 1 §4.2, paper 8 §3).
               k=1 (no accel), k=2 (Aitken), k=3 (Richardson), k=4..6.
  - cube     : max_orders = number of coordinate permutations averaged
               (paper 1 §3.5).  Caps at n!.
  - simplex  : k_pivots = number of independent edge-frame pivots averaged.
               k=1 (default) up to k=6.
  - sparse   : activity_power = exponent applied to support activity weights
               (heuristic).  Power=1 (default), 2, 4, etc.

Test bench: (n, d) = (3, 3), Bezout = 27, 3 seeds.  Use combined orbit
with single-geometry only (no portfolio mixing) to isolate the effect
of dimension.

Each (geom, dim) row reports avg over 3 seeds:
  - coverage
  - wall time
  - F-evaluations
  - paths_ok / paths_fail
  - collisions
"""
from __future__ import annotations
import cmath, itertools, json, math, os, random, sys, time
from itertools import product as iprod, permutations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bench_064_vs_059 import (
    F_eval, eval_poly, is_finite_vec, residual_norm, degree, solve_linear,
    schmidt_path, edge_frame, gen_random_poly_system,
    homotopy_polys, total_degree_start, start_roots, correct_T2K1_portfolio,
)


# ============================================================================
# Geometries with explicit internal dimension parameter
# ============================================================================
def geom_path(polys, anchor, z, F_anchor, dim, counter=None):
    """Path geometry: just the standard telescope (dim=1 baseline)."""
    return schmidt_path(polys, anchor, z, tuple(range(len(z))), counter=counter)


def geom_cube(polys, anchor, z, F_anchor, dim, counter=None):
    """Cube geometry: average max_orders=dim coordinate permutations."""
    n = len(z)
    perms = list(permutations(range(n)))
    if dim < len(perms):
        perms = perms[:dim]
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for order in perms:
        Q = schmidt_path(polys, anchor, z, order, counter=counter)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Q[i][j] / len(perms)
    return Q_avg


def geom_simplex(polys, anchor, z, F_anchor, dim, counter=None):
    """Simplex geometry with k_pivots independent edge-frames averaged.
    For dim=1: standard edge_frame + defect correction (=baseline simplex).
    For dim>=2: average dim edge_frames built from dim distinct anchor
    perturbations."""
    n = len(z)
    if dim <= 1:
        # Standard single-frame simplex with defect correction
        delta = [z[i] - anchor[i] for i in range(n)]
        B = edge_frame(polys, anchor, z, F_anchor, counter=counter)
        if sum(abs(v) ** 2 for v in delta) < 1e-28:
            return B
        F_z = F_eval(polys, z, counter=counter)
        Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
        defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
        w = [v.conjugate() for v in delta]
        denom = sum(w[j] * delta[j] for j in range(n))
        if abs(denom) < 1e-28:
            return B
        Q = [row[:] for row in B]
        for i in range(n):
            for j in range(n):
                Q[i][j] += defect[i] * w[j] / denom
        return Q
    # Multi-frame simplex: average dim independent edge-frames
    rng = random.Random(20260429)
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for k in range(dim):
        # Perturb anchor for each frame
        a_pert = [a + complex(0.05 * rng.gauss(0, 1), 0.05 * rng.gauss(0, 1))
                  for a in anchor]
        F_a = F_eval(polys, a_pert, counter=counter)
        delta = [z[i] - a_pert[i] for i in range(n)]
        B = edge_frame(polys, a_pert, z, F_a, counter=counter)
        if sum(abs(v) ** 2 for v in delta) < 1e-28:
            for i in range(n):
                for j in range(n):
                    Q_avg[i][j] += B[i][j] / dim
            continue
        F_z = F_eval(polys, z, counter=counter)
        Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
        defect = [F_z[i] - F_a[i] - Bd[i] for i in range(n)]
        w = [v.conjugate() for v in delta]
        denom = sum(w[j] * delta[j] for j in range(n))
        if abs(denom) < 1e-28:
            for i in range(n):
                for j in range(n):
                    Q_avg[i][j] += B[i][j] / dim
            continue
        Q = [row[:] for row in B]
        for i in range(n):
            for j in range(n):
                Q[i][j] += defect[i] * w[j] / denom
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Q[i][j] / dim
    return Q_avg


def geom_sparse(polys, anchor, z, F_anchor, dim, counter=None):
    """Sparse geometry: activity_power = dim.  power=1 default."""
    n = len(z)
    activity = [0.0] * n
    for poly in polys:
        for exp, coeff in poly.items():
            mag = abs(coeff)
            for j, e in enumerate(exp):
                activity[j] += mag * (e ** dim)        # power = dim
    m = max(activity) if activity else 0.0
    activity = [a if a > 0 else max(1.0, 0.1 * m) for a in activity]
    delta = [z[i] - anchor[i] for i in range(n)]
    B = edge_frame(polys, anchor, z, F_anchor, counter=counter)
    F_z = F_eval(polys, z, counter=counter)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [activity[j] * delta[j].conjugate() for j in range(n)]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


GEOMS = {"path": geom_path, "cube": geom_cube,
         "simplex": geom_simplex, "sparse": geom_sparse}


# ============================================================================
# T_n accelerator hierarchy (paper 1 §4.2 + paper 8 §3)
# ============================================================================
def pandrosion_h_with_geom(polys, anchor, z, F_anchor, geom_func, dim, counter=None):
    Q = geom_func(polys, anchor, z, F_anchor, dim, counter=counter)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), False
    out = [anchor[i] + step[i] for i in range(len(z))]
    return (out, True) if is_finite_vec(out) else (list(z), False)


def Tk_step(polys, anchor, z, F_anchor, geom_func, dim, k, counter=None):
    """T_k accelerator: produce s_1, ..., s_k via h, then Aitken / Richardson.
    k=1: just z -> s_1 (no acceleration).
    k=2: Aitken Delta-squared.
    k>=3: iterated Richardson with lambda_hat.
    """
    if k <= 1:
        s, ok = pandrosion_h_with_geom(polys, anchor, z, F_anchor, geom_func, dim, counter)
        return s, ok
    # Generate s_1, ..., s_k
    s_list = [list(z)]
    for _ in range(k):
        s_new, ok = pandrosion_h_with_geom(polys, anchor, s_list[-1], F_anchor,
                                            geom_func, dim, counter)
        if not ok:
            return s_list[-1], len(s_list) > 1
        s_list.append(s_new)
    # Aitken on (z=s_0, s_1, s_2)
    n = len(z)
    t2 = []
    for j in range(n):
        d0 = s_list[1][j] - s_list[0][j]
        d2 = s_list[2][j] - 2 * s_list[1][j] + s_list[0][j]
        t2.append(s_list[2][j] if abs(d2) < 1e-300 else s_list[0][j] - d0 * d0 / d2)
    if k == 2:
        return t2, True
    # Richardson higher orders (paper 1 §4.2)
    out = list(t2)
    lam_den = sum(s_list[1][j] - s_list[0][j] for j in range(n))
    if abs(lam_den) < 1e-300:
        return out, True
    lam_hat = sum(s_list[2][j] - s_list[1][j] for j in range(n)) / lam_den
    if abs(lam_hat - 1.0) < 1e-12:
        return s_list[-1], True
    for kk in range(3, k + 1):
        out = [out[j] - (s_list[kk][j] - out[j]) / (lam_hat - 1.0) for j in range(n)]
    return out, True


# ============================================================================
# Combined orbit with explicit (geometry, dim, T_k) parameters
# ============================================================================
def orbit(polys, z_init, geom_name, geom_dim, T_k=2, tol=1e-9,
          max_epochs=60, counter=None):
    n = len(z_init); z = list(z_init)
    anchor = [zi + complex(0.13 * (i + 1), 0.07) for i, zi in enumerate(z_init)]
    geom_func = GEOMS[geom_name]
    for epoch in range(max_epochs):
        F_anchor = F_eval(polys, anchor, counter)
        if is_finite_vec(F_anchor) and max(abs(v) for v in F_anchor) < tol:
            return anchor, True
        r_z = residual_norm(polys, z, counter)
        if r_z < tol:
            return z, True
        z_new, ok = Tk_step(polys, anchor, z, F_anchor, geom_func, geom_dim, T_k, counter)
        r_new = residual_norm(polys, z_new, counter) if ok else float("inf")
        if not ok or not math.isfinite(r_new):
            anchor = [a + complex(0.05 * random.gauss(0, 1), 0.05 * random.gauss(0, 1)) for a in anchor]
            continue
        if r_new < r_z:
            anchor = list(z); z = z_new
        else:
            # half-step retry
            z_half = [z[k] + 0.5 * (z_new[k] - z[k]) for k in range(n)]
            r_half = residual_norm(polys, z_half, counter)
            if r_half < r_z:
                anchor = list(z); z = z_half
            else:
                anchor = [a + complex(0.05 * random.gauss(0, 1), 0.05 * random.gauss(0, 1)) for a in anchor]
    return z, residual_norm(polys, z, counter) < tol


def gen_starts(n, count, seed=20260427):
    rng = random.Random(seed); starts = []
    for k in range(count):
        u = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
        norm = sum(abs(v)**2 for v in u) ** 0.5 or 1.0
        u = [v / norm for v in u]
        rho = 2.0 * (0.7 ** k)
        if rho < 0.05: rho = 0.05 + 0.5 * rng.random()
        starts.append([rho * v for v in u])
    return starts


def is_new(z, found, tol=1e-4):
    return all(max(abs(z[i] - r[i]) for i in range(len(z))) > tol for r in found)


def coverage_run(polys, geom_name, geom_dim, T_k=2, n_factor=3, n_min=12, tol=1e-9):
    n = len(next(iter(polys[0])))
    bez = 1
    for p in polys:
        bez *= max(degree(p), 1)
    n_starts = max(n_min, n_factor * bez)
    starts = gen_starts(n, n_starts)
    counter = [0]; found = []; ok = fail = 0
    t0 = time.time()
    for z0 in starts:
        z, conv = orbit(polys, z0, geom_name, geom_dim, T_k=T_k, tol=tol, counter=counter)
        if conv:
            ok += 1
            if is_new(z, found):
                found.append(z)
        else:
            fail += 1
    elapsed = time.time() - t0
    return {"bezout": bez, "found": len(found), "coverage": len(found)/max(bez, 1),
            "time": elapsed, "F_evals": counter[0],
            "paths_ok": ok, "paths_fail": fail,
            "collisions": ok - len(found)}


# ============================================================================
# Bench harness
# ============================================================================
RESULTS_FILE = "/tmp/bench_geom_dim_results.json"


def load_results():
    if os.path.exists(RESULTS_FILE):
        with open(RESULTS_FILE) as f: return json.load(f)
    return []


def save_results(rs):
    with open(RESULTS_FILE, "w") as f: json.dump(rs, f, indent=1)


def already_done(results, n, d, seed, geom, dim, Tk):
    return any(r["n"] == n and r["d"] == d and r["seed"] == seed
               and r["geom"] == geom and r["dim"] == dim and r["Tk"] == Tk
               for r in results)


def run_chunk(n, d, seeds, geom_dim_pairs, timeout=40):
    results = load_results()
    t_start = time.time()
    print(f"  {'(n,d,seed)':>14} {'geom':>9} {'dim':>3} {'Tk':>3} {'Bez':>4} | "
          f"{'cov%':>6} {'time':>7} {'F-evals':>10} {'ok/fail':>10} {'coll':>4}",
          flush=True)
    for seed in seeds:
        polys = gen_random_poly_system(n, d, seed=seed)
        for (geom, dim, Tk) in geom_dim_pairs:
            if already_done(results, n, d, seed, geom, dim, Tk):
                continue
            if time.time() - t_start > timeout - 5:
                print("  -- chunk budget exhausted --", flush=True)
                save_results(results); return
            r = coverage_run(polys, geom, dim, T_k=Tk, n_factor=3, n_min=12, tol=1e-9)
            row = {"n": n, "d": d, "seed": seed, "geom": geom, "dim": dim, "Tk": Tk,
                   **r}
            results.append(row); save_results(results)
            print(f"  ({n:>2},{d:>2},{seed % 100:>2}) {geom:>9} {dim:>3} {Tk:>3} "
                  f"{r['bezout']:>4} | "
                  f"{100*r['coverage']:>5.1f}% {r['time']:>6.2f}s "
                  f"{r['F_evals']:>10} {r['paths_ok']:>4}/{r['paths_fail']:<3} "
                  f"{r['collisions']:>4}", flush=True)
    save_results(results)


def report():
    results = load_results()
    if not results:
        print("no results"); return
    # Group by (geom, dim, Tk), avg over seeds
    by_key = {}
    for r in results:
        key = (r["geom"], r["dim"], r["Tk"])
        by_key.setdefault(key, []).append(r)
    print()
    print(f"  {'geom':>9} {'dim':>3} {'Tk':>3} {'#s':>3} | "
          f"{'cov':>6} {'time':>7} {'F-evals':>10} {'eff(c/F)':>10} {'collis':>7}")
    print("-" * 80)
    for key in sorted(by_key.keys()):
        rs = by_key[key]
        avg = lambda k: sum(r[k] for r in rs) / len(rs)
        eff = avg("coverage") / max(avg("F_evals"), 1) * 1e6
        print(f"  {key[0]:>9} {key[1]:>3} {key[2]:>3} {len(rs):>3} | "
              f"{100*avg('coverage'):>5.1f}% {avg('time'):>6.2f}s "
              f"{avg('F_evals'):>10.0f} {eff:>9.2f}x {avg('collisions'):>6.1f}")


CHUNKS = {
    # path = T_k accelerator order, T_k in {2,3,4,5,6,8,12}
    "path1": [("path", 1, 2), ("path", 1, 3), ("path", 1, 4)],
    "path2": [("path", 1, 5), ("path", 1, 6)],
    "path3": [("path", 1, 8), ("path", 1, 12)],
    # cube = max_orders, capped at n! (=6 for n=3). On (4,2) test we go higher.
    "cube":  [("cube", k, 2) for k in [2, 3, 4]],
    "cube2": [("cube", k, 2) for k in [5, 6]],
    "cube3": [("cube", k, 2) for k in [12, 24]],
    # simplex = number of pivots averaged (2..7)
    "simp":  [("simplex", k, 2) for k in [2, 3, 4]],
    "simp2": [("simplex", k, 2) for k in [5, 6]],
    "simp3": [("simplex", k, 2) for k in [7]],
    # sparse = activity power
    "sparse":  [("sparse", k, 2) for k in [1, 2, 3]],
    "sparse2": [("sparse", k, 2) for k in [4, 5, 6]],
    "sparse3": [("sparse", k, 2) for k in [7]],
}


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "report"
    if cmd == "reset":
        if os.path.exists(RESULTS_FILE): os.remove(RESULTS_FILE)
        print("reset")
    elif cmd == "report":
        report()
    else:
        # Test on (3,3) Bezout=27, 3 seeds
        seeds = [70033, 70034, 70035]
        n, d = 3, 3
        if cmd in CHUNKS:
            run_chunk(n, d, seeds, CHUNKS[cmd])
        else:
            print(f"unknown chunk {cmd}, options: {list(CHUNKS.keys())}")
