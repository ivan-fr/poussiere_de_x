"""
Rigorous benchmark of flow/065 (gamma-retry homotopy) vs flow/059 (Lairez gamma).
Both use the same Pandrosion T_2 K=1 corrector + 4 universal geometries + Armijo
non-holomorphic fallback.  The only difference:
  - 065: gamma=1 first (straight-line); on collisions, retry with random gamma.
  - 059: gamma=e^{0.73i} fixed (Lairez safety from the start).

Metrics per (n, d, seed) for each method:
  - coverage:        distinct roots found / Bezout
  - wall_time:       total seconds
  - F_evals:         polynomial evaluations
  - paths_ok:        paths that the tracker reports as converged
  - paths_fail:      paths that failed (residual > tol)
  - collisions:      paths_ok - distinct_roots
  - retries:         (065 only) number of additional gamma rotations used

Run mode 'all' may exceed sandbox 45s; use chunk modes (chunk1, chunk2, ...)
and append results to /tmp/bench_065_results.json then aggregate.
"""

from __future__ import annotations
import cmath, itertools, json, math, os, random, sys, time
from itertools import product as iprod

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bench_064_vs_059 import (
    F_eval, eval_poly, is_finite_vec, residual_norm, degree, solve_linear,
    correct_T2K1_portfolio,
    homotopy_polys, total_degree_start, start_roots, gen_random_poly_system
)


# ---------- Path tracker with detailed metrics ----------
def track_path_metric(target, start, gamma, z0, tol=1e-9, counter=None):
    t = 0.0; dt = 0.005; z = list(z0); n = len(z)
    prev_z, prev_t, failures = None, None, 0
    while t < 1.0 - 1e-15:
        trial_dt = min(dt, 1.0 - t); t_next = t + trial_dt
        if prev_z is None:
            pred = list(z)
        else:
            slope = [(z[i] - prev_z[i]) / (t - prev_t) for i in range(n)]
            pred = [z[i] + trial_dt * slope[i] for i in range(n)]
        H = homotopy_polys(target, start, gamma, t_next)
        z_corr, ok = correct_T2K1_portfolio(H, pred, tol=tol, max_epochs=20, counter=counter)
        if ok and residual_norm(H, z_corr, counter) < 10.0 * tol:
            prev_z, prev_t = list(z), t
            z, t = z_corr, t_next
            dt = min(0.02, dt * 1.1)
        else:
            failures += 1; dt *= 0.5
            if dt < 1e-5 or failures > 60:
                return z, False
    return z, residual_norm(target, z, counter) < 1e-7


def run_method_059(target, gamma=cmath.exp(0.73j), tol=1e-9):
    """Single-pass with fixed gamma (Lairez)."""
    n = len(target); degrees = [max(1, degree(p)) for p in target]
    start_sys = total_degree_start(degrees, n); z0_list = start_roots(degrees)
    bez = math.prod(degrees)
    found = []; counter = [0]; ok_count = fail_count = 0
    for z0 in z0_list:
        z, ok = track_path_metric(target, start_sys, gamma, z0, tol=tol, counter=counter)
        if ok:
            ok_count += 1
            if all(max(abs(z[i] - r[i]) for i in range(n)) > 1e-4 for r in found):
                found.append(z)
        else:
            fail_count += 1
    return {"bezout": bez, "found": len(found), "coverage": len(found)/max(bez,1),
            "F_evals": counter[0], "paths_ok": ok_count, "paths_fail": fail_count,
            "collisions": ok_count - len(found), "retries": 0}


def run_method_065(target, tol=1e-9, max_retry=3):
    """Pass 1 gamma=1 straight-line; retry collisions/failures with random gamma."""
    n = len(target); degrees = [max(1, degree(p)) for p in target]
    start_sys = total_degree_start(degrees, n); z0_list = start_roots(degrees)
    bez = math.prod(degrees)
    found = []; counter = [0]; total_ok = total_fail = 0
    # Track which starts produced which root
    path_root = {}    # tuple(z0) -> z if ok else None
    for z0 in z0_list:
        z, ok = track_path_metric(target, start_sys, 1.0, z0, tol=tol, counter=counter)
        if ok:
            total_ok += 1
            path_root[tuple(z0)] = z
            if all(max(abs(z[i] - r[i]) for i in range(n)) > 1e-4 for r in found):
                found.append(z)
        else:
            total_fail += 1
            path_root[tuple(z0)] = None
    retries = 0
    rng = random.Random(20260428)
    for _ in range(max_retry):
        if len(found) >= bez: break
        retries += 1
        gamma_new = cmath.exp(2j * math.pi * rng.random())
        for z0 in z0_list:
            if len(found) >= bez: break
            # Skip starts already at a found root
            z_prev = path_root.get(tuple(z0))
            if z_prev is not None and any(max(abs(z_prev[i] - r[i]) for i in range(n)) <= 1e-4 for r in found):
                # Already accounted for; retry from this start might find a new root
                pass
            z, ok = track_path_metric(target, start_sys, gamma_new, z0, tol=tol, counter=counter)
            if ok and all(max(abs(z[i] - r[i]) for i in range(n)) > 1e-4 for r in found):
                found.append(z); total_ok += 1
            elif not ok:
                total_fail += 1
    return {"bezout": bez, "found": len(found), "coverage": len(found)/max(bez,1),
            "F_evals": counter[0], "paths_ok": total_ok, "paths_fail": total_fail,
            "collisions": total_ok - len(found), "retries": retries}


# ---------- Bench runner with persistence ----------
RESULTS_FILE = "/tmp/bench_065_results.json"


def load_results():
    if os.path.exists(RESULTS_FILE):
        with open(RESULTS_FILE) as f:
            return json.load(f)
    return []


def save_results(results):
    with open(RESULTS_FILE, "w") as f:
        json.dump(results, f, indent=1)


def already_done(results, n, d, seed_offset):
    return any(r["n"] == n and r["d"] == d and r["seed"] == seed_offset for r in results)


def run_chunk(configs, seeds_per_config, timeout_total=40.0):
    results = load_results()
    t0 = time.time()
    print(f"  {'(n,d,seed)':>14} {'Bez':>5} | "
          f"{'cov065':>7} {'cov059':>7} {'t065':>7} {'t059':>7} "
          f"{'fail065':>8} {'fail059':>8} {'coll065':>8} {'coll059':>8} "
          f"{'ret065':>7} {'F065/F059':>14}", flush=True)
    for (n, d) in configs:
        for seed_offset in range(seeds_per_config):
            if already_done(results, n, d, seed_offset):
                continue
            if time.time() - t0 > timeout_total - 5:
                print("  --- chunk budget exhausted, stopping ---", flush=True)
                save_results(results); return False
            seed = 65000 + 1000 * n + 100 * d + seed_offset
            polys = gen_random_poly_system(n, d, seed=seed)
            tt = time.time(); r065 = run_method_065(polys, tol=1e-9); t065 = time.time() - tt
            tt = time.time(); r059 = run_method_059(polys, tol=1e-9); t059 = time.time() - tt
            row = {"n": n, "d": d, "seed": seed_offset, "bezout": r065["bezout"],
                   "cov_065": r065["coverage"], "cov_059": r059["coverage"],
                   "time_065": t065, "time_059": t059,
                   "F_065": r065["F_evals"], "F_059": r059["F_evals"],
                   "ok_065": r065["paths_ok"], "ok_059": r059["paths_ok"],
                   "fail_065": r065["paths_fail"], "fail_059": r059["paths_fail"],
                   "coll_065": r065["collisions"], "coll_059": r059["collisions"],
                   "retries_065": r065["retries"]}
            results.append(row); save_results(results)
            print(f"  ({n:>2},{d:>2},{seed_offset}) {row['bezout']:>5} | "
                  f"{100*row['cov_065']:>6.1f}% {100*row['cov_059']:>6.1f}% "
                  f"{row['time_065']:>6.1f}s {row['time_059']:>6.1f}s "
                  f"{row['fail_065']:>8} {row['fail_059']:>8} "
                  f"{row['coll_065']:>8} {row['coll_059']:>8} "
                  f"{row['retries_065']:>7} {row['F_065']:>6}/{row['F_059']:<6}",
                  flush=True)
    save_results(results)
    return True


def report():
    results = load_results()
    if not results:
        print("No results yet."); return
    by_nd = {}
    for r in results:
        key = (r["n"], r["d"]); by_nd.setdefault(key, []).append(r)
    print()
    print("AGGREGATED RESULTS")
    print(f"  {'(n,d)':>6} {'Bez':>5} {'#s':>3} | "
          f"{'cov065':>8} {'cov059':>8} | "
          f"{'t065':>7} {'t059':>7} ratio | "
          f"{'F065':>10} {'F059':>10} {'F-ratio':>8} | "
          f"{'fail065':>8} {'fail059':>8} | "
          f"{'coll065':>8} {'coll059':>8} | "
          f"{'ret065':>7}")
    print("-" * 162)
    for (n, d) in sorted(by_nd.keys()):
        rs = by_nd[(n, d)]; bez = rs[0]["bezout"]
        avg = lambda key: sum(r[key] for r in rs) / len(rs)
        time_ratio = avg("time_065") / max(avg("time_059"), 1e-9)
        F_ratio = avg("F_065") / max(avg("F_059"), 1.0)
        print(f"  ({n:>2},{d:>2}) {bez:>5} {len(rs):>3} | "
              f"{100*avg('cov_065'):>7.1f}% {100*avg('cov_059'):>7.1f}% | "
              f"{avg('time_065'):>6.1f}s {avg('time_059'):>6.1f}s {time_ratio:>5.2f} | "
              f"{avg('F_065'):>10.0f} {avg('F_059'):>10.0f} {F_ratio:>7.2f}x | "
              f"{avg('fail_065'):>8.1f} {avg('fail_059'):>8.1f} | "
              f"{avg('coll_065'):>8.1f} {avg('coll_059'):>8.1f} | "
              f"{avg('retries_065'):>6.1f}")


CONFIG_CHUNKS = {
    "1": [(2, 3), (3, 2), (4, 2)],
    "2": [(2, 5), (3, 3)],
    "3": [(2, 8), (3, 4)],
    "4": [(4, 3), (5, 2)],
    "5": [(2, 12)],
    "report": [],
}

if __name__ == "__main__":
    chunk = sys.argv[1] if len(sys.argv) > 1 else "1"
    seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    if chunk == "report":
        report()
    elif chunk == "reset":
        if os.path.exists(RESULTS_FILE):
            os.remove(RESULTS_FILE); print("results reset")
    else:
        configs = CONFIG_CHUNKS.get(chunk, [(2, 3)])
        print(f"Chunk {chunk}: configs {configs}, {seeds} seeds each", flush=True)
        run_chunk(configs, seeds)
