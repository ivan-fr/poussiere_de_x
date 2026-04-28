"""
Side-by-side benchmark of flow/064 (Newton-homotopy light, no gamma) vs flow/059
(gamma-homotopy, Lairez-style).  Both use the same Pandrosion T_2 K=1 corrector
with the four universal Q_F geometries.

Metrics per (n, d, seed):
  - coverage %   (distinct roots found / Bezout)
  - wall time
  - paths ok / paths fail
  - F-evaluations
  - paths converging to duplicate roots ("collisions")

Output: CSV-friendly table and average summary.
"""
from __future__ import annotations
import cmath, itertools, math, random, time
from itertools import product as iprod


# -- Polynomial primitives -------------------------------------------------
def eval_poly(poly, z, counter=None):
    if counter is not None: counter[0] += 1
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


def F_eval(polys, z, counter=None):
    return [eval_poly(p, z, counter) for p in polys]


def is_finite_vec(z):
    return all(math.isfinite(complex(v).real) and math.isfinite(complex(v).imag) for v in z)


def residual_norm(polys, z, counter=None):
    if not is_finite_vec(z): return float("inf")
    vals = F_eval(polys, z, counter)
    if any(not math.isfinite(v.real) or not math.isfinite(v.imag) for v in vals):
        return float("inf")
    return max(abs(v) for v in vals)


def degree(poly):
    return max((sum(e) for e in poly), default=0)


# -- Linear algebra --------------------------------------------------------
def solve_linear(A, b):
    n = len(A); M = [list(row) + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-14: return None
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0.0 + 0.0j] * n
    for i in range(n - 1, -1, -1):
        rhs = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = rhs / M[i][i]
    return x


# -- Q_F geometries (4 modes from flow/061-064) ---------------------------
def schmidt_path(polys, anchor, z, order, counter=None):
    n = len(z); Q = [[0.0 + 0.0j] * n for _ in range(n)]
    cur = list(z); F_cur = F_eval(polys, cur, counter)
    for j in order:
        nxt = list(cur); nxt[j] = anchor[j]
        F_next = F_eval(polys, nxt, counter)
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300: continue
        for i in range(n):
            Q[i][j] = (F_cur[i] - F_next[i]) / denom
        cur, F_cur = nxt, F_next
    return Q


def schmidt_cube(polys, anchor, z, max_orders=8, counter=None):
    n = len(z); perms = list(itertools.permutations(range(n)))[:max_orders]
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for order in perms:
        Q = schmidt_path(polys, anchor, z, order, counter)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Q[i][j] / len(perms)
    return Q_avg


def edge_frame(polys, anchor, z, F_anchor, counter=None):
    n = len(z); B = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        edge = list(anchor); edge[j] = z[j]
        F_edge = F_eval(polys, edge, counter)
        denom = z[j] - anchor[j]
        if abs(denom) < 1e-300: continue
        for i in range(n):
            B[i][j] = (F_edge[i] - F_anchor[i]) / denom
    return B


def schmidt_simplex(polys, anchor, z, F_anchor, counter=None):
    n = len(z); delta = [z[i] - anchor[i] for i in range(n)]
    B = edge_frame(polys, anchor, z, F_anchor, counter)
    if sum(abs(v)**2 for v in delta) < 1e-28: return B
    F_z = F_eval(polys, z, counter)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [v.conjugate() for v in delta]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28: return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def schmidt_sparse(polys, anchor, z, F_anchor, counter=None):
    n = len(z); activity = [0.0] * n
    for poly in polys:
        for exp, coeff in poly.items():
            for j, e in enumerate(exp):
                activity[j] += abs(coeff) * e
    m = max(activity) if activity else 0.0
    activity = [a if a > 0 else max(1.0, 0.1 * m) for a in activity]
    delta = [z[i] - anchor[i] for i in range(n)]
    B = edge_frame(polys, anchor, z, F_anchor, counter)
    F_z = F_eval(polys, z, counter)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    w = [activity[j] * delta[j].conjugate() for j in range(n)]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28: return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def Q_geometry(polys, anchor, z, F_anchor, mode, counter=None):
    if mode == "path": return schmidt_path(polys, anchor, z, tuple(range(len(z))), counter)
    if mode == "cube": return schmidt_cube(polys, anchor, z, counter=counter)
    if mode == "simplex": return schmidt_simplex(polys, anchor, z, F_anchor, counter)
    if mode == "sparse": return schmidt_sparse(polys, anchor, z, F_anchor, counter)
    raise ValueError(mode)


# -- T_2 K=1 corrector (paper 8 + 4 geoms) --------------------------------
def pandrosion_h(polys, anchor, z, F_anchor, mode, counter=None):
    Q = Q_geometry(polys, anchor, z, F_anchor, mode, counter)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None: return list(z), False
    out = [anchor[i] + step[i] for i in range(len(z))]
    return (out, True) if is_finite_vec(out) else (list(z), False)


def T2_step(polys, anchor, z, F_anchor, mode, counter=None):
    s1, ok = pandrosion_h(polys, anchor, z, F_anchor, mode, counter)
    if not ok: return list(z), False
    s2, ok = pandrosion_h(polys, anchor, s1, F_anchor, mode, counter)
    if not ok: return s1, True
    n = len(z); out = []
    for k in range(n):
        d0 = s1[k] - z[k]; d2 = s2[k] - 2 * s1[k] + z[k]
        out.append(s2[k] if abs(d2) < 1e-300 else z[k] - d0 * d0 / d2)
    return out, True


def newton_ELS_step(polys, z, counter=None):
    n = len(z); F_z = F_eval(polys, z, counter)
    if not is_finite_vec(F_z): return list(z), False
    h = 1e-7
    J = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        zp = list(z); zp[j] += h; F_p = F_eval(polys, zp, counter)
        for i in range(n):
            J[i][j] = (F_p[i] - F_z[i]) / h
    direction = solve_linear(J, [-v for v in F_z])
    if direction is None: return list(z), False
    r0 = sum(abs(v)**2 for v in F_z); best_z, best_r = list(z), r0
    for k in range(-6, 7):
        tau = 2.0 ** k
        cand = [z[i] + tau * direction[i] for i in range(n)]
        if not is_finite_vec(cand): continue
        F_c = F_eval(polys, cand, counter)
        r = sum(abs(v)**2 for v in F_c) if is_finite_vec(cand) else float("inf")
        if r < best_r: best_z, best_r = cand, r
    return best_z, best_r < r0


def correct_T2K1_portfolio(polys, z_init, tol=1e-9, max_epochs=25, counter=None):
    n = len(z_init); z = list(z_init)
    anchor = [zi + complex(0.13 * (i + 1), 0.07) for i, zi in enumerate(z_init)]
    geom_modes = ["path", "simplex", "sparse"]
    if n <= 4: geom_modes.append("cube")
    for epoch in range(max_epochs):
        F_anchor = F_eval(polys, anchor, counter)
        if is_finite_vec(F_anchor) and max(abs(v) for v in F_anchor) < tol:
            return anchor, True
        r_z = residual_norm(polys, z, counter)
        if r_z < tol: return z, True
        best_r, best_z = r_z, z
        for mode in geom_modes:
            z_t2, ok = T2_step(polys, anchor, z, F_anchor, mode, counter)
            if ok:
                r = residual_norm(polys, z_t2, counter)
                if math.isfinite(r) and r < best_r: best_r, best_z = r, z_t2
        z_n, ok = newton_ELS_step(polys, z, counter)
        if ok:
            r = residual_norm(polys, z_n, counter)
            if math.isfinite(r) and r < best_r: best_r, best_z = r, z_n
        if best_r >= r_z:
            anchor = [a + complex(0.05*random.gauss(0,1), 0.05*random.gauss(0,1)) for a in anchor]
            continue
        anchor = list(z); z = best_z
    return z, residual_norm(polys, z, counter) < tol


# -- Homotopy --------------------------------------------------------------
def add_scaled(out, poly, scale):
    for exp, coeff in poly.items():
        out[exp] = out.get(exp, 0.0 + 0.0j) + scale * coeff


def clean_dict(d, eps=1e-14):
    return {k: v for k, v in d.items() if abs(v) > eps}


def homotopy_polys(target, start, gamma, t):
    polys = []
    for f, g in zip(target, start):
        h = {}; add_scaled(h, g, 1.0 - t); add_scaled(h, f, gamma * t)
        polys.append(clean_dict(h))
    return polys


def total_degree_start(degrees, n):
    polys = []; zero = tuple([0] * n)
    for i, d in enumerate(degrees):
        exp = tuple(d if j == i else 0 for j in range(n))
        polys.append({exp: 1.0 + 0.0j, zero: -1.0 + 0.0j})
    return polys


def roots_of_unity(d):
    return [cmath.exp(2j * math.pi * k / d) for k in range(d)]


def start_roots(degrees):
    return [list(root) for root in iprod(*[roots_of_unity(d) for d in degrees])]


def track_path(target, start, gamma, z0, tol=1e-9, counter=None):
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
                return z, False, failures
    return z, residual_norm(target, z, counter) < 1e-7, failures


def track_all(target, gamma=1.0, tol=1e-9):
    n = len(target); degrees = [max(1, degree(p)) for p in target]
    start_sys = total_degree_start(degrees, n)
    z0_list = start_roots(degrees)
    found = []; paths_ok = paths_fail = 0; total_failures = 0
    counter = [0]
    raw_endpoints = []
    for z0 in z0_list:
        z, ok, fails = track_path(target, start_sys, gamma, z0, tol=tol, counter=counter)
        total_failures += fails
        if ok:
            paths_ok += 1
            raw_endpoints.append(z)
            if all(max(abs(z[i] - r[i]) for i in range(n)) > 1e-4 for r in found):
                found.append(z)
        else:
            paths_fail += 1
    bez = math.prod(degrees)
    # Path collisions = paths that succeeded but collapsed to a duplicate root
    collisions = paths_ok - len(found)
    return {"roots": found, "bezout": bez, "paths_ok": paths_ok,
            "paths_fail": paths_fail, "collisions": collisions,
            "F_evals": counter[0], "tracker_failures": total_failures,
            "coverage": len(found) / max(bez, 1)}


def gen_random_poly_system(n, d, seed):
    rng = random.Random(seed); polys = []
    for i in range(n):
        poly = {}
        for alpha in iprod(*[range(d + 1) for _ in range(n)]):
            if sum(alpha) > d: continue
            if rng.random() < 0.7:
                poly[tuple(alpha)] = complex(rng.gauss(0.0, 1.0), 0.15 * rng.gauss(0.0, 1.0))
        diag = tuple(d if k == i else 0 for k in range(n))
        poly[diag] = poly.get(diag, 0.0 + 0.0j) + 1.0
        m = max(abs(v) for v in poly.values())
        polys.append({k: v / m for k, v in poly.items()})
    return polys


def run_bench(configs, seeds_per_config=3, gamma_064=1.0, gamma_059=cmath.exp(0.73j)):
    """Returns list of dicts with per-(n,d,seed) results for both methods."""
    results = []
    print(f"  {'(n,d,seed)':>14} {'Bez':>5} | "
          f"{'cov064':>7} {'cov059':>7} {'t064':>7} {'t059':>7} "
          f"{'fail064':>8} {'fail059':>8}", flush=True)
    for (n, d) in configs:
        for seed_offset in range(seeds_per_config):
            seed = 65000 + 1000 * n + 100 * d + seed_offset
            polys = gen_random_poly_system(n, d, seed=seed)
            t0 = time.time()
            r064 = track_all(polys, gamma=gamma_064, tol=1e-9)
            t064 = time.time() - t0
            t0 = time.time()
            r059 = track_all(polys, gamma=gamma_059, tol=1e-9)
            t059 = time.time() - t0
            row = {"n": n, "d": d, "seed": seed_offset, "bezout": r064["bezout"],
                   "cov_064": r064["coverage"], "cov_059": r059["coverage"],
                   "time_064": t064, "time_059": t059,
                   "paths_ok_064": r064["paths_ok"], "paths_ok_059": r059["paths_ok"],
                   "paths_fail_064": r064["paths_fail"], "paths_fail_059": r059["paths_fail"],
                   "coll_064": r064["collisions"], "coll_059": r059["collisions"],
                   "Fevals_064": r064["F_evals"], "Fevals_059": r059["F_evals"]}
            results.append(row)
            print(f"  ({n:>2},{d:>2},{seed_offset}) {row['bezout']:>5} | "
                  f"{100*row['cov_064']:>6.1f}% {100*row['cov_059']:>6.1f}% "
                  f"{row['time_064']:>6.1f}s {row['time_059']:>6.1f}s "
                  f"{row['paths_fail_064']:>8} {row['paths_fail_059']:>8}",
                  flush=True)
    return results


def report(results):
    """Group by (n, d) and average."""
    by_nd = {}
    for r in results:
        key = (r["n"], r["d"])
        by_nd.setdefault(key, []).append(r)
    print()
    print(f"  {'(n,d)':>6} {'Bez':>5} {'seeds':>5} | "
          f"{'cov 064':>8} {'cov 059':>8} | "
          f"{'t 064':>7} {'t 059':>7} | "
          f"{'fail 064':>8} {'fail 059':>8} | "
          f"{'coll 064':>8} {'coll 059':>8} | "
          f"{'Fev 064':>8} {'Fev 059':>8}")
    print("-" * 134)
    for (n, d), rs in sorted(by_nd.items()):
        bez = rs[0]["bezout"]
        avg = lambda key: sum(r[key] for r in rs) / len(rs)
        print(f"  ({n:>2},{d:>2}) {bez:>5} {len(rs):>5} | "
              f"{100*avg('cov_064'):>7.1f}% {100*avg('cov_059'):>7.1f}% | "
              f"{avg('time_064'):>6.1f}s {avg('time_059'):>6.1f}s | "
              f"{avg('paths_fail_064'):>8.1f} {avg('paths_fail_059'):>8.1f} | "
              f"{avg('coll_064'):>8.1f} {avg('coll_059'):>8.1f} | "
              f"{avg('Fevals_064'):>8.0f} {avg('Fevals_059'):>8.0f}")


if __name__ == "__main__":
    import sys
    # Allow command-line restriction of configs to fit in 45s sandbox budget
    arg_set = sys.argv[1] if len(sys.argv) > 1 else "small"
    if arg_set == "small":
        configs = [(2, 3), (2, 5), (3, 2), (3, 3), (4, 2)]
    elif arg_set == "medium":
        configs = [(2, 8), (3, 4), (4, 3), (5, 2)]
    elif arg_set == "med2":
        configs = [(3, 4), (4, 3)]
    elif arg_set == "med3":
        configs = [(5, 2), (4, 3)]
    elif arg_set == "large":
        configs = [(2, 12), (3, 5)]
    else:
        configs = [(2, 3), (2, 5), (2, 8), (3, 2), (3, 3), (3, 4), (4, 2), (4, 3), (5, 2)]
    print(f"Bench set: {arg_set}, configs: {configs}")
    results = run_bench(configs, seeds_per_config=3)
    report(results)
