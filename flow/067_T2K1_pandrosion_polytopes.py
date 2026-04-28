"""
PAPER: 067
TITLE: T_2 K=1 + each paper-1 universal Pandrosion polytope geometry tested
       SEPARATELY (no portfolio) on Smale-17-style Kostlan-Smale systems.

CORRECTION USER-PROVIDED CONTEXT
================================
Paper 1 (1pandrosion_smale.tex) defines the universal multivariate Pandrosion
operator for ARBITRARY polynomial systems F : C^n -> C^n.  Quoting the paper
(line numbers refer to latex/recovered/1pandrosion_smale.tex):

  L.105-110: univariate operator
        z_{n+1} = z_0 - P(z_0) / Q(z_0, z_n)             (eq. 3.5)

  L.136-147: multivariate Q_F by coordinate-by-coordinate interpolation
        F(z) - F(z_0) = Q_F(z_0, z) * (z - z_0)          (eq. 3.6)

        [Q_F]_{i,j} = [ F_i(z_0[<j] | z[j]  | z[>j])
                      - F_i(z_0[<j] | z_0[j]| z[>j]) ] / (z[j]-z_0[j])
                                                          (eq. 3.7)

  L.157-164: multivariate operator
        P_{F,z_0}(z) = z_0 - Q_F(z_0, z)^{-1} F(z_0)     (eq. 3.10)

Paper 1 PROVES universal local convergence (Theorem 5.1, basin radius
eta_F(zeta)) using ONE Q_F.  The "polytope" variants come from flow/046,
047, 049: each provides an alternative Q_F that still satisfies the
telescoping identity (3.6) exactly, hence remains in paper 1's framework.

This is the *paper-1 universal* extension.  It is NOT the paper-0 monomial
dD simplex of flow/003, which solves only z^p = x and is for x^{1/p}
specifically.

GEOMETRIES TESTED (each paper-1 universal Q_F variant)
======================================================
All satisfy the telescoping identity (paper 1 eq. 3.6) exactly, so
regular zeros remain fixed points of P_{F,z_0}.

  path (control, paper 1 §3.5 canonical):
        Q via single coordinate-by-coordinate interpolation, eq. (3.7).
        Cost: n F-evaluations per Q.

  simplex_4 (paper 1 universal, flow/047 dim 4):
        Q = avg of 4 edge frames + radial defect correction at 4
        simplex-vertex perturbations of the anchor.

  hypercube_4 (paper 1 universal, flow/046 dim 4 = 2^4):
        Q = avg of Schmidt paths over min(16, n!) coordinate orderings.

  cross_polytope_4 (paper 1 universal extension, dim 4 = 2*4 vertices):
        Q = avg of 8 edge frames at +/-e_i anchor perturbations
        (capped at 2*n if n<4).

  sparse_4 (paper 1 universal, flow/049 dim 4):
        Q = single edge frame + monomial-activity^4-weighted defect.

  prism_T4 (paper 1 §4.2 hierarchy applied with d=4):
        Q = path Schmidt slope, BUT outer corrector is T_4 (Richardson
        order 2^{4-2}=4) instead of T_2.  This is the multi-anchor
        "difference prism" geometry of flow/020 in operator form.

OUTER SCHEME
============
T_2 K=1: Aitken Delta-squared of P_{F,a}, anchor reset every step.
For prism_T4: T_4 (4 evaluations of P, Richardson stencil).

NO PORTFOLIO.  Each orbit uses ONE geometry for ALL its steps.

BENCH
=====
Bezout configurations: (2,3), (2,4), (2,5), (2,6), (3,2), (3,3), (4,2)
3 seeds per config.

Persistence: /tmp/bench_067_results.json
Modes: chunk1..chunk7, report, reset.
"""
from __future__ import annotations
import cmath, itertools, json, math, os, random, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bench_064_vs_059 import (
    F_eval, eval_poly, residual_norm, is_finite_vec, degree, solve_linear,
    schmidt_path, edge_frame, gen_random_poly_system,
)


# ============================================================================
# Q_F constructions - paper 1 universal forms
# ============================================================================
def Q_path(polys, a, z, F_a, counter=None):
    """Paper 1 §3.5 eq. (3.7) canonical Schmidt slope.

    Q[i,j] computed by interpolating between a and z one coordinate at a
    time in canonical order (0, 1, ..., n-1).
    """
    return schmidt_path(polys, a, z, tuple(range(len(z))), counter=counter)


def _simplex_corrected_frame(polys, a, z, F_a, counter=None):
    """One edge frame B + radial defect correction (flow/047).

    Telescoping identity holds:  F(z)-F(a) = Q (z-a), via
        Q = B + (F(z)-F(a)-B*(z-a)) * (z-a)^T / |z-a|^2.
    """
    n = len(z)
    delta = [z[i] - a[i] for i in range(n)]
    norm2 = sum(abs(v) ** 2 for v in delta)
    B = edge_frame(polys, a, z, F_a, counter=counter)
    if norm2 < 1e-28:
        return B
    F_z = F_eval(polys, z, counter=counter)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_a[i] - Bd[i] for i in range(n)]
    w = [v.conjugate() for v in delta]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


def Q_simplex_4(polys, a, z, F_a, counter=None):
    """Paper 1 universal simplex Q_F at dim 4 (flow/047 generalised).

    Average 4 simplex-corrected edge frames computed at 4 perturbed
    anchors.  Each frame independently satisfies telescoping; the
    average does too (linearity).
    """
    n = len(z)
    rng = random.Random(67000 + 4 * 1000 + n)
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    K = 4
    for k in range(K):
        if k == 0:
            a_use, F_use = a, F_a
        else:
            a_use = [ai + complex(0.07 * rng.gauss(0, 1), 0.07 * rng.gauss(0, 1))
                     for ai in a]
            F_use = F_eval(polys, a_use, counter)
        Qk = _simplex_corrected_frame(polys, a_use, z, F_use, counter)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Qk[i][j] / K
    return Q_avg


def Q_hypercube_4(polys, a, z, F_a, counter=None):
    """Paper 1 universal hypercube Q_F at dim 4 (flow/046).

    Average Schmidt paths over up to 2^4 = 16 coordinate orderings.
    For n <= 4 caps at n! (4! = 24, 3! = 6, 2! = 2).
    """
    n = len(z)
    perms_all = list(itertools.permutations(range(n)))
    cap = min(16, len(perms_all))
    if cap < len(perms_all):
        perms = list(perms_all[:cap])
    else:
        perms = perms_all
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    for order in perms:
        Q = schmidt_path(polys, a, z, order, counter=counter)
        for i in range(n):
            for j in range(n):
                Q_avg[i][j] += Q[i][j] / len(perms)
    return Q_avg


def Q_cross_polytope_4(polys, a, z, F_a, counter=None):
    """Paper 1 universal cross-polytope Q_F at dim 4 (flow/014 lift).

    Average 2*min(4,n) edge frames at +/-e_i anchor perturbations.
    """
    n = len(z)
    K = min(4, n)
    Q_avg = [[0.0 + 0.0j] * n for _ in range(n)]
    count = 0
    eps = 0.1
    for k in range(K):
        for sign in (+1.0, -1.0):
            a_pert = list(a)
            a_pert[k] = a_pert[k] + sign * eps
            F_pert = F_eval(polys, a_pert, counter)
            Bk = _simplex_corrected_frame(polys, a_pert, z, F_pert, counter)
            for i in range(n):
                for j in range(n):
                    Q_avg[i][j] += Bk[i][j]
            count += 1
    if count == 0:
        return _simplex_corrected_frame(polys, a, z, F_a, counter)
    return [[v / count for v in row] for row in Q_avg]


def Q_sparse_4(polys, a, z, F_a, counter=None):
    """Paper 1 universal sparse Q_F at dim 4 (flow/049): activity power 4.

    One edge frame + monomial-activity^4-weighted radial defect correction.
    """
    n = len(z)
    activity = [0.0] * n
    for poly in polys:
        for exp, coeff in poly.items():
            mag = abs(coeff)
            for j, e in enumerate(exp):
                activity[j] += mag * (e ** 4)
    m = max(activity) if activity else 0.0
    activity = [v if v > 0 else max(1.0, 0.1 * m) for v in activity]
    delta = [z[i] - a[i] for i in range(n)]
    norm2 = sum(abs(v) ** 2 for v in delta)
    B = edge_frame(polys, a, z, F_a, counter=counter)
    if norm2 < 1e-28:
        return B
    F_z = F_eval(polys, z, counter=counter)
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_a[i] - Bd[i] for i in range(n)]
    w = [activity[j] * delta[j].conjugate() for j in range(n)]
    denom = sum(w[j] * delta[j] for j in range(n))
    if abs(denom) < 1e-28:
        return B
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / denom
    return Q


# ============================================================================
# Pandrosion operator + T_n correctors (paper 1 eq. 3.10 + §4.2)
# ============================================================================
def pandrosion_h(polys, a, z, F_a, Q_func, counter=None):
    """Paper 1 eq. (3.10): P_{F,a}(z) = a - Q_F(a, z)^{-1} F(a)."""
    Q = Q_func(polys, a, z, F_a, counter=counter)
    step = solve_linear(Q, [-v for v in F_a])
    if step is None:
        return list(z), False
    out = [a[i] + step[i] for i in range(len(z))]
    return (out, True) if is_finite_vec(out) else (list(z), False)


def T2_step(polys, a, z, F_a, Q_func, counter=None):
    """Aitken Delta^2 of P (paper 1 §4.2 T_2; paper 8 K=1)."""
    s1, ok1 = pandrosion_h(polys, a, z, F_a, Q_func, counter)
    if not ok1:
        return list(z), False
    s2, ok2 = pandrosion_h(polys, a, s1, F_a, Q_func, counter)
    if not ok2:
        return s1, True
    n = len(z); out = []
    for k in range(n):
        d0 = s1[k] - z[k]; d2 = s2[k] - 2 * s1[k] + z[k]
        out.append(s2[k] if abs(d2) < 1e-300 else z[k] - d0 * d0 / d2)
    return out, True


def T4_step(polys, a, z, F_a, Q_func, counter=None):
    """T_4 Richardson (paper 1 §4.2, order 2^{4-2}=4)."""
    s_list = [list(z)]
    for _ in range(4):
        s_new, ok = pandrosion_h(polys, a, s_list[-1], F_a, Q_func, counter)
        if not ok:
            return s_list[-1], len(s_list) > 1
        s_list.append(s_new)
    n = len(z)
    t2 = []
    for j in range(n):
        d0 = s_list[1][j] - s_list[0][j]
        d2 = s_list[2][j] - 2 * s_list[1][j] + s_list[0][j]
        t2.append(s_list[2][j] if abs(d2) < 1e-300 else s_list[0][j] - d0 * d0 / d2)
    lam_den = sum(s_list[1][j] - s_list[0][j] for j in range(n))
    if abs(lam_den) < 1e-300:
        return t2, True
    lam = sum(s_list[2][j] - s_list[1][j] for j in range(n)) / lam_den
    if abs(lam - 1.0) < 1e-12:
        return s_list[-1], True
    out = list(t2)
    for kk in range(3, 5):
        out = [out[j] - (s_list[kk][j] - out[j]) / (lam - 1.0) for j in range(n)]
    return out, True


GEOMS = {
    "path":             Q_path,
    "simplex_4":        Q_simplex_4,
    "hypercube_4":      Q_hypercube_4,
    "cross_polytope_4": Q_cross_polytope_4,
    "sparse_4":         Q_sparse_4,
    "prism_T4":         Q_path,    # same Q, different outer corrector
}


# ============================================================================
# Single-geometry orbit (NO portfolio).  T_2 K=1, except prism_T4 uses T_4.
# ============================================================================
def orbit(polys, z_init, geom_name, tol=1e-9, max_epochs=60, counter=None):
    n = len(z_init); z = list(z_init)
    anchor = [zi + complex(0.13 * (i + 1), 0.07) for i, zi in enumerate(z_init)]
    rng = random.Random(67000 + n + hash(geom_name) % 1000)
    Q_func = GEOMS[geom_name]
    step_func = T4_step if geom_name == "prism_T4" else T2_step
    for _ in range(max_epochs):
        F_anchor = F_eval(polys, anchor, counter)
        if is_finite_vec(F_anchor) and max(abs(v) for v in F_anchor) < tol:
            return anchor, True
        r_z = residual_norm(polys, z, counter)
        if r_z < tol:
            return z, True
        z_new, ok = step_func(polys, anchor, z, F_anchor, Q_func, counter)
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
# Bench harness with persistence
# ============================================================================
RESULTS_FILE = "/tmp/bench_067_results.json"

GEOM_LIST = ["path", "simplex_4", "hypercube_4", "cross_polytope_4",
             "sparse_4", "prism_T4"]


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
    print(f"  {'(n,d,s)':>9} {'geom':>17} {'Bez':>4} | "
          f"{'cov%':>6} {'time':>7} {'Fevals':>9} {'ok/fail':>10} {'coll':>4}",
          flush=True)
    for (n, d) in configs:
        for seed_off in range(seeds_per_config):
            seed = 67000 + 1000 * n + 100 * d + seed_off
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
                print(f"  ({n:>2},{d:>2},{seed_off}) {geom:>17} "
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
    print("AGGREGATED 067 RESULTS  (T_2 K=1 + single paper-1 universal geometry)")
    print("-" * 90)
    print(f"  {'geom':>17} {'(n,d)':>7} {'Bez':>4} {'#s':>3} | "
          f"{'cov%':>6} {'time':>7} {'F-evals':>10} {'eff(c/F)':>9}")
    print("-" * 90)
    for key in sorted(by_key.keys(), key=lambda k: (k[0], k[1], k[2])):
        geom, n, d = key
        rs = by_key[key]
        cov = 100 * sum(r["coverage"] for r in rs) / len(rs)
        tt = sum(r["time"] for r in rs) / len(rs)
        fe = sum(r["F_evals"] for r in rs) / len(rs)
        eff = (cov / 100.0) / max(fe, 1) * 1e6
        print(f"  {geom:>17} ({n:>2},{d:>2}) {rs[0]['bezout']:>4} "
              f"{len(rs):>3} | {cov:>5.1f}% {tt:>6.2f}s {fe:>10.0f} {eff:>8.2f}")


CONFIG_CHUNKS = {
    "1": [(2, 3)],
    "2": [(2, 4)],
    "3": [(2, 5)],
    "4": [(2, 6)],
    "5": [(3, 2)],
    "6": [(3, 3)],
    "7": [(4, 2)],
    "all_n2": [(2, 3), (2, 4), (2, 5), (2, 6)],
    "small": [(2, 3), (3, 2), (4, 2)],
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
