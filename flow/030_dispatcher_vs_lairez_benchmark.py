"""
PAPER: 030
TITLE: Complete benchmark: 029 geometry dispatcher vs Lairez
STATUS: benchmark of hypercube/simplex/sparse Pandrosion geometries

MISSION
=======

flow/029 introduced a geometry dispatcher:

    box support       -> hypercube Pandrosion
    simplex support   -> simplex Pandrosion
    sparse support    -> exact-support S_A Pandrosion

This file benchmarks that dispatcher against a generic Lairez-style
predictor-corrector homotopy on the same family of polynomial systems:

    P_i(s) = x_i (1-s_i) S_A(s) - (x_i - 1) = 0.

For the Pandrosion side we run the dispatcher with prism dimensions d=3,4,5
and keep the cheapest result.

For Lairez we track x(t)=x_start+t(x_target-x_start), starting from a
symmetric root at x_start=(1.25,...,1.25).  Cost is weighted as

    weighted = F_evals + n * J_evals.

This is generous to Lairez but keeps the comparison consistent with flow/024,
025, and 026.

HONEST EXPECTATION
==================

The dispatcher should win on systems whose support matches a Pandrosion
geometry and whose fixed-point map is contractive.  It is not a theorem that it
wins on every possible S_A.  This benchmark measures that boundary.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_flow(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DISPATCH = load_flow("flow029_dispatcher", "029_sparse_newton_polytope_dispatcher.py")


# -----------------------------------------------------------------------------
# S_A system and Lairez-style homotopy.
# -----------------------------------------------------------------------------
def S_A(s, A):
    total = 0.0
    for alpha in A:
        term = 1.0
        for si, ai in zip(s, alpha):
            term *= si ** ai
        total += term
    return total


def dS_A(s, A, j):
    total = 0.0
    for alpha in A:
        if alpha[j] == 0:
            continue
        term = alpha[j] * (s[j] ** (alpha[j] - 1))
        for k, ak in enumerate(alpha):
            if k != j:
                term *= s[k] ** ak
        total += term
    return total


def F_target(s, x_vec, A):
    S = S_A(s, A)
    return [x_vec[i] * (1.0 - s[i]) * S - (x_vec[i] - 1.0) for i in range(len(s))]


def jacobian(s, x_vec, A):
    n = len(s)
    S = S_A(s, A)
    J = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            term = x_vec[i] * (1.0 - s[i]) * dS_A(s, A, j)
            if i == j:
                term -= x_vec[i] * S
            J[i][j] = term
    return J


def solve_linear(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-300:
            return None
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        rhs = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = rhs / M[i][i]
    return x


def find_symmetric_start(x_const, n, A, tol=1e-15):
    s = 0.5
    for _ in range(1000):
        S = S_A([s] * n, A)
        sn = 1.0 - (x_const - 1.0) / (x_const * S)
        if abs(sn - s) < tol:
            return sn
        s = sn
    return s


def lairez_track(x_start, x_target, A, s0, tol=1e-13, max_steps=1000):
    n = len(x_target)
    s = list(s0)
    t = 0.0
    dt = 0.05
    dt_min = 1e-8
    dt_max = 0.5
    dx = [x_target[i] - x_start[i] for i in range(n)]
    F_evals = 0
    J_evals = 0
    steps = 0

    while t < 1.0 - 1e-15 and steps < max_steps:
        if t + dt > 1.0:
            dt = 1.0 - t
        x_t = [x_start[i] + t * dx[i] for i in range(n)]
        S = S_A(s, A)
        dH_dt = [dx[i] * ((1.0 - s[i]) * S - 1.0) for i in range(n)]
        J = jacobian(s, x_t, A)
        J_evals += 1
        dsdt = solve_linear(J, [-v for v in dH_dt])
        if dsdt is None:
            dt = max(dt / 2.0, dt_min)
            continue

        t_new = t + dt
        s_curr = [s[i] + dt * dsdt[i] for i in range(n)]
        x_new = [x_start[i] + t_new * dx[i] for i in range(n)]
        newton_used = 0
        failed = False

        for _ in range(8):
            Fv = F_target(s_curr, x_new, A)
            F_evals += 1
            if max(abs(v) for v in Fv) < tol:
                break
            Jc = jacobian(s_curr, x_new, A)
            J_evals += 1
            ds = solve_linear(Jc, [-v for v in Fv])
            if ds is None:
                failed = True
                break
            s_curr = [s_curr[i] + ds[i] for i in range(n)]
            newton_used += 1
            if max(abs(v) for v in ds) < 1e-12:
                break

        res = max(abs(v) for v in F_target(s_curr, x_new, A))
        F_evals += 1
        if failed or res > 1e-7 or newton_used > 5:
            dt = max(dt / 2.0, dt_min)
            if dt > dt_min:
                continue

        s = s_curr
        t = t_new
        steps += 1
        if newton_used <= 2:
            dt = min(dt * 1.5, dt_max)

    for _ in range(12):
        Fv = F_target(s, x_target, A)
        F_evals += 1
        if max(abs(v) for v in Fv) < tol:
            break
        Jc = jacobian(s, x_target, A)
        J_evals += 1
        ds = solve_linear(Jc, [-v for v in Fv])
        if ds is None:
            break
        s = [s[i] + ds[i] for i in range(n)]
    return s, F_evals, J_evals, steps


# -----------------------------------------------------------------------------
# Benchmark cases.
# -----------------------------------------------------------------------------
def primes(n):
    pool = [2.0, 3.0, 5.0, 7.0, 11.0, 13.0]
    return pool[:n]


def all_twos(n):
    return [2.0] * n


def partial(n):
    if n == 2:
        return [2.0, 5.0]
    return [2.0] * (n - 1) + [5.0]


def sparse_cycle_A(n):
    A = {tuple([0] * n)}
    for i in range(n):
        e = [0] * n
        e[i] = 1
        A.add(tuple(e))
    for i in range(n - 1):
        e = [0] * n
        e[i] = 1
        e[i + 1] = 1
        A.add(tuple(e))
    return sorted(A)


def sparse_diag_A(n):
    A = {tuple([0] * n)}
    for i in range(n):
        e = [0] * n
        e[i] = 1
        A.add(tuple(e))
    A.add(tuple([1] * n))
    return sorted(A)


def thin_chain_A(n):
    A = {tuple([0] * n)}
    current = [0] * n
    for i in range(n):
        current[i] = 1
        A.add(tuple(current))
    return sorted(A)


def benchmark_cases():
    cases = []
    for n in range(2, 7):
        cases.append(DISPATCH.build_system_from_A(f"box n{n} sym", all_twos(n), DISPATCH.box_support(n, 2)))
        cases.append(DISPATCH.build_system_from_A(f"box n{n} full", primes(n), DISPATCH.box_support(n, 2)))
        cases.append(DISPATCH.build_system_from_A(f"simplex n{n} sym", all_twos(n), DISPATCH.compositions_leq(n, 1)))
        cases.append(DISPATCH.build_system_from_A(f"simplex n{n} full", primes(n), DISPATCH.compositions_leq(n, 1)))
        cases.append(DISPATCH.build_system_from_A(f"sparse-cycle n{n}", primes(n), sparse_cycle_A(n)))
        cases.append(DISPATCH.build_system_from_A(f"sparse-diag n{n}", partial(n), sparse_diag_A(n)))
        cases.append(DISPATCH.build_system_from_A(f"thin-chain n{n}", primes(n), thin_chain_A(n)))
    for n in range(2, 5):
        cases.append(DISPATCH.build_system_from_A(f"box p3 n{n}", partial(n), DISPATCH.box_support(n, 3)))
        cases.append(DISPATCH.build_system_from_A(f"simplex p3 n{n}", partial(n), DISPATCH.compositions_leq(n, 2)))
    return cases


def run_case(system):
    x_vec = DISPATCH.infer_x_vec(system)
    A = set(system["A"])
    n = system["n"]

    prism_results = []
    for d in [3, 4, 5]:
        result = DISPATCH.dispatch_solve(system, d=d)
        prism_results.append((d, result["evals"], result["iters"], result["family"], result["res"]))
    best_prism = min(prism_results, key=lambda item: item[1])

    x_start_const = 1.25
    x_start = [x_start_const] * n
    s0_scalar = find_symmetric_start(x_start_const, n, A)
    s0 = [s0_scalar] * n
    s_l, Fe, Je, steps = lairez_track(x_start, x_vec, A, s0)
    weighted = Fe + n * Je
    l_res = max(abs(v) for v in F_target(s_l, x_vec, A))

    best = "prism" if best_prism[1] <= weighted else "lai"
    speedup = weighted / best_prism[1] if best_prism[1] else float("inf")
    return {
        "name": system["name"],
        "n": n,
        "A_size": len(A),
        "family": best_prism[3],
        "d": best_prism[0],
        "p_evals": best_prism[1],
        "p_iters": best_prism[2],
        "p_res": best_prism[4],
        "l_steps": steps,
        "l_Fe": Fe,
        "l_Je": Je,
        "l_weighted": weighted,
        "l_res": l_res,
        "best": best,
        "speedup": speedup,
    }


def median(values):
    vals = sorted(values)
    return vals[len(vals) // 2]


def main():
    print("=" * 112)
    print("flow/030 -- complete benchmark: 029 dispatcher vs Lairez")
    print("=" * 112)
    print()
    print("Pandrosion side: dispatcher from flow/029, best of prism d=3/4/5.")
    print("Lairez side    : predictor-corrector homotopy, weighted cost F+nJ.")
    print()

    rows = [run_case(system) for system in benchmark_cases()]

    print(
        f"{'system':>20} {'fam':>9} {'n':>2} {'|A|':>4} {'d':>2} | "
        f"{'prism':>10} {'Lairez weighted':>18} {'speed':>8} {'best':>7}"
    )
    print("-" * 112)
    for row in rows:
        print(
            f"{row['name'][:20]:>20} {row['family']:>9} {row['n']:>2} "
            f"{row['A_size']:>4} {row['d']:>2} | "
            f"{row['p_iters']:>3}it/{row['p_evals']:<4} "
            f"{row['l_steps']:>3}st/{row['l_Fe']}F+{row['l_Je']}J={row['l_weighted']:<5} "
            f"{row['speedup']:>7.1f}x {row['best']:>7}"
        )
        if row["p_res"] > 1e-8 or row["l_res"] > 1e-8:
            print(f"  WARNING residuals: prism={row['p_res']:.2e}, lairez={row['l_res']:.2e}")

    wins = {}
    family_stats = {}
    speedups = []
    for row in rows:
        wins[row["best"]] = wins.get(row["best"], 0) + 1
        family_stats.setdefault(row["family"], []).append(row)
        speedups.append(row["speedup"])

    print()
    print("=" * 112)
    print("SUMMARY")
    print("=" * 112)
    print(f"cases                  : {len(rows)}")
    print(f"winners                : {wins}")
    print(f"speedup min/median/max : {min(speedups):.1f}x / {median(speedups):.1f}x / {max(speedups):.1f}x")
    for family, subset in family_stats.items():
        local_wins = {}
        local_speedups = [row["speedup"] for row in subset]
        for row in subset:
            local_wins[row["best"]] = local_wins.get(row["best"], 0) + 1
        print(
            f"{family:>9}: {len(subset):>2} cases, winners={local_wins}, "
            f"median speedup={median(local_speedups):.1f}x"
        )

    print()
    print("Verdict:")
    print("  This benchmark tests the optimistic claim directly.  If Lairez has")
    print("  zero wins above, the dispatcher wins on this structured S_A corpus.")
    print("  The honest caveat remains: the corpus is Pandrosion-structured;")
    print("  arbitrary expanded polynomial systems still need support extraction")
    print("  or lifting before this dispatcher applies.")


if __name__ == "__main__":
    main()
