"""
Extended multivariate benchmark for Pandrosion-T_3 on Smale-17-style systems.
Goes beyond (n, d) = (4, 3) up to (5, 4) and (6, 3) for the Section 3.5 of
Part I.
"""
from __future__ import annotations
import math, time
from itertools import product
import numpy as np

RNG = np.random.default_rng(20260424)


def monomial_exps(n, d):
    return [e for total in range(d + 1)
            for e in product(range(d + 1), repeat=n) if sum(e) == total]


def random_KS_system(n, d, rng):
    """Kostlan-Smale: each monomial coefficient ~ N_C(0, multinomial weight)."""
    exps = monomial_exps(n, d)
    sys = []
    for _ in range(n):
        coefs = {}
        for e in exps:
            # Kostlan weight = sqrt(d! / (e_0! * e_1! * ... * (d - sum_e)!))
            ssum = sum(e)
            mult = math.factorial(d)
            for ek in e:
                mult //= math.factorial(ek)
            mult //= math.factorial(d - ssum)
            std = math.sqrt(mult)
            coefs[e] = (rng.standard_normal() + 1j * rng.standard_normal()) * std / math.sqrt(2)
        sys.append(coefs)
    return sys


def eval_poly(coefs, z):
    acc = 0j
    for e, c in coefs.items():
        t = c
        for k, ek in enumerate(e):
            if ek:
                t *= z[k] ** ek
        acc += t
    return acc


def eval_system(system, z):
    return np.array([eval_poly(coefs, z) for coefs in system], dtype=complex)


def diff_quotient_matrix(system, z0, z, n):
    Q = np.zeros((n, n), dtype=complex)
    z_curr = np.array(z0, dtype=complex)
    for j in range(n):
        delta = z[j] - z0[j]
        z_next = z_curr.copy()
        z_next[j] = z[j]
        if abs(delta) < 1e-13:
            for i, coefs in enumerate(system):
                d_acc = 0j
                for e, c in coefs.items():
                    if e[j] == 0: continue
                    t = c * e[j]
                    for m, em in enumerate(e):
                        if m == j:
                            if em - 1 > 0: t *= z_curr[m] ** (em - 1)
                        elif em > 0:
                            t *= z_curr[m] ** em
                    d_acc += t
                Q[i, j] = d_acc
        else:
            F_next = eval_system(system, z_next)
            F_curr = eval_system(system, z_curr)
            Q[:, j] = (F_next - F_curr) / delta
        z_curr = z_next
    return Q


def adaptive_pandrosion_solve(system, n, z_init, K=3, max_iters=200, tol=1e-8):
    z0 = np.array(z_init, dtype=complex)
    z = z0.copy()
    for it in range(max_iters):
        r = eval_system(system, z)
        res = float(np.max(np.abs(r)))
        if res < tol:
            return z, it, res
        Q = diff_quotient_matrix(system, z0, z, n)
        F0 = eval_system(system, z0)
        try:
            delta = np.linalg.solve(Q, F0)
        except np.linalg.LinAlgError:
            return None, it, float('inf')
        z_new = z0 - delta
        if not np.all(np.isfinite(z_new)):
            return None, it, float('inf')
        z = z_new
        if (it + 1) % K == 0:
            z0 = z.copy()
    return z, max_iters, float(np.max(np.abs(eval_system(system, z))))


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    t0 = time.perf_counter()
    print("=" * 95, flush=True)
    print("Extended multivariate Pandrosion-T_3 (K=3) benchmark on Kostlan-Smale systems", flush=True)
    print("=" * 95, flush=True)
    print(f"{'(n,d)':>7} {'Bezout':>7} {'#monoms':>9} {'trials':>7} "
          f"{'success':>8} {'avg iters':>10} {'avg time/run':>14}", flush=True)
    print("-" * 95, flush=True)

    CONFIGS = [(2, 3), (3, 3), (4, 2), (4, 3), (5, 2), (5, 3), (5, 4), (6, 2), (6, 3)]
    NUM_TRIALS = 5

    summary_rows = []
    for n, d in CONFIGS:
        bezout = d ** n
        n_monoms = len(monomial_exps(n, d))
        rng = np.random.default_rng(42)
        n_ok = 0; iters_sum = 0; time_sum = 0.0
        for _ in range(NUM_TRIALS):
            system = random_KS_system(n, d, rng)
            z0 = rng.standard_normal(n) + 1j * rng.standard_normal(n)
            t_start = time.perf_counter()
            z, it, res = adaptive_pandrosion_solve(system, n, z0, K=3,
                                                   max_iters=300)
            dt = time.perf_counter() - t_start
            time_sum += dt
            if z is not None and res < 1e-6:
                n_ok += 1; iters_sum += it
        avg_it = iters_sum / max(1, n_ok)
        avg_time = time_sum / NUM_TRIALS
        print(f"({n},{d})  {bezout:>7} {n_monoms:>9} {NUM_TRIALS:>7} "
              f"{n_ok}/{NUM_TRIALS}    {avg_it:>10.1f} {avg_time:>14.3f}s",
              flush=True)
        summary_rows.append((n, d, bezout, n_ok, NUM_TRIALS, avg_it))

    print("\n" + "=" * 95, flush=True)
    print("Summary for Section 3.5 of Part I", flush=True)
    print("=" * 95, flush=True)
    largest = max(summary_rows, key=lambda r: r[2])
    print(f"  Largest configuration tested: (n,d) = ({largest[0]}, {largest[1]}), "
          f"Bezout degree D = {largest[2]}", flush=True)
    print(f"  Convergence rates: ", end="", flush=True)
    print(", ".join(f"({n},{d}): {ok}/{tot}" for n, d, _, ok, tot, _ in summary_rows),
          flush=True)
    overall = sum(r[3] for r in summary_rows)
    overall_total = sum(r[4] for r in summary_rows)
    print(f"  Overall: {overall}/{overall_total} = "
          f"{100 * overall / overall_total:.0f}% convergence", flush=True)
    print(f"  Mean iterations across all configs: "
          f"{sum(r[5] for r in summary_rows)/len(summary_rows):.1f}", flush=True)
    print(f"\n[total elapsed: {time.perf_counter() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
