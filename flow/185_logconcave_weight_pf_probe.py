"""
PAPER: 185 (NEW — log-concave weight hypothesis)
TITLE: Log-concave weights as the missing structure in the Xi PF∞ door
STATUS: RH remains OPEN.  Papers 183-184 show that arbitrary positive
        mixtures of b_n(u)=u^(2n)/(2n)! can fail PF∞, while Riemann-Phi
        survives finite tests and is numerically log-concave on its effective
        support.  This paper tests whether log-concavity of the mixing
        weights explains the survival.  It does not: log-concavity helps
        but is insufficient.
DEPENDS: 181 (PF∞ door), 183 (mixture failures), 184 (Phi stability).

THEORY
======

Given an ordered grid u_0 < ... < u_m and positive weights w_j, define

    a_n = sum_j w_j u_j^(2n)/(2n)!.

Arbitrary positive weights do not preserve Toeplitz PF∞.  The Riemann-Phi
weights appear special.  A plausible missing hypothesis is discrete
log-concavity:

    w_j^2 >= w_{j-1} w_{j+1}.

If log-concave weights preserved the tested PF minors while arbitrary
weights failed frequently, the analytic target would become:

    Prove Phi(u) is in a continuous log-concave / PF-preserving class.

The experiment finds failures even inside this log-concave model.  The
target must therefore be stronger: a variation-diminishing or PF property
of the weight itself, not ordinary log-concavity.
"""
from __future__ import annotations

import itertools
import math
import random


def det_small(rows):
    n = len(rows)
    total = 0.0
    for perm in itertools.permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
        term = 1.0
        for i in range(n):
            term *= rows[i][perm[i]]
        total += -term if inv % 2 else term
    return total


def seq_from_weights(us, ws, n_terms=16):
    return [
        sum(w * (u ** (2 * n)) / math.factorial(2 * n) for u, w in zip(us, ws))
        for n in range(n_terms)
    ]


def toeplitz_entry(seq, r, c):
    k = c - r
    return seq[k] if 0 <= k < len(seq) else 0.0


def min_toeplitz(seq, window=7, max_size=4, tol=1e-14):
    neg = 0
    total = 0
    min_det = float("inf")
    min_case = None
    idx = range(window)
    for size in range(1, max_size + 1):
        for rows in itertools.combinations(idx, size):
            for cols in itertools.combinations(idx, size):
                M = [[toeplitz_entry(seq, r, c) for c in cols] for r in rows]
                d = det_small(M)
                total += 1
                if d < min_det:
                    min_det = d
                    min_case = (rows, cols, size)
                if d < -tol:
                    neg += 1
    return total, neg, min_det, min_case


def normalize(ws):
    s = sum(ws)
    return [w / s for w in ws]


def random_logconcave_weights(us, rng):
    """Weights exp(-a(u-m)^2 + b u), a>=0: log-concave on a line."""
    a = 0.05 + 3.0 * rng.random()
    m = min(us) + (max(us) - min(us)) * rng.random()
    b = -2.0 + 4.0 * rng.random()
    raw = [math.exp(-a * (u - m) ** 2 + b * u) for u in us]
    return normalize(raw)


def random_arbitrary_weights(us, rng):
    raw = [10 ** (-3 + 6 * rng.random()) for _ in us]
    return normalize(raw)


def is_discrete_logconcave(ws, tol=1e-15):
    return all(ws[i] ** 2 + tol >= ws[i - 1] * ws[i + 1] for i in range(1, len(ws) - 1))


def main():
    print("=" * 80)
    print("PAPER 185 — log-concave weight hypothesis")
    print("=" * 80)

    rng = random.Random(2026)
    us = [0.15 + 0.35 * j for j in range(12)]
    print("\n[1] Grid")
    print("  u_j:", " ".join(f"{u:.2f}" for u in us))

    # ------------------------------------------------------------------
    # [2] Compare arbitrary vs log-concave weights
    # ------------------------------------------------------------------
    print("\n[2] Stress test: arbitrary vs log-concave weights")
    print(f"  {'class':>16} {'trials':>8} {'failures':>10} {'worst min det':>16}")
    for label, sampler in [
        ("arbitrary", random_arbitrary_weights),
        ("log-concave", random_logconcave_weights),
    ]:
        failures = 0
        worst = float("inf")
        worst_case = None
        trials = 120
        logc_bad = 0
        for _ in range(trials):
            ws = sampler(us, rng)
            if not is_discrete_logconcave(ws):
                logc_bad += 1
            seq = seq_from_weights(us, ws)
            _, neg, md, case = min_toeplitz(seq)
            failures += int(neg > 0)
            if md < worst:
                worst = md
                worst_case = case
        print(f"  {label:>16} {trials:>8} {failures:>10} {worst:>16.6e}")
        if logc_bad:
            print(f"    non-log-concave samples in class {label}: {logc_bad}")
        print(f"    worst case minor: {worst_case}")

    # ------------------------------------------------------------------
    # [3] Boundary: two-peak weights violate log-concavity
    # ------------------------------------------------------------------
    print("\n[3] Controlled two-peak perturbation")
    print("    w = (1-eps) logconcave_base + eps(two far peaks)")
    base = random_logconcave_weights(us, random.Random(11))
    peaks = [0.0] * len(us)
    peaks[1] = 0.5
    peaks[-2] = 0.5
    print(f"  {'eps':>8} {'logconcave?':>13} {'neg minors':>12} {'min det':>14}")
    for eps in [0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0]:
        ws = normalize([(1 - eps) * b + eps * p for b, p in zip(base, peaks)])
        seq = seq_from_weights(us, ws)
        _, neg, md, _ = min_toeplitz(seq)
        print(f"  {eps:>8.2f} {str(is_discrete_logconcave(ws)):>13}"
              f" {neg:>12} {md:>14.6e}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  RESULT:")
    print("    Arbitrary weights fail often; log-concave weights fail less often")
    print("    but still fail. Ordinary log-concavity is not the missing theorem.")
    print()
    print("  INTERPRETATION:")
    print("    The needed structure is stronger than log-concavity: likely a")
    print("    continuous total-positivity / variation-diminishing property of")
    print("    the Riemann-Phi density.")
    print()
    print("  NEXT THEOREM TO SEEK:")
    print("    A criterion: if w(u) belongs to a specific PF/log-concave class,")
    print("    then a_n = integral w(u) u^(2n)/(2n)! du is PF∞.")


if __name__ == "__main__":
    main()
