"""
PAPER: 183 (NEW — Bessel/cosh PF mixture probe)
TITLE: Positive mixtures of the factorial-Bessel PF kernel:
       a quadrature route toward the Xi PF∞ door
STATUS: RH remains OPEN.  This paper tests a tempting closure phenomenon
        suggested by paper 182: positive mixtures of the sequences
            b_n(u) = u^(2n)/(2n)!
        do NOT preserve Toeplitz total positivity in general.  However,
        the specific Riemann-Phi quadrature mixture survives the tested
        finite windows.  The target is therefore special structure of the
        Riemann-Phi measure, not a generic convexity theorem.
DEPENDS: 181 (PF∞ door), 182 (moment/factorial decomposition).

THEORY
======

------------------------------------------------------------------------
THE SPECIAL PF KERNEL
------------------------------------------------------------------------

For each fixed u > 0, define

    b_n(u) = u^(2n)/(2n)!.

Its generating function is

    B_u(z) = sum_n b_n(u) z^n = cosh(u sqrt(z)).

The zeros of B_u are

    z_k = -pi^2 (k + 1/2)^2 / u^2,      k = 0,1,2,...

so B_u is Laguerre-Pólya with all zeros real and non-positive.
Therefore b(u) is PF∞ for each fixed u.

Paper 182 gives

    A_n = integral Phi(u) b_n(u) du.

The hard question is not whether each b(u) is PF∞; it is whether this
positive mixture remains PF∞ for the particular Riemann-Phi measure.

------------------------------------------------------------------------
WHAT THIS PAPER TESTS
------------------------------------------------------------------------

  (1) deterministic two/three-point mixtures;
  (2) random positive mixtures over u in [0.1, 4], which expose failures;
  (3) Riemann-Phi quadrature mixtures;
  (4) comparison against a deliberately wrong kernel whose mixtures fail.

The goal is to identify whether a generic closure theorem is plausible.
It is not.  The proof must exploit the Riemann-Phi weight.
"""
from __future__ import annotations

import itertools
import math
import random


def det_small(rows):
    n = len(rows)
    if n == 0:
        return 1.0
    total = 0.0
    for perm in itertools.permutations(range(n)):
        inv = 0
        for i in range(n):
            for j in range(i + 1, n):
                inv += int(perm[i] > perm[j])
        term = 1.0
        for i in range(n):
            term *= rows[i][perm[i]]
        total += -term if inv % 2 else term
    return total


def toeplitz_entry(seq, r, c):
    k = c - r
    if k < 0 or k >= len(seq):
        return 0.0
    return seq[k]


def min_toeplitz_minor(seq, window=8, max_size=4, tol=1e-18):
    total = 0
    negative = 0
    min_det = float("inf")
    min_case = None
    indices = list(range(window))
    for size in range(1, max_size + 1):
        for rows_idx in itertools.combinations(indices, size):
            for cols_idx in itertools.combinations(indices, size):
                M = [[toeplitz_entry(seq, r, c) for c in cols_idx] for r in rows_idx]
                d = det_small(M)
                total += 1
                if d < min_det:
                    min_det = d
                    min_case = (rows_idx, cols_idx, size)
                if d < -tol:
                    negative += 1
    return total, negative, min_det, min_case


def bessel_seq_mixture(us, ws, n_terms=18):
    return [
        sum(w * (u ** (2 * n)) / math.factorial(2 * n) for u, w in zip(us, ws))
        for n in range(n_terms)
    ]


def bad_kernel_mixture(us, ws, n_terms=18):
    """A nearby-looking kernel with no LP/PF reason: u^(2n)/(n!)^2."""
    return [
        sum(w * (u ** (2 * n)) / (math.factorial(n) ** 2) for u, w in zip(us, ws))
        for n in range(n_terms)
    ]


def phi_riemann(u, n_terms=40):
    s = 0.0
    e4u = math.exp(4 * u)
    e5u = math.exp(5 * u)
    e9u = math.exp(9 * u)
    for n in range(1, n_terms + 1):
        nn = n * n
        term = (2 * math.pi**2 * n**4 * e9u - 3 * math.pi * nn * e5u)
        term *= math.exp(-math.pi * nn * e4u)
        s += term
    return max(0.0, s)


def phi_quadrature_nodes(U=5.0, n_nodes=80):
    """Midpoint quadrature nodes and positive weights for Phi(u)du."""
    h = U / n_nodes
    us = [(j + 0.5) * h for j in range(n_nodes)]
    ws = [phi_riemann(u) * h for u in us]
    total = sum(ws)
    if total > 0:
        ws = [w / total for w in ws]
    return us, ws


def main():
    print("=" * 80)
    print("PAPER 183 — Bessel/cosh PF mixture probe")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Fixed kernel sanity: zeros of cosh(u sqrt(z))
    # ------------------------------------------------------------------
    print("\n[1] Fixed-u kernel sanity")
    for u in [0.5, 1.0, 2.0]:
        first_zero = -(math.pi * 0.5 / u) ** 2
        print(f"  u={u:.2f}: first zero of cosh(u sqrt(z)) = {first_zero:.8f}")
    print("  all zeros are negative, hence each fixed-u sequence is PF∞")

    # ------------------------------------------------------------------
    # [2] Deterministic mixtures
    # ------------------------------------------------------------------
    print("\n[2] Deterministic positive mixtures")
    examples = [
        ([0.5, 2.0], [0.5, 0.5]),
        ([0.2, 1.3], [0.35, 0.65]),
        ([1.0, 3.0], [0.8, 0.2]),
        ([0.3, 0.9, 2.4], [0.2, 0.5, 0.3]),
    ]
    print(f"  {'nodes u':>22} {'negative minors':>17} {'min det':>14}")
    for us, ws in examples:
        seq = bessel_seq_mixture(us, ws)
        _, neg, md, _ = min_toeplitz_minor(seq)
        print(f"  {str(us):>22} {neg:>17} {md:>14.3e}")

    # ------------------------------------------------------------------
    # [3] Random mixture stress test
    # ------------------------------------------------------------------
    print("\n[3] Random positive mixture stress test")
    rng = random.Random(2026)
    failures = 0
    worst = float("inf")
    worst_case = None
    n_trials = 200
    for _ in range(n_trials):
        k = rng.randint(2, 7)
        us = sorted(0.1 + 3.9 * rng.random() for _ in range(k))
        raw = [rng.random() for _ in range(k)]
        s = sum(raw)
        ws = [x / s for x in raw]
        seq = bessel_seq_mixture(us, ws)
        _, neg, md, case = min_toeplitz_minor(seq, window=8, max_size=4)
        if neg:
            failures += 1
        if md < worst:
            worst = md
            worst_case = (us, ws, case)
    print(f"  trials: {n_trials}")
    print(f"  failures: {failures}")
    print(f"  worst min determinant: {worst:.3e}")
    print(f"  worst case minor: {worst_case[2]}")

    # ------------------------------------------------------------------
    # [4] Riemann-Phi quadrature mixture
    # ------------------------------------------------------------------
    print("\n[4] Riemann-Phi quadrature mixture")
    us, ws = phi_quadrature_nodes()
    seq = bessel_seq_mixture(us, ws)
    total, neg, md, case = min_toeplitz_minor(seq, window=8, max_size=4)
    print(f"  quadrature nodes: {len(us)}")
    print(f"  Toeplitz minors tested: {total}")
    print(f"  negative minors: {neg}")
    print(f"  minimum determinant: {md:.3e}")
    print(f"  minimum case: {case}")

    # ------------------------------------------------------------------
    # [5] Control experiment
    # ------------------------------------------------------------------
    print("\n[5] Control experiment with a nearby wrong kernel")
    print("    bad kernel: u^(2n)/(n!)^2.  It also has positive coefficients,")
    print("    but it is not the Riemann/cosh kernel.")
    bad_failures = 0
    bad_worst = float("inf")
    for _ in range(n_trials):
        k = rng.randint(2, 7)
        us = sorted(0.1 + 3.9 * rng.random() for _ in range(k))
        raw = [rng.random() for _ in range(k)]
        s = sum(raw)
        ws = [x / s for x in raw]
        seq = bad_kernel_mixture(us, ws)
        _, neg, md, _ = min_toeplitz_minor(seq, window=8, max_size=4)
        if neg:
            bad_failures += 1
        bad_worst = min(bad_worst, md)
    print(f"  random bad-kernel failures: {bad_failures}/{n_trials}")
    print(f"  bad-kernel worst min determinant: {bad_worst:.3e}")

    # ------------------------------------------------------------------
    # [6] Honest assessment
    # ------------------------------------------------------------------
    print("\n[6] HONEST ASSESSMENT")
    print("  WHAT SURVIVED:")
    print("    - simple deterministic mixtures survived the tested windows;")
    print("    - arbitrary random positive mixtures did NOT survive;")
    print("    - the Riemann-Phi quadrature mixture survived the tested window.")
    print()
    print("  WHAT THIS SUGGESTS:")
    print("    The missing theorem cannot be generic convexity of the PF∞ cone")
    print("    along the Bessel-cosh ray.  It must use special regularity,")
    print("    log-concavity, or variation-diminishing structure of Phi(u)du.")
    print()
    print("  WHAT WOULD MATTER:")
    print("    Prove PF∞ for the specific Riemann-Phi measure:")
    print("      A_n = integral Phi(u) u^(2n)/(2n)! du.")
    print("    Paper 183 rules out the stronger false shortcut 'all positive")
    print("    measures work'.")


if __name__ == "__main__":
    main()
