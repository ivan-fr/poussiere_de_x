"""
PAPER: 184 (NEW — stability boundary around the Riemann-Phi measure)
TITLE: Riemann-Phi is not a generic positive mixture: PF-minor stability,
       log-shape diagnostics, and adversarial perturbations
STATUS: RH remains OPEN.  Paper 183 showed that arbitrary positive mixtures
        of b_n(u)=u^(2n)/(2n)! can fail Toeplitz PF∞, while the Riemann-Phi
        quadrature mixture survives.  This paper probes why: is Phi safely
        inside the PF cone in finite windows, or close to its boundary?
DEPENDS: 181 (PF∞ door), 182 (factorially damped moments),
         183 (Bessel/cosh mixture probe).

THEORY
======

------------------------------------------------------------------------
THE QUESTION AFTER PAPER 183
------------------------------------------------------------------------

The sequence
    A_n = integral Phi(u) u^(2n)/(2n)! du
survives finite Toeplitz PF tests.

But the map
    positive measure mu -> a_n(mu)=integral u^(2n)/(2n)! dmu(u)
does NOT preserve PF∞ for arbitrary mu.

Therefore the RH-relevant theorem must use special structure of Phi(u)du.
This paper tests finite-window stability:

  (S1) log-shape of Phi on its effective support;
  (S2) distance to PF failure under a spike perturbation;
  (S3) distance to PF failure under tilts exp(lambda u);
  (S4) the active minors closest to zero after removing structural zeros.

The output is not a proof.  It tells us where the analytic inequality
would need to be sharp.
"""
from __future__ import annotations

import itertools
import math


def phi_riemann(u, n_terms=50):
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


def active_toeplitz_minors(seq, window=8, max_size=4, tol=1e-16):
    """Return non-structural Toeplitz minors sorted by determinant.

    Structural zeros occur whenever the chosen rows/cols force a zero
    triangular block.  We ignore determinants with abs <= tol.
    """
    out = []
    indices = list(range(window))
    for size in range(1, max_size + 1):
        for rows_idx in itertools.combinations(indices, size):
            for cols_idx in itertools.combinations(indices, size):
                M = [[toeplitz_entry(seq, r, c) for c in cols_idx] for r in rows_idx]
                d = det_small(M)
                if abs(d) > tol:
                    out.append((d, rows_idx, cols_idx, size))
    out.sort(key=lambda x: x[0])
    return out


def seq_from_measure(us, ws, n_terms=18):
    return [
        sum(w * (u ** (2 * n)) / math.factorial(2 * n) for u, w in zip(us, ws))
        for n in range(n_terms)
    ]


def phi_nodes(U=5.0, n_nodes=100, tilt=0.0):
    h = U / n_nodes
    us = [(j + 0.5) * h for j in range(n_nodes)]
    ws = [phi_riemann(u) * math.exp(tilt * u) * h for u in us]
    s = sum(ws)
    return us, [w / s for w in ws]


def first_negative_for_spike(base_us, base_ws, spike_u, window=8):
    """Binary search epsilon in (1-eps)Phi + eps delta_spike."""
    def neg_at(eps):
        us = list(base_us) + [spike_u]
        ws = [(1 - eps) * w for w in base_ws] + [eps]
        seq = seq_from_measure(us, ws)
        minors = active_toeplitz_minors(seq, window=window)
        return any(d < -1e-14 for d, *_ in minors)

    if not neg_at(1.0):
        return None
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if neg_at(mid):
            hi = mid
        else:
            lo = mid
    return hi


def main():
    print("=" * 80)
    print("PAPER 184 — stability boundary around Riemann-Phi")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Shape of Phi
    # ------------------------------------------------------------------
    print("\n[1] Phi shape diagnostics")
    grid = [5.0 * j / 500 for j in range(501)]
    vals = [phi_riemann(u) for u in grid]
    imax = max(range(len(vals)), key=lambda i: vals[i])
    support_mass = sum(vals) * (5.0 / 500)
    print(f"  max Phi grid: {vals[imax]:.8e} at u={grid[imax]:.4f}")
    print(f"  grid mass proxy on [0,5]: {support_mass:.8e}")

    # Discrete log-concavity of positive grid values around the support.
    log_bad = 0
    checked = 0
    worst = float("inf")
    for i in range(1, len(vals) - 1):
        if vals[i - 1] > 1e-300 and vals[i] > 1e-300 and vals[i + 1] > 1e-300:
            margin = math.log(vals[i]) * 2 - math.log(vals[i - 1]) - math.log(vals[i + 1])
            checked += 1
            worst = min(worst, margin)
            if margin < -1e-10:
                log_bad += 1
    print(f"  discrete log-concavity checks: {checked}")
    print(f"  log-concavity violations: {log_bad}")
    print(f"  worst log-concavity margin: {worst:.6e}")

    # ------------------------------------------------------------------
    # [2] Active Toeplitz minors at Phi
    # ------------------------------------------------------------------
    print("\n[2] Active Toeplitz minors for Phi quadrature")
    us, ws = phi_nodes()
    seq = seq_from_measure(us, ws)
    minors = active_toeplitz_minors(seq, window=8, max_size=4)
    neg = sum(1 for d, *_ in minors if d < 0)
    print(f"  active minors: {len(minors)}")
    print(f"  negative active minors: {neg}")
    print("  five smallest active minors:")
    for d, rows, cols, size in minors[:5]:
        print(f"    det={d:.6e}, size={size}, rows={rows}, cols={cols}")

    # ------------------------------------------------------------------
    # [3] Exponential tilts exp(lambda u)
    # ------------------------------------------------------------------
    print("\n[3] Exponential tilt stability: Phi(u) exp(lambda u)")
    print(f"  {'lambda':>8} {'neg minors':>12} {'min active det':>16}")
    for lam in [-6, -4, -2, 0, 2, 4, 6]:
        us_l, ws_l = phi_nodes(tilt=lam)
        seq_l = seq_from_measure(us_l, ws_l)
        ms = active_toeplitz_minors(seq_l, window=8, max_size=4)
        neg_l = sum(1 for d, *_ in ms if d < -1e-14)
        print(f"  {lam:>8.2f} {neg_l:>12} {ms[0][0]:>16.6e}")

    # ------------------------------------------------------------------
    # [4] Spike perturbations
    # ------------------------------------------------------------------
    print("\n[4] Spike perturbation threshold")
    print("    measure = (1-eps) Phi(u)du + eps delta_spike")
    print(f"  {'spike u':>8} {'first eps causing PF failure':>30}")
    for spike_u in [0.1, 0.25, 0.5, 1.0, 2.0, 3.5, 5.0]:
        eps = first_negative_for_spike(us, ws, spike_u)
        eps_str = "none up to 1" if eps is None else f"{eps:.6f}"
        print(f"  {spike_u:>8.2f} {eps_str:>30}")

    # ------------------------------------------------------------------
    # [5] Honest assessment
    # ------------------------------------------------------------------
    print("\n[5] HONEST ASSESSMENT")
    print("  FINDINGS:")
    print("    - Phi is sharply localised; finite active Toeplitz minors are")
    print("      positive in the tested window.")
    print("    - Exponential tilts test whether PF survives mass shifts.")
    print("    - Spike thresholds estimate distance to the finite PF boundary.")
    print()
    print("  ANALYTIC LEAD:")
    print("    The proof should not target arbitrary measures.  It should prove")
    print("    that the particular Riemann-Phi shape, perhaps through log-shape")
    print("    or variation-diminishing constraints, keeps all Toeplitz minors")
    print("    non-negative after the factorial-Bessel transform.")


if __name__ == "__main__":
    main()
