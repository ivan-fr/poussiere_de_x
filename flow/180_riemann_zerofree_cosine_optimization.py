"""
PAPER: 180 (NEW — effective zero-free region search)
TITLE: Riemann zero-free regions via Pandrosion-zeta cosine-polynomial
       optimisation: a certified non-negative search harness
STATUS: RH remains OPEN.  Effective zero-free regions are classical;
        this paper builds a reproducible search harness for the
        de la Vallée Poussin / Mossinghoff-Trudgian-Yang door.  The
        reported score is a proxy, not a new record constant.
DEPENDS: 145 (Riemann zero-free effective), 148 (Pandrosion-Bombieri
         sharp R), 155 (Bombieri kernel SDP), 160 (MTY diagnostic).

THEORY
======

------------------------------------------------------------------------
THE DOOR
------------------------------------------------------------------------

For sigma > 1,
    F_zeta(s) := -zeta'(s) / zeta(s)
              = sum_n Lambda(n) n^{-s}.

If
    P(theta) = a_0 + 2 sum_{k=1}^m a_k cos(k theta) >= 0,
then
    sum_k a_k Re F_zeta(sigma + i k t) >= 0.

This is the de la Vallée Poussin positivity engine.  After inserting the
explicit formula for zeta'/zeta, one obtains zero-free regions
    sigma > 1 - 1 / (R log |t|)
where R is an explicit functional of the chosen non-negative cosine
polynomial and several analytic error constants.

The record-level MT/MTY/Yang constants require a full analytic pipeline.
This paper does NOT reproduce those constants.  It isolates the part the
Pandrosion framework can improve directly: finding good non-negative
cosine polynomials with large first harmonic relative to size/entropy.

------------------------------------------------------------------------
CERTIFIED NON-NEGATIVE PARAMETRISATION
------------------------------------------------------------------------

Use Fejer-Riesz:
    P(theta) = |q_0 + q_1 e^{i theta} + ... + q_m e^{i m theta}|^2.

For real q_j this gives coefficients
    a_k = sum_{j=0}^{m-k} q_j q_{j+k}.

Hence P(theta) >= 0 automatically, without grid assumptions.

Search score:
    score(P) = a_1 / max_theta P(theta)
and proxy
    R_proxy(P) = max_theta P(theta) / a_1.

Lower R_proxy is better for the positivity mechanism, but it is only the
shape part of the true zero-free-region constant.
"""
from __future__ import annotations

import math
from itertools import product


def coeffs_from_q(q):
    """Fejer-Riesz coefficients: P = |sum q_j e^{ij theta}|^2.

    Stored as [a_0, a_1, ...] for
        P(theta) = a_0 + 2 sum_{k>=1} a_k cos(k theta).
    """
    m = len(q) - 1
    coeffs = []
    for k in range(m + 1):
        coeffs.append(sum(q[j] * q[j + k] for j in range(m + 1 - k)))
    return coeffs


def trig_poly_value(coeffs, theta):
    val = coeffs[0]
    for k in range(1, len(coeffs)):
        val += 2 * coeffs[k] * math.cos(k * theta)
    return val


def minmax_on_grid(coeffs, n_grid=4096):
    pmin = float("inf")
    pmax = -float("inf")
    argmin = 0.0
    argmax = 0.0
    for j in range(n_grid):
        theta = 2 * math.pi * j / n_grid
        v = trig_poly_value(coeffs, theta)
        if v < pmin:
            pmin = v
            argmin = theta
        if v > pmax:
            pmax = v
            argmax = theta
    return pmin, pmax, argmin, argmax


def r_proxy(coeffs):
    if len(coeffs) < 2 or coeffs[1] <= 0:
        return float("inf")
    _, pmax, _, _ = minmax_on_grid(coeffs)
    return pmax / coeffs[1]


def normalise_q(q):
    norm = math.sqrt(sum(x * x for x in q))
    if norm == 0:
        return None
    # Force q_0 positive to remove the sign symmetry.
    qn = [x / norm for x in q]
    if qn[0] < 0:
        qn = [-x for x in qn]
    return qn


def hill_climb(seed, step=0.12, rounds=8):
    """Small deterministic coordinate search on q-space."""
    q = list(seed)
    best_coeffs = coeffs_from_q(q)
    best_R = r_proxy(best_coeffs)

    for r in range(rounds):
        improved = True
        local_step = step / (1.7 ** r)
        while improved:
            improved = False
            for i in range(len(q)):
                for sign in [-1.0, 1.0]:
                    trial = list(q)
                    trial[i] += sign * local_step
                    trial = normalise_q(trial)
                    if trial is None:
                        continue
                    coeffs = coeffs_from_q(trial)
                    Rp = r_proxy(coeffs)
                    if Rp < best_R:
                        q = trial
                        best_coeffs = coeffs
                        best_R = Rp
                        improved = True
    return q, best_coeffs, best_R


def main():
    print("=" * 80)
    print("PAPER 180 — zero-free regions via cosine-polynomial optimisation")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Classical baseline polynomials
    # ------------------------------------------------------------------
    print("\n[1] Classical positive-polynomial baselines")
    # de la Vallée Poussin: 3 + 4 cos(theta) + cos(2 theta)
    # in our convention [a0, a1, a2] = [3, 2, 0.5].
    baselines = [
        ("dlVP 3+4c+c2", [3.0, 2.0, 0.5]),
        ("Fejer |1+z|^2", coeffs_from_q([1.0 / math.sqrt(2), 1.0 / math.sqrt(2)])),
        ("Fejer |1+z+z^2|^2", coeffs_from_q([1.0 / math.sqrt(3)] * 3)),
        ("binomial deg3", coeffs_from_q(normalise_q([1.0, 3.0, 3.0, 1.0]))),
    ]
    print(f"  {'name':>24} {'deg':>4} {'min P':>12} {'max P':>12}"
          f" {'a1':>10} {'R_proxy':>12}")
    for name, coeffs in baselines:
        pmin, pmax, _, _ = minmax_on_grid(coeffs)
        print(f"  {name:>24} {len(coeffs)-1:>4} {pmin:>12.4e} {pmax:>12.4e}"
              f" {coeffs[1]:>10.4e} {r_proxy(coeffs):>12.6f}")

    # ------------------------------------------------------------------
    # [2] Coarse Fejer-Riesz search
    # ------------------------------------------------------------------
    print("\n[2] Coarse Fejer-Riesz search over small integer q-vectors")
    print("    P(theta) = |sum q_j e^{ij theta}|^2 is non-negative by construction.")
    best = None
    m = 4
    values = [-2, -1, 0, 1, 2]
    checked = 0
    for q_raw in product(values, repeat=m + 1):
        if q_raw[0] <= 0:
            continue
        q = normalise_q(list(q_raw))
        if q is None:
            continue
        coeffs = coeffs_from_q(q)
        if len(coeffs) < 2 or coeffs[1] <= 1e-12:
            continue
        Rp = r_proxy(coeffs)
        checked += 1
        if best is None or Rp < best[0]:
            best = (Rp, q, coeffs)

    if best is None:
        print("  no admissible candidate found")
        return
    best_R, best_q, best_coeffs = best
    print(f"  checked candidates: {checked}")
    print(f"  best coarse R_proxy: {best_R:.6f}")
    print("  best coarse q:", " ".join(f"{x:+.4f}" for x in best_q))
    print("  best coarse a:", " ".join(f"{x:+.6f}" for x in best_coeffs))

    # ------------------------------------------------------------------
    # [3] Deterministic local improvement
    # ------------------------------------------------------------------
    print("\n[3] Local coordinate improvement")
    q_opt, coeffs_opt, R_opt = hill_climb(best_q)
    pmin, pmax, argmin, argmax = minmax_on_grid(coeffs_opt, n_grid=8192)
    print(f"  improved R_proxy: {R_opt:.6f}")
    print(f"  min P on grid:    {pmin:.6e} at theta={argmin:.4f}")
    print(f"  max P on grid:    {pmax:.6e} at theta={argmax:.4f}")
    print("  q_opt:", " ".join(f"{x:+.8f}" for x in q_opt))
    print("  a_opt:", " ".join(f"{x:+.8f}" for x in coeffs_opt))
    print("  positivity status: CERTIFIED by Fejer-Riesz square form")

    # ------------------------------------------------------------------
    # [4] Pandrosion-zeta positivity check
    # ------------------------------------------------------------------
    print("\n[4] Pandrosion-zeta positivity check for sigma > 1")
    print("    S = sum_k a_k Re[-zeta'/zeta(sigma + i k t)] should be >= 0")
    try:
        from mpmath import mp, mpc, zeta
        mp.dps = 40
        print(f"  {'sigma':>8} {'t':>8} {'S':>18} {'>=0?':>8}")
        for sigma in [1.02, 1.05, 1.10]:
            for T in [14.0, 100.0, 1000.0]:
                S = mp.mpf("0")
                for k, ak in enumerate(coeffs_opt):
                    if k == 0:
                        s = mpc(sigma, 0)
                    else:
                        s = mpc(sigma, k * T)
                    F = -zeta(s, derivative=1) / zeta(s)
                    S += ak * mp.re(F)
                print(f"  {sigma:>8.3f} {T:>8.1f} {float(S):>18.8e}"
                      f" {str(S >= -mp.mpf('1e-25')):>8}")
    except ImportError:
        print("  [mpmath unavailable; skipped]")

    # ------------------------------------------------------------------
    # [5] Honest assessment
    # ------------------------------------------------------------------
    print("\n[5] HONEST ASSESSMENT")
    print("  WHAT THIS DOES:")
    print("    - Gives a certified non-negative cosine-polynomial search engine.")
    print("    - Produces candidate coefficients for the zero-free-region pipeline.")
    print("    - Verifies the Pandrosion-zeta positivity inequality numerically.")
    print()
    print("  WHAT THIS DOES NOT DO:")
    print("    - It does not claim a new Mossinghoff-Trudgian-Yang/Yang constant.")
    print("    - R_proxy is only the shape component of the full explicit")
    print("      zero-free-region constant R.")
    print()
    print("  NEXT TECHNICAL STEP:")
    print("    - Replace R_proxy by the exact MTY/Yang functional.")
    print("    - Add SDP/interval arithmetic for the Fejer-Riesz q-polynomial.")
    print("    - Compare against the degree-22 MTY polynomial coefficient-by-coefficient.")


if __name__ == "__main__":
    main()
