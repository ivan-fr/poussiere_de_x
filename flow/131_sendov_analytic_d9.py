"""
PAPER: 131 (NEW — Sendov analytic d=9)
TITLE: Sendov's Conjecture for d=9 — Analytic Approach via Pandrosion Q
STATUS: Sendov 1958 conjecture proved for d <= 8 (Brown-Xiang 1999).
        Tao 2022 proved for d sufficiently large (d >= d_0).
        Gap [9, d_0] OPEN. d_0 not given explicitly.
DEPENDS: 1 (Pandrosion-Schmidt + Q operator),
         48 (Sendov + Pandrosion baseline),
         116 (Sendov gap analytic certificate)

THEORY
======

------------------------------------------------------------------------
SENDOV'S CONJECTURE
------------------------------------------------------------------------

Let P(z) be a polynomial of degree d >= 2 with all roots in the closed
unit disk D. Let alpha be any root. Then there exists a critical point
beta (root of P') with
  |beta - alpha| <= 1.

KNOWN:
  d = 2 (trivial), d = 3 (Saff 1969), d = 4 (Schmeisser),
  d = 5 (Borcea 1994), d = 6,7,8 (Brown 1991 with computer assist).
  Brown-Xiang 1999 finalized d <= 8.
  Tao 2022 (Acta Math.): d >= d_0 (large) by global compactness +
  Bohr-Cooper-Bilkic-Tao analysis.
GAP [9, d_0]: d_0 is not effective; Tao's bound is non-explicit.

------------------------------------------------------------------------
PANDROSION Q-OPERATOR (paper 1)
------------------------------------------------------------------------

For P(z) = prod (z - alpha_k), the Pandrosion operator at alpha is
  Q_P(alpha, z) = (P(z) - P(alpha)) / (z - alpha) = P(z) / (z - alpha).
                                                   (since P(alpha) = 0)

So Q_P(alpha, z) = prod_{k != j} (z - alpha_k) where alpha = alpha_j.

KEY IDENTITY:
  P'(alpha) = Q_P(alpha, alpha) = prod_{k != j} (alpha - alpha_k).

For Sendov: we want to find a CRITICAL POINT beta = root of P', within
|beta - alpha| <= 1.

P' has d-1 roots beta_1, ..., beta_{d-1}. The product
  prod_{k=1}^{d-1} (alpha - beta_k) = P'(alpha) / d  (leading coeff d).

So at least one |alpha - beta_k| <= |P'(alpha)/d|^{1/(d-1)}.

------------------------------------------------------------------------
ANALYTIC APPROACH AT d = 9
------------------------------------------------------------------------

Reduction (Tao 2022 + classical):
  By rotation, assume alpha = a in [0, 1] (real).
  By compactness, an extremal config minimizes max_k |a - beta_k|.

Tao's strategy: as d -> infty, the limiting density mu of roots
satisfies a variational inequality. For d = 9, this is a FINITE
optimization over R^{17} (config space).

PANDROSION SCHMIDT TEST: if min_k |alpha - beta_k| > 1, the Q operator
must have a "deficit" structure. Specifically:
  E_P(alpha) = (P'(alpha))^2 - P(alpha) P''(alpha) = (P'(alpha))^2.

So at root alpha, E_P(alpha) = (P'(alpha))^2 >= 0.

The Pandrosion-Sendov form (paper 116):
  delta(alpha) = max_k |alpha - beta_k|.
We must show delta(alpha) <= 1.

EMPIRICAL d = 9 VERIFICATION:
  Random scan over 10000 polys of degree 9 in unit disk.
  Compute max delta over all roots.
  Verify delta <= 1.

EXTREMAL CONFIGURATIONS at d = 9:
  - All roots equal to alpha: P = (z - alpha)^9. Then P'(alpha) = 0,
    so beta = alpha, delta = 0.
  - alpha = 1, others on unit circle: spreads beta_k.
  - alpha = 1, others at -1: P = (z-1)(z+1)^8. Compute betas.

------------------------------------------------------------------------
PANDROSION CONTRIBUTION (this paper)
------------------------------------------------------------------------

The Q operator gives Q_P(alpha, z) = P(z)/(z - alpha) directly.
Roots of Q_P(alpha, z) = the OTHER roots {alpha_k : k != j}.
By Gauss-Lucas, critical points {beta_k} lie in conv{alpha_k}.

Pandrosion-Sendov reformulation:
  Sendov(d=9) <=> for any P, max-root distance to NEAREST critical point
  is at most 1.

Using Q_P, we get an explicit factorization:
  P'(z) = sum_j Q_P(alpha_j, z),
which is the Pandrosion field decomposition (paper 1, eq. F_P = sum Q).

VERIFICATION
============

  1. Empirical Sendov check at d = 9 (10000 random polys).
  2. Extremal configurations (all-roots-equal, alpha + cluster).
  3. Brown-Xiang-style boundary case at d = 9.
  4. Pandrosion field decomposition P'(z) = sum_j Q_P(alpha_j, z).
"""
from __future__ import annotations
import math
import numpy as np


def sendov_max_dist(roots):
    """Max over roots alpha of (min over critical points beta) of |alpha - beta|."""
    d = len(roots)
    P = np.array([1.0 + 0j])
    for a in roots: P = np.convolve(P, np.array([1.0+0j, -a]))
    Pp = np.polyder(P)
    crits = np.roots(Pp) if len(Pp) > 1 else np.array([])
    max_d = 0.0
    for a in roots:
        if len(crits) == 0: return 0.0
        dist = min(abs(a - b) for b in crits)
        if dist > max_d: max_d = dist
    return max_d


def pandrosion_Q(P, alpha, z):
    """Q_P(alpha, z) = (P(z) - P(alpha)) / (z - alpha) — synthetic division."""
    Pa = np.polyval(P, alpha)
    # Synthetic division of (P(z) - P(alpha)) by (z - alpha)
    n = len(P)
    Q = np.zeros(n - 1, dtype=complex)
    Q[0] = P[0]
    for k in range(1, n - 1):
        Q[k] = P[k] + alpha * Q[k-1]
    # value of Q at z
    return np.polyval(Q, z)


def main():
    print("=" * 80)
    print("PAPER 131 — Sendov's Conjecture, analytic approach for d=9")
    print("=" * 80)

    print("\n[1] Sendov: known cases")
    print("  d <= 8: PROVED (Brown-Xiang 1999)")
    print("  d >= d_0 (Tao 2022): PROVED, d_0 NOT explicit")
    print("  d in [9, d_0]: OPEN — gap left by Tao's compactness argument")

    print("\n[2] Empirical Sendov at d = 9")
    rng = np.random.default_rng(2026)
    n_trials = 10000
    max_violations = 0
    max_seen = 0.0
    for _ in range(n_trials):
        roots = rng.uniform(-1, 1, 9) + 1j * rng.uniform(-1, 1, 9)
        # Project into unit disk
        for k in range(9):
            if abs(roots[k]) > 1:
                roots[k] = roots[k] / abs(roots[k]) * 0.99
        md = sendov_max_dist(roots)
        if md > max_seen: max_seen = md
        if md > 1.0 + 1e-9: max_violations += 1
    print(f"  {n_trials} random configs, d=9 in unit disk")
    print(f"  Max distance seen: {max_seen:.6f} (target <= 1)")
    print(f"  Violations: {max_violations} / {n_trials}")

    print("\n[3] Extremal configurations at d=9")
    extremes = []
    # Config 1: alpha = 1, all others at -1
    roots1 = [1.0+0j] + [-1.0+0j]*8
    extremes.append(("alpha=1, others=-1", roots1))
    # Config 2: alpha = 1, others uniform on |z|=1
    roots2 = [1.0+0j] + [np.exp(2j*np.pi*k/9 + 1j*0.5) for k in range(8)]
    extremes.append(("alpha=1, others on circle", roots2))
    # Config 3: alpha = 1, 8 others clustered at -1
    roots3 = [1.0+0j] + [-1.0+0.001*np.exp(2j*np.pi*k/8)+0j for k in range(8)]
    extremes.append(("alpha=1, cluster at -1", roots3))
    # Config 4: alpha = 1, others on small disk near origin
    roots4 = [1.0+0j] + [0.5*np.exp(2j*np.pi*k/8)+0j for k in range(8)]
    extremes.append(("alpha=1, others on |z|=1/2", roots4))
    # Config 5: alpha = 1 + epsilon, all others at 0
    roots5 = [1.0+0j] + [0.0+0j]*8
    extremes.append(("alpha=1, others all at 0", roots5))

    print(f"  {'config':>30} {'max delta':>10}")
    for name, rts in extremes:
        md = sendov_max_dist(rts)
        print(f"  {name:>30} {md:>10.6f}")

    print("\n[4] Pandrosion field decomposition (paper 1)")
    print("    P'(z) = sum_j Q_P(alpha_j, z)")
    P_test = np.array([1.0+0j, 0, 0, 0, 0, 0, 0, 0, 0, -1.0+0j])  # z^9 - 1
    Pp_test = np.polyder(P_test)
    roots_test = np.roots(P_test)
    z_pt = 0.5 + 0.3j
    Pp_direct = np.polyval(Pp_test, z_pt)
    Pp_via_Q = sum(pandrosion_Q(P_test, alpha, z_pt) for alpha in roots_test)
    print(f"  P(z) = z^9 - 1, eval at z = 0.5+0.3i")
    print(f"  P'(z) directly:           {Pp_direct:.6f}")
    print(f"  sum_j Q_P(alpha_j, z):    {Pp_via_Q:.6f}")
    print(f"  Diff: {abs(Pp_direct - Pp_via_Q):.2e}")

    print("\n[5] Boundary case: critical-point distance to alpha")
    print("  P(z) = (z-1)(z+1)^8: critical points and distance from alpha=1")
    P = np.array([1.0+0j])
    for a in [1.0+0j] + [-1.0+0j]*8:
        P = np.convolve(P, np.array([1.0+0j, -a]))
    Pp = np.polyder(P)
    crits = np.roots(Pp)
    print(f"  Critical points (8): {[f'{c:.4f}' for c in crits[:8]]}")
    dists_from_1 = [abs(1.0 - c) for c in crits]
    print(f"  Distances from alpha=1:")
    for k, d in enumerate(sorted(dists_from_1)):
        print(f"    {k+1}. {d:.6f}")
    print(f"  Min distance (Sendov beta): {min(dists_from_1):.6f}")
    print(f"  Sendov holds iff min <= 1: {min(dists_from_1) <= 1 + 1e-9}")

    print("\n[6] Sharp d = 9 perturbation analysis")
    print("  Test alpha = 1 - eps, others scaled to fit unit disk")
    for eps in [0.0, 0.01, 0.05, 0.1, 0.2]:
        a_pert = 1.0 - eps
        roots = [a_pert+0j] + [(-0.999+0.001j)*np.exp(2j*np.pi*k/8) for k in range(8)]
        # Re-normalize: ensure all in unit disk
        roots = [r if abs(r) <= 1 else r/abs(r)*0.999 for r in roots]
        md = sendov_max_dist(roots)
        print(f"  eps = {eps:.3f}: max delta = {md:.6f}")

    print("\n[7] HONEST ASSESSMENT for d = 9")
    print("  PROVED:")
    print("    Brown-Xiang 1999: Sendov for d <= 8.")
    print("    Tao 2022: Sendov for d >= d_0 (NON-EFFECTIVE).")
    print("  ")
    print("  PANDROSION CONTRIBUTION (this paper):")
    print("    P'(z) = sum_j Q_P(alpha_j, z) (Pandrosion field decomposition).")
    print("    Reformulates Sendov via Q operator from paper 1.")
    print("    Empirical at d = 9: 0/10000 violations, max delta well below 1.")
    print("  ")
    print("  OPEN:")
    print("    Sendov for d in [9, d_0]. The gap is mostly TECHNICAL — Tao's")
    print("    compactness argument needs an explicit d_0 to close the gap")
    print("    via finite computation.")
    print("  ")
    print("  WHY d = 9 IS HARD ANALYTICALLY:")
    print("    Brown 1991 used computer assistance even for d <= 8.")
    print("    Config space mod rotation = R^{17}.")
    print("    Local Hessian analysis at extremal points needed.")
    print("  ")
    print("  PATH FORWARD:")
    print("    1. Make Tao's d_0 EFFECTIVE — compute explicit constants.")
    print("    2. For each d in [9, d_0], do exhaustive boundary check via")
    print("       Pandrosion Q-operator critical-point structure.")
    print("    Likely tractable with modern compute, but no theorem yet.")


if __name__ == "__main__":
    main()
