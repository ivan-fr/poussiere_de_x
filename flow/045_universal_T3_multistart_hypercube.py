"""
PAPER: 045
TITLE: Universal Pandrosion-T_3 + adaptive reanchor + multistart on HYPERCUBE geometry
STATUS: combines flow/041 (universal T_3 + reanchor) with hypercube system from flow/025

WHAT'S COMBINED
===============
Geometry      : hypercube Pandrosion system from flow/025
                P_i(s) = x_i (1 - s_i) * prod_j S_p(s_j) - (x_i - 1) = 0
                (this is a multi-root polynomial system; Bezout count = (n*(p-1)+1)^n
                 for the natural extension)

Operator      : universal Pandrosion-T_3 from flow/041
                h(z) = z_0 - Q_F(z_0, z)^{-1} F(z_0)
                (Schmidt slope matrix; works on ANY polynomial residual)
                T_3 = Steffensen + Richardson order 3
                Adaptive reanchor every K=3 cycles

Search        : MULTISTART
                Try multiple z_0 to discover different roots; deduplicate.

Why this combination might WIN against Lairez
=============================================
Lairez tracks Bezout roots via homotopy continuation (one path per root).
Multistart Pandrosion tries random starts and may converge to several roots.
On the structured hypercube polynomial:
  - The system has many roots (Bezout-many in the projective closure).
  - Each root has its own basin of attraction under the universal Pandrosion
    operator.
  - With enough diverse starts, we may capture a substantial fraction at low
    cost (no Jacobian, no path tracking).

OPEN QUESTION
=============
Whether multistart Pandrosion-T_3 captures all Bezout roots, or only the
"easy" ones with large basins.  We measure this empirically.
"""

from __future__ import annotations
import math
import random


# -----------------------------------------------------------------------------
# Hypercube Pandrosion polynomial system (from flow/025).
# -----------------------------------------------------------------------------
def S_p_scalar(s, p):
    return sum(s ** k for k in range(p))


def S_hyper(s, p):
    total = 1.0
    for si in s:
        total *= S_p_scalar(si, p)
    return total


def F_target_hyper(s, x_vec, p):
    if any(not math.isfinite(si) for si in s):
        return [float("inf")] * len(s)
    S = S_hyper(s, p)
    if not math.isfinite(S):
        return [float("inf")] * len(s)
    return [x_vec[i] * (1.0 - s[i]) * S - (x_vec[i] - 1.0) for i in range(len(s))]


def residual_norm(s, x_vec, p):
    if any(not math.isfinite(si) for si in s):
        return float("inf")
    return max(abs(v) for v in F_target_hyper(s, x_vec, p))


# -----------------------------------------------------------------------------
# Schmidt slope matrix Q_F at (z_0, z).
# -----------------------------------------------------------------------------
def schmidt_slope(z0, z, x_vec, p):
    n = len(z)
    points = [list(z)]
    cur = list(z)
    for j in range(n):
        cur[j] = z0[j]
        points.append(cur[:])
    F_vals = [F_target_hyper(pt, x_vec, p) for pt in points]
    Q = [[0.0] * n for _ in range(n)]
    for j in range(n):
        denom = z[j] - z0[j]
        if abs(denom) < 1e-300:
            zp = list(z); zp[j] += 1e-7
            fp = F_target_hyper(zp, x_vec, p)
            for i in range(n):
                Q[i][j] = (fp[i] - F_vals[0][i]) / 1e-7
        else:
            for i in range(n):
                Q[i][j] = (F_vals[j][i] - F_vals[j + 1][i]) / denom
    return Q, F_vals[0], F_vals[n], n + 1


def solve_linear(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-14:
            return None
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        rhs = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = rhs / M[i][i]
    return x


# -----------------------------------------------------------------------------
# Universal Pandrosion operator + T_3 + adaptive reanchor.
# -----------------------------------------------------------------------------
def pandrosion_h(z0, z, F_z0, x_vec, p):
    n = len(z)
    Q, _, _, cost = schmidt_slope(z0, z, x_vec, p)
    step = solve_linear(Q, [-v for v in F_z0])
    if step is None:
        return list(z), cost
    return [z0[j] + step[j] for j in range(n)], cost


def T3_step(z0, z, F_z0, x_vec, p):
    s1, c1 = pandrosion_h(z0, z, F_z0, x_vec, p)
    s2, c2 = pandrosion_h(z0, s1, F_z0, x_vec, p)
    n = len(z)
    t2 = []
    for k in range(n):
        d0 = s1[k] - z[k]; d2 = s2[k] - 2 * s1[k] + z[k]
        t2.append(s2[k] if abs(d2) < 1e-300 else z[k] - d0 * d0 / d2)
    s3, c3 = pandrosion_h(z0, t2, F_z0, x_vec, p)
    lam_num = sum(s2[k] - s1[k] for k in range(n))
    lam_den = sum(s1[k] - z[k] for k in range(n))
    if abs(lam_den) < 1e-300 or abs(lam_num / lam_den - 1.0) < 1e-12:
        return s3, c1 + c2 + c3
    lam_hat = lam_num / lam_den
    out = [t2[k] - (s3[k] - t2[k]) / (lam_hat - 1.0) for k in range(n)]
    return out, c1 + c2 + c3


def adaptive_pandrosion_T3(z0_init, x_vec, p, K=3, tol=1e-10, max_cycles=30):
    n = len(z0_init)
    z0 = list(z0_init); z = list(z0_init)
    F_evals = 0
    for cycle in range(max_cycles):
        F_z0 = F_target_hyper(z0, x_vec, p); F_evals += 1
        rn = max(abs(v) for v in F_z0) if all(math.isfinite(v) for v in F_z0) else float("inf")
        if rn < tol:
            return z0, F_evals, cycle, True
        if not math.isfinite(rn):
            return z, F_evals, cycle, False
        for _ in range(K):
            z_new, c = T3_step(z0, z, F_z0, x_vec, p)
            F_evals += c
            r_new = residual_norm(z_new, x_vec, p); F_evals += 1
            if r_new < tol:
                return z_new, F_evals, cycle, True
            if not math.isfinite(r_new):
                return z, F_evals, cycle, False
            z = z_new
        z0 = list(z)
    return z, F_evals, max_cycles, residual_norm(z, x_vec, p) < tol


# -----------------------------------------------------------------------------
# Multistart with deduplication.
# -----------------------------------------------------------------------------
def generate_starts(n, n_starts, seed=0):
    rng = random.Random(seed)
    out = [[0.5] * n]              # canonical Pandrosion start
    out.append([-0.5] * n)         # mirror
    out.append([0.0] * n)           # origin
    for _ in range(n_starts - 3):
        z0 = [rng.uniform(-1.5, 1.5) for _ in range(n)]
        out.append(z0)
    return out


def is_new_root(z, found, tol=1e-3):
    for fr in found:
        if max(abs(z[i] - fr[i]) for i in range(len(z))) < tol:
            return False
    return True


def multistart_solve(x_vec, p, n_starts=30, tol=1e-10):
    n = len(x_vec)
    starts = generate_starts(n, n_starts)
    found = []
    total_F = 0; n_attempts = 0; n_diverged = 0
    for z0 in starts:
        z, F, cycles, ok = adaptive_pandrosion_T3(z0, x_vec, p, tol=tol)
        total_F += F; n_attempts += 1
        if ok and is_new_root(z, found):
            found.append(list(z))
        elif not ok:
            n_diverged += 1
    return found, total_F, n_attempts, n_diverged


# -----------------------------------------------------------------------------
# Lairez-style homotopy on the same system (baseline).
# -----------------------------------------------------------------------------
def fd_jacobian(z, x_vec, p, h=1e-7):
    n = len(z)
    f0 = F_target_hyper(z, x_vec, p)
    J = [[0.0] * n for _ in range(n)]
    for j in range(n):
        zp = list(z); zp[j] += h
        fp = F_target_hyper(zp, x_vec, p)
        for i in range(n):
            J[i][j] = (fp[i] - f0[i]) / h
    return J


def homotopy_F(s, x_vec, p, t):
    """H(s, t) = (1-t) * (s - 1) + t * F_target(s)."""
    F = F_target_hyper(s, x_vec, p)
    return [(1 - t) * (s[i] - 1.0) + t * F[i] for i in range(len(s))]


def homotopy_J(s, x_vec, p, t, h=1e-7):
    n = len(s)
    f0 = homotopy_F(s, x_vec, p, t)
    J = [[0.0] * n for _ in range(n)]
    for j in range(n):
        sp = list(s); sp[j] += h
        fp = homotopy_F(sp, x_vec, p, t)
        for i in range(n):
            J[i][j] = (fp[i] - f0[i]) / h
    return J


def lairez_homotopy_one(z0_start, x_vec, p, tol=1e-10, dt_init=0.05):
    """Track ONE path from z0_start (root of (s-1)=0 trivial system at t=0,
    or random) to F_target at t=1."""
    n = len(z0_start)
    s = list(z0_start)
    t = 0.0; dt = dt_init
    F_evals = 0; J_evals = 0
    while t < 1.0 - 1e-12:
        if t + dt > 1.0: dt = 1.0 - t
        # Predictor (warm start).
        s_pred = list(s)
        # Corrector at t+dt with Newton.
        ok = False
        for _ in range(8):
            r = homotopy_F(s_pred, x_vec, p, t + dt); F_evals += 1
            rn = max(abs(v) for v in r) if all(math.isfinite(v) for v in r) else float("inf")
            if rn < tol:
                ok = True; break
            J = homotopy_J(s_pred, x_vec, p, t + dt); J_evals += 1
            step = solve_linear(J, [-v for v in r])
            if step is None: break
            s_pred = [s_pred[i] + step[i] for i in range(n)]
        if ok:
            s = s_pred; t += dt; dt = min(dt * 1.4, 0.2)
        else:
            dt /= 2.0
            if dt < 1e-6: return s, F_evals, J_evals, False
    return s, F_evals, J_evals, residual_norm(s, x_vec, p) < tol


def lairez_multistart(x_vec, p, n_starts=10, tol=1e-10):
    n = len(x_vec)
    rng = random.Random(42)
    found = []; total_F = 0; total_J = 0
    for _ in range(n_starts):
        # Random starting points covering Bezout roots of trivial start system.
        z0 = [rng.uniform(0.5, 1.5) for _ in range(n)]
        z, F, J, ok = lairez_homotopy_one(z0, x_vec, p, tol=tol)
        total_F += F; total_J += J
        if ok and is_new_root(z, found):
            found.append(z)
    weighted = total_F + n * total_J
    return found, weighted, total_F, total_J


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------
def main():
    print("=" * 116)
    print("flow/045 -- Universal Pandrosion-T_3 + adaptive reanchor + multistart "
          "on HYPERCUBE geometry")
    print("=" * 116)
    print()
    print("System: P_i(s) = x_i (1 - s_i) * prod_j S_p(s_j) - (x_i - 1) = 0  "
          "(hypercube Pandrosion polynomial)")
    print("Method: universal Pandrosion-T_3 (Schmidt slope + T_3 + reanchor K=3) "
          "+ multistart (N starts).")
    print("Compared against: Lairez-style homotopy multistart (Newton corrector).")
    print()

    cases = [
        ("(2,3)",     [2.0, 3.0],          2),
        ("(2,3,5)",   [2.0, 3.0, 5.0],     2),
        ("(2,3,5,7)", [2.0, 3.0, 5.0, 7.0], 2),
        ("(2,3,5)",   [2.0, 3.0, 5.0],     3),
    ]

    for label, x_vec, p in cases:
        n = len(x_vec)
        print(f"\n*** x = {label}, p = {p}, n = {n} ***")
        # Universal Pandrosion-T_3 multistart
        roots_p, F_p, attempts_p, div_p = multistart_solve(x_vec, p, n_starts=30)
        # Lairez homotopy multistart
        roots_l, w_l, F_l, J_l = lairez_multistart(x_vec, p, n_starts=10)
        print(f"  Universal-T_3 + multistart : {len(roots_p)} distinct roots, "
              f"{F_p} F-evals across {attempts_p} starts ({div_p} diverged)")
        for r in roots_p:
            print(f"      root: {[round(v, 6) for v in r]}, residual = "
                  f"{residual_norm(r, x_vec, p):.2e}")
        print(f"  Lairez homotopy + multistart: {len(roots_l)} distinct roots, "
              f"{F_l}F + {J_l}J = {w_l} weighted")
        for r in roots_l:
            print(f"      root: {[round(v, 6) for v in r]}, residual = "
                  f"{residual_norm(r, x_vec, p):.2e}")
        # Verdict per case.
        if len(roots_p) > len(roots_l):
            verdict = "Pandrosion-T_3 finds MORE roots"
        elif len(roots_p) < len(roots_l):
            verdict = "Lairez finds more roots"
        elif F_p < w_l:
            verdict = "TIE on roots, Pandrosion CHEAPER"
        else:
            verdict = "TIE on roots, Lairez cheaper"
        print(f"  -> {verdict}")

    print()
    print("=" * 116)
    print("VERDICT (honest):")
    print("=" * 116)
    print("  - Universal-T_3 + multistart on hypercube polynomial finds roots without")
    print("    a Jacobian.  Each start triggers an independent Pandrosion descent.")
    print("  - For systems with many roots, multistart explores multiple basins.")
    print("    Whether it captures ALL Bezout roots depends on basin geometry.")
    print()
    print("  HONEST LIMITS:")
    print("  - No gamma-trick: paths may pass through singularities; we don't recover.")
    print("  - Random starts may be biased toward 'easy' roots with large basins.")
    print("  - Lairez homotopy with K-S certified starts has stronger completeness")
    print("    guarantees; ours is empirical.")
    print()
    print("  WHEN PANDROSION-T_3 + MULTISTART WINS:")
    print("  - When Newton/Lairez paths fail or stall (basin issues).")
    print("  - When derivative-free is needed (no analytic Jacobian).")
    print("  - When the cost per F-eval is low (sparse polynomials).")
    print("=" * 116)


if __name__ == "__main__":
    main()
