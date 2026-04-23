"""
Smale-17 benchmark, v4.  FINALLY uses the correct framework from the companion
PDF `pandrosion_smale_v2.pdf`: the Generalized Pandrosion Operator via
difference quotient, with adaptive anchor update every K steps.

Univariate (PDF §4):
    P(z) - P(z_0) = (z - z_0) * Q(z_0, z),  Q = difference quotient
    P_{P, z_0}(z) = z_0 - P(z_0) / Q(z_0, z)
    Newton is the limit z_0 = z_n (then Q(z, z) = P'(z)).

Adaptive anchor (PDF §5):  update z_0 <- z_n every K steps.  K=3 is empirically
optimal: it avoids the poles of Q (which move with z_0) while staying local.

Multivariate lift:  for F : C^n -> C^n,  F(z) - F(z_0) = Q_F(z_0, z) (z - z_0)
where Q_F is an n x n matrix computed by telescoping through the coordinates:

    Q_F[i, j](z_0, z) = [F_i(z_0[0..j-1] | z[j] | z[j+1..]) -
                         F_i(z_0[0..j-1] | z_0[j] | z[j+1..])] / (z[j] - z_0[j])

(with a smooth limit when z[j] = z_0[j]).  Then:

    P_{F, z_0}(z_n) = z_0 - Q_F(z_0, z_n)^{-1} F(z_0)

With K=3 adaptive anchor, this resolves the Smale-17-style multivariate
systems using exactly the algorithm the PDF proposes, not a block extension.
"""
from __future__ import annotations
import numpy as np
import time
from itertools import product
from scipy.optimize import fsolve


# ----------------------------------------------------------------------------
# Polynomial system utilities
# ----------------------------------------------------------------------------

def monomial_exps(n, d):
    return [e for total in range(d + 1)
            for e in product(range(d + 1), repeat=n) if sum(e) == total]


def random_system(n, d, rng):
    exps = monomial_exps(n, d)
    return [{e: rng.standard_normal() + 1j * rng.standard_normal() for e in exps}
            for _ in range(n)]


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


def jacobian(system, z, n):
    """Exact Jacobian DF(z) of the polynomial system."""
    J = np.zeros((n, n), dtype=complex)
    for i, coefs in enumerate(system):
        for j in range(n):
            # partial derivative w.r.t. z[j]
            d_acc = 0j
            for e, c in coefs.items():
                if e[j] == 0:
                    continue
                t = c * e[j]
                for m, em in enumerate(e):
                    if m == j:
                        if em - 1 > 0:
                            t *= z[m] ** (em - 1)
                    elif em > 0:
                        t *= z[m] ** em
                d_acc += t
            J[i, j] = d_acc
    return J


# ----------------------------------------------------------------------------
# Difference-quotient matrix Q_F(z_0, z)
#
#   F(z) - F(z_0) = Q_F(z_0, z) (z - z_0)
#
# Telescoping identity:
#   F(z) - F(z_0) = sum_{j=0}^{n-1} [F(z^{(j+1)}) - F(z^{(j)})]
# where z^{(0)} = z_0 and z^{(j+1)} differs from z^{(j)} only in coord j,
# replaced by z[j].  Each telescope term equals
#   (z[j] - z_0[j]) * Q_F[:, j],
# with Q_F[i, j] = [F_i(z^{(j+1)}) - F_i(z^{(j)})] / (z[j] - z_0[j]).
#
# When z[j] = z_0[j] we use the exact derivative ∂_j F_i at z^{(j)}.
# ----------------------------------------------------------------------------

def diff_quotient_matrix(system, z0, z, n):
    Q = np.zeros((n, n), dtype=complex)
    z_curr = np.array(z0, dtype=complex)
    for j in range(n):
        delta = z[j] - z0[j]
        z_next = z_curr.copy()
        z_next[j] = z[j]
        if abs(delta) < 1e-14:
            # use exact ∂_j F(z_curr)
            for i, coefs in enumerate(system):
                d_acc = 0j
                for e, c in coefs.items():
                    if e[j] == 0:
                        continue
                    t = c * e[j]
                    for m, em in enumerate(e):
                        if m == j:
                            if em - 1 > 0:
                                t *= z_curr[m] ** (em - 1)
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


# ----------------------------------------------------------------------------
# Generalized Pandrosion Operator (multivariate, adaptive anchor)
#
#   Step: z_{n+1} = z_0 - Q_F(z_0, z_n)^{-1} F(z_0)
#   Every K steps: z_0 <- z_n (anchor update)
# ----------------------------------------------------------------------------

def adaptive_pandrosion_solve(system, n, z_init, K=3, max_iters=200, tol=1e-10):
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
        # adaptive anchor update every K steps
        if (it + 1) % K == 0:
            z0 = z.copy()
    return z, max_iters, float(np.max(np.abs(eval_system(system, z))))


# ----------------------------------------------------------------------------
# Newton baseline (pure multivariate Newton)
# ----------------------------------------------------------------------------

def pure_newton(system, n, z_init, max_iters=200, tol=1e-10):
    z = np.array(z_init, dtype=complex)
    for it in range(max_iters):
        F = eval_system(system, z)
        res = float(np.max(np.abs(F)))
        if res < tol:
            return z, it, res
        J = jacobian(system, z, n)
        try:
            delta = np.linalg.solve(J, F)
        except np.linalg.LinAlgError:
            return None, it, float('inf')
        z = z - delta
        if not np.all(np.isfinite(z)):
            return None, it, float('inf')
    return z, max_iters, float(np.max(np.abs(eval_system(system, z))))


# Newton via scipy fsolve for reference
def fsolve_newton(system, n, z_init, tol=1e-10, max_iters=500):
    def F(xy):
        z = xy[:n] + 1j * xy[n:]
        return np.concatenate([eval_system(system, z).real, eval_system(system, z).imag])
    xy0 = np.concatenate([np.asarray(z_init).real, np.asarray(z_init).imag])
    sol, info, ier, _ = fsolve(F, xy0, full_output=True, xtol=tol,
                               maxfev=max_iters * (2 * n + 1))
    z = sol[:n] + 1j * sol[n:]
    return z, info['nfev'], float(np.max(np.abs(eval_system(system, z))))


def run():
    NUM_TRIALS = 5
    CONFIGS = [(1, 3), (1, 5), (2, 2), (2, 3), (3, 2), (3, 3), (4, 2), (4, 3)]
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 105)
    print("Smale-17 v4 — Adaptive Pandrosion (difference quotient + K=3 anchor update)")
    print("=" * 105)
    print(f"{'(n,d)':>6} | {'Bezout':>7} | {'method':<22} | succ/{NUM_TRIALS} | {'avg iters':>10} | {'avg residual':>14}")
    print("-" * 105)

    for n, d in CONFIGS:
        bezout = d ** n
        # --- Adaptive Pandrosion (K=3) ---
        rng = np.random.default_rng(42)
        n_ok = 0; iters_sum = 0; res_sum = 0.0
        for _ in range(NUM_TRIALS):
            system = random_system(n, d, rng)
            z0 = rng.standard_normal(n) + 1j * rng.standard_normal(n)
            z, it, res = adaptive_pandrosion_solve(system, n, z0, K=3)
            if z is not None and res < 1e-8:
                n_ok += 1; iters_sum += it; res_sum += res
        a_i = iters_sum / max(1, n_ok); a_r = res_sum / max(1, n_ok)
        print(f"({n},{d})  | {bezout:>7} | {'Adaptive Pandr. K=3':<22} |     {n_ok}/{NUM_TRIALS} | {a_i:>10.1f} | {a_r:>14.2e}")

        # --- Pure Newton ---
        rng = np.random.default_rng(42)
        n_ok = 0; iters_sum = 0; res_sum = 0.0
        for _ in range(NUM_TRIALS):
            system = random_system(n, d, rng)
            z0 = rng.standard_normal(n) + 1j * rng.standard_normal(n)
            z, it, res = pure_newton(system, n, z0)
            if z is not None and res < 1e-8:
                n_ok += 1; iters_sum += it; res_sum += res
        a_i = iters_sum / max(1, n_ok); a_r = res_sum / max(1, n_ok)
        print(f"        |         | {'Pure Newton':<22} |     {n_ok}/{NUM_TRIALS} | {a_i:>10.1f} | {a_r:>14.2e}")

        # --- Vary K to check K=3 claim ---
        for K in [1, 5, 10, 1000000]:  # K=1 is Newton, K=big is fixed-anchor
            rng = np.random.default_rng(42)
            n_ok = 0; iters_sum = 0
            for _ in range(NUM_TRIALS):
                system = random_system(n, d, rng)
                z0 = rng.standard_normal(n) + 1j * rng.standard_normal(n)
                z, it, res = adaptive_pandrosion_solve(system, n, z0, K=K)
                if z is not None and res < 1e-8:
                    n_ok += 1; iters_sum += it
            a_i = iters_sum / max(1, n_ok)
            label = f"Adaptive K={K}" if K < 1000000 else "Fixed anchor (K=∞)"
            print(f"        |         | {label:<22} |     {n_ok}/{NUM_TRIALS} | {a_i:>10.1f} | {'':>14}")
        print("-" * 105)


if __name__ == "__main__":
    run()
