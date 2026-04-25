"""
Multi-precision benchmark for the complex case z^p = X:
  Newton, Halley, multi-start Steffensen-Pandrosion, Grand Master multi-start.

Uses mpmath (MPFR-backed) at 50 and 200 decimal digits. All methods
count a "step" as one call to their update; the seed-setup of multi-
start and Grand Master counts as step 1 since it does nontrivial work
(computes eps_seed, chooses nearest anchor).
"""
from __future__ import annotations
import mpmath as mp
import time

N_REP = 10
MAX_ITERS = 400

# Conservative contraction radius. The empirical R_sigma is usually
# larger (e.g. ~2 for (p=3, x=2)), but we use a uniform-safe value here
# so the seed always lands inside the proven quadratic-contraction disk.
R_SIGMA = mp.mpf("0.5")


# ---------------------------------------------------------------------------
# Iteration primitives
# ---------------------------------------------------------------------------

def h_c(z, x, p):
    Sp = 0
    zk = mp.mpc(1)
    for _ in range(p):
        Sp = Sp + zk
        zk = zk * z
    return 1 - (x - 1) / (x * Sp)


def sigma_c(z, x, p):
    h1 = h_c(z, x, p)
    h2 = h_c(h1, x, p)
    denom = h2 - 2 * h1 + z
    if abs(denom) == 0:
        return h1
    return z - (h1 - z) ** 2 / denom


# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

def newton_c(X, p, z0, tol):
    z = z0
    for n in range(MAX_ITERS):
        if abs(z ** p - X) <= tol * abs(X):
            return z, n
        z = ((p - 1) * z + X / z ** (p - 1)) / p
    return z, MAX_ITERS


def halley_c(X, p, z0, tol):
    z = z0
    for n in range(MAX_ITERS):
        zp = z ** p
        if abs(zp - X) <= tol * abs(X):
            return z, n
        z = z * ((p - 1) * zp + (p + 1) * X) / ((p + 1) * zp + (p - 1) * X)
    return z, MAX_ITERS


def multi_start_pandro_c(X, p, z0, tol):
    """Multi-start seeded at principal anchor alpha = X^{1/p}."""
    alpha = mp.power(X, mp.mpc(1) / p)
    x_P = 1 / X
    dist = abs(z0 - alpha)
    if dist == 0:
        return alpha, 0
    if 2 * dist > R_SIGMA:
        eps_seed = R_SIGMA / (2 * dist)
    else:
        eps_seed = mp.mpc(1)
    z = alpha + eps_seed * (z0 - alpha)
    for n in range(MAX_ITERS):
        if abs(z ** p - X) <= tol * abs(X):
            return z, n + 1
        z = sigma_c(z, x_P, p)
    return z, MAX_ITERS


def grand_master_c(X, p, z0, tol):
    """Grand Master: scaling + Voronoi-nearest selection + multi-start.
    Iterates sigma_{p, 1/X'}, whose fixed points satisfy z^p = X'.
    Then X^{1/p} = beta * X'^{1/p}."""
    X_abs = abs(X)
    beta = int(mp.floor(X_abs ** (mp.mpf(1) / p)))
    if beta < 2:
        beta = 1
    A = mp.mpf(beta) ** p
    X_prime = X / A
    x_P = 1 / X_prime

    # Fixed points of h_{p, x_P} are z with z^p = 1/x_P = X_prime.
    gamma_principal = mp.power(X_prime, mp.mpc(1) / p)
    omega = mp.exp(2j * mp.pi / p)
    gammas = [gamma_principal * omega ** k for k in range(p)]

    # Initial guess for X_prime^{1/p}: w0 = z0/beta.
    if z0 == 0:
        w0 = mp.mpc(1)
    else:
        w0 = z0 / mp.mpf(beta)

    s_star = min(range(p), key=lambda k: abs(gammas[k] - w0))
    target = gammas[s_star]

    dist = abs(w0 - target)
    if dist == 0:
        return mp.mpf(beta) * target, 0
    if 2 * dist > R_SIGMA:
        eps_seed = R_SIGMA / (2 * dist)
    else:
        eps_seed = mp.mpc(1)
    w = target + eps_seed * (w0 - target)

    for n in range(MAX_ITERS):
        z_approx = mp.mpf(beta) * w
        if abs(z_approx ** p - X) <= tol * abs(X):
            return z_approx, n + 1
        w = sigma_c(w, x_P, p)
    return mp.mpf(beta) * w, MAX_ITERS


METHODS = [
    ("Newton",       newton_c),
    ("Halley",       halley_c),
    ("MS-Pandro",    multi_start_pandro_c),
    ("Grand Master", grand_master_c),
]


# ---------------------------------------------------------------------------
# Regimes
# ---------------------------------------------------------------------------

REGIMES = [
    # (p, X, z0, label)
    (3, mp.mpc("2", "3"),          mp.mpc("1", "0"),   "p=3, X=2+3i,       z0=1"),
    (3, mp.mpc("10", "0"),         mp.mpc("2", "1"),   "p=3, X=10,         z0=2+i"),
    (4, mp.mpc("-1", "1"),         mp.mpc("1", "1"),   "p=4, X=-1+i,       z0=1+i"),
    (5, mp.mpc("100", "7"),        mp.mpc("2", "1"),   "p=5, X=100+7i,     z0=2+i"),
    (3, mp.mpc("1000003", "0"),    mp.mpc("50", "0"),  "p=3, X=10^6+3,     z0=50 (far)"),
    (3, mp.mpc("1000003", "17"),   mp.mpc("50", "10"), "p=3, X=10^6+17i,   z0=50+10i"),
    (7, mp.mpc("10", "0"),         mp.mpc("1", "1"),   "p=7, X=10,         z0=1+i"),
]


def time_method(method, X, p, z0, tol):
    best = float("inf")
    last_steps = None
    for _ in range(N_REP):
        t0 = time.perf_counter()
        _, steps = method(X, p, z0, tol)
        dt = time.perf_counter() - t0
        if dt < best:
            best = dt
        last_steps = steps
    return best, last_steps


def run_benchmark():
    for digits in (50, 200):
        mp.mp.dps = digits + 15
        tol = mp.mpf(10) ** (-digits)
        print(f"\n{'=' * 98}")
        print(f"Precision: {digits} decimal digits")
        print(f"{'=' * 98}")
        header = f"{'Regime':<34} |"
        for name, _ in METHODS:
            header += f" {name:>15} |"
        print(header)
        print("-" * len(header))
        for p, X, z0, label in REGIMES:
            row = f"{label:<34} |"
            for name, method in METHODS:
                t, steps = time_method(method, X, p, z0, tol)
                cell = f"{steps:>3}it/{t * 1000:>7.2f}ms"
                row += f" {cell:>15} |"
            print(row)


if __name__ == "__main__":
    print(f"Complex benchmark: Newton / Halley / Multi-start Pandrosion / Grand Master")
    print(f"R_sigma (conservative uniform value): {R_SIGMA}")
    print(f"N_REP = {N_REP}, MAX_ITERS = {MAX_ITERS}")
    run_benchmark()
