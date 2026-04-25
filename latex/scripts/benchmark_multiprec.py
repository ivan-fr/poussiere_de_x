"""
Multi-precision benchmark: Newton vs Halley vs Steffensen-Pandrosion-scaling
for computing x^{1/p}. Uses mpmath (MPFR-backed) at configurable precision.

Iterations are counted with an iterate-independent stopping criterion
  |u^p - x| / x <= tol
so all three methods stop on the same condition. Wall-clock time is
measured over N_REP = 20 repeats to reduce noise.
"""
from __future__ import annotations
import mpmath as mp
import time

PREC_DIGITS = (50, 200, 500)        # decimal digits tested
N_REP = 20                           # repetitions for timing
MAX_ITERS = 300                      # safety cap

REGIMES = [
    # (p, x_str, label)
    (2, "2",        "p=2,   x=2       (baseline)"),
    (3, "2",        "p=3,   x=2       (Pandrosion)"),
    (3, "11",       "p=3,   x=11"),
    (3, "1000003",  "p=3,   x=10^6+3  (large, non-power)"),
    (5, "99",       "p=5,   x=99"),
    (10, "1048579", "p=10,  x=2^20+3  (large, non-power)"),
    (20, "2",       "p=20,  x=2       (large p)"),
]


# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

def newton(x, p, tol):
    """Newton on u -> u^p - x, initial guess from mp.root via Taylor seed."""
    u = mp.mpf(x) ** (mp.mpf(1) / p)  # use mpmath's root as seed? unfair...
    # Fairer: start from floor(x^{1/p}) as integer seed.
    u = mp.mpf(int(mp.floor(mp.mpf(x) ** (mp.mpf(1) / p))))
    if u < 1:
        u = mp.mpf(1)
    steps = 0
    while steps < MAX_ITERS:
        up_p = u ** p
        if mp.fabs(up_p - x) <= tol * x:
            return u, steps
        u = ((p - 1) * u + x / u ** (p - 1)) / p
        steps += 1
    return u, steps


def halley(x, p, tol):
    """Halley on u^p - x = 0; cubic order."""
    u = mp.mpf(int(mp.floor(mp.mpf(x) ** (mp.mpf(1) / p))))
    if u < 1:
        u = mp.mpf(1)
    steps = 0
    while steps < MAX_ITERS:
        up = u ** p
        if mp.fabs(up - x) <= tol * x:
            return u, steps
        u = u * ((p - 1) * up + (p + 1) * x) / ((p + 1) * up + (p - 1) * x)
        steps += 1
    return u, steps


def h_pandro(s, xp, p):
    Sp = mp.mpf(0)
    sk = mp.mpf(1)
    for _ in range(p):
        Sp += sk
        sk *= s
    return 1 - (xp - 1) / (xp * Sp)


def steffensen_pandrosion_scaled(x, p, tol):
    """Scaled Steffensen-Pandrosion: x = A * xp with A = floor(x^{1/p})^p,
    so xp in [1, (1 + 1/floor(x^{1/p}))^p)."""
    a = int(mp.floor(mp.mpf(x) ** (mp.mpf(1) / p)))
    if a < 1:
        a = 1
    A = mp.mpf(a) ** p
    xp = mp.mpf(x) / A
    if xp < 1:
        xp = mp.mpf(x); A = mp.mpf(1); a = 1
    # Initial: s0 = h(1) = 1 - (xp-1)/(xp*p)
    s = 1 - (xp - 1) / (xp * p)
    steps = 0
    while steps < MAX_ITERS:
        v_scaled = a * (xp * s ** (p - 1))  # readout
        u = v_scaled
        if mp.fabs(u ** p - x) <= tol * x:
            return u, steps
        # Steffensen on h
        h1 = h_pandro(s, xp, p)
        h2 = h_pandro(h1, xp, p)
        denom = h2 - 2 * h1 + s
        if mp.fabs(denom) < tol ** 2:
            s = h1
        else:
            s = s - (h1 - s) ** 2 / denom
        steps += 1
    u = a * (xp * s ** (p - 1))
    return u, steps


METHODS = [
    ("Newton",      newton,                  "quadratic, requires f'"),
    ("Halley",      halley,                  "cubic,     requires f', f''"),
    ("S--Pand+scal",steffensen_pandrosion_scaled,"quadratic, derivative-free"),
]


# ---------------------------------------------------------------------------
# Benchmark driver
# ---------------------------------------------------------------------------

def time_method(method, x, p, tol):
    best = float("inf")
    last_steps = None
    for _ in range(N_REP):
        t0 = time.perf_counter()
        _, steps = method(x, p, tol)
        dt = time.perf_counter() - t0
        if dt < best:
            best = dt
        last_steps = steps
    return best, last_steps


def run_benchmark():
    for digits in PREC_DIGITS:
        mp.mp.dps = digits + 15
        tol = mp.mpf(10) ** (-digits)
        print(f"\n{'=' * 78}")
        print(f"Precision: {digits} decimal digits (mp.dps = {mp.mp.dps}, tol = 10^-{digits})")
        print(f"{'=' * 78}")
        header = (f"{'Regime':<25} | " +
                  " | ".join(f"{name:>14}" for name, _, _ in METHODS))
        print(header)
        print("-" * len(header))
        for p, x_str, label in REGIMES:
            x = mp.mpf(x_str)
            row = f"{label:<25} | "
            for name, method, _ in METHODS:
                t, steps = time_method(method, x, p, tol)
                row += f"{steps:>3}it/{t*1000:>6.2f}ms | "
            print(row.rstrip(" |"))


if __name__ == "__main__":
    print(f"mpmath backend: {mp.libmp.BACKEND}")
    print(f"N repetitions for timing: {N_REP}")
    run_benchmark()
