"""
How does the iteration count of multi-start Pandrosion-Steffensen scale
with the polynomial degree p?

We test two families at fixed 50-digit precision:

  (A) Monomial:       P(z) = z^p - 2                 (well-conditioned)
  (B) Trinomial:      P(z) = z^p - z - 1             (Bring-style)

for p in {3, 5, 10, 20, 50, 100, 200}.

Each run: target the principal real root alpha (computed analytically
for the monomial, numerically for the trinomial), seed from z0 = 3+3i,
iterate sigma_P until |P(s)| <= tol * max(1, |s|^p).
"""
from __future__ import annotations
import mpmath as mp
import time


# ---------------------------------------------------------------------------
# Primitives (same as benchmark_general_poly.py, specialised for speed)
# ---------------------------------------------------------------------------

def S_p_fast(p, s):
    """Evaluate S_p(s) = 1 + s + ... + s^(p-1) via geometric formula
    for speed at large p: S_p(s) = (s^p - 1)/(s - 1) when s != 1."""
    if mp.fabs(s - 1) < mp.mpf(10) ** (-30):
        return mp.mpf(p)
    return (s ** p - 1) / (s - 1)


def eval_R(R_coefs, s):
    """Horner evaluation of R(s) = R_coefs[0] + R_coefs[1]*s + ..."""
    acc = mp.mpc(0)
    for c in reversed(R_coefs):
        acc = acc * s + c
    return acc


def h_P(p, R_coefs, s):
    Sp = S_p_fast(p, s)
    if abs(Sp) < mp.mpf(10) ** (-30):
        return mp.inf
    return 1 - (eval_R(R_coefs, s) + 1) / Sp


def sigma_P(p, R_coefs, s):
    h1 = h_P(p, R_coefs, s)
    h2 = h_P(p, R_coefs, h1)
    d = h2 - 2 * h1 + s
    if abs(d) < mp.mpf(10) ** (-30):
        return h1
    return s - (h1 - s) ** 2 / d


def eval_poly(coefs, z):
    """Horner for full polynomial coefficients [a_0, ..., a_d]."""
    acc = mp.mpc(0)
    for c in reversed(coefs):
        acc = acc * z + c
    return acc


def multistart_pandrosion(p, R_coefs, P_coefs, alpha, z, tol, R_sigma,
                          max_iters=500):
    d = abs(z - alpha)
    if d < mp.mpf(10) ** (-40):
        return alpha, 0
    eps = R_sigma / (2 * d) if 2 * d > R_sigma else mp.mpc(1)
    s = alpha + eps * (z - alpha)
    for n in range(max_iters):
        val = eval_poly(P_coefs, s)
        scale = max(1, abs(s) ** p)
        if abs(val) <= tol * scale:
            return s, n + 1
        s = sigma_P(p, R_coefs, s)
    return s, max_iters


# ---------------------------------------------------------------------------
# Empirical R_sigma via bisection on a coarse grid
# ---------------------------------------------------------------------------

def K_sigma_probe(p, R_coefs, alpha, R):
    """Empirical upper bound on |sigma(w) - alpha| / |w - alpha|^2 on the
    circle |w - alpha| = R (8 sample points)."""
    max_ratio = mp.mpf(0)
    for k in range(8):
        theta = mp.mpf(2) * mp.pi * k / 8
        w = alpha + R * mp.exp(1j * theta)
        sw = sigma_P(p, R_coefs, w)
        if abs(w - alpha) < mp.mpf(10) ** (-30):
            continue
        ratio = abs(sw - alpha) / abs(w - alpha) ** 2
        if ratio > max_ratio:
            max_ratio = ratio
    return max_ratio


def find_R_sigma(p, R_coefs, alpha, gap_upper):
    """Find largest R <= gap_upper/2 such that K_sigma(R) * R <= 1/2."""
    R_lo, R_hi = mp.mpf(0), gap_upper / 2
    for _ in range(20):
        R_mid = (R_lo + R_hi) / 2
        if R_mid < mp.mpf(10) ** (-10):
            break
        K = K_sigma_probe(p, R_coefs, alpha, R_mid)
        if K * R_mid <= mp.mpf('0.5'):
            R_lo = R_mid
        else:
            R_hi = R_mid
    return R_lo


# ---------------------------------------------------------------------------
# Main benchmark
# ---------------------------------------------------------------------------

def run():
    mp.mp.dps = 70
    tol = mp.mpf(10) ** (-50)
    z0 = mp.mpc("3", "3")

    p_values = [3, 5, 10, 20, 50, 100, 200]

    print("=" * 80)
    print("Family A: monomial P(z) = z^p - 2    (target alpha = 2^(1/p), real)")
    print("=" * 80)
    print(f"{'p':>5} | {'alpha':>12} | {'R_sigma (est.)':>14} | {'iters':>5} | {'time (s)':>8}")
    print("-" * 80)
    for p in p_values:
        alpha = mp.mpf(2) ** (mp.mpf(1) / p)
        # Monomial R(s) = -2 (constant): h_P = 1 - (-2+1)/S_p = 1 + 1/S_p
        # Wait: for P(z) = z^p - 2, we have R(z) = -2.
        R_coefs = [mp.mpf(-2)] + [mp.mpf(0)] * (p - 1)
        P_coefs = [mp.mpf(-2)] + [mp.mpf(0)] * (p - 1) + [mp.mpf(1)]
        # Estimate R_sigma: other roots are 2^(1/p) * omega^k for k=1..p-1.
        # Gap to nearest: |2^(1/p)(omega - 1)| = 2^(1/p) * 2*sin(pi/p).
        gap = alpha * 2 * mp.sin(mp.pi / p)
        t0 = time.perf_counter()
        R_sig = find_R_sigma(p, R_coefs, alpha, gap)
        _, iters = multistart_pandrosion(p, R_coefs, P_coefs, alpha, z0, tol, R_sig)
        dt = time.perf_counter() - t0
        print(f"{p:>5} | {float(alpha):>12.6f} | {float(R_sig):>14.4e} | {iters:>5} | {dt:>8.2f}")

    print()
    print("=" * 80)
    print("Family B: trinomial P(z) = z^p - z - 1   (Bring-style, one real root)")
    print("=" * 80)
    print(f"{'p':>5} | {'alpha':>12} | {'R_sigma (est.)':>14} | {'iters':>5} | {'time (s)':>8}")
    print("-" * 80)
    for p in p_values:
        R_coefs = [mp.mpf(-1), mp.mpf(-1)] + [mp.mpf(0)] * (p - 2)
        P_coefs = [mp.mpf(-1), mp.mpf(-1)] + [mp.mpf(0)] * (p - 2) + [mp.mpf(1)]
        # Real root of z^p - z - 1 ~ 2^(1/p) for large p (since alpha^p ~ 2).
        seed_guess = mp.mpf(2) ** (mp.mpf(1) / p)
        try:
            alpha = mp.findroot(lambda z: z ** p - z - 1, seed_guess)
            alpha = mp.mpc(alpha)
        except Exception as e:
            print(f"{p:>5} | findroot failed: {e}")
            continue
        # Estimate gap via perturbation: for large p, other roots are near |z| = 1,
        # on the unit circle. Gap from alpha (slightly > 1) to unit circle is small.
        # Use conservative estimate: gap = min(|alpha - 1|, 2*pi/p).
        gap = min(abs(alpha - 1), 2 * mp.pi / p)
        if gap < mp.mpf("1e-5"):
            gap = mp.mpf("0.01")
        t0 = time.perf_counter()
        R_sig = find_R_sigma(p, R_coefs, alpha, gap)
        _, iters = multistart_pandrosion(p, R_coefs, P_coefs, alpha, z0, tol, R_sig)
        dt = time.perf_counter() - t0
        print(f"{p:>5} | {float(alpha.real):>12.6f} | {float(R_sig):>14.4e} | {iters:>5} | {dt:>8.2f}")


if __name__ == "__main__":
    run()
