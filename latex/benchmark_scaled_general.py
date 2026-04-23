"""
Does scaling transport to general polynomials?

For P(z) = z^p + R(z) and target root alpha_k, the substitution z = beta*w
transforms P into  tilde P(w) = w^p + R(beta*w)/A   with A = beta^p.
The target root becomes tilde alpha_k = alpha_k / beta.  If beta is
chosen close to |alpha_k|, then tilde alpha_k is close to 1 and the
iteration should be well-conditioned.

Empirically compare:
  (a) unscaled multi-start Pandrosion-Steffensen on P
  (b) scaled multi-start Pandrosion-Steffensen on tilde P, with beta = round(|alpha_k|)

Test polynomials chosen so that roots are NOT close to 1 in absolute value.
"""
from __future__ import annotations
import mpmath as mp


def S_p_fast(p, s):
    if mp.fabs(s - 1) < mp.mpf(10) ** (-30):
        return mp.mpf(p)
    return (s ** p - 1) / (s - 1)


def eval_R(R_coefs, s):
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
    acc = mp.mpc(0)
    for c in reversed(coefs):
        acc = acc * z + c
    return acc


def K_sigma_probe(p, R_coefs, alpha, R):
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


def multistart(p, R_coefs, P_coefs, alpha, z, tol, R_sigma, max_iters=500):
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


def scale_polynomial(p, R_coefs, beta):
    """Return tilde R such that tilde P(w) = w^p + tilde R(w) satisfies
    P(beta w) / A = tilde P(w), with A = beta^p.
    tilde R(w) = R(beta w) / A, so tilde R[i] = R[i] * beta^i / A."""
    A = mp.mpf(beta) ** p
    return [R_coefs[i] * mp.mpf(beta) ** i / A for i in range(len(R_coefs))]


# ---------------------------------------------------------------------------
# Test cases with roots |alpha| significantly > 1
# ---------------------------------------------------------------------------

def run():
    mp.mp.dps = 70
    tol = mp.mpf(10) ** (-50)
    z0 = mp.mpc("3", "3")

    # Test case 1: z^3 - 2z - 100, real root ~ 4.83
    # Test case 2: z^5 - 3 z - 1000, real root ~ 3.98
    # Test case 3: z^3 - 10, real root ~ 2.154 (monomial, sanity)
    # Test case 4: z^4 + 2 z^3 - 1000, real root ~ 5.03
    CASES = [
        ("z^3 - 2z - 100",   3, [-100, -2, 0], [-100, -2, 0, 1]),
        ("z^5 - 3z - 1000",  5, [-1000, -3, 0, 0, 0], [-1000, -3, 0, 0, 0, 1]),
        ("z^3 - 10 (monomial sanity)", 3, [-10, 0, 0], [-10, 0, 0, 1]),
        ("z^4 + 2 z^3 - 1000", 4, [-1000, 0, 0, 2], [-1000, 0, 0, 2, 1]),
    ]

    for name, p, R_coefs, P_coefs in CASES:
        print(f"\n{'=' * 75}")
        print(f"Polynomial: {name}")
        print(f"{'=' * 75}")

        # Find principal real root
        try:
            alpha = mp.findroot(
                lambda z: eval_poly([mp.mpc(c) for c in P_coefs], z),
                mp.mpc(abs(P_coefs[0]) ** (mp.mpf(1) / p))
            )
            alpha = mp.mpc(alpha)
        except Exception as e:
            print(f"findroot failed: {e}")
            continue
        print(f"Target root alpha     = {alpha}")

        # ---- (a) UNSCALED ----
        gap = max(abs(alpha) * 2 * mp.sin(mp.pi / p), mp.mpf("0.1"))
        R_unscaled = find_R_sigma(p, R_coefs, alpha, gap)
        _, iters_unscaled = multistart(p, R_coefs, P_coefs, alpha, z0, tol, R_unscaled)
        print(f"  (a) Unscaled:  R_sigma = {float(R_unscaled):>10.4e},  iters = {iters_unscaled}")

        # ---- (b) SCALED with beta = round(|alpha|) ----
        beta = max(1, int(mp.nint(abs(alpha))))
        if beta == 1:
            print(f"  (b) Scaled:    beta = 1, no scaling possible")
            continue
        tilde_R = scale_polynomial(p, R_coefs, beta)
        # tilde P(w) = w^p + tilde R(w), so tilde P coefs = tilde R + [0,...,0,1]
        tilde_P = tilde_R + [mp.mpf(1)]
        tilde_alpha = alpha / mp.mpf(beta)
        print(f"  (b) Scaled:    beta = {beta},  tilde_alpha = {tilde_alpha}")

        gap_tilde = max(abs(tilde_alpha) * 2 * mp.sin(mp.pi / p), mp.mpf("0.05"))
        R_scaled = find_R_sigma(p, tilde_R, tilde_alpha, gap_tilde)
        # z_0 in w-space is z_0 / beta:
        w0 = z0 / mp.mpf(beta)
        _, iters_scaled = multistart(p, tilde_R, tilde_P, tilde_alpha, w0, tol, R_scaled)
        speedup = iters_unscaled / iters_scaled if iters_scaled > 0 else float('inf')
        print(f"                R_sigma_scaled = {float(R_scaled):>10.4e},  iters = {iters_scaled}  (speedup {speedup:.2f}x)")


if __name__ == "__main__":
    run()
