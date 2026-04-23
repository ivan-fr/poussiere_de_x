"""
Generalised Pandrosion-Steffensen for degree-p polynomials beyond monomials.

For P(z) = z^p + R(z) with deg R < p, the Pandrosion-style iteration
    h_P(s) := 1 - (R(s) + 1) / S_p(s),       S_p(s) = 1 + s + ... + s^{p-1}
has fixed points exactly at the roots of P:
    h_P(s) = s
  <=> (1 - s) S_p(s) = R(s) + 1
  <=> 1 - s^p = R(s) + 1
  <=> s^p = -R(s)
  <=> P(s) = 0.
This reduces to the monomial Pandrosion h_{p,x}(s) = 1 - (x-1)/(xS_p(s))
when R(s) = -1/x (constant).

Multi-start Pandrosion for general polynomials:
  sigma_P := Aitken(h_P) has super-attracting fixed points at each root
  (provided h_P'(alpha) != 1 and h_P is twice differentiable, which hold
  for simple roots alpha with S_p(alpha) != 0).
  Seed: z_0 = alpha_k + eps * (z - alpha_k) with eps = R_sigma / (2|z - alpha_k|)
  places any z in the quadratic-contraction disk of the chosen target.

The scheme is derivative-free: it uses R (the "low-order part" of P) and
S_p (a fixed universal polynomial), never P' directly.
"""
from __future__ import annotations
import mpmath as mp


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------

def S_p(p, s):
    """Geometric sum 1 + s + s^2 + ... + s^{p-1}."""
    acc = mp.mpc(0); sk = mp.mpc(1)
    for _ in range(p):
        acc = acc + sk; sk = sk * s
    return acc


def h_P(p, R_coefs, s):
    """Pandrosion iteration for P(z) = z^p + R(z), R given by its
    coefficient list R_coefs = [r_0, r_1, ..., r_{p-1}]."""
    R_s = sum(c * s ** i for i, c in enumerate(R_coefs))
    Sp = S_p(p, s)
    if abs(Sp) < mp.mpf(10) ** (-40):
        return mp.inf
    return 1 - (R_s + 1) / Sp


def sigma_P(p, R_coefs, s):
    """Steffensen on h_P: T(s) = s - (h(s) - s)^2 / (h(h(s)) - 2h(s) + s)."""
    h1 = h_P(p, R_coefs, s)
    h2 = h_P(p, R_coefs, h1)
    d = h2 - 2 * h1 + s
    if abs(d) < mp.mpf(10) ** (-40):
        return h1
    return s - (h1 - s) ** 2 / d


def eval_poly(coefs, z):
    """Evaluate polynomial a_0 + a_1 z + ... + a_d z^d at z."""
    acc = mp.mpc(0); zk = mp.mpc(1)
    for a in coefs:
        acc = acc + a * zk; zk = zk * z
    return acc


def plain_pandrosion(p, R_coefs, z0, tol, P_coefs, max_iters=300):
    """Plain Pandrosion (no Steffensen, no multi-start)."""
    s = z0
    for n in range(max_iters):
        if abs(eval_poly(P_coefs, s)) <= tol:
            return s, n
        s = h_P(p, R_coefs, s)
    return s, max_iters


def steff_pandrosion(p, R_coefs, z0, tol, P_coefs, max_iters=300):
    """Plain Steffensen-Pandrosion (no multi-start)."""
    s = z0
    for n in range(max_iters):
        if abs(eval_poly(P_coefs, s)) <= tol:
            return s, n
        s = sigma_P(p, R_coefs, s)
    return s, max_iters


def multistart_pandrosion(p, R_coefs, alpha, z, tol, R_sigma, P_coefs,
                          max_iters=300):
    """Multi-start Steffensen-Pandrosion targeting alpha."""
    d = abs(z - alpha)
    if d < mp.mpf(10) ** (-40):
        return alpha, 0
    eps = R_sigma / (2 * d) if 2 * d > R_sigma else mp.mpc(1)
    s = alpha + eps * (z - alpha)
    for n in range(max_iters):
        if abs(eval_poly(P_coefs, s)) <= tol:
            return s, n + 1
        s = sigma_P(p, R_coefs, s)
    return s, max_iters


# ---------------------------------------------------------------------------
# Test cases.  Format: (name, p, R_coefs, full_P_coefs)
# P(z) = z^p + sum_i R_coefs[i] z^i
# ---------------------------------------------------------------------------

TEST_CASES = [
    # Monomial sanity check: z^3 - 2 (should recover paper's p=3, x=2 case).
    ("z^3 - 2          (monomial, x=2)",
     3, [-2, 0, 0], [-2, 0, 0, 1]),
    # Trinomial / Bring-style:
    ("z^5 - z - 1",
     5, [-1, -1, 0, 0, 0], [-1, -1, 0, 0, 0, 1]),
    # Bombelli's cubic:
    ("z^3 - 2z - 5",
     3, [-5, -2, 0], [-5, -2, 0, 1]),
    # Quartic with small coefficient spread:
    ("z^4 + z + 1",
     4, [1, 1, 0, 0], [1, 1, 0, 0, 1]),
    # Sextic:
    ("z^6 + 3 z^2 - 1",
     6, [-1, 0, 3, 0, 0, 0], [-1, 0, 3, 0, 0, 0, 1]),
]


def all_roots(P_coefs):
    """Return all p = deg(P) roots of P via mpmath.polyroots."""
    return mp.polyroots(list(reversed(P_coefs)), maxsteps=200, extraprec=100)


def K_sigma_est(p, R_coefs, alpha, probe_R):
    """Empirical estimate of K_sigma by measuring |sigma(w) - alpha| /
    |w - alpha|^2 on |w - alpha| = probe_R."""
    angles = [k * mp.mpf(2) * mp.pi / 16 for k in range(16)]
    max_ratio = mp.mpf(0)
    for th in angles:
        w = alpha + probe_R * mp.exp(1j * th)
        sw = sigma_P(p, R_coefs, w)
        if abs(w - alpha) < mp.mpf(10) ** (-40):
            continue
        max_ratio = max(max_ratio, abs(sw - alpha) / abs(w - alpha) ** 2)
    return max_ratio


def find_R_sigma(p, R_coefs, alpha, roots, initial_upper=None):
    """Bisection: find largest R <= gap/2 such that sigma is 1/2-contracting
    on B(alpha, R), i.e. K_sigma(R) * R <= 1/2.  We take K_sigma as the
    empirical maximum of |sigma(w)-alpha|/|w-alpha|^2 on the circle |w-alpha|=R."""
    others = [r for r in roots if abs(r - alpha) > mp.mpf(10) ** (-20)]
    if others:
        gap = min(abs(r - alpha) for r in others) / 2
    else:
        gap = mp.mpf(10)
    if initial_upper is None:
        R_hi = gap
    else:
        R_hi = min(gap, initial_upper)
    R_lo = mp.mpf(0)
    for _ in range(30):
        R_mid = (R_lo + R_hi) / 2
        if R_mid < mp.mpf(10) ** (-15):
            break
        K = K_sigma_est(p, R_coefs, alpha, R_mid)
        if K * R_mid <= mp.mpf('0.5'):
            R_lo = R_mid
        else:
            R_hi = R_mid
    return R_lo


def run():
    mp.mp.dps = 60
    tol = mp.mpf(10) ** (-50)
    z0_external = mp.mpc("3", "3")

    for name, p, R_coefs, P_coefs in TEST_CASES:
        print(f"\n{'=' * 72}")
        print(f"Polynomial: {name}")
        print(f"  P(z) = z^{p} + R(z), deg R = {p - 1}, R_coefs = {R_coefs}")
        print(f"{'=' * 72}")
        roots = all_roots(P_coefs)
        print(f"Roots ({len(roots)} total):")
        for k, r in enumerate(roots):
            R_sig = find_R_sigma(p, R_coefs, r, roots)
            print(f"  root_{k} = {str(r)[:60]:<60}")
            print(f"             R_sigma (empirical) = {float(R_sig):.4f}")

        print(f"\nMulti-start Pandrosion-Steffensen from z_0 = {z0_external}:")
        for k, r in enumerate(roots):
            R_sig = find_R_sigma(p, R_coefs, r, roots)
            if R_sig < mp.mpf(10) ** (-10):
                print(f"  Target root_{k}: R_sigma ~= 0, skipping")
                continue
            _, steps = multistart_pandrosion(
                p, R_coefs, r, z0_external, tol, R_sig, P_coefs
            )
            print(f"  Target root_{k} = {str(r)[:40]:<40} -> {steps} iterations")


if __name__ == "__main__":
    run()
