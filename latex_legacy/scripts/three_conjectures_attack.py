"""
Triple attack: Lehmer + Casas-Alvero + De Bruijn-Newman via Pandrosion tools.

For each conjecture:
  1. Reformulate via Pandrosion-Gram / Pandrosion-spectrum
  2. Search for structural inequality analogous to det G_norm >= det K (Atiyah)
  3. Empirical verification at extreme scale
"""
from __future__ import annotations
import math
import numpy as np
import itertools
from typing import List


# =====================================================================
# PART 1: LEHMER — Pandrosion-Gram reformulation
# =====================================================================

def mahler_measure(coeffs):
    """M(P) = |a_d| * prod max(1, |alpha_j|)."""
    roots = np.roots(coeffs)
    a_d = abs(coeffs[0])
    return float(a_d * np.prod(np.maximum(1.0, np.abs(roots))))


def pandrosion_quotients_on_circle(coeffs, n_pts=128):
    """Compute Q_j(z) = P(z)/(z - alpha_j) sampled on |z|=1.

    Returns matrix Q[j, k] = Q_j(e^{i theta_k}) and the Gram-on-circle.
    """
    roots = np.roots(coeffs)
    d = len(roots)
    thetas = 2 * np.pi * np.arange(n_pts) / n_pts
    z = np.exp(1j * thetas)
    Q = np.zeros((d, n_pts), dtype=complex)
    for j in range(d):
        # Q_j(z) = prod_{k != j} (z - alpha_k)
        prod = np.ones(n_pts, dtype=complex)
        for k in range(d):
            if k == j:
                continue
            prod *= (z - roots[k])
        Q[j] = prod
    # Gram matrix: G_{ij} = (1/n_pts) sum_k Q_i(z_k) bar(Q_j(z_k))
    G = (Q @ Q.conj().T) / n_pts
    return Q, G, roots


def lehmer_pandrosion_gram(coeffs):
    """Compute the Pandrosion-Gram determinant for Lehmer.

    By paper 77: det G^(P) = |Disc(P)|^2 = prod_{i<j} |alpha_i - alpha_j|^2.

    Test: relate log M(P) to log det G (Lehmer in Gram language).
    """
    Q, G, roots = pandrosion_quotients_on_circle(coeffs, n_pts=256)
    # Discriminant
    d = len(roots)
    disc = 1.0
    log_disc = 0.0
    for i in range(d):
        for j in range(i+1, d):
            log_disc += 2 * math.log(max(abs(roots[i] - roots[j]), 1e-300))
    sg, log_det_G = np.linalg.slogdet(G)
    return dict(
        log_M=float(np.log(mahler_measure(coeffs))),
        log_disc=float(log_disc),
        log_det_G=float(np.real(log_det_G)),
        roots_outside=[r for r in roots if abs(r) > 1.0 + 1e-12],
    )


def lehmer_smyth_polynomial():
    """Smyth's polynomial L_0(z) = z^10 + z^9 - z^7 - z^6 - z^5 - z^4 - z^3 + z + 1."""
    return np.array([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], dtype=float)


def lehmer_scan_height_h(d_max, h, n_samples, rng):
    """Scan random monic integer polys of degree <= d_max, height <= h.

    Return min M(P) > 1 found.
    """
    min_M = float('inf')
    n_above_one = 0
    smallest = None
    for _ in range(n_samples):
        d = rng.integers(2, d_max + 1)
        coefs = rng.integers(-h, h + 1, size=d + 1)
        coefs[0] = 1  # monic
        if coefs[-1] == 0:
            continue  # not full degree
        M = mahler_measure(coefs.astype(float))
        if M > 1.0001:
            n_above_one += 1
            if M < min_M:
                min_M = M
                smallest = coefs.copy()
    return min_M, n_above_one, smallest


def lehmer_attack():
    print("=" * 95, flush=True)
    print("LEHMER ATTACK", flush=True)
    print("=" * 95, flush=True)

    # Verify Smyth's polynomial
    L0 = lehmer_smyth_polynomial()
    M_L0 = mahler_measure(L0)
    print(f"\nSmyth's polynomial L_0: M(L_0) = {M_L0:.7f} (lit: 1.17628081)")
    info = lehmer_pandrosion_gram(L0)
    print(f"  log M = {info['log_M']:.6f}")
    print(f"  log disc(P)^2 = {info['log_disc']:.4f}")
    print(f"  log det G (Pandrosion-circle) = {info['log_det_G']:.4f}")
    # By paper 77: det G = |disc|^2 (modulo normalization)
    print(f"  ratio (log det G - log disc^2) = {info['log_det_G'] - info['log_disc']:.4e}")

    # Smyth-bound test: scan
    print("\nSCAN: smallest M(P) > 1 across random monic integer polys")
    print(f"{'d_max':>6} {'height h':>9} {'#samples':>10} {'min M(P)':>12} {'#above 1':>10}",
          flush=True)
    rng = np.random.default_rng(20260601)
    for d_max in [10, 15, 20, 25, 30, 40]:
        for h in [1, 2]:
            n = 5000 if d_max < 30 else 2000
            min_M, n_above, smallest = lehmer_scan_height_h(d_max, h, n, rng)
            print(f"{d_max:>6} {h:>9} {n:>10} {min_M:>12.6f} {n_above:>10}",
                  flush=True)

    # KEY ATTACK: Pandrosion-Gram structural relation
    print("\nKEY ATTACK: relate log M(P) to Pandrosion-Gram structure")
    print(f"{'d':>4} {'M(P)':>10} {'log M':>10} {'log det G':>11} "
          f"{'log disc^2':>12} {'cyclotomic?':>14}", flush=True)
    # Test on cyclotomic, Smyth, plastic, random
    test_polys = [
        ("z^2+1 (cyclotomic)", np.array([1, 0, 1])),
        ("z^3-z-1 (plastic)", np.array([1, 0, -1, -1])),
        ("Smyth L_0", L0),
        ("z^4+z^3+1", np.array([1, 1, 0, 0, 1])),
        ("z^6-z^5-z^3-z+1", np.array([1, -1, 0, -1, 0, -1, 1])),
    ]
    for name, p in test_polys:
        info = lehmer_pandrosion_gram(p.astype(float))
        M = mahler_measure(p.astype(float))
        is_cyclo = "yes" if abs(M - 1.0) < 1e-8 else "no"
        print(f"  {name}", flush=True)
        print(f"    M = {M:.6f}, log M = {info['log_M']:+.6f}, "
              f"log det G = {info['log_det_G']:+.4f}, "
              f"log disc^2 = {info['log_disc']:+.4f}",
              flush=True)


# =====================================================================
# PART 2: CASAS-ALVERO — Pandrosion-spectrum extension
# =====================================================================

def shares_root(P_coeffs, Q_coeffs, tol=1e-6):
    """Test if P and Q share a common root within tol."""
    roots_P = np.roots(P_coeffs)
    roots_Q = np.roots(Q_coeffs)
    for r_p in roots_P:
        for r_q in roots_Q:
            if abs(r_p - r_q) < tol:
                return True, r_p
    return False, None


def derivatives(coeffs):
    """Return list of P', P'', ..., P^(d-1)."""
    d = len(coeffs) - 1
    derivs = []
    P = coeffs.copy()
    for k in range(1, d):
        P = np.polyder(P)
        derivs.append(P)
    return derivs


def casas_alvero_check(coeffs, tol=1e-6):
    """Check Casas-Alvero condition: P shares a root with each P^(k), 1 <= k <= d-1.

    Returns: (satisfies_CA_condition, all_roots_equal)
    """
    d = len(coeffs) - 1
    if d < 2:
        return False, False
    derivs = derivatives(coeffs)
    shared_roots = []
    for k, dP in enumerate(derivs, start=1):
        ok, r = shares_root(coeffs, dP, tol)
        if not ok:
            return False, None
        shared_roots.append(r)
    # If condition holds, check if all shared roots are equal (meaning P = (z-alpha)^d)
    r0 = shared_roots[0]
    all_equal = all(abs(r - r0) < tol for r in shared_roots)
    return True, all_equal


def casas_alvero_search(d, n_trials, rng, perturbation=0.01):
    """Search for a counterexample: random P near (z-1)^d, check if it satisfies CA
    without being (z-1)^d itself."""
    counterexamples = 0
    near_pure_power = 0
    for _ in range(n_trials):
        # Start from (z-1)^d
        base = np.array([1])
        for _ in range(d):
            base = np.convolve(base, np.array([1, -1]))
        base = base.astype(float)
        # Perturb
        perturb = perturbation * rng.standard_normal(len(base))
        perturb[0] = 0  # keep monic
        P = base + perturb
        satisfies, all_equal = casas_alvero_check(P)
        if satisfies:
            near_pure_power += 1
            if not all_equal:
                counterexamples += 1
    return counterexamples, near_pure_power


def casas_alvero_pandrosion_spectrum(coeffs):
    """Compute the Pandrosion-quotient spectrum at all roots and check injectivity.

    Spec_P(z) := {Q(alpha_j, z)}_j where Q(alpha, z) = P(z) - P(alpha) / (z - alpha)
              = P(z) / (z - alpha)  [since alpha is a root]

    For Casas-Alvero: if P shares a root with each P^(k), then certain Spec values
    coincide.
    """
    roots = np.roots(coeffs)
    d = len(roots)
    if d < 2:
        return None

    # Compute Spec at a generic point
    z0 = 0.5 + 0.7j
    spec_at_z0 = []
    for j in range(d):
        # Q(alpha_j, z0) = P(z0) / (z0 - alpha_j)
        Pz = np.polyval(coeffs, z0)
        Q_val = Pz / (z0 - roots[j])
        spec_at_z0.append(Q_val)
    return np.array(spec_at_z0)


def casas_alvero_attack():
    print("\n" + "=" * 95, flush=True)
    print("CASAS-ALVERO ATTACK", flush=True)
    print("=" * 95, flush=True)

    # Verify pure power case
    print("\n[1] Verification: (z-1)^d satisfies Casas-Alvero")
    for d in [3, 5, 8]:
        base = np.array([1.0])
        for _ in range(d):
            base = np.convolve(base, np.array([1, -1]))
        ok, all_eq = casas_alvero_check(base)
        print(f"  d = {d}: satisfies CA condition? {ok}, all roots equal? {all_eq}")

    # Counterexample search at large degrees
    print("\n[2] Search for counterexamples (random polys near (z-1)^d):")
    print(f"{'d':>3} {'#trials':>9} {'#near pure power':>20} {'#counterexamples':>20} "
          f"{'status':>10}",
          flush=True)
    rng = np.random.default_rng(20260602)
    open_degrees = [8, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24]  # non-prime-power, > 7
    for d in open_degrees:
        ctr, near = casas_alvero_search(d, 1000, rng, perturbation=0.001)
        status = "PASSES" if ctr == 0 else f"COUNTEREXAMPLE FOUND ({ctr})"
        print(f"{d:>3} {1000:>9} {near:>20} {ctr:>20} {status:>10}",
              flush=True)

    # Pandrosion-spectrum on test polys
    print("\n[3] Pandrosion spectrum (sample at z_0 = 0.5+0.7i):")
    for d in [4, 5, 6]:
        base = np.array([1.0])
        for _ in range(d):
            base = np.convolve(base, np.array([1, -1]))
        spec = casas_alvero_pandrosion_spectrum(base)
        print(f"  (z-1)^{d}: spec = {[f'{abs(s):.4f}' for s in spec]}",
              flush=True)


# =====================================================================
# PART 3: DE BRUIJN-NEWMAN constant Lambda
# =====================================================================

def riemann_xi(s):
    """Riemann xi function: xi(s) = (1/2) s (s-1) pi^{-s/2} Gamma(s/2) zeta(s).

    For numerical purposes, use Riemann's symmetric form on the critical line:
    xi(1/2 + it) is real and even in t."""
    # Use mpmath if available
    try:
        from mpmath import mp, mpc, zeta, gamma, pi, mpf
        mp.dps = 30
        s_mp = mpc(complex(s))
        return complex(0.5 * s_mp * (s_mp - 1) * pi**(-s_mp/2) *
                      gamma(s_mp/2) * zeta(s_mp))
    except ImportError:
        return None


def Phi_riemann(t, t_max=4.0, n_terms=20):
    """Newman's Phi function: integrand of de Bruijn-Newman.

    Phi(u) = sum_{n=1}^infty (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})

    Used to define H_t(z) = (1/8) int_0^infty Phi(u) e^{i z u} du via heat-flow."""
    s = 0.0
    for n in range(1, n_terms + 1):
        e9u = math.exp(9 * t)
        e5u = math.exp(5 * t)
        e4u = math.exp(4 * t)
        s += (2 * math.pi**2 * n**4 * e9u - 3 * math.pi * n**2 * e5u) * \
             math.exp(-math.pi * n**2 * e4u)
    return s


def H_t_function(z, t):
    """Newman's H_t(z) = (1/8) integral_0^infty Phi(u) cos(z u) e^{t u^2} du.

    Approximate via numerical integration. RH ⟺ H_0 has only real zeros.
    Tao-Rodgers (2018) proved Lambda >= 0 (existence of complex zero of some H_t).
    Newman conjecture: Lambda <= 0 — implies RH.
    """
    # Numerical integration over u in [0, U]
    U = 4.0
    n_pts = 100
    du = U / n_pts
    integral = 0.0 + 0j
    for k in range(n_pts):
        u = (k + 0.5) * du
        phi = Phi_riemann(u)
        integral += phi * np.cos(z * u) * np.exp(t * u**2)
    return integral * du / 8.0


def H_t_real_rooted_test(t, T_max=20.0, n_pts=200, tol=1e-3):
    """Test if H_t(z) has only real zeros for z in [-T_max, T_max].

    Look for sign changes of H_t along the real axis."""
    zs = np.linspace(-T_max, T_max, n_pts)
    vals = [H_t_function(z, t) for z in zs]
    real_parts = [v.real for v in vals]
    # Count sign changes of real part
    n_zeros = 0
    for i in range(len(real_parts) - 1):
        if real_parts[i] * real_parts[i+1] < 0:
            n_zeros += 1
    # All zeros are real iff #(real zeros) >= total zeros (which is harder to count).
    # Simpler test: compute H_t at off-axis z = x + iy for small y.
    # If H_t is real-rooted, |H_t(x + iy)| > 0 for y != 0 (no off-axis zeros nearby).
    return n_zeros, vals


def de_bruijn_newman_attack():
    print("\n" + "=" * 95, flush=True)
    print("DE BRUIJN-NEWMAN LAMBDA ATTACK", flush=True)
    print("=" * 95, flush=True)
    print("Goal: explore the value of Lambda where H_t transitions from real-rooted")
    print("(i.e., RH-like) to having complex zeros. Tao-Rodgers (2018): Lambda >= 0.")
    print("Newman conj: Lambda <= 0. Combined: Lambda = 0 ⇔ RH at boundary.")

    # Test H_t for various t
    print(f"\n{'t':>10} {'#real-axis sign changes':>26} {'H_0(0)':>14} {'H_0(10)':>14}",
          flush=True)
    for t in [-0.1, -0.05, -0.01, 0.0, 0.01, 0.05, 0.1]:
        n_zeros, vals = H_t_real_rooted_test(t, T_max=20.0, n_pts=200)
        H_0 = H_t_function(0.0, t)
        H_10 = H_t_function(10.0, t)
        print(f"{t:>10.3f} {n_zeros:>26} {H_0.real:>14.4f} {H_10.real:>14.4f}",
              flush=True)

    print("\nObservation: at t = 0, H_0 corresponds to the Riemann xi function.")
    print("RH ⟺ H_0 has only real zeros. We see sign changes (real-axis zeros)")
    print("but cannot rule out off-axis zeros from real-axis data alone.")

    # Try real-stability test on shifted polynomials (Pólya approach)
    print("\nPólya-Schur test on polynomial truncations of cosh / Riemann xi:")
    # Approximate xi on critical line by polynomial in t^2
    from mpmath import mp, mpf, mpc, sin, cos, pi, exp, gamma, zeta
    mp.dps = 30
    print("  (mpmath required for Riemann zeros)")

    # Test: does Pandrosion-Pólya-Schur (paper 98) tools say anything about xi?
    print("\n  Riemann xi(1/2 + i t) at t = 14.13... (first zero):")
    try:
        from mpmath import mp, mpc, zeta, gamma, pi
        mp.dps = 30
        s = mpc('0.5', '14.134725')
        xi = 0.5 * s * (s - 1) * pi**(-s/2) * gamma(s/2) * zeta(s)
        print(f"    xi = {complex(xi)}")
    except ImportError:
        print("    mpmath not available")


# =====================================================================
# Main
# =====================================================================

def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    lehmer_attack()
    casas_alvero_attack()
    de_bruijn_newman_attack()


if __name__ == "__main__":
    main()
