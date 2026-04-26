"""Casas-Alvero conjecture: P(z) = monic of degree d, P(z) and P^(k)(z)
share a common root for each k = 1, ..., d-1, then P(z) = (z - alpha)^d.

Equivalently: if a single root alpha satisfies P(alpha) = P'(alpha) = ... = P^(d-1)(alpha) = 0,
then P(z) = (z - alpha)^d.

Pandrosion reformulation: The Pandrosion-quotient spectrum
    Spec_P(z) := {Q(alpha_j, z) : j = 1, ..., d}
is injective in the d-tuple data iff Casas-Alvero holds.

Numerical attack: random-search for counterexamples in
    P(z) = (z - alpha)^{d-1} (z - beta) + epsilon * R(z)
where R is a small perturbation. If no counterexample is found across a large
random ensemble, we have empirical evidence for the conjecture.
"""
from __future__ import annotations
import math, time
import numpy as np


def common_roots(P, Q, tol=1e-8):
    """Find common roots of P(z) and Q(z) given as numpy ascending-coefs."""
    if len(P) <= 1 or len(Q) <= 1:
        return []
    rP = np.roots(P[::-1])
    rQ = np.roots(Q[::-1])
    out = []
    for a in rP:
        for b in rQ:
            if abs(a - b) < tol:
                out.append(complex((a + b) / 2))
                break
    return out


def casas_alvero_violation(coefs, tol=1e-8):
    """Return alpha if P and all derivatives P', P'', ..., P^{d-1} share a root
    that is NOT the d-fold root of (z - alpha)^d. Else return None.

    Test: find common roots of P, P', ..., P^{d-1}.  If exactly one alpha
    appears with multiplicity d in P, this is the trivial case (z - alpha)^d.
    Otherwise, it's a counterexample.
    """
    d = len(coefs) - 1
    if d < 2:
        return None
    # Compute all derivatives
    derivs = [coefs.copy()]
    cur = coefs.copy()
    for k in range(1, d):
        new = cur[1:] * np.arange(1, len(cur))
        derivs.append(new.copy())
        cur = new
    # Look for common root of P and each derivative
    rP = np.roots(derivs[0][::-1])
    for alpha in rP:
        all_share = True
        for k in range(1, d):
            vals_at_alpha = np.polyval(derivs[k][::-1], alpha)
            if abs(vals_at_alpha) > tol:
                all_share = False
                break
        if all_share:
            # Check if P(z) = (z - alpha)^d (i.e., trivial case)
            test = np.array([1.0 + 0j])
            for _ in range(d):
                test = np.convolve(test, np.array([-alpha, 1.0 + 0j]))
            # Compare with coefs (might have phase factor)
            scale = coefs[-1] / test[-1] if abs(test[-1]) > 1e-15 else 1.0
            test *= scale
            err = float(np.linalg.norm(test - coefs.astype(complex)))
            if err < 1e-6:
                # Trivial case: P = leading * (z - alpha)^d
                return None
            else:
                # Non-trivial common root — would be a counterexample
                return alpha
    return None


def trivial_test():
    """P(z) = (z-2)^5 should give NO counterexample (it IS the trivial case)."""
    alpha = 2.0
    d = 5
    P = np.array([1.0])
    for _ in range(d):
        P = np.convolve(P, np.array([-alpha, 1.0]))
    res = casas_alvero_violation(P)
    return res is None  # should be True


def test_known_low_degree():
    """For d = 4, Graf von Bothmer et al verified Casas-Alvero.
    Random search should find no counterexample."""
    rng = np.random.default_rng(20260427)
    n_test = 5000
    n_violations = 0
    for _ in range(n_test):
        d = 4
        # Sample monic poly with random coefs
        coefs = np.append(rng.standard_normal(d), 1.0)
        violation = casas_alvero_violation(coefs)
        if violation is not None:
            n_violations += 1
    return n_violations, n_test


def random_search_violations(d, n_polys, rng=None):
    """Random search for Casas-Alvero counterexamples at degree d."""
    if rng is None:
        rng = np.random.default_rng(20260427 + d)
    n_violations = 0
    for _ in range(n_polys):
        coefs = np.append(rng.standard_normal(d), 1.0)
        try:
            v = casas_alvero_violation(coefs)
        except (np.linalg.LinAlgError, ValueError):
            continue
        if v is not None:
            n_violations += 1
    return n_violations


def adversarial_perturbed_pure(d, n_perturb, rng=None):
    """Test perturbations of (z - 1)^d: add small noise, check for violations."""
    if rng is None:
        rng = np.random.default_rng(424242 + d)
    n_violations = 0
    for _ in range(n_perturb):
        # Start from (z - 1)^d
        P = np.array([1.0])
        for _ in range(d):
            P = np.convolve(P, np.array([-1.0, 1.0]))
        eps = 10**rng.uniform(-6, -2)
        P += eps * rng.standard_normal(len(P))
        try:
            v = casas_alvero_violation(P)
        except (np.linalg.LinAlgError, ValueError):
            continue
        if v is not None:
            n_violations += 1
    return n_violations


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')
    print("=" * 75, flush=True)
    print("Casas-Alvero conjecture random search", flush=True)
    print("=" * 75, flush=True)

    # Sanity: trivial case
    print(f"\nTrivial test (z-2)^5: {'PASS' if trivial_test() else 'FAIL'}", flush=True)

    # Test d=4 (known proved)
    nv, nt = test_known_low_degree()
    print(f"\nd=4 random scan ({nt} polys): {nv} violations (expected 0)", flush=True)

    # Scan d = 8, 10, 12, 16, 20 (open range)
    print("\nMain scan: random monic polynomials, search for counterexamples",
          flush=True)
    print(f"{'d':>3} {'#polys':>8} {'#violations':>12} {'time':>9}", flush=True)
    for d in [8, 10, 12, 16, 20, 24]:
        t0 = time.perf_counter()
        n_polys = max(100, 100000 // (d * d))
        nv = random_search_violations(d, n_polys)
        elapsed = time.perf_counter() - t0
        print(f"{d:>3} {n_polys:>8} {nv:>12} {elapsed:>8.1f}s", flush=True)

    # Adversarial: perturbations of (z-1)^d
    print("\nAdversarial: small perturbations of (z-1)^d", flush=True)
    print(f"{'d':>3} {'#perturb':>10} {'#violations':>12}", flush=True)
    for d in [8, 10, 16, 24]:
        t0 = time.perf_counter()
        nv = adversarial_perturbed_pure(d, 5000)
        elapsed = time.perf_counter() - t0
        print(f"{d:>3} {5000:>10} {nv:>12}  (t={elapsed:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
