"""
PAPER: 004 (canonical: 4pandrosion_smale.pdf)
TITLE: A Structural Pandrosion Reduction for Smale's MVC --
       The Inverse Quotient and the Fiber Vanishing Identity
STATUS: structural reduction (not a proof for d>=5)
DEPENDS: 003

THEORY
======

Setup: P monic of degree d with P(0) = 0 and P'(0) = 1.
Roots: 0 = alpha_0, alpha_1, ..., alpha_{d-1}.
Critical points: c_1, ..., c_{d-1}.

PANDROSION INVERSE (Theorem 2.1):
  Q(c, 0) = P(c)/c, hence S(c) = |Q(c, 0)|.

CO-FIBER (Definition 2.2):
The co-fiber of a critical point c is Fib(c) = {z_1, ..., z_{d-2}},
the roots of Q_2(c, z) = (P(z) - P(c))/(z - c)^2 (i.e., solutions of
P(z) = P(c) other than z = c with multiplicity).

CO-FIBER PRODUCT (Theorem 2.3):
  S(c) = |c| * prod_{j=1}^{d-2} |z_j|.                          (*)

THREE-IDENTITY SYSTEM:
  (Th 3.1) Root vanishing: sum_{j=0}^{d-1} 1/P'(alpha_j) = 0.
  (Th 3.2) Fiber-bridge:  sum_{j=0}^{d-1} 1/((alpha_j - c) P'(alpha_j)) = -1/P(c).
  (Th 3.3) Fiber vanishing: sum_{j=0}^{d-1} 1/((alpha_j - c)^2 P'(alpha_j)) = 0.

PANDROSION REDUCTION (Theorem 4.1):
At each critical point c:
  R(c) = tilde_q(c) / d
where tilde_q(z) = sum_{k=1}^{d-1} k * a_{d-k} * z^{d-1-k},
with P(z) = z^d + a_{d-1} z^{d-1} + ... + a_1 z, so tilde_q(0) = (d-1) * a_1 = d-1.

Therefore: Smale's MVC for d-degree polys reduces to the algebraic question
"does tilde_q not exceed its value at the origin (= d-1) at every root of P'?"

For d=3, this is closed algebraically. For general d, only first-order
perturbation around z^d - z is established.

VERIFICATION
============

This script verifies:
  1. The Pandrosion inverse: S(c) = |Q(c, 0)|.
  2. The three vanishing identities.
  3. The Pandrosion reduction R(c) = tilde_q(c)/d.
  4. tilde_q(0) = d - 1.
"""
from __future__ import annotations
import math
import numpy as np


def normalize_P(P):
    """Force P(0) = 0, P'(0) = 1."""
    P = P.astype(complex).copy()
    P[-1] = 0.0
    P[-2] = 1.0
    return P


def Q(P, a, z):
    """Pandrosion quotient (P(z) - P(a))/(z - a)."""
    if abs(z - a) < 1e-15:
        return np.polyval(np.polyder(P), z)
    return (np.polyval(P, z) - np.polyval(P, a)) / (z - a)


def Q2(P, c, z):
    """Q_2(c, z) = (P(z) - P(c))/(z - c)^2 (Pandrosion inverse for critical c)."""
    if abs(z - c) < 1e-15:
        # Limit: P''(c)/2
        return np.polyval(np.polyder(P, 2), c) / 2.0
    return (np.polyval(P, z) - np.polyval(P, c)) / (z - c)**2


def co_fiber(P, c):
    """Roots of P(z) - P(c) = 0 other than c (with multiplicity).
    Equivalent to roots of Q_2(c, z) (degree d-2)."""
    d = P.shape[0] - 1
    Pc = np.polyval(P, c)
    poly = P.copy()
    poly[-1] = poly[-1] - Pc  # P(z) - P(c)
    roots = np.roots(poly)
    # Remove c (multiplicity 2 since c is critical: P(z) - P(c) has double root at c)
    remaining = []
    found = 0
    for r in sorted(roots, key=lambda x: abs(x - c)):
        if found < 2 and abs(r - c) < 1e-6:
            found += 1
        else:
            remaining.append(r)
    return remaining


def smale_ratio(P, c):
    """S(c) = |R(c)| = |P(c)/c|."""
    if abs(c) < 1e-12: return None
    return abs(np.polyval(P, c) / c)


def co_fiber_product(P, c):
    """|c| * prod |z_j| where z_j are co-fiber points."""
    fib = co_fiber(P, c)
    if not fib: return abs(c)
    return abs(c) * float(np.prod([abs(z) for z in fib]))


def root_vanishing(P):
    """sum 1/P'(alpha_j) over all roots = 0?"""
    roots = np.roots(P)
    Pp = np.polyder(P)
    s = sum(1.0 / np.polyval(Pp, ak) for ak in roots)
    return abs(s)


def fiber_bridge(P, c):
    """sum 1/((alpha_j - c) P'(alpha_j)) = -1/P(c)?"""
    roots = np.roots(P)
    Pp = np.polyder(P)
    s = sum(1.0 / ((ak - c) * np.polyval(Pp, ak)) for ak in roots)
    Pc = np.polyval(P, c)
    if abs(Pc) < 1e-12: return None
    return abs(s + 1.0 / Pc)


def fiber_vanishing(P, c):
    """sum 1/((alpha_j - c)^2 P'(alpha_j)) = 0?"""
    roots = np.roots(P)
    Pp = np.polyder(P)
    s = sum(1.0 / ((ak - c)**2 * np.polyval(Pp, ak)) for ak in roots)
    return abs(s)


def tilde_q(P, c):
    """tilde_q(z) = sum_{k=1}^{d-1} k * a_{d-k} * z^{d-1-k}.

    With P(z) = z^d + a_{d-1} z^{d-1} + ... + a_1 z (so P[0] = 1, P[-2] = a_1),
    coefficients a_{d-k} for k = 1..d-1 are P[1], P[2], ..., P[-2].

    P[i] in numpy convention is coefficient of z^{d-i}. So a_{d-k} = P[k].
    tilde_q(z) = sum_{k=1}^{d-1} k * P[k] * z^{d-1-k}.
    """
    d = P.shape[0] - 1
    s = 0.0 + 0j
    for k in range(1, d):
        s += k * P[k] * c**(d - 1 - k)
    return s


def pandrosion_reduction(P, c):
    """R(c) vs tilde_q(c)/d."""
    if abs(c) < 1e-12: return None
    Rc = np.polyval(P, c) / c
    tq = tilde_q(P, c)
    d = P.shape[0] - 1
    return abs(Rc - tq / d)


def main():
    print("=" * 80)
    print("PAPER 4 — Pandrosion inverse, fiber identities, Pandrosion reduction")
    print("=" * 80)

    rng = np.random.default_rng(0)

    # 1. Pandrosion inverse: S(c) = |Q(c, 0)|
    print("\n[1] Pandrosion inverse: S(c) = |Q(c, 0)| = |P(c)/c|")
    for d in [3, 5, 7]:
        P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
        P[0] = 1.0
        P = normalize_P(P)
        Pp = np.polyder(P)
        crits = np.roots(Pp)
        for c in crits[:2]:
            if abs(c) < 1e-9: continue
            S_via_R = abs(np.polyval(P, c) / c)
            S_via_Q = abs(Q(P, c, 0))
            print(f"  d={d}: c={c:.3f}, S via R(c) = {S_via_R:.4f}, S via Q(c,0) = {S_via_Q:.4f}")

    # 2. Three vanishing identities
    print("\n[2] Three vanishing identities (Theorem 3.1, 3.2, 3.3)")
    print(f"  {'d':>3} {'root vanish':>14} {'fiber bridge':>14} {'fiber vanish':>14}")
    for d in [3, 4, 5, 6, 8]:
        P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
        P[0] = 1.0
        P = normalize_P(P)
        rv = root_vanishing(P)
        Pp = np.polyder(P)
        crits = np.roots(Pp)
        c = crits[0]  # first critical
        fb = fiber_bridge(P, c)
        fv = fiber_vanishing(P, c)
        print(f"  {d:>3} {rv:>14.2e} {fb if fb is None else f'{fb:>14.2e}':>14} {fv:>14.2e}")

    # 3. Co-fiber product: S(c) = |c| * prod |z_j|
    print("\n[3] Co-fiber product: S(c) = |c| * prod |z_j| (Theorem 2.3)")
    for d in [3, 5, 7]:
        P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
        P[0] = 1.0
        P = normalize_P(P)
        Pp = np.polyder(P)
        crits = np.roots(Pp)
        for c in crits[:2]:
            if abs(c) < 1e-9: continue
            S_direct = smale_ratio(P, c)
            S_cofiber = co_fiber_product(P, c)
            print(f"  d={d}: c={c:.3f}, S direct = {S_direct:.6f}, S cofiber = {S_cofiber:.6f}")

    # 4. Pandrosion reduction: R(c) = tilde_q(c)/d
    print("\n[4] Pandrosion reduction (Theorem 4.1): R(c) = tilde_q(c)/d")
    for d in [3, 4, 5, 6]:
        P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
        P[0] = 1.0
        P = normalize_P(P)
        # tilde_q(0) should be d - 1
        tq0 = tilde_q(P, 0).real
        Pp = np.polyder(P)
        crits = np.roots(Pp)
        max_err = 0.0
        for c in crits:
            if abs(c) < 1e-9: continue
            err = pandrosion_reduction(P, c)
            if err is not None and err > max_err: max_err = err
        print(f"  d={d}: tilde_q(0) = {tq0:.4f} (expected {d-1}), max |R(c) - tq(c)/d| = {max_err:.2e}")


if __name__ == "__main__":
    main()
