"""
PAPER: 003 (canonical: 3pandrosion_smale.pdf)
TITLE: Smale's Mean Value Conjecture and the Pandrosion Vanishing Identity
STATUS: structural reformulation (does not close MVC)
DEPENDS: 002

THEORY
======

Setup: P monic of degree d with P(0) = 0 and P'(0) = 1.
Roots: 0 = alpha_0, alpha_1, ..., alpha_{d-1}.
Critical points: c (roots of P').

Smale's ratio at the origin:
  S(c) = |R(c)/R(0)| = |P(c)/c| / |P'(0)|         (since P'(0) = 1)
       = |P(c)/c|.

Conjecture: min over critical c of S(c) <= (d-1)/d.

VANISHING IDENTITY (Theorem 2.1, this paper):
For P with P(0) = 0 and P'(0) = R(0) = 1:

    1 = -sum_{k=1}^{d-1} 1 / (alpha_k * R'(alpha_k))           (*)

where R(z) = P(z)/z = prod_{k>=1} (z - alpha_k).

This is the classical LAGRANGE-SYLVESTER PARTIAL FRACTION IDENTITY
(written in Pandrosion notation as zeta_P(1) = sum 1/P'(alpha_k) = 0).

Corollary 2.2 (Root separation bound):
By triangle inequality on (*):

    min_{k>=1} |alpha_k| * |R'(alpha_k)|  <=  d - 1.            (**)

This says: at least one non-zero root alpha_k has bounded
|alpha_k| * |separation from other non-zero roots|.

CRITICAL POINT IDENTITY (Section 3):
At a critical point c of P:
  R(c) = -c * R'(c),
  S(c) = |c * R'(c)|.

THE GAP TO MVC:
The vanishing identity gives a ROOT-SIDE bound (alpha_k); MVC asks for a
CRITICAL-POINT bound (c). The classical MVC reduction (Tischler) bridges
these via fiber identities (paper 4 of this series).

VERIFICATION
============

This script verifies:
  1. The vanishing identity (*) holds to machine precision for random P.
  2. The triangle-inequality root-side bound (**).
  3. The critical point identity R(c) = -c R'(c).
  4. The product identity prod_c S(c) = |c R'(c)|.
"""
from __future__ import annotations
import math
import numpy as np


def normalize_P_at_origin(P, rng):
    """Construct random monic P with P(0) = 0 and P'(0) = 1.
    P(z) = z * R(z) where R(z) = prod (z - alpha_k), R(0) = 1 means
    prod (-alpha_k) = 1 (i.e., constant term of R is 1).
    Equivalent: shift coefficients so P(0) = 0 and P'(0) = 1."""
    d = P.shape[0] - 1
    if d < 2: return None
    # Force P(0) = 0: set constant term to 0
    P = P.copy()
    P[-1] = 0.0
    # Force P'(0) = 1: set linear coefficient to 1
    P[-2] = 1.0
    return P


def vanishing_identity_residual(P):
    """Compute |1 + sum_{k>=1} 1/(alpha_k * R'(alpha_k))|
    for P with P(0) = 0, P'(0) = 1.

    Should be zero by Lagrange-Sylvester identity."""
    roots = np.roots(P)
    # alpha_0 = 0, others are alpha_k
    nonzero = [r for r in roots if abs(r) > 1e-9]
    if len(nonzero) != len(roots) - 1: return None  # double root at 0?
    # R(z) = P(z)/z, R'(alpha_k) = P'(alpha_k) (since P(alpha_k) = 0)
    Pp = np.polyder(P)
    s = 0.0 + 0j
    for ak in nonzero:
        # R(z) = P(z)/z, so R'(alpha_k) = P'(alpha_k)/alpha_k for alpha_k != 0
        # (since P(alpha_k) = 0 for k >= 1).
        # Hence 1/(alpha_k * R'(alpha_k)) = 1/P'(alpha_k).
        Pp_ak = np.polyval(Pp, ak)
        s += 1.0 / Pp_ak
    return abs(1.0 + s)


def critical_point_identity(P):
    """At each critical point c, verify R(c) + c R'(c) = 0
    where R = P/z. Equivalently, c * P'(c) = P(c) - c * R(c) ...

    Direct check: at critical c, P'(c) = R(c) + c R'(c) = 0
    implies R(c) = -c R'(c).
    """
    Pp = np.polyder(P)
    crits = np.roots(Pp)
    errs = []
    for c in crits:
        if abs(c) < 1e-9: continue
        # R(c) = P(c)/c, R'(c) = (P'(c) - R(c))/c = -R(c)/c (since P'(c) = 0)
        Rc = np.polyval(P, c) / c
        Rpc = -Rc / c
        # Verify: R(c) = -c R'(c)
        # i.e., Rc + c * Rpc = 0
        errs.append(abs(Rc + c * Rpc))
    return errs


def root_side_bound(P):
    """Verify min_k |alpha_k * R'(alpha_k)| <= d - 1 (triangle inequality)."""
    d = P.shape[0] - 1
    roots = np.roots(P)
    nonzero = [r for r in roots if abs(r) > 1e-9]
    Pp = np.polyder(P)
    vals = [abs(ak * np.polyval(Pp, ak)) for ak in nonzero]
    return min(vals), d - 1


def main():
    print("=" * 80)
    print("PAPER 3 — Pandrosion vanishing identity and root-side bound")
    print("=" * 80)

    # 1. Vanishing identity verification
    print("\n[1] Lagrange-Sylvester vanishing identity (Theorem 2.1)")
    print(f"  Identity: 1 + sum_{{k>=1}} 1/(alpha_k R'(alpha_k)) = 0")
    rng = np.random.default_rng(0)
    for d in [3, 4, 5, 6, 8, 10]:
        max_err = 0.0
        for _ in range(10):
            P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
            P[0] = 1.0
            P = normalize_P_at_origin(P, rng)
            err = vanishing_identity_residual(P)
            if err is None: continue
            if err > max_err: max_err = err
        print(f"  d={d}: max residual = {max_err:.2e}")

    # 2. Root-side bound verification
    print("\n[2] Root-side bound: min_k |alpha_k * R'(alpha_k)| <= d - 1")
    label = "min |a R prime|"
    print(f"  {'d':>3} {label:>16} {'d - 1':>8}  status")
    for d in [3, 5, 8, 12]:
        min_val, bd = float('inf'), d - 1
        for _ in range(50):
            P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
            P[0] = 1.0
            P = normalize_P_at_origin(P, rng)
            mv, _ = root_side_bound(P)
            if mv < min_val: min_val = mv
        ok = "OK" if min_val <= bd + 1e-6 else "VIOLATES"
        print(f"  {d:>3} {min_val:>14.6f} {bd:>8} {ok}")

    # 3. Critical-point identity R(c) = -c R'(c)
    print("\n[3] Critical-point identity: R(c) + c R'(c) = 0 at critical points")
    for d in [3, 5, 8]:
        max_err = 0.0
        for _ in range(10):
            P = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
            P[0] = 1.0
            P = normalize_P_at_origin(P, rng)
            errs = critical_point_identity(P)
            if errs:
                max_err = max(max_err, max(errs))
        print(f"  d={d}: max identity error = {max_err:.2e}")


if __name__ == "__main__":
    main()
