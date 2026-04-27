"""
PAPER: 021
TITLE: True Pandrosion 3D via a proportional-chain prism
STATUS: native 3D map, not Steffensen, not Aitken, not h(h)

GOAL
====

The previous "difference prism" is useful, but it is not pure Pandrosion:
it is Aitken-Steffensen applied to the 2D Pandrosion map.

This file tests a stricter 3D construction for the cube root.  The state has
two native coordinates

    (s, t)

representing the two normalized proportional means in the chain

    1, s, t, x^{-1}.

The ideal fixed state is

    beta = x^{-1/3},       (s*, t*) = (beta, beta^2),

and the readout is

    x * t  ->  x^{1/3}.

NATIVE 3D MAP
=============

The 3D Thales sum is not S_3(s)=1+s+s^2.  It is the free chain sum

    C(s,t) = 1 + s + t.

The correction is the same Pandrosion correction applied to the full chain:

    g(s,t) = 1 - (x - 1) / (x C(s,t)).

Then the chain is advanced by one native proportional transport:

    H_3(s,t) = ( g(s,t),  s * g(s,t) ).

This uses one 3D chain sum and one transport.  It does not evaluate h(h), it
does not use finite differences, and it does not apply Aitken/Steffensen.

FIXED POINT
===========

If H_3(s,t)=(s,t), then t=s^2 and

    x (1+s+s^2)(1-s) = x - 1
    x (1-s^3) = x - 1
    s^3 = 1/x.

So the positive fixed point is exactly

    s* = x^{-1/3},       t* = x^{-2/3},

and the readout x*t* is x^{1/3}.

LOCAL STRUCTURE
===============

At the fixed point the Jacobian is generally rank 2:

    J = [[a, a],
         [beta(1+a), beta a]],

where beta=x^{-1/3} and a=(1-beta)/(1+beta+beta^2).

Its determinant is -a*beta, nonzero for x>1.  So this construction avoids the
rank-one symmetric collapse of the simplex/hypercube/cross-polytope attempts.

HONEST VERDICT
==============

This is a true Pandrosion 3D map.  It preserves the cube-root target and has a
native two-coordinate dynamics.  It is not, however, faster than paper-0 on the
scalar cube-root benchmark.  That is useful information: the first strict 3D
Pandrosion geometry exists, but it is not yet the improved algorithm.
"""

from __future__ import annotations

import math


def paper0_h(s: float, x: float) -> float:
    """Paper-0 Pandrosion map for p=3."""
    return 1.0 - (x - 1.0) / (x * (1.0 + s + s * s))


def paper0_readout(s: float, x: float) -> float:
    return x * s * s


def chain_sum_3d(s: float, t: float) -> float:
    return 1.0 + s + t


def pandrosion_3d_step(s: float, t: float, x: float) -> tuple[float, float]:
    """Native 3D Pandrosion chain step H_3(s,t)."""
    g = 1.0 - (x - 1.0) / (x * chain_sum_3d(s, t))
    return g, s * g


def readout_3d(t: float, x: float) -> float:
    return x * t


def canonical_start_3d(x: float) -> tuple[float, float]:
    """
    Start from one uniform first correction and put the chain on t=s^2.
    No root information is used.
    """
    q = 1.0 - (x - 1.0) / (3.0 * x)
    return q, q * q


def iterate_paper0(x: float, tol: float = 1e-12, max_iter: int = 10000) -> tuple[int, float]:
    alpha = x ** (1.0 / 3.0)
    s, _ = canonical_start_3d(x)
    for n in range(max_iter + 1):
        err = abs(paper0_readout(s, x) - alpha)
        if err < tol:
            return n, err
        s = paper0_h(s, x)
    return max_iter, err


def iterate_3d(x: float, tol: float = 1e-12, max_iter: int = 10000) -> tuple[int, float, float, float]:
    alpha = x ** (1.0 / 3.0)
    s, t = canonical_start_3d(x)
    for n in range(max_iter + 1):
        err = abs(readout_3d(t, x) - alpha)
        if err < tol:
            return n, err, s, t
        s, t = pandrosion_3d_step(s, t, x)
    return max_iter, err, s, t


def paper0_lambda(x: float) -> float:
    beta = x ** (-1.0 / 3.0)
    d = 1.0 + beta + beta * beta
    return (x - 1.0) / x * (1.0 + 2.0 * beta) / (d * d)


def jacobian_3d_data(x: float) -> tuple[float, float, float]:
    """Return spectral radius and the two eigenvalues of H_3 at the fixed point."""
    beta = x ** (-1.0 / 3.0)
    a = (1.0 - beta) / (1.0 + beta + beta * beta)
    trace = a * (1.0 + beta)
    disc = trace * trace + 4.0 * a * beta
    lam_plus = 0.5 * (trace + math.sqrt(disc))
    lam_minus = 0.5 * (trace - math.sqrt(disc))
    rho = max(abs(lam_plus), abs(lam_minus))
    return rho, lam_plus, lam_minus


def fixed_point_residual(x: float, s: float, t: float) -> float:
    beta = x ** (-1.0 / 3.0)
    return max(abs(s - beta), abs(t - beta * beta))


def main() -> None:
    print("=" * 86)
    print("flow/021 -- true Pandrosion 3D: proportional-chain prism")
    print("=" * 86)
    print()
    print("Map:")
    print("  C(s,t) = 1+s+t")
    print("  g(s,t) = 1 - (x-1)/(x C(s,t))")
    print("  H_3(s,t) = (g(s,t), s*g(s,t))")
    print()
    print("Fixed point:")
    print("  H_3(s,t)=(s,t) => t=s^2 and x(1+s+s^2)(1-s)=x-1")
    print("  hence s*=x^(-1/3), t*=x^(-2/3), readout x*t*=x^(1/3)")
    print()

    cases = [2.0, 3.0, 5.0, 8.0, 10.0, 100.0]
    print(
        f"{'x':>8} {'paper0 it':>10} {'3D it':>8} {'lambda p0':>11} "
        f"{'rho H3':>11} {'eig +':>11} {'eig -':>11} {'fp err':>11}"
    )
    print("-" * 86)
    for x in cases:
        n0, _ = iterate_paper0(x)
        n3, err3, s3, t3 = iterate_3d(x)
        rho, ep, em = jacobian_3d_data(x)
        fp_err = fixed_point_residual(x, s3, t3)
        print(
            f"{x:8.3g} {n0:10d} {n3:8d} {paper0_lambda(x):11.6f} "
            f"{rho:11.6f} {ep:11.6f} {em:11.6f} {fp_err:11.3e}"
        )

    print()
    print("Verdict:")
    print("  YES: this is a strict 3D Pandrosion map, not Steffensen.")
    print("  YES: the cube-root fixed point is preserved exactly.")
    print("  YES: the local Jacobian has nonzero determinant for x>1, so no rank-1 collapse.")
    print("  NO: on this scalar benchmark it is not faster than paper-0.")
    print("  Next target: find a 3D transport law that keeps the same fixed point")
    print("  but lowers rho(H_3) below the paper-0 multiplier.")


if __name__ == "__main__":
    main()
