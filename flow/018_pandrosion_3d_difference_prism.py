"""
PAPER: 018
TITLE: Pandrosion 3D via the difference prism
STATUS: new 3D geometry; derivative-free quadratic acceleration

GEOMETRY
========

The failed 3D/dD attempts replace the 2D Pandrosion sum S_p by a
multivariate sum/product.  In the symmetric case this either collapses to
rank one, becomes separable, or changes the fixed point.

This file tests a different 3D lift.

Keep the original 2D Thales map

    h(s) = 1 - (x - 1) / (x S_p(s)).

The third dimension is not another symmetric spatial coordinate.  It is a
phase-space coordinate carrying the finite difference between consecutive
Thales layers:

    s0 = s,
    s1 = h(s0),
    s2 = h(s1),
    d0 = s1 - s0,
    d1 = s2 - s1.

In the 3D difference prism the two points

    Q0 = (d0, s0, 0),   Q1 = (d1, s1, 1)

define a ruled chord.  Project this chord to the fixed plane d = 0.  The
intercept is

    Pi_3D(s) = s0 - d0^2 / (d1 - d0)

which is exactly Aitken-Steffensen applied to the Pandrosion map.

This is a genuine 3D geometric construction: two stacked Thales sections plus
one transverse fixed-plane projection.  It preserves the original fixed point
s* = x^(-1/p), avoids symmetric rank-one spatial collapse, and improves the
local order from 1 to 2 without derivatives.
"""

from __future__ import annotations

import math


def S_p(s: float, p: int) -> float:
    return sum(s**k for k in range(p))


def pandrosion_h(s: float, x: float, p: int) -> float:
    return 1.0 - (x - 1.0) / (x * S_p(s, p))


def pandrosion_readout(s: float, x: float, p: int) -> float:
    return x * s ** (p - 1)


def prism_3d_step(s: float, x: float, p: int) -> float:
    """One 3D difference-prism projection step."""
    s1 = pandrosion_h(s, x, p)
    s2 = pandrosion_h(s1, x, p)
    d0 = s1 - s
    d1 = s2 - s1
    denom = d1 - d0
    if abs(denom) < 1e-300:
        return s2
    return s - d0 * d0 / denom


def iterate_raw(x: float, p: int, s0: float, tol: float = 1e-12, max_iter: int = 10000):
    target = x ** (1.0 / p)
    s = s0
    history = []
    for n in range(max_iter + 1):
        err = abs(pandrosion_readout(s, x, p) - target)
        history.append(err)
        if err < tol:
            return n, s, history
        s = pandrosion_h(s, x, p)
    return max_iter, s, history


def iterate_prism(x: float, p: int, s0: float, tol: float = 1e-12, max_iter: int = 100):
    target = x ** (1.0 / p)
    s = s0
    history = []
    for n in range(max_iter + 1):
        err = abs(pandrosion_readout(s, x, p) - target)
        history.append(err)
        if err < tol:
            return n, s, history
        s = prism_3d_step(s, x, p)
    return max_iter, s, history


def estimate_order(errors: list[float]) -> float:
    usable = [e for e in errors if math.isfinite(e) and 1e-16 < e < 1e-2]
    if len(usable) < 3:
        return float("nan")
    e0, e1, e2 = usable[-3:]
    if e0 == 0 or e1 == 0 or e2 == 0 or e0 == e1:
        return float("nan")
    return math.log(e2 / e1) / math.log(e1 / e0)


def canonical_start(x: float, p: int) -> float:
    return pandrosion_h(1.0, x, p)


def main() -> None:
    print("=" * 78)
    print("Pandrosion 3D: difference-prism geometry")
    print("=" * 78)
    print()
    print("3D lift:")
    print("  two Thales layers produce s0, s1=h(s0), s2=h(s1)")
    print("  the transverse coordinate is the finite difference d=s_next-s")
    print("  projecting the chord (d0,s0)->(d1,s1) to d=0 gives")
    print("      Pi_3D(s) = s - (h(s)-s)^2 / (h(h(s))-2h(s)+s)")
    print()

    cases = [
        (2.0, 3),
        (8.0, 3),
        (10.0, 3),
        (100.0, 3),
        (2.0, 5),
        (10.0, 5),
    ]

    print(f"{'x':>10} {'p':>3} {'s0':>13} {'raw iters':>10} {'3D iters':>9} {'order 3D':>10} {'err 3D':>11}")
    print("-" * 78)
    for x, p in cases:
        s0 = canonical_start(x, p)
        n_raw, _, raw_hist = iterate_raw(x, p, s0)
        n_prism, s_prism, prism_hist = iterate_prism(x, p, s0)
        order = estimate_order(prism_hist)
        final_err = abs(pandrosion_readout(s_prism, x, p) - x ** (1.0 / p))
        raw_label = f"{n_raw}" if n_raw < 10000 else ">10000"
        print(
            f"{x:10.3g} {p:3d} {s0:13.9f} {raw_label:>10} {n_prism:9d} "
            f"{order:10.3f} {final_err:11.3e}"
        )

    print()
    print("Verdict:")
    print("  This is not another symmetric polytope.  The new axis stores")
    print("  the residual direction d=s_next-s, so the fixed point is the")
    print("  geometric plane d=0.  The 3D projection preserves x^(1/p) and")
    print("  gives quadratic convergence without derivatives.")


if __name__ == "__main__":
    main()
