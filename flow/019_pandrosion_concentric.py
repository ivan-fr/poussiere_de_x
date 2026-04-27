"""
PAPER: 019 (alternative geometry: concentric circles / spheres / d-spheres)
TITLE: Pandrosion via Thales-on-circles (family eta)
STATUS: empirical exploration

CONSTRUCTION
============
Paper-0 uses parallel LINES to slice a rectangle.  Thales also works with
concentric CIRCLES: a secant line through the center of concentric circles
cuts each circle at a radius proportional to that circle's radius.

  2D annulus:  concentric circles of radius r_0 < r_1 < ... < r_p slice the
               annulus into rings of area A_k = pi (r_k^2 - r_{k-1}^2).
  3D shell:    concentric spheres r_0 < ... < r_p slice into spherical
               shells of volume V_k = (4pi/3)(r_k^3 - r_{k-1}^3).
  d-shell:     concentric (d-1)-spheres slice into shells of volume
               V_k = c_d (r_k^d - r_{k-1}^d).

POLYNOMIAL DERIVATION (RADIAL ONLY):
Imposing geometric progression V_k = s^{k-1} V_1 with u_k = r_k^d / R^d:
    u_k - u_{k-1} = s^{k-1} (u_1 - u_0).
This is LINEAR in u, identical to paper-0 with reparametrisation u = r^d.
So purely radial Thales-on-spheres reduces to paper-0 strict.
==> SAME as family beta (pyramid).  Degenerate.

BUT: spheres carry ANGULAR structure that cubes/pyramids do not.
The d-sphere has area element sin^{d-2}(theta) dtheta dphi_1 ... dphi_{d-2}.
A latitude band at theta has weight sin^{d-2}(theta), giving NON-UNIFORM
weights to angular Thales slices.  This is a TRUE new polynomial.

POLYNOMIAL DERIVATION (RADIAL + ANGULAR, d=3):
Sphere S^2: area element sin(theta) dtheta dphi.
Slice angularly at theta_k = k * pi / p.  Latitude band area:
    A_k = 2 pi (cos theta_{k-1} - cos theta_k).
With u_k = 1 - cos theta_k (in [0, 2]) and geometric progression of A_k:
    A_k / A_1 = (u_k - u_{k-1}) / (u_1 - u_0) = s^{k-1}.
Linear in u again ==> same paper-0 after substitution u = 1 - cos theta.

So the radial AND angular Thales on spheres separately give paper-0.
The TWO can be combined into a hypercube-like structure with two
variables (s_r, s_theta).

POLYNOMIAL DERIVATION (FULL POLAR/SPHERICAL):
2D polar (r, theta):  S_eta,2(s_r, s_theta) = S_p(s_r) * S_p(s_theta).
3D spherical:         S_eta,3(s_r, s_theta, s_phi) = S_p * S_p * S_p.
This is structurally identical to family epsilon (hypercube).

GENUINELY NEW POLYNOMIAL: weighted angular Thales on S^{d-1}.
If we slice angularly at NON-UNIFORM theta to enforce equal-weight rings:
    sin^{d-2}(theta) Delta theta = const,
or equivalently, slice such that V_k = s^{k-1} V_1 with V_k a sin^{d-2}-
weighted measure, the polynomial S becomes:
    S_eta,d^{ang}(s) = sum_k s^k * w_k(d),
where w_k(d) reflects the spherical weighting.  This is the trapezoid
family (zeta) with d-dependent weights: it shifts the FIXED POINT to
a new algebraic value, NOT x^{1/p}.

VERDICT (PRELIMINARY):
- Concentric d-spheres + radial Thales = paper-0 (degenerate).
- Concentric d-spheres + uniform angular Thales = paper-0 again.
- Concentric d-spheres + (radial AND angular) = hypercube product.
- Concentric d-spheres + sin^{d-2}-weighted angular = trapezoid analog
  (different fixed point, different problem).
- No new path to a faster scalar x^{1/p} algorithm.

Verify numerically.
"""
from __future__ import annotations
import math


def S_p(s, p):
    return sum(s**k for k in range(p))


def S_radial_dsphere(s, d, p):
    """Radial-only Thales on d-sphere: u_k = r_k^d, S identical to paper-0."""
    return sum(s**k for k in range(p))


def S_angular_dsphere(s, d, p, n_quad=2000):
    """Angular Thales on S^{d-1} with sin^{d-2} weighting.

    Slice theta in [0, pi] into p bands of equal angular extent pi/p.
    Each band has weight w_k = int_{theta_{k-1}}^{theta_k} sin^{d-2} dtheta.
    Polynomial: S = sum_k s^k * w_k / w_0.
    """
    if d == 2:
        # S^1: weight is uniform, returns paper-0 S_p.
        return sum(s**k for k in range(p))
    # Numerical integration of sin^{d-2} on each band.
    weights = []
    for k in range(p):
        t0 = k * math.pi / p
        t1 = (k + 1) * math.pi / p
        # Trapezoidal rule
        N = n_quad
        h = (t1 - t0) / N
        total = 0.0
        for i in range(N + 1):
            t = t0 + i * h
            val = math.sin(t) ** (d - 2)
            if i == 0 or i == N:
                total += 0.5 * val
            else:
                total += val
        weights.append(total * h)
    w0 = weights[0]
    return sum((complex(s) ** k) * (weights[k] / w0) for k in range(p))


def S_polar_product(s_r, s_t, p):
    """Polar/spherical full product (hypercube-like)."""
    return S_p(complex(s_r), p) * S_p(complex(s_t), p)


def iter_radial(x, d, p, max_iter=200, tol=1e-13):
    """Radial-only Thales on d-sphere: identical to paper-0."""
    s = 0.5
    for n in range(max_iter):
        S = S_p(s, p)
        s_new = 1 - (x - 1) / (x * S)
        if abs(s_new - s) < tol:
            return s_new, n + 1
        s = s_new
    return s, max_iter


def iter_angular(x, d, p, max_iter=200, tol=1e-13):
    """Angular Thales on S^{d-1} with sin-weighting (trapezoid-like)."""
    s = 0.5 + 0.0j
    for n in range(max_iter):
        S = S_angular_dsphere(s, d, p)
        if abs(S) < 1e-300:
            return s, n
        s_new = 1 - (x - 1) / (x * S)
        if abs(s_new - s) < tol:
            return s_new, n + 1
        s = s_new
    return s, max_iter


def iter_polar(x_r, x_t, p, max_iter=200, tol=1e-13):
    """Polar/spherical: 2 variables (s_r, s_theta), hypercube product."""
    s_r, s_t = 0.5 + 0j, 0.5 + 0j
    for n in range(max_iter):
        S = S_polar_product(s_r, s_t, p)
        if abs(S) < 1e-300:
            return (s_r, s_t), n
        s_r_new = 1 - (x_r - 1) / (x_r * S)
        s_t_new = 1 - (x_t - 1) / (x_t * S)
        if abs(s_r_new - s_r) < tol and abs(s_t_new - s_t) < tol:
            return (s_r_new, s_t_new), n + 1
        s_r, s_t = s_r_new, s_t_new
    return (s_r, s_t), max_iter


def main():
    print("=" * 76)
    print("Family eta: Pandrosion via THALES-ON-CIRCLES (concentric d-spheres)")
    print("=" * 76)
    print()
    print("Three sub-constructions:")
    print("  (R) radial slicing of d-shell with parameter u = r^d")
    print("  (A) angular slicing of S^{d-1} with sin^{d-2}-weighted bands")
    print("  (P) polar/spherical product of radial AND angular")
    print()

    # SUB-CONSTRUCTION (R): radial Thales on d-spheres
    print("[R] RADIAL Thales on d-spheres:")
    print(f"  {'d':>3} {'x':>5} {'p':>3} {'s_*':>14} {'paper-0 x^{-1/p}':>20} {'iter':>5}")
    print("-" * 76)
    for d in [2, 3, 4, 5, 6]:
        for p in [2, 3]:
            x = 2.0
            s_star, n_iter = iter_radial(x, d, p)
            s_pred = x ** (-1.0 / p)
            mark = "OK" if abs(s_star - s_pred) < 1e-9 else "DIFF"
            print(f"  {d:>3} {x:>5.1f} {p:>3} {s_star:>14.10f} {s_pred:>20.10f} {n_iter:>5}  {mark}")
    print("  --> radial collapses to paper-0 strict (DEGENERATE).")
    print()

    # SUB-CONSTRUCTION (A): angular Thales with sin-weighting
    print("[A] ANGULAR Thales on S^{d-1} with sin^{d-2}-weighting:")
    print(f"  {'d':>3} {'x':>5} {'p':>3} {'s_*':>14} {'paper-0 x^{-1/p}':>20} {'iter':>5}")
    print("-" * 76)
    for d in [2, 3, 4, 5, 6]:
        for p in [2, 3]:
            x = 2.0
            s_star, n_iter = iter_angular(x, d, p)
            s_pred = x ** (-1.0 / p)
            sr = s_star.real
            diff = abs(sr - s_pred)
            mark = "==paper-0" if diff < 1e-9 else f"shift {diff:.3e}"
            print(f"  {d:>3} {x:>5.1f} {p:>3} {sr:>14.10f} {s_pred:>20.10f} {n_iter:>5}  {mark}")
    print("  --> for d=2, identical to paper-0; for d>=3, sin^{d-2}-weight")
    print("      shifts fixed point ==> trapezoid-like (DIFFERENT problem).")
    print()

    # SUB-CONSTRUCTION (P): polar/spherical product
    print("[P] POLAR product (s_r, s_theta) [hypercube-like]:")
    print(f"  {'p':>3} {'x_r':>5} {'x_t':>5} {'s_r*':>14} {'s_t*':>14} {'iter':>5}")
    print("-" * 76)
    for p in [2, 3]:
        for (xr, xt) in [(2.0, 2.0), (2.0, 3.0), (3.0, 5.0)]:
            (s_r, s_t), n_iter = iter_polar(xr, xt, p)
            print(f"  {p:>3} {xr:>5.1f} {xt:>5.1f} "
                  f"{s_r.real:>14.10f} {s_t.real:>14.10f} {n_iter:>5}")
    print("  --> equivalent to family epsilon (hypercube): rank-1 collapse")
    print("      at symmetric x_r = x_t, genuinely coupled at asymmetric.")
    print()

    # ANGULAR weights for inspection
    print("Angular weight pattern w_k(d) / w_0(d) for p=4:")
    print(f"  {'d':>3}   {'w_0':>10} {'w_1':>10} {'w_2':>10} {'w_3':>10}")
    for d in [2, 3, 4, 5, 6]:
        # raw weights
        p = 4
        weights = []
        for k in range(p):
            t0 = k * math.pi / p
            t1 = (k + 1) * math.pi / p
            N = 2000
            h = (t1 - t0) / N
            total = 0.0
            for i in range(N + 1):
                t = t0 + i * h
                val = math.sin(t) ** (d - 2) if d >= 2 else 1.0
                if i == 0 or i == N:
                    total += 0.5 * val
                else:
                    total += val
            weights.append(total * h)
        w0 = weights[0]
        rs = [w / w0 for w in weights]
        print(f"  {d:>3}   {rs[0]:>10.4f} {rs[1]:>10.4f} {rs[2]:>10.4f} {rs[3]:>10.4f}")
    print("  d=2: uniform (paper-0).  d>=3: weights peak at equator (k=p/2).")
    print()

    print("=" * 76)
    print("VERDICT (family eta = Thales-on-circles/spheres):")
    print("  [R] radial only         ==> paper-0 strict (degenerate)")
    print("  [A] angular sin-weighted==> trapezoid analog (different x^{1/p})")
    print("  [P] radial X angular    ==> hypercube product (already studied)")
    print()
    print("  No new geometric path to a faster scalar x^{1/p} algorithm.")
    print("  However, [A] is interesting in its own right: it gives a")
    print("  PARAMETRISED FAMILY indexed by ambient dimension d, with")
    print("  fixed points that are new algebraic numbers (roots of weighted")
    print("  S_p^{eta,d}(s) * (1 - s) = 1 - 1/x).  These could be useful")
    print("  for solving SPHERICALLY POSED root problems (not the goal here).")
    print("=" * 76)


if __name__ == "__main__":
    main()
