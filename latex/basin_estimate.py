"""
Numerically estimate the radius R_sigma and basin volume lower bound
for the unseeded Steffensen iteration sigma_{p,x} at the principal
anchor gamma_0 = x^{-1/p}.

For (p, x) = (3, 2):
  - Compute the largest R such that every point in B(alpha, R)
    is driven to alpha by plain sigma iteration (no seed).
  - Compare to the analytic lower bound from the quadratic-contraction
    radius construction.
  - Estimate empirical basin fractions on a grid.
"""
from __future__ import annotations
import math
import numpy as np


def h_fun(z: complex, p: int = 3, x: float = 2.0) -> complex:
    Sp = sum(z ** k for k in range(p))
    if Sp == 0:
        return float("nan") + 1j * float("nan")
    return 1 - (x - 1) / (x * Sp)


def sigma_fun(z: complex, p: int = 3, x: float = 2.0) -> complex:
    hz = h_fun(z, p, x)
    hhz = h_fun(hz, p, x)
    denom = hhz - 2 * hz + z
    if abs(denom) < 1e-16:
        return hz  # take h(z) if Aitken denominator vanishes
    return z - (hz - z) ** 2 / denom


def converges_to(z0: complex, target: complex,
                 max_iters: int = 200, tol: float = 1e-10,
                 blowup: float = 1e4) -> bool:
    z = complex(z0)
    for _ in range(max_iters):
        z_new = sigma_fun(z)
        if not np.isfinite(z_new.real) or not np.isfinite(z_new.imag):
            return False
        if abs(z_new) > blowup:
            return False
        if abs(z_new - target) < tol:
            return True
        z = z_new
    return abs(z - target) < 1e-4


def classify(z0: complex, anchors: list[complex],
             max_iters: int = 200, tol: float = 1e-8) -> int:
    """Return the index of the anchor closest to the final iterate,
    or -1 if divergent."""
    z = complex(z0)
    for _ in range(max_iters):
        z_new = sigma_fun(z)
        if not np.isfinite(z_new.real) or not np.isfinite(z_new.imag):
            return -1
        if abs(z_new) > 1e4:
            return -1
        if any(abs(z_new - a) < tol for a in anchors):
            return min(range(len(anchors)), key=lambda k: abs(z_new - anchors[k]))
        z = z_new
    return min(range(len(anchors)), key=lambda k: abs(z - anchors[k]))


def main() -> None:
    p, x = 3, 2.0
    alpha = x ** (-1.0 / p)
    omega = math.e ** (2j * math.pi / p)
    anchors = [alpha * omega ** k for k in range(p)]
    print(f"(p, x) = ({p}, {x})")
    print(f"Anchors: {[f'{a:.4f}' for a in anchors]}")

    # (1) Estimate R_sigma by bisection
    print("\n--- Estimating R_sigma (largest R s.t. every w in B(alpha, R) converges to alpha)")
    N_boundary = 180
    angles = np.linspace(0, 2 * math.pi, N_boundary, endpoint=False)
    R_lo, R_hi = 0.0, 0.5
    for _ in range(30):
        R = (R_lo + R_hi) / 2
        all_conv = all(
            converges_to(alpha + R * np.exp(1j * th), alpha)
            for th in angles
        )
        if all_conv:
            R_lo = R
        else:
            R_hi = R
    R_sigma_emp = R_lo
    print(f"  Empirical R_sigma >= {R_sigma_emp:.4f}")
    print(f"  Lower-bound on basin area: pi R_sigma^2 >= {math.pi * R_sigma_emp**2:.4f}")

    # (2) Grid sweep: estimate basin fractions
    print("\n--- Grid sweep for basin fractions on [-1.5, 1.5]^2")
    N = 300
    xs = np.linspace(-1.5, 1.5, N)
    ys = np.linspace(-1.5, 1.5, N)
    counts = [0, 0, 0, 0]  # 0, 1, 2, -1 (divergent)
    for xi in xs:
        for yi in ys:
            k = classify(xi + 1j * yi, anchors)
            counts[k] += 1
    total = N * N
    print(f"  Fraction -> anchor 0 (principal, real): {counts[0]/total:.4f}")
    print(f"  Fraction -> anchor 1:                   {counts[1]/total:.4f}")
    print(f"  Fraction -> anchor 2:                   {counts[2]/total:.4f}")
    print(f"  Fraction divergent/unclassified:        {counts[3]/total:.4f}")

    square_area = 3.0 * 3.0
    analytic_lb_per_anchor = math.pi * R_sigma_emp ** 2
    total_analytic_lb = p * analytic_lb_per_anchor / square_area
    print(f"\nAnalytic lower bound on total basin fraction in this square:")
    print(f"  p * pi R_sigma^2 / area([-1.5, 1.5]^2) = "
          f"{total_analytic_lb:.4f}  ({100 * total_analytic_lb:.2f}%)")
    print(f"Empirical total basin fraction: "
          f"{(counts[0]+counts[1]+counts[2])/total:.4f}")


if __name__ == "__main__":
    main()
