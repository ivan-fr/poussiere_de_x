"""
PAPER: 196 (NEW — lobe dominance for T_k)
TITLE: Gaussian lobe decomposition for
       T_k = E[Z^4 prod_{m<=k}(1-Z^2/(2*pi*m^2))]
STATUS: RH remains OPEN.  Paper 195 reduced monotonicity of the initial
        block barrier to T_k >= 0.  This paper rewrites T_k as an
        alternating sum over Gaussian lobes separated by the zeros
        sqrt(2*pi)m and tests a lobe-dominance certificate.
DEPENDS: 195 (initial-block recurrence).

THEORY
======

Let
    a_m = sqrt(2*pi) m,
    P_k(z^2) = prod_{m=1}^k (1 - z^2/a_m^2).

Then

    T_k = E[Z^4 P_k(Z^2)]
        = 2/sqrt(2*pi) * integral_0^infty
            z^4 P_k(z^2) exp(-z^2/2) dz.

The sign of P_k(z^2) alternates on intervals

    [0,a_1], [a_1,a_2], ..., [a_k,infty).

Define lobe magnitudes L_j >= 0 by integrating the absolute value on
each interval.  Then

    T_k = C * (L_0 - L_1 + L_2 - ... + (-1)^k L_k).

A possible proof of T_k>=0 is a dominance criterion for the alternating
lobe sum.  This paper computes the lobes and cumulative tail sums.
"""
from __future__ import annotations

import math


def Pk_z2(z, k):
    out = 1.0
    z2 = z * z
    for m in range(1, k + 1):
        out *= 1.0 - z2 / (2 * math.pi * m * m)
    return out


def density_integrand_abs(z, k):
    return abs((z ** 4) * Pk_z2(z, k) * math.exp(-0.5 * z * z))


def simpson(f, a, b, n=800):
    if n % 2:
        n += 1
    h = (b - a) / n
    total = f(a) + f(b)
    for i in range(1, n):
        total += (4 if i % 2 else 2) * f(a + i * h)
    return total * h / 3


def lobe_magnitudes(k):
    roots = [math.sqrt(2 * math.pi) * m for m in range(1, k + 1)]
    bounds = [0.0] + roots
    lobes = []
    for j in range(k):
        a, b = bounds[j], bounds[j + 1]
        lobes.append(simpson(lambda z: density_integrand_abs(z, k), a, b))
    # Tail: Gaussian is tiny; integrate to max(root_k + 12, 20).
    tail_end = max(bounds[-1] + 12.0, 20.0)
    lobes.append(simpson(lambda z: density_integrand_abs(z, k), bounds[-1], tail_end, n=1600))
    return lobes


def Tk_from_lobes(lobes):
    c = 2.0 / math.sqrt(2 * math.pi)
    return c * sum(((-1) ** j) * L for j, L in enumerate(lobes))


def Tk_coeff(k):
    """Coefficient/moment exact-ish computation from paper 195."""
    coeffs = [1.0]
    for m in range(1, k + 1):
        x = 1.0 / (2 * math.pi * m * m)
        coeffs.append(0.0)
        for j in range(len(coeffs) - 1, 0, -1):
            coeffs[j] -= x * coeffs[j - 1]

    def odd_df(n):
        out = 1.0
        for q in range(n, 0, -2):
            out *= q
        return out

    return sum(c * odd_df(2 * j + 3) for j, c in enumerate(coeffs))


def main():
    print("=" * 80)
    print("PAPER 196 — Gaussian lobe dominance for T_k")
    print("=" * 80)

    # ------------------------------------------------------------------
    # [1] Lobe decomposition
    # ------------------------------------------------------------------
    print("\n[1] Alternating lobe decomposition")
    print(f"  {'k':>3} {'T_lobes':>14} {'T_coeff':>14} {'min tail sum':>14}")
    for k in range(1, 13):
        lobes = lobe_magnitudes(k)
        T_lobes = Tk_from_lobes(lobes)
        T_exact = Tk_coeff(k)
        # Tail sums S_j = L_j - L_{j+1}+...
        tail_sums = []
        for j in range(len(lobes)):
            tail_sums.append(sum(((-1) ** (i - j)) * lobes[i]
                                 for i in range(j, len(lobes))))
        print(f"  {k:>3} {T_lobes:>14.8e} {T_exact:>14.8e}"
              f" {min(tail_sums):>14.8e}")

    # ------------------------------------------------------------------
    # [2] Lobe profile examples
    # ------------------------------------------------------------------
    print("\n[2] Lobe profiles")
    for k in [3, 6, 10]:
        lobes = lobe_magnitudes(k)
        print(f"  k={k}:")
        print("    ", " ".join(f"{L:.3e}" for L in lobes[: min(len(lobes), 8)]))
        print(f"    alternating sum raw = {sum(((-1)**j)*L for j,L in enumerate(lobes)):.8e}")

    # ------------------------------------------------------------------
    # [3] Honest assessment
    # ------------------------------------------------------------------
    print("\n[3] HONEST ASSESSMENT")
    print("  PROGRESS:")
    print("    T_k is an alternating Gaussian-lobe sum with explicit zeros at")
    print("    sqrt(2*pi)m.  Numerical lobes reproduce the coefficient T_k.")
    print()
    print("  POSSIBLE PROOF ROUTE:")
    print("    Prove all alternating tail sums of lobe magnitudes are positive.")
    print("    That would imply T_k>=0 directly by a variation-diminishing")
    print("    argument on the Gaussian-weighted polynomial.")
    print()
    print("  LIMITATION:")
    print("    This is still analytic/numeric evidence; a rigorous lobe bound")
    print("    needs explicit estimates on each interval.")


if __name__ == "__main__":
    main()
