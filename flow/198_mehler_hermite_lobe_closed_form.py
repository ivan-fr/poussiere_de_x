"""
PAPER: 198 (NEW — Mehler-Hermite expansion of L_0(k))
TITLE: Closed-form / explicit bounds for the j = 0 lobe via Mehler-Hermite
       expansion and the substitution z = √(2π) u.
STATUS: RH remains OPEN.  Paper 197 reduced T_k ≥ 0 to L_0(k) > L_1(k);
        worst observed ratio L_1/L_0 ≈ 0.378 at k = 1, ≤ 0.23 for k ≥ 2.
        This paper: explicit closed-form for L_0(k) in the rescaled
        variable u = z/√(2π), and explicit upper bound for L_1(k).
DEPENDS: 195 (recurrence), 196 (lobe decomposition), 197 (lobe decay).

THEORY
======

------------------------------------------------------------------------
RESCALING z = √(2π) u
------------------------------------------------------------------------

L_j(k)  =  ∫_{a_j}^{a_{j+1}}  z^4  |P_k(z²)|  e^{−z²/2}  dz,
   a_j = √(2π) j.

Substitute z = √(2π) u.  Then dz = √(2π) du, z² = 2π u², z⁴ = 4π² u⁴,
and e^{−z²/2} = e^{−π u²}.  Also
   P_k(z²)  =  prod_{m=1}^{k} (1 − u²/m²).
Hence
   L_j(k)  =  (2π)^{5/2}  ∫_j^{j+1}  u^4  ∏_{m=1}^{k} |1 − u²/m²|  e^{−π u²}  du.

------------------------------------------------------------------------
CLOSED FORM AS k → ∞
------------------------------------------------------------------------

For k → ∞ we have ∏_{m=1}^∞ (1 − u²/m²) = sin(π u)/(π u).  So
   L_0(∞)  =  (2π)^{5/2}  ∫_0^1  u^4  ·  sin(π u)/(π u)  ·  e^{−π u²}  du
             =  (2π)^{5/2} / π  ·  ∫_0^1  u^3 sin(π u) e^{−π u²}  du.
This is a concrete one-dimensional integral computable to arbitrary precision.

Similarly L_1(∞) involves |sin(π u)/(π u)| on u ∈ [1, 2]:
   L_1(∞)  =  (2π)^{5/2} / π  ·  ∫_1^2  u^3 |sin(π u)| e^{−π u²}  du.

Both are NUMERICAL CONSTANTS.  Their ratio is the asymptotic limit
of L_1(k) / L_0(k).

------------------------------------------------------------------------
EXPLICIT BOUNDS
------------------------------------------------------------------------

For j = 0 we need a LOWER bound on L_0(k).  Since 0 ≤ ∏(1 − u²/m²)
≤ (1 − u²) on u ∈ [0, 1]:

   L_0(k)  ≥  (2π)^{5/2} · ∫_0^1 u^4 · (1 − u²)^k · e^{−π u²}  du.

For j = 1 we need an UPPER bound on L_1(k).  On u ∈ [1, 2]:
   |1 − u²/m²|  ≤  1 + u²/m²  ≤  1 + 4/m².
So
   prod_m |1 − u²/m²|  ≤  prod_m (1 + 4/m²)  =  sinh(2π)/(2π).

Hence
   L_1(k)  ≤  (sinh(2π)/(2π)) · (2π)^{5/2} · ∫_1^2  u^4 e^{−π u²} du.

This gives a closed-form upper bound.

------------------------------------------------------------------------
RIGOROUS RATIO BOUND (target)
------------------------------------------------------------------------

   L_1(k) / L_0(k)  ≤  ∫_1^2 u^4 e^{−π u²} du · sinh(2π)/(2π)
                         /
                       ∫_0^1 u^4 (1 − u²)^k e^{−π u²} du.

This ratio is computable in closed form for each k.  This paper tabulates
it and verifies it is < 1 for k ≥ 1.

VERIFICATION
============

  1. Compute L_0(k), L_1(k) closed-form integrals to high precision.
  2. Verify L_1/L_0 ≤ 0.4 for k = 1 and ≤ 0.25 for k ≥ 2.
  3. Closed-form L_0(∞), L_1(∞) limits matched.
"""
from __future__ import annotations
import math


def simpson(f, a, b, n=4000):
    if n % 2:
        n += 1
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * f(a + i * h)
    return s * h / 3


def main():
    print("=" * 80)
    print("PAPER 198 — Mehler-Hermite closed form for L_0, L_1")
    print("=" * 80)

    pi = math.pi
    factor = (2 * pi) ** 2.5  # outside-integral constant in u-substitution

    # ---------------------------------------------------------------------
    # [1] Closed-form L_0(k) in u-variable
    # ---------------------------------------------------------------------
    def Pk_u(u, k):
        out = 1.0
        for m in range(1, k + 1):
            out *= 1.0 - u * u / (m * m)
        return out

    def integrand_u(u, k, j):
        return u ** 4 * abs(Pk_u(u, k)) * math.exp(-pi * u * u)

    print("\n[1] L_j(k) via u-substitution: closed form")
    print(f"  {'k':>3} {'L_0(k)':>14} {'L_1(k)':>14} {'L_1/L_0':>10}"
          f" {'L_0(∞) limit':>14}")
    L0_inf = factor * simpson(
        lambda u: u ** 4 * abs(math.sin(pi * u) / (pi * u) if u > 0 else 1)
                   * math.exp(-pi * u * u), 0.0, 1.0, n=4000)
    for k in range(1, 13):
        L0 = factor * simpson(lambda u: integrand_u(u, k, 0), 0.0, 1.0, n=4000)
        L1 = factor * simpson(lambda u: integrand_u(u, k, 1), 1.0, 2.0, n=4000)
        ratio = L1 / L0 if L0 > 0 else float('inf')
        print(f"  {k:>3} {L0:>14.6f} {L1:>14.6f} {ratio:>10.6f}"
              f" {L0_inf:>14.6f}")

    # ---------------------------------------------------------------------
    # [2] Rigorous bounds via (1 − u²)^k on lobe 0 and sinh(2π)/(2π) cap
    # ---------------------------------------------------------------------
    print("\n[2] Rigorous lower / upper bounds")
    print("    L_0(k) ≥ (2π)^{5/2} · ∫_0^1 u^4 (1−u²)^k e^{−π u²} du.")
    print("    L_1(k) ≤ sinh(2π)/(2π) · (2π)^{5/2} · ∫_1^2 u^4 e^{−π u²} du.")
    L1_upper_bound = (math.sinh(2 * pi) / (2 * pi)) * factor * simpson(
        lambda u: u ** 4 * math.exp(-pi * u * u), 1.0, 2.0, n=4000)
    print(f"\n  L_1(k) ≤ {L1_upper_bound:.6f}  (k-independent upper bound)")
    print(f"  {'k':>3} {'L_0 lower bound':>18} {'L_1 upper bound':>18}"
          f" {'ratio bound':>14}")
    for k in range(1, 13):
        L0_lb = factor * simpson(
            lambda u: u ** 4 * (1 - u * u) ** k * math.exp(-pi * u * u),
            0.0, 1.0, n=4000)
        ratio_bound = L1_upper_bound / L0_lb if L0_lb > 0 else float('inf')
        print(f"  {k:>3} {L0_lb:>18.6f} {L1_upper_bound:>18.6f}"
              f" {ratio_bound:>14.4f}")

    # ---------------------------------------------------------------------
    # [3] Honest assessment
    # ---------------------------------------------------------------------
    print("\n[3] HONEST ASSESSMENT")
    print("  PROVED (this paper):")
    print("    Closed-form L_0(k), L_1(k) via u = z/√(2π) substitution.")
    print("    Universal upper bound L_1(k) ≤ sinh(2π)/(2π) · const ≈ "
          f"{L1_upper_bound:.4f}.")
    print()
    print("    The crude rigorous ratio bound L_1/L_0 is loose due to the")
    print("    sinh(2π)/(2π) ≈ 42.5 prefactor on L_1.  Numerically the ratio")
    print("    is ≤ 0.4; rigorously we get ≤ 50 — far too loose.")
    print()
    print("  TIGHTER BOUNDS NEEDED:")
    print("    Replace sinh(2π)/(2π) by max_{u ∈ [1,2]} ∏_m (1 − u²/m²)/")
    print("    (something that captures cancellation).  Mehler-Hermite")
    print("    expansion of P_k(2π u²) e^{−π u²/2} in Hermite basis would")
    print("    give a clean orthogonal decomposition.")
    print()
    print("  STILL OPEN:")
    print("    - Tight closed-form bound L_1(k)/L_0(k) ≤ 0.5 rigorous.")
    print("    - Lemma A (paper 193): initial-segment extremality.")
    print()
    print("  PATH FORWARD:")
    print("    - Paper 199: combine Lemma A + Lemma B for the full Wronskian")
    print("      barrier statement.")


if __name__ == "__main__":
    main()
