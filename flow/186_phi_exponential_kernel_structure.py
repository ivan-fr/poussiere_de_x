"""
PAPER: 186 (NEW — exponential-kernel structure of Riemann Phi)
TITLE: Riemann-Phi after y=e^{4u}: exponential kernels, sign regularity,
       and the next variation-diminishing target
STATUS: RH remains OPEN.  Papers 183-185 rule out several shortcuts.
        This paper rewrites Phi in the variable y=e^{4u}, exposing it as
        a positive linear combination of exponential-kernel pieces with a
        polynomial prefactor.  The next plausible theorem is sign-regularity
        of this special exponential mixture, not generic log-concavity.
DEPENDS: 182 (moment/factorial decomposition), 184 (Phi stability),
         185 (log-concavity insufficiency).

THEORY
======

Riemann's Phi kernel in the de Bruijn-Newman integral can be written

    Phi(u) = sum_{m>=1}
      (2 pi^2 m^4 e^{9u} - 3 pi m^2 e^{5u}) exp(-pi m^2 e^{4u}).

Set y = e^{4u}.  Then

    Phi(u) = sum_{m>=1}
      pi m^2 y^(5/4) (2 pi m^2 y - 3) exp(-pi m^2 y).

Each summand is an exponential kernel exp(-lambda y) times a simple
positive polynomial factor for y >= 1.  Since kernels exp(-lambda y)
are totally positive in (lambda,y), this is the first place where a
variation-diminishing theorem could plausibly enter.

The hard part is the factor (2 lambda y - 3).  This paper tests whether
the full kernel

    K(lambda,y) = y^(5/4) (2 lambda y - 3) exp(-lambda y)

is sign-regular / totally positive on lambda=pi m^2, y>=1.
"""
from __future__ import annotations

import itertools
import math


def det_small(rows):
    n = len(rows)
    total = 0.0
    for perm in itertools.permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
        term = 1.0
        for i in range(n):
            term *= rows[i][perm[i]]
        total += -term if inv % 2 else term
    return total


def K_exp(lam, y):
    return math.exp(-lam * y)


def K_phi_piece(lam, y):
    return (y ** 1.25) * (2 * lam * y - 3.0) * math.exp(-lam * y)


def test_kernel(kernel, lambdas, ys, max_size=4):
    neg = 0
    total = 0
    min_det = float("inf")
    min_case = None
    # For exp(-lambda y), total positivity has a sign depending on ordering.
    # Use lambdas increasing and ys increasing, then det is alternating:
    # (-1)^{k(k-1)/2} det >= 0.
    for size in range(1, max_size + 1):
        sign = -1 if (size * (size - 1) // 2) % 2 else 1
        for li in itertools.combinations(range(len(lambdas)), size):
            for yi in itertools.combinations(range(len(ys)), size):
                M = [[kernel(lambdas[i], ys[j]) for j in yi] for i in li]
                d = sign * det_small(M)
                total += 1
                if d < min_det:
                    min_det = d
                    min_case = (li, yi, size)
                if d < -1e-14:
                    neg += 1
    return total, neg, min_det, min_case


def phi_y(y, m_terms=30):
    return sum(math.pi * m * m * (y ** 1.25) * (2 * math.pi * m * m * y - 3.0)
               * math.exp(-math.pi * m * m * y)
               for m in range(1, m_terms + 1))


def main():
    print("=" * 80)
    print("PAPER 186 — exponential-kernel structure of Riemann Phi")
    print("=" * 80)

    lambdas = [math.pi * m * m for m in range(1, 8)]
    ys = [1.0, 1.05, 1.15, 1.35, 1.7, 2.2, 3.0]

    # ------------------------------------------------------------------
    # [1] Kernel sign-regularity
    # ------------------------------------------------------------------
    print("\n[1] Sign-regularity of exp(-lambda y) baseline")
    total, neg, md, case = test_kernel(K_exp, lambdas, ys)
    print(f"  minors tested: {total}")
    print(f"  sign-regular violations: {neg}")
    print(f"  minimum signed determinant: {md:.6e}, case={case}")

    print("\n[2] Sign-regularity of Phi summand kernel")
    total, neg, md, case = test_kernel(K_phi_piece, lambdas, ys)
    print(f"  minors tested: {total}")
    print(f"  sign-regular violations: {neg}")
    print(f"  minimum signed determinant: {md:.6e}, case={case}")

    # ------------------------------------------------------------------
    # [3] Positivity and shape in y
    # ------------------------------------------------------------------
    print("\n[3] Phi(y) positivity/shape on y>=1")
    pmin = float("inf")
    argmin = 1.0
    pmax = -float("inf")
    argmax = 1.0
    for j in range(501):
        y = 1.0 + 4.0 * j / 500
        v = phi_y(y)
        if v < pmin:
            pmin, argmin = v, y
        if v > pmax:
            pmax, argmax = v, y
    print(f"  min Phi_y on [1,5]: {pmin:.12e} at y={argmin:.4f}")
    print(f"  max Phi_y on [1,5]: {pmax:.12e} at y={argmax:.4f}")

    # ------------------------------------------------------------------
    # [4] Honest assessment
    # ------------------------------------------------------------------
    print("\n[4] HONEST ASSESSMENT")
    print("  WHAT THIS ADDS:")
    print("    Phi is not just an arbitrary positive density.  After y=e^{4u},")
    print("    it is built from exp(-lambda y), a classical totally positive")
    print("    kernel, multiplied by y^(5/4)(2 lambda y - 3).")
    print()
    print("  IF THE PHI-PIECE KERNEL PASSES:")
    print("    The next proof attempt should use sign-regularity of this kernel")
    print("    plus Cauchy-Binet / Andreief to transfer positivity to minors.")
    print()
    print("  IF IT FAILS:")
    print("    The proof must use cancellation across the m-sum, not termwise")
    print("    total positivity.  That is much harder but now precisely located.")


if __name__ == "__main__":
    main()
