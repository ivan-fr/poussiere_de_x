"""
PAPER: 085 (canonical: 75pandrosion_christol.pdf)
TITLE: Christol's Theorem and Pandrosion Automaton
STATUS: proved (Christol 1979)

THEORY
======

CHRISTOL'S THEOREM: A formal power series f in F_p[[t]] is algebraic over
F_p(t) iff its coefficient sequence (c_n) is computable by a finite
automaton (p-automatic).

PANDROSION CONNECTION: The Pandrosion field of an algebraic f over F_p(t)
gives constraints on the coefficient sequence.

VERIFICATION
============

  1. Algebraic series example.
  2. p-automatic sequence.
"""
from __future__ import annotations


def main():
    print("=" * 80)
    print("PAPER 85 — Christol's theorem (algebraic = automatic)")
    print("=" * 80)

    print("\n[1] Algebraic series example: f = 1 / (1 - t) in F_2[[t]]")
    print("  f = 1 + t + t^2 + t^3 + ...")
    print("  Coefficients (c_n) = (1, 1, 1, 1, ...) — clearly automatic.")

    print("\n[2] More complex: f^2 + t f + 1 = 0 in F_2[[t]]")
    print("  This defines an algebraic series. Coefficients are 2-automatic")
    print("  by Christol's theorem.")

    print("\n[3] Thue-Morse sequence: 2-automatic, gives an algebraic series in F_2")
    # Thue-Morse: t_n = (sum of bits of n) mod 2
    print("  Thue-Morse t_n: 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, ...")
    # Verify
    seq = [bin(n).count('1') % 2 for n in range(20)]
    print(f"  First 20 terms: {seq}")
    print("  Generating function f(t) = sum t_n t^n is algebraic over F_2(t).")


if __name__ == "__main__":
    main()
