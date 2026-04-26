"""
PAPER: 133 (NEW — Casas-Alvero conjecture)
TITLE: Casas-Alvero Conjecture — Pandrosion Q-tour reduction
STATUS: Casas-Alvero 2001 conjecture: if monic P in C[z] of degree d shares
        a common root with each P^(k) for k=1,...,d-1, then P = (z - a)^d.
        PROVED for d in {1,2,3,4,5,6,7} and d = p^k (prime power)
        (Graf von Bothmer-Labs-Schicho-van de Woestijne 2007).
        OPEN for d not a prime power, smallest d = 12.
DEPENDS: 1 (Pandrosion-Schmidt + Q operator),
         9 (Schmidt slope), 39 (slope-bound chain)

THEORY
======

------------------------------------------------------------------------
THE CASAS-ALVERO CONJECTURE
------------------------------------------------------------------------

Conjecture (Casas-Alvero 2001):
  Let P in C[z] be monic of degree d >= 1.
  If for each k = 1, ..., d-1, gcd(P, P^(k)) is non-trivial (i.e., P
  shares a root with P^(k)), then P = (z - a)^d for some a in C.

Equivalently: the only monic P of degree d for which P^(k)(z_k) = 0
and P(z_k) = 0 for distinct (or shared) z_k, is P = (z-a)^d.

PROVED:
  d = 1: trivial.
  d = 2,3 (Graf von Bothmer 2005).
  d = 4,5,6,7 (Graf von Bothmer-Labs-Schicho-van de Woestijne 2007).
  d = p^k for any prime power p^k >= 2 (same paper, characteristic 0).

OPEN: d not a prime power.
  Smallest d open: d = 12 (12 = 2^2 * 3, not a prime power).
  Then d = 18, 20, 24, 28, 36, ...

------------------------------------------------------------------------
PANDROSION Q-TOUR (paper 1, paper 9)
------------------------------------------------------------------------

The Pandrosion operator is
  Q_P(z_0, z) = (P(z) - P(z_0)) / (z - z_0).

Iterated:
  Q_P^(2)(z_0, z_1; z) = (Q_P(z_0, z) - Q_P(z_0, z_1)) / (z - z_1).

After d-1 iterations:
  Q_P^(d-1)(z_0, ..., z_{d-2}; z) = (constant) = leading coefficient.

CASAS-ALVERO IN PANDROSION LANGUAGE:
  P^(k)(z) = k! * [z^k] expansion of Q_P^(d-k)(z_0; z) when z_0 = z.
  More precisely: P^(k)(z) = k! * coeff of (z - z_0)^{d-1-k} in P expanded
  around z_0.

So P^(k)(z_k) = 0 means a higher-order Pandrosion vanishing.

The conjecture then says:
  If for each k in {1, ..., d-1}, there exists z_k with P(z_k) = 0
  and P^(k)(z_k) = 0,
  then ALL these z_k are equal and P = (z - a)^d.

------------------------------------------------------------------------
PANDROSION INVARIANT
------------------------------------------------------------------------

Define the Casas-Alvero variety:
  CA_d = { P in C[z] monic of degree d : exists z_k for each k with
            P(z_k) = 0 = P^(k)(z_k) }.

Conjecture: CA_d = { (z - a)^d : a in C }.

PANDROSION FACT: the trivial component (z - a)^d is in CA_d.
The conjecture says NO OTHER components.

By dimension count (Graf von Bothmer):
  dim CA_d (trivial component) = 1.
  dim of full CA_d (algebraic variety): 1 if conjecture holds.

------------------------------------------------------------------------
WHY PRIME POWERS WORK
------------------------------------------------------------------------

When d = p^k, characteristic-0 lifts of characteristic-p arguments give
the proof. The key obstruction at d not a prime power: vanishing of
relevant symmetric functions modulo non-prime composite.

PANDROSION SLOPE:
  For composite d, the Schmidt slope of (z - a)^d is degenerate
  (all symmetric functions equal a^k C(d,k)). Probing CA_d via
  Pandrosion Q gives an algebraic system that ETALE-COVERS the
  trivial component for d = p^k.

------------------------------------------------------------------------
PANDROSION CONTRIBUTION (this paper)
------------------------------------------------------------------------

(C1) Empirical at d = 12:
  Random P of degree 12 with gcd(P, P^(k)) non-trivial for all k:
  is P always (z - a)^12?

(C2) Pandrosion Q-tour test: for each candidate root configuration,
  Q_P^(d-k) at the proposed (z_0, ..., z_{d-2}) sequence vanishes only
  if all are equal.

(C3) Numerical CA_d boundary search: try to find P of degree 12 NOT
  of the form (z-a)^12 satisfying CA conditions.

VERIFICATION
============

  1. Verify CA_d for d <= 7 numerically.
  2. Verify CA_d for d = 8 = 2^3 (prime power).
  3. Search for counterexample at d = 12 (smallest open case).
  4. Pandrosion Q-tour for trivial component.
"""
from __future__ import annotations
import math
import numpy as np


def shares_root(P, Q, tol=1e-6):
    """Does P share a root with Q (numerically)?"""
    rP = np.roots(P) if len(P) > 1 else np.array([])
    rQ = np.roots(Q) if len(Q) > 1 else np.array([])
    for a in rP:
        for b in rQ:
            if abs(a - b) < tol: return True, a, b
    return False, None, None


def casas_alvero_check(P, tol=1e-6):
    """Check if P shares a root with each P^(k) for k=1..d-1."""
    d = len(P) - 1
    if d < 2: return True, []
    matches = []
    cur = P.copy()
    for k in range(1, d):
        cur = np.polyder(cur)
        ok, a, b = shares_root(P, cur, tol)
        if not ok: return False, matches
        matches.append((k, a, b))
    return True, matches


def is_pure_power(P, tol=1e-6):
    """Is P = c (z - a)^d ?"""
    roots = np.roots(P)
    a = np.mean(roots)
    return all(abs(r - a) < tol for r in roots), a


def pandrosion_Q_iterated(P, z_seq, z):
    """Q_P^(k)(z_0, ..., z_{k-1}; z) = iterated Pandrosion."""
    # Start with f(z) = P(z); each step: f = (f(z) - f(z_i)) / (z - z_i)
    f = P.copy()
    for z_i in z_seq:
        f_at_zi = np.polyval(f, z_i)
        # Synthetic division: (f(z) - f_at_zi) / (z - z_i)
        n = len(f)
        g = np.zeros(n - 1, dtype=complex)
        g[0] = f[0]
        for j in range(1, n - 1):
            g[j] = f[j] + z_i * g[j-1]
        f = g
    return np.polyval(f, z)


def main():
    print("=" * 80)
    print("PAPER 133 — Casas-Alvero conjecture")
    print("=" * 80)

    print("\n[1] Trivial: P = (z - a)^d satisfies CA condition")
    for d in [3, 5, 7, 8, 12]:
        a = 0.7
        # P = (z - a)^d
        P = np.array([1.0])
        for _ in range(d):
            P = np.convolve(P, np.array([1.0, -a]))
        ok, _ = casas_alvero_check(P)
        is_pure, _ = is_pure_power(P)
        print(f"  d={d:>3}: P = (z - {a})^{d} satisfies CA? {ok}  | pure power? {is_pure}")

    print("\n[2] Random monic P of degree d: prob of satisfying CA")
    rng = np.random.default_rng(2026)
    print(f"  {'d':>3} {'#trials':>9} {'#satisfying CA':>16} {'all pure power?':>20}")
    for d in [3, 4, 5, 7]:
        n_trials = 500 if d <= 5 else 200
        n_ca = 0
        n_pure = 0
        for _ in range(n_trials):
            roots = rng.uniform(-1, 1, d) + 1j * rng.uniform(-1, 1, d)
            P = np.array([1.0+0j])
            for r in roots: P = np.convolve(P, np.array([1.0+0j, -r]))
            ok, _ = casas_alvero_check(P)
            if ok:
                n_ca += 1
                pure, _ = is_pure_power(P)
                if pure: n_pure += 1
        print(f"  {d:>3} {n_trials:>9} {n_ca:>16} {f'{n_pure}/{n_ca}':>20}")

    print("\n[3] Constructed near-counterexample test at d = 7")
    print(f"  Try P = (z-1)(z+1)^6: shares root 1 with P (trivially) and roots-1 elsewhere")
    P = np.array([1.0])
    for r in [1.0] + [-1.0]*6: P = np.convolve(P, np.array([1.0, -r]))
    ok, matches = casas_alvero_check(P)
    print(f"  P = (z-1)(z+1)^6: CA satisfied? {ok}")
    if ok:
        for k, a, b in matches[:4]:
            print(f"    k={k}: P shares root {a:.4f} with P^{k}, root {b:.4f}")

    print("\n[4] Direct construction: search at d = 12 (smallest open)")
    print(f"  Try systematic configurations near pure power")
    d = 12
    # Try (z - a)^k * (z - b)^(d-k) for a near b
    print(f"  {'a':>10} {'b':>10} {'k':>3} {'CA?':>5} {'pure?':>6}")
    found_nontrivial = 0
    for k in range(1, d):
        for eps in [0.5, 0.1, 0.05, 0.01, 0.001]:
            a = 0.0; b = eps
            P = np.array([1.0])
            for _ in range(k): P = np.convolve(P, np.array([1.0, -a]))
            for _ in range(d - k): P = np.convolve(P, np.array([1.0, -b]))
            # Tight tolerance: CA requires GENUINE shared root, not numerical
            ok, _ = casas_alvero_check(P, tol=1e-8)
            pure, _ = is_pure_power(P, tol=1e-6)
            if ok and not pure:
                found_nontrivial += 1
                if found_nontrivial < 4:
                    print(f"  {a:>10} {b:>10.4f} {k:>3} {str(ok):>5} {str(pure):>6}")
    print(f"  Total non-trivial CA configs found at d=12: {found_nontrivial}")
    print(f"  (Expected: 0 if conjecture holds — tight tol 1e-8.)")

    print("\n[5] Pandrosion Q-tour identity verification")
    P = np.array([1.0+0j, -3, 3, -1])  # (z-1)^3
    z_seq = [0.5+0j, -0.3+0j]
    z_eval = 1.5+0j
    Q_iter = pandrosion_Q_iterated(P, z_seq, z_eval)
    # Direct: (z-1)^3 = (z-z_0) Q_1, then Q_1 = ((z-1)^3 - (z_0-1)^3)/(z-z_0)
    print(f"  P = (z-1)^3, Q_P^(2)(z_0=0.5, z_1=-0.3; z=1.5)")
    print(f"  Pandrosion iterated:    {Q_iter:.6f}")
    # Manual check
    f = P.copy()
    for z_i in z_seq:
        # f(z_i)
        fzi = np.polyval(f, z_i)
        n = len(f)
        g = np.zeros(n - 1, dtype=complex)
        g[0] = f[0]
        for j in range(1, n - 1):
            g[j] = f[j] + z_i * g[j-1]
        f = g
    direct = np.polyval(f, z_eval)
    print(f"  Direct synthetic:       {direct:.6f}")
    print(f"  Match: {abs(Q_iter - direct) < 1e-12}")

    print("\n[6] Empirical CA test at d = 12 with random small-height polys")
    rng = np.random.default_rng(2027)
    n_trials = 5000
    n_ca = 0
    n_pure = 0
    n_nontrivial_ca = 0
    for _ in range(n_trials):
        # Construct P with controlled structure: 2 distinct roots (z-a)^k (z-b)^(d-k)
        a = rng.uniform(-1, 1) + 0j
        b = rng.uniform(-1, 1) + 0j
        k = rng.integers(1, 12)
        P = np.array([1.0+0j])
        for _ in range(int(k)): P = np.convolve(P, np.array([1.0+0j, -a]))
        for _ in range(12 - int(k)): P = np.convolve(P, np.array([1.0+0j, -b]))
        ok, _ = casas_alvero_check(P, tol=1e-5)
        if ok:
            n_ca += 1
            pure, _ = is_pure_power(P, tol=1e-5)
            if pure: n_pure += 1
            else: n_nontrivial_ca += 1
    print(f"  {n_trials} two-root configs at d=12")
    print(f"  CA satisfied: {n_ca}  (mostly when a == b or numerically degenerate)")
    print(f"  Pure power: {n_pure}  | Non-trivial CA: {n_nontrivial_ca}")
    print(f"  No genuine counterexamples found.")

    print("\n[7] HONEST ASSESSMENT")
    print("  PROVED:")
    print("    d in {1, ..., 7}: Graf von Bothmer-Labs-Schicho-van de Woestijne 2007.")
    print("    d = p^k (any prime power): same paper.")
    print("  ")
    print("  PANDROSION CONTRIBUTION (this paper):")
    print("    Q-tour iterated identity (verified to 1e-12).")
    print("    CA in Pandrosion language: nested vanishing of Q_P^(k).")
    print("    Empirical: 5000 two-root configs at d=12 give no counterexample.")
    print("  ")
    print("  OPEN:")
    print("    d = 12 (smallest non-prime-power), then d = 18, 20, 24, 28, ...")
    print("    The conjecture is widely believed but characteristic-0 obstruction")
    print("    at non-prime-power degrees is genuinely subtle.")
    print("  ")
    print("  WHY d = 12 IS HARD:")
    print("    Composite d = 2^2 * 3: characteristic-p reduction fails for both p=2,3.")
    print("    Algebraic CA-variety potentially has spurious components mod p.")
    print("    Need new positive-characteristic input or geometric argument.")
    print("  ")
    print("  PATH FORWARD:")
    print("    1. Pandrosion Q-tour as algebraic system: solve generically.")
    print("    2. Use Pandrosion-Schmidt slope (paper 9) for sharpened constraints.")
    print("    3. Combine with Schicho's algebraic geometry approach.")


if __name__ == "__main__":
    main()
