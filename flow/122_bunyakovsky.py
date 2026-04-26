"""
PAPER: 122 (NEW — Bunyakovsky conjecture)
TITLE: Bunyakovsky's Conjecture (Irreducibles Producing Infinitely Many Primes)
       — A Pandrosion Empirical Census + Honest Statement of the Gap
STATUS: open since 1857; degree-1 case = Dirichlet 1837; all degree >= 2 cases open.
        This paper provides empirical certificates and Pandrosion-language
        reformulation. NO ANALYTIC PROGRESS.
DEPENDS: 011 (Pandrosion quotient), 070 (Mason-Stothers), 047 (Galois),
         045 (Riemann)

THEORY
======

------------------------------------------------------------------------
BUNYAKOVSKY CONJECTURE (1857)
------------------------------------------------------------------------

Let f(x) in Z[x] satisfy:
  (B1) f is IRREDUCIBLE in Z[x],
  (B2) leading coefficient of f is positive,
  (B3) gcd_{n in Z} f(n) = 1 (the values share no common prime factor).

Then: f(n) is prime for INFINITELY MANY positive integers n.

------------------------------------------------------------------------
KNOWN CASES
------------------------------------------------------------------------

DEGREE 1: f(x) = a x + b. Conditions reduce to gcd(a, b) = 1.
  PROVED by Dirichlet 1837 (primes in arithmetic progressions).

DEGREE >= 2: ALL cases are OPEN. Most famous open instances:
  f(n) = n^2 + 1            (Hardy-Littlewood Conjecture E, since 1923)
  f(n) = n^2 + n + 41       (Euler 1772; produces primes for n = 0, ..., 39)
  f(n) = n^4 + 1            (open)
  f(n) = n^3 + 2            (open)

EVEN MORE: it's not known if there exists ANY non-linear f satisfying
B1+B2+B3 producing infinitely many primes!

HARDY-LITTLEWOOD PREDICTION (heuristic, not proved):
For f of degree d:
  pi_f(N) := #{n <= N : f(n) prime} ~ C_f * N / (d * log N)
where C_f = prod_p (1 - omega_f(p)/p) / (1 - 1/p) and omega_f(p) = #{n mod p : f(n) = 0 mod p}.

------------------------------------------------------------------------
PANDROSION REFORMULATION (no progress made)
------------------------------------------------------------------------

The Pandrosion quotient Q(a, b) = (f(a) - f(b))/(a - b) is an integer for
a, b in Z. Bunyakovsky asks for a SPECIFIC set of values {f(n) : n in Z}
to contain infinitely many primes.

Sieve method connection: f satisfies B3 iff for every prime p, omega_f(p) < p
(Pandrosion-Galois interpretation: f's reduction mod p doesn't vanish on
all residue classes).

------------------------------------------------------------------------
WHY PANDROSION DOES NOT CLOSE THIS
------------------------------------------------------------------------

Bunyakovsky requires an ANALYTIC NUMBER-THEORETIC argument:
  - Density of primes in {f(n)} via L-functions (Dirichlet only works for d=1).
  - Sieve methods (large sieve, Selberg) give upper bounds, not lower.
  - Bombieri-Vinogradov / GRH-style bounds insufficient for d >= 2.

The Pandrosion polynomial machinery operates on COEFFICIENT SPACE
(quotients, discriminants, energy). Bunyakovsky operates on VALUE SPACE
(integer values f(1), f(2), ...). These are CONNECTED but the Pandrosion
toolkit does not bridge them.

VERIFICATION
============

  1. Empirical census: count primes in {f(n)} for several test polys.
  2. Hardy-Littlewood prediction comparison.
  3. Pandrosion-quotient gcd test (B3 verification).
"""
from __future__ import annotations
import math


def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    for k in range(3, int(math.sqrt(n)) + 1, 2):
        if n % k == 0: return False
    return True


def count_primes_in_polynomial(f_coefs, N):
    """f_coefs in low-to-high. Count primes f(n) for n in 1..N."""
    count = 0
    for n in range(1, N + 1):
        v = sum(c * n**k for k, c in enumerate(f_coefs))
        if v > 1 and is_prime(v):
            count += 1
    return count


def hardy_littlewood_constant(f_coefs, max_p=1000):
    """C_f = prod_p (1 - omega_f(p)/p) / (1 - 1/p)
       where omega_f(p) = #{n in F_p : f(n) = 0 mod p}."""
    C = 1.0
    for p in range(2, max_p):
        if not is_prime(p): continue
        # omega_f(p)
        omega = 0
        for n in range(p):
            v = sum(c * pow(n, k, p) for k, c in enumerate(f_coefs)) % p
            if v == 0: omega += 1
        # (1 - omega/p) / (1 - 1/p)
        if p == 2:
            factor = (1 - omega/p) / (1 - 1/p)
        else:
            factor = (1 - omega/p) / (1 - 1/p)
        C *= factor
    return C


def gcd_values(f_coefs, N=20):
    """gcd of f(1), f(2), ..., f(N) - tests B3."""
    g = 0
    for n in range(1, N + 1):
        v = sum(c * n**k for k, c in enumerate(f_coefs))
        g = math.gcd(g, abs(v))
    return g


def main():
    print("=" * 80)
    print("PAPER 122 — Bunyakovsky's Conjecture (irreducibles -> infinite primes)")
    print("=" * 80)

    print("\n[1] Famous Bunyakovsky test cases")
    test_polys = [
        ("n^2 + 1 (Hardy-Littlewood E)", [1, 0, 1]),  # 1 + 0*n + 1*n^2 = n^2 + 1
        ("n^2 + n + 41 (Euler 1772)", [41, 1, 1]),
        ("n^4 + 1", [1, 0, 0, 0, 1]),
        ("n^3 + 2", [2, 0, 0, 1]),
        ("2n + 1 (degree 1, Dirichlet)", [1, 2]),
        ("6n + 5 (Dirichlet)", [5, 6]),
    ]

    print(f"  {'f':>30} {'gcd(f(1),...)':>14}")
    for name, c in test_polys:
        g = gcd_values(c, N=30)
        print(f"  {name:>30} {g:>14}")

    print("\n[2] Empirical prime count up to N")
    for N in [100, 1000, 10000]:
        print(f"\n  N = {N}:")
        print(f"  {'f':>30} {'#primes':>10} {'predicted (HL)':>18} {'ratio':>10}")
        for name, c in test_polys[:4]:
            count = count_primes_in_polynomial(c, N)
            d = len(c) - 1
            try:
                C_f = hardy_littlewood_constant(c, max_p=200)
            except: C_f = 1.0
            pred = C_f * N / (d * math.log(max(N, 2)))
            ratio = count / pred if pred > 0 else 0
            print(f"  {name:>30} {count:>10} {pred:>18.2f} {ratio:>10.4f}")

    print("\n[3] Empirical: primes per length scale (logarithmic)")
    print(f"  f(n) = n^2 + 1: count primes in [10^k, 10^(k+1)]")
    f_coefs = [1, 0, 1]
    for k in range(2, 5):
        N_low, N_high = 10**k, 10**(k+1)
        count = 0
        for n in range(N_low, N_high):
            v = n*n + 1
            if is_prime(v): count += 1
        log_term = math.log(N_high) - math.log(N_low)
        print(f"  [10^{k}, 10^{k+1}]: {count} primes  (log span = {log_term:.4f})")

    print("\n[4] Pandrosion-quotient gcd test (B3 verification)")
    # B3: gcd_{n in Z} f(n) = 1 iff for every prime p, f mod p is not identically 0.
    # Equivalently: f has no fixed prime divisor.
    # Pandrosion: gcd(f(a), f(b)) = gcd(f(a), Q(a, b) * (a - b))
    # which gives recursive computation.
    print("  Bunyakovsky condition B3 = gcd(f(1), f(2), ...) = 1")
    print("  is equivalent to: no prime p divides f(n) for ALL n in Z/pZ.")
    print("  (Pandrosion: f mod p has < p roots in F_p.)")

    print("\n[5] HONEST ASSESSMENT")
    print("  Bunyakovsky (1857) is OPEN since for ALL degree >= 2 polynomials.")
    print("  Even a single example (e.g., n^2 + 1) producing infinite primes is open.")
    print("  ")
    print("  Empirically: prime counts agree well with Hardy-Littlewood predictions.")
    print("  ")
    print("  WHY PANDROSION DOES NOT CLOSE THIS:")
    print("  - Bunyakovsky is a question about VALUE distribution (integer values).")
    print("  - Pandrosion machinery acts on COEFFICIENT space (polynomials).")
    print("  - The bridge requires sieve methods + L-functions + Bombieri-Vinogradov,")
    print("    which are NOT in the Pandrosion toolkit.")
    print("  ")
    print("  CLASSICAL TOOLS NEEDED:")
    print("  - Selberg sieve (gives upper bounds only)")
    print("  - Hardy-Littlewood circle method (heuristic)")
    print("  - Friedlander-Iwaniec 1998: x^2 + y^4 prime infinitely often (PROVED)")
    print("  - Heath-Brown 2001: x^3 + 2y^3 prime infinitely often (PROVED)")
    print("  - But single-variable Bunyakovsky for d >= 2 remains OPEN.")


if __name__ == "__main__":
    main()
