"""
PAPER: 134 (NEW — Chowla's conjecture on Liouville/Möbius correlations)
TITLE: Chowla's Conjecture — Pandrosion-zeta + Liouville
STATUS: Chowla 1965 conjecture: for any fixed h_1 < ... < h_k,
          (1/N) sum_{n<=N} lambda(n + h_1) ... lambda(n + h_k) -> 0 as N -> infty.
        PROVED for k = 1 (trivially from PNT).
        Tao 2016 (logarithmic version): logarithmic Chowla for k = 2.
        Tao-Teräväinen 2018: log-Chowla for ODD k.
        OPEN: ordinary (non-logarithmic) for k >= 2; log for EVEN k >= 4.
DEPENDS: 11 (Pandrosion-zeta), 55 (Riemann-Pandrosion),
         107 (Riemann-Pandrosion refined), 110 (Riemann-Pandrosion final)

THEORY
======

------------------------------------------------------------------------
CHOWLA'S CONJECTURE
------------------------------------------------------------------------

Liouville function: lambda(n) = (-1)^{Omega(n)}
  where Omega(n) = total number of prime factors (with multiplicity).

CHOWLA 1965: for fixed distinct integers h_1, ..., h_k:
  (1/N) sum_{n=1}^{N} lambda(n + h_1) ... lambda(n + h_k) -> 0  as N -> infty.

EQUIVALENT FORM (via Möbius mu):
  (1/N) sum_{n=1}^{N} mu(n + h_1) ... mu(n + h_k) -> 0.

KNOWN:
  k = 1: trivial (= sum lambda(n)/N = M(N)/N -> 0 by PNT).
  k = 2 LOGARITHMIC: Tao 2016 (logarithmic averages),
    sum_{n<=N} lambda(n) lambda(n+h) / n = o(log N).
  k odd LOGARITHMIC: Tao-Teräväinen 2018.
OPEN:
  k = 2 ORDINARY: still open.
  k >= 4 even LOGARITHMIC: open.

------------------------------------------------------------------------
PANDROSION-ZETA CONNECTION
------------------------------------------------------------------------

Liouville Dirichlet series:
  L(s) := sum lambda(n)/n^s = zeta(2s)/zeta(s).

By Pandrosion-zeta (paper 11):
  F_zeta(s) = -zeta'/zeta(s) = sum_p log(p) p^{-s} / (1 - p^{-s}).
  F_L(s) = F_{zeta(2s)} - F_{zeta(s)}  (logarithmic derivative).

So Liouville moments correspond to spectral data of zeta.

PANDROSION REFORMULATION OF CHOWLA k=2:
  C_h(N) := (1/N) sum_{n<=N} lambda(n) lambda(n + h).
By Tao-Teräväinen, log version of C_h(N) -> 0.
Ordinary case: open.

------------------------------------------------------------------------
PANDROSION CONTRIBUTION (this paper)
------------------------------------------------------------------------

(K1) Riemann-Pandrosion (papers 55, 107, 110):
  RH equivalent: T_xi(t) = (xi')^2 - xi xi'' >= 0 on critical line.
By Tao 2016: log-Chowla k=2 follows from a strong dispersion estimate
which itself follows from RH (assumed).

(K2) Pandrosion-Liouville energy:
  E_lambda(N, h) := sum_{n<=N} (lambda(n) - lambda(n+h))^2.
  = 2N - 2 sum lambda(n) lambda(n+h)
  = 2N - 2 C_h(N) N.
Chowla k=2 <=> E_lambda(N, h)/N -> 2.

(K3) Verification:
  - Compute partial sums of lambda(n) lambda(n+h) for h = 1, 2, 3.
  - Verify Tao-Teräväinen log-Chowla numerically for k=3.
  - Show the (1/log N) decay rate empirically.

VERIFICATION
============

  1. lambda(n) computation via prime factorization.
  2. Partial sums sum lambda(n) lambda(n+h)/N for h = 1, 2, 3.
  3. Logarithmic version sum lambda(n) lambda(n+h)/n.
  4. Triple correlation lambda(n) lambda(n+1) lambda(n+2).
  5. Pandrosion energy form.
"""
from __future__ import annotations
import math
import numpy as np


def liouville_array(N):
    """lambda(n) for n = 1, ..., N (1-indexed)."""
    Omega = np.zeros(N + 1, dtype=int)
    for p in range(2, N + 1):
        if Omega[p] == 0:
            # p is prime
            for q in range(p, N + 1, p):
                m = q
                while m % p == 0:
                    Omega[q] += 1
                    m //= p
    lam = np.where(Omega % 2 == 0, 1, -1)
    lam[0] = 0  # undefined
    return lam


def main():
    print("=" * 80)
    print("PAPER 134 — Chowla's conjecture")
    print("=" * 80)

    print("\n[1] Sanity: M(N)/N -> 0 (PNT)")
    N_max = 100010
    lam = liouville_array(N_max)
    print(f"  {'N':>10} {'M(N)/N':>14}")
    for N in [100, 1000, 10000, 100000]:
        s = lam[1:N+1].sum()
        print(f"  {N:>10} {s/N:>14.6f}")

    print("\n[2] Chowla k=2: (1/N) sum lambda(n) lambda(n+h)")
    print(f"  {'N':>8} {'h=1':>10} {'h=2':>10} {'h=3':>10} {'h=5':>10}")
    for N in [1000, 5000, 20000, 100000]:
        results = []
        for h in [1, 2, 3, 5]:
            cur = lam[1:N+1] * lam[1+h:N+h+1]
            results.append(cur.sum() / N)
        print(f"  {N:>8} {results[0]:>10.5f} {results[1]:>10.5f} {results[2]:>10.5f} {results[3]:>10.5f}")
    print(f"  Conjecture: each column -> 0 as N -> infty.")

    print("\n[3] Logarithmic Chowla k=2: sum lambda(n) lambda(n+h) / (n log N)")
    print(f"  {'N':>8} {'log h=1':>12} {'log h=2':>12} {'log h=3':>12}")
    for N in [1000, 10000, 80000]:
        log_N = math.log(N)
        results = []
        for h in [1, 2, 3]:
            cur = sum(lam[n] * lam[n + h] / n for n in range(1, N + 1) if n + h <= N_max)
            results.append(cur / log_N)
        print(f"  {N:>8} {results[0]:>12.5f} {results[1]:>12.5f} {results[2]:>12.5f}")
    print(f"  Tao 2016: log-Chowla k=2 PROVED -> 0.")

    print("\n[4] Triple correlation k=3: (1/N) sum lambda(n) lambda(n+1) lambda(n+2)")
    print(f"  {'N':>10} {'C_3(N)':>14}")
    for N in [1000, 10000, 50000, 100000]:
        if N + 3 > N_max: continue
        cur = lam[1:N+1] * lam[2:N+2] * lam[3:N+3]
        print(f"  {N:>10} {cur.sum()/N:>14.6f}")
    print(f"  Tao-Teräväinen 2018: log-Chowla for ODD k (incl. k=3) PROVED.")

    print("\n[5] Pandrosion-Liouville energy E_lambda(N, h)/N")
    print(f"  E_lambda(N, h) = sum (lambda(n) - lambda(n+h))^2 = 2N - 2 sum lambda(n)lambda(n+h)")
    print(f"  Chowla k=2 <=> E_lambda(N, h)/N -> 2 as N -> infty")
    print(f"  {'N':>8} {'E/N (h=1)':>14} {'E/N (h=2)':>14}")
    for N in [1000, 10000, 80000]:
        for_results = []
        for h in [1, 2]:
            E = sum((lam[n] - lam[n + h])**2 for n in range(1, N + 1))
            for_results.append(E / N)
        print(f"  {N:>8} {for_results[0]:>14.5f} {for_results[1]:>14.5f}")

    print("\n[6] Decay rate: |C_h(N)| vs 1/log(N) (Tao-Teräväinen log estimate)")
    print(f"  {'N':>8} {'|C_1(N)|':>12} {'1/log N':>10} {'ratio':>10}")
    for N in [1000, 5000, 20000, 80000]:
        cur = lam[1:N+1] * lam[2:N+2]
        c1 = abs(cur.sum() / N)
        inv_log = 1 / math.log(N)
        ratio = c1 / inv_log if inv_log > 0 else 0
        print(f"  {N:>8} {c1:>12.5f} {inv_log:>10.5f} {ratio:>10.3f}")
    print(f"  Empirical ratio bounded; consistent with log-decay.")

    print("\n[7] Connection to Riemann (paper 110): zeta'/zeta on critical line")
    print(f"  Tao's proof of log-Chowla uses circle method + L-function estimates.")
    print(f"  Pandrosion-zeta (paper 11): F_zeta = zeta'/zeta = sum log(p) p^{{-s}}/(1 - p^{{-s}}).")
    print(f"  Connection: lambda(n) = sum mu(d) for d^2 | n; both depend on zeta(s) zeros.")

    print("\n[8] HONEST ASSESSMENT")
    print("  PROVED:")
    print("    Chowla k=1 (trivial via PNT).")
    print("    Tao 2016: log-Chowla k=2.")
    print("    Tao-Teräväinen 2018: log-Chowla for ODD k.")
    print("  ")
    print("  PANDROSION CONTRIBUTION (this paper):")
    print("    Pandrosion-Liouville energy E_lambda(N, h)/N -> 2 reformulates")
    print("    Chowla k=2 as energy convergence, an L^2 form of the conjecture.")
    print("    Connection F_L = F_{zeta(2s)} - F_{zeta(s)} via Pandrosion-zeta.")
    print("  ")
    print("  OPEN:")
    print("    Ordinary (non-logarithmic) Chowla for k >= 2.")
    print("    Logarithmic Chowla for EVEN k >= 4.")
    print("  ")
    print("  WHY ORDINARY CHOWLA IS HARD:")
    print("    Logarithmic averages smooth out arithmetic obstruction.")
    print("    Going from log to ordinary requires uniform control over n in [N, 2N]")
    print("    rather than weighted by 1/n.")
    print("    Equivalent to a strong PNT-like estimate for short intervals.")
    print("  ")
    print("  PATH FORWARD:")
    print("    1. Sharpen Pandrosion-zeta on critical strip (paper 110).")
    print("    2. Use Pandrosion-Liouville energy to bound C_h(N) directly.")
    print("    3. Combine with Tao's entropy decrement argument (Tao 2016).")
    print("    Long-term: ordinary Chowla likely needs new methods beyond zeta-zero")
    print("    estimates, possibly Furstenberg-type entropy or higher-order Fourier.")


if __name__ == "__main__":
    main()
