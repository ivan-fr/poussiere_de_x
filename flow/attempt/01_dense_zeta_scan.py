"""
ATTEMPT 01 — Dense empirical scan of |ζ(1+it)| / (log|t|)^{2/3}.

Goal: produce a much tighter EMPIRICAL upper bound on
    C_emp(t) := |ζ(1+it)| / (log|t|)^{2/3}
than what paper 154 reports (max ~ 1 over coarse log samples up to 10^7).

This is NOT rigorous — but it sets the realistic target for any
rigorous reduction of the V-K constant C_VK = 76.2 (MTY 2022).

Strategy
--------
- High mpmath precision (mp.dps = 30).
- Dense uniform sampling over t ∈ [10, 10^4] (where |ζ(1+it)| has its
  largest oscillations relative to the smooth growth).
- Coarser log-spaced sampling over t ∈ [10^4, 10^9].
- Track sup, near-sup, and the ratio sup / 76.2.
"""
from __future__ import annotations
import math
import time

from mpmath import mp, mpc, mpf, log as mlog, zeta as mzeta

mp.dps = 30

MTY_C = 76.2


def C_emp(t: float) -> float:
    z = abs(mzeta(mpc(1, t)))
    return float(z / mlog(t) ** (mpf("2") / 3))


def dense_scan(t_lo: float, t_hi: float, n: int):
    step = (t_hi - t_lo) / n
    sup = 0.0
    sup_t = t_lo
    top = []  # (ratio, t)
    for k in range(n + 1):
        t = t_lo + k * step
        if t < 3:
            continue
        r = C_emp(t)
        if r > sup:
            sup = r
            sup_t = t
        top.append((r, t))
    top.sort(reverse=True)
    return sup, sup_t, top[:5]


def log_scan(t_lo: float, t_hi: float, n: int):
    a, b = math.log10(t_lo), math.log10(t_hi)
    step = (b - a) / n
    sup = 0.0
    sup_t = t_lo
    top = []
    for k in range(n + 1):
        t = 10 ** (a + k * step)
        r = C_emp(t)
        if r > sup:
            sup = r
            sup_t = t
        top.append((r, t))
    top.sort(reverse=True)
    return sup, sup_t, top[:5]


def main():
    print("=" * 80)
    print("ATTEMPT 01 — Dense empirical scan of C_emp(t) = |ζ(1+it)|/(log t)^{2/3}")
    print("=" * 80)

    bands = [
        ("dense uniform t in [10, 100]",       10.0,    100.0, 2000, "uniform"),
        ("dense uniform t in [100, 1000]",     100.0,   1000.0, 4000, "uniform"),
        ("dense uniform t in [1000, 10000]",   1000.0,  10000.0, 5000, "uniform"),
        ("log-spaced t in [10^4, 10^6]",       1e4,     1e6, 1500, "log"),
        ("log-spaced t in [10^6, 10^9]",       1e6,     1e9, 1500, "log"),
    ]

    overall_sup = 0.0
    overall_t = 0.0

    for label, lo, hi, n, kind in bands:
        t0 = time.time()
        if kind == "uniform":
            sup, sup_t, top = dense_scan(lo, hi, n)
        else:
            sup, sup_t, top = log_scan(lo, hi, n)
        dt = time.time() - t0
        print(f"\n[{label}]  n={n}  ({dt:.1f}s)")
        print(f"  sup C_emp = {sup:.6f}  at t ≈ {sup_t:.4f}")
        print(f"  top 5: " + ", ".join(f"({r:.4f}@{t:.2f})" for r, t in top))
        if sup > overall_sup:
            overall_sup = sup
            overall_t = sup_t

    print("\n" + "=" * 80)
    print(f"  OVERALL EMPIRICAL SUP   :  {overall_sup:.6f}  at t ≈ {overall_t:.4f}")
    print(f"  MTY 2022 rigorous C_VK  :  {MTY_C}")
    print(f"  rigorous / empirical    :  {MTY_C / overall_sup:.1f} ×")
    print()
    print("  -> empirical C is well below 1, while MTY's rigorous C is 76.2.")
    print("     Rigorously squeezing C below 76.2 would already beat MTY 2022.")
    print("     This file PROVES NOTHING — it only sets the empirical target.")

    # ----------------------------------------------------------------
    # Heuristic projection of R given a hypothetical rigorous C
    # ----------------------------------------------------------------
    V0 = math.log(2 * math.pi) - 0.5772156649  # log 2π - γ_Euler
    alpha = (1.5587 - V0) / math.log(MTY_C)    # calibrate to MTY 2022
    print()
    print("[heuristic R projection]   R = 4 + V0 + α log C,"
          f" V0={V0:.4f}, α={alpha:.4f}")
    print(f"  {'C':>10} {'V':>10} {'R':>10} {'gain vs MTY':>14}")
    for C in [76.2, 70.0, 50.0, 30.0, 10.0, 5.0, 2.0, 1.0, overall_sup]:
        V = V0 + alpha * math.log(C)
        R = 4 + V
        print(f"  {C:>10.4f} {V:>10.4f} {R:>10.4f} {5.5587 - R:>14.4f}")


if __name__ == "__main__":
    main()
