"""
Rigorous audit of mason_stothers_landscape.py.
==============================================

Checks:
  (A) Sympy-symbolic verification of slack = 0 on the line y = -x,
      x != 0, for the family a = t^6, c = (t-x)^3 (t-y)^3, b = c - a.

  (B) Numerical scan for points (x, y) NOT on y = -x with slack = 0
      (i.e., other extremal loci), to confirm or refute that y = -x
      is the entire real-locus extremal.

  (C) Hand-check of slack and quality at 6+ specific (x, y) values.

  (D) Identification of NaN-handling artifacts in the figures.

  (E) Verification of the log|disc(b)| computation via independent
      sympy evaluation.
"""

from __future__ import annotations
import numpy as np
import sympy as sp
from mason_stothers_landscape import slack_3_3, quality, log_disc_b


# =====================================================================
# (A) Symbolic verification: y = -x extremal
# =====================================================================

def symbolic_extremal_check():
    print("=" * 72)
    print("(A) SYMBOLIC: slack = 0 on y = -x for the family a=t^6, c=(t-x)^3(t+x)^3")
    print("=" * 72)
    t, x = sp.symbols('t x', real=True, nonzero=True)
    a = t**6
    c = (t - x)**3 * (t + x)**3
    b = sp.expand(c - a)
    print(f"  c(t) = (t-x)^3 (t+x)^3 = {sp.expand(c)}")
    print(f"  b(t) = c - t^6        = {b}")
    print(f"  deg_t(b) = {sp.Poly(b, t).degree()}")
    # Roots of b in t:
    b_poly = sp.Poly(b, t)
    print(f"  b as poly in t (coefficients high->low): {b_poly.all_coeffs()}")
    # Reduce: pull out -x^2
    b_simplified = sp.expand(b / (-x**2))
    print(f"  b / (-x^2) = {b_simplified}")
    # Treat as polynomial in u = t^2
    u = sp.symbols('u', real=True)
    b_in_u = b_simplified.subs(t**2, u).expand()
    print(f"  with u = t^2:  3u^2 - 3 x^2 u + x^4 = {b_in_u}")
    # Discriminant of 3u^2 - 3x^2 u + x^4
    disc_u = (3 * x**2)**2 - 4 * 3 * x**4
    print(f"  discriminant in u: {sp.simplify(disc_u)}  (= -3 x^4, always negative for x real ≠ 0)")
    print(f"  ⟹ two complex-conjugate values of u, four distinct t-roots.")
    print()
    # Distinct roots count:
    # a contributes {0}: 1 distinct
    # c contributes {x, -x}: 2 distinct
    # b contributes 4 distinct roots, none of which are 0, x, or -x
    # (because b(0) = -x^6, b(x) = -x^6, b(-x) = -x^6, all nonzero for x != 0)
    print("  b(0) = c(0) - 0 =", sp.expand(c.subs(t, 0)))   # = x^3 (-x)^3 = -x^6
    # Actually c(0) = (-x)^3 (x)^3 = -x^3 · x^3 = -x^6. So b(0) = -x^6.
    print("  b(x) =", sp.expand(b.subs(t, x)))
    print("  b(-x) =", sp.expand(b.subs(t, -x)))
    print()
    print("  ⟹ Distinct roots of abc: {0, x, -x} ∪ {4 complex roots of b} = 7 distinct.")
    print("  ⟹ deg rad(abc) = 7;  max deg = 6;  slack = 7 - 1 - 6 = 0.  EXTREMAL ✓")
    print()


# =====================================================================
# (B) Numerical scan for OTHER extremals (slack = 0 outside y = -x)
# =====================================================================

def scan_other_extremals():
    print("=" * 72)
    print("(B) NUMERICAL SCAN for slack = 0 OUTSIDE the line y = -x")
    print("=" * 72)
    found = []
    N = 300
    span = 3.0
    for x in np.linspace(-span, span, N):
        for y in np.linspace(-span, span, N):
            if abs(x) < 0.1 or abs(y) < 0.1 or abs(x - y) < 0.1:
                continue
            if abs(x + y) < 0.05:
                continue   # exclude y = -x line
            s = slack_3_3(x, y)
            if np.isfinite(s) and s == 0:
                found.append((x, y, s))
    print(f"  Scanned grid: {N}x{N}, span [-{span}, {span}], excluding the y = -x band.")
    print(f"  Points with slack = 0 OUTSIDE y = -x: {len(found)}")
    if found:
        for x, y, s in found[:10]:
            print(f"    (x={x:+.4f}, y={y:+.4f}) slack={s}, x+y={x+y:.6f}")
    else:
        print("  ⟹ No other extremals found on this grid. y = -x appears to be")
        print("     the entire real extremal locus for this family.")
    # Explicit small-margin search near specific candidates
    candidates = [(1.0, 1.0), (1.5, 0.5), (2.0, 0.5), (-1.0, 0.5)]
    print()
    print(f"  Quick sanity check at non-extremal points:")
    for x, y in candidates:
        s = slack_3_3(x, y)
        print(f"    ({x}, {y}): slack = {s}")


# =====================================================================
# (C) Hand-check of specific values
# =====================================================================

def hand_check_specific():
    print()
    print("=" * 72)
    print("(C) HAND-CHECK of specific (x, y) values")
    print("=" * 72)
    cases = [
        ("(1, -1) extremal expected",     1.0, -1.0,  0,  1.0000),
        ("(2, -2) extremal expected",     2.0, -2.0,  0,  1.0000),
        ("(0.7, -0.7) extremal expected", 0.7, -0.7,  0,  1.0000),
        ("(1, 2) generic expected",       1.0,  2.0,  1,  6/7),
        ("(2, 1) generic expected",       2.0,  1.0,  1,  6/7),
        ("(-2, 1) generic expected",     -2.0,  1.0,  1,  6/7),
        ("(1.5, -0.5) generic expected",  1.5, -0.5,  1,  6/7),
    ]
    print(f"  {'description':<35s}  {'slack':>6s}  {'q':>9s}  {'pass'}")
    for label, x, y, exp_s, exp_q in cases:
        s = slack_3_3(x, y)
        q = quality(x, y)
        ok = (s == exp_s) and (abs(q - exp_q) < 1e-3)
        flag = "OK" if ok else "FAIL"
        print(f"  {label:<35s}  {s:>6.0f}  {q:>9.4f}  {flag}")


# =====================================================================
# (D) NaN-handling audit in the figures
# =====================================================================

def nan_audit():
    print()
    print("=" * 72)
    print("(D) NaN-HANDLING AUDIT in the three figures")
    print("=" * 72)
    print()
    # Recompute the slack and quality grids
    N = 110
    span = 3.0
    xs = np.linspace(-span, span, N)
    ys = np.linspace(-span, span, N)
    S = np.zeros((N, N))
    Q = np.zeros((N, N))
    LD = np.zeros((N, N))
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            S[j, i] = slack_3_3(x, y)
            Q[j, i] = quality(x, y)
            LD[j, i] = log_disc_b(x, y)
    n_total = N * N
    n_nan_S = int(np.sum(np.isnan(S)))
    n_nan_Q = int(np.sum(np.isnan(Q)))
    n_nan_LD = int(np.sum(np.isnan(LD)))
    print(f"  Grid size: {N}x{N} = {n_total} cells")
    print(f"  NaN cells in S (slack):     {n_nan_S} ({100*n_nan_S/n_total:.2f}%)")
    print(f"  NaN cells in Q (quality):   {n_nan_Q} ({100*n_nan_Q/n_total:.2f}%)")
    print(f"  NaN cells in LD (log disc): {n_nan_LD} ({100*n_nan_LD/n_total:.2f}%)")
    print()
    print("  Figure 1 (slack):  NaN ↦ nanmax(S) = max value seen.")
    print("    This treats degenerate configs as 'maximally non-extremal'.")
    print("    Effect: bands along x=0, y=0, x=y appear at the MAX height.")
    print("    Visual: defensive but doesn't pretend they are 'low quality'. OK.")
    print()
    print("  Figure 2 (quality):  NaN ↦ 0.")
    print("    BUG: cells along the degeneracy bands (x=0, y=0, x=y) are")
    print("    plotted as q=0 (deep dark colour), giving the illusion of")
    print("    'genuine quality 0' regions.  These are actually undefined.")
    print("    Should be masked or set to NaN-aware plotting.")
    print()
    print("  Figure 3 (log|disc|):  NaN ↦ 0, then clipped.")
    print("    Similar issue to Figure 2.  The 'central depression' visible")
    print("    in the rendered figure is partly an artefact of NaN→0 along")
    print("    the degeneracy bands.")
    print()
    distinct_S = sorted(set(S[~np.isnan(S)].astype(int).flatten().tolist()))
    distinct_Q_set = sorted(set(np.round(Q[~np.isnan(Q)], 4).flatten().tolist()))[:10]
    print(f"  Distinct slack values seen: {distinct_S}")
    print(f"  Distinct quality values (first 10): {distinct_Q_set}")


# =====================================================================
# (E) Independent sympy evaluation of log|disc(b)|
# =====================================================================

def disc_audit():
    print()
    print("=" * 72)
    print("(E) Independent sympy verification of log|disc(b)|")
    print("=" * 72)
    t, x_s, y_s = sp.symbols('t x y', real=True)
    a = t**6
    c = (t - x_s)**3 * (t - y_s)**3
    b = sp.expand(c - a)
    b_poly = sp.Poly(b, t)
    disc_sym = b_poly.discriminant()
    print(f"  disc(b) (symbolic, simplified) =")
    disc_simplified = sp.factor(disc_sym)
    print(f"    {disc_simplified}")
    test_pts = [(1.0, -1.0), (2.0, -2.0), (1.0, 2.0), (1.5, -0.5)]
    print()
    print(f"  {'(x, y)':>20s}  {'sympy disc(b)':>30s}  {'log10|disc|':>15s}  "
          f"{'numpy log_disc_b':>20s}")
    for x_v, y_v in test_pts:
        d_sym = float(disc_sym.subs({x_s: x_v, y_s: y_v}))
        ld_sym = np.log10(abs(d_sym) + 1e-15)
        ld_num = log_disc_b(x_v, y_v)
        match = "OK" if abs(ld_sym - ld_num) < 1e-6 else "MISMATCH"
        print(f"  ({x_v:+.2f}, {y_v:+.2f})  {d_sym:>+30.6e}  "
              f"{ld_sym:>15.4f}  {ld_num:>20.4f}  [{match}]")


def main():
    symbolic_extremal_check()
    scan_other_extremals()
    hand_check_specific()
    nan_audit()
    disc_audit()


if __name__ == "__main__":
    main()
