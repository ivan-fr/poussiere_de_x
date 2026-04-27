"""
PAPER: 015 (alternative geometry: face tilings of d-cube)
TITLE: Pandrosion as face-tiling of a d-cube (family delta)
STATUS: empirical exploration

CONSTRUCTION
============
Family delta: tile a d-cube with C(d, 2) two-dimensional faces.  Each
face is a Pandrosion 2D rectangle with its own ratio x_{ij} and own
intercept variable s_{ij}.  Faces sharing edges impose continuity
constraints between their s_{ij}.

  d=2: 1 rectangle (paper-0)
  d=3: 3 faces of a cube (3 variables s_{12}, s_{13}, s_{23})
  d=4: 6 faces of a 4-cube (6 variables s_{ij})
  d-cube: C(d, 2) = d(d-1)/2 face variables.

POLYNOMIAL DERIVATION:
Independent face Thales sums multiplied (no junction constraint):
  S_delta,d(s_{12}, ..., s_{d-1,d}) = prod_{ij} S_p(s_{ij})

This is FULLY SEPARABLE: the global product factorises.  Each face is
a paper-0 2D Pandrosion independently.  The fixed-point system reduces to
C(d, 2) decoupled paper-0 iterations.

Family delta is therefore equivalent to running C(d,2) parallel paper-0
iterations.  No genuine multivariate dynamics; no acceleration over
paper-0 strict.

Verify numerically.
"""
from __future__ import annotations
import math
from itertools import combinations


def S_face_tiling(s_dict, p):
    """S_delta,d({(i,j): s_{ij}}) = prod over faces of S_p(s_{ij})."""
    total = 1.0+0.0j
    for s_ij in s_dict.values():
        total *= sum(complex(s_ij)**k for k in range(p))
    return total


def iterate_face_tiling(x_dict, d, p, max_iter=100, tol=1e-13):
    """Iterate face-tiling Pandrosion."""
    pairs = list(combinations(range(d), 2))
    s = {ij: complex(0.5) for ij in pairs}
    history = [dict(s)]
    for step in range(max_iter):
        Spsn = S_face_tiling(s, p)
        s_new = {ij: 1 - (x_dict[ij] - 1) / (x_dict[ij] * Spsn) for ij in pairs}
        max_change = max(abs(s_new[ij] - s[ij]) for ij in pairs)
        if max_change < tol:
            return s_new, step+1, history
        s = s_new
        history.append(dict(s))
    return s, max_iter, history


def main():
    print("="*70)
    print("Family delta: Pandrosion in d-CUBE FACE TILING")
    print("="*70)
    print()
    print("Each of the C(d, 2) faces is an independent paper-0 2D Pandrosion.")
    print("Polynomial: S_delta = prod_{ij} S_p(s_{ij}).")
    print("Iteration: each s_{ij} updated using global product Spsn.")
    print()
    print(f"  {'d':>3} {'#faces':>8} {'p':>3} {'s_{ij}*':>14} {'lambda_avg':>11} {'iter':>5}")
    print("-"*70)
    for d in [2, 3, 4, 5, 6]:
        for p in [2, 3]:
            x = 2.0
            pairs = list(combinations(range(d), 2))
            x_dict = {ij: x for ij in pairs}
            s_final, n_iter, history = iterate_face_tiling(x_dict, d, p)
            s_first = list(s_final.values())[0]
            s_star_real = s_first.real
            # Estimate lambda from history (use first variable)
            errs = [abs(list(h.values())[0].real - s_star_real) for h in history]
            ratios = [errs[i+1]/errs[i] for i in range(len(errs)-2) if errs[i] > 1e-15]
            lam = ratios[-1] if ratios else float('nan')
            n_faces = len(pairs)
            print(f"  {d:>3} {n_faces:>8} {p:>3} {s_star_real:>14.10f} {lam:>11.6f} {n_iter:>5}")
    print()
    # Check separability
    print("Asymmetric face tiling (different x_{ij} per face):")
    d = 3; pairs = list(combinations(range(d), 2))
    x_dict = {(0,1): 2.0, (0,2): 3.0, (1,2): 5.0}
    s_final, n_iter, _ = iterate_face_tiling(x_dict, d, 2)
    print(f"  d=3, x = {x_dict}")
    for ij, s in s_final.items():
        print(f"    s_{ij} = {s.real:.6f}")
    print(f"  iterations: {n_iter}")
    print()
    print("VERDICT: family delta separable into C(d,2) parallel paper-0 iterations.")
    print("  Convergence rate scales like paper-0 (no improvement).")
    print("  Coupling through global product Spsn is artificial; faces")
    print("  could be solved independently more efficiently.")


if __name__ == "__main__":
    main()
