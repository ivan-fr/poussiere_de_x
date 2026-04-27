"""
PAPER: 014 (alternative geometry: cross-polytopes / octahedra)
TITLE: Pandrosion in d-cross-polytope (family gamma)
STATUS: empirical exploration

CONSTRUCTION
============
Cross-polytope (orthoplex) in dimension d: 2d vertices at +/- e_i,
faces are (d-1)-simplices on the orthants.
  d=2: square (4 vertices)
  d=3: octahedron (6 vertices, 8 faces)
  d=4: 16-cell (8 vertices, 16 faces)

Thales construction: insert hyperplanes parallel to one face triangle,
intersecting (d) edges in (d) intercept positions.  Each intercept is
parameterised by a separate variable s_i.

POLYNOMIAL DERIVATION (d=3 octahedron, p=2):
A face of the octahedron is the simplex {(x,y,z) : x+y+z=1, x,y,z>=0}.
The hyperplane parallel to this face at level h cuts 3 edges (toward
+e_1, +e_2, +e_3).  Three intercepts s_1, s_2, s_3 are coupled by the
single hyperplane equation.

Imposing geometric progression of slice volumes gives a system mixing
the 3 intercepts.  The natural polynomial:
  S_gamma,3(s_1, s_2, s_3) = sum_{cyclic} (1 + s_i + s_i*s_j)
where the sum is over cyclic triples (i, j, k).

SYMMETRY: cyclic but NOT permutation-invariant (oriented faces).
Should produce non-rank-1 collapse.
"""
from __future__ import annotations
import math


def S_octahedron(s_vec, p):
    """Cross-polytope-style polynomial: cyclic sum of geometric progressions
    along edge orientations.
    For p=2: S_gamma,d(s_1,...,s_n) = sum_i (1 + s_i + s_i*s_{i+1 mod n}).
    """
    s_vec = [complex(s) for s in s_vec]
    n = len(s_vec)
    total = 0.0+0.0j
    for i in range(n):
        s_i = s_vec[i]
        s_next = s_vec[(i+1) % n]
        # Contribution of edge i: 1 (start) + s_i + s_i * s_next (after junction)
        for k in range(p):
            total += s_i**k
        # Junction term: s_i * s_next
        total += s_i * s_next
    return total / n  # normalised by number of edges


def iterate_cross(x_vec, d, p, max_iter=200, tol=1e-13):
    n = d
    s = [complex(0.5)] * n
    history = [tuple(s)]
    for step in range(max_iter):
        Spsn = S_octahedron(s, p)
        s_new = [1 - (x_vec[i]-1)/(x_vec[i]*Spsn) for i in range(n)]
        if max(abs(s_new[i]-s[i]) for i in range(n)) < tol:
            return tuple(s_new), step+1, history
        s = s_new
        history.append(tuple(s))
    return tuple(s), max_iter, history


def estimate_lambda(history, s_star_real):
    errs = [abs(complex(h[0]).real - s_star_real) for h in history if h]
    ratios = [errs[i+1]/errs[i] for i in range(len(errs)-2) if errs[i] > 1e-15]
    return ratios[-1] if ratios else float('nan')


def main():
    print("="*70)
    print("Family gamma: Pandrosion in d-CROSS-POLYTOPE (orthoplex)")
    print("="*70)
    print()
    print("Polynomial: S_gamma,d(s_1, ..., s_d) = (1/d) * sum_i [S_p(s_i) + s_i*s_{i+1}]")
    print("Cyclic structure: not symmetric under arbitrary permutation.")
    print()
    print(f"  {'d':>3} {'x':>5} {'p':>3} {'s_*[0]':>14} {'lambda':>11} {'iter':>5}")
    print("-"*70)
    for d in [2, 3, 4, 5, 6]:
        for p in [2, 3]:
            x = 2.0
            n = d
            s_final, n_iter, history = iterate_cross((x,)*n, d, p)
            s_star_real = s_final[0].real
            lam = estimate_lambda(history, s_star_real)
            print(f"  {d:>3} {x:>5.1f} {p:>3} {s_star_real:>14.10f} {lam:>11.6f} {n_iter:>5}")
    print()
    # Check rank of Jacobian -- is rank-1 collapse breaking?
    print("Check: does iteration land on the diagonal in 1 step?")
    print("Setup d=3, p=2, off-diagonal start s_0 = (0.1, 0.5, 0.9):")
    s0 = [0.1+0j, 0.5+0j, 0.9+0j]
    Spsn = S_octahedron(s0, 2)
    s1 = [1 - (2-1)/(2*Spsn) for _ in range(3)]
    print(f"  s_1 = {[c.real for c in s1]}")
    print(f"  std(s_1) = {abs(s1[0]-s1[1]):.2e}")
    print()
    print("VERDICT: when all x_i are equal, iteration still goes to diagonal")
    print("(because H_i depends only on S which is the same scalar).")
    print("BUT: with asymmetric x_i (different ratios per direction), the")
    print("octahedral structure produces genuinely coupled non-diagonal dynamics.")
    print()
    print("Asymmetric octahedron test (x = (2, 3, 5), d=3, p=2):")
    s_final, n_iter, _ = iterate_cross((2.0, 3.0, 5.0), 3, 2)
    print(f"  fixed point: ({s_final[0].real:.6f}, {s_final[1].real:.6f}, {s_final[2].real:.6f})")
    print(f"  iterations: {n_iter}")
    Sp = S_octahedron(s_final, 2)
    for i, x in enumerate([2, 3, 5]):
        res = abs(x*Sp*(1-s_final[i]) - (x-1))
        print(f"    eq_{i+1} residual = {res:.2e}")


if __name__ == "__main__":
    main()
