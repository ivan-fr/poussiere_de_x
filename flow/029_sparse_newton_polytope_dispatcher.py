"""
PAPER: 029
TITLE: Sparse/Newton-polytope Pandrosion dispatcher
STATUS: extends flow/028 with an exact-support geometry

MISSION
=======

flow/028 dispatches between:

    hypercube geometry  S_box     = prod_j S_p(s_j)
    simplex geometry    S_simplex = sum_{|alpha|<=p-1} s^alpha

flow/029 adds a third geometry:

    sparse geometry     S_A       = sum_{alpha in A} s^alpha

where A is the detected support of the geometric sum inside the system.

This is the Newton-polytope direction: if the support is not dense enough to
justify a box or a total-degree simplex, use the actual support.

SUPPORTED SYSTEM FORM
=====================

For a chosen support A, the Pandrosion system is

    P_i(s) = x_i (1-s_i) S_A(s) - (x_i - 1) = 0.

The fixed-point map is

    F_i(s) = 1 - (x_i - 1)/(x_i S_A(s)).

The dispatcher extracts A by looking for the shared positive support across
equations after the pattern

    +x_i s^alpha - x_i s^(alpha+e_i).

If the detected A is box-like or simplex-like, flow/028 would choose those
families.  If A is sparse, flow/029 chooses "sparse".
"""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_flow(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


HYPER = load_flow("flow025_hypercube", "025_symmetry_reduced_prism_vs_lairez.py")
SIMPLEX = load_flow("flow027_simplex", "027_symmetry_reduced_simplex_vs_lairez.py")


# -----------------------------------------------------------------------------
# Support utilities.
# -----------------------------------------------------------------------------
def clean(poly, eps=1e-12):
    return {exp: coeff for exp, coeff in poly.items() if abs(coeff) > eps}


def add_term(poly, exp, coeff):
    poly[exp] = poly.get(exp, 0.0) + coeff


def shift(exp, i):
    values = list(exp)
    values[i] += 1
    return tuple(values)


def unshift(exp, i):
    if exp[i] == 0:
        return None
    values = list(exp)
    values[i] -= 1
    return tuple(values)


def compositions_leq(n, max_total):
    out = []

    def rec(prefix, slots, remaining):
        if slots == 1:
            for k in range(remaining + 1):
                out.append(tuple(prefix + [k]))
            return
        for k in range(remaining + 1):
            rec(prefix + [k], slots - 1, remaining - k)

    rec([], n, max_total)
    return out


def box_support(n, p):
    out = []

    def rec(prefix, slots):
        if slots == 0:
            out.append(tuple(prefix))
            return
        for k in range(p):
            rec(prefix + [k], slots - 1)

    rec([], n)
    return out


def binom(n, k):
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    num = 1
    den = 1
    for i in range(1, k + 1):
        num *= n - k + i
        den *= i
    return num // den


def support(system):
    out = set()
    for poly in system["polys"]:
        for exp, coeff in poly.items():
            if abs(coeff) > 1e-12:
                out.add(exp)
    return out


def infer_x_vec(system):
    zero = tuple([0] * system["n"])
    x_vec = []
    for poly in system["polys"]:
        vals = [abs(coeff) for exp, coeff in poly.items() if exp != zero and abs(coeff) > 1e-12]
        x_vec.append(max(vals) if vals else 1.0)
    return x_vec


def detect_A_support(system, x_vec):
    """
    Detect A from negative -x_i terms of the form -x_i*s^(alpha+e_i).
    This is more reliable than reading positive terms, because the constant
    term x_i from x_i*S_A cancels with -(x_i-1), leaving coefficient 1.
    """
    detected = None
    for i, poly in enumerate(system["polys"]):
        xi = x_vec[i]
        Ai = set()
        for exp, coeff in poly.items():
            if abs(coeff + xi) > 1e-8:
                continue
            alpha = unshift(exp, i)
            if alpha is not None:
                Ai.add(alpha)
        if detected is None:
            detected = Ai
        else:
            detected = detected.intersection(Ai)
    return detected or set()


def is_box(A, n):
    if not A:
        return False
    max_exp = [max(exp[i] for exp in A) for i in range(n)]
    if n > 1 and any(m == 0 for m in max_exp):
        return False
    expected = set(box_support(n, max(max_exp) + 1))
    rectangular = set()
    def rec(prefix, slots):
        if slots == n:
            rectangular.add(tuple(prefix))
            return
        for k in range(max_exp[slots] + 1):
            rec(prefix + [k], slots + 1)
    rec([], 0)
    return set(A) == rectangular or set(A) == expected


def is_simplex(A, n):
    if not A:
        return False
    max_total = max(sum(exp) for exp in A)
    return set(A) == set(compositions_leq(n, max_total))


def geometry_scores(system, A):
    n = system["n"]
    if not A:
        return {"hypercube": 0.0, "simplex": 0.0, "sparse": 0.0}
    max_exp = [max(exp[i] for exp in A) for i in range(n)]
    box_volume = 1
    for m in max_exp:
        box_volume *= m + 1
    max_total = max(sum(exp) for exp in A)
    simplex_volume = binom(n + max_total, n)
    box_density = len(A) / box_volume
    simplex_density = len(A) / simplex_volume
    sparse_bonus = 1.0 - max(box_density, simplex_density)
    return {
        "hypercube": box_density,
        "simplex": simplex_density,
        "sparse": sparse_bonus,
        "box_density": box_density,
        "simplex_density": simplex_density,
        "A_size": len(A),
        "box_volume": box_volume,
        "simplex_volume": simplex_volume,
    }


def choose_geometry(system):
    x_vec = infer_x_vec(system)
    A = set(system["A"]) if "A" in system else detect_A_support(system, x_vec)
    if is_box(A, system["n"]):
        family = "hypercube"
    elif is_simplex(A, system["n"]):
        family = "simplex"
    else:
        family = "sparse"
    return family, A, x_vec, geometry_scores(system, A)


# -----------------------------------------------------------------------------
# Sparse Pandrosion solver.
# -----------------------------------------------------------------------------
def S_sparse(s, A):
    total = 0.0
    for alpha in A:
        term = 1.0
        for si, ai in zip(s, alpha):
            term *= si ** ai
        total += term
    return total


def dS_sparse(s, A, j):
    total = 0.0
    for alpha in A:
        if alpha[j] == 0:
            continue
        term = alpha[j] * (s[j] ** (alpha[j] - 1))
        for k, ak in enumerate(alpha):
            if k != j:
                term *= s[k] ** ak
        total += term
    return total


def F_target_sparse(s, x_vec, A):
    S = S_sparse(s, A)
    return [x_vec[i] * (1.0 - s[i]) * S - (x_vec[i] - 1.0) for i in range(len(s))]


def F_iter_sparse(s, args):
    x_vec, A = args
    S = S_sparse(s, A)
    return [1.0 - (x_vec[i] - 1.0) / (x_vec[i] * S) for i in range(len(s))]


def symmetry_groups(x_vec, tol=0.0):
    values = []
    groups = []
    for i, x in enumerate(x_vec):
        found = None
        for g, value in enumerate(values):
            if abs(x - value) <= tol:
                found = g
                break
        if found is None:
            values.append(x)
            groups.append([i])
        else:
            groups[found].append(i)
    return values, groups


def expand_from_quotient(y, groups, n):
    s = [0.0] * n
    for g, idxs in enumerate(groups):
        for i in idxs:
            s[i] = y[g]
    return s


def collapse_A_to_quotient(A, groups):
    collapsed = {}
    for alpha in A:
        beta = []
        for idxs in groups:
            beta.append(sum(alpha[i] for i in idxs))
        beta = tuple(beta)
        collapsed[beta] = collapsed.get(beta, 0) + 1
    return collapsed


def S_sparse_quotient(y, collapsed_A):
    total = 0.0
    for beta, coeff in collapsed_A.items():
        term = coeff
        for yi, bi in zip(y, beta):
            term *= yi ** bi
        total += term
    return total


def F_iter_sparse_quotient(y, args):
    x_unique, collapsed_A = args
    S = S_sparse_quotient(y, collapsed_A)
    if abs(S) < 1e-300:
        return list(y)
    return [1.0 - (x_unique[i] - 1.0) / (x_unique[i] * S) for i in range(len(y))]


def diag_steffensen(T, x0, args):
    x1 = T(x0, args)
    x2 = T(x1, args)
    out = []
    for a, b, c in zip(x0, x1, x2):
        d0 = b - a
        d1 = c - b
        d2 = d1 - d0
        out.append(c if abs(d2) < 1e-300 else a - d0 * d0 / d2)
    return out


def make_T_d(base, d):
    if d == 2:
        return base
    prev = make_T_d(base, d - 1)

    def T(x0, args):
        return diag_steffensen(prev, x0, args)
    return T


def prism_solve_sparse_quotient(x_vec, A, d=3, tol=1e-13, max_iter=200):
    x_unique, groups = symmetry_groups(x_vec)
    collapsed_A = collapse_A_to_quotient(A, groups)
    T = make_T_d(F_iter_sparse_quotient, d)
    y = [0.5] * len(x_unique)
    args = (x_unique, collapsed_A)
    cost = 2 ** (d - 2)
    evals = 0
    for it in range(max_iter):
        yn = T(y, args)
        evals += cost
        if max(abs(yn[i] - y[i]) for i in range(len(y))) < tol:
            s = expand_from_quotient(yn, groups, len(x_vec))
            return s, evals, it + 1, len(x_unique)
        y = yn
    return expand_from_quotient(y, groups, len(x_vec)), evals, max_iter, len(x_unique)


# -----------------------------------------------------------------------------
# Build systems.
# -----------------------------------------------------------------------------
def build_system_from_A(name, x_vec, A):
    n = len(x_vec)
    zero = tuple([0] * n)
    polys = []
    for i, x in enumerate(x_vec):
        poly = {}
        for alpha in A:
            add_term(poly, alpha, x)
            add_term(poly, shift(alpha, i), -x)
        add_term(poly, zero, -(x - 1.0))
        polys.append(clean(poly))
    return {"name": name, "n": n, "polys": polys, "A": set(A)}


def build_cases():
    return [
        build_system_from_A("box p2 n3", [2.0, 3.0, 5.0], box_support(3, 2)),
        build_system_from_A("simplex p2 n3", [2.0, 3.0, 5.0], compositions_leq(3, 1)),
        build_system_from_A("sparse cycle n3", [2.0, 3.0, 5.0], [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1)]),
        build_system_from_A("sparse diagonal n3", [2.0, 2.0, 5.0], [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]),
        build_system_from_A("thin chain n4", [2.0, 3.0, 5.0, 7.0], [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0), (1, 1, 1, 1)]),
    ]


def dispatch_solve(system, d=3):
    family, A, x_vec, scores = choose_geometry(system)
    if family == "hypercube":
        # Infer p from box support max exponent + 1.
        p = max(max(exp) for exp in A) + 1
        s, evals, iters, qdim, _ = HYPER.prism_solve_quotient(x_vec, p, d=d)
        res = max(abs(v) for v in HYPER.F_target(s, x_vec, p))
    elif family == "simplex":
        p = max(sum(exp) for exp in A) + 1
        s, evals, iters, qdim, _ = SIMPLEX.prism_solve_quotient(x_vec, p, d=d)
        res = max(abs(v) for v in SIMPLEX.F_target(s, x_vec, p))
    else:
        s, evals, iters, qdim = prism_solve_sparse_quotient(x_vec, A, d=d)
        res = max(abs(v) for v in F_target_sparse(s, x_vec, A))
        p = "-"
    return {
        "family": family,
        "p": p,
        "A_size": len(A),
        "x_vec": x_vec,
        "iters": iters,
        "evals": evals,
        "qdim": qdim,
        "res": res,
        "scores": scores,
    }


def main():
    print("=" * 104)
    print("flow/029 -- sparse/Newton-polytope Pandrosion dispatcher")
    print("=" * 104)
    print()
    print("Dispatcher families: hypercube | simplex | sparse exact-support S_A")
    print()

    print(
        f"{'system':>22} | {'chosen':>9} {'|A|':>4} {'p':>3} {'qdim':>4} "
        f"{'iters/evals':>12} {'res':>10} | {'box dens':>8} {'simp dens':>9}"
    )
    print("-" * 104)
    for system in build_cases():
        result = dispatch_solve(system, d=3)
        scores = result["scores"]
        print(
            f"{system['name']:>22} | "
            f"{result['family']:>9} {result['A_size']:>4} {str(result['p']):>3} "
            f"{result['qdim']:>4} {result['iters']:>3}it/{result['evals']:<5} "
            f"{result['res']:>10.1e} | "
            f"{scores.get('box_density', 0):>8.3f} {scores.get('simplex_density', 0):>9.3f}"
        )

    print()
    print("Verdict:")
    print("  The dispatcher no longer has to force everything into hypercube or")
    print("  simplex.  Sparse supports get their own Pandrosion sum S_A.")
    print("  This is the Newton-polytope direction: geometry follows support.")


if __name__ == "__main__":
    main()
