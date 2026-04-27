"""
PAPER: 028
TITLE: Pandrosion geometry dispatcher
STATUS: first automatic selector for hypercube vs simplex geometry

MISSION
=======

Given a polynomial fixed-point system, decide which Pandrosion geometry should
be used:

    hypercube geometry  S(s)=prod_j S_p(s_j)
    simplex geometry    S(s)=sum_{|alpha|<=p-1} s^alpha

Then run the matching symmetry-reduced prism solver:

    hypercube -> flow/025 engine
    simplex   -> flow/027 engine

This is the algorithmic answer to:

    "on adapte la geometrie au systeme qu'on nous donne ?"

Yes.  The dispatcher looks at the monomial support / Newton-polytope shape.
If the support fills a box better than a total-degree simplex, use hypercube.
If it fills a total-degree simplex better than a box, use simplex.

This is not yet a universal solver.  It is the first layer of a Pandrosion
geometry dispatcher.  More geometries can be added later:

    chain, face-tiling, sparse toric/Newton polytope, lifted auxiliary systems.
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
# Sparse polynomial helpers.
# A system is {"name": str, "n": int, "polys": list[dict[exp_tuple,float]]}.
# -----------------------------------------------------------------------------
def clean(poly, eps=1e-12):
    return {exp: coeff for exp, coeff in poly.items() if abs(coeff) > eps}


def add_term(poly, exp, coeff):
    poly[exp] = poly.get(exp, 0.0) + coeff


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


def shift(exp, i):
    values = list(exp)
    values[i] += 1
    return tuple(values)


def support(system):
    supp = set()
    for poly in system["polys"]:
        for exp, coeff in poly.items():
            if abs(coeff) > 1e-12:
                supp.add(exp)
    return supp


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


def geometry_scores(system):
    supp = support(system)
    n = system["n"]
    if not supp:
        return {"hypercube": 0.0, "simplex": 0.0}

    max_exp = [max(exp[i] for exp in supp) for i in range(n)]
    box_volume = 1
    for m in max_exp:
        box_volume *= m + 1
    box_density = len(supp) / box_volume

    max_total = max(sum(exp) for exp in supp)
    simplex_volume = binom(n + max_total, n)
    simplex_density = len(supp) / simplex_volume

    # Boundary bonus: simplex systems tend to have all monomials under a
    # total-degree cap; hypercube systems tend to use a rectangular cap.
    balanced_box = min(max_exp) / max(max_exp) if max(max_exp) > 0 else 1.0
    hypercube_score = box_density * (0.75 + 0.25 * balanced_box)
    simplex_score = simplex_density

    return {
        "hypercube": hypercube_score,
        "simplex": simplex_score,
        "box_density": box_density,
        "simplex_density": simplex_density,
        "max_exp": max_exp,
        "max_total": max_total,
    }


def choose_geometry(system):
    scores = geometry_scores(system)
    family = "hypercube" if scores["hypercube"] >= scores["simplex"] else "simplex"
    return family, scores


def infer_p(system, family):
    supp = support(system)
    if family == "simplex":
        return max(sum(exp) for exp in supp)
    return max(max(exp) for exp in supp)


def infer_x_vec(system):
    """
    Infer x_i for Pandrosion systems P_i=x_i(1-s_i)S-(x_i-1).
    Nonconstant coefficients are typically +/-x_i, while the constant term
    cancels to 1.
    """
    x_vec = []
    zero = tuple([0] * system["n"])
    for poly in system["polys"]:
        vals = [abs(coeff) for exp, coeff in poly.items() if exp != zero and abs(coeff) > 1e-12]
        x_vec.append(max(vals) if vals else 1.0)
    return x_vec


# -----------------------------------------------------------------------------
# Build benchmark systems with known geometry.
# -----------------------------------------------------------------------------
def build_hypercube_system(x_vec, p):
    n = len(x_vec)
    S_supp = box_support(n, p)
    polys = []
    zero = tuple([0] * n)
    for i, x in enumerate(x_vec):
        poly = {}
        for exp in S_supp:
            add_term(poly, exp, x)
            add_term(poly, shift(exp, i), -x)
        add_term(poly, zero, -(x - 1.0))
        polys.append(clean(poly))
    return {"name": f"hypercube x={x_vec}, p={p}", "n": n, "polys": polys}


def build_simplex_system(x_vec, p):
    n = len(x_vec)
    S_supp = compositions_leq(n, p - 1)
    polys = []
    zero = tuple([0] * n)
    for i, x in enumerate(x_vec):
        poly = {}
        for exp in S_supp:
            add_term(poly, exp, x)
            add_term(poly, shift(exp, i), -x)
        add_term(poly, zero, -(x - 1.0))
        polys.append(clean(poly))
    return {"name": f"simplex x={x_vec}, p={p}", "n": n, "polys": polys}


def dispatch_solve(system, d=3):
    family, scores = choose_geometry(system)
    p = infer_p(system, family)
    x_vec = infer_x_vec(system)
    if family == "hypercube":
        s, evals, iters, qdim, counts = HYPER.prism_solve_quotient(x_vec, p, d=d)
        residual = max(abs(v) for v in HYPER.F_target(s, x_vec, p))
    else:
        s, evals, iters, qdim, counts = SIMPLEX.prism_solve_quotient(x_vec, p, d=d)
        residual = max(abs(v) for v in SIMPLEX.F_target(s, x_vec, p))
    return {
        "family": family,
        "p": p,
        "x_vec": x_vec,
        "iters": iters,
        "evals": evals,
        "qdim": qdim,
        "counts": counts,
        "residual": residual,
        "scores": scores,
    }


def main():
    print("=" * 100)
    print("flow/028 -- Pandrosion geometry dispatcher")
    print("=" * 100)
    print()
    print("Rule: score monomial support as box-like vs total-degree-simplex-like,")
    print("      choose hypercube or simplex, reduce symmetries, then prism-solve.")
    print()

    cases = [
        build_hypercube_system([2.0, 2.0, 2.0], 2),
        build_hypercube_system([2.0, 3.0, 5.0], 2),
        build_hypercube_system([2.0, 2.0, 5.0], 3),
        build_simplex_system([2.0, 2.0, 2.0], 2),
        build_simplex_system([2.0, 3.0, 5.0], 2),
        build_simplex_system([2.0, 2.0, 5.0], 3),
        build_simplex_system([2.0, 3.0, 5.0, 7.0], 2),
        build_hypercube_system([2.0, 3.0, 5.0, 7.0], 2),
    ]

    print(
        f"{'true system':>34} | {'chosen':>9} {'p':>2} {'qdim':>4} "
        f"{'iters/evals':>12} {'res':>10} | {'score H':>8} {'score S':>8}"
    )
    print("-" * 100)
    ok = 0
    for system in cases:
        result = dispatch_solve(system, d=3)
        true_family = "hypercube" if system["name"].startswith("hypercube") else "simplex"
        if result["family"] == true_family:
            ok += 1
        scores = result["scores"]
        print(
            f"{system['name'][:34]:>34} | "
            f"{result['family']:>9} {result['p']:>2} {result['qdim']:>4} "
            f"{result['iters']:>3}it/{result['evals']:<5} {result['residual']:>10.1e} | "
            f"{scores['hypercube']:>8.3f} {scores['simplex']:>8.3f}"
        )

    print()
    print(f"classification accuracy on benchmark: {ok}/{len(cases)}")
    print()
    print("Verdict:")
    print("  This is the first practical dispatcher: inspect support, choose")
    print("  hypercube or simplex geometry, then call the specialised solver.")
    print("  The next version should add sparse toric/Newton-polytope geometries")
    print("  instead of forcing every system into only these two families.")


if __name__ == "__main__":
    main()
