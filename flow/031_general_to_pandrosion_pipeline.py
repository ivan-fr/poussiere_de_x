"""
PAPER: 031
TITLE: General polynomial system -> Pandrosion geometry pipeline
STATUS: prototype

MISSION
=======

Implement the full pipeline:

    general polynomial system
      -> pandrosionisation / coordinate change / auxiliary lifting
      -> S_A Pandrosion form
      -> dispatcher hypercube/simplex/sparse

The target Pandrosion form is

    P_i(s) = x_i (1-s_i) S_A(s) - (x_i - 1).

If we can represent or lift a system into this form, flow/029 can dispatch the
geometry:

    hypercube / simplex / sparse exact support.

IMPORTANT
=========

This is a prototype, not a theorem.  It has three routes:

1. direct:
   The system is already in Pandrosion S_A form.

2. linear:
   Try small linear coordinate changes z=A y, then test if P(Ay) has
   Pandrosion S_A form.

3. auxiliary:
   Introduce auxiliary variables for mixed monomials, then test if the lifted
   system has Pandrosion S_A form.

If none of these routes succeeds, the system is not handled by this prototype.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_flow(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DISPATCH = load_flow("flow029_dispatcher", "029_sparse_newton_polytope_dispatcher.py")


# -----------------------------------------------------------------------------
# Sparse polynomial utilities.
# -----------------------------------------------------------------------------
def clean(poly, eps=1e-12):
    return {exp: coeff for exp, coeff in poly.items() if abs(coeff) > eps}


def add_term(poly, exp, coeff):
    poly[exp] = poly.get(exp, 0.0) + coeff


def eval_poly(poly, z):
    total = 0.0
    try:
        for exp, coeff in poly.items():
            term = coeff
            for zi, power in zip(z, exp):
                term *= zi ** power
            total += term
    except OverflowError:
        return float("inf")
    return total if math.isfinite(total) else float("inf")


def residual(system, z):
    return [eval_poly(poly, z) for poly in system["polys"]]


def residual_norm(system, z):
    values = residual(system, z)
    return max(abs(v) for v in values)


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


def degree(exp):
    return sum(exp)


def support_size(exp):
    return sum(1 for v in exp if v != 0)


def pure_exp(n, i, power):
    return tuple(power if j == i else 0 for j in range(n))


def extend_exp(exp, new_n):
    return tuple(list(exp) + [0] * (new_n - len(exp)))


# -----------------------------------------------------------------------------
# Linear composition P(Ay).
# -----------------------------------------------------------------------------
def poly_add(a, b):
    out = dict(a)
    for exp, coeff in b.items():
        out[exp] = out.get(exp, 0.0) + coeff
    return clean(out)


def poly_mul(a, b):
    if not a or not b:
        return {}
    n = len(next(iter(a)))
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            exp = tuple(ea[i] + eb[i] for i in range(n))
            out[exp] = out.get(exp, 0.0) + ca * cb
    return clean(out)


def poly_pow(poly, k, n):
    out = {tuple([0] * n): 1.0}
    for _ in range(k):
        out = poly_mul(out, poly)
    return out


def linear_form(row):
    n = len(row)
    return clean({tuple(1 if i == j else 0 for i in range(n)): row[j] for j in range(n)})


def compose_linear(poly, A):
    n = len(A)
    forms = [linear_form(A[j]) for j in range(n)]
    out = {}
    one = {tuple([0] * n): 1.0}
    for exp, coeff in poly.items():
        term = one
        for j, power in enumerate(exp):
            term = poly_mul(term, poly_pow(forms[j], power, n))
        term = {e: coeff * c for e, c in term.items()}
        out = poly_add(out, term)
    return clean(out)


def transform_system(system, A):
    z0 = system.get("z0", [0.5] * system["n"])
    return {
        "name": system["name"] + " | linear",
        "n": system["n"],
        "polys": [compose_linear(poly, A) for poly in system["polys"]],
        "z0": solve_linear(A, z0) or list(z0),
        "A_linear": A,
        "source": system,
    }


def identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def candidate_matrices(n):
    if n == 1:
        return [identity(1)]
    if n == 2:
        return [
            identity(2),
            [[1.0, 1.0], [1.0, -1.0]],
            [[1.0, 1.0], [1.0, 2.0]],
            [[1.0, 2.0], [1.0, -1.0]],
        ]
    if n == 3:
        return [
            identity(3),
            [[1.0, 1.0, 1.0], [1.0, -1.0, 0.0], [1.0, 0.0, -1.0]],
            [[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]],
        ]
    return [identity(n)]


def solve_linear(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-14:
            return None
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        rhs = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = rhs / M[i][i]
    return x


def mat_vec(A, y):
    return [sum(A[i][j] * y[j] for j in range(len(y))) for i in range(len(A))]


# -----------------------------------------------------------------------------
# Auxiliary lifting.
# -----------------------------------------------------------------------------
def collect_mixed_monomials(system):
    mixed = []
    seen = set()
    for poly in system["polys"]:
        for exp in poly:
            if degree(exp) >= 2 and support_size(exp) >= 2 and exp not in seen:
                seen.add(exp)
                mixed.append(exp)
    return mixed


def monomial_value(exp, z):
    value = 1.0
    for zi, power in zip(z, exp):
        value *= zi ** power
    return value


def lift_system(system):
    n0 = system["n"]
    mixed = collect_mixed_monomials(system)
    aux_index = {exp: n0 + k for k, exp in enumerate(mixed)}
    n = n0 + len(mixed)

    lifted_polys = []
    for poly in system["polys"]:
        out = {}
        for exp, coeff in poly.items():
            if exp in aux_index:
                new_exp = pure_exp(n, aux_index[exp], 1)
            else:
                new_exp = extend_exp(exp, n)
            add_term(out, new_exp, coeff)
        lifted_polys.append(clean(out))

    for exp, aux in aux_index.items():
        aux_poly = {pure_exp(n, aux, 1): 1.0}
        add_term(aux_poly, extend_exp(exp, n), -1.0)
        lifted_polys.append(clean(aux_poly))

    z0 = list(system.get("z0", [0.5] * n0))
    base_z0 = system.get("z0", [0.5] * n0)
    for exp in mixed:
        z0.append(monomial_value(exp, base_z0))

    return {
        "name": system["name"] + " | lift",
        "n": n,
        "n_original": n0,
        "polys": lifted_polys,
        "z0": z0,
        "mixed": mixed,
        "source": system,
    }


# -----------------------------------------------------------------------------
# Detect exact Pandrosion S_A form.
# -----------------------------------------------------------------------------
def infer_x_from_poly(poly, n):
    zero = tuple([0] * n)
    vals = [abs(c) for exp, c in poly.items() if exp != zero and abs(c) > 1e-12]
    return max(vals) if vals else 1.0


def detect_pandrosion_form(system):
    n = system["n"]
    x_vec = [infer_x_from_poly(poly, n) for poly in system["polys"]]
    if "A" in system:
        A = set(system["A"])
    else:
        A_common = None
        for i, poly in enumerate(system["polys"]):
            xi = x_vec[i]
            Ai = set()
            for exp, coeff in poly.items():
                if abs(coeff + xi) < 1e-8:
                    alpha = unshift(exp, i)
                    if alpha is not None:
                        Ai.add(alpha)
            if A_common is None:
                A_common = Ai
            else:
                A_common = A_common.intersection(Ai)
        A = A_common or set()
    if not A:
        return None

    # Rebuild the expected Pandrosion system and compare exactly enough.
    expected = DISPATCH.build_system_from_A("expected", x_vec, sorted(A))
    for poly, exp_poly in zip(system["polys"], expected["polys"]):
        keys = set(poly) | set(exp_poly)
        for exp in keys:
            if abs(poly.get(exp, 0.0) - exp_poly.get(exp, 0.0)) > 1e-8:
                return None
    return {"x_vec": x_vec, "A": sorted(A)}


def solve_pandrosion_form(system, form, d=3):
    structured = DISPATCH.build_system_from_A(system["name"] + " | S_A", form["x_vec"], form["A"])
    result = DISPATCH.dispatch_solve(structured, d=d)
    return result, structured


def result_residual(result):
    return result.get("residual", result.get("res", float("inf")))


def attempt_pipeline(system):
    attempts = []

    direct = detect_pandrosion_form(system)
    if direct is not None:
        result, structured = solve_pandrosion_form(system, direct, d=3)
        attempts.append(("direct", result, structured, None))

    for A in candidate_matrices(system["n"]):
        transformed = transform_system(system, A)
        form = detect_pandrosion_form(transformed)
        if form is not None:
            result, structured = solve_pandrosion_form(transformed, form, d=3)
            attempts.append(("linear", result, structured, A))

    lifted = lift_system(system)
    form = detect_pandrosion_form(lifted)
    if form is not None:
        result, structured = solve_pandrosion_form(lifted, form, d=3)
        attempts.append(("lift", result, structured, None))

    if not attempts:
        return None
    return min(attempts, key=lambda item: (result_residual(item[1]), item[1]["evals"]))


# -----------------------------------------------------------------------------
# Test systems.
# -----------------------------------------------------------------------------
def systems():
    return [
        DISPATCH.build_system_from_A("already hypercube", [2.0, 3.0, 5.0], DISPATCH.box_support(3, 2)),
        DISPATCH.build_system_from_A("already simplex", [2.0, 3.0, 5.0], DISPATCH.compositions_leq(3, 1)),
        DISPATCH.build_system_from_A("already sparse", [2.0, 3.0, 5.0], [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)]),
        {
            "name": "rotated 2D Pandrosion",
            "n": 2,
            "z0": [1.0, 0.2],
            "polys": [
                # This is intentionally not exact Pandrosion in z; current
                # prototype may or may not recover it through candidate A.
                {(2, 0): 1.0, (1, 1): 2.0, (0, 2): 1.0, (1, 0): -2.0, (0, 1): -2.0, (0, 0): 1.0},
                {(2, 0): 1.0, (1, 1): -2.0, (0, 2): 1.0, (1, 0): -3.0, (0, 1): 3.0, (0, 0): 1.0},
            ],
        },
        {
            "name": "generic not handled",
            "n": 2,
            "z0": [0.8, 0.8],
            "polys": [
                {(2, 0): 1.0, (1, 1): 0.37, (0, 1): -0.2, (0, 0): -1.0},
                {(0, 2): 1.0, (1, 0): 0.11, (1, 1): -0.41, (0, 0): -0.7},
            ],
        },
    ]


def main():
    print("=" * 104)
    print("flow/031 -- general system -> Pandrosion geometry pipeline")
    print("=" * 104)
    print()
    print("Routes: direct S_A detection | linear coordinate change | auxiliary lift")
    print()
    print(f"{'system':>24} | {'route':>8} {'geom':>9} {'|A|':>4} {'qdim':>4} {'cost':>10} {'res':>10}")
    print("-" * 104)
    for system in systems():
        attempt = attempt_pipeline(system)
        if attempt is None:
            print(f"{system['name']:>24} | {'NONE':>8} {'-':>9} {'-':>4} {'-':>4} {'-':>10} {'not handled':>10}")
            continue
        route, result, structured, A = attempt
        print(
            f"{system['name']:>24} | {route:>8} {result['family']:>9} "
            f"{result['A_size']:>4} {result['qdim']:>4} "
            f"{result['iters']:>3}it/{result['evals']:<4} {result_residual(result):>10.1e}"
        )

    print()
    print("Verdict:")
    print("  031 is the first full pipeline skeleton.  It successfully routes")
    print("  systems already in Pandrosion S_A form into the dispatcher.")
    print("  The hard part is still automatic pandrosionisation of arbitrary")
    print("  coefficients; when detection fails, the system is honestly marked")
    print("  not handled instead of pretending universality.")


if __name__ == "__main__":
    main()
