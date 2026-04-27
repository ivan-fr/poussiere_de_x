"""
PAPER: 027
TITLE: Pandrosionization of general polynomial systems
STATUS: experimental bridge from special Pandrosion systems to general systems

QUESTION
========

Can we keep the method Pandrosion-based for general polynomial systems,
instead of saying "025 only works because the system already has a Pandrosion
fixed-point structure"?

CORE IDEA
=========

Yes, if each equation can be put in coordinate-monic form

    z_i^p + R_i(z_1,...,z_n) = 0.

Then the Pandrosion map is

    H_i(z) = 1 - (R_i(z)+1) / S_p(z_i),
    S_p(u) = 1 + u + ... + u^{p-1}.

Indeed H_i(z)=z_i iff

    S_p(z_i)(1-z_i) = R_i(z)+1
    1 - z_i^p = R_i(z)+1
    z_i^p + R_i(z) = 0.

So every fixed point of H is a root of the original system, and conversely.
No derivative is used in H.  This is the exact multivariate analogue of the
"beyond monomials" construction in paper-0.

HOW GENERAL?
============

Univariate: every nonconstant polynomial can be normalized to monic, so this
is exact for all univariate polynomials.

Multivariate: it is exact when each equation has a usable pure pivot monomial
z_i^p with nonzero coefficient.  More general systems can often be moved into
this form by a coordinate change or by lifting to auxiliary variables
arithmetic-circuit style.  That lifting increases dimension, so the caveat
does not disappear: it moves from "not applicable" to "convergence and size
of the lift must be controlled".

This flow tests the coordinate-monic path and compares:

    Pandrosionized map + vector prism + Armijo fallback
    vs.
    finite-difference Newton + Armijo

The Newton baseline is not Lairez.  It is included only as a sanity check on
whether the Pandrosionized map is numerically viable beyond the hypercube
system of flow/025/026.
"""

from __future__ import annotations


# -----------------------------------------------------------------------------
# Sparse polynomial utilities.
# A polynomial is dict[exponent_tuple] = coefficient.
# -----------------------------------------------------------------------------
def eval_poly(poly, z):
    total = 0.0
    for exp, coeff in poly.items():
        term = coeff
        for zi, power in zip(z, exp):
            term *= zi ** power
        total += term
    return total


def poly_without_term(poly, term_exp):
    return {exp: coeff for exp, coeff in poly.items() if exp != term_exp and abs(coeff) > 0.0}


def scale_poly(poly, scale):
    return {exp: coeff * scale for exp, coeff in poly.items()}


def S_p(u, p):
    if p == 1:
        return 1.0
    return sum(u ** k for k in range(p))


def norm_inf(v):
    return max(abs(x) for x in v)


def residual(system, z):
    return [eval_poly(poly, z) for poly in system["polys"]]


def residual_norm(system, z):
    return norm_inf(residual(system, z))


# -----------------------------------------------------------------------------
# Pandrosionization.
# -----------------------------------------------------------------------------
def find_pivot(poly, n, pivot_index):
    """Find the highest pure power c*z_i^p present in one equation."""
    best = None
    for exp, coeff in poly.items():
        if abs(coeff) == 0.0:
            continue
        pure = all(power == 0 for j, power in enumerate(exp) if j != pivot_index)
        if pure and exp[pivot_index] >= 1:
            if best is None or exp[pivot_index] > best[0]:
                best = (exp[pivot_index], coeff, exp)
    return best


def pandrosionize(system):
    """
    Convert polynomials P_i(z)=0 to z_i^p + R_i(z)=0 by dividing each equation
    by the coefficient of a pure pivot monomial z_i^p.
    """
    polys = system["polys"]
    n = system["n"]
    data = []
    for i, poly in enumerate(polys):
        pivot = find_pivot(poly, n, i)
        if pivot is None:
            raise ValueError(f"equation {i} has no pure pivot monomial in variable {i}")
        p, coeff, exp = pivot
        normalized = scale_poly(poly, 1.0 / coeff)
        R = poly_without_term(normalized, exp)
        data.append({"p": p, "R": R})
    return {"n": n, "data": data, "source": system}


def pandrosion_map(z, pdata):
    out = []
    for i, item in enumerate(pdata["data"]):
        p = item["p"]
        Rz = eval_poly(item["R"], z)
        out.append(1.0 - (Rz + 1.0) / S_p(z[i], p))
    return out


# -----------------------------------------------------------------------------
# Vector prism with residual-aware Armijo fallback.
# -----------------------------------------------------------------------------
def diag_steffensen_step(T, z):
    z1 = T(z)
    z2 = T(z1)
    out = []
    for a, b, c in zip(z, z1, z2):
        d0 = b - a
        d1 = c - b
        d2 = d1 - d0
        out.append(c if abs(d2) < 1e-300 else a - d0 * d0 / d2)
    return out, 2


def damped_step(T, system, z, base_res):
    Hz = T(z)
    evals = 1
    lam = 1.0
    for _ in range(24):
        cand = [(1.0 - lam) * zi + lam * hi for zi, hi in zip(z, Hz)]
        r = residual_norm(system, cand)
        if r < base_res:
            return cand, evals, r
        lam *= 0.5
    return Hz, evals, residual_norm(system, Hz)


def solve_pandrosionized(system, z0, tol=1e-10, max_iter=400):
    pdata = pandrosionize(system)
    T = lambda z: pandrosion_map(z, pdata)
    z = list(z0)
    evals = 0
    for it in range(max_iter + 1):
        r0 = residual_norm(system, z)
        if r0 < tol:
            return z, it, evals, r0, True

        cand, used = diag_steffensen_step(T, z)
        evals += used
        rc = residual_norm(system, cand)
        if rc < r0:
            z = cand
            continue

        z, used, _ = damped_step(T, system, z, r0)
        evals += used
    return z, max_iter, evals, residual_norm(system, z), False


# -----------------------------------------------------------------------------
# Baseline: finite-difference Newton + Armijo.
# -----------------------------------------------------------------------------
def finite_difference_jacobian(system, z, h=1e-6):
    n = system["n"]
    f0 = residual(system, z)
    J = [[0.0] * n for _ in range(n)]
    for j in range(n):
        zp = list(z)
        zp[j] += h
        fp = residual(system, zp)
        for i in range(n):
            J[i][j] = (fp[i] - f0[i]) / h
    return J


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


def solve_fd_newton(system, z0, tol=1e-10, max_iter=40):
    z = list(z0)
    evals = 0
    n = system["n"]
    for it in range(max_iter + 1):
        r0 = residual_norm(system, z)
        evals += 1
        if r0 < tol:
            return z, it, evals, r0, True
        J = finite_difference_jacobian(system, z)
        evals += n + 1
        f = residual(system, z)
        step = solve_linear(J, [-v for v in f])
        if step is None:
            return z, it, evals, r0, False
        lam = 1.0
        accepted = False
        for _ in range(24):
            cand = [zi + lam * si for zi, si in zip(z, step)]
            rc = residual_norm(system, cand)
            evals += 1
            if rc < r0:
                z = cand
                accepted = True
                break
            lam *= 0.5
        if not accepted:
            return z, it, evals, r0, False
    return z, max_iter, evals, residual_norm(system, z), False


# -----------------------------------------------------------------------------
# Test systems.
# -----------------------------------------------------------------------------
def systems():
    return [
        {
            "name": "univariate cubic  z^3 - z - 1",
            "n": 1,
            "z0": [0.7],
            "polys": [
                {(3,): 1.0, (1,): -1.0, (0,): -1.0},
            ],
        },
        {
            "name": "coupled quadratic 2D",
            "n": 2,
            "z0": [0.8, 0.8],
            "polys": [
                {(2, 0): 1.0, (1, 1): 0.10, (0, 1): 0.20, (0, 0): -2.0},
                {(0, 2): 1.0, (1, 1): 0.05, (1, 0): -0.30, (0, 0): -1.5},
            ],
        },
        {
            "name": "general-looking 2D with linear pivot",
            "n": 2,
            "z0": [1.0, 1.0],
            "polys": [
                {(2, 0): 2.0, (1, 1): 1.0, (0, 2): 1.0, (0, 0): -3.0},
                {(0, 1): 1.0, (1, 0): -1.0, (0, 0): 0.25},
            ],
        },
        {
            "name": "coupled quadratic 3D",
            "n": 3,
            "z0": [0.9, 0.9, 0.9],
            "polys": [
                {(2, 0, 0): 1.0, (0, 1, 0): 0.30, (0, 1, 1): 0.05, (0, 0, 0): -1.7},
                {(0, 2, 0): 1.0, (1, 0, 0): -0.20, (1, 0, 1): 0.04, (0, 0, 0): -1.2},
                {(0, 0, 2): 1.0, (1, 0, 0): 0.25, (1, 1, 0): -0.03, (0, 0, 0): -1.6},
            ],
        },
    ]


def main():
    print("=" * 92)
    print("flow/027 -- Pandrosionization of general polynomial systems")
    print("=" * 92)
    print()
    print("Rule: rewrite each equation as z_i^p + R_i(z)=0, then use")
    print("      H_i(z)=1-(R_i(z)+1)/S_p(z_i).")
    print("Acceleration: vector difference-prism; fallback: damped Pandrosion step.")
    print()

    print(f"{'system':>38} | {'Pandrosionized prism':>28} | {'FD Newton':>24} | {'winner':>8}")
    print("-" * 108)
    for system in systems():
        try:
            pz, pit, pe, pr, pok = solve_pandrosionized(system, system["z0"])
        except ValueError as exc:
            pz, pit, pe, pr, pok = [], 0, 0, float("inf"), False
            print(f"  cannot pandrosionize {system['name']}: {exc}")

        nz, nit, ne, nr, nok = solve_fd_newton(system, system["z0"])
        ptag = f"{pit:>2}it/{pe:<3} eval r={pr:.1e} {'OK' if pok else 'NO'}"
        ntag = f"{nit:>2}it/{ne:<3} eval r={nr:.1e} {'OK' if nok else 'NO'}"
        if pok and (not nok or pe <= ne):
            winner = "pand"
        elif nok:
            winner = "newton"
        else:
            winner = "none"
        print(f"{system['name']:>38} | {ptag:>28} | {ntag:>24} | {winner:>8}")

    print()
    print("Verdict:")
    print("  Univariate: yes, every polynomial can be Pandrosionized after monic normalization.")
    print("  Multivariate: yes for coordinate-monic systems z_i^p+R_i=0.")
    print("  Fully general systems need a coordinate change or auxiliary-variable lift.")
    print("  This keeps the method Pandrosion-based, but the open problem becomes")
    print("  convergence/control of the lifted map, not fixed-point existence.")


if __name__ == "__main__":
    main()
