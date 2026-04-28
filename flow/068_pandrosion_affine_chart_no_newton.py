"""
FLOW 068 -- dD polynomial Pandrosion with affine-chart recovery, no Newton-ELS.

Core identity
-------------
All local geometries keep the polynomial Pandrosion identity from
1pandrosion_smale.tex:

    F(z) - F(a) = Q_G(a,z) (z-a)
    P_G(a,z)   = a - Q_G(a,z)^(-1) F(a)

so every regular zero of F is a fixed point of P_G.  The main slope is

    Q_int,ij(a,z) = int_0^1 d_j F_i(a+t(z-a)) dt,

which reduces in one variable to x*S_p(s) for F(s)=x*s^p-1 and a=1.
The portfolio also includes exactified dD simplex, sparse, support,
hypercube-covector, path and cube geometries.

What is new vs 067
------------------
067 could replace Newton-ELS locally, but the hard (4,2) straight homotopy case
still coalesced at 15/16 roots.  068 adds affine-chart recovery: solve
F(Ay+b)=0 under deterministic near-identity dense charts and map back z=Ay+b.
This changes the path geometry without using Lairez gamma rotation and without
using Newton-ELS.  The random test generator is now the same as flow/064 for
fair A/B comparison.
"""
from __future__ import annotations

import cmath
import itertools
import math
import time
from functools import lru_cache
from itertools import product as iprod
from math import comb
from typing import Dict, Iterable, List, Sequence, Tuple

Complex = complex
Exponent = Tuple[int, ...]
Poly = Dict[Exponent, Complex]
System = List[Poly]
Vector = List[Complex]
Matrix = List[List[Complex]]

MAX_Z = 1.0e7
MAX_TERM = 1.0e220


def finite(v: Complex) -> bool:
    v = complex(v)
    return math.isfinite(v.real) and math.isfinite(v.imag)


def finite_vec(z: Sequence[Complex]) -> bool:
    return all(finite(v) for v in z)


def safe_z(z: Sequence[Complex], bound: float = MAX_Z) -> bool:
    return finite_vec(z) and all(abs(v) <= bound for v in z)


def degree(poly: Poly) -> int:
    return max((sum(e) for e in poly), default=0)


def system_degree(polys: System) -> int:
    return max((degree(p) for p in polys), default=1)


def eval_poly(poly: Poly, z: Sequence[Complex]) -> Complex:
    if not safe_z(z):
        return complex("inf")
    total = 0.0 + 0.0j
    try:
        for exp, coeff in poly.items():
            term = complex(coeff)
            for zj, aj in zip(z, exp):
                if aj:
                    term *= zj ** aj
                    if not finite(term) or abs(term) > MAX_TERM:
                        return complex("inf")
            total += term
            if not finite(total) or abs(total) > MAX_TERM:
                return complex("inf")
    except (OverflowError, ZeroDivisionError):
        return complex("inf")
    return total if finite(total) else complex("inf")


def F_eval(polys: System, z: Sequence[Complex]) -> Vector:
    return [eval_poly(p, z) for p in polys]


def residual_norm(polys: System, z: Sequence[Complex]) -> float:
    if not safe_z(z):
        return float("inf")
    vals = F_eval(polys, z)
    if not finite_vec(vals):
        return float("inf")
    return max(abs(v) for v in vals)


def residual_norm2(polys: System, z: Sequence[Complex]) -> float:
    vals = F_eval(polys, z)
    if not finite_vec(vals):
        return float("inf")
    return sum(abs(v) ** 2 for v in vals)


def solve_linear(A: Matrix, b: Sequence[Complex], tol: float = 1e-13) -> Vector | None:
    n = len(A)
    M = [list(row) + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        piv = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[piv] = M[piv], M[k]
        if abs(M[k][k]) < tol:
            return None
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            if f == 0:
                continue
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0.0 + 0.0j] * n
    for i in range(n - 1, -1, -1):
        rhs = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        if abs(M[i][i]) < tol:
            return None
        x[i] = rhs / M[i][i]
    return x if safe_z(x) else None


# -----------------------------------------------------------------------------
# Exact integral-simplex slope: polynomial version of S_p.
# -----------------------------------------------------------------------------

def conv(p: Sequence[Complex], q: Sequence[Complex]) -> List[Complex]:
    out = [0.0 + 0.0j] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        if pi == 0:
            continue
        for j, qj in enumerate(q):
            out[i + j] += pi * qj
    return out


def linear_power(a: Complex, d: Complex, m: int) -> List[Complex]:
    if m <= 0:
        return [1.0 + 0.0j]
    return [comb(m, q) * (a ** (m - q)) * (d ** q) for q in range(m + 1)]


def integral_product_linear(a: Sequence[Complex], delta: Sequence[Complex], powers: Sequence[int]) -> Complex:
    coeffs = [1.0 + 0.0j]
    for ak, dk, mk in zip(a, delta, powers):
        coeffs = conv(coeffs, linear_power(ak, dk, mk))
    return sum(c / (q + 1) for q, c in enumerate(coeffs))


def Q_integral(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], F_anchor: Sequence[Complex] | None = None) -> Matrix:
    n = len(z)
    a = [complex(x) for x in anchor]
    delta = [complex(z[i]) - a[i] for i in range(n)]
    Q = [[0.0 + 0.0j] * n for _ in range(n)]
    for i, poly in enumerate(polys):
        for alpha, coeff in poly.items():
            c = complex(coeff)
            if c == 0:
                continue
            for j, aj in enumerate(alpha):
                if aj == 0:
                    continue
                powers = list(alpha)
                powers[j] -= 1
                Q[i][j] += c * aj * integral_product_linear(a, delta, powers)
    return Q


# -----------------------------------------------------------------------------
# dD lattice covectors: simplex, sparse, hypercube.
# -----------------------------------------------------------------------------

def multi_indices_le(total_degree: int, n: int) -> Iterable[Tuple[int, ...]]:
    if n == 0:
        yield ()
        return
    if n == 1:
        for k in range(total_degree + 1):
            yield (k,)
        return
    for k in range(total_degree + 1):
        for rest in multi_indices_le(total_degree - k, n - 1):
            yield (k,) + rest


@lru_cache(maxsize=512)
def cached_indices(total_degree: int, n: int) -> Tuple[Tuple[int, ...], ...]:
    return tuple(multi_indices_le(total_degree, n))


def simplex_grad(r: Sequence[Complex], p: int) -> Vector:
    n = len(r)
    grad = [0.0 + 0.0j] * n
    for alpha in cached_indices(max(1, p) - 1, n):
        for j, aj in enumerate(alpha):
            if aj == 0:
                continue
            term = aj * (r[j] ** (aj - 1) if aj > 1 else 1.0 + 0.0j)
            for k, ak in enumerate(alpha):
                if k != j and ak:
                    term *= r[k] ** ak
            grad[j] += term
    return grad


def geom_sum_grad_1d(x: Complex, p: int) -> Tuple[Complex, Complex]:
    S = 0.0 + 0.0j
    G = 0.0 + 0.0j
    term = 1.0 + 0.0j
    for k in range(max(1, p)):
        S += term
        if k >= 1:
            G += k * (x ** (k - 1))
        term *= x
    return S, G


def hypercube_grad(r: Sequence[Complex], p: int) -> Vector:
    sg = [geom_sum_grad_1d(x, p) for x in r]
    grad = [0.0 + 0.0j] * len(r)
    for j in range(len(r)):
        v = sg[j][1]
        for k in range(len(r)):
            if k != j:
                v *= sg[k][0]
        grad[j] = v
    return grad


def coefficient_activity(polys: System, n: int) -> List[float]:
    act = [0.0] * n
    for poly in polys:
        for alpha, coeff in poly.items():
            c = abs(coeff)
            for j, aj in enumerate(alpha):
                act[j] += c * aj
    m = max(act) if act else 0.0
    return [a / m if m > 0 and a > 0 else 0.05 for a in act]


def support_grad(polys: System, r: Sequence[Complex]) -> Vector:
    n = len(r)
    grad = [0.0 + 0.0j] * n
    maxc = max((abs(c) for poly in polys for c in poly.values()), default=1.0) or 1.0
    for poly in polys:
        for alpha, coeff in poly.items():
            w = abs(coeff) / maxc
            if w == 0:
                continue
            for j, aj in enumerate(alpha):
                if aj == 0:
                    continue
                term = w * aj
                for k, ak in enumerate(alpha):
                    power = ak - 1 if k == j else ak
                    if power:
                        term *= r[k] ** power
                grad[j] += term
    return grad


def normalized_delta(anchor: Sequence[Complex], z: Sequence[Complex]) -> Tuple[Vector, float]:
    d = [z[j] - anchor[j] for j in range(len(z))]
    s = max(1e-14, max(abs(x) for x in d))
    return [x / s for x in d], s


def covector(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], mode: str) -> Vector:
    n = len(z)
    delta = [z[j] - anchor[j] for j in range(n)]
    r, _ = normalized_delta(anchor, z)
    act = coefficient_activity(polys, n)
    p = system_degree(polys)
    if mode == "simplex":
        grad = simplex_grad(tuple(r), p)
    elif mode == "hypercube_cov":
        grad = hypercube_grad(tuple(r), p)
    elif mode == "support":
        grad = support_grad(polys, tuple(r))
    elif mode == "sparse":
        scores = sorted([(act[j] * (abs(delta[j]) + 1e-12), j) for j in range(n)], reverse=True)
        active = {j for _, j in scores[:max(1, int(math.ceil(math.sqrt(n))))]}
        return [delta[j].conjugate() * (1.0 + act[j]) if j in active else 0.0 + 0.0j for j in range(n)]
    else:
        grad = [0.0 + 0.0j] * n
    return [delta[j].conjugate() * (1.0 + abs(grad[j])) * (0.25 + act[j]) for j in range(n)]


def edge_frame(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], F_anchor: Sequence[Complex]) -> Matrix:
    n = len(z)
    B = [[0.0 + 0.0j] * n for _ in range(n)]
    for j in range(n):
        den = z[j] - anchor[j]
        if abs(den) < 1e-300:
            continue
        e = list(anchor)
        e[j] = z[j]
        Fe = F_eval(polys, e)
        if not finite_vec(Fe):
            continue
        for i in range(n):
            B[i][j] = (Fe[i] - F_anchor[i]) / den
    return B


def exactify(B: Matrix, F_anchor: Sequence[Complex], F_z: Sequence[Complex], anchor: Sequence[Complex], z: Sequence[Complex], w: Sequence[Complex]) -> Matrix:
    n = len(z)
    delta = [z[j] - anchor[j] for j in range(n)]
    Bd = [sum(B[i][j] * delta[j] for j in range(n)) for i in range(n)]
    defect = [F_z[i] - F_anchor[i] - Bd[i] for i in range(n)]
    den = sum(w[j] * delta[j] for j in range(n))
    if abs(den) < 1e-28:
        w = [d.conjugate() for d in delta]
        den = sum(w[j] * delta[j] for j in range(n))
    if abs(den) < 1e-28:
        return [row[:] for row in B]
    Q = [row[:] for row in B]
    for i in range(n):
        for j in range(n):
            Q[i][j] += defect[i] * w[j] / den
    return Q


def Q_edge_exact(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], F_anchor: Sequence[Complex], mode: str) -> Matrix:
    B = edge_frame(polys, anchor, z, F_anchor)
    Fz = F_eval(polys, z)
    if not finite_vec(Fz):
        return B
    return exactify(B, F_anchor, Fz, anchor, z, covector(polys, anchor, z, mode))


def Q_path(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], order: Sequence[int]) -> Matrix:
    n = len(z)
    Q = [[0.0 + 0.0j] * n for _ in range(n)]
    cur = list(z)
    Fcur = F_eval(polys, cur)
    for j in order:
        nxt = list(cur)
        nxt[j] = anchor[j]
        Fnxt = F_eval(polys, nxt)
        den = z[j] - anchor[j]
        if abs(den) >= 1e-300 and finite_vec(Fcur) and finite_vec(Fnxt):
            for i in range(n):
                Q[i][j] = (Fcur[i] - Fnxt[i]) / den
        cur, Fcur = nxt, Fnxt
    return Q


def Q_cube(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], max_orders: int = 8) -> Matrix:
    n = len(z)
    if n <= 4:
        orders = list(itertools.permutations(range(n)))[:max_orders]
    else:
        orders = [tuple(range(n)), tuple(reversed(range(n)))]
    out = [[0.0 + 0.0j] * n for _ in range(n)]
    inv = 1.0 / len(orders)
    for order in orders:
        Q = Q_path(polys, anchor, z, order)
        for i in range(n):
            for j in range(n):
                out[i][j] += inv * Q[i][j]
    return out


def Q_geometry(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], F_anchor: Sequence[Complex], mode: str) -> Matrix:
    if mode == "integral":
        return Q_integral(polys, anchor, z, F_anchor)
    if mode in {"simplex", "sparse", "support", "hypercube_cov"}:
        return Q_edge_exact(polys, anchor, z, F_anchor, mode)
    if mode == "path":
        return Q_path(polys, anchor, z, tuple(range(len(z))))
    if mode == "cube":
        return Q_cube(polys, anchor, z)
    raise ValueError(mode)


def pandrosion_h(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], F_anchor: Sequence[Complex], mode: str) -> Tuple[Vector, bool]:
    if not safe_z(anchor) or not safe_z(z) or not finite_vec(F_anchor):
        return list(z), False
    Q = Q_geometry(polys, anchor, z, F_anchor, mode)
    step = solve_linear(Q, [-v for v in F_anchor])
    if step is None:
        return list(z), False
    out = [anchor[i] + step[i] for i in range(len(z))]
    return (out, True) if safe_z(out) else (list(z), False)


def T2_step(polys: System, anchor: Sequence[Complex], z: Sequence[Complex], F_anchor: Sequence[Complex], mode: str) -> Tuple[Vector, bool]:
    s1, ok = pandrosion_h(polys, anchor, z, F_anchor, mode)
    if not ok:
        return list(z), False
    s2, ok = pandrosion_h(polys, anchor, s1, F_anchor, mode)
    if not ok:
        return s1, True
    out = []
    for k in range(len(z)):
        d0 = s1[k] - z[k]
        d2 = s2[k] - 2 * s1[k] + z[k]
        cand = s2[k] if abs(d2) < 1e-260 else z[k] - d0 * d0 / d2
        if abs(cand) > MAX_Z or abs(cand - z[k]) > 80.0 * max(1.0, abs(d0)):
            cand = s2[k]
        out.append(cand)
    return (out, True) if safe_z(out) else (s2, True)


def armijo_fallback(polys: System, z: Sequence[Complex]) -> Tuple[Vector, bool]:
    n = len(z)
    F0 = F_eval(polys, z)
    if not finite_vec(F0):
        return list(z), False
    r0 = sum(abs(v) ** 2 for v in F0)
    if r0 == 0:
        return list(z), True
    eps = min(1e-3, max(1e-7, math.sqrt(r0) * 1e-3))
    best = None
    best_pred = -1.0
    for j in range(n):
        for u in (1.0 + 0j, -1.0 + 0j, 1j, -1j):
            zp = list(z)
            zp[j] += eps * u
            Fp = F_eval(polys, zp)
            if not finite_vec(Fp):
                continue
            col = [(Fp[i] - F0[i]) / eps for i in range(n)]
            cn2 = sum(abs(c) ** 2 for c in col)
            if cn2 < 1e-24:
                continue
            ip = sum(col[i].conjugate() * F0[i] for i in range(n))
            pred = abs(ip) ** 2 / cn2
            alpha = -ip / cn2
            if pred > best_pred:
                best_pred = pred
                best = (j, u, alpha)
    if best is None:
        return list(z), False
    j, u, alpha = best
    for k in range(10):
        b = 0.5 ** k
        cand = list(z)
        cand[j] += b * alpha * u
        r = residual_norm2(polys, cand)
        if math.isfinite(r) and r < r0:
            return cand, True
    return list(z), False


def deterministic_anchor(z: Sequence[Complex], epoch: int) -> Vector:
    scale = 0.12 * max(1.0, max(abs(x) for x in z)) / (1.0 + 0.12 * epoch)
    out = []
    for j, zj in enumerate(z):
        th = 2.3999632297 * (epoch + 1) * (j + 1)
        out.append(zj + scale * complex(math.cos(th), math.sin(1.6180339887 * th)) + complex(0.01 * (j + 1), -0.007))
    return out


def accept_damped(polys: System, z: Sequence[Complex], cand: Sequence[Complex], best_r: float, best_z: Vector) -> Tuple[float, Vector]:
    if not safe_z(cand):
        return best_r, best_z
    direction = [cand[i] - z[i] for i in range(len(z))]
    for tau in (1.0, 0.5, 0.25, 0.125):
        trial = [z[i] + tau * direction[i] for i in range(len(z))]
        r = residual_norm(polys, trial)
        if math.isfinite(r) and r < best_r:
            best_r, best_z = r, trial
    return best_r, best_z


def corrector(polys: System, z_init: Sequence[Complex], tol: float = 1e-9, max_epochs: int = 20,
              modes: Sequence[str] = ("integral", "simplex", "sparse", "path", "cube", "support")) -> Tuple[Vector, bool, int]:
    z = list(z_init)
    anchor = deterministic_anchor(z, 0)
    for epoch in range(max_epochs):
        rz = residual_norm(polys, z)
        if rz < tol:
            return z, True, epoch
        Fa = F_eval(polys, anchor)
        if finite_vec(Fa) and max(abs(v) for v in Fa) < tol:
            return list(anchor), True, epoch
        best_r, best_z = rz, list(z)
        for mode in modes:
            # integral first usually gives the strongest mathematical Pandrosion step.
            for stepper in (pandrosion_h, T2_step):
                cand, ok = stepper(polys, anchor, z, Fa, mode)
                if ok:
                    best_r, best_z = accept_damped(polys, z, cand, best_r, best_z)
            if best_r < 0.35 * rz:
                break
        if best_r > 0.995 * rz:
            cand, ok = armijo_fallback(polys, z)
            if ok:
                best_r, best_z = accept_damped(polys, z, cand, best_r, best_z)
        if best_r >= rz:
            anchor = deterministic_anchor(z, epoch + 1)
            continue
        anchor = list(z)
        z = list(best_z)
    return z, residual_norm(polys, z) < tol, max_epochs


# -----------------------------------------------------------------------------
# Homotopy shell, copied in spirit from 064 but without Newton-ELS.
# -----------------------------------------------------------------------------

def add_scaled(out: Poly, poly: Poly, scale: Complex) -> None:
    for exp, coeff in poly.items():
        out[exp] = out.get(exp, 0.0 + 0.0j) + scale * coeff


def clean(poly: Poly, eps: float = 1e-14) -> Poly:
    return {k: v for k, v in poly.items() if abs(v) > eps}


def homotopy(target: System, start: System, t: float) -> System:
    out = []
    for f, g in zip(target, start):
        h: Poly = {}
        add_scaled(h, g, 1.0 - t)
        add_scaled(h, f, t)
        out.append(clean(h))
    return out


def start_system(degrees: Sequence[int], n: int) -> System:
    zero = tuple([0] * n)
    out = []
    for i, d in enumerate(degrees):
        exp = tuple(d if j == i else 0 for j in range(n))
        out.append({exp: 1.0 + 0.0j, zero: -1.0 + 0.0j})
    return out


def roots_unity(d: int) -> List[Complex]:
    return [cmath.exp(2j * math.pi * k / d) for k in range(d)]


def start_roots(degrees: Sequence[int]) -> List[Vector]:
    return [list(x) for x in iprod(*[roots_unity(d) for d in degrees])]


def phase_start_system(degrees: Sequence[int], n: int, phases: Sequence[Complex]) -> System:
    zero = tuple([0] * n)
    out = []
    for i, d in enumerate(degrees):
        exp = tuple(d if j == i else 0 for j in range(n))
        out.append({exp: 1.0 + 0.0j, zero: -complex(phases[i])})
    return out


def roots_of_phase(d: int, phase: Complex) -> List[Complex]:
    radius = abs(phase) ** (1.0 / d)
    theta = cmath.phase(phase) / d
    base = radius * cmath.exp(1j * theta)
    return [base * cmath.exp(2j * math.pi * k / d) for k in range(d)]


def phase_start_roots(degrees: Sequence[int], phases: Sequence[Complex]) -> List[Vector]:
    return [list(x) for x in iprod(*[roots_of_phase(d, ph) for d, ph in zip(degrees, phases)])]


def deterministic_phases(n: int, pass_index: int, seed: int) -> List[Complex]:
    if pass_index == 0:
        return [1.0 + 0.0j] * n
    # deterministic low-discrepancy-ish phases, no random module
    out = []
    for j in range(n):
        theta = (seed * 0.000123 + (pass_index + 1) * (j + 1) * 2.3999632297) % (2.0 * math.pi)
        radius = 0.75 + 0.45 * (0.5 + 0.5 * math.sin(seed * 0.017 + (j + 1) * (pass_index + 3)))
        out.append(radius * cmath.exp(1j * theta))
    return out


def track_one(target: System, start: System, z0: Sequence[Complex], tol: float = 1e-9,
              max_steps: int = 220) -> Tuple[Vector, bool, int, int]:
    z = list(z0)
    n = len(z)
    t = 0.0
    dt = 0.006
    prev_z = None
    prev_t = None
    fails = 0
    epochs = 0
    steps = 0
    while t < 1.0 - 1e-15 and steps < max_steps:
        steps += 1
        tnext = min(1.0, t + dt)
        if prev_z is None or prev_t is None or t == prev_t:
            pred = list(z)
        else:
            slope = [(z[i] - prev_z[i]) / (t - prev_t) for i in range(n)]
            pred = [z[i] + (tnext - t) * slope[i] for i in range(n)]
        H = homotopy(target, start, tnext)
        zc, ok, ep = corrector(H, pred, tol=tol, max_epochs=18)
        epochs += ep
        rh = residual_norm(H, zc)
        if ok and rh < 25.0 * tol:
            prev_z, prev_t = list(z), t
            z, t = list(zc), tnext
            dt = min(0.03, dt * 1.16)
            fails = max(0, fails - 1)
        else:
            fails += 1
            dt *= 0.5
            if dt < 8e-6 or fails > 70:
                return z, False, epochs, steps
    return z, residual_norm(target, z) < 1e-7, epochs, steps


def track_all(target: System, tol: float = 1e-9, passes: int = 1, seed: int = 67000) -> dict:
    n = len(target)
    degrees = [max(1, degree(p)) for p in target]
    bez = math.prod(degrees)
    found: List[Vector] = []
    okp = failp = epochs = steps = paths = 0
    for pass_index in range(max(1, passes)):
        phases = deterministic_phases(n, pass_index, seed + 17 * bez)
        start = phase_start_system(degrees, n, phases)
        roots0 = phase_start_roots(degrees, phases)
        for z0 in roots0:
            paths += 1
            z, ok, ep, st = track_one(target, start, z0, tol=tol)
            epochs += ep
            steps += st
            if ok:
                okp += 1
                if all(max(abs(z[i] - r[i]) for i in range(n)) > 1e-4 for r in found):
                    found.append(z)
                    if len(found) >= bez:
                        return {"roots": found, "bezout": bez, "coverage": 1.0,
                                "paths_ok": okp, "paths_fail": failp, "paths": paths,
                                "passes_used": pass_index + 1, "epochs": epochs, "steps": steps}
            else:
                failp += 1
    return {"roots": found, "bezout": bez, "coverage": len(found) / max(1, bez),
            "paths_ok": okp, "paths_fail": failp, "paths": paths,
            "passes_used": max(1, passes), "epochs": epochs, "steps": steps}



# -----------------------------------------------------------------------------
# 068 additions: fair 064-compatible systems + affine-chart recovery.
# -----------------------------------------------------------------------------

import argparse
import random


def gen_random_poly_system(n: int, d: int, seed: int) -> System:
    """Same random-system generator as flow/064 for fair A/B benchmarks."""
    rng = random.Random(seed)
    polys: System = []
    for i in range(n):
        poly: Poly = {}
        for alpha in iprod(*[range(d + 1) for _ in range(n)]):
            if sum(alpha) > d:
                continue
            if rng.random() < 0.7:
                poly[tuple(alpha)] = complex(rng.gauss(0.0, 1.0), 0.15 * rng.gauss(0.0, 1.0))
        diag = tuple(d if k == i else 0 for k in range(n))
        poly[diag] = poly.get(diag, 0.0 + 0.0j) + 1.0
        m = max(abs(v) for v in poly.values()) or 1.0
        polys.append({k: v / m for k, v in poly.items()})
    return polys


def dict_poly_mul(p: Poly, q: Poly, n: int, eps: float = 1e-14) -> Poly:
    out: Poly = {}
    for a, ca in p.items():
        if abs(ca) <= eps:
            continue
        for b, cb in q.items():
            if abs(cb) <= eps:
                continue
            e = tuple(a[i] + b[i] for i in range(n))
            out[e] = out.get(e, 0.0 + 0.0j) + ca * cb
    return {e: c for e, c in out.items() if abs(c) > eps}


def linear_form_power(const: Complex, coeffs: Sequence[Complex], power: int) -> Poly:
    """Return (const + sum_j coeffs[j] y_j)^power as a sparse dict poly."""
    n = len(coeffs)
    one = tuple([0] * n)
    out: Poly = {one: 1.0 + 0.0j}
    if power <= 0:
        return out
    lin: Poly = {one: complex(const)}
    for j, c in enumerate(coeffs):
        if abs(c) > 0:
            e = [0] * n
            e[j] = 1
            lin[tuple(e)] = lin.get(tuple(e), 0.0 + 0.0j) + complex(c)
    for _ in range(power):
        out = dict_poly_mul(out, lin, n)
    return out


def compose_affine_poly(poly: Poly, A: Matrix, b: Sequence[Complex]) -> Poly:
    """Compose f(z) by z=A*y+b.  Used for affine chart recovery."""
    n = len(b)
    one = tuple([0] * n)
    out: Poly = {}
    for alpha, coeff in poly.items():
        term: Poly = {one: complex(coeff)}
        for k, power in enumerate(alpha):
            if power:
                term = dict_poly_mul(term, linear_form_power(b[k], A[k], power), n)
        for e, c in term.items():
            out[e] = out.get(e, 0.0 + 0.0j) + c
    return {e: c for e, c in out.items() if abs(c) > 1e-12}


def compose_affine_system(polys: System, A: Matrix, b: Sequence[Complex]) -> System:
    return [compose_affine_poly(p, A, b) for p in polys]


def affine_map(A: Matrix, b: Sequence[Complex], y: Sequence[Complex]) -> Vector:
    n = len(b)
    return [b[i] + sum(A[i][j] * y[j] for j in range(n)) for i in range(n)]


def deterministic_affine_chart(n: int, chart_index: int, seed: int) -> Tuple[Matrix, Vector]:
    """Near-identity dense affine charts.  chart_index=0 is the identity chart."""
    if chart_index == 0:
        return [[1.0 + 0.0j if i == j else 0.0 + 0.0j for j in range(n)] for i in range(n)], [0.0 + 0.0j] * n
    rng = random.Random(seed + 997 * chart_index + 13 * n)
    eps = 0.25 if chart_index <= 3 else 0.45
    A: Matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            noise = eps * complex(rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)) / math.sqrt(2.0 * max(1, n))
            row.append((1.0 + 0.0j if i == j else 0.0 + 0.0j) + noise)
        A.append(row)
    b = [0.25 * complex(rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)) / math.sqrt(2.0) for _ in range(n)]
    return A, b


def unique_append_root(found: List[Vector], z: Sequence[Complex], sep: float = 1e-4) -> bool:
    n = len(z)
    if not safe_z(z):
        return False
    if all(max(abs(z[i] - r[i]) for i in range(n)) > sep for r in found):
        found.append(list(z))
        return True
    return False


def polish_original_root(polys: System, z: Sequence[Complex], tol: float = 1e-10) -> Vector | None:
    """Pandrosion-only polish after mapping back from an affine chart."""
    r = residual_norm(polys, z)
    if math.isfinite(r) and r < 1e-8:
        return list(z)
    zp, ok, _ = corrector(polys, z, tol=tol, max_epochs=24)
    if ok and residual_norm(polys, zp) < 1e-7:
        return zp
    if math.isfinite(r) and r < 1e-7:
        return list(z)
    return None


def track_all_affine_charts(target: System, tol: float = 1e-9, identity_passes: int = 2,
                            affine_charts: int = 4, chart_passes: int = 1,
                            seed: int = 68000, stop_at_bezout: bool = True) -> dict:
    """Track identity phase starts first, then affine-coordinate recovery charts.

    No Newton-ELS is used.  Each affine chart solves F(Ay+b)=0 and maps y back
    to z=Ay+b.  This changes the straight homotopy path geometry without using
    Lairez gamma rotation.
    """
    n = len(target)
    degrees = [max(1, degree(p)) for p in target]
    bez = math.prod(degrees)
    found: List[Vector] = []
    paths_ok = paths_fail = paths_total = steps = epochs = 0
    charts_used = 0
    chart_log: List[Tuple[int, int, int, int, float]] = []

    # Chart 0: original coordinates, with phase starts.
    t0 = time.time()
    res0 = track_all(target, tol=tol, passes=max(1, identity_passes), seed=seed)
    charts_used = 1
    paths_ok += res0["paths_ok"]
    paths_fail += res0["paths_fail"]
    paths_total += res0.get("paths", res0["paths_ok"] + res0["paths_fail"])
    steps += res0.get("steps", 0)
    epochs += res0.get("epochs", 0)
    added = 0
    for z in res0["roots"]:
        zp = polish_original_root(target, z, tol=tol)
        if zp is not None and unique_append_root(found, zp):
            added += 1
    chart_log.append((0, len(res0["roots"]), added, len(found), time.time() - t0))
    if stop_at_bezout and len(found) >= bez:
        return {"roots": found, "bezout": bez, "coverage": 1.0,
                "paths_ok": paths_ok, "paths_fail": paths_fail, "paths": paths_total,
                "steps": steps, "epochs": epochs, "charts_used": charts_used,
                "chart_log": chart_log}

    for chart_index in range(1, max(1, affine_charts) + 1):
        tc = time.time()
        A, b = deterministic_affine_chart(n, chart_index, 80000 + 101 * (n + bez))
        transformed = compose_affine_system(target, A, b)
        res = track_all(transformed, tol=tol, passes=max(1, chart_passes), seed=seed + 1009 * chart_index)
        charts_used = chart_index + 1
        paths_ok += res["paths_ok"]
        paths_fail += res["paths_fail"]
        paths_total += res.get("paths", res["paths_ok"] + res["paths_fail"])
        steps += res.get("steps", 0)
        epochs += res.get("epochs", 0)
        added = 0
        for y in res["roots"]:
            z = affine_map(A, b, y)
            zp = polish_original_root(target, z, tol=tol)
            if zp is not None and unique_append_root(found, zp):
                added += 1
        chart_log.append((chart_index, len(res["roots"]), added, len(found), time.time() - tc))
        if stop_at_bezout and len(found) >= bez:
            break

    return {"roots": found, "bezout": bez, "coverage": len(found) / max(1, bez),
            "paths_ok": paths_ok, "paths_fail": paths_fail, "paths": paths_total,
            "steps": steps, "epochs": epochs, "charts_used": charts_used,
            "chart_log": chart_log}


def telescope_error(polys: System, mode: str, anchor: Sequence[Complex], z: Sequence[Complex]) -> float:
    Fa = F_eval(polys, anchor)
    Fz = F_eval(polys, z)
    Q = Q_geometry(polys, anchor, z, Fa, mode)
    d = [z[j] - anchor[j] for j in range(len(z))]
    Qd = [sum(Q[i][j] * d[j] for j in range(len(z))) for i in range(len(z))]
    return max(abs(Qd[i] - (Fz[i] - Fa[i])) for i in range(len(z)))


def smoke() -> None:
    print("[068 smoke] telescope identity for exact dD polynomial Pandrosion slopes")
    modes = ["integral", "simplex", "sparse", "hypercube_cov", "support", "path", "cube"]
    for n, d in [(2, 3), (3, 2), (4, 2)]:
        polys = gen_random_poly_system(n, d, 67000 + 10 * n + d)
        a = [complex(math.sin(j + 1), math.cos(2 * j + 1)) for j in range(n)]
        z = [complex(math.cos(3 * j + 1), math.sin(4 * j + 2)) for j in range(n)]
        errs = [telescope_error(polys, m, a, z) for m in modes]
        print(f"  n={n}, d={d}, worst={max(errs):.2e} | " + ", ".join(f"{m}:{e:.1e}" for m, e in zip(modes, errs)))


def parse_case(s: str) -> Tuple[int, int]:
    if "," in s:
        a, b = s.split(",", 1)
    elif "x" in s:
        a, b = s.split("x", 1)
    else:
        raise argparse.ArgumentTypeError("case must look like 4,2 or 4x2")
    return int(a), int(b)


def benchmark(cases: Sequence[Tuple[int, int]], identity_passes: int = 2,
              affine_charts: int = 4, chart_passes: int = 1,
              show_chart_log: bool = False) -> None:
    print("=" * 132)
    print("068 -- dD polynomial Pandrosion + phase starts + affine-chart recovery; no Newton-ELS")
    print("=" * 132)
    print("Same random system generator/seed family as 064.  Q modes: integral, simplex, sparse, hypercube/support, path/cube.")
    print(f"  {'(n,d)':>8} {'Bez':>5} | {'roots':>5} {'cov%':>6} {'paths ok/fail':>14} {'paths':>6} {'charts':>6} {'steps':>7} {'epochs':>8} {'time':>7}")
    print("-" * 132)
    for n, d in cases:
        polys = gen_random_poly_system(n, d, seed=61000 + 100 * n + d)
        t0 = time.time()
        result = track_all_affine_charts(polys, identity_passes=identity_passes,
                                         affine_charts=affine_charts,
                                         chart_passes=chart_passes,
                                         seed=68000 + 100 * n + d)
        elapsed = time.time() - t0
        print(f"  ({n:>2},{d:>2}) {result['bezout']:>5} | {len(result['roots']):>5} {100*result['coverage']:>5.1f}% "
              f"{result['paths_ok']:>4}/{result['paths_fail']:<4} {result['paths']:>6} "
              f"{result['charts_used']:>6} {result['steps']:>7} {result['epochs']:>8} {elapsed:>6.1f}s")
        if show_chart_log:
            for ci, chart_roots, added, total, dt in result["chart_log"]:
                print(f"      chart {ci:>2}: chart_roots={chart_roots:>3}, added={added:>2}, total={total:>3}, time={dt:>5.1f}s")
    print("=" * 132)


def main() -> None:
    parser = argparse.ArgumentParser(description="flow/068 polynomial Pandrosion affine-chart benchmark")
    parser.add_argument("--cases", nargs="*", type=parse_case,
                        default=[(2, 2), (2, 3), (3, 2), (2, 4), (4, 2)],
                        help="cases like 2,3 or 4x2")
    parser.add_argument("--passes", type=int, default=2, help="phase passes in the original coordinates")
    parser.add_argument("--chart-passes", type=int, default=1, help="phase passes per affine recovery chart")
    parser.add_argument("--charts", type=int, default=4, help="number of affine recovery charts after the identity chart")
    parser.add_argument("--no-smoke", action="store_true")
    parser.add_argument("--chart-log", action="store_true")
    args = parser.parse_args()
    if not args.no_smoke:
        smoke()
        print()
    benchmark(args.cases, identity_passes=args.passes, affine_charts=args.charts,
              chart_passes=args.chart_passes, show_chart_log=args.chart_log)


if __name__ == "__main__":
    main()
