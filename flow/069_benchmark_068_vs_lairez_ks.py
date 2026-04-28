"""
FLOW 069 -- benchmark 068 polynomial Pandrosion vs Lairez-style gamma homotopy

Purpose
-------
This file is a reusable benchmark harness for comparing:

  1. flow/068: dD polynomial Pandrosion with phase starts + affine-chart
     recovery, no Newton-ELS.

  2. Lairez-style baseline: total-degree start system + random gamma homotopy
     + analytic Newton corrector.  This is not the official HomotopyContinuation
     or lairez library; it is a compact local implementation of the classical
     randomized gamma total-degree continuation used as the reference baseline.

System families
---------------
  dense064   : same dense random generator/seed convention as flow/064/068.
  ks         : affine Kostlan-Shub-Smale/Kostlan systems, coefficients scaled by
               sqrt(d!/(alpha!*(d-|alpha|)!)).
  sparse_ks  : sparse support sampled from the same Kostlan scaling, with
               diagonal degree-d and constant terms forced in each equation.
  high_ks    : shorthand suite of higher-degree 2D/3D Kostlan cases.

All reported coverage is against the total-degree Bezout number prod_i deg(F_i).
For generic dense/Kostlan systems this is the expected number of isolated affine
complex roots.  Sparse systems can have structure-dependent degeneracies, so
coverage there is best read as a total-degree stress score.

Examples
--------
Quick validation:
  python 069_benchmark_068_vs_lairez_ks.py --suite quick --charts 1 --passes 1

Fuller benchmark with CSV/JSON logs:
  python 069_benchmark_068_vs_lairez_ks.py --suite full --seeds 3 --charts 2 --passes 2 \
    --csv /mnt/data/069_results.csv --json /mnt/data/069_results.json

High-degree KS only:
  python 069_benchmark_068_vs_lairez_ks.py --family ks --cases 2,8 2,10 2,12 3,4 3,5 \
    --seeds 2 --charts 2 --passes 2
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import random
import statistics
import sys
import time
import signal
import multiprocessing as mp
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from itertools import product as iprod
from math import comb, factorial
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

Complex = complex
Exponent = Tuple[int, ...]
Poly = Dict[Exponent, Complex]
System = List[Poly]
Vector = List[Complex]
Matrix = List[List[Complex]]

HERE = Path(__file__).resolve().parent
M068_PATH = HERE / "068_pandrosion_affine_chart_no_newton.py"


def load_068():
    spec = importlib.util.spec_from_file_location("flow068", str(M068_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import 068 from {M068_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["flow068"] = module
    spec.loader.exec_module(module)
    return module


m068 = load_068()


# -----------------------------------------------------------------------------
# Basic polynomial helpers, independent of 068 so the Lairez-style baseline is
# explicit and auditable.
# -----------------------------------------------------------------------------

def finite(v: Complex) -> bool:
    z = complex(v)
    return math.isfinite(z.real) and math.isfinite(z.imag)


def finite_vec(z: Sequence[Complex]) -> bool:
    return all(finite(v) for v in z)


def degree(poly: Poly) -> int:
    return max((sum(e) for e in poly), default=0)


def degrees(polys: System) -> List[int]:
    return [max(1, degree(p)) for p in polys]


def bezout(polys: System) -> int:
    out = 1
    for d in degrees(polys):
        out *= d
    return out


def eval_poly(poly: Poly, z: Sequence[Complex]) -> Complex:
    total = 0.0 + 0.0j
    try:
        for exp, coeff in poly.items():
            term = complex(coeff)
            for zj, aj in zip(z, exp):
                if aj:
                    term *= zj ** aj
            total += term
    except (OverflowError, ZeroDivisionError):
        return complex("inf")
    return total if finite(total) else complex("inf")


def F_eval(polys: System, z: Sequence[Complex]) -> Vector:
    return [eval_poly(p, z) for p in polys]


def residual_norm(polys: System, z: Sequence[Complex]) -> float:
    if not finite_vec(z):
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
    M = [list(row) + [complex(b[i])] for i, row in enumerate(A)]
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
    return x if finite_vec(x) else None


def jacobian_eval(polys: System, z: Sequence[Complex]) -> Matrix:
    n = len(z)
    J = [[0.0 + 0.0j for _ in range(n)] for _ in range(n)]
    try:
        for i, poly in enumerate(polys):
            for alpha, coeff in poly.items():
                c = complex(coeff)
                if c == 0:
                    continue
                for j, aj in enumerate(alpha):
                    if aj == 0:
                        continue
                    term = c * aj
                    for k, ak in enumerate(alpha):
                        pk = ak - 1 if k == j else ak
                        if pk:
                            term *= z[k] ** pk
                    J[i][j] += term
    except (OverflowError, ZeroDivisionError):
        return [[complex("inf") for _ in range(n)] for _ in range(n)]
    return J


def add_scaled(out: Poly, poly: Poly, scale: Complex) -> None:
    for e, c in poly.items():
        out[e] = out.get(e, 0.0 + 0.0j) + scale * c


def clean(poly: Poly, eps: float = 1e-14) -> Poly:
    return {e: c for e, c in poly.items() if abs(c) > eps}


def scale_system(polys: System, gammas: Sequence[Complex]) -> System:
    return [{e: gammas[i] * c for e, c in p.items()} for i, p in enumerate(polys)]


def homotopy_polys(target_gamma: System, start: System, t: float) -> System:
    out: System = []
    for f, g in zip(target_gamma, start):
        h: Poly = {}
        add_scaled(h, g, 1.0 - t)
        add_scaled(h, f, t)
        out.append(clean(h))
    return out


def total_degree_start_system(degs: Sequence[int], n: int) -> System:
    zero = tuple([0] * n)
    out: System = []
    for i, d in enumerate(degs):
        exp = tuple(d if j == i else 0 for j in range(n))
        out.append({exp: 1.0 + 0.0j, zero: -1.0 + 0.0j})
    return out


def roots_unity(d: int) -> List[Complex]:
    return [complex(math.cos(2 * math.pi * k / d), math.sin(2 * math.pi * k / d)) for k in range(d)]


def start_roots(degs: Sequence[int]) -> List[Vector]:
    return [list(x) for x in iprod(*[roots_unity(d) for d in degs])]


def unique_append(found: List[Vector], z: Sequence[Complex], sep: float = 1e-4) -> bool:
    n = len(z)
    if not finite_vec(z):
        return False
    if all(max(abs(z[i] - r[i]) for i in range(n)) > sep for r in found):
        found.append(list(z))
        return True
    return False


# -----------------------------------------------------------------------------
# System generators.
# -----------------------------------------------------------------------------

def multi_indices_le(d: int, n: int) -> Iterable[Exponent]:
    if n == 0:
        yield ()
        return
    if n == 1:
        for k in range(d + 1):
            yield (k,)
        return
    for k in range(d + 1):
        for rest in multi_indices_le(d - k, n - 1):
            yield (k,) + rest


def complex_gaussian(rng: random.Random, sigma: float = 1.0) -> Complex:
    # E|z|^2 = sigma^2
    return sigma * complex(rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)) / math.sqrt(2.0)


def normalize_system(polys: System, mode: str = "rms") -> System:
    out: System = []
    for p in polys:
        if not p:
            out.append(p)
            continue
        if mode == "max":
            s = max(abs(c) for c in p.values()) or 1.0
        else:
            s = math.sqrt(sum(abs(c) ** 2 for c in p.values()) / max(1, len(p))) or 1.0
        out.append({e: c / s for e, c in p.items()})
    return out


def kostlan_scale(alpha: Exponent, d: int) -> float:
    # Affine chart of homogeneous degree-d Kostlan polynomial in n+1 variables:
    # scale^2 = d! / ((d-|alpha|)! prod_j alpha_j!).
    rem = d - sum(alpha)
    denom = factorial(rem)
    for a in alpha:
        denom *= factorial(a)
    return math.sqrt(factorial(d) / denom)


def gen_ks_system(n: int, d: int, seed: int) -> System:
    rng = random.Random(seed)
    polys: System = []
    alphas = list(multi_indices_le(d, n))
    for _ in range(n):
        p: Poly = {}
        for alpha in alphas:
            p[alpha] = complex_gaussian(rng, kostlan_scale(alpha, d))
        polys.append(p)
    return normalize_system(polys, "rms")


def gen_sparse_ks_system(n: int, d: int, seed: int, density: float | None = None) -> System:
    rng = random.Random(seed)
    if density is None:
        # Decrease density as support grows, but keep enough mixed terms to stress the solver.
        density = min(0.45, max(0.12, 8.0 / max(1, comb(n + d, d))))
    zero = tuple([0] * n)
    polys: System = []
    for i in range(n):
        p: Poly = {}
        for alpha in multi_indices_le(d, n):
            force = alpha == zero or alpha == tuple(d if j == i else 0 for j in range(n))
            if force or rng.random() < density:
                p[alpha] = complex_gaussian(rng, kostlan_scale(alpha, d))
        # Guarantee actual degree d in every polynomial.
        diag = tuple(d if j == i else 0 for j in range(n))
        p[diag] = p.get(diag, 0.0 + 0.0j) + complex_gaussian(rng, kostlan_scale(diag, d))
        if zero not in p:
            p[zero] = complex_gaussian(rng, kostlan_scale(zero, d))
        polys.append(p)
    return normalize_system(polys, "rms")


def gen_dense064_system(n: int, d: int, seed: int) -> System:
    return m068.gen_random_poly_system(n, d, seed)


def gen_system(family: str, n: int, d: int, seed: int) -> System:
    family = family.lower()
    if family in {"dense", "dense064", "064"}:
        return gen_dense064_system(n, d, seed)
    if family in {"ks", "kostlan", "kss"}:
        return gen_ks_system(n, d, seed)
    if family in {"sparse", "sparse_ks", "sparse-kostlan"}:
        return gen_sparse_ks_system(n, d, seed)
    if family in {"high_ks", "high-kostlan"}:
        return gen_ks_system(n, d, seed)
    raise ValueError(f"unknown family: {family}")


def term_count(polys: System) -> int:
    return sum(len(p) for p in polys)


# -----------------------------------------------------------------------------
# Lairez-style total-degree gamma homotopy + analytic Newton.
# -----------------------------------------------------------------------------

@dataclass
class TrackResult:
    roots: List[Vector]
    bezout: int
    coverage: float
    paths_ok: int
    paths_fail: int
    paths: int
    steps: int
    newton_iters: int
    max_residual: float


def deterministic_gamma_vector(n: int, seed: int, retry_index: int = 0) -> List[Complex]:
    rng = random.Random(seed + 104729 * retry_index)
    # Nonzero random phases.  Vector phases are at least as generic as a scalar
    # phase and keep the zero set unchanged equation-wise.
    return [complex(math.cos(rng.uniform(0.0, 2.0 * math.pi)), math.sin(rng.uniform(0.0, 2.0 * math.pi))) for _ in range(n)]


def newton_correct_analytic(polys: System, z_init: Sequence[Complex], tol: float = 1e-9,
                            max_iter: int = 12) -> Tuple[Vector, bool, int]:
    z = list(z_init)
    r = residual_norm(polys, z)
    if r < tol:
        return z, True, 0
    if not math.isfinite(r):
        return z, False, 0
    it_used = 0
    for it in range(1, max_iter + 1):
        it_used = it
        Fz = F_eval(polys, z)
        J = jacobian_eval(polys, z)
        if not finite_vec(Fz) or any(not finite_vec(row) for row in J):
            return z, False, it_used
        step = solve_linear(J, [-v for v in Fz])
        if step is None:
            return z, False, it_used
        r2 = sum(abs(v) ** 2 for v in Fz)
        best_z = None
        best_r2 = r2
        # Simple damped Newton / line search.  This is intentionally standard,
        # not a Pandrosion step.
        for k in range(12):
            tau = 0.5 ** k
            cand = [z[j] + tau * step[j] for j in range(len(z))]
            if not finite_vec(cand):
                continue
            rc2 = residual_norm2(polys, cand)
            if math.isfinite(rc2) and rc2 < best_r2:
                best_r2 = rc2
                best_z = cand
                if rc2 < tol * tol:
                    return cand, True, it_used
        if best_z is None:
            return z, False, it_used
        z = best_z
        if math.sqrt(best_r2) < tol:
            return z, True, it_used
        if best_r2 > 0.999999 * r2:
            return z, False, it_used
    return z, residual_norm(polys, z) < tol, it_used


def track_one_lairez(target: System, start: System, target_gamma: System, z0: Sequence[Complex],
                     tol: float = 1e-9, max_steps: int = 420,
                     max_newton_iter: int = 12) -> Tuple[Vector, bool, int, int]:
    z = list(z0)
    n = len(z)
    t = 0.0
    dt = 0.01
    prev_z = None
    prev_t = None
    fails = 0
    steps = 0
    newton_iters = 0
    while t < 1.0 - 1e-15 and steps < max_steps:
        steps += 1
        tnext = min(1.0, t + dt)
        # Euler predictor from dH/dz * dz/dt + dH/dt = 0.
        Hcur = homotopy_polys(target_gamma, start, t)
        Jcur = jacobian_eval(Hcur, z)
        dHdt = [fg - gs for fg, gs in zip(F_eval(target_gamma, z), F_eval(start, z))]
        dzdt = solve_linear(Jcur, [-v for v in dHdt])
        if dzdt is not None and finite_vec(dzdt) and max(abs(v) for v in dzdt) * (tnext - t) < 8.0 * max(1.0, max(abs(x) for x in z)):
            pred = [z[j] + (tnext - t) * dzdt[j] for j in range(n)]
        elif prev_z is None or prev_t is None or t == prev_t:
            pred = list(z)
        else:
            slope = [(z[j] - prev_z[j]) / (t - prev_t) for j in range(n)]
            pred = [z[j] + (tnext - t) * slope[j] for j in range(n)]
        H = homotopy_polys(target_gamma, start, tnext)
        zc, ok, nit = newton_correct_analytic(H, pred, tol=tol, max_iter=max_newton_iter)
        newton_iters += nit
        rh = residual_norm(H, zc)
        if ok and rh < 20.0 * tol:
            prev_z, prev_t = list(z), t
            z, t = list(zc), tnext
            dt = min(0.025, dt * 1.15)
            fails = max(0, fails - 1)
        else:
            fails += 1
            dt *= 0.5
            if dt < 2e-6 or fails > 80:
                return z, False, steps, newton_iters
    return z, residual_norm(target, z) < 1e-7, steps, newton_iters


def track_all_lairez(target: System, seed: int, tol: float = 1e-9,
                     max_steps: int = 420, max_newton_iter: int = 12,
                     retries: int = 3) -> TrackResult:
    """Lairez-style gamma homotopy with deterministic gamma retries.

    A production continuation code should not need repeated retries on generic
    systems; the retries here compensate for this compact local tracker's
    simple predictor/corrector and make the baseline harder to beat.
    """
    n = len(target)
    degs = degrees(target)
    B = math.prod(degs)
    start = total_degree_start_system(degs, n)
    roots0 = start_roots(degs)
    found: List[Vector] = []
    okp = failp = steps = nit = 0
    for retry in range(max(1, retries)):
        if len(found) >= B:
            break
        gammas = deterministic_gamma_vector(n, seed, retry)
        target_gamma = scale_system(target, gammas)
        for z0 in roots0:
            if len(found) >= B:
                break
            z, ok, st, ni = track_one_lairez(target, start, target_gamma, z0, tol=tol,
                                             max_steps=max_steps, max_newton_iter=max_newton_iter)
            steps += st
            nit += ni
            if ok:
                okp += 1
                unique_append(found, z)
            else:
                failp += 1
    if found:
        max_res = max(residual_norm(target, z) for z in found)
    else:
        max_res = float("inf")
    return TrackResult(roots=found, bezout=B, coverage=len(found) / max(1, B), paths_ok=okp,
                       paths_fail=failp, paths=okp + failp, steps=steps,
                       newton_iters=nit, max_residual=max_res)



# -----------------------------------------------------------------------------
# Benchmark orchestration.
# -----------------------------------------------------------------------------

@dataclass
class BenchRow:
    family: str
    n: int
    d: int
    seed_index: int
    seed: int
    terms: int
    bezout: int
    alg: str
    roots: int
    coverage: float
    paths_ok: int
    paths_fail: int
    paths: int
    charts: int
    steps: int
    epochs_or_newton: int
    max_residual: float
    seconds: float




class TimeLimitExpired(Exception):
    pass


@contextmanager
def time_limit(seconds: float | None):
    """Unix wall-clock time limit for one path/corrector call."""
    if seconds is None or seconds <= 0:
        yield
        return
    old_handler = signal.getsignal(signal.SIGALRM)
    def handler(signum, frame):
        raise TimeLimitExpired()
    signal.signal(signal.SIGALRM, handler)
    signal.setitimer(signal.ITIMER_REAL, float(seconds))
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)



def track_one_068_fast(target: System, start: System, z0: Sequence[Complex], tol: float = 1e-9,
                       max_steps: int = 180,
                       modes: Sequence[str] = ("integral", "simplex", "sparse", "path")) -> Tuple[Vector, bool, int, int]:
    """Budgeted 068 path tracker.

    This is the same homotopy shell as flow/068, but the local Pandrosion
    corrector is restricted to the fast exact slopes.  It keeps the benchmark
    usable on high-degree KS systems and avoids spending benchmark time in the
    heavier support/cube portfolio unless the standalone 068 file is desired.
    """
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
        H = m068.homotopy(target, start, tnext)
        zc, ok, ep = m068.corrector(H, pred, tol=tol, max_epochs=12, modes=modes)
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
            if dt < 8e-6 or fails > 50:
                return z, False, epochs, steps
    return z, residual_norm(target, z) < 1e-7, epochs, steps

def _track_phase_pass_068(target: System, pass_index: int, seed: int, tol: float,
                          max_steps: int, path_timeout: float,
                          found: List[Vector], stop_at: int) -> Tuple[int, int, int, int, int]:
    """Manual version of 068.track_all to avoid unbounded suite-level stalls."""
    n = len(target)
    degs = degrees(target)
    phases = m068.deterministic_phases(n, pass_index, seed + 17 * math.prod(degs))
    start = m068.phase_start_system(degs, n, phases)
    roots0 = m068.phase_start_roots(degs, phases)
    okp = failp = paths = steps = epochs = 0
    for z0 in roots0:
        if len(found) >= stop_at:
            break
        paths += 1
        try:
            with time_limit(path_timeout):
                z, ok, ep, st = track_one_068_fast(target, start, z0, tol=tol, max_steps=max_steps)
        except TimeLimitExpired:
            failp += 1
            continue
        steps += st
        epochs += ep
        if ok and residual_norm(target, z) < 1e-7:
            okp += 1
            unique_append(found, z)
        else:
            failp += 1
    return okp, failp, paths, steps, epochs


def _polish_068(target: System, z: Sequence[Complex], tol: float, timeout: float) -> Vector | None:
    r = residual_norm(target, z)
    if math.isfinite(r) and r < 1e-7:
        return list(z)
    try:
        with time_limit(timeout):
            zp, ok, _ = m068.corrector(target, z, tol=tol, max_epochs=18)
    except TimeLimitExpired:
        return None
    if ok and residual_norm(target, zp) < 1e-7:
        return zp
    return None


def run_068(target: System, seed: int, passes: int, charts: int, chart_passes: int,
            tol: float, max_steps: int = 220, path_timeout: float = 8.0) -> Tuple[dict, float]:
    """Run the 068 method, but with explicit per-path time limits.

    This keeps long high-degree stress tests honest: a stalled Pandrosion path is
    counted as a failed path rather than freezing the whole benchmark.
    """
    t0 = time.time()
    n = len(target)
    B = bezout(target)
    found: List[Vector] = []
    okp = failp = paths = steps = epochs = 0
    charts_used = 1

    for pi in range(max(1, passes)):
        a, b, c, dsteps, e = _track_phase_pass_068(target, pi, seed, tol, max_steps, path_timeout, found, B)
        okp += a; failp += b; paths += c; steps += dsteps; epochs += e
        if len(found) >= B:
            break

    # Affine recovery charts, same geometric idea as 068, but budgeted.
    for ci in range(1, max(0, charts) + 1):
        if len(found) >= B:
            break
        charts_used = ci + 1
        A, bvec = m068.deterministic_affine_chart(n, ci, 80000 + 101 * (n + B))
        transformed = m068.compose_affine_system(target, A, bvec)
        chart_found: List[Vector] = []
        for pi in range(max(1, chart_passes)):
            a, b, c, dsteps, e = _track_phase_pass_068(transformed, pi, seed + 1009 * ci, tol,
                                                       max_steps, path_timeout, chart_found, B)
            okp += a; failp += b; paths += c; steps += dsteps; epochs += e
        for y in chart_found:
            z = m068.affine_map(A, bvec, y)
            zp = _polish_068(target, z, tol, timeout=max(2.0, path_timeout / 2.0))
            if zp is not None:
                unique_append(found, zp)
            if len(found) >= B:
                break

    elapsed = time.time() - t0
    max_res = max((residual_norm(target, z) for z in found), default=float("inf"))
    return {
        "roots": found,
        "bezout": B,
        "coverage": len(found) / max(1, B),
        "paths_ok": okp,
        "paths_fail": failp,
        "paths": paths,
        "charts_used": charts_used,
        "steps": steps,
        "epochs": epochs,
        "max_residual": max_res,
    }, elapsed


def run_lairez(target: System, seed: int, tol: float, max_steps: int,
               max_newton_iter: int, retries: int = 3) -> Tuple[TrackResult, float]:
    t0 = time.time()
    result = track_all_lairez(target, seed=seed, tol=tol,
                              max_steps=max_steps, max_newton_iter=max_newton_iter, retries=retries)
    return result, time.time() - t0




def _method_worker(conn, alg: str, target: System, kwargs: dict) -> None:
    try:
        if alg == "068":
            result, seconds = run_068(target, **kwargs)
            conn.send({"ok": True, "alg": alg, "result": result, "seconds": seconds})
        elif alg == "lairez":
            result, seconds = run_lairez(target, **kwargs)
            conn.send({"ok": True, "alg": alg, "result": result, "seconds": seconds})
        else:
            conn.send({"ok": False, "alg": alg, "error": f"unknown algorithm {alg}"})
    except BaseException as exc:
        conn.send({"ok": False, "alg": alg, "error": repr(exc)})
    finally:
        try:
            conn.close()
        except Exception:
            pass


def run_method_with_timeout(alg: str, target: System, kwargs: dict, timeout: float) -> Tuple[object | None, float, str]:
    """Run one method in a child process and kill it on timeout.

    This is intentionally process-based rather than signal-based: in this
    execution environment SIGALRM is unreliable, while Process.terminate/kill
    works.  Use timeout<=0 to run directly in the current process.
    """
    if timeout <= 0:
        if alg == "068":
            result, seconds = run_068(target, **kwargs)
        elif alg == "lairez":
            result, seconds = run_lairez(target, **kwargs)
        else:
            raise ValueError(alg)
        return result, seconds, "ok"

    ctx = mp.get_context("fork") if hasattr(mp, "get_context") else mp
    parent, child = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=_method_worker, args=(child, alg, target, kwargs))
    t0 = time.time()
    proc.start()
    child.close()
    if parent.poll(timeout):
        payload = parent.recv()
        proc.join(timeout=1.0)
        elapsed = payload.get("seconds", time.time() - t0)
        if payload.get("ok"):
            return payload.get("result"), elapsed, "ok"
        return None, elapsed, payload.get("error", "error")
    # Timed out: terminate, then kill if needed.
    proc.terminate()
    proc.join(timeout=2.0)
    if proc.is_alive():
        try:
            proc.kill()
        except AttributeError:
            pass
        proc.join(timeout=2.0)
    return None, time.time() - t0, "timeout"

def seed_for(family: str, n: int, d: int, seed_index: int) -> int:
    base = {
        "dense064": 61000,
        "dense": 61000,
        "ks": 69000,
        "kostlan": 69000,
        "sparse_ks": 70000,
        "sparse": 70000,
        "high_ks": 71000,
    }.get(family, 72000)
    return base + 1000 * seed_index + 100 * n + d


def parse_case(s: str) -> Tuple[int, int]:
    if "," in s:
        a, b = s.split(",", 1)
    elif "x" in s:
        a, b = s.split("x", 1)
    else:
        raise argparse.ArgumentTypeError("case must look like 4,2 or 4x2")
    return int(a), int(b)


def suite_cases(name: str) -> List[Tuple[str, Tuple[int, int]]]:
    name = name.lower()
    if name == "quick":
        return [
            ("dense064", (2, 3)),
            ("dense064", (3, 2)),
            ("ks", (2, 4)),
            ("ks", (2, 6)),
            ("ks", (3, 3)),
            ("sparse_ks", (2, 8)),
        ]
    if name in {"full", "complete"}:
        return [
            ("dense064", (2, 2)),
            ("dense064", (2, 3)),
            ("dense064", (2, 4)),
            ("dense064", (2, 6)),
            ("dense064", (3, 2)),
            ("dense064", (3, 3)),
            ("dense064", (4, 2)),
            ("dense064", (5, 2)),
            ("ks", (2, 4)),
            ("ks", (2, 6)),
            ("ks", (2, 8)),
            ("ks", (2, 10)),
            ("ks", (3, 3)),
            ("ks", (3, 4)),
            ("ks", (4, 2)),
            ("ks", (5, 2)),
            ("sparse_ks", (2, 8)),
            ("sparse_ks", (2, 12)),
            ("sparse_ks", (3, 4)),
            ("sparse_ks", (4, 3)),
        ]
    if name in {"high", "high_ks"}:
        return [
            ("ks", (2, 8)),
            ("ks", (2, 10)),
            ("ks", (2, 12)),
            ("ks", (2, 14)),
            ("ks", (3, 4)),
            ("ks", (3, 5)),
            ("sparse_ks", (2, 16)),
            ("sparse_ks", (3, 5)),
        ]
    raise ValueError(f"unknown suite: {name}")


def format_pct(x: float) -> str:
    return f"{100.0 * x:5.1f}%"


def print_row(row: BenchRow) -> None:
    print(
        f"{row.family:>9} ({row.n:>2},{row.d:<2}) seed{row.seed_index:<2} "
        f"{row.alg:>12} | Bez={row.bezout:<4} terms={row.terms:<4} "
        f"roots={row.roots:<4} cov={format_pct(row.coverage)} "
        f"paths={row.paths_ok:>4}/{row.paths_fail:<4} tot={row.paths:<4} "
        f"charts={row.charts:<2} steps={row.steps:<6} work={row.epochs_or_newton:<6} "
        f"res={row.max_residual:.1e} time={row.seconds:6.2f}s"
    )


def summarize(rows: List[BenchRow]) -> None:
    print("\n" + "=" * 132)
    print("SUMMARY BY FAMILY / ALGORITHM")
    print("=" * 132)
    groups: Dict[Tuple[str, str], List[BenchRow]] = {}
    for r in rows:
        groups.setdefault((r.family, r.alg), []).append(r)
    print(f"{'family':>10} {'alg':>12} {'cases':>5} {'avg cov':>8} {'full%':>7} {'avg time':>9} {'avg paths/Bez':>13}")
    print("-" * 132)
    for (fam, alg), rs in sorted(groups.items()):
        avg_cov = statistics.mean(r.coverage for r in rs)
        full = sum(1 for r in rs if r.coverage >= 0.999999) / len(rs)
        avg_time = statistics.mean(r.seconds for r in rs)
        avg_pbez = statistics.mean(r.paths / max(1, r.bezout) for r in rs)
        print(f"{fam:>10} {alg:>12} {len(rs):>5} {100*avg_cov:7.2f}% {100*full:6.1f}% {avg_time:8.2f}s {avg_pbez:12.2f}")

    # Head-to-head rows by exact case/seed.
    print("\n" + "=" * 132)
    print("HEAD-TO-HEAD DELTAS: 068 minus Lairez-style")
    print("=" * 132)
    by_key: Dict[Tuple[str, int, int, int], Dict[str, BenchRow]] = {}
    for r in rows:
        by_key.setdefault((r.family, r.n, r.d, r.seed_index), {})[r.alg] = r
    wins_cov = wins_time = ties_cov = compared = 0
    for key, dct in sorted(by_key.items()):
        if "068" not in dct or "lairez" not in dct:
            continue
        compared += 1
        a, l = dct["068"], dct["lairez"]
        dcov = a.coverage - l.coverage
        ratio = a.seconds / l.seconds if l.seconds > 0 else float("inf")
        if dcov > 1e-12:
            wins_cov += 1
        elif abs(dcov) <= 1e-12:
            ties_cov += 1
        if a.seconds < l.seconds:
            wins_time += 1
        fam, n, dd, si = key
        print(f"{fam:>9} ({n},{dd}) seed{si:<2}: Δcov={100*dcov:+6.1f}%  time_ratio_068/lairez={ratio:6.2f}  paths_ratio={a.paths/max(1,l.paths):5.2f}")
    if compared:
        print("-" * 132)
        print(f"coverage: 068 wins {wins_cov}/{compared}, ties {ties_cov}/{compared}; time: 068 faster {wins_time}/{compared}")


def write_csv(rows: List[BenchRow], path: str) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


def write_json(rows: List[BenchRow], path: str) -> None:
    with open(path, "w") as f:
        json.dump([asdict(r) for r in rows], f, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(description="flow/069 benchmark: 068 vs Lairez-style gamma homotopy on KS/high-degree systems")
    parser.add_argument("--suite", choices=["quick", "full", "complete", "high", "high_ks"], default="quick")
    parser.add_argument("--family", default=None, help="Override suite with one family: dense064, ks, sparse_ks")
    parser.add_argument("--cases", nargs="*", type=parse_case, default=None, help="Cases like 2,8 3,4")
    parser.add_argument("--seeds", type=int, default=1, help="Number of deterministic seeds per case")
    parser.add_argument("--only", choices=["both", "068", "lairez"], default="both")
    parser.add_argument("--passes", type=int, default=1, help="068 phase passes in identity chart")
    parser.add_argument("--chart-passes", type=int, default=1, help="068 phase passes per affine chart")
    parser.add_argument("--charts", type=int, default=1, help="068 affine recovery charts after identity")
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--path-timeout", type=float, default=0.0, help="Per-path wall-clock timeout for 068 paths; 0 disables")
    parser.add_argument("--068-max-steps", type=int, default=220, help="Max homotopy steps per 068 path")
    parser.add_argument("--lairez-max-steps", type=int, default=420)
    parser.add_argument("--lairez-newton-iters", type=int, default=12)
    parser.add_argument("--lairez-retries", type=int, default=3, help="gamma retries for the compact Lairez-style baseline")
    parser.add_argument("--method-timeout", type=float, default=0.0, help="Per algorithm/case timeout in seconds; 0 disables. In this sandbox, keep 0 because Python timed waits are unreliable.")
    parser.add_argument("--max-bezout", type=int, default=300, help="Skip cases with Bezout above this unless set to 0")
    parser.add_argument("--csv", default=None)
    parser.add_argument("--json", default=None)
    args = parser.parse_args()

    if args.family:
        if not args.cases:
            raise SystemExit("--family requires --cases")
        tasks = [(args.family, c) for c in args.cases]
    else:
        tasks = suite_cases(args.suite)

    print("=" * 132)
    print("069 -- 068 polynomial Pandrosion vs Lairez-style gamma total-degree continuation")
    print("=" * 132)
    print(f"suite={args.suite}, seeds={args.seeds}, 068 passes={args.passes}, charts={args.charts}, chart_passes={args.chart_passes}, tol={args.tol:g}")
    print("Lairez baseline = gamma total-degree homotopy + analytic Newton corrector (compact local implementation).")
    print("=" * 132)

    rows: List[BenchRow] = []
    for family, (n, d) in tasks:
        B = d ** n
        if args.max_bezout and B > args.max_bezout:
            print(f"SKIP {family:>9} ({n},{d}) Bez={B} > --max-bezout={args.max_bezout}")
            continue
        for si in range(args.seeds):
            seed = seed_for(family, n, d, si)
            target = gen_system(family, n, d, seed)
            terms = term_count(target)
            B = bezout(target)

            if args.only in {"both", "068"}:
                try:
                    kwargs068 = dict(seed=68000 + seed, passes=args.passes, charts=args.charts,
                                     chart_passes=args.chart_passes, tol=args.tol,
                                     max_steps=args.__dict__["068_max_steps"],
                                     path_timeout=args.path_timeout)
                    res068, sec068, status068 = run_method_with_timeout("068", target, kwargs068, args.method_timeout)
                    if status068 != "ok" or res068 is None:
                        row = BenchRow(family, n, d, si, seed, terms, B, "068", 0, 0.0, 0, B, B,
                                       args.charts + 1, 0, 0, float("inf"), sec068)
                        print(f"WARN 068 {family} ({n},{d}) seed{si}: {status068}")
                    else:
                        row = BenchRow(
                            family=family, n=n, d=d, seed_index=si, seed=seed,
                            terms=terms, bezout=B, alg="068",
                            roots=len(res068.get("roots", [])), coverage=res068.get("coverage", 0.0),
                            paths_ok=res068.get("paths_ok", 0), paths_fail=res068.get("paths_fail", 0),
                            paths=res068.get("paths", 0), charts=res068.get("charts_used", 0),
                            steps=res068.get("steps", 0), epochs_or_newton=res068.get("epochs", 0),
                            max_residual=res068.get("max_residual", float("inf")), seconds=sec068,
                        )
                except Exception as exc:  # keep the benchmark moving
                    row = BenchRow(family, n, d, si, seed, terms, B, "068", 0, 0.0, 0, B, B, args.charts + 1, 0, 0, float("inf"), 0.0)
                    print(f"ERROR 068 {family} ({n},{d}) seed{si}: {exc}")
                rows.append(row)
                print_row(row)

            if args.only in {"both", "lairez"}:
                try:
                    kwargs_l = dict(seed=90000 + seed, tol=args.tol,
                                    max_steps=args.lairez_max_steps,
                                    max_newton_iter=args.lairez_newton_iters, retries=args.lairez_retries)
                    resl, secl, statusl = run_method_with_timeout("lairez", target, kwargs_l, args.method_timeout)
                    if statusl != "ok" or resl is None:
                        row = BenchRow(family, n, d, si, seed, terms, B, "lairez", 0, 0.0, 0, B, B,
                                       0, 0, 0, float("inf"), secl)
                        print(f"WARN lairez {family} ({n},{d}) seed{si}: {statusl}")
                    else:
                        row = BenchRow(
                            family=family, n=n, d=d, seed_index=si, seed=seed,
                            terms=terms, bezout=B, alg="lairez",
                            roots=len(resl.roots), coverage=resl.coverage,
                            paths_ok=resl.paths_ok, paths_fail=resl.paths_fail, paths=resl.paths,
                            charts=0, steps=resl.steps, epochs_or_newton=resl.newton_iters,
                            max_residual=resl.max_residual, seconds=secl,
                        )
                except Exception as exc:
                    row = BenchRow(family, n, d, si, seed, terms, B, "lairez", 0, 0.0, 0, B, B, 0, 0, 0, float("inf"), 0.0)
                    print(f"ERROR lairez {family} ({n},{d}) seed{si}: {exc}")
                rows.append(row)
                print_row(row)

    if rows:
        summarize(rows)
        if args.csv:
            write_csv(rows, args.csv)
            print(f"\nCSV written to {args.csv}")
        if args.json:
            write_json(rows, args.json)
            print(f"JSON written to {args.json}")
    print("=" * 132)


if __name__ == "__main__":
    main()
