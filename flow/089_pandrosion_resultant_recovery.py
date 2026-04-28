"""
FLOW 089 -- branch-safe Pandrosion with coordinate-resultant recovery.

089 uses the best branch-safe tracker policy found in 087
(`polygonS-b0.5-hsoft`) and then performs a small deterministic bivariate
coordinate-resultant recovery on the target system:

    R_1(z_1) = Res_{z_2}(F_1,F_2),   R_2(z_2) = Res_{z_1}(F_1,F_2).

The resultants are sampled on a circle and reconstructed by a DFT coefficient
formula.  Their roots are used only as Pandrosion polish starts.  No Lairez,
Newton-ELS, external dependency, or extra unknown is used.
"""
from __future__ import annotations

import cmath
import importlib.util
import math
import os
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Sequence

Complex = complex
Matrix = list[list[Complex]]

HERE = Path(__file__).resolve().parent
FLOW087_PATH = HERE / "087_pandrosion_branch_safe_tracker.py"


def _load_087():
    spec = importlib.util.spec_from_file_location("flow087_for_089", str(FLOW087_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW087_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow087_for_089"] = mod
    spec.loader.exec_module(mod)
    return mod


def _consume_089_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--resultant-radius" and i + 1 < len(argv):
            os.environ["PANDROSION_089_RESULTANT_RADIUS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--resultant-samples" and i + 1 < len(argv):
            os.environ["PANDROSION_089_RESULTANT_SAMPLES"] = argv[i + 1]
            i += 2
            continue
        if arg == "--no-resultant-recovery":
            os.environ["PANDROSION_089_RESULTANT_OFF"] = "1"
            i += 1
            continue
        out.append(arg)
        i += 1
    return out


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(os.environ.get(name, str(default))))
    except Exception:
        return default


def _resultant_radii() -> list[float]:
    raw = os.environ.get("PANDROSION_089_RESULTANT_RADIUS", "")
    if raw:
        vals = []
        for part in raw.replace(";", ",").split(","):
            part = part.strip()
            if part:
                vals.append(max(1e-6, float(part)))
        if vals:
            return vals
    return [4.0, 6.0]


def _det_complex(A: Matrix) -> Complex:
    n = len(A)
    M = [list(row) for row in A]
    det = 1.0 + 0.0j
    for k in range(n):
        piv = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[piv][k]) < 1e-220:
            return 0.0 + 0.0j
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
            det = -det
        pv = M[k][k]
        det *= pv
        for i in range(k + 1, n):
            factor = M[i][k] / pv
            if factor:
                for j in range(k + 1, n):
                    M[i][j] -= factor * M[k][j]
    return det


def _sylvester_two_ascending(f: Sequence[Complex], g: Sequence[Complex]) -> Matrix:
    def trim(p: Sequence[Complex]) -> list[Complex]:
        q = list(p)
        while len(q) > 1 and abs(q[-1]) < 1e-13:
            q.pop()
        return q

    f = trim(f)
    g = trim(g)
    mdeg = len(f) - 1
    ndeg = len(g) - 1
    size = mdeg + ndeg
    if size <= 0:
        return [[1.0 + 0.0j]]
    M = [[0.0 + 0.0j for _ in range(size)] for __ in range(size)]
    fd = list(reversed(f))
    gd = list(reversed(g))
    for row in range(ndeg):
        for j, val in enumerate(fd):
            M[row][row + j] = val
    for row in range(mdeg):
        for j, val in enumerate(gd):
            M[ndeg + row][row + j] = val
    return M


def _resultant_value(polys, coord: int, zval: Complex) -> Complex:
    elim = 1 - int(coord)
    coeffs = []
    for poly in polys:
        maxdeg = max((alpha[elim] for alpha in poly), default=0)
        arr = [0.0 + 0.0j for _ in range(maxdeg + 1)]
        for alpha, coeff in poly.items():
            arr[alpha[elim]] += coeff * (zval ** alpha[coord] if alpha[coord] else 1.0)
        coeffs.append(arr)
    return _det_complex(_sylvester_two_ascending(coeffs[0], coeffs[1]))


def _resultant_coeffs(polys, coord: int, radius: float, samples: int, degree: int) -> list[Complex]:
    values = [
        _resultant_value(polys, coord, radius * cmath.exp(2j * math.pi * j / samples))
        for j in range(samples)
    ]
    coeffs = []
    for k in range(degree + 1):
        s = 0.0 + 0.0j
        for j, value in enumerate(values):
            s += value * cmath.exp(-2j * math.pi * j * k / samples)
        coeffs.append(s / samples / (radius ** k))
    return coeffs


def _poly_eval(coeffs: Sequence[Complex], z: Complex) -> Complex:
    out = 0.0 + 0.0j
    for coeff in reversed(coeffs):
        out = out * z + coeff
    return out


def _durand_kerner(coeffs: Sequence[Complex], radius: float, max_iter: int = 220) -> list[Complex]:
    n = len(coeffs) - 1
    if n <= 0 or abs(coeffs[-1]) == 0:
        return []
    mono = [c / coeffs[-1] for c in coeffs]
    roots = [radius * cmath.exp(2j * math.pi * (i + 0.37) / n) for i in range(n)]
    for _ in range(max_iter):
        max_delta = 0.0
        for i, z in enumerate(roots):
            den = 1.0 + 0.0j
            for j, w in enumerate(roots):
                if i != j:
                    den *= z - w if z != w else 1e-12
            delta = _poly_eval(mono, z) / den
            if abs(delta) > 10.0 * radius:
                delta = delta / abs(delta) * (10.0 * radius)
            roots[i] = z - delta
            max_delta = max(max_delta, abs(delta))
        if max_delta < 1e-9:
            break
    return roots


def _poly_in_other_variable(poly, coord: int, fixed_value: Complex) -> list[Complex]:
    other = 1 - int(coord)
    maxdeg = max((alpha[other] for alpha in poly), default=0)
    out = [0.0 + 0.0j for _ in range(maxdeg + 1)]
    for alpha, coeff in poly.items():
        out[alpha[other]] += coeff * (fixed_value ** alpha[coord] if alpha[coord] else 1.0)
    return out


def _unique_append(m, roots: list[list[Complex]], z: Sequence[Complex], sep: float) -> bool:
    if z is None or not m.safe_z(z):
        return False
    for root in roots:
        scale = max(1.0, max(abs(v) for v in z), max(abs(v) for v in root))
        if max(abs(z[i] - root[i]) for i in range(len(z))) <= sep * scale:
            return False
    roots.append(list(z))
    return True


def _resultant_recovery(mod, target, roots: list[list[Complex]], args) -> tuple[int, int, float]:
    if os.environ.get("PANDROSION_089_RESULTANT_OFF") == "1" or len(target) != 2:
        return 0, 0, 0.0
    m = mod.m
    f076 = mod.f076
    radii = _resultant_radii()
    samples = max(128, _env_int("PANDROSION_089_RESULTANT_SAMPLES", 256))
    degree = m.bezout(target)
    attempts = 0
    hits = 0
    sep = float(getattr(args, "cluster_sep", 1e-6))
    t0 = time.time()

    # Keep the final accounting strict: roots with residual above the acceptance
    # threshold are not counted as solved roots for the resultant phase.
    valid_roots = [z for z in roots if m.residual_norm(target, z) < args.residual_accept]
    roots[:] = valid_roots

    for radius in radii:
        for coord in (0, 1):
            coeffs = _resultant_coeffs(target, coord, radius, samples, degree)
            candidates = _durand_kerner(coeffs, radius=max(18.0, 3.0 * radius), max_iter=220)
            for fixed in candidates:
                if not (math.isfinite(fixed.real) and math.isfinite(fixed.imag) and abs(fixed) < 80.0):
                    continue
                other_poly = _poly_in_other_variable(target[0], coord, fixed)
                other_roots = _durand_kerner(other_poly, radius=max(3.0, abs(fixed) + 3.0), max_iter=100)
                best = None
                best_res = float("inf")
                for other in other_roots:
                    z = [0.0 + 0.0j, 0.0 + 0.0j]
                    z[coord] = fixed
                    z[1 - coord] = other
                    res = m.residual_norm(target, z)
                    if math.isfinite(res) and res < best_res:
                        best = z
                        best_res = res
                attempts += 1
                if best is None:
                    continue
                polished = m.polish_070(
                    target,
                    best,
                    args.tol,
                    max(20, int(args.quad_cap)),
                    f076.DEFAULT_MODES,
                    f076.DEFAULT_RESCUE_MODES,
                    None,
                )
                if polished is not None and m.residual_norm(target, polished) < args.residual_accept:
                    if _unique_append(m, roots, polished, sep):
                        hits += 1
                        if len(roots) >= degree:
                            return attempts, hits, time.time() - t0
    return attempts, hits, time.time() - t0


def install_resultant_recovery(mod) -> None:
    original_run_case = mod.run_case

    def run_case_089(args, case: str):
        summaries, scores, batches, roots = original_run_case(args, case)
        n, d = mod.parse_case(case)
        seed = mod.m.seed_for(args.family, n, d, args.seed_index)
        target = mod.m.gen_system(args.family, n, d, seed)
        B = mod.m.bezout(target)
        attempts, hits, sec = _resultant_recovery(mod, target, roots, args) if len(roots) < B else (0, 0, 0.0)

        if summaries:
            primary = summaries[0]
            residuals = [mod.m.residual_norm(target, z) for z in roots]
            maxres = max(residuals) if residuals else float("inf")
            status = "ok" if len(roots) >= B and maxres < args.residual_accept else "partial"
            summaries[0] = replace(
                primary,
                alg="089-branch-safe-resultant-recovery",
                roots=len(roots),
                coverage=len(roots) / max(1, B),
                path_rows=primary.path_rows + attempts,
                candidates=primary.candidates + hits,
                max_residual=maxres,
                seconds_observed=primary.seconds_observed + sec,
                status=status,
                notes=(
                    f"{primary.notes}; resultant_attempts={attempts}; resultant_hits={hits}; "
                    f"resultant_seconds={sec:.2f}; coordinate-resultant Pandrosion recovery"
                ),
            )
            print(f"resultant-recovery attempts={attempts} hits={hits} roots={len(roots)}/{B} sec={sec:.2f}", flush=True)
            print(
                f"    089-branch-safe-resultant-recovery roots={len(roots)}/{B} "
                f"cov={100.0 * len(roots) / max(1, B):.1f}% paths={primary.path_rows + attempts} "
                f"maxres={maxres:.2e} sec={primary.seconds_observed + sec:.2f} status={status}",
                flush=True,
            )
        return summaries, scores, batches, roots

    mod.run_case = run_case_089


def _rewrite_argv_for_089(w087, argv: list[str]) -> list[str]:
    out = list(argv)
    if "--batch-worker" not in out and "--outdir" not in out:
        out.extend(["--outdir", "089_batches"])
    if "--batch-worker" not in out and "--force-policy" not in out:
        out.extend(["--force-policy", "polygonS-b0.5-hsoft"])
    return w087.rewrite_argv(out)


def main() -> None:
    sys.argv = _consume_089_args(sys.argv)
    w087 = _load_087()
    sys.argv = w087._consume_branch_args(sys.argv)
    mod = w087.load_081()
    w087.install_branch_safe_tracker(mod)
    w087.install_compat_and_single_policy(mod)
    install_resultant_recovery(mod)
    sys.argv = _rewrite_argv_for_089(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
