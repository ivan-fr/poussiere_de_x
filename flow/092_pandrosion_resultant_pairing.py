"""
FLOW 092 -- bivariate resultant-pair Pandrosion recovery.

092 keeps the full 091 pipeline and tightens root clustering by default.  This
matters on sparse bivariate systems where valid nearby roots were previously
merged by the generic 1e-6 cluster radius.

It also provides an opt-in bivariate recovery layer.  Instead of using one
coordinate resultant only as a fixed-coordinate start, it reconstructs both
coordinate resultants

    R_x(x) = Res_y(F_1, F_2),   R_y(y) = Res_x(F_1, F_2),

clusters their roots, builds candidate pairs from:

  * roots of F_i(x, y) after fixing one coordinate,
  * low-residual pairs from R_x x R_y,

and accepts a candidate only after Pandrosion polish on the original system.

No Lairez/Newton-ELS fallback is used.
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
Vector = list[Complex]

HERE = Path(__file__).resolve().parent
FLOW091_PATH = HERE / "091_pandrosion_projective_riemann.py"

BEST_DEFAULT_VALUE_ARGS = {
    "--batch-timeout": "20",
    "--parallel-batches": "4",
    "--batch-size": "16",
    "--micro-batch": "4",
}
BEST_DEFAULT_FLAGS = ("--equation-normalize", "--stop-at-bezout")


def _load_091():
    spec = importlib.util.spec_from_file_location("flow091_for_092", str(FLOW091_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW091_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow091_for_092"] = mod
    spec.loader.exec_module(mod)
    return mod


def _consume_092_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    pair_seen = False
    if "--cluster-sep" in argv:
        os.environ["PANDROSION_092_USER_CLUSTER_SEP"] = "1"
    while i < len(argv):
        arg = argv[i]
        if arg == "--pair-radius" and i + 1 < len(argv):
            os.environ["PANDROSION_092_PAIR_RADIUS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--pair-samples" and i + 1 < len(argv):
            os.environ["PANDROSION_092_PAIR_SAMPLES"] = argv[i + 1]
            i += 2
            continue
        if arg == "--pair-top" and i + 1 < len(argv):
            os.environ["PANDROSION_092_PAIR_TOP"] = argv[i + 1]
            i += 2
            continue
        if arg == "--pair-budget" and i + 1 < len(argv):
            os.environ["PANDROSION_092_PAIR_BUDGET"] = argv[i + 1]
            i += 2
            continue
        if arg == "--pair-max-abs" and i + 1 < len(argv):
            os.environ["PANDROSION_092_PAIR_MAX_ABS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--no-pair-recovery":
            os.environ["PANDROSION_092_PAIR_OFF"] = "1"
            os.environ.pop("PANDROSION_092_PAIR_ON", None)
            pair_seen = True
            i += 1
            continue
        if arg == "--pair-recovery":
            os.environ["PANDROSION_092_PAIR_ON"] = "1"
            os.environ.pop("PANDROSION_092_PAIR_OFF", None)
            pair_seen = True
            i += 1
            continue
        if arg in {"--legacy-defaults", "--no-best-defaults"}:
            os.environ["PANDROSION_092_BEST_DEFAULTS_OFF"] = "1"
            i += 1
            continue
        if arg == "--no-equation-normalize":
            os.environ["PANDROSION_092_NO_EQUATION_NORMALIZE"] = "1"
            i += 1
            continue
        if arg == "--no-stop-at-bezout":
            os.environ["PANDROSION_092_NO_STOP_AT_BEZOUT"] = "1"
            i += 1
            continue
        out.append(arg)
        i += 1
    if not pair_seen and os.environ.get("PANDROSION_092_PAIR_OFF") != "1":
        os.environ["PANDROSION_092_PAIR_ON"] = "1"
    return out


def _has_option(argv: Sequence[str], opt: str) -> bool:
    prefix = opt + "="
    return any(arg == opt or arg.startswith(prefix) for arg in argv[1:])


def _apply_best_defaults(argv: list[str]) -> list[str]:
    if os.environ.get("PANDROSION_092_BEST_DEFAULTS_OFF") == "1":
        return argv
    out = list(argv)
    for opt, value in BEST_DEFAULT_VALUE_ARGS.items():
        if not _has_option(out, opt):
            out.extend([opt, value])
    if os.environ.get("PANDROSION_092_NO_EQUATION_NORMALIZE") != "1" and not _has_option(out, "--equation-normalize"):
        out.append("--equation-normalize")
    if os.environ.get("PANDROSION_092_NO_STOP_AT_BEZOUT") != "1" and not _has_option(out, "--stop-at-bezout"):
        out.append("--stop-at-bezout")
    return out


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(os.environ.get(name, str(default))))
    except Exception:
        return default


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def _pair_radii() -> list[float]:
    raw = os.environ.get("PANDROSION_092_PAIR_RADIUS", "")
    if raw:
        vals = []
        for part in raw.replace(";", ",").split(","):
            part = part.strip()
            if part:
                vals.append(max(1e-6, float(part)))
        if vals:
            return vals
    return [3.0, 4.0, 6.0, 8.0]


def _finite(v: Complex) -> bool:
    z = complex(v)
    return math.isfinite(z.real) and math.isfinite(z.imag)


def _finite_vec(z: Sequence[Complex]) -> bool:
    return all(_finite(v) for v in z)


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


def _durand_kerner(coeffs: Sequence[Complex], radius: float, max_iter: int = 240) -> list[Complex]:
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
        if max_delta < 1e-10:
            break
    return roots


def _cluster_scalars(vals: Sequence[Complex], sep: float = 1e-6, max_abs: float = 120.0) -> list[Complex]:
    out: list[Complex] = []
    for v in vals:
        if not _finite(v) or abs(v) > max_abs:
            continue
        dup = False
        for w in out:
            scale = max(1.0, abs(v), abs(w))
            if abs(v - w) <= sep * scale:
                dup = True
                break
        if not dup:
            out.append(complex(v))
    return out


def _coord_roots(polys, coord: int, degree: int, samples: int, radii: Sequence[float], max_abs: float) -> list[Complex]:
    vals: list[Complex] = []
    for radius in radii:
        coeffs = _resultant_coeffs(polys, coord, radius, samples, degree)
        vals.extend(_durand_kerner(coeffs, radius=max(18.0, 3.0 * radius), max_iter=240))
    return _cluster_scalars(vals, sep=2e-6, max_abs=max_abs)


def _poly_in_other_variable(poly, coord: int, fixed_value: Complex) -> list[Complex]:
    other = 1 - int(coord)
    maxdeg = max((alpha[other] for alpha in poly), default=0)
    out = [0.0 + 0.0j for _ in range(maxdeg + 1)]
    for alpha, coeff in poly.items():
        out[alpha[other]] += coeff * (fixed_value ** alpha[coord] if alpha[coord] else 1.0)
    return out


def _residual_candidates_from_fixed(mod, target, coord: int, fixed: Complex, max_abs: float) -> list[tuple[float, Vector]]:
    other = 1 - int(coord)
    candidates: list[Complex] = []
    for poly in target:
        coeffs = _poly_in_other_variable(poly, coord, fixed)
        candidates.extend(_durand_kerner(coeffs, radius=max(3.0, abs(fixed) + 3.0), max_iter=120))
    vals = _cluster_scalars(candidates, sep=1e-6, max_abs=max_abs)
    out: list[tuple[float, Vector]] = []
    for v in vals:
        z = [0.0 + 0.0j, 0.0 + 0.0j]
        z[coord] = fixed
        z[other] = v
        res = mod.m.residual_norm(target, z)
        if math.isfinite(res):
            out.append((res, z))
    out.sort(key=lambda item: item[0])
    return out[:4]


def _unique_pair(cands: list[tuple[float, Vector]], z: Sequence[Complex], score: float, sep: float) -> None:
    if not _finite_vec(z):
        return
    for _s, w in cands:
        scale = max(1.0, max(abs(v) for v in z), max(abs(v) for v in w))
        if max(abs(z[i] - w[i]) for i in range(2)) <= sep * scale:
            return
    cands.append((float(score), [complex(z[0]), complex(z[1])]))


def _unique_append_root(mod, roots: list[Vector], z: Sequence[Complex] | None, sep: float) -> bool:
    if z is None or not _finite_vec(z) or not mod.m.safe_z(z):
        return False
    for root in roots:
        scale = max(1.0, max(abs(v) for v in z), max(abs(v) for v in root))
        if max(abs(z[i] - root[i]) for i in range(len(z))) <= sep * scale:
            return False
    roots.append(list(z))
    return True


def _pair_recovery(mod, target, roots: list[Vector], args) -> tuple[int, int, float]:
    if os.environ.get("PANDROSION_092_PAIR_OFF") == "1" or len(target) != 2:
        return 0, 0, 0.0
    m = mod.m
    degree = m.bezout(target)
    if len(roots) >= degree:
        return 0, 0, 0.0

    t0 = time.time()
    budget = max(0.0, _env_float("PANDROSION_092_PAIR_BUDGET", 14.0))
    samples = max(128, _env_int("PANDROSION_092_PAIR_SAMPLES", max(256, 2 * degree + 64)))
    top = max(1, _env_int("PANDROSION_092_PAIR_TOP", 640))
    max_abs = max(2.0, _env_float("PANDROSION_092_PAIR_MAX_ABS", 120.0))
    radii = _pair_radii()
    sep = float(getattr(args, "cluster_sep", 1e-6))

    valid_roots = [z for z in roots if m.residual_norm(target, z) < args.residual_accept]
    roots[:] = valid_roots

    xs = _coord_roots(target, 0, degree, samples, radii, max_abs)
    ys = _coord_roots(target, 1, degree, samples, radii, max_abs)
    candidates: list[tuple[float, Vector]] = []

    for x in xs:
        for score, z in _residual_candidates_from_fixed(mod, target, 0, x, max_abs):
            _unique_pair(candidates, z, score, sep=2e-6)
    for y in ys:
        for score, z in _residual_candidates_from_fixed(mod, target, 1, y, max_abs):
            _unique_pair(candidates, z, score, sep=2e-6)

    # Cross-match both coordinate resultants.  This is cheap for B <= a few
    # hundred and catches cases where fixing only one coordinate chooses the
    # wrong branch of F_i(fixed, other).
    cross_cap = min(len(xs) * len(ys), max(top * 4, degree * 16))
    scored_cross: list[tuple[float, Vector]] = []
    for x in xs:
        if len(scored_cross) >= cross_cap:
            break
        local: list[tuple[float, Vector]] = []
        for y in ys:
            z = [x, y]
            res = m.residual_norm(target, z)
            if math.isfinite(res):
                local.append((res, z))
        local.sort(key=lambda item: item[0])
        for item in local[:4]:
            _unique_pair(scored_cross, item[1], item[0], sep=2e-6)
    for score, z in scored_cross:
        _unique_pair(candidates, z, score, sep=2e-6)

    candidates.sort(key=lambda item: item[0])
    # The resultant reconstruction is the candidate generator.  The user-facing
    # budget controls the potentially unbounded polish phase, otherwise large
    # bivariate cases can spend all their budget before trying any candidate.
    deadline = time.time() + budget if budget > 0.0 else 1e100
    attempts = hits = 0
    for _score, z in candidates[:top]:
        if time.time() >= deadline:
            break
        attempts += 1
        polished = mod.f076.polish_070_compat(target, z, args.tol, max(20, int(args.quad_cap)))
        if polished is not None and m.residual_norm(target, polished) < args.residual_accept:
            if _unique_append_root(mod, roots, polished, sep):
                hits += 1
                if len(roots) >= degree:
                    break
    return attempts, hits, time.time() - t0


def install_pair_recovery(mod) -> None:
    original_run_case = mod.run_case

    def run_case_092(args, case: str):
        n, d = mod.parse_case(case)
        old_sep = getattr(args, "cluster_sep", 1e-6)
        auto_tight = os.environ.get("PANDROSION_092_USER_CLUSTER_SEP") != "1" and n == 2
        if auto_tight:
            args.cluster_sep = float(os.environ.get("PANDROSION_092_CLUSTER_SEP", "1e-10"))
        summaries, scores, batches, roots = original_run_case(args, case)
        seed = mod.m.seed_for(args.family, n, d, args.seed_index)
        target = mod.m.gen_system(args.family, n, d, seed)
        B = mod.m.bezout(target)
        attempts, hits, sec = _pair_recovery(mod, target, roots, args) if len(roots) < B else (0, 0, 0.0)

        if summaries:
            primary = summaries[0]
            residuals = [mod.m.residual_norm(target, z) for z in roots]
            maxres = max(residuals) if residuals else float("inf")
            status = "ok" if len(roots) >= B and maxres < args.residual_accept else "partial"
            summaries[0] = replace(
                primary,
                alg="092-pandrosion-resultant-pairing",
                roots=len(roots),
                coverage=len(roots) / max(1, B),
                path_rows=primary.path_rows + attempts,
                candidates=primary.candidates + hits,
                max_residual=maxres,
                seconds_observed=primary.seconds_observed + sec,
                status=status,
                notes=(
                    f"{primary.notes}; 092_pair_attempts={attempts}; 092_pair_hits={hits}; "
                    f"092_pair_seconds={sec:.2f}; coordinate-resultant pair matching"
                ),
            )
            print(f"092-pair-recovery attempts={attempts} hits={hits} roots={len(roots)}/{B} sec={sec:.2f}", flush=True)
            print(
                f"    092-pandrosion-resultant-pairing roots={len(roots)}/{B} "
                f"cov={100.0 * len(roots) / max(1, B):.1f}% paths={primary.path_rows + attempts} "
                f"maxres={maxres:.2e} sec={primary.seconds_observed + sec:.2f} status={status}",
                flush=True,
            )
        if auto_tight:
            args.cluster_sep = old_sep
        return summaries, scores, batches, roots

    mod.run_case = run_case_092


def main() -> None:
    sys.argv = _consume_092_args(sys.argv)
    sys.argv = _apply_best_defaults(sys.argv)
    f091 = _load_091()
    sys.argv = f091._consume_091_args(sys.argv)
    f090 = f091._load_090()
    sys.argv = f090._consume_090_args(sys.argv)
    f089 = f090._load_089()
    sys.argv = f089._consume_089_args(sys.argv)
    w087 = f089._load_087()
    sys.argv = w087._consume_branch_args(sys.argv)
    mod = w087.load_081()
    w087.install_branch_safe_tracker(mod)
    w087.install_compat_and_single_policy(mod)
    f089.install_resultant_recovery(mod)
    f090.install_multivariate_completion(mod)
    f091.install_projective_completion(mod)
    install_pair_recovery(mod)
    sys.argv = f090._rewrite_argv_for_090(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
