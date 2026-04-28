"""
FLOW 080 -- per-geometry homothety with polygon start optimization.

080 takes the useful validation shape from 079 but fixes the level at which the
homothety is applied.  The path is tracked in the original z coordinates.  Each
time 070 receives a current system H_t and generates a Pandrosion geometry, the
corrector first builds

    G_t(y) = D_t H_t(S_t y),      z = S_t y,

with S_t computed from that H_t.  Thus the scaling is attached to the generated
geometry, not imposed once globally on the whole run.

The optional polygon start optimization is adapted from 079: for bivariate
systems it samples Jensen polygons of coordinate resultants and uses the
paper-0 h(1) radius only for the total-degree start roots.  No Lairez fallback or Newton-ELS is used by this script.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import signal
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

HERE = Path(__file__).resolve().parent
CORE_PATH = HERE / "076_pandrosion_homothetic_geometry.py"

spec = importlib.util.spec_from_file_location("flow076_core_for_080", str(CORE_PATH))
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {CORE_PATH}")
f076 = importlib.util.module_from_spec(spec)
sys.modules["flow076_core_for_080"] = f076
spec.loader.exec_module(f076)

Complex = complex
Vector = List[Complex]
Matrix = List[List[Complex]]


@dataclass
class PathRow080:
    family: str
    n: int
    d: int
    seed: int
    retry: int
    idx: int
    ok_track: bool
    t: float
    residual: float
    steps: int
    epochs: int
    seconds: float
    status: str
    has_root: bool


@dataclass
class SummaryRow080:
    family: str
    n: int
    d: int
    seed: int
    terms: int
    bezout: int
    alg: str
    roots: int
    coverage: float
    path_rows: int
    candidates: int
    retries_used: int
    max_residual: float
    seconds_observed: float
    status: str
    notes: str


@dataclass
class GeometryScaleStats:
    calls: int = 0
    scaled_calls: int = 0
    min_scale: float = float("inf")
    max_scale: float = 0.0
    max_anisotropy: float = 1.0

    def observe(self, scales: Sequence[float]) -> None:
        vals = [max(1e-12, float(s)) for s in scales]
        self.calls += 1
        if any(abs(math.log(v)) > 1e-12 for v in vals):
            self.scaled_calls += 1
        self.min_scale = min(self.min_scale, min(vals))
        self.max_scale = max(self.max_scale, max(vals))
        self.max_anisotropy = max(self.max_anisotropy, max(vals) / max(1e-12, min(vals)))

    def note(self) -> str:
        if self.calls <= 0:
            return "geometry_scales=none"
        return (
            f"geometry_scale_calls={self.calls}; scaled_calls={self.scaled_calls}; "
            f"scale_range=[{self.min_scale:.4g},{self.max_scale:.4g}]; "
            f"max_anisotropy={self.max_anisotropy:.4g}"
        )


SCALE_STATS = GeometryScaleStats()
ORIGINAL_CORRECTOR = f076.m.corrector


class _PathTimeout(Exception):
    pass


def _alarm_handler(_signum, _frame):
    raise _PathTimeout()


class path_timer:
    def __init__(self, seconds: float):
        self.seconds = float(seconds)
        self.old_handler = None

    def __enter__(self):
        if self.seconds > 0:
            self.old_handler = signal.getsignal(signal.SIGALRM)
            signal.signal(signal.SIGALRM, _alarm_handler)
            signal.setitimer(signal.ITIMER_REAL, self.seconds)

    def __exit__(self, exc_type, exc, tb):
        if self.seconds > 0:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
            signal.signal(signal.SIGALRM, self.old_handler)
        return False


def parse_case(s: str) -> Tuple[int, int]:
    a, b = str(s).split(",")
    return int(a), int(b)


def root_to_json(z: Sequence[Complex]):
    return [[float(v.real), float(v.imag)] for v in z]


def unique_append(roots: List[Vector], z: Sequence[Complex], sep: float = 1e-6) -> bool:
    if not f076.m.safe_z(z):
        return False
    n = len(z)
    for r in roots:
        scale = max(1.0, max(abs(z[i]) for i in range(n)), max(abs(r[i]) for i in range(n)))
        if max(abs(z[i] - r[i]) for i in range(n)) <= sep * scale:
            return False
    roots.append(list(z))
    return True


def det_logabs_complex(A: Matrix) -> float:
    n = len(A)
    if n == 0:
        return 0.0
    M = [list(map(complex, row)) for row in A]
    logabs = 0.0
    for k in range(n):
        piv = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[piv][k]) < 1e-260 or not math.isfinite(abs(M[piv][k])):
            return -700.0
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
        pv = M[k][k]
        ap = abs(pv)
        if ap < 1e-260 or not math.isfinite(ap):
            return -700.0
        logabs += math.log(ap)
        inv = 1.0 / pv
        for j in range(k, n):
            M[k][j] *= inv
        for i in range(k + 1, n):
            fac = M[i][k]
            if fac:
                for j in range(k, n):
                    M[i][j] -= fac * M[k][j]
    return max(-700.0, min(700.0, logabs))


def sylvester_two_ascending(f: Sequence[Complex], g: Sequence[Complex]) -> Matrix:
    def trim(p: Sequence[Complex]) -> List[Complex]:
        q = list(p)
        while len(q) > 1 and abs(q[-1]) < 1e-14:
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


def coord_resultant_logabs_2d(polys, coord: int, zval: Complex, maxdeg: int) -> float:
    elim = 1 - int(coord)
    coeffs = []
    for poly in polys:
        arr = [0.0 + 0.0j for _ in range(maxdeg + 1)]
        for alpha, c in poly.items():
            ae = alpha[elim]
            af = alpha[coord]
            arr[ae] += c * (zval ** af if af else 1.0)
        coeffs.append(arr)
    return det_logabs_complex(sylvester_two_ascending(coeffs[0], coeffs[1]))


def coordinate_resultant_polygon_2d(
    polys,
    coord: int,
    n_radii: int = 44,
    n_phases: int = 6,
    log_span: float = 3.0,
) -> Dict[str, float]:
    B = f076.m.bezout(polys)
    maxdeg = max((sum(a) for p0 in polys for a in p0), default=1)
    log_rs = [-log_span + 2.0 * log_span * i / max(1, n_radii - 1) for i in range(n_radii)]
    phases = [2.0 * math.pi * k / max(1, n_phases) for k in range(max(1, n_phases))]
    vals: List[float] = []
    for lr in log_rs:
        r = math.exp(lr)
        acc = 0.0
        for ph in phases:
            z = r * complex(math.cos(ph), math.sin(ph))
            acc += coord_resultant_logabs_2d(polys, coord, z, maxdeg)
        vals.append(acc / len(phases))

    slopes = []
    for i in range(n_radii - 1):
        dr = log_rs[i + 1] - log_rs[i]
        sv = (vals[i + 1] - vals[i]) / dr if dr else 0.0
        slopes.append(max(0, min(B, int(round(sv)))))
    for i in range(1, len(slopes)):
        if slopes[i] < slopes[i - 1]:
            slopes[i] = slopes[i - 1]

    def radius_at(q: float) -> float:
        target = q * B
        for i, count in enumerate(slopes):
            if count >= target:
                return math.exp(0.5 * (log_rs[i] + log_rs[i + 1]))
        return math.exp(log_rs[-1])

    r10 = radius_at(0.10)
    r50 = radius_at(0.50)
    r90 = radius_at(0.90)
    k_raw = max(1.0, r90 / max(r10, 1e-300))
    k_dy = 1.0 if k_raw <= 1.0 else 2.0 ** round(math.log(k_raw, 2.0))
    return {"r10": r10, "r50": r50, "r90": r90, "k_raw": k_raw, "k_dyadic": k_dy, "maxdeg": float(maxdeg)}


def h1_start_radius(k: float, degree: int) -> float:
    k = max(1.0, float(k))
    degree = max(2, int(degree))
    return max(0.55, min(1.20, 1.0 - (k - 1.0) / (k * degree)))


def polygon_start_radii(polys, mode: str) -> Tuple[List[float], str]:
    n = len(polys)
    if mode == "off" or n != 2:
        return [1.0 for _ in range(n)], "polygon_start=off" if mode == "off" else "polygon_start=unsupported-n"
    stats = [coordinate_resultant_polygon_2d(polys, j) for j in range(2)]
    degs = f076.m.degrees(polys)
    radii: List[float] = []
    for j, st in enumerate(stats):
        k = st["k_dyadic"] if mode == "dyadic" else st["k_raw"]
        radii.append(h1_start_radius(k, degs[j] if j < len(degs) else int(st["maxdeg"])))
    note = ";".join(
        f"j{j}:r10={stats[j]['r10']:.4g},r50={stats[j]['r50']:.4g},"
        f"r90={stats[j]['r90']:.4g},kraw={stats[j]['k_raw']:.4g},"
        f"kdy={stats[j]['k_dyadic']:.4g},h1={radii[j]:.4g}"
        for j in range(2)
    )
    return radii, note


def geometry_scales(polys, args) -> List[float]:
    n = len(polys)
    if args.geometry_homothety == "off":
        scales = [1.0 for _ in range(n)]
    else:
        scales = f076.system_diagonal_scales(
            polys,
            min_scale=args.scale_min,
            max_scale=args.scale_max,
            strength=args.system_scale_strength,
        )
    SCALE_STATS.observe(scales)
    return scales


def install_per_geometry_homothety(args) -> None:
    def corrector_080(polys, z_init: Sequence[Complex], tol: float = 1e-9,
                      max_epochs: int = 12, quad_cap: int = 16,
                      modes: Sequence[str] = ("system", "integral_gl", "blend"),
                      rescue_modes: Sequence[str] = (),
                      deadline: float | None = None):
        scales = geometry_scales(polys, args)
        use_scaled = any(abs(math.log(max(1e-12, s))) > 1e-12 for s in scales)
        if not use_scaled and not args.equation_normalize:
            return ORIGINAL_CORRECTOR(polys, z_init, tol=tol, max_epochs=max_epochs,
                                      quad_cap=quad_cap, modes=modes,
                                      rescue_modes=rescue_modes, deadline=deadline)
        scaled_polys = f076.variable_scale_system(
            polys,
            scales,
            normalize_equations=bool(args.equation_normalize),
        )
        y_init = f076.vec_to_scaled(z_init, scales)
        y, ok_scaled, epochs = ORIGINAL_CORRECTOR(
            scaled_polys,
            y_init,
            tol=tol,
            max_epochs=max_epochs,
            quad_cap=quad_cap,
            modes=modes,
            rescue_modes=rescue_modes,
            deadline=deadline,
        )
        z = f076.vec_from_scaled(y, scales)
        rz = f076.m.residual_norm(polys, z)
        ok = bool(ok_scaled and math.isfinite(rz) and rz < 60.0 * tol)
        return z, ok, epochs

    f076.m.corrector = corrector_080


def build_retry_data(args, family: str, n: int, d: int, retry: int, start_radii: Sequence[float]):
    seed = f076.m.seed_for(family, n, d, args.seed_index)
    target = f076.m.gen_system(family, n, d, seed)
    degs = f076.m.degrees(target)
    B = math.prod(degs)
    phases = f076.m.m068.deterministic_phases(n, 0, seed + 17 * B + 10007 * retry)
    radii = list(start_radii)
    if len(radii) < n:
        radii.extend([1.0] * (n - len(radii)))
    phases = [complex(phases[j]) * (max(1e-12, float(radii[j])) ** degs[j]) for j in range(n)]
    start = f076.m.m068.phase_start_system(degs, n, phases)
    roots0 = f076.m.m068.phase_start_roots(degs, phases)
    gammas = f076.m.deterministic_gamma_vector(n, seed + 991 * B, retry)
    target_gamma = f076.m.scale_system(target, gammas)
    return seed, target, start, target_gamma, roots0, B


def track_path(args, family: str, n: int, d: int, seed: int, retry: int, idx: int,
               target, start, target_gamma, roots0, deadline: float | None) -> Tuple[PathRow080, Vector | None]:
    p0 = time.time()
    remaining = max(0.001, deadline - p0) if deadline is not None else 0.0
    guard = min(float(args.path_timeout), remaining) if args.path_timeout > 0 and deadline is not None else float(args.path_timeout)
    modes = f076.parse_modes_arg(args.modes, f076.DEFAULT_MODES)
    rescue_modes = f076.parse_modes_arg(args.rescue_modes, f076.DEFAULT_RESCUE_MODES)
    try:
        with path_timer(guard):
            tr = f076.m.track_one_070(
                target,
                start,
                target_gamma,
                roots0[idx],
                tol=args.tol,
                max_steps=args.max_steps,
                max_epochs=args.max_epochs,
                quad_cap=args.quad_cap,
                modes=modes,
                rescue_modes=rescue_modes,
                deadline=deadline,
            )
            z_polished = None
            if tr.ok or tr.t > args.accept_t:
                zp = f076.m.polish_070(
                    target,
                    tr.z,
                    args.tol,
                    args.quad_cap,
                    modes,
                    rescue_modes,
                    deadline,
                )
                if zp is not None and f076.m.residual_norm(target, zp) < args.residual_accept:
                    z_polished = list(zp)
            residual = f076.m.residual_norm(target, z_polished) if z_polished is not None else float("inf")
            row = PathRow080(
                family=family, n=n, d=d, seed=seed, retry=retry, idx=idx,
                ok_track=bool(tr.ok), t=float(tr.t), residual=float(residual),
                steps=int(tr.steps), epochs=int(tr.epochs),
                seconds=float(time.time() - p0),
                status="ok" if z_polished is not None else tr.status,
                has_root=z_polished is not None,
            )
            return row, z_polished
    except _PathTimeout:
        return PathRow080(family, n, d, seed, retry, idx, False, 0.0, float("inf"), 0, 0,
                          float(time.time() - p0), "timeout", False), None
    except BaseException as exc:
        return PathRow080(family, n, d, seed, retry, idx, False, 0.0, float("inf"), 0, 0,
                          float(time.time() - p0), f"error:{type(exc).__name__}", False), None


def retry_indices(B: int, retry: int, path_rows: Sequence[PathRow080], retry_limit: int) -> List[int]:
    if retry == 0:
        return list(range(B))
    failed = [r.idx for r in path_rows if r.retry == 0 and not r.has_root]
    order = f076.golden_order(B)
    out: List[int] = []
    seen = set()
    for idx in [*failed, *order]:
        if 0 <= idx < B and idx not in seen:
            seen.add(idx)
            out.append(idx)
        if len(out) >= retry_limit:
            break
    return out


def run_case(args, case: str) -> Tuple[List[SummaryRow080], List[PathRow080], List[Vector]]:
    family = args.family
    n, d = parse_case(case)
    seed = f076.m.seed_for(family, n, d, args.seed_index)
    target0 = f076.m.gen_system(family, n, d, seed)
    B = f076.m.bezout(target0)
    terms = f076.m.term_count(target0)
    start_radii, polygon_note = polygon_start_radii(target0, args.polygon_start)
    deadline = time.time() + args.time_budget if args.time_budget > 0 else None

    print("=" * 120, flush=True)
    print("080 -- per-geometry homothety Pandrosion", flush=True)
    print("=" * 120, flush=True)
    print(f"family={family}, case=({n},{d}), seed={seed}, terms={terms}, Bezout={B}", flush=True)
    print(
        f"geometry_homothety={args.geometry_homothety}, eq_norm={args.equation_normalize}, "
        f"polygon_start={args.polygon_start}, start_radii={f076.encode_scales(start_radii)}, "
        f"time_budget={args.time_budget:g}, path_timeout={args.path_timeout:g}; no Lairez fallback; no Newton-ELS",
        flush=True,
    )
    if polygon_note:
        print(f"polygon_note={polygon_note}", flush=True)

    roots: List[Vector] = []
    path_rows: List[PathRow080] = []
    candidates = 0
    retries_used = 0
    t0 = time.time()

    for retry in range(max(1, args.retries)):
        if deadline is not None and time.time() >= deadline:
            break
        if args.stop_at_bezout and len(roots) >= B:
            break
        retries_used = retry + 1
        seed2, target, start, target_gamma, roots0, B2 = build_retry_data(args, family, n, d, retry, start_radii)
        assert B2 == B
        indices = retry_indices(B, retry, path_rows, max(0, int(args.retry_limit)))
        for k, idx in enumerate(indices, 1):
            if deadline is not None and time.time() >= deadline:
                break
            if args.stop_at_bezout and len(roots) >= B:
                break
            row, z = track_path(args, family, n, d, seed2, retry, idx, target, start, target_gamma, roots0, deadline)
            path_rows.append(row)
            if z is not None:
                candidates += 1
                unique_append(roots, z, sep=args.cluster_sep)
            if args.progress_every and (k % args.progress_every == 0 or k == len(indices)):
                print(
                    f"retry={retry} progress={k}/{len(indices)} roots={len(roots)}/{B} "
                    f"last_idx={idx} status={row.status} sec_path={row.seconds:.2f}",
                    flush=True,
                )

    residuals = [f076.m.residual_norm(target0, z) for z in roots]
    max_res = max(residuals) if residuals else float("inf")
    sec = time.time() - t0
    budget_hit = deadline is not None and time.time() >= deadline and len(roots) < B
    status = "ok" if len(roots) >= B and max_res < args.residual_accept else ("budget" if budget_hit else "partial")
    notes = (
        f"per-geometry homothety; geometry_homothety={args.geometry_homothety}; "
        f"eq_norm={args.equation_normalize}; polygon_start={args.polygon_start}; "
        f"start_radii={f076.encode_scales(start_radii)}; {polygon_note}; "
        f"{SCALE_STATS.note()}; time_budget={args.time_budget:g}; path_timeout={args.path_timeout:g}; "
        f"no Lairez fallback; no Newton-ELS"
    )
    summary = [SummaryRow080(
        family, n, d, seed, terms, B, "080-per-geometry-homothety-startopt",
        len(roots), len(roots) / max(1, B), len(path_rows), candidates,
        retries_used, max_res, sec, status, notes,
    )]
    print("-" * 120, flush=True)
    for r in summary:
        print(
            f"{r.alg:>40} roots={r.roots}/{r.bezout} cov={100*r.coverage:.1f}% "
            f"paths={r.path_rows} maxres={r.max_residual:.2e} sec={r.seconds_observed:.2f} status={r.status}",
            flush=True,
        )
    return summary, path_rows, roots


def write_outputs(all_summary: List[SummaryRow080], all_paths: List[PathRow080],
                  roots_by_case: Dict[str, List[Vector]], args) -> None:
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(asdict(all_summary[0]).keys()))
            w.writeheader()
            for row in all_summary:
                w.writerow(asdict(row))
    if args.path_csv:
        fields = list(asdict(all_paths[0]).keys()) if all_paths else list(PathRow080("", 0, 0, 0, 0, 0, False, 0, 0, 0, 0, 0, "", False).__dict__.keys())
        with open(args.path_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for row in all_paths:
                w.writerow(asdict(row))
    if args.roots_json:
        payload = {case: [root_to_json(z) for z in roots] for case, roots in roots_by_case.items()}
        Path(args.roots_json).write_text(json.dumps(payload, indent=2))
    if args.md:
        lines = ["# 080 benchmark: per-geometry homothety\n\n"]
        lines.append("| family | case | method | Bezout | roots | coverage | paths | max residual | seconds | status |\n")
        lines.append("|---|---:|---|---:|---:|---:|---:|---:|---:|---|\n")
        for r in all_summary:
            lines.append(
                f"| {r.family} | ({r.n},{r.d}) | {r.alg} | {r.bezout} | {r.roots} | "
                f"{100*r.coverage:.1f}% | {r.path_rows} | {r.max_residual:.2e} | "
                f"{r.seconds_observed:.2f} | {r.status} |\n"
            )
        lines.append(
            "\n080 applies the homothety at the geometry level: each current system H_t is "
            "scaled before 070 generates its Pandrosion slope. The optional start radius "
            "uses the Jensen/resultant polygon from 079 only for the total-degree starts. "
            "No Lairez fallback or Newton-ELS is used.\n"
        )
        Path(args.md).write_text("".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description="080 per-geometry homothety Pandrosion")
    ap.add_argument("--family", default="ks")
    ap.add_argument("--cases", nargs="+", default=["2,8"])
    ap.add_argument("--case", default=None, help="single-case alias, e.g. 2,8")
    ap.add_argument("--seed-index", type=int, default=0)
    ap.add_argument("--retries", type=int, default=2)
    ap.add_argument("--retry-limit", type=int, default=16)
    ap.add_argument("--stop-at-bezout", action="store_true")
    ap.add_argument("--cluster-sep", type=float, default=1e-6)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--residual-accept", type=float, default=1e-7)
    ap.add_argument("--max-steps", type=int, default=120)
    ap.add_argument("--max-epochs", type=int, default=4)
    ap.add_argument("--quad-cap", type=int, default=12)
    ap.add_argument("--accept-t", type=float, default=0.90)
    ap.add_argument("--time-budget", type=float, default=20.0)
    ap.add_argument("--path-timeout", type=float, default=1.5)
    ap.add_argument("--progress-every", type=int, default=16)
    ap.add_argument("--geometry-homothety", choices=["off", "system"], default="system")
    ap.add_argument("--scale-min", type=float, default=0.25)
    ap.add_argument("--scale-max", type=float, default=4.0)
    ap.add_argument("--system-scale-strength", type=float, default=1.0)
    ap.add_argument("--equation-normalize", dest="equation_normalize", action="store_true", default=True)
    ap.add_argument("--no-equation-normalize", dest="equation_normalize", action="store_false")
    ap.add_argument("--polygon-start", choices=["off", "startopt", "dyadic"], default="startopt")
    ap.add_argument("--modes", default=",".join(f076.DEFAULT_MODES))
    ap.add_argument("--rescue-modes", default=",".join(f076.DEFAULT_RESCUE_MODES))
    ap.add_argument("--csv", default=None)
    ap.add_argument("--path-csv", default=None)
    ap.add_argument("--md", default=None)
    ap.add_argument("--roots-json", default=None)
    args = ap.parse_args()

    install_per_geometry_homothety(args)
    cases = [args.case] if args.case else list(args.cases)
    all_summary: List[SummaryRow080] = []
    all_paths: List[PathRow080] = []
    roots_by_case: Dict[str, List[Vector]] = {}
    for case in cases:
        summary, paths, roots = run_case(args, case)
        all_summary.extend(summary)
        all_paths.extend(paths)
        roots_by_case[case] = roots
    write_outputs(all_summary, all_paths, roots_by_case, args)


if __name__ == "__main__":
    main()
