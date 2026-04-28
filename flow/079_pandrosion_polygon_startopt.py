"""
FLOW 079 -- polygon-scaled start optimized Pandrosion benchmark for KS systems.

079 keeps the corrected 076 scaling principle and adds an adapted paper-0 Pandrosion polygon / starting-point optimization

    z = S y,        G(y) = D F(S y),

where S is generated from the polynomial system and D optionally normalizes the
equations.  The local tracker/corrector is still the dD polynomial Pandrosion
tracker from 070/076; no Newton-ELS fallback is used.

The difference from 076/077 orchestration is pragmatic: this file can run in a
single Python process with an optional per-path SIGALRM guard.  It avoids the
sandbox failure mode where isolated workers write their JSON output but the
parent remains blocked during interpreter shutdown.  It is therefore meant as a
stable numerical validation script for KS(2,8), KS(2,10), etc.  The new polygon mode samples Jensen polygons of coordinate resultants, chooses a dyadic shell width k as in pandrosion_vs_lairez.py, and uses the paper-0 h(1) start principle to set non-unit start radii without forcing the reduced ratio to 1.
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
CORE_CANDIDATES = [
    HERE / "076_pandrosion_homothetic_geometry_fixed.py",
    HERE / "076_pandrosion_homothetic_geometry_complete_fixed.py",
    HERE / "076_pandrosion_homothetic_geometry.py",
]
CORE_PATH = next((p for p in CORE_CANDIDATES if p.exists()), CORE_CANDIDATES[0])
if not CORE_PATH.exists():
    raise RuntimeError("079 requires 076_pandrosion_homothetic_geometry_fixed.py next to it")

spec = importlib.util.spec_from_file_location("flow076_core_for_079_direct", str(CORE_PATH))
f076 = importlib.util.module_from_spec(spec)
sys.modules["flow076_core_for_079_direct"] = f076
assert spec.loader is not None
spec.loader.exec_module(f076)

Complex = complex
Vector = List[Complex]


@dataclass
class PathRow079:
    family: str
    n: int
    d: int
    seed: int
    retry: int
    idx: int
    ok_track: bool
    t: float
    residual_scaled: float
    residual: float
    steps: int
    epochs: int
    seconds: float
    status: str
    has_root: bool


@dataclass
class SummaryRow079:
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


def parse_case(s: str) -> Tuple[int, int]:
    a, b = str(s).split(",")
    return int(a), int(b)


def root_to_json(z: Sequence[Complex]):
    return [[float(v.real), float(v.imag)] for v in z]


def unique_append(roots: List[Vector], z: Sequence[Complex], sep: float = 1e-6) -> bool:
    n = len(z)
    for r in roots:
        scale = max(1.0, max(abs(z[i]) for i in range(n)), max(abs(r[i]) for i in range(n)))
        if max(abs(z[i] - r[i]) for i in range(n)) <= sep * scale:
            return False
    roots.append(list(z))
    return True


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


Matrix = List[List[Complex]]


def _det_logabs_complex(A: Matrix) -> float:
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
            f = M[i][k]
            if f:
                for j in range(k, n):
                    M[i][j] -= f * M[k][j]
    return max(-700.0, min(700.0, logabs))


def _sylvester_two_ascending(f: Sequence[Complex], g: Sequence[Complex]) -> Matrix:
    def trim(p):
        q = list(p)
        while len(q) > 1 and abs(q[-1]) < 1e-14:
            q.pop()
        return q
    f = trim(f); g = trim(g)
    mdeg = len(f) - 1; ndeg = len(g) - 1
    size = mdeg + ndeg
    if size <= 0:
        return [[1.0 + 0.0j]]
    M = [[0.0 + 0.0j for _ in range(size)] for __ in range(size)]
    fd = list(reversed(f)); gd = list(reversed(g))
    for row in range(ndeg):
        for j, val in enumerate(fd):
            M[row][row + j] = val
    for row in range(mdeg):
        for j, val in enumerate(gd):
            M[ndeg + row][row + j] = val
    return M


def _coord_resultant_logabs_2d(polys, coord: int, zval: Complex, maxdeg: int) -> float:
    elim = 1 - int(coord)
    coeffs = []
    for poly in polys:
        arr = [0.0 + 0.0j for _ in range(maxdeg + 1)]
        for alpha, c in poly.items():
            ae = alpha[elim]
            af = alpha[coord]
            arr[ae] += c * (zval ** af if af else 1.0)
        coeffs.append(arr)
    return _det_logabs_complex(_sylvester_two_ascending(coeffs[0], coeffs[1]))


def coordinate_resultant_polygon_2d(polys, coord: int, n_radii: int = 44, n_phases: int = 6,
                                    log_span: float = 3.0) -> Dict[str, float]:
    """Jensen/Pandrosion polygon for one coordinate of a bivariate system.

    This adapts pandrosion_vs_lairez.py:k_optimal_dyadic to systems: eliminate
    the other variable by a Sylvester resultant and sample Jensen averages of
    log|R_j(r e^{i phi})|.  It reads the coordinate root-shell geometry without
    solving the system.
    """
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
            acc += _coord_resultant_logabs_2d(polys, coord, z, maxdeg)
        vals.append(acc / len(phases))
    slopes = []
    for i in range(n_radii - 1):
        dr = log_rs[i + 1] - log_rs[i]
        sv = (vals[i + 1] - vals[i]) / dr if dr else 0.0
        slopes.append(max(0, min(B, int(round(sv)))))
    # Jensen cumulative root count is monotone; enforce it against sampling noise.
    for i in range(1, len(slopes)):
        if slopes[i] < slopes[i - 1]:
            slopes[i] = slopes[i - 1]
    def radius_at(q: float) -> float:
        target = q * B
        for i, c in enumerate(slopes):
            if c >= target:
                return math.exp(0.5 * (log_rs[i] + log_rs[i + 1]))
        return math.exp(log_rs[-1])
    r10, r50, r90 = radius_at(0.10), radius_at(0.50), radius_at(0.90)
    k_raw = max(1.0, r90 / max(r10, 1e-300))
    k_dy = 1.0 if k_raw <= 1.0 else 2.0 ** round(math.log(k_raw, 2.0))
    return {"r10": r10, "r50": r50, "r90": r90, "k_raw": k_raw, "k_dyadic": k_dy, "count_max": float(max(slopes) if slopes else 0), "maxdeg": float(maxdeg)}


def h1_start_radius(k: float, p: int) -> float:
    """Paper-0 h(1)=1-(x-1)/(x p), with x replaced by shell width k.

    Unlike the naive r10/r50 start, this is conservative.  It does not force the
    reduced ratio to 1; k may be 4, 5, 8, ... from the Pandrosion polygon.
    """
    k = max(1.0, float(k))
    p = max(2, int(p))
    return max(0.55, min(1.20, 1.0 - (k - 1.0) / (k * p)))


def polygon_geometry_2d(polys, scale_min: float, scale_max: float, policy: str = "startopt"):
    if len(polys) != 2:
        return [1.0 for _ in range(len(polys))], [1.0 for _ in range(len(polys))], "polygon:unsupported-n"
    stats = [coordinate_resultant_polygon_2d(polys, j) for j in range(2)]
    degs = f076.m.degrees(polys)
    scales = [1.0, 1.0]
    starts = [1.0, 1.0]
    policy = (policy or "startopt").lower()
    for j, st in enumerate(stats):
        r10 = max(1e-12, st["r10"]); r50 = max(1e-12, st["r50"]); r90 = max(1e-12, st["r90"])
        k_raw = max(1.0, st["k_raw"])
        k_dy = max(1.0, st["k_dyadic"])
        p_eff = degs[j] if j < len(degs) else max(2, int(st["maxdeg"]))
        if policy == "off":
            s, a = 1.0, 1.0
        elif policy == "median":
            s, a = r50, 1.0
        elif policy == "shell":
            # Places the lower shell around 1 and leaves width k in y-space.
            s, a = r10, h1_start_radius(k_raw, p_eff)
        elif policy == "dyadic":
            # Dyadic constructive homothety: puts a k-dyadic shell in y-space.
            s, a = r90 / k_dy, h1_start_radius(k_dy, p_eff)
        else:  # startopt default: median scale + h(1)-style anchor radius.
            s, a = r50, h1_start_radius(k_raw, p_eff)
        scales[j] = max(scale_min, min(scale_max, float(s)))
        starts[j] = max(0.55, min(1.20, float(a)))
    note = ";".join(
        f"j{j}:r10={stats[j]['r10']:.4g},r50={stats[j]['r50']:.4g},r90={stats[j]['r90']:.4g},kraw={stats[j]['k_raw']:.4g},kdy={stats[j]['k_dyadic']:.4g},h1={starts[j]:.4g}"
        for j in range(2))
    return scales, starts, note


def compute_geometry(args, target, n: int):
    system_scales = f076.system_diagonal_scales(
        target, min_scale=args.scale_min, max_scale=args.scale_max,
        strength=args.system_scale_strength,
    )
    base_scales = f076.combine_scales(system_scales, [1.0] * n, args.homothety, args.scale_min, args.scale_max)
    if args.homothety == "roots":
        base_scales = [1.0] * n
    start_radii = [1.0] * n
    polygon_note = ""
    if args.polygon_geometry != "off" and n == 2:
        try:
            poly_scales, poly_starts, polygon_note = polygon_geometry_2d(target, args.scale_min, args.scale_max, args.polygon_geometry)
            start_radii = poly_starts
            if args.polygon_scale_blend > 0:
                blend = max(0.0, min(1.0, float(args.polygon_scale_blend)))
                out = []
                for s0, sp in zip(base_scales, poly_scales):
                    x = (1.0 - blend) * math.log(max(1e-12, s0)) + blend * math.log(max(1e-12, sp))
                    out.append(max(args.scale_min, min(args.scale_max, math.exp(x))))
                base_scales = out
        except Exception as exc:
            polygon_note = f"polygon-error:{type(exc).__name__}:{exc}"
    return base_scales, start_radii, polygon_note


def compute_scales(args, target, n: int) -> List[float]:
    scales, _starts, _note = compute_geometry(args, target, n)
    return scales




def build_retry_data(family: str, n: int, d: int, seed_index: int, retry: int, scales: Sequence[float], equation_normalize: bool, start_radii: Sequence[float] | None = None):
    seed = f076.m.seed_for(family, n, d, seed_index)
    target = f076.m.gen_system(family, n, d, seed)
    degs = f076.m.degrees(target)
    B = math.prod(degs)
    phases = f076.m.m068.deterministic_phases(n, 0, seed + 17 * B + 10007 * retry)
    radii = list(start_radii or [1.0] * n)
    if len(radii) < n:
        radii += [1.0] * (n - len(radii))
    phases = [complex(phases[j]) * (float(radii[j]) ** degs[j]) for j in range(n)]
    start_z = f076.m.m068.phase_start_system(degs, n, phases)
    roots0_z = f076.m.m068.phase_start_roots(degs, phases)
    use_scaled = any(abs(math.log(max(1e-12, s))) > 1e-12 for s in scales)
    if use_scaled or equation_normalize:
        target_track = f076.variable_scale_system(target, scales, normalize_equations=equation_normalize)
        start_track = f076.variable_scale_system(start_z, scales, normalize_equations=equation_normalize)
        roots0 = [f076.vec_to_scaled(r, scales) for r in roots0_z]
    else:
        target_track = target
        start_track = start_z
        roots0 = roots0_z
    gammas = f076.m.deterministic_gamma_vector(n, seed + 991 * B, retry)
    target_gamma = f076.m.scale_system(target_track, gammas)
    return seed, target, target_track, start_track, target_gamma, roots0, B


def track_path(args, family: str, n: int, d: int, seed: int, target, target_track, start_track, target_gamma, roots0, scales, retry: int, idx: int) -> Tuple[PathRow079, Vector | None]:
    p0 = time.time()
    try:
        with path_timer(args.path_timeout):
            tr = f076.m.track_one_070(
                target_track, start_track, target_gamma, roots0[idx],
                tol=args.tol, max_steps=args.max_steps,
                max_epochs=args.max_epochs, quad_cap=args.quad_cap,
            )
            z_polished = None
            if tr.ok or tr.t > args.accept_t:
                y = f076.polish_070_compat(target_track, tr.z, args.tol, args.quad_cap)
                if y is not None:
                    z_polished = f076.vec_from_scaled(y, scales)
                    rz = f076.m.residual_norm(target, z_polished)
                    if not (math.isfinite(rz) and rz < args.residual_accept):
                        zp2 = f076.polish_070_compat(target, z_polished, args.tol, args.quad_cap)
                        if zp2 is not None:
                            z_polished = zp2
            residual = f076.m.residual_norm(target, z_polished) if z_polished is not None else float("inf")
            row = PathRow079(
                family=family, n=n, d=d, seed=seed, retry=retry, idx=idx,
                ok_track=bool(tr.ok), t=float(tr.t), residual_scaled=float(tr.residual),
                residual=float(residual), steps=int(tr.steps), epochs=int(tr.epochs),
                seconds=float(time.time() - p0), status="ok" if z_polished is not None else "no-candidate",
                has_root=z_polished is not None,
            )
            return row, list(z_polished) if z_polished is not None else None
    except _PathTimeout:
        return PathRow079(family,n,d,seed,retry,idx,False,0.0,float("inf"),float("inf"),0,0,float(time.time()-p0),"timeout",False), None
    except BaseException as exc:
        return PathRow079(family,n,d,seed,retry,idx,False,0.0,float("inf"),float("inf"),0,0,float(time.time()-p0),f"error:{type(exc).__name__}",False), None


def low_discrepancy_order(B: int) -> List[int]:
    return f076.golden_order(B)


def run_case(args, case: str) -> Tuple[List[SummaryRow079], List[PathRow079], List[Vector]]:
    family = args.family
    n, d = parse_case(case)
    seed = f076.m.seed_for(family, n, d, args.seed_index)
    target0 = f076.m.gen_system(family, n, d, seed)
    B = f076.m.bezout(target0)
    terms = f076.m.term_count(target0)
    scales, start_radii, polygon_note = compute_geometry(args, target0, n)
    print("=" * 120, flush=True)
    print("079 -- direct homothetic Pandrosion validation", flush=True)
    print("=" * 120, flush=True)
    print(f"family={family}, case=({n},{d}), seed={seed}, terms={terms}, Bezout={B}", flush=True)
    print(f"homothety={args.homothety}, scales={f076.encode_scales(scales)}, start_radii={f076.encode_scales(start_radii)}, polygon={args.polygon_geometry}, blend={args.polygon_scale_blend}, eq_norm={args.equation_normalize}, path_timeout={args.path_timeout:g}; no Newton-ELS", flush=True)
    if polygon_note:
        print(f"polygon_note={polygon_note}", flush=True)

    roots: List[Vector] = []
    path_rows: List[PathRow079] = []
    candidates = 0
    retries_used = 1
    t0 = time.time()

    for retry in range(max(1, args.retries)):
        if retry == 0:
            indices = list(range(B))
        else:
            if args.stop_at_bezout and len(roots) >= B:
                break
            order = low_discrepancy_order(B)
            # Retry failed/no-candidate paths first, then golden sentinels.
            failed = [r.idx for r in path_rows if r.retry == 0 and not r.has_root]
            seen = set()
            indices = []
            for i in [*failed, *order]:
                if 0 <= i < B and i not in seen:
                    seen.add(i); indices.append(i)
                if len(indices) >= args.retry_limit:
                    break
        retries_used = retry + 1
        seed2, target, target_track, start_track, target_gamma, roots0, B2 = build_retry_data(
            family, n, d, args.seed_index, retry, scales, args.equation_normalize, start_radii)
        assert B2 == B
        for k, idx in enumerate(indices, 1):
            row, z = track_path(args, family, n, d, seed2, target, target_track, start_track, target_gamma, roots0, scales, retry, idx)
            path_rows.append(row)
            if z is not None:
                candidates += 1
                unique_append(roots, z, sep=args.cluster_sep)
            if args.progress_every and (k % args.progress_every == 0 or k == len(indices)):
                print(f"retry={retry} progress={k}/{len(indices)} roots={len(roots)}/{B} last_idx={idx} status={row.status} sec_path={row.seconds:.2f}", flush=True)
            if args.stop_at_bezout and len(roots) >= B:
                break
        if args.stop_at_bezout and len(roots) >= B:
            break

    residuals = [f076.m.residual_norm(target0, z) for z in roots]
    max_res = max(residuals) if residuals else float("inf")
    sec = time.time() - t0
    status = "ok" if len(roots) >= B and max_res < args.residual_accept else "partial"
    notes = (
        f"direct in-process validation; homothety={args.homothety}; scales={f076.encode_scales(scales)}; start_radii={f076.encode_scales(start_radii)}; polygon={args.polygon_geometry}; blend={args.polygon_scale_blend}; {polygon_note}; "
        f"eq_norm={args.equation_normalize}; path_timeout={args.path_timeout:g}; retry_limit={args.retry_limit}; "
        f"no Newton-ELS; worker-shutdown issue avoided"
    )
    summary: List[SummaryRow079] = [SummaryRow079(
        family, n, d, seed, terms, B, "079-polygon-startopt-pandrosion",
        len(roots), len(roots)/max(1,B), len(path_rows), candidates, retries_used,
        max_res, sec, status, notes,
    )]
    if args.include_lairez_reference:
        ref = f076.load_lairez_reference(family, n, d, seed, B)
        if ref is not None:
            summary.append(SummaryRow079(ref.family, ref.n, ref.d, ref.seed, ref.terms, ref.bezout, "lairez-style-reference", ref.roots, ref.coverage, ref.path_rows, ref.candidates, ref.retries_used, ref.max_residual, ref.seconds_observed, ref.status, ref.notes))
    if args.run_lairez:
        ns = argparse.Namespace(**vars(args))
        ns.case = f"{n},{d}"
        lr = f076.run_lairez_now(ns, target0, seed, terms, B)
        summary.append(SummaryRow079(lr.family, lr.n, lr.d, lr.seed, lr.terms, lr.bezout, "lairez-style-run", lr.roots, lr.coverage, lr.path_rows, lr.candidates, lr.retries_used, lr.max_residual, lr.seconds_observed, lr.status, lr.notes))

    print("-" * 120, flush=True)
    for r in summary:
        print(f"{r.alg:>34} roots={r.roots}/{r.bezout} cov={100*r.coverage:.1f}% paths={r.path_rows} maxres={r.max_residual:.2e} sec={r.seconds_observed:.2f} status={r.status}", flush=True)
    return summary, path_rows, roots


def write_outputs(all_summary: List[SummaryRow079], all_paths: List[PathRow079], roots_by_case: Dict[str, List[Vector]], args) -> None:
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(asdict(all_summary[0]).keys()))
            w.writeheader(); [w.writerow(asdict(r)) for r in all_summary]
    if args.path_csv:
        with open(args.path_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(asdict(all_paths[0]).keys()) if all_paths else list(PathRow079("",0,0,0,0,0,False,0,0,0,0,0,0,"",False).__dict__.keys()))
            w.writeheader(); [w.writerow(asdict(r)) for r in all_paths]
    if args.roots_json:
        payload = {case: [root_to_json(z) for z in roots] for case, roots in roots_by_case.items()}
        Path(args.roots_json).write_text(json.dumps(payload, indent=2))
    if args.md:
        lines = ["# 079 benchmark: polygon-scaled start-optimized Pandrosion\n\n"]
        lines.append("| family | case | method | Bezout | roots | coverage | paths | max residual | seconds | status |\n")
        lines.append("|---|---:|---|---:|---:|---:|---:|---:|---:|---|\n")
        for r in all_summary:
            lines.append(f"| {r.family} | ({r.n},{r.d}) | {r.alg} | {r.bezout} | {r.roots} | {100*r.coverage:.1f}% | {r.path_rows} | {r.max_residual:.2e} | {r.seconds_observed:.2f} | {r.status} |\n")
        lines.append("\n079 keeps the corrected 076 homothetic system geometry and adds a Jensen/Pandrosion-polygon start optimization adapted from pandrosion_vs_lairez.py:k_optimal_dyadic. It avoids forcing x' close to 1: the dyadic shell width k can be 4, 5, 8, etc., and h(1) gives conservative anchor radii. No Newton-ELS is used by the Pandrosion method.\n")
        Path(args.md).write_text("".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description="079 direct homothetic Pandrosion validation")
    ap.add_argument("--family", default="ks")
    ap.add_argument("--cases", nargs="+", default=["2,8", "2,10"])
    ap.add_argument("--seed-index", type=int, default=0)
    ap.add_argument("--retries", type=int, default=2)
    ap.add_argument("--retry-limit", type=int, default=40)
    ap.add_argument("--stop-at-bezout", action="store_true")
    ap.add_argument("--cluster-sep", type=float, default=1e-6)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--residual-accept", type=float, default=1e-7)
    ap.add_argument("--max-steps", type=int, default=120)
    ap.add_argument("--max-epochs", type=int, default=4)
    ap.add_argument("--quad-cap", type=int, default=12)
    ap.add_argument("--accept-t", type=float, default=0.90)
    ap.add_argument("--path-timeout", type=float, default=0.0, help="0 disables SIGALRM per-path timeout")
    ap.add_argument("--progress-every", type=int, default=10)
    ap.add_argument("--homothety", choices=["none", "system", "roots", "hybrid"], default="system")
    ap.add_argument("--polygon-geometry", choices=["off", "startopt", "median", "shell", "dyadic"], default="off",
                    help="adapted pandrosion_vs_lairez k_optimal/Jensen polygon and paper-0 h(1) start optimization")
    ap.add_argument("--polygon-scale-blend", type=float, default=0.0,
                    help="0 keeps 076 system S; 1 uses polygon S; intermediate log-blend")
    ap.add_argument("--scale-min", type=float, default=0.25)
    ap.add_argument("--scale-max", type=float, default=4.0)
    ap.add_argument("--system-scale-strength", type=float, default=1.0)
    ap.add_argument("--equation-normalize", action="store_true")
    ap.add_argument("--include-lairez-reference", action="store_true")
    ap.add_argument("--run-lairez", action="store_true")
    ap.add_argument("--lairez-max-steps", type=int, default=420)
    ap.add_argument("--lairez-newton-iters", type=int, default=12)
    ap.add_argument("--lairez-retries", type=int, default=2)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--path-csv", default=None)
    ap.add_argument("--md", default=None)
    ap.add_argument("--roots-json", default=None)
    args = ap.parse_args()

    all_summary: List[SummaryRow079] = []
    all_paths: List[PathRow079] = []
    roots_by_case: Dict[str, List[Vector]] = {}
    for case in args.cases:
        summary, paths, roots = run_case(args, case)
        all_summary.extend(summary)
        all_paths.extend(paths)
        roots_by_case[case] = roots
    write_outputs(all_summary, all_paths, roots_by_case, args)


if __name__ == "__main__":
    main()
