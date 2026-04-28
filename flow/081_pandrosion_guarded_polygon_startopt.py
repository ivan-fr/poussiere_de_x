"""
FLOW 080 -- batch-guarded Pandrosion polygon/start optimization for homothetic dD systems.

080 keeps the corrected 076 dD scaling principle

    z = S y,      G(y) = D F(S y),

and adapts the paper-0 / pandrosion_vs_lairez.py Jensen polygon idea safely.
It builds a small set of system-generated homotheties and Pandrosion h(1)
start radii, scores them on a few sentinel paths, then runs the full homotopy
with the selected geometry.  Tracking is done in batch workers: one subprocess
handles many paths, which avoids the sandbox shutdown/path-granularity issues.
No Newton-ELS fallback is used by the Pandrosion method.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

Complex = complex
Vector = List[Complex]
Matrix = List[List[Complex]]

HERE = Path(__file__).resolve().parent
CORE_CANDIDATES = [
    HERE / "076_pandrosion_homothetic_geometry_fixed.py",
    HERE / "076_pandrosion_homothetic_geometry_complete_fixed.py",
    HERE / "076_pandrosion_homothetic_geometry.py",
]
CORE_PATH = next((p for p in CORE_CANDIDATES if p.exists()), CORE_CANDIDATES[0])
if not CORE_PATH.exists():
    raise RuntimeError("080 requires 076_pandrosion_homothetic_geometry_fixed.py next to it")

spec = importlib.util.spec_from_file_location("flow076_core_for_080", str(CORE_PATH))
f076 = importlib.util.module_from_spec(spec)
sys.modules["flow076_core_for_080"] = f076
assert spec.loader is not None
spec.loader.exec_module(f076)
m = f076.m


def parse_case(s: str) -> Tuple[int, int]:
    a, b = str(s).split(",")
    return int(a), int(b)


def parse_indices(s: str) -> List[int]:
    out: List[int] = []
    for part in str(s).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            a, b = part.split(":", 1)
            out.extend(range(int(a), int(b)))
        else:
            out.append(int(part))
    seen = set(); uniq = []
    for i in out:
        if i not in seen:
            seen.add(i); uniq.append(i)
    return uniq


def encode_indices(indices: Sequence[int]) -> str:
    return ",".join(str(int(i)) for i in indices)


def parse_floats(s: str | None, n: int, default: float = 1.0) -> List[float]:
    if not s:
        return [float(default) for _ in range(n)]
    vals: List[float] = []
    for part in str(s).replace(";", ",").split(","):
        part = part.strip()
        if part:
            vals.append(float(part))
    if len(vals) == 1 and n > 1:
        vals *= n
    if len(vals) != n:
        raise ValueError(f"expected {n} floats, got {len(vals)} from {s!r}")
    return vals


def encode_floats(vals: Sequence[float]) -> str:
    return ",".join(f"{float(v):.12g}" for v in vals)


def root_to_json(z: Sequence[Complex]):
    return [[float(v.real), float(v.imag)] for v in z]


def root_from_json(zjson) -> Vector | None:
    if zjson is None:
        return None
    return [complex(float(a), float(b)) for a, b in zjson]


def unique_append(found: List[Vector], z: Sequence[Complex] | None, sep: float = 1e-6) -> bool:
    if z is None:
        return False
    try:
        if not all(math.isfinite(v.real) and math.isfinite(v.imag) and abs(v) < 1e100 for v in z):
            return False
    except Exception:
        return False
    n = len(z)
    for r in found:
        scale = max(1.0, max(abs(z[i]) for i in range(n)), max(abs(r[i]) for i in range(n)))
        if max(abs(z[i] - r[i]) for i in range(n)) <= sep * scale:
            return False
    found.append(list(z))
    return True


def cluster_roots(rows_or_roots: Sequence, sep: float = 1e-6) -> List[Vector]:
    roots: List[Vector] = []
    for item in rows_or_roots:
        if isinstance(item, dict):
            z = root_from_json(item.get("z")) if item.get("z") is not None else None
        else:
            z = item
        unique_append(roots, z, sep)
    return roots


def golden_order(B: int) -> List[int]:
    if hasattr(f076, "golden_order"):
        return list(f076.golden_order(B))
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    out: List[int] = []; seen = set(); k = 0
    while len(out) < B and k < 10 * B + 10:
        idx = int(math.floor(B * ((k * phi) % 1.0))) if k else 0
        idx = max(0, min(B - 1, idx))
        if idx not in seen:
            seen.add(idx); out.append(idx)
        k += 1
    out.extend(i for i in range(B) if i not in seen)
    return out[:B]


# -----------------------------------------------------------------------------
# Jensen/Pandrosion polygon for bivariate coordinate resultants.
# -----------------------------------------------------------------------------


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


def coordinate_resultant_polygon_2d(polys, coord: int, n_radii: int = 36, n_phases: int = 6,
                                    log_span: float = 3.0) -> Dict[str, float]:
    B = m.bezout(polys)
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
    return {"r10": r10, "r50": r50, "r90": r90, "k_raw": k_raw, "k_dyadic": k_dy,
            "maxdeg": float(maxdeg), "count_max": float(max(slopes) if slopes else 0)}


def h1_start_radius(k: float, p: int) -> float:
    k = max(1.0, float(k))
    p = max(2, int(p))
    return max(0.55, min(1.20, 1.0 - (k - 1.0) / (k * p)))


def log_blend(a: Sequence[float], b: Sequence[float], blend: float, lo: float, hi: float) -> List[float]:
    out = []
    blend = max(0.0, min(1.0, float(blend)))
    llo, lhi = math.log(lo), math.log(hi)
    for x, y in zip(a, b):
        v = (1.0 - blend) * math.log(max(1e-12, float(x))) + blend * math.log(max(1e-12, float(y)))
        if abs(v) < 0.06:
            v = 0.0
        out.append(math.exp(max(llo, min(lhi, v))))
    return out


@dataclass
class GeometryPolicy:
    name: str
    scales: List[float]
    start_radii: List[float]
    note: str


@dataclass
class BatchResult080:
    stage: str
    retry: int
    policy: str
    indices: str
    paths: int
    candidates: int
    unique_roots: int
    max_residual: float
    seconds: float
    status: str
    path_json: str


@dataclass
class PolicyScore080:
    policy: str
    sentinel_paths: int
    candidates: int
    unique_roots: int
    mean_seconds: float
    max_residual: float
    score: float
    status: str
    note: str


@dataclass
class Summary080:
    family: str
    n: int
    d: int
    seed: int
    terms: int
    bezout: int
    alg: str
    selected_policy: str
    roots: int
    coverage: float
    path_rows: int
    candidates: int
    batches: int
    retries_used: int
    max_residual: float
    seconds_observed: float
    status: str
    notes: str


def polygon_base_2d(polys, scale_min: float, scale_max: float) -> Tuple[List[float], List[float], List[float], str]:
    n = len(polys)
    if n != 2:
        return [1.0] * n, [1.0] * n, [1.0] * n, "polygon:unsupported-n"
    degs = m.degrees(polys)
    stats = [coordinate_resultant_polygon_2d(polys, j) for j in range(2)]
    poly_scales: List[float] = []
    starts_raw: List[float] = []
    starts_dy: List[float] = []
    for j, st in enumerate(stats):
        r50 = max(1e-12, st["r50"])
        k_raw = max(1.0, st["k_raw"])
        k_dy = max(1.0, st["k_dyadic"])
        p_eff = degs[j] if j < len(degs) else max(2, int(st["maxdeg"]))
        poly_scales.append(max(scale_min, min(scale_max, float(r50))))
        starts_raw.append(h1_start_radius(k_raw, p_eff))
        starts_dy.append(h1_start_radius(k_dy, p_eff))
    note = ";".join(
        f"j{j}:r10={stats[j]['r10']:.4g},r50={stats[j]['r50']:.4g},r90={stats[j]['r90']:.4g},kraw={stats[j]['k_raw']:.4g},kdy={stats[j]['k_dyadic']:.4g},hraw={starts_raw[j]:.4g},hdy={starts_dy[j]:.4g}"
        for j in range(2)
    )
    return poly_scales, starts_raw, starts_dy, note


def build_policies(polys, args) -> Tuple[List[GeometryPolicy], str]:
    n = len(polys)
    sys0 = f076.system_diagonal_scales(polys, min_scale=args.scale_min, max_scale=args.scale_max,
                                       strength=args.system_scale_strength)
    system_scales = f076.combine_scales(sys0, [1.0] * n, args.homothety, args.scale_min, args.scale_max)
    ones = [1.0] * n
    policies: List[GeometryPolicy] = [
        GeometryPolicy("system-unit", list(system_scales), list(ones), "076 baseline: system homothety, unit starts"),
    ]
    if args.include_unit_policy:
        policies.append(GeometryPolicy("unit-unit", list(ones), list(ones), "no variable homothety, unit starts"))
    polygon_note = "polygon:off"
    if n == 2 and args.polygon:
        try:
            poly_scales, hraw, hdy, polygon_note = polygon_base_2d(polys, args.scale_min, args.scale_max)
            for alpha in args.start_blends:
                soft_raw = log_blend(ones, hraw, alpha, 0.55, 1.20)
                soft_dy = log_blend(ones, hdy, alpha, 0.55, 1.20)
                policies.append(GeometryPolicy(f"system-hraw-a{alpha:g}", list(system_scales), soft_raw,
                                               f"paper-0 h(1) raw start, log blend alpha={alpha:g}"))
                policies.append(GeometryPolicy(f"system-hdy-a{alpha:g}", list(system_scales), soft_dy,
                                               f"dyadic k h(1) start, log blend alpha={alpha:g}"))
            for beta in args.scale_blends:
                ps = log_blend(system_scales, poly_scales, beta, args.scale_min, args.scale_max)
                soft = log_blend(ones, hraw, min(0.50, max(args.start_blends or [0.5])), 0.55, 1.20)
                policies.append(GeometryPolicy(f"polygonS-b{beta:g}-hsoft", ps, soft,
                                               f"Jensen/resultant median scale log-blend beta={beta:g} + guarded h(1)"))
        except Exception as exc:
            polygon_note = f"polygon-error:{type(exc).__name__}:{exc}"
    seen = set(); uniq: List[GeometryPolicy] = []
    for p in policies:
        key = (tuple(round(x, 10) for x in p.scales), tuple(round(x, 10) for x in p.start_radii))
        if key in seen:
            continue
        seen.add(key); uniq.append(p)
    if args.max_policies > 0:
        uniq = uniq[:args.max_policies]
    return uniq, polygon_note


# -----------------------------------------------------------------------------
# Batch worker
# -----------------------------------------------------------------------------


def track_one_index(args, target, target_track, start_track, target_gamma, roots0, scales: Sequence[float], idx: int, use_scaled: bool) -> dict:
    p0 = time.time()
    z_polished = None
    try:
        tr = m.track_one_070(target_track, start_track, target_gamma, roots0[idx],
                             tol=args.tol, max_steps=args.max_steps, max_epochs=args.max_epochs,
                             quad_cap=args.quad_cap)
        if tr.ok or tr.t > args.accept_t:
            y = f076.polish_070_compat(target_track, tr.z, args.tol, args.quad_cap)
            if y is not None:
                z_polished = f076.vec_from_scaled(y, scales) if use_scaled else list(y)
                rz = m.residual_norm(target, z_polished)
                if not (math.isfinite(rz) and rz < args.residual_accept):
                    zp2 = f076.polish_070_compat(target, z_polished, args.tol, args.quad_cap)
                    if zp2 is not None:
                        z_polished = zp2
        residual = m.residual_norm(target, z_polished) if z_polished is not None else float("inf")
        return {
            "idx": int(idx), "ok": bool(tr.ok), "has_root": z_polished is not None,
            "t": float(tr.t), "residual_scaled": float(tr.residual), "residual": float(residual),
            "steps": int(tr.steps), "epochs": int(tr.epochs), "seconds": float(time.time() - p0),
            "status": "ok" if z_polished is not None and residual < args.residual_accept else "no-candidate",
            "z": root_to_json(z_polished) if z_polished is not None else None,
        }
    except BaseException as exc:
        return {"idx": int(idx), "ok": False, "has_root": False, "t": 0.0,
                "residual_scaled": float("inf"), "residual": float("inf"), "steps": 0, "epochs": 0,
                "seconds": float(time.time() - p0), "status": f"error:{type(exc).__name__}",
                "error": repr(exc), "z": None}


def batch_worker(args: argparse.Namespace) -> None:
    n, d = parse_case(args.case)
    seed = m.seed_for(args.family, n, d, args.seed_index)
    target = m.gen_system(args.family, n, d, seed)
    degs = m.degrees(target)
    B = math.prod(degs)
    indices = [i for i in parse_indices(args.indices) if 0 <= i < B]
    retry = int(args.retry)
    scales = parse_floats(args.scales, n, 1.0)
    start_radii = parse_floats(args.start_radii, n, 1.0)
    t0 = time.time()
    phases = m.m068.deterministic_phases(n, 0, seed + 17 * B + 10007 * retry)
    phases = [complex(phases[j]) * (float(start_radii[j]) ** int(degs[j])) for j in range(n)]
    start_z = m.m068.phase_start_system(degs, n, phases)
    roots0_z = m.m068.phase_start_roots(degs, phases)
    use_scaled = any(abs(math.log(max(1e-12, s))) > 1e-12 for s in scales) or bool(args.equation_normalize)
    if use_scaled:
        target_track = f076.variable_scale_system(target, scales, normalize_equations=args.equation_normalize)
        start_track = f076.variable_scale_system(start_z, scales, normalize_equations=args.equation_normalize)
        roots0 = [f076.vec_to_scaled(r, scales) for r in roots0_z]
    else:
        target_track = target
        start_track = start_z
        roots0 = roots0_z
    gammas = m.deterministic_gamma_vector(n, seed + 991 * B, retry)
    target_gamma = m.scale_system(target_track, gammas)
    rows = []
    for idx in indices:
        r = track_one_index(args, target, target_track, start_track, target_gamma, roots0, scales, idx, use_scaled)
        r.update({"family": args.family, "n": n, "d": d, "seed": seed, "retry": retry,
                  "policy": args.policy, "scales": scales, "start_radii": start_radii})
        rows.append(r)
    payload = {"family": args.family, "case": args.case, "seed_index": args.seed_index, "seed": seed,
               "retry": retry, "policy": args.policy, "indices": indices, "rows": rows,
               "scales": scales, "start_radii": start_radii, "seconds": float(time.time() - t0)}
    Path(args.out).write_text(json.dumps(payload, indent=2))
    os._exit(0)


def summarize_rows(rows: Sequence[dict], sep: float) -> Tuple[int, int, float, float, int]:
    roots = cluster_roots(rows, sep=sep)
    cand = sum(1 for r in rows if r.get("z") is not None)
    residuals = [float(r.get("residual", float("inf")) or float("inf")) for r in rows if r.get("z") is not None]
    maxres = max(residuals) if residuals else float("inf")
    mean_sec = sum(float(r.get("seconds", 0.0) or 0.0) for r in rows) / max(1, len(rows))
    return len(roots), cand, maxres, mean_sec, len(rows)


def launch_batch(args: argparse.Namespace, policy: GeometryPolicy, retry: int, indices: Sequence[int],
                 outdir: Path, stage: str) -> Tuple[BatchResult080, List[dict]]:
    outdir.mkdir(parents=True, exist_ok=True)
    idx_tag = f"{min(indices)}_{max(indices)+1}_n{len(indices)}" if indices else "empty"
    safe_pol = ''.join(ch if ch.isalnum() or ch in '-_' else '_' for ch in policy.name)
    path = outdir / f"080_{args.family}_{args.case.replace(',', 'x')}_{stage}_r{retry}_{idx_tag}_{safe_pol}.json"
    cmd = [
        sys.executable, "-S", str(Path(__file__).resolve()),
        "--batch-worker", "--family", args.family, "--case", args.case,
        "--seed-index", str(args.seed_index), "--retry", str(retry),
        "--indices", encode_indices(indices), "--policy", policy.name,
        "--scales", encode_floats(policy.scales), "--start-radii", encode_floats(policy.start_radii),
        "--out", str(path), "--tol", str(args.tol), "--max-steps", str(args.max_steps),
        "--max-epochs", str(args.max_epochs), "--quad-cap", str(args.quad_cap),
        "--accept-t", str(args.accept_t), "--residual-accept", str(args.residual_accept),
    ]
    if args.equation_normalize:
        cmd.append("--equation-normalize")
    t0 = time.time(); status = "ok"
    rows: List[dict] = []
    # Poll for the JSON artifact rather than waiting for interpreter shutdown.
    # Some sandbox workers finish numerically and write JSON but hang during final
    # cleanup.  Once the artifact exists, the result is usable; terminate the worker
    # if it has not exited promptly.
    stdout = None if args.verbose else subprocess.DEVNULL
    stderr = None if args.verbose else subprocess.DEVNULL
    proc = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
    deadline = time.time() + (float(args.batch_timeout) if args.batch_timeout and args.batch_timeout > 0 else 1e9)
    while True:
        if path.exists():
            status = "ok"
            try:
                proc.wait(timeout=0.15)
            except subprocess.TimeoutExpired:
                try:
                    proc.kill()
                except Exception:
                    pass
            break
        rc = proc.poll()
        if rc is not None:
            status = "ok" if rc == 0 else f"error:{rc}"
            break
        if time.time() >= deadline:
            status = "timeout"
            try:
                proc.kill()
            except Exception:
                pass
            break
        time.sleep(0.05)
    seconds = time.time() - t0
    if path.exists():
        try:
            data = json.loads(path.read_text())
            rows = list(data.get("rows", []))
            seconds = float(data.get("seconds", seconds))
        except Exception:
            pass
    unique, cand, maxres, _mean, paths = summarize_rows(rows, args.cluster_sep)
    br = BatchResult080(stage, retry, policy.name, encode_indices(indices), paths if paths else len(indices),
                        cand, unique, maxres, seconds, status, str(path))
    return br, rows


def score_policy(policy: GeometryPolicy, rows: Sequence[dict], sep: float) -> PolicyScore080:
    unique, cand, maxres, mean_sec, paths = summarize_rows(rows, sep)
    timeouts = sum(1 for r in rows if str(r.get("status", "")).startswith("timeout"))
    resid_pen = 0.0 if (math.isfinite(maxres) and maxres < 1e-7) else 50.0
    radius_shift = sum(abs(math.log(max(1e-12, r))) for r in policy.start_radii)
    score = 10000.0 * unique + 100.0 * cand - 10.0 * mean_sec - 500.0 * timeouts - resid_pen - 2.0 * radius_shift
    status = "ok" if cand == paths and paths > 0 and timeouts == 0 else "mixed"
    return PolicyScore080(policy.name, paths, cand, unique, mean_sec, maxres, score, status, policy.note)


def run_batches_parallel(args, policy: GeometryPolicy, retry: int, chunks: Sequence[Sequence[int]], outdir: Path, stage: str) -> Tuple[List[BatchResult080], List[dict]]:
    out_batches: List[BatchResult080] = []
    out_rows: List[dict] = []
    if args.parallel_batches <= 1 or len(chunks) <= 1:
        for ch in chunks:
            br, rows = launch_batch(args, policy, retry, ch, outdir, stage)
            out_batches.append(br); out_rows.extend(rows)
        return out_batches, out_rows
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, int(args.parallel_batches))) as ex:
        futs = [(pos, ch, ex.submit(launch_batch, args, policy, retry, ch, outdir, stage)) for pos, ch in enumerate(chunks)]
        tmp = []
        for pos, ch, fut in futs:
            br, rows = fut.result()
            tmp.append((pos, br, rows))
        tmp.sort(key=lambda x: x[0])
        for _pos, br, rows in tmp:
            out_batches.append(br); out_rows.extend(rows)
    return out_batches, out_rows


# -----------------------------------------------------------------------------
# Case orchestration
# -----------------------------------------------------------------------------


def run_case(args: argparse.Namespace, case: str) -> Tuple[List[Summary080], List[PolicyScore080], List[BatchResult080], List[Vector]]:
    args.case = case
    n, d = parse_case(case)
    seed = m.seed_for(args.family, n, d, args.seed_index)
    target = m.gen_system(args.family, n, d, seed)
    B = m.bezout(target)
    terms = m.term_count(target)
    policies, polygon_note = build_policies(target, args)
    outdir = Path(args.outdir) / f"{args.family}_{n}x{d}_seed{args.seed_index}"
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 124, flush=True)
    print("080 -- batch-guarded polygon/startopt homothetic Pandrosion", flush=True)
    print("=" * 124, flush=True)
    print(f"family={args.family}, case=({n},{d}), seed={seed}, terms={terms}, Bezout={B}", flush=True)
    print(f"policies={len(policies)}, sentinel={args.sentinel_count}, batch={args.batch_size}, parallel_batches={args.parallel_batches}, timeout={args.batch_timeout:g}s, eq_norm={args.equation_normalize}; no Newton-ELS", flush=True)
    print(f"polygon_note={polygon_note}", flush=True)
    for p in policies:
        print(f"  policy {p.name}: S={encode_floats(p.scales)} start={encode_floats(p.start_radii)}", flush=True)

    t0 = time.time()
    order = golden_order(B)
    sentinel = order[:min(B, max(1, int(args.sentinel_count)))]
    scores: List[PolicyScore080] = []
    guard_batches: List[BatchResult080] = []
    for pol in policies:
        br, rows = launch_batch(args, pol, 0, sentinel, outdir, "sentinel")
        guard_batches.append(br)
        sc = score_policy(pol, rows, args.cluster_sep)
        scores.append(sc)
        print(f"sentinel {pol.name:>24}: unique={sc.unique_roots}/{len(sentinel)} cand={sc.candidates}/{sc.sentinel_paths} mean={sc.mean_seconds:.2f}s score={sc.score:.1f} status={sc.status}", flush=True)

    scores_sorted = sorted(scores, key=lambda s: (-s.score, s.mean_seconds, s.policy))
    selected = next(p for p in policies if p.name == scores_sorted[0].policy)
    sys_score = next((s for s in scores if s.policy == "system-unit"), None)
    if sys_score is not None and (scores_sorted[0].score - sys_score.score) <= args.system_tie_margin:
        selected = next((p for p in policies if p.name == "system-unit"), selected)
    if args.force_policy:
        selected = next((p for p in policies if p.name == args.force_policy), selected)
    print(f"selected_policy={selected.name}: S={encode_floats(selected.scales)} start={encode_floats(selected.start_radii)}", flush=True)

    all_batches: List[BatchResult080] = []
    if args.count_guard_paths:
        all_batches.extend([b for b in guard_batches if b.policy == selected.name])
    all_rows: List[dict] = []
    for b in all_batches:
        try:
            all_rows.extend(json.loads(Path(b.path_json).read_text()).get("rows", []))
        except Exception:
            pass

    chunks = [list(range(i, min(B, i + args.batch_size))) for i in range(0, B, args.batch_size)]
    base_batches, base_rows = run_batches_parallel(args, selected, 0, chunks, outdir, "base")
    all_batches.extend(base_batches); all_rows.extend(base_rows)
    roots = cluster_roots(all_rows, sep=args.cluster_sep)
    candidates = sum(1 for r in all_rows if r.get("z") is not None)
    print(f"base selected={selected.name} candidates={candidates}/{len(all_rows)} roots={len(roots)}/{B}", flush=True)

    retries_used = 1
    retry_rows: List[dict] = []
    if len(roots) < B and args.retries > 1:
        # Failed base paths first, then low-discrepancy order.
        failed = [int(r.get("idx", -1)) for r in base_rows if r.get("z") is None]
        retry_order = list(dict.fromkeys([*failed, *order]))[:max(0, int(args.retry_limit))]
        retry_policies = [selected]
        # Try second-best policy if genuinely different.
        for sc in scores_sorted:
            p = next(p for p in policies if p.name == sc.policy)
            if p.name != selected.name:
                retry_policies.append(p); break
        for retry in range(1, int(args.retries)):
            retries_used = retry + 1
            for pol in retry_policies:
                if args.stop_at_bezout and len(roots) >= B:
                    break
                chunks_r = [retry_order[i:i + args.micro_batch] for i in range(0, len(retry_order), args.micro_batch)]
                for ch in chunks_r:
                    br, rows = launch_batch(args, pol, retry, ch, outdir, f"retry{retry}")
                    all_batches.append(br); retry_rows.extend(rows); all_rows.extend(rows)
                    roots = cluster_roots(all_rows, sep=args.cluster_sep)
                    candidates = sum(1 for r in all_rows if r.get("z") is not None)
                    print(f"retry={retry} policy={pol.name} indices={encode_indices(ch)} roots={len(roots)}/{B} candidates={candidates}", flush=True)
                    if args.stop_at_bezout and len(roots) >= B:
                        break
            if args.stop_at_bezout and len(roots) >= B:
                break

    roots = cluster_roots(all_rows, sep=args.cluster_sep)
    residuals = [m.residual_norm(target, z) for z in roots]
    maxres = max(residuals) if residuals else float("inf")
    candidates = sum(1 for r in all_rows if r.get("z") is not None)
    status = "ok" if len(roots) >= B and maxres < args.residual_accept else "partial"
    sec = time.time() - t0
    notes = (f"batch-guarded polygon/startopt; selected={selected.name}; S={encode_floats(selected.scales)}; "
             f"start={encode_floats(selected.start_radii)}; {polygon_note}; sentinel={encode_indices(sentinel)}; "
             f"guard_scores={[asdict(s) for s in scores]}; no Newton-ELS")
    summaries: List[Summary080] = [Summary080(args.family, n, d, seed, terms, B,
        "080-batch-guarded-polygon-startopt", selected.name, len(roots), len(roots)/max(1,B),
        len(all_rows), candidates, len(all_batches), retries_used, maxres, sec, status, notes)]

    if args.include_lairez_reference:
        ref = f076.load_lairez_reference(args.family, n, d, seed, B)
        if ref is not None:
            summaries.append(Summary080(ref.family, ref.n, ref.d, ref.seed, ref.terms, ref.bezout,
                "lairez-style-reference", "n/a", ref.roots, ref.coverage, ref.path_rows,
                ref.candidates, 0, ref.retries_used, ref.max_residual, ref.seconds_observed,
                ref.status, ref.notes))
    if args.run_lairez:
        ns = argparse.Namespace(**vars(args)); ns.case = f"{n},{d}"
        lr = f076.run_lairez_now(ns, target, seed, terms, B)
        summaries.append(Summary080(lr.family, lr.n, lr.d, lr.seed, lr.terms, lr.bezout,
            "lairez-style-run", "n/a", lr.roots, lr.coverage, lr.path_rows,
            lr.candidates, 0, lr.retries_used, lr.max_residual, lr.seconds_observed,
            lr.status, lr.notes))

    for r in summaries:
        print(f"{r.alg:>38} roots={r.roots}/{r.bezout} cov={100*r.coverage:.1f}% paths={r.path_rows} batches={r.batches} maxres={r.max_residual:.2e} sec={r.seconds_observed:.2f} status={r.status}", flush=True)
    return summaries, scores, all_batches, roots


def write_outputs(summaries: List[Summary080], scores: List[PolicyScore080], batches: List[BatchResult080], roots_by_case: Dict[str, List[Vector]], args) -> None:
    if args.csv and summaries:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(asdict(summaries[0]).keys()))
            w.writeheader(); [w.writerow(asdict(r)) for r in summaries]
    if args.score_csv and scores:
        with open(args.score_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(asdict(scores[0]).keys()))
            w.writeheader(); [w.writerow(asdict(r)) for r in scores]
    if args.batch_csv and batches:
        with open(args.batch_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(asdict(batches[0]).keys()))
            w.writeheader(); [w.writerow(asdict(r)) for r in batches]
    if args.roots_json:
        payload = {case: [root_to_json(z) for z in roots] for case, roots in roots_by_case.items()}
        Path(args.roots_json).write_text(json.dumps(payload, indent=2))
    if args.md and summaries:
        lines = ["# 080 benchmark: batch-guarded Pandrosion polygon/startopt\n\n"]
        lines.append("| family | case | method | policy | Bezout | roots | coverage | paths | batches | max residual | seconds | status |\n")
        lines.append("|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for r in summaries:
            lines.append(f"| {r.family} | ({r.n},{r.d}) | {r.alg} | {r.selected_policy} | {r.bezout} | {r.roots} | {100*r.coverage:.1f}% | {r.path_rows} | {r.batches} | {r.max_residual:.2e} | {r.seconds_observed:.2f} | {r.status} |\n")
        lines.append("\n`080` adapts the Pandrosion/Jensen polygon and the paper-0 `h(1)` starting point to dD systems, but applies them only after sentinel scoring. The selected geometry is then used with the corrected 076 homothetic Pandrosion tracker in batch workers. No Newton-ELS fallback is used by the Pandrosion method.\n")
        Path(args.md).write_text("".join(lines))
    sys.stdout.flush(); sys.stderr.flush()
    if os.environ.get("PANDROSION_NO_FAST_EXIT", "") != "1" and not getattr(args, "batch_worker", False):
        os._exit(0)


def main() -> None:
    ap = argparse.ArgumentParser(description="080 batch-guarded polygon/startopt homothetic Pandrosion")
    ap.add_argument("--family", default="ks")
    ap.add_argument("--cases", nargs="+", default=["2,8", "2,10"])
    ap.add_argument("--seed-index", type=int, default=0)
    ap.add_argument("--homothety", choices=["none", "system", "roots", "hybrid"], default="system")
    ap.add_argument("--equation-normalize", action="store_true")
    ap.add_argument("--scale-min", type=float, default=0.25)
    ap.add_argument("--scale-max", type=float, default=4.0)
    ap.add_argument("--system-scale-strength", type=float, default=1.0)
    ap.add_argument("--polygon", action="store_true", default=True)
    ap.add_argument("--no-polygon", dest="polygon", action="store_false")
    ap.add_argument("--include-unit-policy", action="store_true", default=False)
    ap.add_argument("--start-blends", type=float, nargs="*", default=[0.35, 0.65])
    ap.add_argument("--scale-blends", type=float, nargs="*", default=[0.25, 0.50])
    ap.add_argument("--max-policies", type=int, default=6)
    ap.add_argument("--sentinel-count", type=int, default=6)
    ap.add_argument("--system-tie-margin", type=float, default=5.0)
    ap.add_argument("--force-policy", default="")
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--parallel-batches", type=int, default=4)
    ap.add_argument("--batch-timeout", type=float, default=90.0)
    ap.add_argument("--retries", type=int, default=2)
    ap.add_argument("--retry-limit", type=int, default=24)
    ap.add_argument("--micro-batch", type=int, default=2)
    ap.add_argument("--count-guard-paths", action="store_true", default=False)
    ap.add_argument("--stop-at-bezout", action="store_true")
    ap.add_argument("--cluster-sep", type=float, default=1e-6)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--residual-accept", type=float, default=1e-7)
    ap.add_argument("--max-steps", type=int, default=120)
    ap.add_argument("--max-epochs", type=int, default=4)
    ap.add_argument("--quad-cap", type=int, default=12)
    ap.add_argument("--accept-t", type=float, default=0.90)
    ap.add_argument("--outdir", default="080_batches")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--score-csv", default=None)
    ap.add_argument("--batch-csv", default=None)
    ap.add_argument("--md", default=None)
    ap.add_argument("--roots-json", default=None)
    ap.add_argument("--include-lairez-reference", action="store_true")
    ap.add_argument("--run-lairez", action="store_true")
    ap.add_argument("--lairez-max-steps", type=int, default=420)
    ap.add_argument("--lairez-newton-iters", type=int, default=12)
    ap.add_argument("--lairez-retries", type=int, default=2)
    ap.add_argument("--verbose", action="store_true")
    # worker-only
    ap.add_argument("--batch-worker", action="store_true")
    ap.add_argument("--case", default="")
    ap.add_argument("--indices", default="")
    ap.add_argument("--idx", type=int, default=0)  # ignored, for compatibility
    ap.add_argument("--retry", type=int, default=0)
    ap.add_argument("--policy", default="")
    ap.add_argument("--scales", default="")
    ap.add_argument("--start-radii", default="")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    args.batch_size = max(1, int(args.batch_size))
    args.parallel_batches = max(1, int(args.parallel_batches))
    args.micro_batch = max(1, int(args.micro_batch))
    args.scale_min = max(1e-6, float(args.scale_min))
    args.scale_max = max(args.scale_min, float(args.scale_max))
    if args.batch_worker:
        if not args.out:
            raise SystemExit("--batch-worker requires --out")
        batch_worker(args)
        return
    all_summary: List[Summary080] = []
    all_scores: List[PolicyScore080] = []
    all_batches: List[BatchResult080] = []
    roots_by_case: Dict[str, List[Vector]] = {}
    for case in args.cases:
        summaries, scores, batches, roots = run_case(args, case)
        all_summary.extend(summaries); all_scores.extend(scores); all_batches.extend(batches); roots_by_case[case] = roots
    write_outputs(all_summary, all_scores, all_batches, roots_by_case, args)
    sys.stdout.flush(); sys.stderr.flush()
    if os.environ.get("PANDROSION_NO_FAST_EXIT", "") != "1":
        os._exit(0)


if __name__ == "__main__":
    main()
