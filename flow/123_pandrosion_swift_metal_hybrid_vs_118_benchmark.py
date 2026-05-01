#!/usr/bin/env python3
"""
123_pandrosion_swift_metal_hybrid_vs_118_benchmark.py

Next step toward a full Swift + Metal port of 118:

  1. Generate many 118 Mobius/Thales starts.
  2. Score them in one Swift + Metal batch by ||F(y)||.
  3. Refine the best candidates with the original 118 complex128 corrector.
  4. Compare against the unmodified 118 solver on the same degree > 33 case.

This is deliberately hybrid.  It keeps the final solve/validation in 118's
complex128 path while moving the broad candidate search onto Metal, which is the
part that can be batched efficiently on Apple Silicon.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
FLOW = Path(__file__).resolve().parent


def _load_118() -> Any:
    path = FLOW / "118_pandrosion_probe_aware_pure_thales_engine.py"
    spec = importlib.util.spec_from_file_location("pandrosion_118_for_123", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load 118 engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ENGINE118 = _load_118()
np = ENGINE118.np


def now() -> float:
    return time.perf_counter()


def parse_cases(raw: str) -> list[str]:
    return [c.strip() for c in str(raw).replace("|", ";").split(";") if c.strip()]


def fmt_float(value: Any, digits: int = 4) -> str:
    if value is None:
        return ""
    try:
        x = float(value)
    except Exception:
        return ""
    if not math.isfinite(x):
        return ""
    return f"{x:.{digits}f}"


def fmt_sci(value: Any) -> str:
    if value is None:
        return ""
    try:
        x = float(value)
    except Exception:
        return ""
    if not math.isfinite(x):
        return ""
    return f"{x:.3e}"


def build_swift_helper(source: Path, binary: Path, force: bool = False) -> dict[str, Any]:
    binary.parent.mkdir(parents=True, exist_ok=True)
    needs_build = force or not binary.exists() or source.stat().st_mtime > binary.stat().st_mtime
    out: dict[str, Any] = {"source": str(source), "binary": str(binary), "rebuilt": bool(needs_build)}
    if not needs_build:
        return out
    cmd = [
        "swiftc",
        "-O",
        str(source),
        "-o",
        str(binary),
        "-framework",
        "Foundation",
        "-framework",
        "Metal",
    ]
    t0 = now()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True)
    out["build_seconds"] = now() - t0
    out["build_stdout"] = proc.stdout
    out["build_stderr"] = proc.stderr
    if proc.returncode != 0:
        raise RuntimeError(f"swiftc failed with code {proc.returncode}: {proc.stderr}")
    return out


def make_engine_args(args: argparse.Namespace) -> argparse.Namespace:
    ns = ENGINE118.build_parser().parse_args([])
    for name in [
        "seed_index",
        "equation_normalize",
        "linear_scale",
        "count",
        "pool",
        "epochs",
        "tol",
        "accept",
        "cluster_sep",
        "trial_timeout",
        "line_search",
        "probe_scale",
        "probe_candidates",
        "probe_radii",
        "probe_self",
        "powers",
        "power_cap",
        "angles",
        "rays",
        "startopt_steps",
        "startopt_candidates",
        "startopt_gains",
        "startopt_micro_epochs",
        "keep_trials",
        "verbose_trials",
        "outdir",
    ]:
        if hasattr(ns, name) and hasattr(args, name):
            setattr(ns, name, getattr(args, name))
    ns.out = None
    ns.cases = str(args.cases)
    return ns


def write_selector_input(path: Path, system: Any, points: Any, trial_ids: Sequence[int]) -> None:
    coeff = np.asarray(system.coeff, dtype=np.complex64)
    pts = np.asarray(points, dtype=np.complex64)
    payload = {
        "n": int(system.n),
        "equations": int(system.n),
        "terms": int(system.terms_per_poly),
        "points": int(pts.shape[0]),
        "exps": [int(x) for x in np.asarray(system.exps, dtype=np.int32).reshape(-1).tolist()],
        "coeffRe": [float(x) for x in coeff.real.reshape(-1).tolist()],
        "coeffIm": [float(x) for x in coeff.imag.reshape(-1).tolist()],
        "pointsRe": [float(x) for x in pts.real.reshape(-1).tolist()],
        "pointsIm": [float(x) for x in pts.imag.reshape(-1).tolist()],
        "trialIds": [int(x) for x in trial_ids],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")


def write_start_selector_input(path: Path, args: argparse.Namespace, system: Any) -> None:
    powers = sorted(set(round(float(x), 16) for x in ENGINE118.parse_float_list(args.powers, ENGINE118.DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in ENGINE118.parse_float_list(args.angles, ENGINE118.DEFAULT_ANGLES_DEG)]
    radii = ENGINE118.parse_float_list(args.rays, ENGINE118.DEFAULT_RADII, positive=True)
    coeff = np.asarray(system.coeff, dtype=np.complex64)
    payload = {
        "n": int(system.n),
        "equations": int(system.n),
        "terms": int(system.terms_per_poly),
        "candidates": int(args.metal_candidates),
        "targetCount": int(args.count),
        "seed": int(system.seed) + 0x113000,
        "powerCap": float(args.power_cap),
        "exps": [int(x) for x in np.asarray(system.exps, dtype=np.int32).reshape(-1).tolist()],
        "coeffRe": [float(x) for x in coeff.real.reshape(-1).tolist()],
        "coeffIm": [float(x) for x in coeff.imag.reshape(-1).tolist()],
        "powers": [float(x) for x in powers],
        "angles": [float(x) for x in angles],
        "radii": [float(x) for x in radii],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")


def start_selector_params(args: argparse.Namespace, system: Any) -> dict[str, Any]:
    powers = sorted(set(round(float(x), 16) for x in ENGINE118.parse_float_list(args.powers, ENGINE118.DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in ENGINE118.parse_float_list(args.angles, ENGINE118.DEFAULT_ANGLES_DEG)]
    radii = ENGINE118.parse_float_list(args.rays, ENGINE118.DEFAULT_RADII, positive=True)
    return {
        "op": "start_select",
        "candidates": int(args.metal_candidates),
        "targetCount": int(args.count),
        "topK": int(args.selector_top) if int(args.selector_top) > 0 else int(args.refine_top) * 8,
        "seed": int(system.seed) + 0x113000,
        "powerCap": float(args.power_cap),
        "powers": [float(x) for x in powers],
        "angles": [float(x) for x in angles],
        "radii": [float(x) for x in radii],
    }


def load_selector_system(proc: Any, system: Any) -> None:
    if proc.stdin is None or proc.stdout is None:
        raise RuntimeError("selector server has no stdin/stdout pipe")
    coeff = np.asarray(system.coeff, dtype=np.complex64)
    payload = {
        "op": "load",
        "n": int(system.n),
        "equations": int(system.n),
        "terms": int(system.terms_per_poly),
        "exps": [int(x) for x in np.asarray(system.exps, dtype=np.int32).reshape(-1).tolist()],
        "coeffRe": [float(x) for x in coeff.real.reshape(-1).tolist()],
        "coeffIm": [float(x) for x in coeff.imag.reshape(-1).tolist()],
    }
    proc.stdin.write(json.dumps(payload, separators=(",", ":")) + "\n")
    proc.stdin.flush()
    response = proc.stdout.readline()
    if not response:
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        raise RuntimeError(f"selector server load failed: {stderr}")
    result = json.loads(response)
    if result.get("error") or not result.get("ok"):
        raise RuntimeError(f"selector server load error: {result}")


def generate_batch_starts(args: argparse.Namespace, system: Any, case_raw: str) -> tuple[Any, list[int], list[dict[str, Any]], float]:
    n = int(system.n)
    powers = sorted(set(round(float(x), 16) for x in ENGINE118.parse_float_list(args.powers, ENGINE118.DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in ENGINE118.parse_float_list(args.angles, ENGINE118.DEFAULT_ANGLES_DEG)]
    radii = ENGINE118.parse_float_list(args.rays, ENGINE118.DEFAULT_RADII, positive=True)

    points = []
    trials = []
    meta = []
    t0 = now()
    pressure_span = max(1, int(args.metal_candidates) // max(1, int(args.count)))
    for trial in range(int(args.metal_candidates)):
        if bool(args.synthetic_pressure):
            roots_found = min(max(0, int(args.count) - 1), trial // pressure_span)
            duplicates = (trial // 17 + 3 * roots_found) % max(2, 2 * int(args.count) + 3)
            failures = (trial // 29) % max(2, int(args.count) + 5)
        else:
            roots_found = 0
            duplicates = 0
            failures = 0
        y_raw, geom = ENGINE118.mobius_homothety_start(
            n,
            trial,
            int(system.seed) + 0x113000,
            powers,
            angles,
            radii,
            float(args.power_cap),
            roots_found=roots_found,
            duplicates=duplicates,
            failures=failures,
            target_count=int(args.count),
        )
        points.append(np.asarray(y_raw, dtype=np.complex128))
        trials.append(int(trial))
        meta.append(geom)
    return np.asarray(points, dtype=np.complex128), trials, meta, now() - t0


def candidate_point(cand: dict[str, Any]) -> Any:
    return np.asarray(cand["re"], dtype=np.float64) + 1j * np.asarray(cand["im"], dtype=np.float64)


def diversify_selected(selected: Sequence[dict[str, Any]], limit: int, sep: float) -> list[dict[str, Any]]:
    wanted = max(1, int(limit))
    sep = max(0.0, float(sep))
    if sep <= 0:
        return [dict(c) for c in selected[:wanted]]

    chosen: list[dict[str, Any]] = []
    chosen_points: list[Any] = []
    deferred: list[dict[str, Any]] = []
    for cand in selected:
        z = candidate_point(cand)
        zn = max(1.0, float(np.linalg.norm(z)))
        near = False
        for prev in chosen_points:
            dist = float(np.linalg.norm(z - prev)) / max(zn, float(np.linalg.norm(prev)), 1.0)
            if dist < sep:
                near = True
                break
        if near:
            deferred.append(dict(cand))
            continue
        chosen.append(dict(cand))
        chosen_points.append(z)
        if len(chosen) >= wanted:
            return chosen

    seen = {int(c.get("index", -1)) for c in chosen}
    for cand in deferred:
        idx = int(cand.get("index", -1))
        if idx in seen:
            continue
        chosen.append(dict(cand))
        seen.add(idx)
        if len(chosen) >= wanted:
            break
    return chosen


def run_selector(binary: Path, input_path: Path, output_path: Path, top_k: int) -> dict[str, Any]:
    t0 = now()
    proc = subprocess.run(
        [str(binary), str(input_path), str(max(1, int(top_k))), str(output_path)],
        cwd=str(ROOT),
        text=True,
        capture_output=True,
    )
    wall = now() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Swift Metal selector failed with code {proc.returncode}: {proc.stderr}")
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    payload["process_wall_seconds"] = float(wall)
    payload["stderr"] = proc.stderr
    return payload


def run_selector_server(proc: Any, input_path: Path, output_path: Path, top_k: int) -> dict[str, Any]:
    if proc.stdin is None or proc.stdout is None:
        raise RuntimeError("selector server has no stdin/stdout pipe")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    payload["topK"] = int(top_k)
    line = json.dumps(payload, separators=(",", ":"))
    t0 = now()
    proc.stdin.write(line + "\n")
    proc.stdin.flush()
    response = proc.stdout.readline()
    wall = now() - t0
    if not response:
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        raise RuntimeError(f"selector server stopped without response: {stderr}")
    result = json.loads(response)
    if result.get("error"):
        raise RuntimeError(f"selector server error: {result['error']}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    result["process_wall_seconds"] = float(wall)
    result["stderr"] = ""
    return result


def run_stateful_selector_server(proc: Any, payload: dict[str, Any], output_path: Path) -> dict[str, Any]:
    if proc.stdin is None or proc.stdout is None:
        raise RuntimeError("selector server has no stdin/stdout pipe")
    t0 = now()
    proc.stdin.write(json.dumps(payload, separators=(",", ":")) + "\n")
    proc.stdin.flush()
    response = proc.stdout.readline()
    wall = now() - t0
    if not response:
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        raise RuntimeError(f"selector server stopped without response: {stderr}")
    result = json.loads(response)
    if result.get("error"):
        raise RuntimeError(f"selector server error: {result['error']}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    result["process_wall_seconds"] = float(wall)
    result["stderr"] = ""
    return result


def warm_selector_server(proc: Any) -> None:
    if proc.stdin is None or proc.stdout is None:
        raise RuntimeError("selector server has no stdin/stdout pipe")
    payload = {
        "n": 1,
        "equations": 1,
        "terms": 1,
        "candidates": 1,
        "targetCount": 1,
        "topK": 1,
        "seed": 1,
        "powerCap": 1.0,
        "exps": [0],
        "coeffRe": [1.0],
        "coeffIm": [0.0],
        "powers": [1.0],
        "angles": [0.0],
        "radii": [0.1],
    }
    proc.stdin.write(json.dumps(payload, separators=(",", ":")) + "\n")
    proc.stdin.flush()
    response = proc.stdout.readline()
    if not response:
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        raise RuntimeError(f"selector server warmup failed: {stderr}")
    result = json.loads(response)
    if result.get("error"):
        raise RuntimeError(f"selector server warmup error: {result['error']}")


def residual_stats(values: Sequence[float]) -> dict[str, Optional[float]]:
    finite = [float(v) for v in values if math.isfinite(float(v))]
    if not finite:
        return {"residual_min": None, "residual_mean": None, "residual_max": None}
    return {
        "residual_min": float(min(finite)),
        "residual_mean": float(sum(finite) / len(finite)),
        "residual_max": float(max(finite)),
    }


def refine_selected(args: argparse.Namespace, system: Any, selected: Sequence[dict[str, Any]], geom_meta: Sequence[dict[str, Any]]) -> dict[str, Any]:
    n = int(system.n)
    chart = ENGINE118.LinearChart.identity(n, scale=float(args.linear_scale))
    target = ENGINE118.TargetTrack(system, chart)
    gains = ENGINE118.parse_float_list(args.startopt_gains, ENGINE118.DEFAULT_GAINS, positive=True)
    probe_radii = ENGINE118.parse_float_list(args.probe_radii, [0.0, 0.35, 0.7, 1.0, 1.6, 2.6, 4.2], positive=False)
    probe_radii = [r for r in probe_radii if r >= 0] or [0.0, 1.0]

    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    failures = 0
    duplicates = 0
    t0 = now()
    for cand in selected:
        if len(roots) >= int(args.count):
            break
        trial = int(cand.get("trial", cand.get("index", 0)))
        y_raw = np.asarray(cand["re"], dtype=np.float64) + 1j * np.asarray(cand["im"], dtype=np.float64)
        y_raw = y_raw.astype(np.complex128)
        geom = dict(geom_meta[trial]) if 0 <= trial < len(geom_meta) else {}
        y0, smeta = ENGINE118.startopt(
            target,
            y_raw,
            trial,
            int(system.seed) + 0x112555,
            int(args.startopt_steps),
            int(args.startopt_candidates),
            gains,
            int(args.startopt_micro_epochs),
        )
        loc = ENGINE118.pandrosion_corrector(
            target,
            y0,
            max_epochs=int(args.epochs),
            tol=float(args.tol),
            accept=float(args.accept),
            trial_timeout=float(args.trial_timeout),
            line_search=int(args.line_search),
            probe_scale=float(getattr(args, "probe_scale", 0.035)),
            direction_seed=int(system.seed) + 7919 * trial,
            probe_candidates=int(args.probe_candidates),
            probe_radii=probe_radii,
            include_self_probe=bool(args.probe_self),
        )
        z = chart.z_from_y(loc["y"])
        r_orig = float(np.linalg.norm(system.eval(z)))
        accepted = bool(math.isfinite(r_orig) and r_orig < float(args.accept))
        rec = {
            "trial": int(trial),
            "selector_rank": int(cand.get("rank", -1)),
            "selector_residual": float(cand.get("residual", float("inf"))),
            "accepted": accepted,
            "status": loc.get("status"),
            "r1": r_orig,
            "epochs": int(loc.get("epochs", 0)),
            "seconds": float(loc.get("seconds", 0.0)),
            "probe_total_evals": int(loc.get("probe_total_evals", 0)),
            **geom,
            **smeta,
        }
        if not accepted:
            failures += 1
            trials.append(rec)
            continue
        dup = ENGINE118.cluster_index(roots, z, float(args.cluster_sep))
        if dup is not None:
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = int(dup)
            trials.append(rec)
            continue
        rid = len(roots)
        cond = ENGINE118.slope_condition_from_corrector(loc)
        realv = ENGINE118.realness(z)
        roots.append(
            {
                "id": rid,
                "source": "123-swift-metal-selector/118-complex128-refine",
                "trial": int(trial),
                "selector_rank": int(cand.get("rank", -1)),
                "selector_residual": float(cand.get("residual", float("inf"))),
                "z_complex": np.asarray(z, dtype=np.complex128).copy(),
                "y_complex": np.asarray(loc["y"], dtype=np.complex128).copy(),
                "residual": float(r_orig),
                "realness": float(realv),
                "cond": cond,
                "score": ENGINE118.score_root(float(r_orig), realv, cond),
                "epochs": int(loc.get("epochs", 0)),
                "seconds": float(loc.get("seconds", 0.0)),
                **geom,
                **smeta,
            }
        )
        rec["status"] = "new-root"
        rec["root_id"] = rid
        trials.append(rec)

    encoded_roots = []
    for r in sorted(roots, key=lambda q: (float(q.get("score", float("inf"))), int(q.get("id", 0)))):
        rr = dict(r)
        zc = rr.pop("z_complex")
        yc = rr.pop("y_complex")
        rr["z"] = ENGINE118.root_to_json(zc)
        rr["y"] = ENGINE118.root_to_json(yc)
        encoded_roots.append(rr)
    seconds = now() - t0
    return {
        "roots": encoded_roots,
        "trials": trials,
        "summary": {
            "requested_roots": int(args.count),
            "unique_roots": len(roots),
            "success": bool(len(roots) >= int(args.count)),
            "trials_used": len(trials),
            "duplicates": int(duplicates),
            "failures": int(failures),
            "refine_seconds": float(seconds),
            "eval_stats": system.stats(),
        },
    }


def run_123_case(args: argparse.Namespace, case_raw: str, helper_build: dict[str, Any], binary: Path) -> dict[str, Any]:
    t_case = now()
    n, d = ENGINE118.parse_case(case_raw)
    if d <= 33:
        raise ValueError(f"case {case_raw!r} is not degree > 33")
    system = ENGINE118.DenseKostlanSystem.make(
        n,
        d,
        seed_index=int(args.seed_index),
        equation_normalize=bool(args.equation_normalize),
    )
    input_path = Path(args.outdir) / "inputs" / f"123_{n}x{d}_starts{args.metal_candidates}.json"
    output_path = Path(args.outdir) / "inputs" / f"123_{n}x{d}_top{args.refine_top}.json"
    if bool(args.metal_generate_starts):
        t_gen = now()
        geom_meta = [{} for _ in range(int(args.metal_candidates))]
        if bool(getattr(args, "stateful_selector_server", False)):
            input_path.parent.mkdir(parents=True, exist_ok=True)
            stateful_payload = start_selector_params(args, system)
            input_path.write_text(json.dumps(stateful_payload, indent=2), encoding="utf-8")
        else:
            stateful_payload = None
            write_start_selector_input(input_path, args, system)
        generate_seconds = now() - t_gen
    else:
        stateful_payload = None
        points, trial_ids, geom_meta, generate_seconds = generate_batch_starts(args, system, case_raw)
        write_selector_input(input_path, system, points, trial_ids)
    selector_top = int(args.selector_top) if int(args.selector_top) > 0 else max(int(args.refine_top), int(args.refine_top) * 8)
    server_proc = getattr(args, "_selector_server_proc", None)
    if server_proc is not None and bool(getattr(args, "stateful_selector_server", False)):
        load_selector_system(server_proc, system)
        selector = run_stateful_selector_server(server_proc, stateful_payload, output_path)
    elif server_proc is not None:
        selector = run_selector_server(server_proc, input_path, output_path, selector_top)
    else:
        selector = run_selector(binary, input_path, output_path, selector_top)
    selected = diversify_selected(selector.get("selected", []), int(args.refine_top), float(args.diversity_sep))
    selector["selected_for_refine"] = selected
    selector["selected_for_refine_count"] = int(len(selected))
    selector["diversity_sep"] = float(args.diversity_sep)
    refined = refine_selected(args, system, selected, geom_meta)
    summary = refined["summary"]
    total_seconds = now() - t_case
    summary.update(
        {
            "generation_seconds": float(system.generation_seconds + generate_seconds),
            "metal_select_kernel_seconds": float(selector.get("kernel_seconds", 0.0)),
            "metal_select_process_seconds": float(selector.get("process_wall_seconds", 0.0)),
            "total_seconds": float(total_seconds),
        }
    )
    return {
        "script": "123_pandrosion_swift_metal_hybrid_vs_118_benchmark.py",
        "case": f"{n},{d}",
        "seed_index": int(args.seed_index),
        "seed": int(system.seed),
        "n": int(n),
        "degree": int(d),
        "terms_per_poly": int(system.terms_per_poly),
        "terms": int(system.total_terms),
        "bezout": int(system.bezout),
        "mode": "swift-metal-batch-selector/118-complex128-refine",
        "helper_build": helper_build,
        "selector": selector,
        "roots": refined["roots"],
        "trials": refined["trials"] if bool(args.verbose_trials) else refined["trials"][: min(len(refined["trials"]), int(args.keep_trials))],
        "summary": summary,
    }


def summarize_result(engine: str, case_raw: str, result: dict[str, Any], seconds: float) -> dict[str, Any]:
    summary = result.get("summary", {})
    eval_stats = summary.get("eval_stats", {})
    residuals = [float(r.get("residual", float("nan"))) for r in result.get("roots", [])]
    out = {
        "case": case_raw,
        "engine": engine,
        "status": "ok",
        "success": bool(summary.get("success")),
        "requested_roots": int(summary.get("requested_roots", 0)),
        "unique_roots": int(summary.get("unique_roots", 0)),
        "trials_used": int(summary.get("trials_used", 0)),
        "duplicates": int(summary.get("duplicates", 0)),
        "failures": int(summary.get("failures", 0)),
        "seconds": float(seconds),
        "reported_seconds": float(summary.get("total_seconds", seconds)),
        "generation_seconds": float(summary.get("generation_seconds", 0.0)),
        "refine_seconds": float(summary.get("refine_seconds", summary.get("extract_seconds", 0.0))),
        "metal_select_kernel_seconds": float(summary.get("metal_select_kernel_seconds", 0.0)),
        "metal_select_process_seconds": float(summary.get("metal_select_process_seconds", 0.0)),
        "eval_used": int(eval_stats.get("eval_count", 0)) if eval_stats.get("eval_count") is not None else None,
        "slope_calls": int(eval_stats.get("slope_count", 0)) if eval_stats.get("slope_count") is not None else None,
        "seconds_eval": float(eval_stats.get("seconds_eval", 0.0)),
        "seconds_slope": float(eval_stats.get("seconds_slope", 0.0)),
        "seed": int(result.get("seed", 0)),
        "terms": int(result.get("terms", 0)),
        "bezout": int(result.get("bezout", 0)),
        "error": None,
    }
    out.update(residual_stats(residuals))
    return out


def error_row(engine: str, case_raw: str, exc: BaseException, seconds: float) -> dict[str, Any]:
    return {
        "case": case_raw,
        "engine": engine,
        "status": "error",
        "success": False,
        "requested_roots": None,
        "unique_roots": 0,
        "trials_used": None,
        "duplicates": None,
        "failures": None,
        "seconds": float(seconds),
        "reported_seconds": None,
        "generation_seconds": None,
        "refine_seconds": None,
        "metal_select_kernel_seconds": None,
        "metal_select_process_seconds": None,
        "eval_used": None,
        "slope_calls": None,
        "seconds_eval": None,
        "seconds_slope": None,
        "seed": None,
        "terms": None,
        "bezout": None,
        "residual_min": None,
        "residual_mean": None,
        "residual_max": None,
        "error": f"{type(exc).__name__}: {exc}",
    }


def run_engine(engine: str, fn: Any, args: argparse.Namespace, case_raw: str, helper_build: Optional[dict[str, Any]] = None, binary: Optional[Path] = None) -> tuple[dict[str, Any], dict[str, Any]]:
    t0 = now()
    try:
        if helper_build is not None and binary is not None:
            result = fn(args, case_raw, helper_build, binary)
        else:
            result = fn(args, case_raw)
    except Exception as exc:
        return error_row(engine, case_raw, exc, now() - t0), {}
    return summarize_result(engine, case_raw, result, now() - t0), result


def compare_pair(row118: dict[str, Any], row123: dict[str, Any]) -> dict[str, Any]:
    roots118 = int(row118.get("unique_roots") or 0)
    roots123 = int(row123.get("unique_roots") or 0)
    seconds118 = float(row118.get("seconds") or 0.0)
    seconds123 = float(row123.get("seconds") or 0.0)
    if roots123 > roots118:
        winner = "123(hybrid)"
    elif roots118 > roots123:
        winner = "118"
    elif bool(row123.get("success")) and not bool(row118.get("success")):
        winner = "123(hybrid)"
    elif bool(row118.get("success")) and not bool(row123.get("success")):
        winner = "118"
    elif seconds123 > 0 and seconds118 > 0:
        winner = "123(hybrid)" if seconds123 < seconds118 else "118"
    else:
        winner = "tie"
    return {
        "case": row118.get("case"),
        "roots_118": roots118,
        "roots_123": roots123,
        "root_delta_123_minus_118": roots123 - roots118,
        "success_118": bool(row118.get("success")),
        "success_123": bool(row123.get("success")),
        "seconds_118": seconds118,
        "seconds_123": seconds123,
        "speedup_123_vs_118": (seconds118 / seconds123 if seconds123 > 0 else None),
        "evals_118": row118.get("eval_used"),
        "evals_123": row123.get("eval_used"),
        "slopes_118": row118.get("slope_calls"),
        "slopes_123": row123.get("slope_calls"),
        "residual_min_118": row118.get("residual_min"),
        "residual_min_123": row123.get("residual_min"),
        "winner": winner,
    }


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    fields = [
        "case",
        "engine",
        "status",
        "success",
        "requested_roots",
        "unique_roots",
        "trials_used",
        "duplicates",
        "failures",
        "seconds",
        "reported_seconds",
        "generation_seconds",
        "refine_seconds",
        "metal_select_kernel_seconds",
        "metal_select_process_seconds",
        "eval_used",
        "slope_calls",
        "seconds_eval",
        "seconds_slope",
        "seed",
        "terms",
        "bezout",
        "residual_min",
        "residual_mean",
        "residual_max",
        "error",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_pair_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    fields = [
        "case",
        "roots_118",
        "roots_123",
        "root_delta_123_minus_118",
        "success_118",
        "success_123",
        "seconds_118",
        "seconds_123",
        "speedup_123_vs_118",
        "evals_118",
        "evals_123",
        "slopes_118",
        "slopes_123",
        "residual_min_118",
        "residual_min_123",
        "winner",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    p = payload["parameters"]
    s = payload["summary"]
    lines = [
        "# 123 Swift Metal Hybrid vs 118 Benchmark",
        "",
        f"- cases: `{', '.join(payload['cases'])}`",
        f"- metal candidates/selector top/refine top: `{p['metal_candidates']}/{p['selector_top']}/{p['refine_top']}`",
        f"- count/pool/epochs: `{p['count']}/{p['pool']}/{p['epochs']}`",
        "- scope: Metal batch start selection, 118 complex128 refinement",
        "",
        "## Summary",
        "",
        f"- 118 roots: `{s['roots_118']}/{s['requested_118']}`; complete cases: `{s['successes_118']}/{s['pairs']}`",
        f"- 123 roots: `{s['roots_123']}/{s['requested_123']}`; complete cases: `{s['successes_123']}/{s['pairs']}`",
        f"- wall seconds: 118=`{s['seconds_118']:.4f}`, 123=`{s['seconds_123']:.4f}`, total=`{s['seconds']:.4f}`",
        f"- root delta 123-118: `{s['root_delta_123_minus_118']}`",
        "",
        "## Pair Results",
        "",
        "| case | 118 roots | 123 roots | delta | 118 sec | 123 sec | 123 speedup | 118 evals | 123 evals | winner |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in payload["pairs"]:
        lines.append(
            f"| `{row.get('case')}` | {row.get('roots_118')} | {row.get('roots_123')} | "
            f"{row.get('root_delta_123_minus_118')} | {fmt_float(row.get('seconds_118'), 4)} | "
            f"{fmt_float(row.get('seconds_123'), 4)} | {fmt_float(row.get('speedup_123_vs_118'), 3)} | "
            f"{row.get('evals_118') or ''} | {row.get('evals_123') or ''} | `{row.get('winner')}` |"
        )
    lines.extend(
        [
            "",
            "## Timing Breakdown",
            "",
            "| case | engine | generation | metal select process | refine/extract | eval sec | slope sec |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in payload["runs"]:
        lines.append(
            f"| `{row.get('case')}` | `{row.get('engine')}` | {fmt_float(row.get('generation_seconds'), 4)} | "
            f"{fmt_float(row.get('metal_select_process_seconds'), 4)} | {fmt_float(row.get('refine_seconds'), 4)} | "
            f"{fmt_float(row.get('seconds_eval'), 4)} | {fmt_float(row.get('seconds_slope'), 4)} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = ENGINE118.build_parser()
    p.description = "123 benchmark: Swift + Metal batch selector with 118 complex128 refine vs 118."
    p.set_defaults(
        cases="2,34",
        count=8,
        pool=512,
        epochs=24,
        outdir="verification/123_swift_metal_hybrid_vs_118_benchmark",
    )
    p.add_argument("--metal-candidates", type=int, default=8192)
    p.add_argument("--selector-top", type=int, default=0, help="0 means 8 * --refine-top")
    p.add_argument("--refine-top", type=int, default=192)
    p.add_argument("--diversity-sep", type=float, default=0.08)
    p.add_argument("--synthetic-pressure", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--metal-generate-starts", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--selector-server", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--stateful-selector-server", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--swift-source", default=str(FLOW / "123_swift_metal_select.swift"))
    p.add_argument("--swift-bin", default=None)
    p.add_argument("--rebuild-swift", action="store_true")
    p.add_argument("--run-order", choices=["118-first", "123-first"], default="118-first")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ENGINE118.ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    cases = parse_cases(args.cases)
    out_json = Path(args.out) if args.out else Path(args.outdir) / "123_swift_metal_hybrid_vs_118_benchmark.json"
    swift_source = Path(args.swift_source)
    if bool(args.stateful_selector_server):
        args.selector_server = True
        args.metal_generate_starts = True
    if bool(args.selector_server):
        if not bool(args.metal_generate_starts):
            raise ValueError("--selector-server currently requires --metal-generate-starts")
        if swift_source.name == "123_swift_metal_select.swift":
            swift_source = FLOW / "125_swift_metal_start_select_server.swift"
    binary_name = "125_swift_metal_start_select_server" if bool(args.selector_server) else "123_swift_metal_select"
    binary = Path(args.swift_bin) if args.swift_bin else Path(args.outdir) / "bin" / binary_name
    helper_build = build_swift_helper(swift_source, binary, force=bool(args.rebuild_swift))

    server_proc = None
    if bool(args.selector_server):
        server_proc = subprocess.Popen(
            [str(binary)],
            cwd=str(ROOT),
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        setattr(args, "_selector_server_proc", server_proc)
        warm_selector_server(server_proc)

    args118 = make_engine_args(args)
    rows: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    full_results: dict[str, Any] = {}
    t_all = now()

    print("=" * 120, flush=True)
    print("123 benchmark: Swift + Metal batch selector / 118 complex128 refine vs 118", flush=True)
    print(
        f"cases={cases} count={args.count} pool={args.pool} epochs={args.epochs} "
        f"metal_candidates={args.metal_candidates} selector_top={args.selector_top or int(args.refine_top) * 8} "
        f"refine_top={args.refine_top} metal_generate_starts={args.metal_generate_starts}",
        flush=True,
    )
    print("=" * 120, flush=True)

    try:
        for case_raw in cases:
            print(f"case={case_raw}", flush=True)
            if args.run_order == "123-first":
                row123, result123 = run_engine("123(hybrid)", run_123_case, args, case_raw, helper_build, binary)
                row118, result118 = run_engine("118", ENGINE118.run_case, args118, case_raw)
            else:
                row118, result118 = run_engine("118", ENGINE118.run_case, args118, case_raw)
                row123, result123 = run_engine("123(hybrid)", run_123_case, args, case_raw, helper_build, binary)

            rows.extend([row118, row123])
            print(
                f"  118: status={row118.get('status')} roots={row118.get('unique_roots')}/{row118.get('requested_roots')} "
                f"evals={row118.get('eval_used')} slopes={row118.get('slope_calls')} seconds={float(row118.get('seconds') or 0.0):.4f}",
                flush=True,
            )
            print(
                f"  123(hybrid): status={row123.get('status')} roots={row123.get('unique_roots')}/{row123.get('requested_roots')} "
                f"evals={row123.get('eval_used')} slopes={row123.get('slope_calls')} seconds={float(row123.get('seconds') or 0.0):.4f}",
                flush=True,
            )
            pair = compare_pair(row118, row123)
            pairs.append(pair)
            print(
                f"  winner={pair['winner']} delta123-118={pair['root_delta_123_minus_118']} "
                f"speedup123={fmt_float(pair.get('speedup_123_vs_118'), 3)}",
                flush=True,
            )
            full_results[case_raw.replace(",", "x")] = {"118": result118, "123(hybrid)": result123}
    finally:
        if server_proc is not None:
            try:
                if server_proc.stdin is not None:
                    server_proc.stdin.write("__quit__\n")
                    server_proc.stdin.flush()
            except Exception:
                pass
            try:
                server_proc.wait(timeout=5)
            except Exception:
                server_proc.kill()

    seconds_all = now() - t_all
    summary = {
        "pairs": len(pairs),
        "runs": len(rows),
        "successes_118": int(sum(1 for row in rows if row.get("engine") == "118" and row.get("success") is True)),
        "successes_123": int(sum(1 for row in rows if row.get("engine") == "123(hybrid)" and row.get("success") is True)),
        "roots_118": int(sum(int(row.get("unique_roots") or 0) for row in rows if row.get("engine") == "118")),
        "roots_123": int(sum(int(row.get("unique_roots") or 0) for row in rows if row.get("engine") == "123(hybrid)")),
        "requested_118": int(sum(int(row.get("requested_roots") or 0) for row in rows if row.get("engine") == "118")),
        "requested_123": int(sum(int(row.get("requested_roots") or 0) for row in rows if row.get("engine") == "123(hybrid)")),
        "seconds_118": float(sum(float(row.get("seconds") or 0.0) for row in rows if row.get("engine") == "118")),
        "seconds_123": float(sum(float(row.get("seconds") or 0.0) for row in rows if row.get("engine") == "123(hybrid)")),
        "root_delta_123_minus_118": int(sum(int(pair.get("root_delta_123_minus_118") or 0) for pair in pairs)),
        "wins_118": int(sum(1 for pair in pairs if pair.get("winner") == "118")),
        "wins_123": int(sum(1 for pair in pairs if pair.get("winner") == "123(hybrid)")),
        "ties": int(sum(1 for pair in pairs if pair.get("winner") == "tie")),
        "seconds": float(seconds_all),
    }

    payload = {
        "script": "123_pandrosion_swift_metal_hybrid_vs_118_benchmark.py",
        "engines": {
            "118": "flow/118_pandrosion_probe_aware_pure_thales_engine.py",
            "123(hybrid)": "Swift+Metal batch selector over 118 starts, then 118 complex128 refinement",
        },
        "cases": cases,
        "parameters": {
            "count": int(args.count),
            "pool": int(args.pool),
            "epochs": int(args.epochs),
            "accept": float(args.accept),
            "tol": float(args.tol),
            "cluster_sep": float(args.cluster_sep),
            "line_search": int(args.line_search),
            "probe_scale": float(args.probe_scale),
            "probe_candidates": int(args.probe_candidates),
            "probe_radii": str(args.probe_radii),
            "probe_self": bool(args.probe_self),
            "seed_index": int(args.seed_index),
            "metal_candidates": int(args.metal_candidates),
            "selector_top": int(args.selector_top) if int(args.selector_top) > 0 else int(args.refine_top) * 8,
            "refine_top": int(args.refine_top),
            "diversity_sep": float(args.diversity_sep),
            "synthetic_pressure": bool(args.synthetic_pressure),
            "metal_generate_starts": bool(args.metal_generate_starts),
            "selector_server": bool(args.selector_server),
            "stateful_selector_server": bool(args.stateful_selector_server),
            "run_order": str(args.run_order),
        },
        "summary": summary,
        "pairs": pairs,
        "runs": rows,
        "full_results": full_results if bool(args.verbose_trials) else {},
    }

    out_csv = out_json.with_suffix(".csv")
    out_pairs_csv = out_json.with_name(out_json.stem + "_pairs.csv")
    out_md = out_json.with_suffix(".md")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_csv(out_csv, rows)
    write_pair_csv(out_pairs_csv, pairs)
    write_markdown(out_md, payload)

    print("=" * 120, flush=True)
    print(
        f"summary: 118 roots={summary['roots_118']}/{summary['requested_118']} "
        f"successes={summary['successes_118']}/{summary['pairs']} seconds={summary['seconds_118']:.4f}",
        flush=True,
    )
    print(
        f"summary: 123 roots={summary['roots_123']}/{summary['requested_123']} "
        f"successes={summary['successes_123']}/{summary['pairs']} seconds={summary['seconds_123']:.4f}",
        flush=True,
    )
    print(f"root_delta_123_minus_118={summary['root_delta_123_minus_118']}", flush=True)
    print(f"json={out_json}", flush=True)
    print(f"csv={out_csv}", flush=True)
    print(f"pairs_csv={out_pairs_csv}", flush=True)
    print(f"report={out_md}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    finally:
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except Exception:
            pass
