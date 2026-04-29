#!/usr/bin/env python3
"""
117_pandrosion_solver_benchmark.py

Benchmark harness for comparing script 116 against classical/external solvers.

The script is intentionally conservative:
  - it reuses flow/116 for the exact same multifamily polynomial generators;
  - it runs Pandrosion-116 directly;
  - it runs a local SciPy multistart baseline when SciPy is installed;
  - it exports the exact same systems to Bertini, PHCpack, and Julia
    HomotopyContinuation.jl formats;
  - it runs those external solvers only when their commands are available;
  - it provides a configurable lairez_custom hook because Lairez's Smale-17
    algorithm is not a standardized command-line executable in this repo.

This makes the comparison reproducible and honest: unavailable solvers are
reported as skipped, not replaced by invented numbers.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]


def _load_engine_116() -> Any:
    path = Path(__file__).resolve().with_name("116_pandrosion_multifamily_vectorized_pure_pandrosion.py")
    spec = importlib.util.spec_from_file_location("pandrosion_116_multifamily", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load 116 engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ENGINE116 = _load_engine_116()
ENGINE115 = ENGINE116.ENGINE
np = ENGINE116.np


SOLVER_GROUPS: dict[str, list[str]] = {
    "all": ["pandrosion116", "scipy_multistart", "bertini", "phcpack", "julia_hc", "lairez_custom"],
    "local": ["pandrosion116", "scipy_multistart"],
    "external": ["bertini", "phcpack", "julia_hc", "lairez_custom"],
    "homotopy": ["bertini", "phcpack", "julia_hc", "lairez_custom"],
}

SOLVER_ALIASES = {
    "116": "pandrosion116",
    "pandrosion": "pandrosion116",
    "pandrosion_116": "pandrosion116",
    "scipy": "scipy_multistart",
    "scipy_root": "scipy_multistart",
    "bertini2": "bertini",
    "bertini_classic": "bertini",
    "phc": "phcpack",
    "phcpack": "phcpack",
    "hc": "julia_hc",
    "homotopycontinuation": "julia_hc",
    "homotopycontinuation_jl": "julia_hc",
    "julia": "julia_hc",
    "lairez": "lairez_custom",
    "lairez17": "lairez_custom",
}


def now() -> float:
    return time.perf_counter()


def ensure_numpy() -> None:
    ENGINE116.ensure_numpy()


def parse_solvers(raw: Optional[str]) -> list[str]:
    if raw is None or str(raw).strip() == "":
        return list(SOLVER_GROUPS["all"])
    out: list[str] = []
    for part in str(raw).replace("|", ",").replace(";", ",").split(","):
        token = part.strip().lower().replace("-", "_")
        if not token:
            continue
        names = SOLVER_GROUPS.get(token)
        if names is None:
            name = SOLVER_ALIASES.get(token, token)
            if name not in SOLVER_GROUPS["all"]:
                valid = ", ".join(SOLVER_GROUPS["all"])
                groups = ", ".join(sorted(SOLVER_GROUPS))
                raise ValueError(f"unknown solver {part!r}; valid solvers: {valid}; groups: {groups}")
            names = [name]
        for name in names:
            if name not in out:
                out.append(name)
    return out or list(SOLVER_GROUPS["all"])


def command_path(command: str) -> Optional[str]:
    parts = shlex.split(str(command))
    if not parts:
        return None
    return shutil.which(parts[0])


def has_scipy() -> bool:
    return importlib.util.find_spec("scipy") is not None


def system_id(case_raw: str, family: str) -> str:
    return f"{case_raw.replace(',', 'x')}_{family}"


def cjson(z: complex) -> list[float]:
    return [float(complex(z).real), float(complex(z).imag)]


def root_to_json(z: Sequence[complex]) -> list[list[float]]:
    return [cjson(v) for v in z]


def complex_literal(c: complex, style: str) -> str:
    z = complex(c)
    r = 0.0 if abs(z.real) < 5e-17 else z.real
    im = 0.0 if abs(z.imag) < 5e-17 else z.imag
    if style == "julia":
        return f"({r:.17g}{im:+.17g}im)"
    if style == "phc":
        return f"({r:.17g}{im:+.17g}*i)"
    return f"({r:.17g}{im:+.17g}*I)"


def monomial_expr(exps: Sequence[int], variables: Sequence[str], power_op: str = "^") -> str:
    factors = []
    for var, exp in zip(variables, exps):
        e = int(exp)
        if e <= 0:
            continue
        if e == 1:
            factors.append(var)
        else:
            factors.append(f"{var}{power_op}{e}")
    return "*".join(factors)


def polynomial_expr(system: Any, row: int, style: str) -> str:
    n = int(system.n)
    if style == "julia":
        variables = [f"x[{i + 1}]" for i in range(n)]
    else:
        variables = [f"x{i + 1}" for i in range(n)]
    pieces = []
    for coeff, exp in zip(system.coeff[row, :], system.exps):
        if abs(complex(coeff)) <= 1e-15:
            continue
        c = complex_literal(complex(coeff), style)
        mon = monomial_expr(exp, variables)
        pieces.append(c if not mon else f"{c}*{mon}")
    return " + ".join(pieces) if pieces else "0"


def write_system_json(system: Any, path: Path, case_raw: str, family: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "case": case_raw,
        "family": family,
        "n": int(system.n),
        "degree": int(system.d),
        "seed": int(system.seed),
        "terms_per_poly": int(system.terms_per_poly),
        "active_terms": int(getattr(system, "active_terms", system.total_terms)),
        "exponents": [[int(v) for v in row] for row in system.exps.tolist()],
        "coefficients": [[cjson(v) for v in row] for row in system.coeff.tolist()],
        "generator_metadata": dict(getattr(system, "generator_metadata", {})),
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_bertini_input(system: Any, path: Path) -> None:
    n = int(system.n)
    variables = ", ".join(f"x{i + 1}" for i in range(n))
    functions = ", ".join(f"f{i + 1}" for i in range(n))
    lines = [
        "CONFIG",
        "TrackType: 1;",
        "MPTYPE: 0;",
        "END;",
        "",
        "INPUT",
        f"variable_group {variables};",
        f"function {functions};",
    ]
    for i in range(n):
        lines.append(f"f{i + 1} = {polynomial_expr(system, i, 'bertini')};")
    lines.append("END;")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_phc_input(system: Any, path: Path) -> None:
    n = int(system.n)
    lines = [str(n)]
    for i in range(n):
        lines.append(polynomial_expr(system, i, "phc") + ";")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_julia_hc_script(system: Any, path: Path) -> None:
    n = int(system.n)
    polynomials = ",\n    ".join(polynomial_expr(system, i, "julia") for i in range(n))
    lines = [
        "using HomotopyContinuation",
        f"@var x[1:{n}]",
        f"F = System([{polynomials}])",
        "t0 = time()",
        "result = solve(F; show_progress=false)",
        "sols = solutions(result)",
        'println("HC_SOLUTIONS=", length(sols))',
        'println("HC_TIME=", time() - t0)',
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_solver_inputs(system: Any, base: Path, case_raw: str, family: str) -> dict[str, str]:
    sid = system_id(case_raw, family)
    folder = base / sid
    folder.mkdir(parents=True, exist_ok=True)
    paths = {
        "system_json": folder / f"{sid}.json",
        "bertini": folder / "input.bertini",
        "phcpack": folder / "input.phc",
        "julia_hc": folder / "solve_hc.jl",
    }
    write_system_json(system, paths["system_json"], case_raw, family)
    write_bertini_input(system, paths["bertini"])
    write_phc_input(system, paths["phcpack"])
    write_julia_hc_script(system, paths["julia_hc"])
    return {k: str(v) for k, v in paths.items()}


def make_system(args: argparse.Namespace, case_raw: str, family: str) -> Any:
    n, d = ENGINE115.parse_case(case_raw)
    return ENGINE116.MultiFamilySystem.make(
        n,
        d,
        seed_index=int(args.seed_index),
        equation_normalize=bool(args.equation_normalize),
        family=family,
        sparse_terms=int(args.sparse_terms),
        sparse_frac=float(args.sparse_frac),
    )


def residual_stats(values: Sequence[float]) -> dict[str, Optional[float]]:
    finite = [float(v) for v in values if math.isfinite(float(v))]
    if not finite:
        return {"residual_min": None, "residual_mean": None, "residual_max": None}
    return {
        "residual_min": float(min(finite)),
        "residual_mean": float(sum(finite) / len(finite)),
        "residual_max": float(max(finite)),
    }


def summarize_pandrosion_result(case_raw: str, family: str, result: dict[str, Any], seconds: float) -> dict[str, Any]:
    summary = result["summary"]
    eval_stats = summary.get("eval_stats", {})
    residuals = [float(r.get("residual", float("nan"))) for r in result.get("roots", [])]
    out = {
        "case": case_raw,
        "family": family,
        "solver": "pandrosion116",
        "status": "ok",
        "success": bool(summary.get("success")),
        "requested_roots": int(summary.get("requested_roots", 0)),
        "unique_roots": int(summary.get("unique_roots", 0)),
        "trials_used": int(summary.get("trials_used", 0)),
        "duplicates": int(summary.get("duplicates", 0)),
        "failures": int(summary.get("failures", 0)),
        "seconds": float(seconds),
        "solver_seconds_reported": float(summary.get("total_seconds", seconds)),
        "eval_used": int(eval_stats.get("eval_count", 0)) if eval_stats.get("eval_count") is not None else None,
        "slope_calls": int(eval_stats.get("slope_count", 0)) if eval_stats.get("slope_count") is not None else None,
        "eval_stats": eval_stats,
        "skip_reason": None,
    }
    out.update(residual_stats(residuals))
    return out


def run_pandrosion116(args: argparse.Namespace, case_raw: str, family: str) -> dict[str, Any]:
    t0 = now()
    try:
        result = ENGINE116.run_case(args, case_raw, family)
    except Exception as exc:
        return {
            "case": case_raw,
            "family": family,
            "solver": "pandrosion116",
            "status": "error",
            "success": False,
            "seconds": float(now() - t0),
            "error": f"{type(exc).__name__}: {exc}",
        }
    return summarize_pandrosion_result(case_raw, family, result, now() - t0)


def pack_complex(z: Any) -> Any:
    zz = np.asarray(z, dtype=np.complex128)
    return np.concatenate([zz.real, zz.imag])


def unpack_complex(x: Any) -> Any:
    xx = np.asarray(x, dtype=np.float64)
    n = xx.size // 2
    return xx[:n] + 1j * xx[n:]


def run_scipy_multistart(args: argparse.Namespace, case_raw: str, family: str, system: Any) -> dict[str, Any]:
    class ScipyEvalBudgetExhausted(RuntimeError):
        pass

    t0 = now()
    if not has_scipy():
        return skipped(case_raw, family, "scipy_multistart", "SciPy is not installed in this Python environment")
    try:
        import scipy.optimize as opt
    except Exception as exc:
        return skipped(case_raw, family, "scipy_multistart", f"SciPy import failed: {type(exc).__name__}: {exc}")

    n = int(system.n)
    chart = ENGINE115.LinearChart.identity(n, scale=float(args.linear_scale))
    roots: list[dict[str, Any]] = []
    residuals: list[float] = []
    duplicates = 0
    failures = 0
    starts = int(args.scipy_starts)
    maxfev = int(args.scipy_maxfev)
    raw_eval_budget = int(args.scipy_eval_budget)
    eval_budget = int(args.pool) if raw_eval_budget == 0 else (None if raw_eval_budget < 0 else max(0, raw_eval_budget))
    time_budget = float(args.scipy_time_budget)
    deadline = t0 + time_budget if time_budget > 0 else None
    eval_used = 0
    nfev_reported = 0
    trials_started = 0
    stop_reason = "starts_exhausted"

    powers = sorted(set(round(float(x), 16) for x in ENGINE115.parse_float_list(args.powers, ENGINE115.DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in ENGINE115.parse_float_list(args.angles, ENGINE115.DEFAULT_ANGLES_DEG)]
    radii = ENGINE115.parse_float_list(args.rays, ENGINE115.DEFAULT_RADII, positive=True)

    def eval_with_budget(z: Any) -> Any:
        nonlocal eval_used
        if eval_budget is not None and eval_used >= eval_budget:
            raise ScipyEvalBudgetExhausted("SciPy global F-evaluation budget exhausted")
        f = system.eval(z)
        eval_used += 1
        return f

    def fun(x: Any) -> Any:
        z = unpack_complex(x)
        f = eval_with_budget(z)
        return np.concatenate([f.real, f.imag])

    for trial in range(max(0, starts)):
        if len(roots) >= int(args.count):
            stop_reason = "target_count_reached"
            break
        if deadline is not None and now() >= deadline:
            stop_reason = "time_budget_exhausted"
            break
        if eval_budget is not None and eval_used >= eval_budget:
            stop_reason = "eval_budget_exhausted"
            break
        try:
            y0, _ = ENGINE115.mobius_homothety_start(
                n,
                trial,
                int(system.seed) + 0x117000,
                powers,
                angles,
                radii,
                float(args.power_cap),
                roots_found=len(roots),
                duplicates=duplicates,
                failures=failures,
                target_count=int(args.count),
            )
            x0 = pack_complex(chart.z_from_y(y0))
            trials_started += 1
            if eval_budget is None:
                maxfev_this = maxfev
            else:
                remaining = max(0, eval_budget - eval_used)
                if remaining <= 1:
                    stop_reason = "eval_budget_exhausted"
                    break
                # Reserve one evaluation for the final validation residual.
                maxfev_this = max(1, min(maxfev, remaining - 1))
            sol = opt.root(fun, x0, method=str(args.scipy_method), options={"maxfev": maxfev_this})
            try:
                nfev_reported += int(getattr(sol, "nfev", 0))
            except Exception:
                pass
            z = unpack_complex(sol.x)
            r = float(np.linalg.norm(eval_with_budget(z)))
            residuals.append(r)
        except ScipyEvalBudgetExhausted:
            stop_reason = "eval_budget_exhausted"
            break
        except Exception:
            failures += 1
            continue
        if not (math.isfinite(r) and r < float(args.accept)):
            failures += 1
            continue
        dup = ENGINE115.cluster_index(roots, z, float(args.cluster_sep))
        if dup is not None:
            duplicates += 1
            continue
        roots.append({"z_complex": np.asarray(z, dtype=np.complex128).copy(), "residual": r, "trial": trial})
    else:
        if len(roots) >= int(args.count):
            stop_reason = "target_count_reached"

    accepted_residuals = [float(r["residual"]) for r in roots]
    out = {
        "case": case_raw,
        "family": family,
        "solver": "scipy_multistart",
        "status": "ok",
        "success": bool(len(roots) >= int(args.count)),
        "requested_roots": int(args.count),
        "unique_roots": int(len(roots)),
        "trials_used": int(trials_started),
        "duplicates": int(duplicates),
        "failures": int(failures),
        "seconds": float(now() - t0),
        "method": str(args.scipy_method),
        "maxfev": int(maxfev),
        "eval_budget": None if eval_budget is None else int(eval_budget),
        "eval_used": int(eval_used),
        "nfev_reported": int(nfev_reported),
        "starts_budget": int(starts),
        "time_budget": None if time_budget <= 0 else float(time_budget),
        "stop_reason": stop_reason,
        "skip_reason": None,
    }
    out.update(residual_stats(accepted_residuals))
    if bool(args.keep_solver_roots):
        out["roots"] = [
            {"trial": int(r["trial"]), "residual": float(r["residual"]), "z": root_to_json(r["z_complex"])}
            for r in roots
        ]
    return out


def skipped(case_raw: str, family: str, solver: str, reason: str, extra: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    out = {
        "case": case_raw,
        "family": family,
        "solver": solver,
        "status": "skipped",
        "success": False,
        "requested_roots": None,
        "unique_roots": None,
        "trials_used": None,
        "duplicates": None,
        "failures": None,
        "seconds": 0.0,
        "skip_reason": reason,
    }
    if extra:
        out.update(extra)
    return out


def error_result(case_raw: str, family: str, solver: str, seconds: float, exc: Exception) -> dict[str, Any]:
    return {
        "case": case_raw,
        "family": family,
        "solver": solver,
        "status": "error",
        "success": False,
        "seconds": float(seconds),
        "error": f"{type(exc).__name__}: {exc}",
        "skip_reason": None,
    }


def parse_external_output(text: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key in ("HC_SOLUTIONS", "HC_TIME"):
        m = re.search(rf"{key}\s*=\s*([0-9.eE+-]+)", text)
        if m:
            val = m.group(1)
            out[key.lower()] = float(val) if "." in val or "e" in val.lower() else int(val)
    generic = re.search(r"(?i)(?:number of )?(?:finite )?solutions?\D{0,30}(\d+)", text)
    if generic and "external_solutions" not in out:
        out["external_solutions"] = int(generic.group(1))
    return out


def run_subprocess(command: Sequence[str], cwd: Path, timeout: float) -> tuple[int, str, str, float]:
    t0 = now()
    proc = subprocess.run(
        list(command),
        cwd=str(cwd),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=float(timeout) if timeout and timeout > 0 else None,
    )
    return int(proc.returncode), proc.stdout, proc.stderr, float(now() - t0)


def run_bertini(args: argparse.Namespace, case_raw: str, family: str, paths: dict[str, str]) -> dict[str, Any]:
    solver = "bertini"
    cmd_path = command_path(args.bertini_command)
    if cmd_path is None:
        return skipped(case_raw, family, solver, f"command not found: {args.bertini_command}", {"input_file": paths.get("bertini")})
    if not bool(args.run_external):
        return skipped(case_raw, family, solver, "external execution disabled; input file generated", {"input_file": paths.get("bertini")})
    input_path = Path(paths["bertini"])
    command = shlex.split(args.bertini_command) + [input_path.name]
    try:
        code, stdout, stderr, seconds = run_subprocess(command, input_path.parent, float(args.external_timeout))
    except Exception as exc:
        return error_result(case_raw, family, solver, 0.0, exc)
    parsed = parse_external_output(stdout + "\n" + stderr)
    return {
        "case": case_raw,
        "family": family,
        "solver": solver,
        "status": "ok" if code == 0 else "nonzero-exit",
        "success": code == 0,
        "seconds": seconds,
        "returncode": code,
        "input_file": str(input_path),
        "stdout_tail": stdout[-4000:],
        "stderr_tail": stderr[-4000:],
        "skip_reason": None,
        **parsed,
    }


def run_phcpack(args: argparse.Namespace, case_raw: str, family: str, paths: dict[str, str]) -> dict[str, Any]:
    solver = "phcpack"
    cmd_path = command_path(args.phc_command)
    if cmd_path is None:
        return skipped(case_raw, family, solver, f"command not found: {args.phc_command}", {"input_file": paths.get("phcpack")})
    if not bool(args.run_external):
        return skipped(case_raw, family, solver, "external execution disabled; input file generated", {"input_file": paths.get("phcpack")})
    input_path = Path(paths["phcpack"])
    output_path = input_path.parent / "phc_output.txt"
    command = shlex.split(args.phc_command) + ["-b", str(input_path), str(output_path)]
    try:
        code, stdout, stderr, seconds = run_subprocess(command, input_path.parent, float(args.external_timeout))
    except Exception as exc:
        return error_result(case_raw, family, solver, 0.0, exc)
    output_text = output_path.read_text(encoding="utf-8", errors="replace") if output_path.exists() else ""
    parsed = parse_external_output(stdout + "\n" + stderr + "\n" + output_text)
    return {
        "case": case_raw,
        "family": family,
        "solver": solver,
        "status": "ok" if code == 0 else "nonzero-exit",
        "success": code == 0,
        "seconds": seconds,
        "returncode": code,
        "input_file": str(input_path),
        "output_file": str(output_path),
        "stdout_tail": stdout[-4000:],
        "stderr_tail": stderr[-4000:],
        "skip_reason": None,
        **parsed,
    }


def run_julia_hc(args: argparse.Namespace, case_raw: str, family: str, paths: dict[str, str]) -> dict[str, Any]:
    solver = "julia_hc"
    cmd_path = command_path(args.julia_command)
    if cmd_path is None:
        return skipped(case_raw, family, solver, f"command not found: {args.julia_command}", {"input_file": paths.get("julia_hc")})
    if not bool(args.run_external):
        return skipped(case_raw, family, solver, "external execution disabled; Julia script generated", {"input_file": paths.get("julia_hc")})
    script_path = Path(paths["julia_hc"])
    command = shlex.split(args.julia_command) + [str(script_path)]
    try:
        code, stdout, stderr, seconds = run_subprocess(command, script_path.parent, float(args.external_timeout))
    except Exception as exc:
        return error_result(case_raw, family, solver, 0.0, exc)
    parsed = parse_external_output(stdout + "\n" + stderr)
    return {
        "case": case_raw,
        "family": family,
        "solver": solver,
        "status": "ok" if code == 0 else "nonzero-exit",
        "success": code == 0,
        "seconds": seconds,
        "returncode": code,
        "input_file": str(script_path),
        "stdout_tail": stdout[-4000:],
        "stderr_tail": stderr[-4000:],
        "skip_reason": None,
        **parsed,
    }


def run_lairez_custom(args: argparse.Namespace, case_raw: str, family: str, paths: dict[str, str]) -> dict[str, Any]:
    solver = "lairez_custom"
    if not str(args.lairez_command or "").strip():
        return skipped(
            case_raw,
            family,
            solver,
            "no standard Lairez CLI is configured; pass --lairez-command with placeholders such as {system_json}",
            {"system_json": paths.get("system_json")},
        )
    command_text = str(args.lairez_command).format(
        system_json=paths.get("system_json"),
        bertini_input=paths.get("bertini"),
        phc_input=paths.get("phcpack"),
        julia_script=paths.get("julia_hc"),
        case=case_raw,
        family=family,
    )
    command = shlex.split(command_text)
    if not command:
        return skipped(case_raw, family, solver, "empty --lairez-command after formatting", {"system_json": paths.get("system_json")})
    if command_path(command[0]) is None and not Path(command[0]).exists():
        return skipped(case_raw, family, solver, f"command not found: {command[0]}", {"system_json": paths.get("system_json")})
    try:
        code, stdout, stderr, seconds = run_subprocess(command, Path(paths["system_json"]).parent, float(args.external_timeout))
    except Exception as exc:
        return error_result(case_raw, family, solver, 0.0, exc)
    parsed = parse_external_output(stdout + "\n" + stderr)
    return {
        "case": case_raw,
        "family": family,
        "solver": solver,
        "status": "ok" if code == 0 else "nonzero-exit",
        "success": code == 0,
        "seconds": seconds,
        "returncode": code,
        "system_json": paths.get("system_json"),
        "stdout_tail": stdout[-4000:],
        "stderr_tail": stderr[-4000:],
        "skip_reason": None,
        **parsed,
    }


def run_solver(args: argparse.Namespace, solver: str, case_raw: str, family: str, system: Any, paths: dict[str, str]) -> dict[str, Any]:
    if solver == "pandrosion116":
        return run_pandrosion116(args, case_raw, family)
    if solver == "scipy_multistart":
        return run_scipy_multistart(args, case_raw, family, system)
    if solver == "bertini":
        return run_bertini(args, case_raw, family, paths)
    if solver == "phcpack":
        return run_phcpack(args, case_raw, family, paths)
    if solver == "julia_hc":
        return run_julia_hc(args, case_raw, family, paths)
    if solver == "lairez_custom":
        return run_lairez_custom(args, case_raw, family, paths)
    return skipped(case_raw, family, solver, "unsupported solver")


def availability(args: argparse.Namespace, solvers: Sequence[str]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for solver in solvers:
        if solver == "pandrosion116":
            out[solver] = {"available": True, "detail": "local Python engine"}
        elif solver == "scipy_multistart":
            out[solver] = {"available": has_scipy(), "detail": "SciPy local baseline"}
        elif solver == "bertini":
            out[solver] = {"available": command_path(args.bertini_command) is not None, "command": args.bertini_command}
        elif solver == "phcpack":
            out[solver] = {"available": command_path(args.phc_command) is not None, "command": args.phc_command}
        elif solver == "julia_hc":
            out[solver] = {"available": command_path(args.julia_command) is not None, "command": args.julia_command}
        elif solver == "lairez_custom":
            out[solver] = {"available": bool(str(args.lairez_command or "").strip()), "command": args.lairez_command}
    return out


def write_csv(path: Path, runs: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case",
        "family",
        "solver",
        "status",
        "success",
        "requested_roots",
        "unique_roots",
        "trials_used",
        "duplicates",
        "failures",
        "eval_budget",
        "eval_used",
        "slope_calls",
        "stop_reason",
        "seconds",
        "residual_min",
        "residual_mean",
        "residual_max",
        "skip_reason",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in runs:
            writer.writerow(row)


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    rows = payload["runs"]
    lines = [
        "# 117 solver benchmark",
        "",
        f"- cases: `{', '.join(payload['cases'])}`",
        f"- families: `{', '.join(payload['families'])}`",
        f"- solvers: `{', '.join(payload['solvers'])}`",
        f"- run_external: `{payload['parameters']['run_external']}`",
        "",
        "## Availability",
        "",
        "| solver | available | detail |",
        "|---|---:|---|",
    ]
    for name, item in payload["availability"].items():
        detail = item.get("detail") or item.get("command") or ""
        lines.append(f"| `{name}` | `{item.get('available')}` | {detail} |")
    lines.extend(["", "## Results", "", "| case | family | solver | status | roots | evals | seconds | residual min | note |", "|---|---|---|---|---:|---:|---:|---:|---|"])
    for r in rows:
        roots = "" if r.get("unique_roots") is None else f"{r.get('unique_roots')}/{r.get('requested_roots')}"
        if r.get("eval_used") is None:
            evals = ""
        elif r.get("eval_budget") is None:
            evals = str(r.get("eval_used"))
        else:
            evals = f"{r.get('eval_used')}/{r.get('eval_budget')}"
        sec = "" if r.get("seconds") is None else f"{float(r.get('seconds', 0.0)):.4f}"
        rmin = "" if r.get("residual_min") is None else f"{float(r.get('residual_min')):.3e}"
        note = r.get("skip_reason") or r.get("error") or r.get("stop_reason") or ""
        lines.append(f"| `{r.get('case')}` | `{r.get('family')}` | `{r.get('solver')}` | `{r.get('status')}` | {roots} | {evals} | {sec} | {rmin} | {note} |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = ENGINE116.build_parser()
    p.description = "117 benchmark: compare Pandrosion-116 against optional external homotopy solvers and local baselines."
    p.set_defaults(cases="3,4", families="all", count=4, pool=512, outdir="verification/117_solver_benchmark")
    p.add_argument("--solvers", default="all", help="solver list or group: all/local/external/homotopy")
    p.add_argument("--list-solvers", action="store_true", help="print available solver adapters and exit")
    p.add_argument("--scipy-starts", type=int, default=128, help="number of local SciPy starts per case/family")
    p.add_argument("--scipy-method", default="hybr", help="scipy.optimize.root method")
    p.add_argument("--scipy-maxfev", type=int, default=2000)
    p.add_argument("--scipy-eval-budget", type=int, default=0, help="global SciPy F-evaluation budget; 0 means use --pool, -1 means unlimited")
    p.add_argument("--scipy-time-budget", type=float, default=0.0, help="optional wall-time budget in seconds for SciPy per case/family; 0 disables")
    p.add_argument("--run-external", action=argparse.BooleanOptionalAction, default=True, help="run external solvers when commands exist")
    p.add_argument("--external-timeout", type=float, default=120.0)
    p.add_argument("--bertini-command", default="bertini")
    p.add_argument("--phc-command", default="phc")
    p.add_argument("--julia-command", default="julia")
    p.add_argument("--lairez-command", default="", help="custom command template; placeholders: {system_json}, {case}, {family}")
    p.add_argument("--inputs-dir", default=None)
    p.add_argument("--keep-solver-roots", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(getattr(args, "list_families", False)):
        print("117 benchmark families inherited from 116")
        for name in ENGINE116.DEFAULT_FAMILIES:
            print(f"{name}: {ENGINE116.FAMILY_DESCRIPTIONS[name]}")
        print("groups: " + ", ".join(f"{k}={','.join(v)}" for k, v in ENGINE116.FAMILY_GROUPS.items()))
        return 0
    if bool(getattr(args, "list_solvers", False)):
        print("117 solver adapters")
        print("pandrosion116: local flow/116 engine")
        print("scipy_multistart: local scipy.optimize.root multistart baseline")
        print("bertini: Bertini input export and optional command execution")
        print("phcpack: PHCpack input export and optional phc execution")
        print("julia_hc: Julia HomotopyContinuation.jl script export and optional execution")
        print("lairez_custom: configurable command hook; no standard Lairez CLI is assumed")
        print("groups: " + ", ".join(f"{k}={','.join(v)}" for k, v in SOLVER_GROUPS.items()))
        return 0
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    families = ENGINE116.parse_families(args.families)
    solvers = parse_solvers(args.solvers)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs_dir = Path(args.inputs_dir) if args.inputs_dir else outdir / "solver_inputs"

    runs: list[dict[str, Any]] = []
    t_all = now()
    print("=" * 120, flush=True)
    print("117 Pandrosion-116 solver benchmark", flush=True)
    print("Same generated systems are exported to optional external solvers; unavailable solvers are skipped.", flush=True)
    print("=" * 120, flush=True)
    print(f"cases={cases} families={families} solvers={solvers}", flush=True)

    for case_raw in cases:
        for family in families:
            system = make_system(args, case_raw, family)
            paths = export_solver_inputs(system, inputs_dir, case_raw, family)
            print(f"case={case_raw} family={family} n={system.n} d={system.d} active_terms={getattr(system, 'active_terms', system.total_terms)}", flush=True)
            for solver in solvers:
                result = run_solver(args, solver, case_raw, family, system, paths)
                runs.append(result)
                roots = result.get("unique_roots")
                requested = result.get("requested_roots")
                roots_s = "skip" if roots is None else f"{roots}/{requested}"
                evals = ""
                if result.get("eval_used") is not None:
                    if result.get("eval_budget") is None:
                        evals = f" evals={result.get('eval_used')}"
                    else:
                        evals = f" evals={result.get('eval_used')}/{result.get('eval_budget')}"
                note = result.get("skip_reason") or result.get("error") or ""
                print(f"  {solver}: status={result.get('status')} roots={roots_s}{evals} seconds={float(result.get('seconds', 0.0)):.3f} {note}", flush=True)

    payload = {
        "script": "117_pandrosion_solver_benchmark.py",
        "cases": cases,
        "families": families,
        "solvers": solvers,
        "availability": availability(args, solvers),
        "parameters": {
            "count": int(args.count),
            "pool": int(args.pool),
            "accept": float(args.accept),
            "tol": float(args.tol),
            "cluster_sep": float(args.cluster_sep),
            "epochs": int(args.epochs),
            "scipy_starts": int(args.scipy_starts),
            "scipy_eval_budget": int(args.scipy_eval_budget),
            "scipy_time_budget": float(args.scipy_time_budget),
            "run_external": bool(args.run_external),
            "external_timeout": float(args.external_timeout),
            "inputs_dir": str(inputs_dir),
        },
        "summary": {
            "runs": len(runs),
            "ok": int(sum(1 for r in runs if r.get("status") == "ok")),
            "skipped": int(sum(1 for r in runs if r.get("status") == "skipped")),
            "errors": int(sum(1 for r in runs if r.get("status") == "error")),
            "successes": int(sum(1 for r in runs if r.get("success") is True)),
            "seconds": float(now() - t_all),
        },
        "runs": runs,
    }

    out_json = Path(args.out) if args.out else outdir / "117_solver_benchmark.json"
    out_csv = out_json.with_suffix(".csv")
    out_md = out_json.with_suffix(".md")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_csv(out_csv, runs)
    write_markdown(out_md, payload)
    print(f"json={out_json}", flush=True)
    print(f"csv={out_csv}", flush=True)
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
