"""
FLOW 082 -- guarded polygon/startopt runner targeting KS(2,10).

082 reuses the fast batch-guarded orchestration from 081, but fixes the local
compatibility gap with the current 076: 081 expects f076.polish_070_compat,
which is not present in this checkout.  This wrapper installs the missing
Pandrosion-only polish helper in both parent and worker processes, then makes
081 relaunch this file for its workers.

No Lairez fallback and no Newton-ELS are used.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FLOW081_PATH = HERE / "081_pandrosion_guarded_polygon_startopt.py"


def load_081():
    spec = importlib.util.spec_from_file_location("flow081_for_082", str(FLOW081_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW081_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow081_for_082"] = mod
    spec.loader.exec_module(mod)
    return mod


def install_compat(mod) -> None:
    def polish_070_compat(target, z, tol, quad_cap):
        return mod.f076.m.polish_070(
            target,
            z,
            tol,
            quad_cap,
            mod.f076.DEFAULT_MODES,
            mod.f076.DEFAULT_RESCUE_MODES,
            None,
        )

    mod.f076.polish_070_compat = polish_070_compat
    # 081.launch_batch uses Path(__file__) for worker subprocesses.  Point it at
    # this wrapper so every worker gets the same compatibility patch.
    mod.__file__ = str(Path(__file__).resolve())

    original_build_policies = mod.build_policies

    def build_policies_no_fallback(polys, args):
        policies, note = original_build_policies(polys, args)
        chosen = None
        if getattr(args, "force_policy", ""):
            chosen = next((p for p in policies if p.name == args.force_policy), None)
        if chosen is None:
            chosen = next((p for p in policies if p.name == "system-unit"), policies[0])
        return [chosen], f"{note}; no-fallback-single-policy={chosen.name}"

    mod.build_policies = build_policies_no_fallback


def rewrite_argv(argv: list[str]) -> list[str]:
    out = list(argv)
    if "--case" in out and "--batch-worker" not in out and "--cases" not in out:
        i = out.index("--case")
        if i + 1 < len(out):
            case = out[i + 1]
            del out[i:i + 2]
            out.extend(["--cases", case])
    if "--batch-worker" not in out and "--cases" not in out:
        out.extend(["--cases", "2,10"])
    if "--batch-worker" not in out and "--equation-normalize" not in out:
        out.append("--equation-normalize")
    if "--batch-worker" not in out and "--batch-timeout" not in out:
        out.extend(["--batch-timeout", "20"])
    if "--batch-worker" not in out and "--parallel-batches" not in out:
        out.extend(["--parallel-batches", "4"])
    if "--batch-worker" not in out and "--batch-size" not in out:
        out.extend(["--batch-size", "16"])
    if "--batch-worker" not in out and "--micro-batch" not in out:
        out.extend(["--micro-batch", "2"])
    if "--batch-worker" not in out and "--retry-limit" not in out:
        out.extend(["--retry-limit", "32"])
    if "--batch-worker" not in out and "--retries" not in out:
        out.extend(["--retries", "3"])
    if "--batch-worker" not in out and "--stop-at-bezout" not in out:
        out.append("--stop-at-bezout")
    return out


def main() -> None:
    mod = load_081()
    install_compat(mod)
    sys.argv = rewrite_argv(sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
