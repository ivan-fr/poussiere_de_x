"""
FLOW 087 -- branch-safe residual-gated Pandrosion tracker.

087 extends 085 with a continuation guard: a corrected point must have a low
residual and remain compatible with the tangent predictor.  If the correction is
too large compared with both the local scale and the predictor move, the step is
treated as ambiguous and the homotopy step is reduced instead of accepting a
possible branch jump.

This is still Pandrosion-only: no Lairez/Newton fallback and no extra unknowns.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path
from typing import Sequence

HERE = Path(__file__).resolve().parent
FLOW081_PATH = HERE / "081_pandrosion_guarded_polygon_startopt.py"


def _consume_branch_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--branch-scale" and i + 1 < len(argv):
            os.environ["PANDROSION_BRANCH_SCALE"] = argv[i + 1]
            i += 2
            continue
        if arg == "--branch-move" and i + 1 < len(argv):
            os.environ["PANDROSION_BRANCH_MOVE"] = argv[i + 1]
            i += 2
            continue
        out.append(arg)
        i += 1
    return out


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def load_081():
    spec = importlib.util.spec_from_file_location("flow081_for_087", str(FLOW081_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW081_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow081_for_087"] = mod
    spec.loader.exec_module(mod)
    return mod


def install_branch_safe_tracker(mod) -> None:
    m = mod.m

    def correction_is_ambiguous(z, pred, zc) -> bool:
        corr = max((abs(zc[i] - pred[i]) for i in range(len(z))), default=0.0)
        move = max((abs(pred[i] - z[i]) for i in range(len(z))), default=0.0)
        scale = max(
            1.0,
            max((abs(v) for v in z), default=0.0),
            max((abs(v) for v in pred), default=0.0),
            max((abs(v) for v in zc), default=0.0),
        )
        scale_gate = max(0.02, _env_float("PANDROSION_BRANCH_SCALE", 0.10))
        move_gate = max(0.50, _env_float("PANDROSION_BRANCH_MOVE", 1.25))
        return corr > max(scale_gate * scale, move_gate * max(move, 1e-12))

    def track_one_branch_safe(target, start, target_gamma, z0,
                              tol: float = 1e-9, max_steps: int = 420,
                              max_epochs: int = 12, quad_cap: int = 16,
                              modes: Sequence[str] = ("system", "integral_gl", "blend"),
                              rescue_modes: Sequence[str] = (),
                              deadline: float | None = None):
        z = list(z0)
        D = max(m.degrees(target))
        dt = min(0.010, max(0.0015, 0.060 / max(2, D * D)))
        dtmax = min(0.030, max(0.006, 0.18 / max(2, D)))
        t = 0.0
        steps = epochs = fails = guarded = 0
        prev_z = None
        prev_t = None
        status = "ok"
        accept_residual = 80.0 * tol
        while t < 1.0 - 1e-15 and steps < max_steps:
            if m.timed_out(deadline):
                status = "budget"
                break
            steps += 1
            tnext = min(1.0, t + dt)
            pred = m.tangent_predictor(target_gamma, start, z, t, tnext - t, prev_z=prev_z, prev_t=prev_t)
            Hnext = m.homotopy_polys(target_gamma, start, tnext)
            zc, _ok, ep = m.corrector(Hnext, pred, tol=tol, max_epochs=max_epochs,
                                      quad_cap=quad_cap, modes=modes,
                                      rescue_modes=rescue_modes, deadline=deadline)
            epochs += ep
            rh = m.residual_norm(Hnext, zc)
            ambiguous = correction_is_ambiguous(z, pred, zc)
            if math.isfinite(rh) and rh < accept_residual and not ambiguous:
                prev_z, prev_t = list(z), t
                z, t = list(zc), tnext
                dt = min(dtmax, dt * (1.20 if ep <= max(3, max_epochs // 3) else 1.06))
                fails = max(0, fails - 1)
                continue
            if math.isfinite(rh) and rh < 1e-6 and not ambiguous and not m.timed_out(deadline):
                zr, _okr, epr = m.corrector(Hnext, zc, tol=tol, max_epochs=max_epochs + 6,
                                            quad_cap=quad_cap,
                                            modes=("blend", "system", "integral_gl"),
                                            rescue_modes=rescue_modes,
                                            deadline=deadline)
                epochs += epr
                rr = m.residual_norm(Hnext, zr)
                if math.isfinite(rr) and rr < accept_residual and not correction_is_ambiguous(z, pred, zr):
                    prev_z, prev_t = list(z), t
                    z, t = list(zr), tnext
                    dt = min(dtmax, dt * 1.03)
                    fails = max(0, fails - 1)
                    continue
            if ambiguous:
                guarded += 1
            fails += 1
            dt *= 0.5
            if dt < 5e-7 or fails > 80:
                status = "guarded" if guarded else "fail"
                break
        if t >= 1.0 - 1e-12 and not m.timed_out(deadline):
            zp, _okp, ep = m.corrector(target, z, tol=tol, max_epochs=max_epochs + 10,
                                       quad_cap=quad_cap,
                                       modes=("system", "integral_gl", "blend"),
                                       rescue_modes=rescue_modes,
                                       deadline=deadline)
            epochs += ep
            rp = m.residual_norm(target, zp)
            if math.isfinite(rp) and rp < 1e-7:
                z = zp
        res = m.residual_norm(target, z)
        ok_final = t >= 1.0 - 1e-12 and res < 1e-7
        return m.PathResult(z=z, ok=ok_final, steps=steps, epochs=epochs, t=t,
                            residual=res, status=status if not ok_final else "ok")

    m.track_one_070 = track_one_branch_safe


def install_compat_and_single_policy(mod) -> None:
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
    mod.__file__ = str(Path(__file__).resolve())

    original_build_policies = mod.build_policies

    def build_policies_no_fallback(polys, args):
        policies, note = original_build_policies(polys, args)
        chosen = None
        if getattr(args, "force_policy", ""):
            chosen = next((p for p in policies if p.name == args.force_policy), None)
        if chosen is None:
            chosen = next((p for p in policies if p.name == "system-unit"), policies[0])
        scale_gate = _env_float("PANDROSION_BRANCH_SCALE", 0.10)
        move_gate = _env_float("PANDROSION_BRANCH_MOVE", 1.25)
        return [chosen], (
            f"{note}; branch-safe scale_gate={scale_gate:.4g}, move_gate={move_gate:.4g}; "
            f"no-fallback-single-policy={chosen.name}"
        )

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
    defaults = {
        "--batch-timeout": "20",
        "--outdir": "087_batches",
        "--parallel-batches": "4",
        "--batch-size": "16",
        "--micro-batch": "2",
        "--retry-limit": "32",
        "--retries": "1",
        "--accept-t": "0.999999",
        "--max-epochs": "16",
        "--max-steps": "240",
    }
    if "--batch-worker" not in out:
        if "--equation-normalize" not in out:
            out.append("--equation-normalize")
        if "--stop-at-bezout" not in out:
            out.append("--stop-at-bezout")
        for flag, value in defaults.items():
            if flag not in out:
                out.extend([flag, value])
    return out


def main() -> None:
    sys.argv = _consume_branch_args(sys.argv)
    mod = load_081()
    install_branch_safe_tracker(mod)
    install_compat_and_single_policy(mod)
    sys.argv = rewrite_argv(sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
