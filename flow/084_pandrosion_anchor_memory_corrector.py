"""
FLOW 084 -- anchor-memory Pandrosion corrector.

This wrapper keeps the fast batch runner from 081 but patches the local
Pandrosion corrector.  The important change is small: after a successful
Pandrosion step, the next anchor is the previous point, not the new point.

If the anchor is set equal to the new point, the next local geometry has
delta = z - a = 0 and the Pandrosion slope degenerates.  Keeping the previous
point as anchor preserves a non-zero local secant geometry without adding
Newton/Lairez fallback or extra unknowns.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import os
import sys
from pathlib import Path
from typing import Sequence

HERE = Path(__file__).resolve().parent
FLOW081_PATH = HERE / "081_pandrosion_guarded_polygon_startopt.py"


def load_081():
    spec = importlib.util.spec_from_file_location("flow081_for_084", str(FLOW081_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW081_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow081_for_084"] = mod
    spec.loader.exec_module(mod)
    return mod


def install_anchor_memory_corrector(mod) -> None:
    m = mod.m

    def corrector_anchor_memory(polys, z_init: Sequence[complex], tol: float = 1e-9,
                                max_epochs: int = 12, quad_cap: int = 16,
                                modes: Sequence[str] = ("system", "integral_gl", "blend"),
                                rescue_modes: Sequence[str] = (),
                                deadline: float | None = None):
        z = list(z_init)
        epochs_used = 0
        anchor = m.deterministic_anchor(z, 0)
        for epoch in range(max_epochs):
            if m.timed_out(deadline):
                return z, False, epochs_used
            epochs_used = epoch + 1
            rz = m.residual_norm(polys, z)
            if rz < tol:
                return z, True, epoch
            if not math.isfinite(rz):
                return z, False, epoch

            # Guard against accidental anchor collapse.
            if max((abs(z[i] - anchor[i]) for i in range(len(z))), default=0.0) <= 1e-14 * max(
                1.0, max((abs(v) for v in z), default=0.0)
            ):
                anchor = m.deterministic_anchor(z, epoch + 1, strength=0.060)

            Fa = m.F_eval(polys, anchor)
            if m.finite_vec(Fa) and max(abs(v) for v in Fa) < tol:
                return list(anchor), True, epoch

            best_r, best_z = rz, list(z)
            tried = set()
            for mode in modes:
                if m.timed_out(deadline):
                    break
                tried.add(mode)
                local_before = best_r
                cand, ok = m.pandrosion_h(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                if ok:
                    best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
                if best_r < 0.92 * local_before:
                    cand, ok = m.T2_step(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                    if ok:
                        best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
                if best_r < 0.22 * rz:
                    break

            if best_r > 0.60 * rz:
                for mode in rescue_modes:
                    if mode in tried or m.timed_out(deadline):
                        continue
                    local_before = best_r
                    cand, ok = m.pandrosion_h(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                    if ok:
                        best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
                    if best_r < 0.85 * local_before:
                        cand, ok = m.T2_step(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                        if ok:
                            best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
                    if best_r < 0.38 * rz:
                        break

            # No fallback: if this anchor did not expose progress, regenerate a
            # Pandrosion anchor and continue with the same current point.
            if best_r >= rz:
                anchor = m.deterministic_anchor(z, epoch + 1, strength=0.050)
                continue

            previous_z = list(z)
            z = list(best_z)
            anchor = previous_z
        return z, m.residual_norm(polys, z) < tol, epochs_used

    m.corrector = corrector_anchor_memory


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
        return [chosen], f"{note}; anchor-memory-corrector; no-fallback-single-policy={chosen.name}"

    mod.build_policies = build_policies_no_fallback


def install_residual_gated_tracker(mod) -> None:
    m = mod.m

    def track_one_residual_gate(target, start, target_gamma, z0,
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
        steps = epochs = fails = 0
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
            zc, ok, ep = m.corrector(Hnext, pred, tol=tol, max_epochs=max_epochs,
                                     quad_cap=quad_cap, modes=modes,
                                     rescue_modes=rescue_modes, deadline=deadline)
            epochs += ep
            rh = m.residual_norm(Hnext, zc)
            if math.isfinite(rh) and rh < accept_residual:
                prev_z, prev_t = list(z), t
                z, t = list(zc), tnext
                dt = min(dtmax, dt * (1.22 if ep <= max(3, max_epochs // 3) else 1.08))
                fails = max(0, fails - 1)
                continue
            if math.isfinite(rh) and rh < 1e-6 and not m.timed_out(deadline):
                zr, okr, epr = m.corrector(Hnext, zc, tol=tol, max_epochs=max_epochs + 6,
                                           quad_cap=quad_cap,
                                           modes=("blend", "system", "integral_gl"),
                                           rescue_modes=rescue_modes,
                                           deadline=deadline)
                epochs += epr
                rr = m.residual_norm(Hnext, zr)
                if math.isfinite(rr) and rr < accept_residual:
                    prev_z, prev_t = list(z), t
                    z, t = list(zr), tnext
                    dt = min(dtmax, dt * 1.04)
                    fails = max(0, fails - 1)
                    continue
            fails += 1
            dt *= 0.5
            if dt < 5e-7 or fails > 80:
                status = "fail"
                break
        if t >= 1.0 - 1e-12 and not m.timed_out(deadline):
            zp, okp, ep = m.corrector(target, z, tol=tol, max_epochs=max_epochs + 10,
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

    m.track_one_070 = track_one_residual_gate


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
        "--outdir": "084_batches",
        "--parallel-batches": "4",
        "--batch-size": "16",
        "--micro-batch": "2",
        "--retry-limit": "32",
        "--retries": "1",
        "--accept-t": "0.999999",
        "--max-epochs": "8",
        "--max-steps": "180",
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
    mod = load_081()
    install_anchor_memory_corrector(mod)
    install_residual_gated_tracker(mod)
    install_compat_and_single_policy(mod)
    sys.argv = rewrite_argv(sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
