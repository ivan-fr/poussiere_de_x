"""
FLOW 086 -- projective-chart Pandrosion corrector.

086 keeps the fast 081 batch runner and the residual-gated tracker from 085,
then enriches the local corrector with projective charts.  When a predictor has
large coordinates, the same polynomial system is corrected in the reciprocal
chart y_j = 1 / z_j for those large coordinates, and the result is projected
back to the original C^n coordinates.

No unknowns are added and no Lairez/Newton fallback is used as a separate
algorithm; this is a local Pandrosion chart choice inside the corrector.
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path
from typing import Sequence

HERE = Path(__file__).resolve().parent
FLOW081_PATH = HERE / "081_pandrosion_guarded_polygon_startopt.py"


def load_081():
    spec = importlib.util.spec_from_file_location("flow081_for_086", str(FLOW081_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW081_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow081_for_086"] = mod
    spec.loader.exec_module(mod)
    return mod


def invert_chart_system(polys, mask: Sequence[bool]):
    out = []
    for poly in polys:
        powers = [
            max((alpha[j] for alpha in poly), default=0) if mask[j] else 0
            for j in range(len(mask))
        ]
        q = {}
        for alpha, coeff in poly.items():
            beta = tuple(powers[j] - alpha[j] if mask[j] else alpha[j] for j in range(len(mask)))
            q[beta] = q.get(beta, 0.0 + 0.0j) + coeff
        out.append({a: c for a, c in q.items() if abs(c) > 1e-14})
    return out


def z_to_chart(z: Sequence[complex], mask: Sequence[bool]):
    y = []
    for zj, inv in zip(z, mask):
        if inv:
            if abs(zj) < 1e-14:
                return None
            y.append(1.0 / zj)
        else:
            y.append(zj)
    return y


def chart_to_z(y: Sequence[complex], mask: Sequence[bool]):
    z = []
    for yj, inv in zip(y, mask):
        if inv:
            if abs(yj) < 1e-14:
                return None
            z.append(1.0 / yj)
        else:
            z.append(yj)
    return z


def install_projective_corrector(mod) -> None:
    m = mod.m
    base_corrector = m.corrector

    def projective_corrector(polys, z_init: Sequence[complex], tol: float = 1e-9,
                             max_epochs: int = 12, quad_cap: int = 16,
                             modes: Sequence[str] = ("system", "integral_gl", "blend"),
                             rescue_modes: Sequence[str] = (),
                             deadline: float | None = None):
        z_best, ok_best, ep_best = base_corrector(
            polys,
            z_init,
            tol=tol,
            max_epochs=max_epochs,
            quad_cap=quad_cap,
            modes=modes,
            rescue_modes=rescue_modes,
            deadline=deadline,
        )
        r_best = m.residual_norm(polys, z_best)
        if m.timed_out(deadline):
            return z_best, ok_best, ep_best

        scale = max((abs(v) for v in z_init), default=0.0)
        if scale < 3.25:
            return z_best, ok_best, ep_best

        masks = []
        large_mask = tuple(abs(v) >= 2.75 for v in z_init)
        if any(large_mask):
            masks.append(large_mask)
        all_mask = tuple(True for _ in z_init)
        if all_mask not in masks:
            masks.append(all_mask)

        total_ep = ep_best
        for mask in masks:
            if m.timed_out(deadline):
                break
            y0 = z_to_chart(z_init, mask)
            if y0 is None or not m.safe_z(y0):
                continue
            chart_polys = invert_chart_system(polys, mask)
            yc, okc, epc = base_corrector(
                chart_polys,
                y0,
                tol=tol,
                max_epochs=max(6, max_epochs),
                quad_cap=quad_cap,
                modes=modes,
                rescue_modes=rescue_modes,
                deadline=deadline,
            )
            total_ep += epc
            zc = chart_to_z(yc, mask)
            if zc is None or not m.safe_z(zc):
                continue
            rz = m.residual_norm(polys, zc)
            if math.isfinite(rz) and rz < r_best:
                z_best, r_best, ok_best = list(zc), rz, bool(okc and rz < 80.0 * tol)
                if rz < tol:
                    break
        return z_best, bool(ok_best or (math.isfinite(r_best) and r_best < tol)), total_ep

    m.corrector = projective_corrector


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
            zc, _ok, ep = m.corrector(Hnext, pred, tol=tol, max_epochs=max_epochs,
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
                zr, _okr, epr = m.corrector(Hnext, zc, tol=tol, max_epochs=max_epochs + 6,
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

    m.track_one_070 = track_one_residual_gate


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
        return [chosen], f"{note}; projective-chart-corrector; no-fallback-single-policy={chosen.name}"

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
        "--outdir": "086_batches",
        "--parallel-batches": "4",
        "--batch-size": "16",
        "--micro-batch": "2",
        "--retry-limit": "32",
        "--retries": "1",
        "--accept-t": "0.999999",
        "--max-epochs": "16",
        "--max-steps": "220",
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
    install_projective_corrector(mod)
    install_residual_gated_tracker(mod)
    install_compat_and_single_policy(mod)
    sys.argv = rewrite_argv(sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
