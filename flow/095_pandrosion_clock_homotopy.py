"""
FLOW 095 -- Pandrosion-clock homotopy experiment.

This file does not modify 092.  It loads the complete 092 solver stack and
overrides only the path tracker, replacing the linear homotopy time t by a
Pandrosion-style clock lambda(t):

    H_t(z) = (1 - lambda(t)) G(z) + lambda(t) Gamma F(z).

The zero set path is the same as the usual gamma total-degree homotopy when
lambda is monotone from 0 to 1, but the geometry is traversed with a nonlinear
clock.  The default clock is the normalized geometric-sum clock

    lambda_p(t) = (t + t^2 + ... + t^p) / p,

which is deliberately conservative near the start and faster near the target.
All correction/polish steps remain Pandrosion-only through the 070/076/087
geometry portfolio.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path
from typing import Sequence

HERE = Path(__file__).resolve().parent
FLOW092_PATH = HERE / "092_pandrosion_resultant_pairing.py"


def _load_092():
    spec = importlib.util.spec_from_file_location("flow092_for_095", str(FLOW092_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW092_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow092_for_095"] = mod
    spec.loader.exec_module(mod)
    return mod


def _consume_095_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--pan-clock" and i + 1 < len(argv):
            os.environ["PANDROSION_095_CLOCK"] = argv[i + 1]
            i += 2
            continue
        if arg == "--pan-clock-power" and i + 1 < len(argv):
            os.environ["PANDROSION_095_POWER"] = argv[i + 1]
            i += 2
            continue
        if arg == "--pan-kappa" and i + 1 < len(argv):
            os.environ["PANDROSION_095_KAPPA"] = argv[i + 1]
            i += 2
            continue
        if arg == "--pan-dt-scale" and i + 1 < len(argv):
            os.environ["PANDROSION_095_DT_SCALE"] = argv[i + 1]
            i += 2
            continue
        out.append(arg)
        i += 1
    return out


def _apply_095_defaults(argv: list[str]) -> list[str]:
    out = list(argv)
    pair_flags = {"--pair-recovery", "--no-pair-recovery"}
    if not any(arg in pair_flags for arg in out[1:]):
        out.append("--pair-recovery")
    return out


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(os.environ.get(name, str(default))))
    except Exception:
        return default


def _clock(t: float) -> float:
    t = max(0.0, min(1.0, float(t)))
    kind = os.environ.get("PANDROSION_095_CLOCK", "geosum").strip().lower()
    if kind in {"linear", "identity"}:
        return t
    if kind in {"smooth", "smoothstep"}:
        return t * t * (3.0 - 2.0 * t)
    if kind in {"rational", "homographic"}:
        k = max(1e-6, _env_float("PANDROSION_095_KAPPA", 1.6))
        den = t + k * (1.0 - t)
        return 0.0 if den <= 0.0 else t / den
    if kind in {"front", "fast"}:
        k = max(1e-6, _env_float("PANDROSION_095_KAPPA", 0.55))
        den = t + k * (1.0 - t)
        return 0.0 if den <= 0.0 else t / den
    # Normalized geometric sum: (t + ... + t^p) / p.
    p = max(1, _env_int("PANDROSION_095_POWER", 3))
    s = 0.0
    tp = t
    for _ in range(p):
        s += tp
        tp *= t
    return max(0.0, min(1.0, s / p))


def _clock_delta(t0: float, t1: float) -> float:
    return max(0.0, _clock(t1) - _clock(t0))


def install_pandrosion_clock_tracker(mod) -> None:
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

    def track_one_clock(target, start, target_gamma, z0,
                        tol: float = 1e-9, max_steps: int = 420,
                        max_epochs: int = 12, quad_cap: int = 16,
                        modes: Sequence[str] = ("system", "integral_gl", "blend"),
                        rescue_modes: Sequence[str] = (),
                        deadline: float | None = None):
        z = list(z0)
        D = max(m.degrees(target))
        dt_scale = max(0.10, _env_float("PANDROSION_095_DT_SCALE", 1.0))
        dt = dt_scale * min(0.010, max(0.0015, 0.060 / max(2, D * D)))
        dtmax = dt_scale * min(0.030, max(0.006, 0.18 / max(2, D)))
        t = 0.0
        lam = _clock(0.0)
        steps = epochs = fails = guarded = 0
        prev_z = None
        prev_lam = None
        status = "ok"
        accept_residual = 80.0 * tol
        while t < 1.0 - 1e-15 and steps < max_steps:
            if m.timed_out(deadline):
                status = "budget"
                break
            steps += 1
            tnext = min(1.0, t + dt)
            lam_next = _clock(tnext)
            dlam = max(1e-15, lam_next - lam)
            pred = m.tangent_predictor(
                target_gamma,
                start,
                z,
                lam,
                dlam,
                prev_z=prev_z,
                prev_t=prev_lam,
            )
            Hnext = m.homotopy_polys(target_gamma, start, lam_next)
            zc, _ok, ep = m.corrector(
                Hnext,
                pred,
                tol=tol,
                max_epochs=max_epochs,
                quad_cap=quad_cap,
                modes=modes,
                rescue_modes=rescue_modes,
                deadline=deadline,
            )
            epochs += ep
            rh = m.residual_norm(Hnext, zc)
            ambiguous = correction_is_ambiguous(z, pred, zc)
            if math.isfinite(rh) and rh < accept_residual and not ambiguous:
                prev_z, prev_lam = list(z), lam
                z, t, lam = list(zc), tnext, lam_next
                dt = min(dtmax, dt * (1.20 if ep <= max(3, max_epochs // 3) else 1.06))
                fails = max(0, fails - 1)
                continue
            if math.isfinite(rh) and rh < 1e-6 and not ambiguous and not m.timed_out(deadline):
                zr, _okr, epr = m.corrector(
                    Hnext,
                    zc,
                    tol=tol,
                    max_epochs=max_epochs + 6,
                    quad_cap=quad_cap,
                    modes=("blend", "system", "integral_gl"),
                    rescue_modes=rescue_modes,
                    deadline=deadline,
                )
                epochs += epr
                rr = m.residual_norm(Hnext, zr)
                if math.isfinite(rr) and rr < accept_residual and not correction_is_ambiguous(z, pred, zr):
                    prev_z, prev_lam = list(z), lam
                    z, t, lam = list(zr), tnext, lam_next
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
            zp, _okp, ep = m.corrector(
                target,
                z,
                tol=tol,
                max_epochs=max_epochs + 10,
                quad_cap=quad_cap,
                modes=("system", "integral_gl", "blend"),
                rescue_modes=rescue_modes,
                deadline=deadline,
            )
            epochs += ep
            rp = m.residual_norm(target, zp)
            if math.isfinite(rp) and rp < 1e-7:
                z = zp
        res = m.residual_norm(target, z)
        ok_final = t >= 1.0 - 1e-12 and res < 1e-7
        return m.PathResult(z=z, ok=ok_final, steps=steps, epochs=epochs, t=t,
                            residual=res, status=status if not ok_final else "ok")

    m.track_one_070 = track_one_clock
    print(
        "095-pandrosion-clock tracker installed "
        f"clock={os.environ.get('PANDROSION_095_CLOCK', 'geosum')} "
        f"power={_env_int('PANDROSION_095_POWER', 3)} "
        f"kappa={_env_float('PANDROSION_095_KAPPA', 1.6):g}",
        flush=True,
    )


def main() -> None:
    sys.argv = _consume_095_args(sys.argv)
    sys.argv = _apply_095_defaults(sys.argv)
    f092 = _load_092()
    sys.argv = f092._consume_092_args(sys.argv)
    sys.argv = f092._apply_best_defaults(sys.argv)
    f091 = f092._load_091()
    sys.argv = f091._consume_091_args(sys.argv)
    f090 = f091._load_090()
    sys.argv = f090._consume_090_args(sys.argv)
    f089 = f090._load_089()
    sys.argv = f089._consume_089_args(sys.argv)
    w087 = f089._load_087()
    sys.argv = w087._consume_branch_args(sys.argv)
    mod = w087.load_081()
    w087.install_branch_safe_tracker(mod)
    install_pandrosion_clock_tracker(mod)
    w087.install_compat_and_single_policy(mod)
    # 081 launches batch workers through its module-level __file__.  Point it
    # back to this wrapper so subprocesses install the 095 tracker as well.
    mod.__file__ = str(Path(__file__).resolve())
    f089.install_resultant_recovery(mod)
    f090.install_multivariate_completion(mod)
    f091.install_projective_completion(mod)
    f092.install_pair_recovery(mod)
    sys.argv = f090._rewrite_argv_for_090(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
