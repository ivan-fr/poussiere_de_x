"""
FLOW 083 -- lifted/reprojected Pandrosion geometry.

083 keeps the algebraic system unchanged, but lets each local Pandrosion
geometry live in a larger generated frame.  For a current pair (a,z), it samples
more local secant directions than the n original coordinate directions, solves a
least-squares lifted slope in that overcomplete frame, and reprojects it back to
an exact n x n Pandrosion slope before the step is taken.

No extra unknown is added: paths still live in C^n and Bezout is unchanged.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path
from typing import List, Sequence

Complex = complex
Vector = List[Complex]
Matrix = List[List[Complex]]

HERE = Path(__file__).resolve().parent
FLOW081_PATH = HERE / "081_pandrosion_guarded_polygon_startopt.py"


def _consume_lift_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--lift-dim" and i + 1 < len(argv):
            os.environ["PANDROSION_LIFT_DIM"] = argv[i + 1]
            i += 2
            continue
        if arg == "--lift-radius" and i + 1 < len(argv):
            os.environ["PANDROSION_LIFT_RADIUS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--lift-blend" and i + 1 < len(argv):
            os.environ["PANDROSION_LIFT_BLEND"] = argv[i + 1]
            i += 2
            continue
        if arg == "--lift-ridge" and i + 1 < len(argv):
            os.environ["PANDROSION_LIFT_RIDGE"] = argv[i + 1]
            i += 2
            continue
        if arg == "--lift-off":
            os.environ["PANDROSION_LIFT_DIM"] = "0"
            i += 1
            continue
        out.append(arg)
        i += 1
    return out


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(os.environ.get(name, str(default))))
    except Exception:
        return default


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def load_081():
    spec = importlib.util.spec_from_file_location("flow081_for_083", str(FLOW081_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW081_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow081_for_083"] = mod
    spec.loader.exec_module(mod)
    return mod


def _lift_directions(anchor: Sequence[Complex], z: Sequence[Complex], lift_dim: int,
                     radius: float, safe_z) -> list[Vector]:
    n = len(z)
    delta = [z[j] - anchor[j] for j in range(n)]
    dirs: list[Vector] = []

    def add(d: Sequence[Complex]) -> None:
        if max((abs(v) for v in d), default=0.0) <= 1e-300:
            return
        point = [anchor[j] + d[j] for j in range(n)]
        if safe_z(point):
            dirs.append([complex(v) for v in d])

    # Original coordinate-edge frame.
    for j in range(n):
        d = [0.0 + 0.0j for _ in range(n)]
        d[j] = delta[j]
        add(d)

    if len(dirs) >= lift_dim:
        return dirs[:lift_dim]

    # Full chord and deterministic mixed chords.  These add cross-term geometry
    # but are always reprojected to the original C^n coordinates.
    add(delta)
    q = 0
    golden = 2.399963229728653
    radius = max(0.15, min(1.35, float(radius)))
    while len(dirs) < lift_dim and q < 4 * max(1, lift_dim + n):
        d: Vector = []
        for j in range(n):
            amp = radius * (0.45 + 0.18 * ((q + 2 * j) % 4))
            th = golden * (q + 1) * (j + 1)
            d.append(delta[j] * amp * complex(math.cos(th), math.sin(th)))
        add(d)
        q += 1
    return dirs[:lift_dim]


def install_lifted_geometry(mod) -> None:
    original_system = mod.m.Q_system_adaptive

    def lifted_system(polys, anchor: Sequence[Complex], z: Sequence[Complex],
                      F_anchor: Sequence[Complex]) -> Matrix:
        base = original_system(polys, anchor, z, F_anchor)
        n = len(z)
        lift_dim = _env_int("PANDROSION_LIFT_DIM", max(n + 2, 2 * n))
        if lift_dim <= n:
            return base
        radius = _env_float("PANDROSION_LIFT_RADIUS", 0.82)
        blend = max(0.0, min(1.0, _env_float("PANDROSION_LIFT_BLEND", 0.65)))
        ridge_factor = max(0.0, _env_float("PANDROSION_LIFT_RIDGE", 1e-11))

        dirs = _lift_directions(anchor, z, lift_dim, radius, mod.m.safe_z)
        if len(dirs) <= n:
            return base

        ys: list[Vector] = []
        for d in dirs:
            p = [anchor[j] + d[j] for j in range(n)]
            Fp = mod.m.F_eval(polys, p)
            if not mod.m.finite_vec(Fp):
                continue
            ys.append([Fp[i] - F_anchor[i] for i in range(n)])
        if len(ys) != len(dirs) or len(dirs) <= n:
            return base

        A = [[0.0 + 0.0j for _ in range(n)] for _ in range(n)]
        b_rows = [[0.0 + 0.0j for _ in range(n)] for _ in range(n)]
        for d, y in zip(dirs, ys):
            norm2 = max(1e-300, sum(abs(v) ** 2 for v in d))
            wt = 1.0 / norm2
            for j in range(n):
                cdj = d[j].conjugate() * wt
                for k in range(n):
                    A[j][k] += cdj * d[k]
                for row in range(n):
                    b_rows[row][j] += cdj * y[row]

        trace = sum(abs(A[j][j]) for j in range(n)) / max(1, n)
        ridge = ridge_factor * max(1.0, trace)
        for j in range(n):
            A[j][j] += ridge

        lifted: Matrix = []
        for row in range(n):
            q = mod.m.solve_linear(A, b_rows[row])
            if q is None or not mod.m.finite_vec(q):
                return base
            lifted.append(q)

        Fz = mod.m.F_eval(polys, z)
        if not mod.m.finite_vec(Fz):
            return base
        lifted = mod.m.exactify(
            lifted,
            F_anchor,
            Fz,
            anchor,
            z,
            mod.m.system_covector(polys, anchor, z, F_anchor),
        )
        if any(not mod.m.finite_vec(row) for row in lifted):
            return base
        return mod.m.blend_matrices(lifted, base, theta=blend)

    mod.m.Q_system_adaptive = lifted_system


def install_compat_and_policy(mod) -> None:
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
        lift_dim = _env_int("PANDROSION_LIFT_DIM", max(len(polys) + 2, 2 * len(polys)))
        lift_note = (
            f"lifted-geometry dim={lift_dim}, radius={_env_float('PANDROSION_LIFT_RADIUS', 0.82):.4g}, "
            f"blend={_env_float('PANDROSION_LIFT_BLEND', 0.65):.4g}; "
            f"no-fallback-single-policy={chosen.name}"
        )
        return [chosen], f"{note}; {lift_note}"

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
    if "--batch-worker" not in out and "--outdir" not in out:
        out.extend(["--outdir", "083_batches"])
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
    sys.argv = _consume_lift_args(sys.argv)
    mod = load_081()
    install_lifted_geometry(mod)
    install_compat_and_policy(mod)
    sys.argv = rewrite_argv(sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
