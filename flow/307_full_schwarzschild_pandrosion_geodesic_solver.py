#!/usr/bin/env python3
"""
307_full_schwarzschild_pandrosion_geodesic_solver.py

Standalone NumPy solver for equatorial Schwarzschild geodesics using the
Pandrosion-relativity protocol of Papers 027--028 and the autonomous style of
flow/304.

This is not a SciPy integrator wrapper and it does not import previous flow
scripts.  It contains its own:
  - NumPy bootstrap and JSON export;
  - Schwarzschild and Painleve-Gullstrand chart formulas;
  - Hamiltonian constraint projection;
  - finite causal/proper-time probe scoring;
  - dynamic atlas conditioning diagnostics;
  - RK4 reference integration for verification.

Core 307 flow:

    state y_k = (T, r, phi, p_r)
      -> finite radial/phase probe atlas near a Hamiltonian predictor
      -> project each probe to the Schwarzschild Hamiltonian shell
      -> score by proper time, finite-slope geodesic defect, constraint,
         tidal stiffness, and dynamic chart condition
      -> select y_{k+1}
      -> record horizon crossing, conservation errors, and perihelion advance

The physical equations are classical.  The Pandrosion contribution is the
finite-probe organization: causal endpoint candidates, proper-time scoring,
Thales-like local scale control, and dynamic atlas selection.

Dependencies: Python stdlib + NumPy only.
"""
from __future__ import annotations

import argparse
import dataclasses
import glob
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


def _bootstrap_numpy_path() -> None:
    """Expose usual site-packages directories when launched with python -S."""
    here = Path(__file__).resolve()
    candidates: list[str] = []
    for pat in (
        str(here.parents[1] / ".venv/lib/python*/site-packages"),
        str(here.parents[1] / "venv/lib/python*/site-packages"),
        "/mnt/data/venv/lib/python*/site-packages",
        "/usr/local/lib/python*/site-packages",
        "/usr/lib/python*/dist-packages",
        "/usr/lib/python*/site-packages",
    ):
        candidates.extend(glob.glob(pat))
    for path in candidates:
        if path not in sys.path:
            sys.path.append(path)


try:
    import numpy as np
except Exception as exc:  # pragma: no cover
    _bootstrap_numpy_path()
    try:
        import numpy as np
    except Exception as exc2:  # pragma: no cover
        np = None
        _NUMPY_IMPORT_ERROR = exc2
    else:
        _NUMPY_IMPORT_ERROR = None
else:
    _NUMPY_IMPORT_ERROR = None


def ensure_numpy() -> None:
    if np is None:
        raise RuntimeError(f"NumPy is required by 307. Import error: {_NUMPY_IMPORT_ERROR!r}")


def now() -> float:
    return time.time()


def finite_float(x: Any, default: float = float("nan")) -> float:
    try:
        y = float(x)
        return y if math.isfinite(y) else default
    except Exception:
        return default


def array_to_rows(a: Any, max_rows: Optional[int] = None) -> list[list[float]]:
    ensure_numpy()
    arr = np.asarray(a, dtype=float)
    if max_rows is not None:
        arr = arr[: max(0, int(max_rows))]
    return [[finite_float(v) for v in row] for row in arr.reshape((-1, arr.shape[-1]))]


def json_safe(obj: Any) -> Any:
    ensure_numpy()
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    return obj


@dataclasses.dataclass
class SchwarzschildConfig:
    mass: float = 1.0
    kappa: float = 1.0
    energy: float = 1.0
    angular_momentum: float = 0.0
    min_radius: float = 1.0e-3

    def f(self, r: float) -> float:
        return 1.0 - 2.0 * self.mass / max(float(r), self.min_radius)

    def a_pg(self, r: float) -> float:
        return math.sqrt(max(0.0, 2.0 * self.mass / max(float(r), self.min_radius)))

    def fprime(self, r: float) -> float:
        rr = max(float(r), self.min_radius)
        return 2.0 * self.mass / (rr * rr)

    def aprime_pg(self, r: float) -> float:
        rr = max(float(r), self.min_radius)
        a = self.a_pg(rr)
        return -0.5 * a / rr


@dataclasses.dataclass
class SolverWeights:
    defect: float = 1000.0
    constraint: float = 25.0
    atlas: float = 0.002
    tidal: float = 0.0
    pr_smooth: float = 0.0005


@dataclasses.dataclass
class SolverParams:
    scenario: str
    step: float
    steps: int
    stop_radius: float
    probe_candidates: int
    probe_radius: float
    record_stride: int
    dynamic_atlas: bool
    reference: bool
    weights: SolverWeights


def pg_hamiltonian(y: Sequence[float], cfg: SchwarzschildConfig) -> float:
    """Hamiltonian H = 1/2 g^{mu nu} p_mu p_nu in PG coordinates."""
    _, r, _, pr = [float(v) for v in y]
    f = cfg.f(r)
    a = cfg.a_pg(r)
    e = cfg.energy
    ell = cfg.angular_momentum
    return 0.5 * (-e * e - 2.0 * a * e * pr + f * pr * pr + (ell * ell) / (r * r))


def constraint_error(y: Sequence[float], cfg: SchwarzschildConfig) -> float:
    return float(pg_hamiltonian(y, cfg) + 0.5 * cfg.kappa)


def radial_potential(r: float, cfg: SchwarzschildConfig) -> float:
    f = cfg.f(r)
    ell = cfg.angular_momentum
    return float(cfg.energy * cfg.energy - f * (cfg.kappa + ell * ell / (r * r)))


def project_pr_to_constraint(r: float, pr_guess: float, cfg: SchwarzschildConfig) -> Optional[float]:
    """Project p_r to the Hamiltonian shell H=-kappa/2 in PG coordinates."""
    r = float(r)
    if r <= cfg.min_radius:
        return None
    f = cfg.f(r)
    a = cfg.a_pg(r)
    e = cfg.energy
    ell = cfg.angular_momentum
    # f pr^2 - 2 a E pr + (L^2/r^2 + kappa - E^2) = 0.
    A = f
    B = -2.0 * a * e
    C = ell * ell / (r * r) + cfg.kappa - e * e
    if abs(A) < 1.0e-10:
        if abs(B) < 1.0e-14:
            return None
        val = -C / B
        return float(val) if math.isfinite(val) else None
    disc = B * B - 4.0 * A * C
    if disc < -1.0e-12:
        return None
    disc = max(0.0, disc)
    root = math.sqrt(disc)
    p1 = (-B + root) / (2.0 * A)
    p2 = (-B - root) / (2.0 * A)
    if not math.isfinite(p1) and not math.isfinite(p2):
        return None
    if not math.isfinite(p1):
        return float(p2)
    if not math.isfinite(p2):
        return float(p1)
    return float(p1 if abs(p1 - pr_guess) <= abs(p2 - pr_guess) else p2)


def initial_pr_from_direction(r0: float, sign: float, cfg: SchwarzschildConfig) -> float:
    """Build p_r from desired radial direction sign(d r/d lambda)."""
    r0 = float(r0)
    R = radial_potential(r0, cfg)
    if R < -1.0e-12:
        raise ValueError(f"inadmissible initial radius: radial potential {R:.6g} < 0")
    dr = math.copysign(math.sqrt(max(0.0, R)), float(sign) if sign != 0 else -1.0)
    f = cfg.f(r0)
    a = cfg.a_pg(r0)
    guess = (dr + a * cfg.energy) / f if abs(f) > 1.0e-10 else 0.0
    pr = project_pr_to_constraint(r0, guess, cfg)
    if pr is None:
        raise ValueError("could not project initial p_r to Hamiltonian shell")
    return float(pr)


def pg_rhs(y: Sequence[float], cfg: SchwarzschildConfig) -> "np.ndarray":
    ensure_numpy()
    T, r, phi, pr = [float(v) for v in y]
    r = max(r, cfg.min_radius)
    f = cfg.f(r)
    a = cfg.a_pg(r)
    fp = cfg.fprime(r)
    ap = cfg.aprime_pg(r)
    e = cfg.energy
    ell = cfg.angular_momentum
    dT = e + a * pr
    dr = -a * e + f * pr
    dphi = ell / (r * r)
    dpr = ap * e * pr - 0.5 * fp * pr * pr + ell * ell / (r * r * r)
    return np.asarray([dT, dr, dphi, dpr], dtype=float)


def rk4_step(y: Sequence[float], h: float, cfg: SchwarzschildConfig, project: bool = True) -> "np.ndarray":
    ensure_numpy()
    yy = np.asarray(y, dtype=float)
    k1 = pg_rhs(yy, cfg)
    k2 = pg_rhs(yy + 0.5 * h * k1, cfg)
    k3 = pg_rhs(yy + 0.5 * h * k2, cfg)
    k4 = pg_rhs(yy + h * k3, cfg)
    out = yy + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    out[1] = max(float(out[1]), cfg.min_radius)
    if project:
        pr = project_pr_to_constraint(float(out[1]), float(out[3]), cfg)
        if pr is not None:
            out[3] = pr
    return out


def pg_segment_ds2(y0: Sequence[float], y1: Sequence[float], cfg: SchwarzschildConfig) -> float:
    T0, r0, phi0, _ = [float(v) for v in y0]
    T1, r1, phi1, _ = [float(v) for v in y1]
    rm = max(0.5 * (r0 + r1), cfg.min_radius)
    f = cfg.f(rm)
    a = cfg.a_pg(rm)
    dT = T1 - T0
    dr = r1 - r0
    dphi = phi1 - phi0
    return float(-f * dT * dT + 2.0 * a * dT * dr + dr * dr + rm * rm * dphi * dphi)


def segment_score_tau(y0: Sequence[float], y1: Sequence[float], cfg: SchwarzschildConfig) -> tuple[float, float]:
    ds2 = pg_segment_ds2(y0, y1, cfg)
    if cfg.kappa > 0.5:
        tau2 = -ds2
        tau = math.sqrt(max(0.0, tau2))
        return float(tau), float(tau2)
    # Null mode: best score is closest to ds2=0.
    return float(-math.sqrt(abs(ds2))), float(-abs(ds2))


def chart_costs(r: float, cfg: SchwarzschildConfig) -> dict[str, float]:
    f = cfg.f(r)
    eps = 1.0e-4
    sch = 1.0 + 0.02 / max(abs(f), eps) ** 2
    if r <= 2.0 * cfg.mass:
        sch += 1.0e6
    pg = 1.28 + 0.08 * (2.0 * cfg.mass / max(r, cfg.min_radius)) ** 2
    return {"schwarzschild": float(sch), "painleve_gullstrand": float(pg)}


def choose_chart(r: float, cfg: SchwarzschildConfig, dynamic: bool = True) -> tuple[str, float, dict[str, float]]:
    costs = chart_costs(r, cfg)
    if not dynamic:
        return "painleve_gullstrand", costs["painleve_gullstrand"], costs
    if costs["schwarzschild"] <= costs["painleve_gullstrand"]:
        return "schwarzschild", costs["schwarzschild"], costs
    return "painleve_gullstrand", costs["painleve_gullstrand"], costs


def finite_slope_defect(y0: "np.ndarray", y1: "np.ndarray", h: float, cfg: SchwarzschildConfig) -> float:
    ensure_numpy()
    mid = 0.5 * (np.asarray(y0, dtype=float) + np.asarray(y1, dtype=float))
    vel = (np.asarray(y1, dtype=float) - np.asarray(y0, dtype=float)) / max(float(h), 1.0e-300)
    rhs = pg_rhs(mid, cfg)
    scale = np.maximum(1.0, np.abs(rhs))
    return float(np.linalg.norm((vel - rhs) / scale))


def probe_offsets(count: int) -> list[float]:
    base = [0.0, -0.5, 0.5, -1.0, 1.0, -1.6, 1.6, -2.4, 2.4, -3.2, 3.2]
    if count <= len(base):
        return base[: max(1, count)]
    out = list(base)
    k = 4.0
    while len(out) < count:
        out.extend([-k, k])
        k *= 1.35
    return out[:count]


def pandrosion_probe_step(
    y: Sequence[float],
    h: float,
    cfg: SchwarzschildConfig,
    params: SolverParams,
) -> tuple["np.ndarray", dict[str, Any]]:
    """One finite-probe Pandrosion geodesic step in PG coordinates."""
    ensure_numpy()
    yy = np.asarray(y, dtype=float)
    rhs0 = pg_rhs(yy, cfg)
    h = float(h)
    predictor = rk4_step(yy, h, cfg, project=True)
    euler_probe = yy + h * rhs0
    rhs_probe = pg_rhs(euler_probe, cfg)
    acc_probe = (rhs_probe - rhs0) / max(h, 1.0e-300)
    euler = predictor.copy()
    radial_scale = max(abs(h * rhs0[1]), abs(0.5 * h * h * acc_probe[1]), 0.03 * abs(h), 1.0e-8)
    pr_scale = max(abs(h * rhs0[3]), abs(0.5 * h * h * acc_probe[3]), 0.03 * abs(h), 1.0e-8)
    best: Optional[tuple[float, np.ndarray, dict[str, Any]]] = None
    tested = 0
    admissible = 0

    for q in probe_offsets(int(params.probe_candidates)):
        tested += 1
        if q == 0.0:
            r1 = float(predictor[1])
            pr_guess = float(predictor[3])
        else:
            r1 = float(euler[1] + float(q) * float(params.probe_radius) * radial_scale)
            pr_guess = float(euler[3] + 0.25 * float(q) * float(params.probe_radius) * pr_scale)
        if r1 <= cfg.min_radius or (params.stop_radius > 0 and r1 < 0.5 * params.stop_radius):
            continue
        pr1 = project_pr_to_constraint(r1, pr_guess, cfg)
        if pr1 is None:
            continue
        trial = euler.copy()
        trial[1] = r1
        trial[3] = pr1
        rhs1 = pg_rhs(trial, cfg)
        trial[0] = yy[0] + 0.5 * h * (rhs0[0] + rhs1[0])
        trial[2] = yy[2] + 0.5 * h * (rhs0[2] + rhs1[2])
        if not np.all(np.isfinite(trial)):
            continue
        tau, causal_margin = segment_score_tau(yy, trial, cfg)
        if cfg.kappa > 0.5 and causal_margin < -1.0e-10:
            continue
        rm = max(0.5 * (float(yy[1]) + float(trial[1])), cfg.min_radius)
        chart, chart_cost, costs = choose_chart(rm, cfg, params.dynamic_atlas)
        defect = finite_slope_defect(yy, trial, h, cfg)
        con = abs(constraint_error(trial, cfg))
        tidal = cfg.mass / max(rm**3, cfg.min_radius**3)
        smooth = abs(float(trial[3]) - float(yy[3]))
        score = (
            tau
            - params.weights.defect * h * defect * defect
            - params.weights.constraint * con * con
            - params.weights.atlas * h * chart_cost
            - params.weights.tidal * h * tidal
            - params.weights.pr_smooth * smooth * smooth
        )
        admissible += 1
        meta = {
            "score": float(score),
            "tau": float(tau),
            "causal_margin": float(causal_margin),
            "finite_slope_defect": float(defect),
            "constraint_error": float(con),
            "tidal": float(tidal),
            "chart": chart,
            "chart_cost": float(chart_cost),
            "chart_costs": costs,
            "probe_offset": float(q),
        }
        if best is None or score > best[0]:
            best = (float(score), trial.copy(), meta)

    if best is None:
        # Last-resort deterministic finite step: use RK4 and project.  It is
        # flagged in metadata and should not normally activate.
        fallback = rk4_step(yy, h, cfg, project=True)
        chart, chart_cost, costs = choose_chart(float(fallback[1]), cfg, params.dynamic_atlas)
        tau, causal_margin = segment_score_tau(yy, fallback, cfg)
        meta = {
            "score": float(tau),
            "tau": float(tau),
            "causal_margin": float(causal_margin),
            "finite_slope_defect": finite_slope_defect(yy, fallback, h, cfg),
            "constraint_error": abs(constraint_error(fallback, cfg)),
            "tidal": cfg.mass / max(float(fallback[1]) ** 3, cfg.min_radius**3),
            "chart": chart,
            "chart_cost": float(chart_cost),
            "chart_costs": costs,
            "probe_offset": None,
            "fallback": "rk4-projected-after-no-admissible-probe",
        }
        return fallback, {"tested": int(tested), "admissible": int(admissible), **meta}

    return best[1], {"tested": int(tested), "admissible": int(admissible), **best[2]}


def bound_orbit_energy_angular_momentum(p: float, e: float, mass: float) -> tuple[float, float, float]:
    """Schwarzschild bound-orbit constants for semi-latus p and eccentricity e.

    Formula uses p in units of M.  Returns (E,L,r_periastron).
    """
    pp = float(p)
    ee = float(e)
    den = pp * (pp - 3.0 - ee * ee)
    if den <= 0:
        raise ValueError("bound orbit needs p > 3 + e^2")
    E2 = ((pp - 2.0 - 2.0 * ee) * (pp - 2.0 + 2.0 * ee)) / den
    L2 = (mass * mass * pp * pp) / (pp - 3.0 - ee * ee)
    if E2 <= 0 or L2 <= 0:
        raise ValueError("invalid p,e produced nonpositive E^2 or L^2")
    rp = mass * pp / (1.0 + ee)
    return math.sqrt(E2), math.sqrt(L2), rp


def scenario_defaults(args: argparse.Namespace, scenario: str) -> tuple[SchwarzschildConfig, np.ndarray, SolverParams, dict[str, Any]]:
    ensure_numpy()
    mass = float(args.mass)
    kappa = float(args.kappa)
    meta: dict[str, Any] = {"scenario": scenario}

    if scenario == "radial-plunge":
        energy = float(args.energy) if args.energy is not None else 1.0
        ell = float(args.angular_momentum) if args.angular_momentum is not None else 0.0
        r0 = float(args.r0) if args.r0 is not None else 8.0 * mass
        sign = -1.0
        step = float(args.step) if args.step > 0 else 0.01 * mass
        steps = int(args.steps) if args.steps > 0 else 1100
        stop_radius = float(args.stop_radius) if args.stop_radius > 0 else 1.65 * mass
    elif scenario == "precession":
        p = float(args.p)
        ecc = float(args.eccentricity)
        energy, ell, r0 = bound_orbit_energy_angular_momentum(p, ecc, mass)
        if args.energy is not None:
            energy = float(args.energy)
        if args.angular_momentum is not None:
            ell = float(args.angular_momentum)
        if args.r0 is not None:
            r0 = float(args.r0)
        sign = 1.0
        step = float(args.step) if args.step > 0 else 0.02 * mass
        steps = int(args.steps) if args.steps > 0 else 19000
        stop_radius = float(args.stop_radius) if args.stop_radius > 0 else 2.05 * mass
        meta.update({"p": p, "eccentricity": ecc})
    elif scenario == "custom":
        energy = float(args.energy) if args.energy is not None else 1.0
        ell = float(args.angular_momentum) if args.angular_momentum is not None else 3.8 * mass
        r0 = float(args.r0) if args.r0 is not None else 10.0 * mass
        sign = -1.0 if str(args.radial_direction).lower().startswith("in") else 1.0
        step = float(args.step) if args.step > 0 else 0.01 * mass
        steps = int(args.steps) if args.steps > 0 else 2000
        stop_radius = float(args.stop_radius) if args.stop_radius > 0 else 2.05 * mass
    else:
        raise ValueError(f"unknown scenario {scenario!r}")

    cfg = SchwarzschildConfig(
        mass=mass,
        kappa=kappa,
        energy=energy,
        angular_momentum=ell,
        min_radius=max(1.0e-6, float(args.min_radius)),
    )
    pr0 = initial_pr_from_direction(r0, sign, cfg)
    y0 = np.asarray([0.0, r0, float(args.phi0), pr0], dtype=float)
    params = SolverParams(
        scenario=scenario,
        step=step,
        steps=steps,
        stop_radius=stop_radius,
        probe_candidates=int(args.probe_candidates),
        probe_radius=float(args.probe_radius),
        record_stride=max(1, int(args.record_stride)),
        dynamic_atlas=bool(args.dynamic_atlas),
        reference=bool(args.reference),
        weights=SolverWeights(
            defect=float(args.defect_weight),
            constraint=float(args.constraint_weight),
            atlas=float(args.atlas_weight),
            tidal=float(args.tidal_weight),
            pr_smooth=float(args.pr_smooth_weight),
        ),
    )
    meta.update({
        "mass": mass,
        "kappa": kappa,
        "energy": energy,
        "angular_momentum": ell,
        "r0": float(r0),
        "phi0": float(args.phi0),
        "radial_direction": "inward" if sign < 0 else "outward",
        "initial_pr": float(pr0),
    })
    return cfg, y0, params, meta


def perihelion_diagnostics(states: "np.ndarray") -> dict[str, Any]:
    ensure_numpy()
    if states.shape[0] < 5:
        return {"perihelia": [], "advance_mean": None}
    r = states[:, 1]
    phi = states[:, 2]
    idx: list[int] = []
    span = max(1.0e-12, float(np.max(r) - np.min(r)))
    min_separation = max(8, int(0.08 * len(r)))
    if r[1] > r[0]:
        idx.append(0)
    for i in range(1, len(r) - 1):
        if r[i] <= r[i - 1] and r[i] < r[i + 1] and (not idx or i - idx[-1] >= min_separation):
            if r[i] <= float(np.min(r)) + 0.15 * span:
                idx.append(i)
    if len(idx) >= 2 and phi[idx[-1]] - phi[idx[-2]] < math.pi:
        # Drop terminal numerical wiggles that are not a full radial cycle.
        idx.pop()
    cleaned: list[int] = []
    for i in idx:
        if not cleaned or phi[i] - phi[cleaned[-1]] >= math.pi:
            cleaned.append(i)
    idx = cleaned
    peri = [{"index": int(i), "r": float(r[i]), "phi": float(phi[i])} for i in idx]
    advances: list[float] = []
    if len(idx) >= 2:
        for a, b in zip(idx[:-1], idx[1:]):
            advances.append(float(phi[b] - phi[a] - 2.0 * math.pi))
    return {
        "perihelia": peri[:8],
        "advance_values": advances[:8],
        "advance_mean": (float(sum(advances) / len(advances)) if advances else None),
        "advance_count": int(len(advances)),
    }


def solve_case(args: argparse.Namespace, scenario: str) -> dict[str, Any]:
    ensure_numpy()
    cfg, y0, params, init_meta = scenario_defaults(args, scenario)
    t0 = now()
    states: list[np.ndarray] = [np.asarray(y0, dtype=float)]
    records: list[dict[str, Any]] = []
    metas: list[dict[str, Any]] = []
    chart_counts = {"schwarzschild": 0, "painleve_gullstrand": 0}
    horizon_cross_step: Optional[int] = None
    stop_reason = "max-steps"

    y = np.asarray(y0, dtype=float)
    for k in range(int(params.steps)):
        if params.stop_radius > 0 and float(y[1]) <= float(params.stop_radius):
            stop_reason = "stop-radius"
            break
        if float(y[1]) <= cfg.min_radius * 1.01:
            stop_reason = "min-radius"
            break
        y_next, meta = pandrosion_probe_step(y, params.step, cfg, params)
        if float(y[1]) > 2.0 * cfg.mass >= float(y_next[1]) and horizon_cross_step is None:
            horizon_cross_step = int(k + 1)
        chart_counts[str(meta.get("chart", "painleve_gullstrand"))] = chart_counts.get(str(meta.get("chart")), 0) + 1
        metas.append(meta)
        y = y_next
        states.append(y.copy())
        if (k + 1) % params.record_stride == 0 or k == 0:
            records.append({
                "step": int(k + 1),
                "lambda": float((k + 1) * params.step),
                "T": float(y[0]),
                "r": float(y[1]),
                "phi": float(y[2]),
                "p_r": float(y[3]),
                "H_error": float(constraint_error(y, cfg)),
                "chart": str(meta.get("chart")),
                "finite_slope_defect": float(meta.get("finite_slope_defect", float("nan"))),
                "tau": float(meta.get("tau", float("nan"))),
            })

    arr = np.asarray(states, dtype=float)
    h_errors = np.asarray([abs(constraint_error(row, cfg)) for row in arr], dtype=float)
    rhs_norm_errors = np.asarray([abs(2.0 * pg_hamiltonian(row, cfg) + cfg.kappa) for row in arr], dtype=float)
    r_values = arr[:, 1]
    phi_values = arr[:, 2]
    probe_defects = np.asarray([finite_float(m.get("finite_slope_defect"), 0.0) for m in metas], dtype=float)
    tau_sum = float(sum(finite_float(m.get("tau"), 0.0) for m in metas))

    reference_summary: dict[str, Any] = {}
    if params.reference:
        ref = np.asarray(y0, dtype=float)
        ref_states = [ref.copy()]
        for _ in range(len(arr) - 1):
            ref = rk4_step(ref, params.step, cfg, project=True)
            ref_states.append(ref.copy())
        ref_arr = np.asarray(ref_states, dtype=float)
        n = min(len(ref_arr), len(arr))
        reference_summary = {
            "enabled": True,
            "max_abs_r_difference": float(np.max(np.abs(arr[:n, 1] - ref_arr[:n, 1]))),
            "max_abs_phi_difference": float(np.max(np.abs(arr[:n, 2] - ref_arr[:n, 2]))),
            "final_r_reference": float(ref_arr[n - 1, 1]),
            "final_phi_reference": float(ref_arr[n - 1, 2]),
        }
    else:
        reference_summary = {"enabled": False}

    peri = perihelion_diagnostics(arr)
    final = arr[-1]
    summary = {
        "scenario": scenario,
        "success": bool(len(arr) > 2 and np.all(np.isfinite(arr[-1]))),
        "stop_reason": stop_reason,
        "steps_completed": int(len(arr) - 1),
        "lambda_final": float((len(arr) - 1) * params.step),
        "final_T": float(final[0]),
        "final_r": float(final[1]),
        "final_phi": float(final[2]),
        "final_p_r": float(final[3]),
        "min_r": float(np.min(r_values)),
        "max_r": float(np.max(r_values)),
        "delta_phi": float(phi_values[-1] - phi_values[0]),
        "horizon_crossed": bool(np.any(r_values <= 2.0 * cfg.mass)),
        "horizon_cross_step": horizon_cross_step,
        "max_abs_hamiltonian_error": float(np.max(h_errors)) if h_errors.size else None,
        "rms_hamiltonian_error": float(math.sqrt(float(np.mean(h_errors * h_errors)))) if h_errors.size else None,
        "max_abs_normalization_error": float(np.max(rhs_norm_errors)) if rhs_norm_errors.size else None,
        "mean_finite_slope_defect": float(np.mean(probe_defects)) if probe_defects.size else None,
        "max_finite_slope_defect": float(np.max(probe_defects)) if probe_defects.size else None,
        "proper_time_score_sum": tau_sum,
        "chart_counts": chart_counts,
        "perihelion": peri,
        "reference": reference_summary,
        "seconds": float(now() - t0),
    }

    return {
        "script": "307_full_schwarzschild_pandrosion_geodesic_solver.py",
        "autonomous": True,
        "dependencies": {"python_scripts": [], "numpy": bool(np is not None), "scipy": False},
        "mode": "307-full-schwarzschild-pandrosion-geodesic-solver",
        "papers": {
            "027": "Pandrosion Relativity: finite-slope causal probes and dynamic atlases",
            "028": "Pandrosion Relativity Benchmarks: Minkowski, Newtonian shadow, horizon atlas",
            "304": "Universal atlas autonomous NumPy style reused here for CLI/JSON discipline",
        },
        "flow_formula": "PG Hamiltonian state -> finite radial probe atlas around a Hamiltonian predictor -> Hamiltonian-shell projection -> proper-time + finite-slope-defect + dynamic-atlas score -> selected causal event; separate RK4 track is verification",
        "initial": init_meta,
        "parameters": {
            "step": float(params.step),
            "steps": int(params.steps),
            "stop_radius": float(params.stop_radius),
            "probe_candidates": int(params.probe_candidates),
            "probe_radius": float(params.probe_radius),
            "record_stride": int(params.record_stride),
            "dynamic_atlas": bool(params.dynamic_atlas),
            "weights": dataclasses.asdict(params.weights),
        },
        "summary": summary,
        "records": records,
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="307 standalone NumPy Full Schwarzschild Pandrosion Geodesic Solver")
    p.add_argument("--scenario", choices=["all", "radial-plunge", "precession", "custom"], default="all")
    p.add_argument("--mass", type=float, default=1.0)
    p.add_argument("--kappa", type=float, default=1.0, help="1 for timelike, 0 for null-style scoring")
    p.add_argument("--energy", type=float, default=None)
    p.add_argument("--angular-momentum", type=float, default=None)
    p.add_argument("--r0", type=float, default=None)
    p.add_argument("--phi0", type=float, default=0.0)
    p.add_argument("--radial-direction", choices=["inward", "outward"], default="inward")
    p.add_argument("--p", type=float, default=12.0, help="bound-orbit semi-latus rectum in units of M")
    p.add_argument("--eccentricity", type=float, default=0.25)
    p.add_argument("--step", type=float, default=0.0)
    p.add_argument("--steps", type=int, default=0)
    p.add_argument("--stop-radius", type=float, default=0.0)
    p.add_argument("--min-radius", type=float, default=1.0e-5)
    p.add_argument("--probe-candidates", type=int, default=9)
    p.add_argument("--probe-radius", type=float, default=0.75)
    p.add_argument("--defect-weight", type=float, default=1000.0)
    p.add_argument("--constraint-weight", type=float, default=25.0)
    p.add_argument("--atlas-weight", type=float, default=0.002)
    p.add_argument("--tidal-weight", type=float, default=0.0)
    p.add_argument("--pr-smooth-weight", type=float, default=0.0005)
    p.add_argument("--dynamic-atlas", action="store_true", default=True)
    p.add_argument("--no-dynamic-atlas", dest="dynamic_atlas", action="store_false")
    p.add_argument("--reference", action="store_true", default=True)
    p.add_argument("--no-reference", dest="reference", action="store_false")
    p.add_argument("--record-stride", type=int, default=20)
    p.add_argument("--self-test", action="store_true", help="run smaller radial+precession smoke tests")
    p.add_argument("--out", default=None)
    p.add_argument("--outdir", default="/private/tmp/307_schwarzschild_out")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)

    if bool(args.self_test):
        args.scenario = "all"
        args.record_stride = max(int(args.record_stride), 10)
        # Scenario defaults still apply, but keep the smoke test quick.
        if int(args.steps) <= 0:
            args.steps = 500
        if float(args.step) <= 0:
            args.step = 0.015
        if float(args.stop_radius) <= 0:
            args.stop_radius = 1.75

    scenarios = ["radial-plunge", "precession"] if args.scenario == "all" else [str(args.scenario)]
    outputs = [solve_case(args, sc) for sc in scenarios]
    final: dict[str, Any]
    if len(outputs) == 1:
        final = outputs[0]
    else:
        final = {
            "script": "307_full_schwarzschild_pandrosion_geodesic_solver.py",
            "autonomous": True,
            "mode": "307-full-schwarzschild-pandrosion-geodesic-solver",
            "cases": outputs,
            "summary": {
                "success": all(bool(o.get("summary", {}).get("success")) for o in outputs),
                "case_count": len(outputs),
                "scenarios": scenarios,
            },
        }

    if args.out:
        out = Path(args.out)
    else:
        name = "self_test_307.json" if bool(args.self_test) else f"307_{args.scenario}.json"
        out = Path(args.outdir) / name
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(json_safe(final), indent=2), encoding="utf-8")

    print("=" * 118, flush=True)
    print("307 autonomous FULL SCHWARZSCHILD PANDROSION GEODESIC SOLVER NumPy", flush=True)
    print("No SciPy and no imports from previous flow scripts; a separate RK4 track is used for verification.", flush=True)
    print("=" * 118, flush=True)
    for out_case in outputs:
        s = out_case["summary"]
        init = out_case["initial"]
        print(
            f"scenario={s['scenario']} E={init['energy']:.8g} L={init['angular_momentum']:.8g} "
            f"steps={s['steps_completed']} stop={s['stop_reason']}",
            flush=True,
        )
        print(
            f"final: r={s['final_r']:.8g} phi={s['final_phi']:.8g} "
            f"horizon_crossed={s['horizon_crossed']} chart_counts={s['chart_counts']}",
            flush=True,
        )
        print(
            f"errors: max|H+k/2|={s['max_abs_hamiltonian_error']:.3e} "
            f"max_defect={s['max_finite_slope_defect']:.3e} seconds={s['seconds']:.2f}",
            flush=True,
        )
        ref = s.get("reference", {})
        if ref.get("enabled"):
            print(
                f"reference: max_dr={ref['max_abs_r_difference']:.3e} "
                f"max_dphi={ref['max_abs_phi_difference']:.3e}",
                flush=True,
            )
        peri = s.get("perihelion", {})
        if peri.get("advance_mean") is not None:
            print(f"perihelion_advance_mean={float(peri['advance_mean']):.8g} rad", flush=True)
    print(f"out={out}", flush=True)
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
