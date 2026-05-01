"""Runtime wrapper around Pandrosion probe-aware correctors.

The bot now prefers the persistent Swift/Metal 129 server for the correction
step and falls back to the copied 118 Python engine when Swift is unavailable.
The copied engine filenames start with digits, so normal Python imports cannot
load them directly.
"""

from __future__ import annotations

import atexit
import importlib.util
import json
import math
import os
import subprocess
import sys
import threading
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence


def _finite_float(value: Any, default: float = 0.0) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return default
    return x if math.isfinite(x) else default


def _clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def _signed_cuberoot(value: float) -> float:
    if value == 0.0:
        return 0.0
    return math.copysign(abs(value) ** (1.0 / 3.0), value)


def _signed_power_root(value: float, degree: int) -> float:
    if value == 0.0:
        return 0.0
    return math.copysign(abs(value) ** (1.0 / max(1, int(degree))), value)


_SWIFT129_LAST_ERROR = ""


def _swift129_enabled() -> bool:
    disabled = os.getenv("PANDROSION_DISABLE_SWIFT_129", "").strip().lower()
    if disabled in {"1", "true", "yes", "on"}:
        return False
    engine = os.getenv("PANDROSION_CORE_ENGINE", "129").strip().lower()
    return engine in {"129", "swift129", "swift-129", "swift_metal", "swift-metal"}


def _record_swift129_error(message: str) -> None:
    global _SWIFT129_LAST_ERROR
    _SWIFT129_LAST_ERROR = message[:500]


@lru_cache(maxsize=1)
def _swift129_binary() -> Path | None:
    source = Path(__file__).resolve().with_name("129_pandrosion_swift_metal_full.swift")
    configured = os.getenv("PANDROSION_129_BINARY", "").strip()
    if configured:
        binary = Path(configured).expanduser()
        return binary if binary.exists() else None
    if not source.exists():
        _record_swift129_error(f"missing Swift 129 source: {source}")
        return None
    binary = Path(__file__).resolve().with_name("bin") / "129_pandrosion_swift_metal_full"
    try:
        needs_build = not binary.exists() or binary.stat().st_mtime < source.stat().st_mtime
        if needs_build:
            binary.parent.mkdir(parents=True, exist_ok=True)
            result = subprocess.run(
                [
                    "swiftc",
                    "-O",
                    str(source),
                    "-o",
                    str(binary),
                    "-framework",
                    "Foundation",
                    "-framework",
                    "Metal",
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=90.0,
            )
            if result.returncode != 0:
                _record_swift129_error(result.stderr or result.stdout or f"swiftc exited {result.returncode}")
                return None
    except Exception as exc:
        _record_swift129_error(f"Swift 129 build failed: {exc}")
        return None
    return binary if binary.exists() else None


class _Swift129Server:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._proc: subprocess.Popen[str] | None = None
        self._next_id = 1

    def close(self) -> None:
        proc = self._proc
        self._proc = None
        if proc is None:
            return
        try:
            proc.terminate()
            proc.wait(timeout=1.0)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

    def _ensure(self) -> subprocess.Popen[str]:
        proc = self._proc
        if proc is not None and proc.poll() is None and proc.stdin and proc.stdout:
            return proc
        binary = _swift129_binary()
        if binary is None:
            raise RuntimeError(_SWIFT129_LAST_ERROR or "Swift 129 binary is unavailable")
        proc = subprocess.Popen(
            [str(binary), "--server"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            bufsize=1,
        )
        self._proc = proc
        return proc

    def solve(self, payload: dict[str, Any]) -> dict[str, Any]:
        with self._lock:
            proc = self._ensure()
            job_id = self._next_id
            self._next_id += 1
            payload = dict(payload)
            payload["id"] = job_id
            line = json.dumps(payload, separators=(",", ":"), allow_nan=False)
            assert proc.stdin is not None and proc.stdout is not None
            try:
                proc.stdin.write(line + "\n")
                proc.stdin.flush()
                raw = proc.stdout.readline()
            except Exception:
                self.close()
                raise
            if not raw:
                code = proc.poll()
                self.close()
                raise RuntimeError(f"Swift 129 server exited before response, code={code}")
            response = json.loads(raw)
            if not response.get("ok", False):
                raise RuntimeError(str(response.get("error", "Swift 129 solve failed")))
            return response


_SWIFT129_SERVER = _Swift129Server()
atexit.register(_SWIFT129_SERVER.close)


def _swift129_corrector(
    engine: Any,
    *,
    n: int,
    degree: int,
    exps: Any,
    coeff: Any,
    y0: Any,
    max_epochs: int,
    tol: float,
    accept: float,
    line_search: int,
    probe_scale: float,
    direction_seed: int,
    probe_candidates: int,
    probe_radii: Sequence[float],
    include_self_probe: bool = True,
) -> dict[str, Any] | None:
    if not _swift129_enabled():
        return None
    try:
        np = engine.np
        exps_arr = np.asarray(exps, dtype=np.int64)
        if exps_arr.ndim == 1:
            exps_arr = exps_arr.reshape((-1, int(n)))
        coeff_arr = np.asarray(coeff, dtype=np.complex128)
        y0_arr = np.asarray(y0, dtype=np.complex128)
        if coeff_arr.shape[0] != int(n) or coeff_arr.shape[1] != exps_arr.shape[0]:
            raise ValueError("coefficient shape does not match exponent table")
        if y0_arr.shape[0] != int(n):
            raise ValueError("start point length does not match n")
        nz = np.any(np.abs(coeff_arr) > 0.0, axis=0)
        if not bool(np.any(nz)):
            return None
        exps_used = np.ascontiguousarray(exps_arr[nz], dtype=np.int64)
        coeff_used = np.ascontiguousarray(coeff_arr[:, nz], dtype=np.complex128)
        payload = {
            "op": "solve_custom",
            "n": int(n),
            "d": int(degree),
            "exps": [int(v) for v in exps_used.reshape(-1).tolist()],
            "coeff_re": [float(v) for v in coeff_used.real.reshape(-1).tolist()],
            "coeff_im": [float(v) for v in coeff_used.imag.reshape(-1).tolist()],
            "y0_re": [float(v) for v in y0_arr.real.reshape(-1).tolist()],
            "y0_im": [float(v) for v in y0_arr.imag.reshape(-1).tolist()],
            "direction_seed": int(direction_seed) & ((1 << 63) - 1),
            "max_epochs": int(max_epochs),
            "tol": float(tol),
            "accept": float(accept),
            "line_search": int(line_search),
            "probe_scale": float(probe_scale),
            "probe_candidates": int(probe_candidates),
            "probe_radii": [float(v) for v in probe_radii],
            "include_self_probe": bool(include_self_probe),
            "equation_normalize": False,
        }
        response = _SWIFT129_SERVER.solve(payload)
        y_re = response.get("y_re")
        y_im = response.get("y_im")
        if not isinstance(y_re, list) or not isinstance(y_im, list) or len(y_re) != int(n) or len(y_im) != int(n):
            raise ValueError("Swift 129 response has invalid y vector")
        y = np.asarray([complex(float(re), float(im)) for re, im in zip(y_re, y_im)], dtype=np.complex128)
        residual = _finite_float(response.get("residual"), float("inf"))
        return {
            "accepted": bool(response.get("accepted", False)),
            "status": response.get("status", "swift129"),
            "y": y,
            "residual": residual,
            "epochs": int(response.get("epochs", 0)),
            "seconds": _finite_float(response.get("seconds"), 0.0),
            "probe_name": "swift129",
            "probe_total_evals": int(response.get("probe_total_evals", 0)),
            "source": "129",
            "swift_terms_per_poly": int(response.get("terms_per_poly", 0)),
            "swift_total_terms": int(response.get("total_terms", 0)),
        }
    except Exception as exc:
        _record_swift129_error(str(exc))
        return None


@lru_cache(maxsize=1)
def _load_engine_118() -> Any:
    path = Path(__file__).resolve().with_name("118_pandrosion_probe_aware_pure_thales_engine.py")
    spec = importlib.util.spec_from_file_location("bot_pandrosion_118_probe_aware", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load Pandrosion 118 engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.ensure_numpy()
    return module


@lru_cache(maxsize=1)
def _market_consensus_target() -> Any:
    engine = _load_engine_118()
    np = engine.np
    n = 3
    d = 2
    exps = engine.compositions_leq(d, n)
    index = {tuple(int(v) for v in row): i for i, row in enumerate(exps.tolist())}
    coeff = np.zeros((n, exps.shape[0]), dtype=np.complex128)

    # f0 = x0^2 - 1, f1 = x1 - x0, f2 = x2 - x1.
    # The roots (+1,+1,+1) and (-1,-1,-1) are the two coherent trend basins.
    coeff[0, index[(0, 0, 0)]] = -1.0
    coeff[0, index[(2, 0, 0)]] = 1.0
    coeff[1, index[(1, 0, 0)]] = -1.0
    coeff[1, index[(0, 1, 0)]] = 1.0
    coeff[2, index[(0, 1, 0)]] = -1.0
    coeff[2, index[(0, 0, 1)]] = 1.0

    row_norm = np.linalg.norm(coeff, axis=1)
    row_norm = np.where(row_norm > 0, row_norm, 1.0)
    coeff = coeff / row_norm[:, None]

    system = engine.DenseKostlanSystem(
        n=n,
        d=d,
        seed=0x119B07,
        exps=exps,
        coeff=coeff,
        weights=np.ones(exps.shape[0], dtype=np.float64),
        equation_normalize=True,
    )
    chart = engine.LinearChart.identity(n)
    return engine.TargetTrack(system, chart)


@lru_cache(maxsize=16)
def _univariate_template(p: int) -> tuple[Any, Any, dict[int, int]]:
    engine = _load_engine_118()
    np = engine.np
    degree = int(max(1, p))
    exps = engine.compositions_leq(degree, 1)
    index = {int(row[0]): i for i, row in enumerate(exps.tolist())}
    weights = np.ones(exps.shape[0], dtype=np.float64)
    return engine, (exps, weights), index


def _univariate_start_terms_118(x_c: complex, p: int) -> tuple[complex, complex, complex]:
    """Deterministic first-order start for the 118 univariate solve.

    This is only a seed, not a fallback solver.  The returned signal is derived
    from the final 118 root in `solve_univariate_core_118`.
    """
    if x_c == 0j:
        return 0j, 0j, 0j
    if x_c == 1.0 + 0j:
        return 1.0 + 0j, 1.0 + 0j, 0j
    p_float = float(p)
    s_0 = 1.0 - (x_c - 1.0) / (x_c * p_float)
    return complex(s_0), complex(s_0), 0j


def _univariate_contraction_witness_118(x_c: complex, p: int, s_0: complex) -> complex:
    """Local 118 basin-rate witness for x*s**p - 1.

    The final endpoint still comes only from the 118 corrector.  This witness is
    used as a signed structural feature so the strategy sees contraction rate,
    not only raw root displacement.
    """
    if x_c == 0j or x_c == 1.0 + 0j:
        return 0j

    p_float = float(p)

    def S_p(s: complex) -> complex:
        if s == 1.0 + 0j:
            return complex(p_float, 0.0)
        return (1.0 - s**p) / (1.0 - s)

    def h(s: complex) -> complex:
        sp = S_p(s)
        if sp == 0j:
            return s
        return 1.0 - (x_c - 1.0) / (x_c * sp)

    try:
        s_1 = h(s_0)
        s_2 = h(s_1)
    except (ZeroDivisionError, OverflowError, ValueError):
        return 0j
    if not all(math.isfinite(z.real) and math.isfinite(z.imag) for z in (s_1, s_2)):
        return 0j
    diff = s_1 - s_0
    if abs(diff) <= 1e-12:
        return 0j
    return complex((s_2 - s_1) / diff)


def _principal_root_seed_118(x_c: complex, p: int) -> complex:
    """Principal analytic seed for multivariate 118 solves.

    This is only the starting chart point for the 5D/12D coupled systems; the
    accepted endpoint and quality still come from `pandrosion_corrector`.
    """
    if x_c == 0j:
        return 0j
    try:
        root = x_c ** (-1.0 / float(max(1, p)))
    except (ZeroDivisionError, OverflowError, ValueError):
        return _univariate_start_terms_118(x_c, p)[0]
    if math.isfinite(root.real) and math.isfinite(root.imag):
        return complex(root)
    return _univariate_start_terms_118(x_c, p)[0]


@lru_cache(maxsize=1)
def _multiscale_template_5d() -> tuple[Any, Any, dict[tuple[int, int, int, int, int], int]]:
    engine = _load_engine_118()
    np = engine.np
    n = 5
    degree = 3
    exps = engine.compositions_leq(degree, n)
    index = {tuple(int(v) for v in row): i for i, row in enumerate(exps.tolist())}
    weights = np.ones(exps.shape[0], dtype=np.float64)
    return engine, (exps, weights), index


@lru_cache(maxsize=1)
def _trade_decision_template_10d() -> tuple[Any, Any, dict[tuple[int, ...], int]]:
    engine = _load_engine_118()
    np = engine.np
    n = 10
    degree = 3
    exps = engine.compositions_leq(degree, n)
    index = {tuple(int(v) for v in row): i for i, row in enumerate(exps.tolist())}
    weights = np.ones(exps.shape[0], dtype=np.float64)
    return engine, (exps, weights), index


@lru_cache(maxsize=1)
def _trade_decision_template_12d() -> tuple[Any, Any, dict[tuple[int, ...], int]]:
    engine = _load_engine_118()
    np = engine.np
    n = 12
    degree = 3
    exps = engine.compositions_leq(degree, n)
    index = {tuple(int(v) for v in row): i for i, row in enumerate(exps.tolist())}
    weights = np.ones(exps.shape[0], dtype=np.float64)
    return engine, (exps, weights), index


@lru_cache(maxsize=1)
def _trade_management_template_8d() -> tuple[Any, Any, dict[tuple[int, ...], int]]:
    engine = _load_engine_118()
    np = engine.np
    n = 8
    degree = 3
    exps = engine.compositions_leq(degree, n)
    index = {tuple(int(v) for v in row): i for i, row in enumerate(exps.tolist())}
    weights = np.ones(exps.shape[0], dtype=np.float64)
    return engine, (exps, weights), index


@lru_cache(maxsize=1)
def _entry_allocation_template_6d() -> tuple[Any, Any, dict[tuple[int, ...], int]]:
    engine = _load_engine_118()
    np = engine.np
    n = 6
    degree = 3
    exps = engine.compositions_leq(degree, n)
    index = {tuple(int(v) for v in row): i for i, row in enumerate(exps.tolist())}
    weights = np.ones(exps.shape[0], dtype=np.float64)
    return engine, (exps, weights), index


@lru_cache(maxsize=1)
def _curve_regime_template_33d() -> tuple[Any, Any, dict[tuple[int, ...], int]]:
    engine = _load_engine_118()
    np = engine.np
    n = 4
    degree = 33
    exps = engine.compositions_leq(degree, n)
    index = {tuple(int(v) for v in row): i for i, row in enumerate(exps.tolist())}
    weights = np.ones(exps.shape[0], dtype=np.float64)
    return engine, (exps, weights), index


def _univariate_start_seed_118(x_c: complex, p: int) -> complex:
    return _univariate_start_terms_118(x_c, p)[0]


def solve_univariate_core_118(x_c: complex, p: int, direction_seed: int = 0) -> dict[str, Any]:
    """Use the 118 probe-aware finite-slope corrector for x*s**p - 1 = 0.

    The bot consumes a stable `(s0, s1, T2, T4, hat_lambda, gap)` contract.
    Here `T4` and `gap` are produced only by the 118 corrector; `hat_lambda`
    is a directional root-displacement witness (`1 - T4`) on the same scale
    used by strategy features.
    """
    x_c = complex(x_c)
    p = int(max(1, p))
    if x_c == 1.0 + 0j:
        return {
            "accepted": True,
            "status": "degenerate_identity",
            "s0": 1.0 + 0j,
            "s1": 1.0 + 0j,
            "T2": 1.0 + 0j,
            "T4": 1.0 + 0j,
            "hat_lambda": 0j,
            "gap": 0j,
            "residual": 0.0,
            "epochs": 0,
        }
    if x_c == 0j:
        return {
            "accepted": True,
            "status": "degenerate_zero",
            "s0": 0j,
            "s1": 0j,
            "T2": 0j,
            "T4": 0j,
            "hat_lambda": 0j,
            "gap": 0j,
            "residual": 0.0,
            "epochs": 0,
        }

    s_0, s_1, T_2 = _univariate_start_terms_118(x_c, p)
    engine, template, index = _univariate_template(p)
    exps, weights = template
    np = engine.np

    coeff = np.zeros((1, exps.shape[0]), dtype=np.complex128)
    coeff[0, index[0]] = -1.0
    coeff[0, index[p]] = x_c

    system = engine.DenseKostlanSystem(
        n=1,
        d=p,
        seed=0x118C0DE ^ int(direction_seed),
        exps=exps,
        coeff=coeff,
        weights=weights,
        equation_normalize=False,
    )
    target = engine.TargetTrack(system, engine.LinearChart.identity(1))
    y0 = np.asarray([s_0], dtype=np.complex128)

    loc = _swift129_corrector(
        engine,
        n=1,
        degree=p,
        exps=exps,
        coeff=coeff,
        y0=y0,
        max_epochs=4,
        tol=1e-12,
        accept=1e-10,
        line_search=8,
        probe_scale=0.05,
        direction_seed=int(direction_seed),
        probe_candidates=5,
        probe_radii=(0.0, 0.25, 0.5, 1.0, 1.8),
        include_self_probe=True,
    )
    if loc is None:
        loc = engine.pandrosion_corrector(
            target,
            y0,
            max_epochs=4,
            tol=1e-12,
            accept=1e-10,
            trial_timeout=0.0,
            line_search=8,
            probe_scale=0.05,
            direction_seed=int(direction_seed),
            probe_candidates=5,
            probe_radii=(0.0, 0.25, 0.5, 1.0, 1.8),
            include_self_probe=True,
        )

    root = complex(np.asarray(loc.get("y", y0), dtype=np.complex128)[0])
    residual = _finite_float(loc.get("residual"), float("inf"))
    accepted = bool(loc.get("accepted", False)) and math.isfinite(residual) and residual < 1e-8
    if not (math.isfinite(root.real) and math.isfinite(root.imag)):
        root = s_0
    hat_lambda = _univariate_contraction_witness_118(x_c, p, s_0)
    if hat_lambda == 0j:
        hat_lambda = 1.0 - root

    return {
        "accepted": bool(accepted),
        "status": loc.get("status"),
        "s0": s_0,
        "s1": s_1,
        "T2": root if accepted else T_2,
        "T4": root,
        "hat_lambda": hat_lambda,
        "gap": s_0 - root,
        "residual": float(residual),
        "epochs": int(loc.get("epochs", 0)),
        "probe_total_evals": int(loc.get("probe_total_evals", 0)),
    }


def solve_multiscale_consensus_118(
    ratios: Sequence[complex],
    liquidity_bias: float = 0.0,
    direction_seed: int = 0,
    coupling_strength: float = 0.08,
) -> dict[str, Any]:
    """Solve a coupled 5-timeframe cubic system with the 118 corrector.

    Variables are the five scale roots `(s3, s5, s10, s20, s38)`.  Each row keeps
    the local market equation `x_k*s_k**3 - 1 = 0`, while a Laplacian term couples
    neighbouring scales.  This is the multivariate problem the old local bot
    engine could not solve as one system.
    """
    if len(ratios) != 5:
        raise ValueError("ratios must contain exactly five complex scale ratios")

    engine, template, index = _multiscale_template_5d()
    exps, weights = template
    np = engine.np
    n = 5
    p = 3
    coupling = max(0.0, min(0.25, float(coupling_strength)))
    liq = max(-2.0, min(2.0, _finite_float(liquidity_bias))) * 0.015

    coeff = np.zeros((n, exps.shape[0]), dtype=np.complex128)
    zero = (0, 0, 0, 0, 0)

    def add(row: int, powers: tuple[int, int, int, int, int], value: complex) -> None:
        coeff[row, index[powers]] += complex(value)

    for i, raw_ratio in enumerate(ratios):
        ratio = complex(raw_ratio)
        cubic = [0, 0, 0, 0, 0]
        cubic[i] = p
        linear_i = [0, 0, 0, 0, 0]
        linear_i[i] = 1

        add(i, zero, -1.0)
        add(i, tuple(cubic), ratio)

        if coupling > 0.0:
            if i == 0:
                add(i, tuple(linear_i), coupling)
                neighbour = [0, 0, 0, 0, 0]
                neighbour[1] = 1
                add(i, tuple(neighbour), -coupling)
            elif i == n - 1:
                add(i, tuple(linear_i), coupling)
                neighbour = [0, 0, 0, 0, 0]
                neighbour[n - 2] = 1
                add(i, tuple(neighbour), -coupling)
            else:
                add(i, tuple(linear_i), 2.0 * coupling)
                left = [0, 0, 0, 0, 0]
                right = [0, 0, 0, 0, 0]
                left[i - 1] = 1
                right[i + 1] = 1
                add(i, tuple(left), -coupling)
                add(i, tuple(right), -coupling)

    if liq != 0.0:
        mid = [0, 0, 1, 0, 0]
        add(2, tuple(mid), 1j * liq)
        add(2, zero, -1j * liq)

    system = engine.DenseKostlanSystem(
        n=n,
        d=p,
        seed=0x1185CA1E ^ int(direction_seed),
        exps=exps,
        coeff=coeff,
        weights=weights,
        equation_normalize=False,
    )
    target = engine.TargetTrack(system, engine.LinearChart.identity(n))
    y0 = np.asarray([_principal_root_seed_118(complex(r), p) for r in ratios], dtype=np.complex128)
    initial_residual = target.residual(y0)

    loc = _swift129_corrector(
        engine,
        n=n,
        degree=p,
        exps=exps,
        coeff=coeff,
        y0=y0,
        max_epochs=4,
        tol=1e-10,
        accept=5e-8,
        line_search=8,
        probe_scale=0.045,
        direction_seed=int(direction_seed),
        probe_candidates=6,
        probe_radii=(0.0, 0.25, 0.5, 1.0, 1.8, 3.0),
        include_self_probe=True,
    )
    if loc is None:
        loc = engine.pandrosion_corrector(
            target,
            y0,
            max_epochs=4,
            tol=1e-10,
            accept=5e-8,
            trial_timeout=0.0,
            line_search=8,
            probe_scale=0.045,
            direction_seed=int(direction_seed),
            probe_candidates=6,
            probe_radii=(0.0, 0.25, 0.5, 1.0, 1.8, 3.0),
            include_self_probe=True,
        )

    y = np.asarray(loc.get("y", y0), dtype=np.complex128)
    residual = _finite_float(loc.get("residual"), float("inf"))
    accepted = bool(loc.get("accepted", False)) and math.isfinite(residual) and residual < 1e-6

    lambda_proxy = 1.0 - y.real
    abs_lambda = np.abs(lambda_proxy)
    mean_abs = float(np.mean(abs_lambda))
    direction_score = float(np.mean(lambda_proxy))
    direction = "BUY" if direction_score >= 0.0 else "SELL"
    consensus = abs(direction_score) / max(1e-12, mean_abs)
    displacement = float(np.linalg.norm(y - y0) / math.sqrt(n))
    tension = displacement / max(1e-8, mean_abs)

    residual_quality = math.exp(-min(50.0, 2.5e6 * max(0.0, residual)))
    tension_quality = math.exp(-min(6.0, 0.6 * max(0.0, tension)))
    quality = residual_quality * (0.35 + 0.65 * max(0.0, min(1.0, consensus))) * (0.55 + 0.45 * tension_quality)
    if not accepted:
        quality = 0.0

    return {
        "accepted": bool(accepted),
        "status": loc.get("status"),
        "direction": direction,
        "direction_score": direction_score,
        "quality": float(max(0.0, min(1.0, quality))),
        "consensus": float(max(0.0, min(1.0, consensus))),
        "initial_residual": float(initial_residual),
        "residual": float(residual),
        "displacement": float(displacement),
        "tension": float(max(0.0, tension)),
        "epochs": int(loc.get("epochs", 0)),
        "probe_name": loc.get("probe_name"),
        "probe_total_evals": int(loc.get("probe_total_evals", 0)),
    }


def solve_trade_decision_system_118(
    ratios: Sequence[complex],
    *,
    volatility_bps: float = 5.0,
    volatility_expansion_ratio: float = 1.0,
    momentum_3_bps: float = 0.0,
    momentum_10_bps: float = 0.0,
    z_score: float = 0.0,
    z_score_20: float = 0.0,
    acceleration_bps: float = 0.0,
    liquidity_momentum_bps: float = 0.0,
    stretch_bps: float = 0.0,
    swing_compression: float = 0.0,
    lambda_consensus: float = 0.5,
    margin_pressure: float = 0.0,
    funding_rate: float = 0.0,
    global_ema_trend: float = 0.0,
    global_rsi: float | None = None,
    direction_seed: int = 0,
) -> dict[str, Any]:
    """Solve a coupled 12D trade-decision system with the 118 corrector.

    Variables:
      d, s3, s5, s10, s20, s38, liquidity, volatility, stretch, risk,
      runner, capital

    The five `s*` variables retain the exact market equations
    `x_k*s_k**3 - 1 = 0`; the additional variables couple direction, liquidity,
    volatility expansion, mean stretch, margin/risk pressure, runner quality,
    and capital efficiency into one square polynomial system.  The old local bot
    engine could not express this as one solve.
    """
    if len(ratios) != 5:
        raise ValueError("ratios must contain exactly five complex scale ratios")

    engine, template, index = _trade_decision_template_12d()
    exps, weights = template
    np = engine.np
    n = 12
    p = 3
    zero = (0,) * n

    def unit(axis: int, power: int = 1) -> tuple[int, ...]:
        powers = [0] * n
        powers[axis] = power
        return tuple(powers)

    def mixed(a: int, b: int) -> tuple[int, ...]:
        powers = [0] * n
        powers[a] += 1
        powers[b] += 1
        return tuple(powers)

    def add(row: int, powers: tuple[int, ...], value: complex) -> None:
        coeff[row, index[powers]] += complex(value)

    vol = max(0.1, _finite_float(volatility_bps, 5.0))
    mom3_norm = _clamp(_finite_float(momentum_3_bps) / vol, -6.0, 6.0)
    mom10_norm = _clamp(_finite_float(momentum_10_bps) / vol, -6.0, 6.0)
    z_norm = _clamp(_finite_float(z_score), -6.0, 6.0)
    z20_norm = _clamp(_finite_float(z_score_20), -6.0, 6.0)
    accel_norm = _clamp(_finite_float(acceleration_bps) / vol, -6.0, 6.0)
    liq_norm = _clamp(_finite_float(liquidity_momentum_bps) / vol, -6.0, 6.0)
    stretch_norm = _clamp(_finite_float(stretch_bps) / max(1.0, vol * 6.0), -6.0, 6.0)
    expansion_norm = _clamp(_finite_float(volatility_expansion_ratio, 1.0) - 1.0, -3.0, 3.0)
    consensus = _clamp(_finite_float(lambda_consensus, 0.5), 0.0, 1.0)
    compression = _clamp(_finite_float(swing_compression), 0.0, 1.0)
    margin = _clamp(_finite_float(margin_pressure), 0.0, 1.5)
    funding_norm = _clamp(_finite_float(funding_rate) * -2500.0, -2.0, 2.0)
    ema_norm = _clamp(_finite_float(global_ema_trend) / 2.0, -3.0, 3.0)
    rsi = _finite_float(global_rsi, 50.0)
    rsi_extreme = _clamp(abs(rsi - 50.0) / 50.0, 0.0, 1.0)

    direction_anchor = math.tanh(
        0.30 * z_norm
        + 0.22 * z20_norm
        + 0.24 * mom3_norm
        + 0.18 * mom10_norm
        + 0.13 * accel_norm
        + 0.06 * liq_norm
        + 0.05 * ema_norm
        + 0.03 * funding_norm
    )
    liquidity_anchor = math.tanh(0.70 * liq_norm)
    volatility_anchor = math.tanh(1.35 * expansion_norm)
    stretch_anchor = math.tanh(stretch_norm)
    risk_anchor = _clamp(
        0.45 * min(1.0, margin)
        + 0.18 * max(0.0, abs(z_norm) - 2.1) / 3.9
        + 0.14 * min(1.0, abs(stretch_norm) / 2.5)
        + 0.13 * (1.0 - consensus)
        + 0.08 * rsi_extreme
        + 0.06 * compression,
        0.0,
        0.96,
    )
    directional_strength = _clamp(
        0.32 * abs(direction_anchor)
        + 0.16 * min(1.0, abs(mom3_norm) / 3.0)
        + 0.12 * min(1.0, abs(mom10_norm) / 3.0)
        + 0.10 * min(1.0, abs(accel_norm) / 3.0),
        0.0,
        1.0,
    )
    runner_anchor = _clamp(
        0.30 * directional_strength
        + 0.22 * consensus
        + 0.18 * (1.0 - risk_anchor)
        + 0.12 * compression
        + 0.10 * max(0.0, liquidity_anchor)
        + 0.08 * max(0.0, 1.0 - abs(expansion_norm) / 2.5)
        - 0.10 * min(1.0, abs(stretch_norm) / 3.0),
        0.0,
        1.0,
    )
    capital_anchor = _clamp(
        0.30 * (1.0 - risk_anchor)
        + 0.22 * consensus
        + 0.18 * directional_strength
        + 0.14 * (1.0 - min(1.0, margin))
        + 0.10 * max(0.0, liquidity_anchor)
        + 0.06 * (1.0 - rsi_extreme),
        0.0,
        1.0,
    )

    coeff = np.zeros((n, exps.shape[0]), dtype=np.complex128)
    alpha = 0.16
    direction_coupling = 0.08
    root_coupling = 0.025
    root_gain = _clamp(30000.0 / max(18.0, vol), 80.0, 650.0)
    scale_weights = (0.23, 0.22, 0.23, 0.18, 0.14)

    # f0: latent direction.  It is anchored by observable normalized features
    # and coupled back to the five exact ratio roots through root_gain*(1-s_k).
    add(0, unit(0, 3), 1.0)
    add(0, unit(0), alpha + direction_coupling)
    add(0, zero, -direction_anchor - direction_coupling * root_gain * sum(scale_weights))
    for axis, weight in zip(range(1, 6), scale_weights):
        add(0, unit(axis), direction_coupling * root_gain * weight)
    add(0, unit(6), -0.045)
    add(0, unit(8), 0.060)
    add(0, unit(10), -0.035)
    add(0, unit(11), -0.025)
    add(0, mixed(0, 7), -0.040)
    add(0, mixed(0, 9), 0.180)
    add(0, mixed(0, 10), -0.030)
    add(0, mixed(0, 11), -0.020)

    # f1..f5: exact cubic scale equations, weakly coupled to direction.
    for row, raw_ratio in enumerate(ratios, start=1):
        ratio = complex(raw_ratio)
        add(row, unit(row, p), ratio)
        add(row, zero, -1.0 - root_coupling)
        add(row, unit(row), root_coupling)
        add(row, unit(0), root_coupling / root_gain)

    # f6..f11: liquidity, volatility expansion, stretch, risk pressure,
    # runner quality, and capital efficiency.
    anchors = (
        liquidity_anchor,
        volatility_anchor,
        stretch_anchor,
        risk_anchor,
        runner_anchor,
        capital_anchor,
    )
    for offset, anchor in enumerate(anchors, start=6):
        add(offset, unit(offset, 3), 1.0)
        add(offset, unit(offset), alpha)
        add(offset, zero, -anchor)
    add(9, mixed(0, 0), 0.030 * risk_anchor)
    add(10, mixed(0, 10), -0.045)
    add(10, unit(9), 0.070)
    add(11, unit(10), -0.050)
    add(11, unit(9), 0.060)

    system = engine.DenseKostlanSystem(
        n=n,
        d=p,
        seed=0x118DEC15 ^ int(direction_seed),
        exps=exps,
        coeff=coeff,
        weights=weights,
        equation_normalize=False,
    )
    target = engine.TargetTrack(system, engine.LinearChart.identity(n))

    y0 = np.zeros(n, dtype=np.complex128)
    y0[0] = complex(_signed_cuberoot(direction_anchor), 0.0)
    for idx, raw_ratio in enumerate(ratios, start=1):
        y0[idx] = _principal_root_seed_118(complex(raw_ratio), p)
    y0[6] = complex(_signed_cuberoot(liquidity_anchor), 0.0)
    y0[7] = complex(_signed_cuberoot(volatility_anchor), 0.0)
    y0[8] = complex(_signed_cuberoot(stretch_anchor), 0.0)
    y0[9] = complex(_signed_cuberoot(risk_anchor), 0.0)
    y0[10] = complex(_signed_cuberoot(runner_anchor), 0.0)
    y0[11] = complex(_signed_cuberoot(capital_anchor), 0.0)
    initial_residual = target.residual(y0)

    loc = _swift129_corrector(
        engine,
        n=n,
        degree=p,
        exps=exps,
        coeff=coeff,
        y0=y0,
        max_epochs=5,
        tol=1e-10,
        accept=2e-7,
        line_search=8,
        probe_scale=0.040,
        direction_seed=int(direction_seed),
        probe_candidates=6,
        probe_radii=(0.0, 0.18, 0.35, 0.7, 1.2, 2.0),
        include_self_probe=True,
    )
    if loc is None:
        loc = engine.pandrosion_corrector(
            target,
            y0,
            max_epochs=5,
            tol=1e-10,
            accept=2e-7,
            trial_timeout=0.0,
            line_search=8,
            probe_scale=0.040,
            direction_seed=int(direction_seed),
            probe_candidates=6,
            probe_radii=(0.0, 0.18, 0.35, 0.7, 1.2, 2.0),
            include_self_probe=True,
        )

    y = np.asarray(loc.get("y", y0), dtype=np.complex128)
    residual = _finite_float(loc.get("residual"), float("inf"))
    accepted = bool(loc.get("accepted", False)) and math.isfinite(residual) and residual < 2e-5

    root_signal = float(root_gain * sum(w * (1.0 - y[i].real) for i, w in zip(range(1, 6), scale_weights)))
    direction_score = float(_clamp(0.55 * y[0].real + 0.45 * math.tanh(root_signal), -1.0, 1.0))
    direction = "BUY" if direction_score >= 0.0 else "SELL"
    sign = 1.0 if direction == "BUY" else -1.0
    directional_features = (mom3_norm, mom10_norm, z_norm, z20_norm, accel_norm, liq_norm)
    aligned = sum(1 for value in directional_features if value * sign >= 0.0)
    alignment = aligned / float(len(directional_features))
    risk_pressure = float(_clamp(abs(y[9].real), 0.0, 1.0))
    runner_score = float(_clamp(y[10].real, 0.0, 1.0))
    capital_efficiency = float(_clamp(y[11].real, 0.0, 1.0))
    residual_quality = math.exp(-min(50.0, 2.0e5 * max(0.0, residual)))
    basin_quality = max(0.0, min(1.0, abs(direction_score)))
    risk_quality = 1.0 - 0.55 * risk_pressure
    runner_quality = 0.72 + 0.28 * runner_score
    capital_quality = 0.82 + 0.18 * capital_efficiency
    quality = (
        residual_quality
        * (0.25 + 0.75 * basin_quality)
        * (0.45 + 0.55 * alignment)
        * risk_quality
        * runner_quality
        * capital_quality
    )
    if not accepted:
        quality = 0.0

    edge_score = quality * (2.0 + 5.2 * abs(direction_score) + 1.1 * consensus + 0.8 * runner_score)
    size_scale = _clamp(
        (0.38 + 0.86 * quality + 0.36 * capital_efficiency + 0.22 * runner_score)
        * (1.0 - 0.48 * risk_pressure),
        0.10,
        1.45,
    )
    exit_bias = _clamp(
        0.55 * risk_pressure + 0.25 * (1.0 - alignment) + 0.20 * min(1.0, margin) - 0.18 * runner_score,
        0.0,
        1.0,
    )
    late_entry_risk = _clamp(
        0.55 * risk_pressure
        + 0.35 * min(1.0, abs(stretch_norm) / 3.0)
        + 0.10 * max(0.0, expansion_norm)
        - 0.15 * runner_score,
        0.0,
        1.0,
    )

    return {
        "accepted": bool(accepted),
        "status": loc.get("status"),
        "direction": direction,
        "direction_score": direction_score,
        "quality": float(_clamp(quality, 0.0, 1.0)),
        "consensus": float(_clamp(consensus, 0.0, 1.0)),
        "alignment": float(_clamp(alignment, 0.0, 1.0)),
        "edge_score": float(max(0.0, edge_score)),
        "size_scale": float(size_scale),
        "exit_bias": float(exit_bias),
        "risk_pressure": float(risk_pressure),
        "runner_score": float(runner_score),
        "capital_efficiency": float(capital_efficiency),
        "late_entry_risk": float(late_entry_risk),
        "margin_pressure": float(min(1.0, margin)),
        "root_signal": float(root_signal),
        "initial_residual": float(initial_residual),
        "residual": float(residual),
        "epochs": int(loc.get("epochs", 0)),
        "probe_name": loc.get("probe_name"),
        "probe_total_evals": int(loc.get("probe_total_evals", 0)),
    }


@lru_cache(maxsize=8192)
def solve_trade_management_system_118(
    *,
    runner_score: float = 0.0,
    exit_bias: float = 0.0,
    risk_pressure: float = 0.0,
    late_entry_risk: float = 0.0,
    capital_efficiency: float = 0.0,
    decision_quality: float = 0.0,
    multiscale_quality: float = 0.0,
    multiscale_consensus: float = 0.0,
    opposition: float = 0.0,
    peak_profit_pct: float = 0.0,
    volatility_bps: float = 5.0,
    fee_round_trip_pct: float = 0.0007,
    direction_seed: int = 0,
) -> dict[str, Any]:
    """Solve an 8D 118 trade-management system.

    Variables are protect, retain, activation, urgency, breathe, stop budget,
    consistency, and capital preservation.  This is deliberately separate from
    entry selection: it solves the exit-management geometry for an open trade
    using the same probe-aware PURE Thales corrector loaded from engine 118.
    """
    engine, template, index = _trade_management_template_8d()
    exps, weights = template
    np = engine.np
    n = 8
    p = 3
    zero = (0,) * n

    def unit(axis: int, power: int = 1) -> tuple[int, ...]:
        powers = [0] * n
        powers[axis] = power
        return tuple(powers)

    def mixed(a: int, b: int) -> tuple[int, ...]:
        powers = [0] * n
        powers[a] = 1
        powers[b] = 1
        return tuple(powers)

    def add(row: int, powers: tuple[int, ...], value: complex) -> None:
        coeff[row, index[powers]] += complex(value)

    runner = _clamp(_finite_float(runner_score), 0.0, 1.0)
    exit_p = _clamp(_finite_float(exit_bias), 0.0, 1.0)
    risk_p = _clamp(_finite_float(risk_pressure), 0.0, 1.0)
    late = _clamp(_finite_float(late_entry_risk), 0.0, 1.0)
    capital = _clamp(_finite_float(capital_efficiency), 0.0, 1.0)
    decision_q = _clamp(_finite_float(decision_quality), 0.0, 1.0)
    ms_q = _clamp(_finite_float(multiscale_quality), 0.0, 1.0)
    ms_consensus = _clamp(_finite_float(multiscale_consensus), 0.0, 1.0)
    oppose = _clamp(_finite_float(opposition), 0.0, 1.0)
    peak_units = _clamp(_finite_float(peak_profit_pct) / max(1e-6, _finite_float(fee_round_trip_pct, 0.0007)), 0.0, 12.0)
    vol_norm = _clamp(_finite_float(volatility_bps, 5.0) / 35.0, 0.0, 3.0)

    quality_anchor = _clamp(0.48 * decision_q + 0.34 * ms_q * ms_consensus + 0.18 * capital, 0.0, 1.0)
    urgency_anchor = _clamp(
        0.36 * risk_p
        + 0.28 * exit_p
        + 0.16 * late
        + 0.16 * oppose
        + 0.08 * min(1.0, vol_norm)
        - 0.14 * runner,
        0.0,
        1.0,
    )
    protect_anchor = _clamp(
        0.18
        + 0.34 * urgency_anchor
        + 0.20 * oppose
        + 0.12 * min(1.0, peak_units / 2.5)
        - 0.14 * runner
        + 0.08 * quality_anchor,
        0.04,
        0.94,
    )
    retain_anchor = _clamp(
        0.22
        + 0.30 * protect_anchor
        + 0.18 * urgency_anchor
        + 0.14 * oppose
        + 0.10 * capital
        - 0.12 * runner,
        0.10,
        0.88,
    )
    activation_anchor = _clamp(
        0.92
        + 0.20 * runner
        + 0.10 * capital
        - 0.30 * urgency_anchor
        - 0.16 * oppose,
        0.36,
        1.28,
    )
    breathe_anchor = _clamp(0.18 + 0.62 * runner + 0.16 * quality_anchor - 0.22 * urgency_anchor, 0.0, 1.0)
    stop_anchor = _clamp(
        0.16
        + 0.30 * runner
        + 0.18 * capital
        + 0.10 * quality_anchor
        - 0.34 * urgency_anchor
        - 0.22 * oppose,
        0.06,
        0.72,
    )
    consistency_anchor = _clamp(quality_anchor * (1.0 - 0.35 * urgency_anchor) + 0.15 * runner, 0.0, 1.0)

    coeff = np.zeros((n, exps.shape[0]), dtype=np.complex128)
    alpha = 0.15
    anchors = (
        protect_anchor,
        retain_anchor,
        activation_anchor,
        urgency_anchor,
        breathe_anchor,
        stop_anchor,
        consistency_anchor,
        capital,
    )
    for row, anchor in enumerate(anchors):
        add(row, unit(row, p), 1.0)
        add(row, unit(row), alpha)
        add(row, zero, -anchor)

    # Couplings: protection and retention respond to urgency/opposition;
    # breathe and stop budget oppose urgency; consistency ties the solve together.
    add(0, unit(3), 0.08)
    add(0, unit(4), -0.06)
    add(0, mixed(1, 3), 0.035)
    add(1, unit(0), -0.05)
    add(1, unit(4), -0.035)
    add(1, mixed(3, 6), 0.045)
    add(2, unit(3), -0.055)
    add(2, unit(4), 0.040)
    add(3, unit(0), 0.045)
    add(3, unit(6), -0.040)
    add(4, unit(3), -0.080)
    add(4, unit(7), 0.045)
    add(5, unit(0), -0.075)
    add(5, unit(3), -0.090)
    add(5, unit(4), 0.065)
    add(6, unit(3), -0.045)
    add(6, mixed(4, 7), 0.035)
    add(7, unit(3), -0.050)
    add(7, unit(6), 0.035)

    system = engine.DenseKostlanSystem(
        n=n,
        d=p,
        seed=0x118E117 ^ int(direction_seed),
        exps=exps,
        coeff=coeff,
        weights=weights,
        equation_normalize=False,
    )
    target = engine.TargetTrack(system, engine.LinearChart.identity(n))

    y0 = np.asarray([complex(_signed_cuberoot(anchor), 0.0) for anchor in anchors], dtype=np.complex128)
    initial_residual = target.residual(y0)

    loc = _swift129_corrector(
        engine,
        n=n,
        degree=p,
        exps=exps,
        coeff=coeff,
        y0=y0,
        max_epochs=5,
        tol=1e-10,
        accept=2e-7,
        line_search=8,
        probe_scale=0.035,
        direction_seed=int(direction_seed),
        probe_candidates=6,
        probe_radii=(0.0, 0.16, 0.32, 0.64, 1.1, 1.8),
        include_self_probe=True,
    )
    if loc is None:
        loc = engine.pandrosion_corrector(
            target,
            y0,
            max_epochs=5,
            tol=1e-10,
            accept=2e-7,
            trial_timeout=0.0,
            line_search=8,
            probe_scale=0.035,
            direction_seed=int(direction_seed),
            probe_candidates=6,
            probe_radii=(0.0, 0.16, 0.32, 0.64, 1.1, 1.8),
            include_self_probe=True,
        )

    y = np.asarray(loc.get("y", y0), dtype=np.complex128)
    residual = _finite_float(loc.get("residual"), float("inf"))
    accepted = bool(loc.get("accepted", False)) and math.isfinite(residual) and residual < 2e-5

    protect = _clamp(y[0].real, 0.0, 1.0)
    retain = _clamp(y[1].real, 0.0, 1.0)
    activation = _clamp(y[2].real, 0.30, 1.35)
    urgency = _clamp(y[3].real, 0.0, 1.0)
    breathe = _clamp(y[4].real, 0.0, 1.0)
    stop_budget = _clamp(y[5].real, 0.04, 0.82)
    consistency = _clamp(y[6].real, 0.0, 1.0)
    capital_guard = _clamp(y[7].real, 0.0, 1.0)

    residual_quality = math.exp(-min(50.0, 2.0e5 * max(0.0, residual)))
    shape_quality = (
        0.40
        + 0.24 * consistency
        + 0.16 * capital_guard
        + 0.12 * breathe
        + 0.08 * (1.0 - urgency)
    )
    quality = residual_quality * _clamp(shape_quality, 0.0, 1.0)
    if not accepted:
        quality = 0.0

    return {
        "accepted": bool(accepted),
        "status": loc.get("status"),
        "quality": float(_clamp(quality, 0.0, 1.0)),
        "protect_tightness": float(protect),
        "retain_fraction": float(retain),
        "activation_scale": float(activation),
        "urgency": float(urgency),
        "breathe_score": float(breathe),
        "stop_budget_scale": float(stop_budget),
        "consistency": float(consistency),
        "capital_guard": float(capital_guard),
        "initial_residual": float(initial_residual),
        "residual": float(residual),
        "epochs": int(loc.get("epochs", 0)),
        "probe_name": loc.get("probe_name"),
        "probe_total_evals": int(loc.get("probe_total_evals", 0)),
    }


@lru_cache(maxsize=8192)
def solve_entry_allocation_system_118(
    *,
    edge_to_floor_ratio: float = 1.0,
    signed_chase_bps: float = 0.0,
    volatility_bps: float = 5.0,
    swing_compression: float = 0.0,
    multiscale_tension: float = 0.0,
    risk_pressure: float = 0.0,
    late_entry_risk: float = 0.0,
    runner_score: float = 0.0,
    capital_efficiency: float = 0.0,
    decision_quality: float = 0.0,
    margin_pressure: float = 0.0,
    consecutive_losses: int = 0,
    direction_seed: int = 0,
) -> dict[str, Any]:
    """Solve a 6D 118 entry-allocation system.

    Variables are conviction, fragility, allocation, chase guard, cooldown
    pressure, and capital throttle.  The solve turns the 12D decision outputs
    plus live path-shape data into one coherent admission/sizing surface.
    """
    engine, template, index = _entry_allocation_template_6d()
    exps, weights = template
    np = engine.np
    n = 6
    p = 3
    zero = (0,) * n

    def unit(axis: int, power: int = 1) -> tuple[int, ...]:
        powers = [0] * n
        powers[axis] = power
        return tuple(powers)

    def mixed(a: int, b: int) -> tuple[int, ...]:
        powers = [0] * n
        powers[a] = 1
        powers[b] = 1
        return tuple(powers)

    def add(row: int, powers: tuple[int, ...], value: complex) -> None:
        coeff[row, index[powers]] += complex(value)

    edge_ratio = _clamp(_finite_float(edge_to_floor_ratio, 1.0), 0.0, 12.0)
    vol = max(0.1, _finite_float(volatility_bps, 5.0))
    chase_norm = _clamp(_finite_float(signed_chase_bps) / max(10.0, vol * 0.65), -4.0, 4.0)
    chase_pressure = _clamp(max(0.0, chase_norm) / 3.0, 0.0, 1.0)
    compression = _clamp(_finite_float(swing_compression), 0.0, 1.0)
    tension_pressure = _clamp(_finite_float(multiscale_tension) / 0.07, 0.0, 1.0)
    risk_p = _clamp(_finite_float(risk_pressure), 0.0, 1.0)
    late = _clamp(_finite_float(late_entry_risk), 0.0, 1.0)
    runner = _clamp(_finite_float(runner_score), 0.0, 1.0)
    capital = _clamp(_finite_float(capital_efficiency), 0.0, 1.0)
    decision_q = _clamp(_finite_float(decision_quality), 0.0, 1.0)
    margin = _clamp(_finite_float(margin_pressure), 0.0, 1.5)
    loss_pressure = _clamp(float(max(0, int(consecutive_losses))) / 4.0, 0.0, 1.0)
    edge_strength = _clamp((edge_ratio - 5.0) / 3.0, 0.0, 1.0)

    fragility_anchor = _clamp(
        0.28 * risk_p
        + 0.22 * compression
        + 0.18 * late
        + 0.16 * chase_pressure
        + 0.18 * tension_pressure
        + 0.12 * loss_pressure
        + 0.08 * min(1.0, margin)
        - 0.14 * runner
        - 0.08 * capital,
        0.0,
        1.0,
    )
    conviction_anchor = _clamp(
        0.30 * decision_q
        + 0.22 * runner
        + 0.18 * capital
        + 0.17 * edge_strength
        + 0.13 * (1.0 - risk_p),
        0.0,
        1.0,
    )
    chase_guard_anchor = _clamp(
        0.10
        + 0.54 * chase_pressure
        + 0.16 * late
        + 0.12 * risk_p
        + 0.10 * compression
        + 0.08 * tension_pressure
        - 0.16 * runner,
        0.0,
        1.0,
    )
    cooldown_anchor = _clamp(
        0.16
        + 0.28 * loss_pressure
        + 0.22 * risk_p
        + 0.16 * late
        + 0.10 * tension_pressure
        + 0.18 * fragility_anchor,
        0.0,
        1.0,
    )
    capital_throttle_anchor = _clamp(
        0.18
        + 0.38 * capital
        + 0.24 * (1.0 - risk_p)
        + 0.14 * runner
        + 0.06 * edge_strength
        - 0.30 * fragility_anchor
        - 0.10 * chase_guard_anchor
        - 0.12 * tension_pressure,
        0.05,
        1.0,
    )
    allocation_anchor = _clamp(
        0.20
        + 0.92 * conviction_anchor
        + 0.14 * edge_strength
        + 0.10 * capital_throttle_anchor
        - 0.72 * fragility_anchor
        - 0.18 * chase_guard_anchor
        - 0.14 * cooldown_anchor,
        0.04,
        1.40,
    )

    coeff = np.zeros((n, exps.shape[0]), dtype=np.complex128)
    alpha = 0.14
    anchors = (
        conviction_anchor,
        fragility_anchor,
        allocation_anchor,
        chase_guard_anchor,
        cooldown_anchor,
        capital_throttle_anchor,
    )
    for row, anchor in enumerate(anchors):
        add(row, unit(row, p), 1.0)
        add(row, unit(row), alpha)
        add(row, zero, -anchor)

    # Coupled admission geometry: fragility and chase suppress allocation;
    # conviction and capital throttle stabilize it, while cooldown pressure
    # reacts to both fragility and recent loss state.
    add(0, unit(1), 0.045)
    add(0, unit(5), -0.035)
    add(1, unit(0), -0.055)
    add(1, unit(3), 0.060)
    add(1, mixed(3, 4), 0.040)
    add(2, unit(0), -0.090)
    add(2, unit(1), 0.120)
    add(2, unit(3), 0.075)
    add(2, unit(4), 0.065)
    add(2, unit(5), -0.080)
    add(3, unit(0), -0.035)
    add(3, unit(1), 0.055)
    add(4, unit(1), 0.060)
    add(4, unit(5), -0.040)
    add(5, unit(1), 0.050)
    add(5, unit(0), -0.045)
    add(5, mixed(0, 2), -0.030)

    system = engine.DenseKostlanSystem(
        n=n,
        d=p,
        seed=0x118A110C ^ int(direction_seed),
        exps=exps,
        coeff=coeff,
        weights=weights,
        equation_normalize=False,
    )
    target = engine.TargetTrack(system, engine.LinearChart.identity(n))

    y0 = np.asarray([complex(_signed_cuberoot(anchor), 0.0) for anchor in anchors], dtype=np.complex128)
    initial_residual = target.residual(y0)

    loc = _swift129_corrector(
        engine,
        n=n,
        degree=p,
        exps=exps,
        coeff=coeff,
        y0=y0,
        max_epochs=5,
        tol=1e-10,
        accept=2e-7,
        line_search=8,
        probe_scale=0.035,
        direction_seed=int(direction_seed),
        probe_candidates=6,
        probe_radii=(0.0, 0.16, 0.32, 0.64, 1.1, 1.8),
        include_self_probe=True,
    )
    if loc is None:
        loc = engine.pandrosion_corrector(
            target,
            y0,
            max_epochs=5,
            tol=1e-10,
            accept=2e-7,
            trial_timeout=0.0,
            line_search=8,
            probe_scale=0.035,
            direction_seed=int(direction_seed),
            probe_candidates=6,
            probe_radii=(0.0, 0.16, 0.32, 0.64, 1.1, 1.8),
            include_self_probe=True,
        )

    y = np.asarray(loc.get("y", y0), dtype=np.complex128)
    residual = _finite_float(loc.get("residual"), float("inf"))
    accepted = bool(loc.get("accepted", False)) and math.isfinite(residual) and residual < 2e-5

    conviction = _clamp(y[0].real, 0.0, 1.0)
    fragility = _clamp(y[1].real, 0.0, 1.0)
    allocation = _clamp(y[2].real, 0.0, 1.55)
    chase_guard = _clamp(y[3].real, 0.0, 1.0)
    cooldown_pressure = _clamp(y[4].real, 0.0, 1.0)
    capital_throttle = _clamp(y[5].real, 0.0, 1.0)

    entry_gate = _clamp(
        conviction
        * (1.0 - 0.52 * fragility)
        * (1.0 - 0.34 * chase_guard)
        * (0.62 + 0.38 * capital_throttle)
        * (1.0 - 0.22 * cooldown_pressure),
        0.0,
        1.0,
    )
    allocation_scale = _clamp(
        0.18
        + 1.58 * allocation * (0.78 + 0.22 * capital_throttle)
        + 0.18 * conviction
        - 0.42 * fragility
        - 0.22 * chase_guard,
        0.08,
        1.82,
    )
    residual_quality = math.exp(-min(50.0, 2.0e5 * max(0.0, residual)))
    shape_quality = _clamp(
        0.30
        + 0.24 * conviction
        + 0.18 * capital_throttle
        + 0.16 * (1.0 - fragility)
        + 0.12 * (1.0 - chase_guard),
        0.0,
        1.0,
    )
    quality = residual_quality * shape_quality if accepted else 0.0

    return {
        "accepted": bool(accepted),
        "status": loc.get("status"),
        "quality": float(_clamp(quality, 0.0, 1.0)),
        "conviction": float(conviction),
        "fragility": float(fragility),
        "allocation": float(allocation),
        "allocation_scale": float(allocation_scale),
        "entry_gate": float(entry_gate),
        "chase_guard": float(chase_guard),
        "cooldown_pressure": float(cooldown_pressure),
        "capital_throttle": float(capital_throttle),
        "chase_pressure": float(chase_pressure),
        "tension_pressure": float(tension_pressure),
        "initial_residual": float(initial_residual),
        "residual": float(residual),
        "epochs": int(loc.get("epochs", 0)),
        "probe_name": loc.get("probe_name"),
        "probe_total_evals": int(loc.get("probe_total_evals", 0)),
    }


def solve_curve_regime_system_118(
    scale_returns_bps: Sequence[float],
    *,
    volatility_bps: float = 5.0,
    volatility_expansion_ratio: float = 1.0,
    acceleration_bps: float = 0.0,
    stretch_bps: float = 0.0,
    swing_compression: float = 0.0,
    lambda_consensus: float = 0.5,
    path_efficiency: float = 0.5,
    drawdown_pressure: float = 0.0,
    direction_seed: int = 0,
) -> dict[str, Any]:
    """Solve a degree-33 118 curve-regime governor.

    This is deliberately curve-only: it consumes multi-timeframe returns and
    path-shape features, not per-coin PnL. Variables are direction, persistence,
    fragility, and harvest quality.
    """
    if len(scale_returns_bps) < 5:
        raise ValueError("scale_returns_bps must contain at least five timeframes")

    engine, template, index = _curve_regime_template_33d()
    exps, weights = template
    np = engine.np
    n = 4
    degree = 33
    zero = (0,) * n

    def unit(axis: int, power: int = 1) -> tuple[int, ...]:
        powers = [0] * n
        powers[axis] = power
        return tuple(powers)

    def mixed(a: int, b: int) -> tuple[int, ...]:
        powers = [0] * n
        powers[a] = 1
        powers[b] = 1
        return tuple(powers)

    def add(row: int, powers: tuple[int, ...], value: complex) -> None:
        coeff[row, index[powers]] += complex(value)

    vol = max(0.1, _finite_float(volatility_bps, 5.0))
    returns = [_clamp(_finite_float(v) / max(1.0, vol), -8.0, 8.0) for v in scale_returns_bps[:9]]
    if len(returns) <= 5:
        weights_by_scale = (0.12, 0.16, 0.24, 0.24, 0.24)
        macro_weight = 0.0
    else:
        raw_weights = (0.16, 0.15, 0.14, 0.13, 0.12, 0.10, 0.08, 0.07, 0.05)[: len(returns)]
        weight_sum = max(1e-9, sum(raw_weights))
        weights_by_scale = tuple(w / weight_sum for w in raw_weights)
        macro_weight = 1.0
    weighted_return = sum(w * r for w, r in zip(weights_by_scale, returns))
    fast_return = 0.58 * returns[0] + 0.42 * returns[1]
    slow_return = 0.46 * returns[2] + 0.30 * returns[3] + 0.24 * returns[4]
    long_tail = returns[5:] if len(returns) > 5 else [slow_return]
    long_return = sum(long_tail) / max(1, len(long_tail))
    sign_pos = sum(1 for r in returns if r > 0.0)
    sign_neg = sum(1 for r in returns if r < 0.0)
    sign_agreement = max(sign_pos, sign_neg) / max(1, len(returns))
    dominant_sign = 1.0 if sign_pos >= sign_neg else -1.0
    sign_conflict = 1.0 - sign_agreement
    bend_pressure = _clamp(abs(fast_return - slow_return) / 3.25, 0.0, 1.0)
    macro_bend_pressure = macro_weight * _clamp(abs(slow_return - long_return) / 3.25, 0.0, 1.0)
    rotation_conflict = (
        1.0
        if fast_return * slow_return < -0.18
        else _clamp(abs(fast_return - slow_return) / 5.0, 0.0, 0.65)
    )
    macro_rotation_conflict = macro_weight * (
        1.0
        if slow_return * long_return < -0.18
        else _clamp(abs(slow_return - long_return) / 5.0, 0.0, 0.65)
    )

    accel_norm = _clamp(_finite_float(acceleration_bps) / vol, -6.0, 6.0)
    stretch_norm = _clamp(_finite_float(stretch_bps) / max(1.0, vol * 5.0), -6.0, 6.0)
    expansion = _clamp(_finite_float(volatility_expansion_ratio, 1.0) - 1.0, -3.0, 3.0)
    compression = _clamp(_finite_float(swing_compression), 0.0, 1.0)
    consensus = _clamp(_finite_float(lambda_consensus, 0.5), 0.0, 1.0)
    efficiency = _clamp(_finite_float(path_efficiency, 0.5), 0.0, 1.0)
    drawdown = _clamp(_finite_float(drawdown_pressure), 0.0, 1.0)
    accel_conflict = 1.0 if accel_norm * dominant_sign < -0.12 else 0.0
    stretch_chase = _clamp(abs(stretch_norm) / 1.8, 0.0, 1.0)
    expansion_pressure = _clamp(max(0.0, expansion) / 1.8, 0.0, 1.0)

    direction_anchor = math.tanh(
        0.50 * weighted_return
        + 0.26 * slow_return
        + 0.16 * fast_return
        + 0.06 * macro_weight * long_return
        + 0.08 * accel_norm
    )
    fragility_anchor = _clamp(
        0.24 * expansion_pressure
        + 0.20 * sign_conflict
        + 0.18 * (1.0 - efficiency)
        + 0.15 * bend_pressure
        + 0.05 * macro_bend_pressure
        + 0.16 * stretch_chase
        + 0.14 * drawdown
        + 0.09 * rotation_conflict
        + 0.04 * macro_rotation_conflict
        + 0.08 * compression
        + 0.12 * accel_conflict
        - 0.14 * consensus,
        0.0,
        1.0,
    )
    persistence_anchor = _clamp(
        0.28 * abs(direction_anchor)
        + 0.24 * sign_agreement
        + 0.18 * consensus
        + 0.16 * efficiency
        + 0.14 * (1.0 - fragility_anchor)
        - 0.08 * rotation_conflict
        - 0.03 * macro_rotation_conflict,
        0.0,
        1.0,
    )
    harvest_anchor = _clamp(
        0.24
        + 0.38 * persistence_anchor
        + 0.18 * abs(direction_anchor)
        + 0.12 * efficiency
        + 0.08 * max(0.0, expansion)
        - 0.30 * fragility_anchor,
        0.0,
        1.0,
    )

    coeff = np.zeros((n, exps.shape[0]), dtype=np.complex128)
    anchors = (direction_anchor, persistence_anchor, fragility_anchor, harvest_anchor)
    for row, anchor in enumerate(anchors):
        add(row, unit(row, degree), 0.010)
        add(row, unit(row), 1.0)
        add(row, zero, -anchor)

    add(0, unit(1), -0.070 * dominant_sign)
    add(0, unit(2), 0.090 * dominant_sign)
    add(0, mixed(0, 2), 0.050)
    add(1, unit(2), 0.170)
    add(1, unit(3), -0.100)
    add(1, mixed(0, 0), -0.040)
    add(2, unit(1), -0.100)
    add(2, unit(3), -0.060)
    add(2, mixed(0, 1), -0.045 * dominant_sign)
    add(2, mixed(0, 2), 0.030 + 0.050 * max(bend_pressure, macro_bend_pressure))
    add(3, unit(2), 0.140)
    add(3, unit(1), -0.110)
    add(3, mixed(0, 0), -0.035)
    add(3, mixed(1, 2), -0.045)

    system = engine.DenseKostlanSystem(
        n=n,
        d=degree,
        seed=0x118C021D ^ int(direction_seed),
        exps=exps,
        coeff=coeff,
        weights=weights,
        equation_normalize=False,
    )
    target = engine.TargetTrack(system, engine.LinearChart.identity(n))

    y0 = np.asarray([complex(_signed_power_root(anchor, degree), 0.0) for anchor in anchors], dtype=np.complex128)
    initial_residual = target.residual(y0)

    loc = _swift129_corrector(
        engine,
        n=n,
        degree=degree,
        exps=exps,
        coeff=coeff,
        y0=y0,
        max_epochs=4,
        tol=1e-10,
        accept=4e-7,
        line_search=6,
        probe_scale=0.026,
        direction_seed=int(direction_seed),
        probe_candidates=4,
        probe_radii=(0.0, 0.12, 0.30, 0.70),
        include_self_probe=True,
    )
    if loc is None:
        loc = engine.pandrosion_corrector(
            target,
            y0,
            max_epochs=4,
            tol=1e-10,
            accept=4e-7,
            trial_timeout=0.0,
            line_search=6,
            probe_scale=0.026,
            direction_seed=int(direction_seed),
            probe_candidates=4,
            probe_radii=(0.0, 0.12, 0.30, 0.70),
            include_self_probe=True,
        )

    y = np.asarray(loc.get("y", y0), dtype=np.complex128)
    residual = _finite_float(loc.get("residual"), float("inf"))
    accepted = bool(loc.get("accepted", False)) and math.isfinite(residual) and residual < 8e-5

    direction_score = _clamp(float(y[0].real), -1.0, 1.0)
    persistence = _clamp(float(y[1].real), 0.0, 1.0)
    fragility = _clamp(float(y[2].real), 0.0, 1.0)
    harvest = _clamp(float(y[3].real), 0.0, 1.0)
    direction = "BUY" if direction_score >= 0.0 else "SELL"
    residual_quality = math.exp(-min(50.0, 1.2e5 * max(0.0, residual)))
    shape_quality = _clamp(
        0.18
        + 0.24 * abs(direction_score)
        + 0.22 * persistence
        + 0.18 * harvest
        + 0.18 * (1.0 - fragility),
        0.0,
        1.0,
    )
    quality = residual_quality * shape_quality if accepted else 0.0
    reversal_pressure = _clamp(
        0.42 * fragility
        + 0.22 * sign_conflict
        + 0.18 * accel_conflict
        + 0.18 * drawdown
        + 0.20 * rotation_conflict
        + 0.06 * macro_rotation_conflict
        + 0.14 * bend_pressure
        + 0.04 * macro_bend_pressure
        - 0.22 * persistence,
        0.0,
        1.0,
    )
    geometry_cleanliness = _clamp(
        0.26 * persistence
        + 0.24 * harvest
        + 0.20 * efficiency
        + 0.16 * sign_agreement
        + 0.14 * consensus
        - 0.34 * fragility
        - 0.24 * reversal_pressure
        - 0.10 * bend_pressure
        - 0.03 * macro_bend_pressure,
        0.0,
        1.0,
    )
    exposure_scale = _clamp(
        0.18
        + 1.34 * quality * geometry_cleanliness
        + 0.24 * harvest
        + 0.22 * persistence
        - 1.08 * fragility
        - 0.58 * reversal_pressure
        - 0.18 * bend_pressure,
        0.05,
        1.92,
    )

    return {
        "accepted": bool(accepted),
        "status": loc.get("status"),
        "degree": degree,
        "quality": float(_clamp(quality, 0.0, 1.0)),
        "direction": direction,
        "direction_score": float(direction_score),
        "persistence": float(persistence),
        "fragility": float(fragility),
        "harvest_score": float(harvest),
        "reversal_pressure": float(reversal_pressure),
        "exposure_scale": float(exposure_scale),
        "geometry_cleanliness": float(geometry_cleanliness),
        "bend_pressure": float(bend_pressure),
        "rotation_conflict": float(rotation_conflict),
        "macro_bend_pressure": float(macro_bend_pressure),
        "macro_rotation_conflict": float(macro_rotation_conflict),
        "long_return_norm": float(long_return),
        "scale_count": int(len(returns)),
        "sign_agreement": float(sign_agreement),
        "path_efficiency": float(efficiency),
        "drawdown_pressure": float(drawdown),
        "initial_residual": float(initial_residual),
        "residual": float(residual),
        "epochs": int(loc.get("epochs", 0)),
        "probe_name": loc.get("probe_name"),
        "probe_total_evals": int(loc.get("probe_total_evals", 0)),
    }


def solve_market_consensus_probe(
    scale_scores: Sequence[float],
    liquidity_bias: float = 0.0,
    direction_seed: int = 0,
) -> dict[str, Any]:
    """Project short/mid/long market scores into a probe-aware trend basin.

    `scale_scores` should contain three standardized signed scores ordered from
    fast to slow.  The 118 corrector then decides whether the point lies in the
    coherent positive or negative basin and reports the residual of that pure
    derivative-free solve.
    """
    if len(scale_scores) != 3:
        raise ValueError("scale_scores must contain exactly three values")

    engine = _load_engine_118()
    np = engine.np
    target = _market_consensus_target()

    imag = max(-0.12, min(0.12, _finite_float(liquidity_bias) * 0.05))
    start = []
    for raw in scale_scores:
        score = max(-4.0, min(4.0, _finite_float(raw)))
        start.append(complex(math.tanh(score), imag))
    y0 = np.asarray(start, dtype=np.complex128)
    initial_residual = target.residual(y0)

    target_system = getattr(target, "system", None)
    loc = None
    if target_system is not None:
        loc = _swift129_corrector(
            engine,
            n=3,
            degree=2,
            exps=getattr(target_system, "exps", None),
            coeff=getattr(target_system, "coeff", None),
            y0=y0,
            max_epochs=4,
            tol=1e-10,
            accept=1e-8,
            line_search=8,
            probe_scale=0.08,
            direction_seed=int(direction_seed),
            probe_candidates=6,
            probe_radii=(0.0, 0.35, 0.7, 1.0, 1.6),
            include_self_probe=True,
        )
    if loc is None:
        loc = engine.pandrosion_corrector(
            target,
            y0,
            max_epochs=4,
            tol=1e-10,
            accept=1e-8,
            trial_timeout=0.0,
            line_search=8,
            probe_scale=0.08,
            direction_seed=int(direction_seed),
            probe_candidates=6,
            probe_radii=(0.0, 0.35, 0.7, 1.0, 1.6),
            include_self_probe=True,
        )

    y = np.asarray(loc.get("y", y0), dtype=np.complex128)
    residual = _finite_float(loc.get("residual"), float("inf"))
    mean_real = float(np.mean(y.real))
    displacement = float(np.linalg.norm(y - y0) / math.sqrt(3.0))
    direction = "BUY" if mean_real >= 0.0 else "SELL"
    residual_quality = math.exp(-min(50.0, 25.0 * max(0.0, residual)))
    displacement_quality = math.exp(-0.35 * min(5.0, max(0.0, displacement)))
    accepted = bool(loc.get("accepted", False)) and residual < 1e-6
    quality = residual_quality * displacement_quality if accepted else 0.0

    return {
        "accepted": accepted,
        "status": loc.get("status"),
        "direction": direction,
        "direction_sign": 1.0 if direction == "BUY" else -1.0,
        "quality": float(max(0.0, min(1.0, quality))),
        "initial_residual": float(initial_residual),
        "residual": float(residual),
        "displacement": float(displacement),
        "epochs": int(loc.get("epochs", 0)),
        "probe_name": loc.get("probe_name"),
        "probe_residual": _finite_float(loc.get("probe_residual"), float("inf")),
        "probe_total_evals": int(loc.get("probe_total_evals", 0)),
    }
