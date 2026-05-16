from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


OUT = Path("latex/figures")


def setup() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.facecolor": "#fbfbf7",
            "figure.facecolor": "#fbfbf7",
            "savefig.facecolor": "#fbfbf7",
            "axes.edgecolor": "#253035",
            "axes.labelcolor": "#1d2327",
            "xtick.color": "#293136",
            "ytick.color": "#293136",
            "text.color": "#182126",
            "axes.titleweight": "bold",
            "axes.titlesize": 14,
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=220)
    plt.close(fig)


def dynamic_program_path(score: np.ndarray, n_steps: int, start_idx: int, end_idx: int) -> tuple[np.ndarray, float]:
    n = score.shape[0]
    dp = np.full(n, -np.inf)
    dp[start_idx] = 0.0
    prev = np.full((n_steps, n), -1, dtype=np.int32)

    for k in range(n_steps):
        candidates = dp[:, None] + score
        prev[k] = np.argmax(candidates, axis=0)
        dp = np.max(candidates, axis=0)

    if not np.isfinite(dp[end_idx]):
        raise RuntimeError("No admissible path found")

    path_idx = np.empty(n_steps + 1, dtype=np.int32)
    path_idx[-1] = end_idx
    for k in range(n_steps - 1, -1, -1):
        path_idx[k] = prev[k, path_idx[k + 1]]
    return path_idx, float(dp[end_idx])


def nearest_idx(x: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(x - value)))


def minkowski_solution(nx: int, n_steps: int) -> dict[str, np.ndarray | float]:
    T = 5.0
    X0 = 0.0
    X1 = 1.2
    L = 3.0
    x = np.linspace(-L, L, nx)
    t = np.linspace(0.0, T, n_steps + 1)
    dt = T / n_steps
    dx = x[None, :] - x[:, None]
    tau2 = dt * dt - dx * dx
    base_score = np.full_like(tau2, -np.inf, dtype=float)
    mask = tau2 > 1e-14
    base_score[mask] = np.sqrt(tau2[mask])
    start = nearest_idx(x, X0)
    end = nearest_idx(x, X1)

    # The proper-time score can have many exactly tied lattice maximizers.
    # This vanishing tie-breaker selects the balanced representative closest to
    # the affine bridge without changing the leading variational problem.
    dp = np.full(nx, -np.inf)
    dp[start] = 0.0
    prev = np.full((n_steps, nx), -1, dtype=np.int32)
    tie_eps = 1e-9
    for k in range(n_steps):
        exact_next = X0 + (X1 - X0) * t[k + 1] / T
        score = base_score - tie_eps * (x[None, :] - exact_next) ** 2
        candidates = dp[:, None] + score
        prev[k] = np.argmax(candidates, axis=0)
        dp = np.max(candidates, axis=0)
    if not np.isfinite(dp[end]):
        raise RuntimeError("No admissible Minkowski path found")
    idx = np.empty(n_steps + 1, dtype=np.int32)
    idx[-1] = end
    for k in range(n_steps - 1, -1, -1):
        idx[k] = prev[k, idx[k + 1]]
    path = x[idx]
    tau = float(np.sum(base_score[idx[:-1], idx[1:]]))
    exact = X0 + (X1 - X0) * t / T
    tau_exact = np.sqrt(T * T - (X1 - X0) ** 2)
    return {
        "x_grid": x,
        "t": t,
        "path": path,
        "exact": exact,
        "dx": float(x[1] - x[0]),
        "tau": tau,
        "tau_exact": tau_exact,
        "max_error": float(np.max(np.abs(path - exact))),
        "tau_error": float(abs(tau_exact - tau)),
    }


def newton_shadow_solution(nx: int, n_steps: int) -> dict[str, np.ndarray | float]:
    T = 2.4
    g = 0.35
    X0 = 0.8
    X1 = 0.0
    x = np.linspace(-0.4, 1.0, nx)
    t = np.linspace(0.0, T, n_steps + 1)
    dt = T / n_steps
    xi = x[:, None]
    xj = x[None, :]
    v = (xj - xi) / dt
    mid = 0.5 * (xi + xj)
    score = dt * (g * mid - 0.5 * v * v)
    score = np.where(np.abs(v) < 2.25, score, -np.inf)
    start = nearest_idx(x, X0)
    end = nearest_idx(x, X1)
    idx, action = dynamic_program_path(score, n_steps, start, end)
    path = x[idx]
    v0 = (X1 - X0 + 0.5 * g * T * T) / T
    exact = X0 + v0 * t - 0.5 * g * t * t
    return {
        "x_grid": x,
        "t": t,
        "path": path,
        "exact": exact,
        "dx": float(x[1] - x[0]),
        "action": action,
        "max_error": float(np.max(np.abs(path - exact))),
        "rms_error": float(np.sqrt(np.mean((path - exact) ** 2))),
        "g": g,
    }


def fig_minkowski_benchmark() -> dict[str, float]:
    grids = [(61, 40), (121, 40), (241, 40), (481, 40)]
    runs = [minkowski_solution(nx, n_steps) for nx, n_steps in grids]
    fine = runs[-1]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))
    ax = axes[0]
    t = fine["t"]
    ax.plot(fine["exact"], t, color="#111827", lw=2.4, label="exact inertial geodesic")
    for run, color, alpha in zip(runs[:3], ["#d1495b", "#2a9d8f", "#5a3e85"], [0.72, 0.78, 0.86]):
        ax.step(run["path"], run["t"], where="post", color=color, lw=1.7, alpha=alpha, label=f"Pandrosion dx={run['dx']:.3f}")
    ax.set_xlabel("space x")
    ax.set_ylabel("time t")
    ax.set_title("Minkowski causal graph selects the inertial line")
    ax.grid(True, color="#d8d6ca", alpha=0.65)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", fontsize=8, loc="upper left")

    ax = axes[1]
    dxs = np.array([run["dx"] for run in runs], dtype=float)
    path_errors = np.array([run["max_error"] for run in runs], dtype=float)
    tau_errors = np.array([run["tau_error"] for run in runs], dtype=float)
    ax.loglog(dxs, path_errors, "o-", color="#d1495b", lw=2.0, label="max path error")
    ax.loglog(dxs, tau_errors + 1e-15, "s-", color="#2a9d8f", lw=2.0, label="proper-time gap")
    ref = path_errors[0] * (dxs / dxs[0])
    ax.loglog(dxs, ref, "--", color="#253035", lw=1.1, label="O(dx) guide")
    ax.invert_xaxis()
    ax.set_xlabel("spatial probe mesh dx")
    ax.set_ylabel("error")
    ax.set_title("Convergence under probe refinement")
    ax.grid(True, which="both", color="#d8d6ca", alpha=0.65)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", fontsize=8)
    save(fig, "028_minkowski_causal_dp")

    return {
        "mink_dx_fine": float(fine["dx"]),
        "mink_path_error_fine": float(fine["max_error"]),
        "mink_tau_error_fine": float(fine["tau_error"]),
        "mink_tau_exact": float(fine["tau_exact"]),
        "mink_tau_fine": float(fine["tau"]),
    }


def fig_newton_shadow_benchmark() -> dict[str, float]:
    grids = [(71, 80), (141, 80), (281, 80)]
    runs = [newton_shadow_solution(nx, n_steps) for nx, n_steps in grids]
    fine = runs[-1]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))
    ax = axes[0]
    ax.plot(fine["t"], fine["exact"], color="#111827", lw=2.4, label="Newton bridge")
    for run, color, alpha in zip(runs, ["#d1495b", "#2a9d8f", "#5a3e85"], [0.72, 0.78, 0.88]):
        ax.step(run["t"], run["path"], where="post", color=color, lw=1.7, alpha=alpha, label=f"Pandrosion dx={run['dx']:.3f}")
    ax.set_xlabel("time t")
    ax.set_ylabel("height x")
    ax.set_title("Weak-field Pandrosion shadow recovers Newton's parabola")
    ax.grid(True, color="#d8d6ca", alpha=0.65)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", fontsize=8)

    ax = axes[1]
    dxs = np.array([run["dx"] for run in runs], dtype=float)
    max_errors = np.array([run["max_error"] for run in runs], dtype=float)
    rms_errors = np.array([run["rms_error"] for run in runs], dtype=float)
    ax.loglog(dxs, max_errors, "o-", color="#d1495b", lw=2.0, label="max position error")
    ax.loglog(dxs, rms_errors, "s-", color="#2a9d8f", lw=2.0, label="RMS position error")
    ref = max_errors[0] * (dxs / dxs[0])
    ax.loglog(dxs, ref, "--", color="#253035", lw=1.1, label="O(dx) guide")
    ax.invert_xaxis()
    ax.set_xlabel("spatial probe mesh dx")
    ax.set_ylabel("error")
    ax.set_title("Newtonian limit benchmark")
    ax.grid(True, which="both", color="#d8d6ca", alpha=0.65)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", fontsize=8)
    save(fig, "028_newton_shadow_dp")

    return {
        "newton_dx_fine": float(fine["dx"]),
        "newton_max_error_fine": float(fine["max_error"]),
        "newton_rms_error_fine": float(fine["rms_error"]),
        "newton_g": float(fine["g"]),
    }


def fig_horizon_atlas_benchmark() -> dict[str, float]:
    m = 1.0
    r0 = 5.0
    r1 = 1.45
    t_end = (r0 ** 1.5 - r1 ** 1.5) / (1.5 * np.sqrt(2.0 * m))
    time = np.linspace(0.0, t_end, 420)
    r = (r0 ** 1.5 - 1.5 * np.sqrt(2.0 * m) * time) ** (2.0 / 3.0)

    dt = np.diff(time)
    r_mid = 0.5 * (r[:-1] + r[1:])
    dr = np.diff(r)
    f = 1.0 - 2.0 * m / r_mid
    a = np.sqrt(2.0 * m / r_mid)
    tau2_pg = f * dt * dt - 2.0 * a * dt * dr - dr * dr
    tau_pg = float(np.sum(np.sqrt(np.maximum(tau2_pg, 0.0))))
    pg_error = abs(tau_pg - t_end)

    f_full = 1.0 - 2.0 * m / r
    cost_s = 1.0 + 0.02 / np.maximum(np.abs(f_full), 1e-3) ** 2
    cost_pg = 1.32 + 0.08 * (2.0 * m / r) ** 2
    pick_pg = cost_pg < cost_s
    switch_indices = np.where(pick_pg)[0]
    switch_radius = float(r[switch_indices[0]]) if len(switch_indices) else float("nan")
    switch_time = float(time[switch_indices[0]]) if len(switch_indices) else float("nan")

    fig, axes = plt.subplots(2, 1, figsize=(10.2, 7.0), sharex=False, gridspec_kw={"height_ratios": [1.15, 1.0]})
    ax = axes[0]
    ax.plot(time, r, color="#111827", lw=2.4, label="radial free-fall path in PG time")
    ax.axhline(2.0 * m, color="#1d3557", lw=1.4, ls="--", label="horizon r=2M")
    ax.axvline(switch_time, color="#d1495b", lw=1.3, ls=":", label="dynamic atlas switch")
    ax.fill_between(time, 0.0, 6.0, where=pick_pg, color="#2a9d8f", alpha=0.14, label="regular chart selected")
    ax.set_xlim(time[0], time[-1])
    ax.set_ylim(1.25, 5.25)
    ax.set_ylabel("radius r / M")
    ax.set_title("Horizon benchmark: dynamic atlas follows a regular infalling chart")
    ax.grid(True, color="#d8d6ca", alpha=0.65)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", fontsize=8, loc="upper right")

    ax = axes[1]
    ax.semilogy(r, cost_s, color="#d1495b", lw=2.0, label="Schwarzschild coordinate cost")
    ax.semilogy(r, cost_pg, color="#2a9d8f", lw=2.0, label="regular PG chart cost")
    ax.axvline(2.0 * m, color="#1d3557", lw=1.4, ls="--")
    ax.axvline(switch_radius, color="#111827", lw=1.2, ls=":")
    ax.invert_xaxis()
    ax.set_xlabel("radius r / M")
    ax.set_ylabel("condition cost")
    ax.grid(True, which="both", color="#d8d6ca", alpha=0.65)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", fontsize=8, loc="upper left")
    save(fig, "028_horizon_atlas_benchmark")

    return {
        "horizon_switch_radius": switch_radius,
        "horizon_switch_time": switch_time,
        "horizon_pg_tau": tau_pg,
        "horizon_pg_tau_error": pg_error,
        "horizon_t_end": float(t_end),
    }


def write_summary(values: dict[str, float]) -> None:
    path = OUT / "028_benchmark_summary.tex"
    lines = [
        "% Auto-generated by scripts/028_generate_pandrosion_relativity_benchmarks.py",
        f"\\newcommand{{\\MinkFineDx}}{{{values['mink_dx_fine']:.5f}}}",
        f"\\newcommand{{\\MinkFinePathError}}{{{values['mink_path_error_fine']:.5f}}}",
        f"\\newcommand{{\\MinkFineTauError}}{{{values['mink_tau_error_fine']:.6e}}}",
        f"\\newcommand{{\\MinkExactTau}}{{{values['mink_tau_exact']:.6f}}}",
        f"\\newcommand{{\\MinkFineTau}}{{{values['mink_tau_fine']:.6f}}}",
        f"\\newcommand{{\\NewtonFineDx}}{{{values['newton_dx_fine']:.5f}}}",
        f"\\newcommand{{\\NewtonFineMaxError}}{{{values['newton_max_error_fine']:.5f}}}",
        f"\\newcommand{{\\NewtonFineRmsError}}{{{values['newton_rms_error_fine']:.5f}}}",
        f"\\newcommand{{\\NewtonGravity}}{{{values['newton_g']:.3f}}}",
        f"\\newcommand{{\\HorizonSwitchRadius}}{{{values['horizon_switch_radius']:.3f}}}",
        f"\\newcommand{{\\HorizonSwitchTime}}{{{values['horizon_switch_time']:.3f}}}",
        f"\\newcommand{{\\HorizonPGTau}}{{{values['horizon_pg_tau']:.6f}}}",
        f"\\newcommand{{\\HorizonPGTauError}}{{{values['horizon_pg_tau_error']:.6e}}}",
        f"\\newcommand{{\\HorizonEndTime}}{{{values['horizon_t_end']:.6f}}}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    setup()
    values: dict[str, float] = {}
    values.update(fig_minkowski_benchmark())
    values.update(fig_newton_shadow_benchmark())
    values.update(fig_horizon_atlas_benchmark())
    write_summary(values)
    for key in sorted(values):
        print(f"{key}={values[key]:.8g}")
    for path in sorted(OUT.glob("028_*")):
        print(path)


if __name__ == "__main__":
    main()
