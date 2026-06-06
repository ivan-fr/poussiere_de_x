from __future__ import annotations

import csv
import functools
import json
import math
import random
import statistics
import time
from pathlib import Path
from typing import Callable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


OUT_BENCH = Path("benchmarks/029_fused_irp")
OUT_FIG = Path("latex/figures")

P_VALUES = [2, 3, 5, 8, 12, 20]
TOL = 2.0e-13
MAX_ITERS = 14
SAMPLES = 7
MIN_TIME = 0.18
VECTOR_N = 200_000
LOG2 = math.log(2.0)


def setup() -> None:
    OUT_BENCH.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.facecolor": "#fbfbf7",
            "figure.facecolor": "#fbfbf7",
            "savefig.facecolor": "#fbfbf7",
            "axes.edgecolor": "#263238",
            "axes.labelcolor": "#1f2933",
            "xtick.color": "#263238",
            "ytick.color": "#263238",
            "text.color": "#172026",
            "axes.titleweight": "bold",
            "axes.titlesize": 13,
        }
    )


def binary_normalize(x: float, p: int) -> tuple[int, float]:
    """Return m,y with x=2^(pm)y and 1 <= y < 2^p."""
    _, exp2 = math.frexp(x)
    m = exp2 // p
    y = math.ldexp(x, -p * m)
    if y < 1.0:
        m -= 1
        y = math.ldexp(y, p)
    elif y >= 2.0**p:
        m += 1
        y = math.ldexp(y, -p)
    return m, y


def newton_scaled_fixed(x: float, p: int, steps: int) -> float:
    m, y = binary_normalize(x, p)
    u = 1.0
    invp = 1.0 / p
    pm1 = p - 1
    for _ in range(steps):
        u = (pm1 * u + y / (u ** pm1)) * invp
    return math.ldexp(u, m)


def palette_k(p: int) -> int:
    return max(4, p * p)


def nearest_palette_shift(eta: float, p: int) -> tuple[float, float]:
    root_step = LOG2 / palette_k(p)
    residual_period = p * root_step
    j = int(round(eta / residual_period))
    root_shift = j * root_step
    return root_shift, eta - p * root_shift


def newton_palette_fixed(x: float, p: int, steps: int) -> float:
    m, y = binary_normalize(x, p)
    eta = math.log(y)
    root_shift, local_eta = nearest_palette_shift(eta, p)
    local_y = math.exp(local_eta)
    u = 1.0
    invp = 1.0 / p
    pm1 = p - 1
    for _ in range(steps):
        u = (pm1 * u + local_y / (u ** pm1)) * invp
    return math.ldexp(math.exp(root_shift) * u, m)


def theta_exact(p: int, d: float) -> float:
    if d == 0.0:
        return 0.0
    return d + (p - 1) * math.log1p(math.expm1(-d) / p)


@functools.lru_cache(maxsize=None)
def theta_poly6_coeffs(p: int) -> tuple[float, float, float, float, float, float]:
    pp = float(p)
    c1 = 1.0 / pp
    c2 = (pp - 1.0) ** 2 / (2.0 * pp**2)
    c3 = -(pp - 1.0) ** 2 * (pp - 2.0) / (6.0 * pp**3)
    c4 = (pp**4 - 8.0 * pp**3 + 19.0 * pp**2 - 18.0 * pp + 6.0) / (24.0 * pp**4)
    c5 = (-pp**5 + 16.0 * pp**4 - 65.0 * pp**3 + 110.0 * pp**2 - 84.0 * pp + 24.0) / (120.0 * pp**5)
    c6 = (
        pp**6
        - 32.0 * pp**5
        + 211.0 * pp**4
        - 570.0 * pp**3
        + 750.0 * pp**2
        - 480.0 * pp
        + 120.0
    ) / (720.0 * pp**6)
    return c1, c2, c3, c4, c5, c6


def theta_poly6(p: int, d: float) -> float:
    c1, c2, c3, c4, c5, c6 = theta_poly6_coeffs(p)
    return d * (c1 + d * (c2 + d * (c3 + d * (c4 + d * (c5 + d * c6)))))


def theta_hybrid(p: int, d: float) -> float:
    if d < 0.05:
        return theta_poly6(p, d)
    return theta_exact(p, d)


def fused_irp_fixed(x: float, p: int, steps: int, theta: Callable[[int, float], float]) -> float:
    m, y = binary_normalize(x, p)
    eta = math.log(y)
    corr = 0.0
    for _ in range(steps):
        root_shift, local_eta = nearest_palette_shift(eta, p)
        corr += root_shift
        d = abs(local_eta)
        if d == 0.0:
            eta = 0.0
            break
        sgn = 1.0 if local_eta > 0.0 else -1.0
        th = theta(p, d)
        corr += sgn * th
        eta = local_eta - p * sgn * th
    return math.ldexp(math.exp(corr), m)


def fused_irp_exact_fixed(x: float, p: int, steps: int) -> float:
    return fused_irp_fixed(x, p, steps, theta_exact)


def fused_irp_hybrid_fixed(x: float, p: int, steps: int) -> float:
    return fused_irp_fixed(x, p, steps, theta_hybrid)


def fused_irp_quad_fixed(x: float, p: int, steps: int) -> float:
    m, y = binary_normalize(x, p)
    eta = math.log(y)
    corr = 0.0
    pp = float(p)
    c1 = 1.0 / pp
    c2 = (pp - 1.0) ** 2 / (2.0 * pp**2)
    root_step = LOG2 / palette_k(p)
    residual_period = p * root_step
    for _ in range(steps):
        j = int(round(eta / residual_period))
        root_shift = j * root_step
        local_eta = eta - p * root_shift
        corr += root_shift
        d = abs(local_eta)
        if d == 0.0:
            eta = 0.0
            break
        sgn = 1.0 if local_eta > 0.0 else -1.0
        th = d * (c1 + c2 * d)
        corr += sgn * th
        eta = local_eta - p * sgn * th
    return math.ldexp(math.exp(corr), m)


def math_pow_root(x: float, p: int, steps: int = 0) -> float:
    _ = steps
    return math.pow(x, 1.0 / p)


def make_cases() -> list[float]:
    rng = random.Random(29029)
    exponents = [-900, -700, -500, -300, -120, -40, -12, -3, 0, 3, 12, 40, 120, 300, 500, 700, 900]
    cases: list[float] = []
    for exp2 in exponents:
        for _ in range(6):
            mant = 1.0 + rng.random()
            cases.append(math.ldexp(mant, exp2))
    rng.shuffle(cases)
    return cases


def rel_error(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0e-300)


def max_error(method: Callable[[float, int, int], float], cases: list[float], p: int, steps: int) -> float:
    worst = 0.0
    for x in cases:
        ref = math_pow_root(x, p)
        worst = max(worst, rel_error(method(x, p, steps), ref))
    return worst


def calibrate(method: Callable[[float, int, int], float], cases: list[float], p: int) -> tuple[int, float]:
    for steps in range(MAX_ITERS + 1):
        err = max_error(method, cases, p, steps)
        if err <= TOL:
            return steps, err
    return MAX_ITERS, max_error(method, cases, p, MAX_ITERS)


def bench_method(
    method: Callable[[float, int, int], float],
    cases: list[float],
    p: int,
    steps: int,
) -> tuple[float, float]:
    checksum = 0.0
    for x in cases:
        checksum += method(x, p, steps)

    loops = 1
    while True:
        start = time.perf_counter_ns()
        acc = 0.0
        for _ in range(loops):
            for x in cases:
                acc += method(x, p, steps)
        elapsed = (time.perf_counter_ns() - start) / 1.0e9
        checksum += acc
        if elapsed >= MIN_TIME:
            break
        loops *= 2

    samples: list[float] = []
    for _ in range(SAMPLES):
        start = time.perf_counter_ns()
        acc = 0.0
        for _ in range(loops):
            for x in cases:
                acc += method(x, p, steps)
        elapsed_ns = time.perf_counter_ns() - start
        checksum += acc
        samples.append(elapsed_ns / (loops * len(cases)))
    return statistics.median(samples), checksum


def write_table(rows: list[dict[str, object]]) -> None:
    table_path = OUT_FIG / "029_fused_irp_benchmark_summary.tex"
    by_p = {int(row["p"]): row for row in rows}
    with table_path.open("w", encoding="utf-8") as f:
        f.write("\\begin{tabular}{rrrrrrrrr}\n")
        f.write("\\toprule\n")
        f.write(
            "$p$ & Newton-K it. & Newton-K ns & IRP exact it. & IRP exact ns & "
            "IRP quad it. & IRP quad ns & quad/Newton-K & pow ns \\\\\n"
        )
        f.write("\\midrule\n")
        for p in P_VALUES:
            row = by_p[p]
            newton = row["methods"]["Newton-K"]  # type: ignore[index]
            exact = row["methods"]["IRP-F exact"]  # type: ignore[index]
            quad = row["methods"]["IRP-F quad"]  # type: ignore[index]
            pow_row = row["methods"]["math.pow"]  # type: ignore[index]
            ratio = float(quad["ns_per_root"]) / float(newton["ns_per_root"])
            f.write(
                f"{p} & {newton['steps']} & {float(newton['ns_per_root']):.1f} & "
                f"{exact['steps']} & {float(exact['ns_per_root']):.1f} & "
                f"{quad['steps']} & {float(quad['ns_per_root']):.1f} & "
                f"{ratio:.2f} & {float(pow_row['ns_per_root']):.1f} \\\\\n"
            )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")


def write_plot(rows: list[dict[str, object]]) -> None:
    ps = [int(row["p"]) for row in rows]
    speedups = []
    exact_speedups = []
    for row in rows:
        methods = row["methods"]  # type: ignore[assignment]
        n = float(methods["Newton-K"]["ns_per_root"])  # type: ignore[index]
        speedups.append(n / float(methods["IRP-F quad"]["ns_per_root"]))  # type: ignore[index]
        exact_speedups.append(n / float(methods["IRP-F exact"]["ns_per_root"]))  # type: ignore[index]

    fig, ax = plt.subplots(figsize=(7.2, 3.7))
    width = 0.34
    xs = list(range(len(ps)))
    ax.bar([x - width / 2 for x in xs], exact_speedups, width, label="IRP-F exact", color="#3a6ea5")
    ax.bar([x + width / 2 for x in xs], speedups, width, label="IRP-F quad", color="#8f5a83")
    ax.axhline(1.0, color="#222222", linewidth=1.0)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(p) for p in ps])
    ax.set_xlabel("monomial degree p")
    ax.set_ylabel("speedup vs Newton-K")
    ax.set_title("Fused IRP wall-clock ratio on calibrated scalar roots")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#d9d4c7", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "029_fused_irp_speedup_vs_newton.pdf", bbox_inches="tight")
    fig.savefig(OUT_FIG / "029_fused_irp_speedup_vs_newton.png", bbox_inches="tight", dpi=220)
    plt.close(fig)


def make_vector_cases() -> np.ndarray:
    rng = np.random.default_rng(29029)
    exponents = rng.integers(-900, 901, size=VECTOR_N, dtype=np.int32)
    mantissas = 1.0 + rng.random(VECTOR_N)
    return np.ldexp(mantissas, exponents)


def vector_binary_normalize(x: np.ndarray, p: int) -> tuple[np.ndarray, np.ndarray]:
    _, exp2 = np.frexp(x)
    m = np.floor_divide(exp2, p).astype(np.int32)
    y = np.ldexp(x, -p * m)
    mask = y < 1.0
    if np.any(mask):
        m = np.where(mask, m - 1, m).astype(np.int32)
        y = np.where(mask, np.ldexp(y, p), y)
    return m, y


def vector_newton_palette(x: np.ndarray, p: int, steps: int) -> np.ndarray:
    m, y = vector_binary_normalize(x, p)
    eta = np.log(y)
    root_step = LOG2 / palette_k(p)
    residual_period = p * root_step
    j = np.rint(eta / residual_period)
    root_shift = j * root_step
    local_eta = eta - p * root_shift
    local_y = np.exp(local_eta)
    u = np.ones_like(x)
    invp = 1.0 / p
    pm1 = p - 1
    for _ in range(steps):
        u = (pm1 * u + local_y / (u**pm1)) * invp
    return np.ldexp(np.exp(root_shift) * u, m)


def vector_irp_quad(x: np.ndarray, p: int, steps: int) -> np.ndarray:
    m, y = vector_binary_normalize(x, p)
    eta = np.log(y)
    corr = np.zeros_like(x)
    pp = float(p)
    c1 = 1.0 / pp
    c2 = (pp - 1.0) ** 2 / (2.0 * pp**2)
    root_step = LOG2 / palette_k(p)
    residual_period = p * root_step
    for _ in range(steps):
        j = np.rint(eta / residual_period)
        root_shift = j * root_step
        local_eta = eta - p * root_shift
        corr += root_shift
        d = np.abs(local_eta)
        sgn = np.where(local_eta >= 0.0, 1.0, -1.0)
        th = d * (c1 + c2 * d)
        corr += sgn * th
        eta = local_eta - p * sgn * th
    return np.ldexp(np.exp(corr), m)


def vector_max_error(values: np.ndarray, refs: np.ndarray) -> float:
    return float(np.max(np.abs(values - refs) / np.maximum(np.abs(refs), 1.0e-300)))


def bench_vector(fn: Callable[[np.ndarray, int, int], np.ndarray], x: np.ndarray, p: int, steps: int) -> tuple[float, float]:
    checksum = float(np.sum(fn(x, p, steps)[:16]))
    samples: list[float] = []
    for _ in range(SAMPLES):
        start = time.perf_counter_ns()
        values = fn(x, p, steps)
        elapsed = time.perf_counter_ns() - start
        checksum += float(np.sum(values[:16]))
        samples.append(elapsed / len(x))
    return statistics.median(samples), checksum


def run_vector_benchmark(scalar_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    x = make_vector_cases()
    by_p = {int(row["p"]): row for row in scalar_rows}
    rows: list[dict[str, object]] = []
    for p in P_VALUES:
        refs = np.power(x, 1.0 / p)
        newton_steps = int(by_p[p]["methods"]["Newton-K"]["steps"])  # type: ignore[index]
        irp_steps = int(by_p[p]["methods"]["IRP-F quad"]["steps"])  # type: ignore[index]
        newton_values = vector_newton_palette(x, p, newton_steps)
        irp_values = vector_irp_quad(x, p, irp_steps)
        newton_err = vector_max_error(newton_values, refs)
        irp_err = vector_max_error(irp_values, refs)
        newton_ns, newton_checksum = bench_vector(vector_newton_palette, x, p, newton_steps)
        irp_ns, irp_checksum = bench_vector(vector_irp_quad, x, p, irp_steps)
        speedup = newton_ns / irp_ns
        rows.append(
            {
                "p": p,
                "case_count": len(x),
                "Newton-K": {
                    "steps": newton_steps,
                    "max_rel_error": newton_err,
                    "ns_per_root": newton_ns,
                    "checksum": newton_checksum,
                },
                "IRP-F quad": {
                    "steps": irp_steps,
                    "max_rel_error": irp_err,
                    "ns_per_root": irp_ns,
                    "checksum": irp_checksum,
                },
                "speedup_vs_newton": speedup,
            }
        )
        print(
            f"vector p={p:2d} Newton-K={newton_ns:6.2f} ns/root "
            f"IRP-F quad={irp_ns:6.2f} ns/root speedup={speedup:.2f}"
        )
    return rows


def write_vector_table(rows: list[dict[str, object]]) -> None:
    table_path = OUT_FIG / "029_fused_irp_vector_benchmark_summary.tex"
    with table_path.open("w", encoding="utf-8") as f:
        f.write("\\begin{tabular}{rrrrrrr}\n")
        f.write("\\toprule\n")
        f.write("$p$ & Newton-K it. & Newton-K ns & IRP quad it. & IRP quad ns & speedup & IRP max err \\\\\n")
        f.write("\\midrule\n")
        for row in rows:
            newton = row["Newton-K"]  # type: ignore[index]
            irp = row["IRP-F quad"]  # type: ignore[index]
            f.write(
                f"{row['p']} & {newton['steps']} & {float(newton['ns_per_root']):.2f} & "
                f"{irp['steps']} & {float(irp['ns_per_root']):.2f} & "
                f"{float(row['speedup_vs_newton']):.2f} & {float(irp['max_rel_error']):.1e} \\\\\n"
            )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")


def write_vector_plot(rows: list[dict[str, object]]) -> None:
    ps = [int(row["p"]) for row in rows]
    speedups = [float(row["speedup_vs_newton"]) for row in rows]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.bar([str(p) for p in ps], speedups, color="#2f7f6f")
    ax.axhline(1.0, color="#222222", linewidth=1.0)
    ax.set_xlabel("monomial degree p")
    ax.set_ylabel("IRP-F quad speedup vs Newton-K")
    ax.set_title("Vectorized NumPy batch benchmark")
    ax.grid(axis="y", color="#d9d4c7", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "029_fused_irp_vector_speedup.pdf", bbox_inches="tight")
    fig.savefig(OUT_FIG / "029_fused_irp_vector_speedup.png", bbox_inches="tight", dpi=220)
    plt.close(fig)


def main() -> None:
    setup()
    cases = make_cases()
    methods: list[tuple[str, Callable[[float, int, int], float]]] = [
        ("Newton-K", newton_palette_fixed),
        ("IRP-F exact", fused_irp_exact_fixed),
        ("IRP-F quad", fused_irp_quad_fixed),
        ("math.pow", math_pow_root),
    ]
    rows: list[dict[str, object]] = []
    checksum_total = 0.0
    for p in P_VALUES:
        method_rows: dict[str, dict[str, float | int]] = {}
        for name, fn in methods:
            if name == "math.pow":
                steps = 0
                err = max_error(fn, cases, p, steps)
            else:
                steps, err = calibrate(fn, cases, p)
            ns, checksum = bench_method(fn, cases, p, steps)
            checksum_total += checksum
            method_rows[name] = {
                "steps": steps,
                "max_rel_error": err,
                "ns_per_root": ns,
            }
            print(f"p={p:2d} {name:12s} steps={steps:2d} ns/root={ns:8.1f} err={err:.3e}")
        rows.append({"p": p, "methods": method_rows})

    summary = {
        "description": "Calibrated fixed-budget scalar benchmark for fused Pandrosion IRP versus palette-scaled Newton.",
        "tolerance": TOL,
        "max_iters": MAX_ITERS,
        "samples": SAMPLES,
        "min_time_seconds": MIN_TIME,
        "case_count": len(cases),
        "p_values": P_VALUES,
        "checksum": checksum_total,
        "rows": rows,
    }
    (OUT_BENCH / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (OUT_BENCH / "summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["p", "method", "steps", "max_rel_error", "ns_per_root"])
        for row in rows:
            for name, data in row["methods"].items():  # type: ignore[union-attr]
                writer.writerow([row["p"], name, data["steps"], data["max_rel_error"], data["ns_per_root"]])
    write_table(rows)
    write_plot(rows)
    vector_rows = run_vector_benchmark(rows)
    write_vector_table(vector_rows)
    write_vector_plot(vector_rows)
    summary["vector_rows"] = vector_rows
    (OUT_BENCH / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {OUT_BENCH / 'summary.json'}")
    print(f"wrote {OUT_FIG / '029_fused_irp_benchmark_summary.tex'}")
    print(f"wrote {OUT_FIG / '029_fused_irp_vector_benchmark_summary.tex'}")


if __name__ == "__main__":
    main()
