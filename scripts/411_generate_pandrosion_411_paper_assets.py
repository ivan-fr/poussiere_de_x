#!/usr/bin/env python3
"""Generate 411-only figures and LaTeX fragments for the research paper."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BENCH = Path("/private/tmp/411_pandrosion_cached_general_complete_benchmark.json")
FIG = ROOT / "latex" / "figures"


FAMILY_LABELS = {
    "affine_linear_polynomial": "affine",
    "near_linear_quadratic": "quadratic",
    "near_linear_ring_quadratic": "ring quadratic",
}


def fmt(x: Any, nd: int = 3) -> str:
    try:
        y = float(x)
    except Exception:
        return "--"
    if not math.isfinite(y):
        return "--"
    return f"{y:.{nd}f}"


def sci(x: Any, nd: int = 2) -> str:
    try:
        y = float(x)
    except Exception:
        return "--"
    if not math.isfinite(y):
        return "--"
    return f"{y:.{nd}e}"


def sci_tex(x: Any, nd: int = 2) -> str:
    raw = sci(x, nd)
    if raw == "--" or "e" not in raw:
        return raw
    mant, exp = raw.split("e", 1)
    return f"{mant}\\times 10^{{{int(exp)}}}"


def load_rows() -> list[dict[str, Any]]:
    if not BENCH.exists():
        raise FileNotFoundError(f"missing benchmark JSON: {BENCH}")
    data = json.loads(BENCH.read_text(encoding="utf-8"))
    return list(data["polynomial_qpath"])


def family_rows(rows: list[dict[str, Any]], family: str) -> list[dict[str, Any]]:
    out = [r for r in rows if r["family"] == family]
    return sorted(out, key=lambda r: int(r["n"]))


def plot_runtime(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colors = {
        "affine_linear_polynomial": "#1f77b4",
        "near_linear_quadratic": "#2ca02c",
        "near_linear_ring_quadratic": "#9467bd",
    }
    for family in FAMILY_LABELS:
        fr = family_rows(rows, family)
        if not fr:
            continue
        n = [r["n"] for r in fr]
        t = [r["411_cached_ms"] for r in fr]
        ax.plot(n, t, marker="o", linewidth=2.0, color=colors[family], label=f"411 {FAMILY_LABELS[family]}")
    base = family_rows(rows, "near_linear_quadratic")
    ax.plot([r["n"] for r in base], [r["matmul_float16_ms"] for r in base], "--", marker="s", color="#d62728", linewidth=2.0, label="torch.matmul fp16")
    ax.plot([r["n"] for r in base], [r["matmul_float32_ms"] for r in base], ":", marker="^", color="#111111", linewidth=2.0, label="torch.matmul fp32")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("wall time (ms, lower is better)")
    ax.set_title("411 q-path runtime against dense AI matmul")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "411_runtime_vs_matmul.pdf")
    fig.savefig(FIG / "411_runtime_vs_matmul.png", dpi=200)
    plt.close(fig)


def plot_speedup(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for family in FAMILY_LABELS:
        fr = family_rows(rows, family)
        if not fr:
            continue
        ax.plot(
            [r["n"] for r in fr],
            [r.get("speedup_411_cached_vs_matmul_fp16", 0.0) for r in fr],
            marker="o",
            linewidth=2.0,
            label=FAMILY_LABELS[family],
        )
    ax.axhline(1.0, color="#222222", linewidth=1.2, linestyle="--")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("dimension n")
    ax.set_ylabel("speedup over torch.matmul fp16")
    ax.set_title("Where 411 crosses dense fp16 matmul")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "411_speedup_vs_fp16.pdf")
    fig.savefig(FIG / "411_speedup_vs_fp16.png", dpi=200)
    plt.close(fig)


def plot_cache_ablation(rows: list[dict[str, Any]]) -> None:
    fr = family_rows(rows, "near_linear_quadratic")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot([r["n"] for r in fr], [r["411_cached_ms"] for r in fr], marker="o", linewidth=2.0, label="411 cached q-path")
    ax.plot([r["n"] for r in fr], [r["411_no_cache_no_reuse_ms"] for r in fr], marker="s", linewidth=2.0, label="411 no cache/no reuse")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("wall time (ms)")
    ax.set_title("411 internal cache removes the fixed setup cost")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "411_cache_ablation.pdf")
    fig.savefig(FIG / "411_cache_ablation.png", dpi=200)
    plt.close(fig)


def plot_accuracy(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for family in FAMILY_LABELS:
        fr = family_rows(rows, family)
        if not fr:
            continue
        ax.plot(
            [r["n"] for r in fr],
            [r["411_cached"]["residual"] for r in fr],
            marker="o",
            linewidth=2.0,
            label=FAMILY_LABELS[family],
        )
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("faithful residual")
    ax.set_title("411 residuals on q-active polynomial systems")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "411_residuals.pdf")
    fig.savefig(FIG / "411_residuals.png", dpi=200)
    plt.close(fig)


def write_summary_tex(rows: list[dict[str, Any]]) -> None:
    q_active = sum(1 for r in rows if r["411_cached"]["q_path_active"])
    wins = sum(1 for r in rows if r.get("beats_matmul_fp16"))
    best_speed = max(float(r.get("speedup_411_cached_vs_matmul_fp16", 0.0)) for r in rows)
    best_fp32_speed = max(float(r.get("speedup_411_cached_vs_matmul_fp32", 0.0)) for r in rows)
    best = max(rows, key=lambda r: float(r.get("speedup_411_cached_vs_matmul_fp16", 0.0)))
    max_res = max(float(r["411_cached"]["residual"]) for r in rows)
    max_root_err = max(float(r["411_cached"]["root_relative_error"]) for r in rows)
    selected = [r for r in rows if int(r["n"]) >= 1024]

    lines = [
        "% Auto-generated by scripts/411_generate_pandrosion_411_paper_assets.py",
        f"\\newcommand{{\\RowsCount}}{{{len(rows)}}}",
        f"\\newcommand{{\\QActiveRows}}{{{q_active}}}",
        f"\\newcommand{{\\FpSixteenWins}}{{{wins}}}",
        f"\\newcommand{{\\BestFpSixteenSpeedup}}{{{fmt(best_speed, 2)}}}",
        f"\\newcommand{{\\BestFpThirtyTwoSpeedup}}{{{fmt(best_fp32_speed, 2)}}}",
        f"\\newcommand{{\\BestSpeedFamily}}{{{FAMILY_LABELS[best['family']]}}}",
        f"\\newcommand{{\\BestSpeedN}}{{{best['n']}}}",
        f"\\newcommand{{\\BestSpeedQK}}{{{best['q']}/{best['k']}}}",
        f"\\newcommand{{\\BestSpeedTime}}{{{fmt(best['411_cached_ms'], 3)}}}",
        f"\\newcommand{{\\BestSpeedMatmul}}{{{fmt(best['matmul_float16_ms'], 3)}}}",
        f"\\newcommand{{\\MaxResidual}}{{{sci_tex(max_res)}}}",
        f"\\newcommand{{\\MaxRootError}}{{{sci_tex(max_root_err)}}}",
        "",
        "\\newcommand{\\BenchmarkTableRows}{%",
    ]
    for r in selected:
        c = r["411_cached"]
        lines.append(
            f"{FAMILY_LABELS[r['family']]} & {r['n']} & {r['q']}/{r['k']} & "
            f"{fmt(r['411_cached_ms'])} & {fmt(r['matmul_float16_ms'])} & "
            f"{fmt(r.get('speedup_411_cached_vs_matmul_fp16', 0.0), 2)} & "
            f"{fmt(r['matmul_float32_ms'])} & {fmt(r.get('speedup_411_cached_vs_matmul_fp32', 0.0), 2)} & "
            f"{sci(c['residual'])} & {sci(c['root_relative_error'])} \\\\"
        )
    lines.append("}")
    lines.append("")
    lines.append("\\newcommand{\\FullBenchmarkRows}{%")
    for r in rows:
        c = r["411_cached"]
        active = "yes" if c["q_path_active"] else "no"
        lines.append(
            f"{FAMILY_LABELS[r['family']]} & {r['n']} & {r['q']}/{r['k']} & "
            f"{fmt(r['411_cached_ms'])} & {fmt(r['matmul_float16_ms'])} & "
            f"{fmt(r.get('speedup_411_cached_vs_matmul_fp16', 0.0), 2)} & "
            f"{sci(c['residual'])} & {active} \\\\"
        )
    lines.append("}")
    (FIG / "411_benchmark_summary.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    plot_runtime(rows)
    plot_speedup(rows)
    plot_cache_ablation(rows)
    plot_accuracy(rows)
    write_summary_tex(rows)
    print(f"generated 411 paper assets in {FIG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
