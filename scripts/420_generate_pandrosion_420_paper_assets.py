#!/usr/bin/env python3
"""Generate 420 figures and LaTeX fragments for the research paper."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FULL = Path("/private/tmp/420_full_irp_no_fallback_vs_torch_matmul_benchmark.json")
PARTIAL = Path("/private/tmp/420_partial_cascade_vs_torch_matmul_benchmark.json")
FIG = ROOT / "latex" / "figures"


LABELS = {
    "q_active": "q active",
    "partial_stage2": "partial stage 2",
    "partial_stage3": "partial stage 3",
    "random": "random full n",
}
COLORS = {
    "q_active": "#1f77b4",
    "partial_stage2": "#2ca02c",
    "partial_stage3": "#9467bd",
    "random": "#d62728",
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing benchmark JSON: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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


def rows_for(rows: list[dict[str, Any]], family: str) -> list[dict[str, Any]]:
    return sorted([r for r in rows if r["family"] == family], key=lambda r: int(r["n"]))


def val(row: dict[str, Any], key: str, field: str = "mean") -> float:
    return float(row[key][field])


def plot_runtime(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for fam in ("q_active", "partial_stage2", "partial_stage3", "random"):
        fr = rows_for(rows, fam)
        ax.plot(
            [r["n"] for r in fr],
            [val(r, "ms_420") for r in fr],
            marker="o",
            linewidth=2.0,
            color=COLORS[fam],
            label=f"420 {LABELS[fam]}",
        )
    base = rows_for(rows, "q_active")
    ax.plot([r["n"] for r in base], [val(r, "square_matmul_float16_ms") for r in base], "--", color="#111111", marker="s", label="torch.matmul fp16 GEMM")
    ax.plot([r["n"] for r in base], [val(r, "matrix_vector_float16_ms") for r in base], ":", color="#555555", marker="^", label="torch.matmul fp16 matvec")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("wall time (ms, lower is better)")
    ax.set_title("420 runtime: q, partial cascade, and full-n IRP")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "420_runtime_cascade_vs_matmul.pdf")
    fig.savefig(FIG / "420_runtime_cascade_vs_matmul.png", dpi=220)
    plt.close(fig)


def plot_speedup(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for fam in ("q_active", "partial_stage2", "partial_stage3", "random"):
        fr = rows_for(rows, fam)
        ax.plot(
            [r["n"] for r in fr],
            [val(r, "speedup_vs_square_matmul_float16") for r in fr],
            marker="o",
            linewidth=2.0,
            color=COLORS[fam],
            label=LABELS[fam],
        )
    ax.axhline(1.0, color="#222222", linestyle="--", linewidth=1.2)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("speedup over torch.matmul fp16 GEMM")
    ax.set_title("420 speedup appears when the cascade stops before n")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "420_speedup_cascade_vs_fp16_gemm.pdf")
    fig.savefig(FIG / "420_speedup_cascade_vs_fp16_gemm.png", dpi=220)
    plt.close(fig)


def plot_accepted_dim(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for fam in ("q_active", "partial_stage2", "partial_stage3", "random"):
        fr = rows_for(rows, fam)
        ax.plot(
            [r["n"] for r in fr],
            [val(r, "accepted_dim") for r in fr],
            marker="o",
            linewidth=2.0,
            color=COLORS[fam],
            label=LABELS[fam],
        )
    base = rows_for(rows, "q_active")
    ax.plot([r["n"] for r in base], [r["n"] for r in base], "--", color="#222222", linewidth=1.1, label="full n")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xlabel("dimension n")
    ax.set_ylabel("accepted IRP dimension")
    ax.set_title("Cascade depth selected by 420")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "420_accepted_dimension.pdf")
    fig.savefig(FIG / "420_accepted_dimension.png", dpi=220)
    plt.close(fig)


def plot_error(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    for fam in ("q_active", "partial_stage2", "partial_stage3", "random"):
        fr = rows_for(rows, fam)
        err = [max(val(r, "relative_error", "max"), 1e-18) for r in fr]
        ax.plot([r["n"] for r in fr], err, marker="o", linewidth=2.0, color=COLORS[fam], label=LABELS[fam])
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("max relative error")
    ax.set_title("420 accuracy across cascade regimes")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "420_relative_error.pdf")
    fig.savefig(FIG / "420_relative_error.png", dpi=220)
    plt.close(fig)


def write_summary(rows: list[dict[str, Any]], full: dict[str, Any], partial: dict[str, Any]) -> None:
    q_rows = rows_for(rows, "q_active")
    p2_rows = rows_for(rows, "partial_stage2")
    p3_rows = rows_for(rows, "partial_stage3")
    random_rows = rows_for(rows, "random")
    best_q = max(q_rows, key=lambda r: val(r, "speedup_vs_square_matmul_float16"))
    best_p2 = max(p2_rows, key=lambda r: val(r, "speedup_vs_square_matmul_float16"))
    best_p3 = max(p3_rows, key=lambda r: val(r, "speedup_vs_square_matmul_float16"))
    max_err = max(val(r, "relative_error", "max") for r in rows)
    max_p2_err = max(val(r, "relative_error", "max") for r in p2_rows)
    max_p3_err = max(val(r, "relative_error", "max") for r in p3_rows)
    selected = [r for r in rows if int(r["n"]) in {512, 1024, 2048, 4096}]

    lines = [
        "% Auto-generated by scripts/420_generate_pandrosion_420_paper_assets.py",
        f"\\newcommand{{\\RowsCount}}{{{len(rows)}}}",
        f"\\newcommand{{\\BenchDevice}}{{{partial.get('device', 'unknown')}}}",
        f"\\newcommand{{\\BenchTorch}}{{{partial.get('torch', 'unknown')}}}",
        f"\\newcommand{{\\BenchNumpy}}{{{partial.get('numpy', 'unknown')}}}",
        f"\\newcommand{{\\QBestSpeedup}}{{{fmt(val(best_q, 'speedup_vs_square_matmul_float16'), 2)}}}",
        f"\\newcommand{{\\QBestN}}{{{best_q['n']}}}",
        f"\\newcommand{{\\QBestDim}}{{{fmt(val(best_q, 'accepted_dim'), 0)}}}",
        f"\\newcommand{{\\PtwoBestSpeedup}}{{{fmt(val(best_p2, 'speedup_vs_square_matmul_float16'), 2)}}}",
        f"\\newcommand{{\\PtwoBestN}}{{{best_p2['n']}}}",
        f"\\newcommand{{\\PtwoBestDim}}{{{fmt(val(best_p2, 'accepted_dim'), 0)}}}",
        f"\\newcommand{{\\PthreeBestSpeedup}}{{{fmt(val(best_p3, 'speedup_vs_square_matmul_float16'), 2)}}}",
        f"\\newcommand{{\\PthreeBestN}}{{{best_p3['n']}}}",
        f"\\newcommand{{\\PthreeBestDim}}{{{fmt(val(best_p3, 'accepted_dim'), 0)}}}",
        f"\\newcommand{{\\MaxRelError}}{{{sci_tex(max_err)}}}",
        f"\\newcommand{{\\PtwoMaxRelError}}{{{sci_tex(max_p2_err)}}}",
        f"\\newcommand{{\\PthreeMaxRelError}}{{{sci_tex(max_p3_err)}}}",
        "",
        "\\newcommand{\\CascadeTableRows}{%",
    ]
    for r in selected:
        lines.append(
            f"{LABELS[r['family']]} & {r['n']} & {r['q']}/{r['k']} & "
            f"{fmt(val(r, 'accepted_dim'), 0)} & {fmt(val(r, 'ms_420'))} & "
            f"{fmt(val(r, 'square_matmul_float16_ms'))} & "
            f"{fmt(val(r, 'speedup_vs_square_matmul_float16'), 2)} & "
            f"{fmt(val(r, 'matrix_vector_float16_ms'))} & "
            f"{fmt(val(r, 'speedup_vs_matrix_vector_float16'), 2)} & "
            f"{sci(val(r, 'relative_error', 'max'))} \\\\"
        )
    lines.append("}")
    lines.append("")
    lines.append("\\newcommand{\\MiddleCascadeRows}{%")
    for r in p2_rows + p3_rows:
        lines.append(
            f"{LABELS[r['family']]} & {r['n']} & {r['q']}/{r['k']} & "
            f"{fmt(val(r, 'expected_dim'), 0)} & {fmt(val(r, 'accepted_dim'), 0)} & "
            f"{fmt(val(r, 'ms_420'))} & {fmt(val(r, 'speedup_vs_square_matmul_float16'), 2)} & "
            f"{fmt(val(r, 'speedup_vs_matrix_vector_float16'), 2)} & {sci(val(r, 'relative_error', 'max'))} \\\\"
        )
    lines.append("}")
    lines.append("")
    lines.append("\\newcommand{\\RandomControlRows}{%")
    for r in random_rows:
        lines.append(
            f"random & {r['n']} & {r['q']}/{r['k']} & {fmt(val(r, 'accepted_dim'), 0)} & "
            f"{fmt(val(r, 'ms_420'))} & {fmt(val(r, 'speedup_vs_square_matmul_float16'), 2)} & "
            f"{fmt(val(r, 'speedup_vs_matrix_vector_float16'), 2)} & {sci(val(r, 'relative_error', 'max'))} \\\\"
        )
    lines.append("}")
    (FIG / "420_benchmark_summary.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    full = load_json(FULL)
    partial = load_json(PARTIAL)
    rows = list(partial["rows"])
    plot_runtime(rows)
    plot_speedup(rows)
    plot_accepted_dim(rows)
    plot_error(rows)
    write_summary(rows, full, partial)
    print(f"wrote 420 figures and summary to {FIG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
