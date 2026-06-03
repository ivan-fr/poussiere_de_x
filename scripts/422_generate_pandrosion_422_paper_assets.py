#!/usr/bin/env python3
"""Generate 422 figures and LaTeX fragments for the research paper."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BENCH = Path("/private/tmp/422_complete_vs_torch_matmul_benchmark.json")
ABLT = Path("/private/tmp/422_adaptive_dense_jump_vs_421_benchmark.json")
FIG = ROOT / "latex" / "figures"


LABELS = {
    "q_active": "q active",
    "partial_stage2": "partial stage 2",
    "partial_stage3": "partial stage 3",
    "random_full_n": "random full n",
}
COLORS = {
    "q_active": "#1f77b4",
    "partial_stage2": "#2ca02c",
    "partial_stage3": "#9467bd",
    "random_full_n": "#d62728",
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing benchmark JSON: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rows_for(rows: list[dict[str, Any]], family: str) -> list[dict[str, Any]]:
    return sorted([r for r in rows if r["family"] == family], key=lambda r: int(r["n"]))


def val(row: dict[str, Any], key: str, field: str = "mean") -> float:
    return float(row[key][field])


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


def plot_runtime(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for fam in LABELS:
        fr = rows_for(rows, fam)
        ax.plot([r["n"] for r in fr], [val(r, "ms_422") for r in fr], marker="o", linewidth=2.0, color=COLORS[fam], label=f"422 {LABELS[fam]}")
    base = rows_for(rows, "q_active")
    ax.plot([r["n"] for r in base], [val(r, "square_matmul_float16_ms") for r in base], "--", marker="s", color="#111111", linewidth=2.0, label="torch.matmul fp16 GEMM")
    ax.plot([r["n"] for r in base], [val(r, "matrix_vector_float16_ms") for r in base], ":", marker="^", color="#555555", linewidth=2.0, label="torch.matmul fp16 matvec")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("wall time (ms, lower is better)")
    ax.set_title("422 runtime against torch.matmul up to n=8192")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "422_runtime_vs_matmul.pdf")
    fig.savefig(FIG / "422_runtime_vs_matmul.png", dpi=220)
    plt.close(fig)


def plot_speedup(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for fam in LABELS:
        fr = rows_for(rows, fam)
        ax.plot([r["n"] for r in fr], [val(r, "speedup_vs_square_matmul_float16") for r in fr], marker="o", linewidth=2.0, color=COLORS[fam], label=LABELS[fam])
    ax.axhline(1.0, color="#222222", linestyle="--", linewidth=1.2)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("speedup over torch.matmul fp16 GEMM")
    ax.set_title("422 speedup over dense fp16 GEMM")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "422_speedup_vs_fp16_gemm.pdf")
    fig.savefig(FIG / "422_speedup_vs_fp16_gemm.png", dpi=220)
    plt.close(fig)


def plot_matvec_speedup(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for fam in LABELS:
        fr = rows_for(rows, fam)
        ax.plot([r["n"] for r in fr], [val(r, "speedup_vs_matrix_vector_float16") for r in fr], marker="o", linewidth=2.0, color=COLORS[fam], label=LABELS[fam])
    ax.axhline(1.0, color="#222222", linestyle="--", linewidth=1.2)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("dimension n")
    ax.set_ylabel("speedup over torch.matmul fp16 matvec")
    ax.set_title("422 against the stricter vector-shaped baseline")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "422_speedup_vs_fp16_matvec.pdf")
    fig.savefig(FIG / "422_speedup_vs_fp16_matvec.png", dpi=220)
    plt.close(fig)


def plot_accepted_dim(rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for fam in LABELS:
        fr = rows_for(rows, fam)
        ax.plot([r["n"] for r in fr], [val(r, "accepted_dim") for r in fr], marker="o", linewidth=2.0, color=COLORS[fam], label=LABELS[fam])
    base = rows_for(rows, "q_active")
    ax.plot([r["n"] for r in base], [r["n"] for r in base], "--", color="#222222", linewidth=1.1, label="full n")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xlabel("dimension n")
    ax.set_ylabel("accepted IRP dimension")
    ax.set_title("422 accepted cascade dimension")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIG / "422_accepted_dimension.pdf")
    fig.savefig(FIG / "422_accepted_dimension.png", dpi=220)
    plt.close(fig)


def plot_jump_vs_421(rows: list[dict[str, Any]]) -> None:
    random_rows = rows_for(rows, "random")
    if not random_rows:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot([r["n"] for r in random_rows], [val(r, "speedup_422_vs_421") for r in random_rows], marker="o", linewidth=2.2, color="#d62728", label="random full n")
    ax.axhline(1.0, color="#222222", linestyle="--", linewidth=1.2)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("dimension n")
    ax.set_ylabel("speedup of 422 over 421")
    ax.set_title("Adaptive dense jump improvement")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "422_adaptive_jump_vs_421.pdf")
    fig.savefig(FIG / "422_adaptive_jump_vs_421.png", dpi=220)
    plt.close(fig)


def write_summary(bench: dict[str, Any], ablt: dict[str, Any]) -> None:
    rows = list(bench["rows"])
    q = rows_for(rows, "q_active")
    p2 = rows_for(rows, "partial_stage2")
    p3 = rows_for(rows, "partial_stage3")
    rnd = rows_for(rows, "random_full_n")
    best_gemm = max(rows, key=lambda r: val(r, "speedup_vs_square_matmul_float16"))
    best_matvec = max(rows, key=lambda r: val(r, "speedup_vs_matrix_vector_float16"))
    max_err = max(val(r, "relative_error", "max") for r in rows)
    random_8192 = [r for r in rnd if int(r["n"]) == 8192][0]
    q_8192 = [r for r in q if int(r["n"]) == 8192][0]
    p2_8192 = [r for r in p2 if int(r["n"]) == 8192][0]
    p3_8192 = [r for r in p3 if int(r["n"]) == 8192][0]
    ab_rows = list(ablt.get("rows", []))
    ab_random = [r for r in ab_rows if r.get("family") == "random"]
    best_jump = max(ab_random, key=lambda r: val(r, "speedup_422_vs_421")) if ab_random else None

    lines = [
        "% Auto-generated by scripts/422_generate_pandrosion_422_paper_assets.py",
        f"\\newcommand{{\\RowsCount}}{{{len(rows)}}}",
        f"\\newcommand{{\\BenchDevice}}{{{bench.get('device', 'unknown')}}}",
        f"\\newcommand{{\\BenchTorch}}{{{bench.get('torch', 'unknown')}}}",
        f"\\newcommand{{\\BenchNumpy}}{{{bench.get('numpy', 'unknown')}}}",
        f"\\newcommand{{\\MaxDimension}}{{{max(int(r['n']) for r in rows)}}}",
        f"\\newcommand{{\\BestGemmSpeedup}}{{{fmt(val(best_gemm, 'speedup_vs_square_matmul_float16'), 2)}}}",
        f"\\newcommand{{\\BestGemmFamily}}{{{LABELS[best_gemm['family']]}}}",
        f"\\newcommand{{\\BestGemmN}}{{{best_gemm['n']}}}",
        f"\\newcommand{{\\BestMatvecSpeedup}}{{{fmt(val(best_matvec, 'speedup_vs_matrix_vector_float16'), 2)}}}",
        f"\\newcommand{{\\BestMatvecFamily}}{{{LABELS[best_matvec['family']]}}}",
        f"\\newcommand{{\\BestMatvecN}}{{{best_matvec['n']}}}",
        f"\\newcommand{{\\MaxRelError}}{{{sci_tex(max_err)}}}",
        f"\\newcommand{{\\RandomFullNTime}}{{{fmt(val(random_8192, 'ms_422'))}}}",
        f"\\newcommand{{\\RandomFullNGemmSpeedup}}{{{fmt(val(random_8192, 'speedup_vs_square_matmul_float16'), 2)}}}",
        f"\\newcommand{{\\RandomFullNMatvecSpeedup}}{{{fmt(val(random_8192, 'speedup_vs_matrix_vector_float16'), 2)}}}",
        f"\\newcommand{{\\QActiveTimeEightK}}{{{fmt(val(q_8192, 'ms_422'))}}}",
        f"\\newcommand{{\\PartialTwoTimeEightK}}{{{fmt(val(p2_8192, 'ms_422'))}}}",
        f"\\newcommand{{\\PartialThreeTimeEightK}}{{{fmt(val(p3_8192, 'ms_422'))}}}",
        f"\\newcommand{{\\QActiveSpeedEightK}}{{{fmt(val(q_8192, 'speedup_vs_square_matmul_float16'), 2)}}}",
        f"\\newcommand{{\\PartialTwoSpeedEightK}}{{{fmt(val(p2_8192, 'speedup_vs_square_matmul_float16'), 2)}}}",
        f"\\newcommand{{\\PartialThreeSpeedEightK}}{{{fmt(val(p3_8192, 'speedup_vs_square_matmul_float16'), 2)}}}",
    ]
    if best_jump is not None:
        lines.extend([
            f"\\newcommand{{\\BestJumpSpeedup}}{{{fmt(val(best_jump, 'speedup_422_vs_421'), 2)}}}",
            f"\\newcommand{{\\BestJumpN}}{{{best_jump['n']}}}",
        ])
    else:
        lines.extend(["\\newcommand{\\BestJumpSpeedup}{--}", "\\newcommand{\\BestJumpN}{--}"])
    lines.extend(["", "\\newcommand{\\MainRows}{%"])
    for r in [x for x in rows if int(x["n"]) in {1024, 2048, 4096, 8192}]:
        lines.append(
            f"{LABELS[r['family']]} & {r['n']} & {r['q']}/{r['k']} & "
            f"{fmt(val(r, 'accepted_dim'), 0)} & {fmt(val(r, 'ms_422'))} & "
            f"{fmt(val(r, 'square_matmul_float16_ms'))} & {fmt(val(r, 'speedup_vs_square_matmul_float16'), 2)} & "
            f"{fmt(val(r, 'matrix_vector_float16_ms'))} & {fmt(val(r, 'speedup_vs_matrix_vector_float16'), 2)} & "
            f"{sci(val(r, 'relative_error', 'max'))} \\\\"
        )
    lines.append("}")
    lines.append("")
    lines.append("\\newcommand{\\EightKRows}{%")
    for r in [q_8192, p2_8192, p3_8192, random_8192]:
        lines.append(
            f"{LABELS[r['family']]} & {r['q']}/{r['k']} & {fmt(val(r, 'accepted_dim'), 0)} & "
            f"{fmt(val(r, 'ms_422'))} & {fmt(val(r, 'speedup_vs_square_matmul_float16'), 2)} & "
            f"{fmt(val(r, 'speedup_vs_square_matmul_float32'), 2)} & "
            f"{fmt(val(r, 'speedup_vs_matrix_vector_float16'), 2)} & "
            f"{fmt(val(r, 'speedup_vs_matrix_vector_float32'), 2)} & "
            f"{sci(val(r, 'relative_error', 'max'))} \\\\"
        )
    lines.append("}")
    lines.append("")
    lines.append("\\newcommand{\\RandomRows}{%")
    for r in rnd:
        lines.append(
            f"{r['n']} & {r['q']}/{r['k']} & {fmt(val(r, 'accepted_dim'), 0)} & "
            f"{fmt(val(r, 'ms_422'))} & {fmt(r['adaptive_dense_jump_rate'], 2)} & "
            f"{fmt(val(r, 'speedup_vs_square_matmul_float16'), 2)} & "
            f"{fmt(val(r, 'speedup_vs_matrix_vector_float16'), 2)} & "
            f"{sci(val(r, 'relative_error', 'max'))} \\\\"
        )
    lines.append("}")
    (FIG / "422_benchmark_summary.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    bench = load_json(BENCH)
    ablt = load_json(ABLT)
    rows = list(bench["rows"])
    plot_runtime(rows)
    plot_speedup(rows)
    plot_matvec_speedup(rows)
    plot_accepted_dim(rows)
    plot_jump_vs_421(list(ablt.get("rows", [])))
    write_summary(bench, ablt)
    print(f"wrote 422 figures and summary to {FIG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
