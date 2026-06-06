from __future__ import annotations

import math
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
            "axes.edgecolor": "#263238",
            "axes.labelcolor": "#1f2933",
            "xtick.color": "#263238",
            "ytick.color": "#263238",
            "text.color": "#172026",
            "axes.titleweight": "bold",
            "axes.titlesize": 13,
        }
    )


def exp0(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    scale = np.where(n > 0.0, np.tanh(0.5 * n) / n, 0.5)
    return scale * v


def log0(x: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(x, axis=-1, keepdims=True)
    scale = np.where(n > 0.0, 2.0 * np.arctanh(n) / n, 2.0)
    return scale * x


def theta_exact(p: int, d: np.ndarray) -> np.ndarray:
    return d + (p - 1) * np.log1p(np.expm1(-d) / p)


def irp_layer_log(U: np.ndarray, X: np.ndarray, p: int) -> np.ndarray:
    A = X - p * U
    d = np.linalg.norm(A, axis=-1, keepdims=True)
    direction = np.divide(A, d, out=np.zeros_like(A), where=d > 0.0)
    th = theta_exact(p, d)
    return U + th * direction


def run_layers(U: np.ndarray, X: np.ndarray, p: int, layers: int) -> np.ndarray:
    out = U.copy()
    for _ in range(layers):
        out = irp_layer_log(out, X, p)
    return out


def local_order(xs: np.ndarray, ys: np.ndarray) -> float:
    mask = ys > 1.0e-300
    coeff = np.polyfit(np.log(xs[mask]), np.log(ys[mask]), 1)
    return float(coeff[0])


def main() -> None:
    setup()
    p = 3
    X = np.array([1.0, -0.55])
    root_log = X / p
    direction = np.array([0.45, 0.89])
    direction = direction / np.linalg.norm(direction)
    eps = np.logspace(-2, -0.5, 25)
    U0 = root_log[None, :] + eps[:, None] * direction[None, :]

    # Exercise the actual Poincare Exp/Log roundtrip before applying the log theorem.
    points = exp0(U0)
    U0_roundtrip = log0(points)

    errors: dict[int, np.ndarray] = {}
    for layers in [1, 2, 3]:
        Uk = run_layers(U0_roundtrip, X, p, layers)
        errors[layers] = np.linalg.norm(Uk - root_log[None, :], axis=1)

    orders = {layers: local_order(eps, err) for layers, err in errors.items()}

    with (OUT / "031_hadamard_irp_order_table.tex").open("w", encoding="utf-8") as f:
        f.write("\\begin{tabular}{rrrr}\n")
        f.write("\\toprule\n")
        f.write("layers & predicted order & fitted order & error at $10^{-2}$ \\\\\n")
        f.write("\\midrule\n")
        for layers in [1, 2, 3]:
            predicted = 2**layers
            fitted = orders[layers]
            final_err = errors[layers][0]
            f.write(f"{layers} & {predicted} & {fitted:.3f} & {final_err:.2e} \\\\\n")
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    colors = {1: "#2f6f9f", 2: "#8f5a83", 3: "#2f7f6f"}
    for layers in [1, 2, 3]:
        ax.loglog(eps, errors[layers], marker="o", markersize=3.5, linewidth=1.6, color=colors[layers], label=f"{layers} layer(s)")
    ax.set_xlabel("initial tangent error norm")
    ax.set_ylabel("post-IRP tangent error norm")
    ax.set_title("Poincare-ball realization of Hadamard IRP local order")
    ax.grid(True, which="both", color="#d9d4c7", linewidth=0.7)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "031_hadamard_irp_order.pdf", bbox_inches="tight")
    fig.savefig(OUT / "031_hadamard_irp_order.png", bbox_inches="tight", dpi=220)
    plt.close(fig)
    print("wrote latex/figures/031_hadamard_irp_order_table.tex")
    print("wrote latex/figures/031_hadamard_irp_order.pdf")


if __name__ == "__main__":
    main()
