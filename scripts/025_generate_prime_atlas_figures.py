from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle


OUT = Path("latex/figures")


def primes_upto(n: int) -> list[int]:
    primes: list[int] = []
    for k in range(2, n + 1):
        ok = True
        for p in primes:
            if p * p > k:
                break
            if k % p == 0:
                ok = False
                break
        if ok:
            primes.append(k)
    return primes


def setup() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.facecolor": "#fbfbf7",
            "figure.facecolor": "#fbfbf7",
            "savefig.facecolor": "#fbfbf7",
            "axes.edgecolor": "#2f3437",
            "axes.labelcolor": "#1d2327",
            "xtick.color": "#293136",
            "ytick.color": "#293136",
            "text.color": "#182126",
            "axes.titleweight": "bold",
            "axes.titlesize": 15,
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=220)
    plt.close(fig)


def fig_prime_atlas() -> None:
    primes = primes_upto(37)
    fig, ax = plt.subplots(figsize=(10.5, 4.3))
    ax.set_xlim(-0.3, len(primes) - 0.7)
    ax.set_ylim(-1.3, 1.45)
    ax.axis("off")

    colors = ["#2f6f9f", "#2a9d8f", "#e9c46a", "#e76f51", "#8d5a97"]
    for i, p in enumerate(primes):
        h = np.log(p)
        width = 0.58
        direct = Rectangle((i - width / 2, 0.05), width, 0.82, facecolor=colors[i % len(colors)], alpha=0.82, edgecolor="#1e292e", lw=0.8)
        inv = Rectangle((i - width / 2, -0.87), width, 0.82, facecolor=colors[i % len(colors)], alpha=0.32, edgecolor="#1e292e", lw=0.8)
        ax.add_patch(direct)
        ax.add_patch(inv)
        ax.text(i, 0.47, rf"$\ell={p}$", ha="center", va="center", fontsize=8.5, color="white", weight="bold")
        ax.text(i, -0.47, rf"$1-s$", ha="center", va="center", fontsize=8, color="#162026", weight="bold")
        ax.plot([i, i], [-0.02, 0.02], color="#1e292e", lw=1)
        ax.text(i, 1.05, rf"$\log\ell={h:.2f}$", ha="center", va="bottom", fontsize=7.5)

    ax.text(-0.22, 0.47, r"direct chart  $\ell^{-s}$", ha="right", va="center", fontsize=12, weight="bold")
    ax.text(-0.22, -0.47, r"inverse chart  $\ell^{-(1-s)}$", ha="right", va="center", fontsize=12, weight="bold")
    ax.text((len(primes) - 1) / 2, 1.34, "Prime-by-prime Pandrosion atlas: each Euler factor has two logarithmic charts", ha="center", va="center", fontsize=15, weight="bold")
    ax.text((len(primes) - 1) / 2, -1.16, r"Euler factor as a geometric cell:  $1+\ell^{-s}+\ell^{-2s}+\cdots$", ha="center", va="center", fontsize=12)
    save(fig, "025_prime_atlas_cells")


def fig_critical_balance() -> None:
    sig = np.linspace(-0.2, 1.2, 500)
    primes = np.array(primes_upto(97), dtype=float)
    weights = np.exp(-np.log(primes) / np.log(97))
    c = np.sum(weights * np.log(primes) ** 2)
    e_direct = c * sig**2
    e_inverse = c * (1 - sig) ** 2
    e_min = np.minimum(e_direct, e_inverse)
    imbalance = e_direct - e_inverse

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(sig, e_direct, color="#2563a5", lw=2.6, label=r"direct cost  $\sum w_\ell(\sigma\log\ell)^2$")
    ax.plot(sig, e_inverse, color="#c14953", lw=2.6, label=r"inverse cost  $\sum w_\ell((1-\sigma)\log\ell)^2$")
    ax.fill_between(sig, 0, e_min, color="#2a9d8f", alpha=0.18, label="selected chart cost")
    ax.axvline(0.5, color="#151d22", ls="--", lw=1.8)
    ax.scatter([0.5], [c * 0.25], s=90, color="#f4a261", edgecolor="#151d22", zorder=5)
    ax.text(0.515, c * 0.25 * 1.06, r"critical balance  $\sigma=\frac{1}{2}$", fontsize=12, weight="bold")
    ax.set_xlabel(r"real part  $\sigma=\Re(s)$")
    ax.set_ylabel("truncated logarithmic chart cost")
    ax.set_title("The critical line as Thales--Riemann chart equilibrium")
    ax.grid(True, color="#d8d6ca", lw=0.8, alpha=0.65)
    ax.legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", loc="upper center", fontsize=9)
    ax2 = ax.twinx()
    ax2.plot(sig, imbalance / np.max(np.abs(imbalance)), color="#5a3e85", lw=1.4, alpha=0.65)
    ax2.set_ylabel("normalized direct-minus-inverse imbalance", color="#5a3e85")
    ax2.tick_params(axis="y", colors="#5a3e85")
    save(fig, "025_critical_line_balance")


def truncated_zeta(sigma: float, t: np.ndarray, prime_limit: int = 97, terms: int = 9) -> np.ndarray:
    primes = np.array(primes_upto(prime_limit), dtype=float)
    out = np.ones_like(t, dtype=np.complex128)
    s = sigma + 1j * t
    for p in primes:
        q = np.exp(-s * np.log(p))
        geom = (1 - q**terms) / (1 - q)
        out *= geom
    return out


def fig_phase_interference() -> None:
    t = np.linspace(0, 90, 1400)
    sigmas = [0.35, 0.50, 0.65]
    colors = ["#2563a5", "#111827", "#c14953"]

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.6), sharex=True, gridspec_kw={"height_ratios": [1.05, 1.0]})
    for sigma, color in zip(sigmas, colors):
        z = truncated_zeta(sigma, t)
        amp = np.log1p(np.abs(z))
        axes[0].plot(t, amp, color=color, lw=1.45, label=rf"$\sigma={sigma:.2f}$")
    axes[0].set_ylabel(r"$\log(1+|\prod_{\ell\leq 97}E_{\ell,9}(\sigma+it)|)$")
    axes[0].set_title("Truncated Euler-Pandrosion product: phase interference is most balanced on the critical line")
    axes[0].grid(True, color="#d8d6ca", alpha=0.6)
    axes[0].legend(ncol=3, frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", loc="upper right")

    zcrit = truncated_zeta(0.5, t)
    phase = np.unwrap(np.angle(zcrit))
    phase_rate = np.gradient(phase, t)
    axes[1].plot(t, phase_rate, color="#2a9d8f", lw=1.35)
    axes[1].fill_between(t, 0, phase_rate, where=phase_rate >= 0, color="#2a9d8f", alpha=0.2)
    axes[1].fill_between(t, 0, phase_rate, where=phase_rate < 0, color="#e76f51", alpha=0.2)
    axes[1].set_xlabel(r"height  $t$")
    axes[1].set_ylabel("phase velocity on " + r"$\sigma=1/2$")
    axes[1].grid(True, color="#d8d6ca", alpha=0.6)

    save(fig, "025_euler_pandrosion_interference")


def main() -> None:
    setup()
    fig_prime_atlas()
    fig_critical_balance()
    fig_phase_interference()
    for path in sorted(OUT.glob("025_*")):
        print(path)


if __name__ == "__main__":
    main()
