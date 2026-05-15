from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.special import digamma, loggamma


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


def finite_cell_log_and_derivative(s: np.ndarray, prime_limit: int = 97, depth: int = 9) -> tuple[np.ndarray, np.ndarray]:
    primes = np.array(primes_upto(prime_limit), dtype=float)
    log_prod = np.zeros_like(s, dtype=np.complex128)
    deriv = np.zeros_like(s, dtype=np.complex128)
    for ell in primes:
        L = np.log(ell)
        q = np.exp(-s * L)
        powers = np.vstack([q**r for r in range(depth)])
        E = powers.sum(axis=0)
        dE = -L * np.vstack([r * q**r for r in range(depth)]).sum(axis=0)
        log_prod += np.log(E)
        deriv += dE / E
    return log_prod, deriv


def arch_log_and_derivative(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    log_c = np.log(0.5) + np.log(s) + np.log(s - 1) - 0.5 * s * np.log(np.pi) + loggamma(s / 2)
    deriv = 1 / s + 1 / (s - 1) - 0.5 * np.log(np.pi) + 0.5 * digamma(s / 2)
    return log_c, deriv


def fig_completed_atlas_architecture() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    def box(x, y, w, h, text, color, alpha=0.9, fs=11):
        rect = Rectangle((x, y), w, h, facecolor=color, edgecolor="#162026", lw=1.2, alpha=alpha)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, weight="bold")
        return rect

    box(0.4, 4.35, 2.1, 0.9, "archimedean\nGamma chart", "#e9c46a")
    box(0.4, 2.55, 2.1, 0.9, "prime atlas\nEuler cells", "#2a9d8f", 0.82)
    box(0.4, 0.75, 2.1, 0.9, "reciprocal\n1-s charts", "#8d5a97", 0.72)
    box(3.5, 3.35, 2.5, 1.05, "completed\nphase sum", "#f4a261")
    box(6.9, 3.35, 2.45, 1.05, "Hardy/Pandrosion\nreal slice", "#5b8db8", 0.88)
    box(6.9, 1.45, 2.45, 1.05, "zero candidates:\nphase + amplitude", "#c14953", 0.82)

    arrows = [
        ((2.55, 4.8), (3.45, 4.0)),
        ((2.55, 3.0), (3.45, 3.85)),
        ((2.55, 1.2), (3.45, 3.55)),
        ((6.05, 3.88), (6.85, 3.88)),
        ((8.1, 3.3), (8.1, 2.55)),
    ]
    for a, b in arrows:
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="->", mutation_scale=16, lw=1.6, color="#162026"))

    ax.text(5.15, 5.35, r"$\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$", ha="center", fontsize=15)
    ax.text(5.12, 2.78, r"$\Theta_{P,N}(t)=\arg C(1/2+it)+\arg\prod_{\ell\in P}E_{\ell,N}(1/2+it)$", ha="center", fontsize=11)
    ax.text(5.08, 0.45, "026 adds the archimedean chart and studies phase velocity, not only critical-line balance.", ha="center", fontsize=11)
    ax.set_title("Completed Pandrosion prime atlas: archimedean chart + prime geometric cells")
    save(fig, "026_completed_atlas_architecture")


def fig_phase_velocity() -> None:
    t = np.linspace(1.0, 90.0, 1200)
    s = 0.5 + 1j * t
    _, dprime = finite_cell_log_and_derivative(s)
    _, darch = arch_log_and_derivative(s)
    prime_vel = np.real(dprime)
    arch_vel = np.real(darch)
    total_vel = prime_vel + arch_vel
    theta_prime = 0.5 * np.real(digamma(0.25 + 0.5j * t)) - 0.5 * np.log(np.pi)

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.7), sharex=True)
    axes[0].plot(t, theta_prime, color="#111827", lw=2.0, label=r"Riemann--Siegel $\vartheta'(t)$")
    axes[0].plot(t, arch_vel, color="#e76f51", lw=1.25, ls="--", label="archimedean phase velocity")
    axes[0].set_ylabel("archimedean velocity")
    axes[0].grid(True, color="#d8d6ca", alpha=0.65)
    axes[0].legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", loc="upper left")

    axes[1].plot(t, prime_vel, color="#2a9d8f", lw=1.0, alpha=0.8, label="prime-cell velocity")
    axes[1].plot(t, total_vel, color="#5a3e85", lw=1.35, label="completed finite-atlas velocity")
    axes[1].axhline(0, color="#162026", lw=0.8)
    axes[1].set_xlabel(r"height $t$ on $s=1/2+it$")
    axes[1].set_ylabel(r"$d\arg/dt$")
    axes[1].grid(True, color="#d8d6ca", alpha=0.65)
    axes[1].legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", loc="upper right")
    fig.suptitle("Phase velocity separates into archimedean and prime Pandrosion contributions", y=0.98, fontsize=16, weight="bold")
    save(fig, "026_completed_phase_velocity")


def fig_phase_amplitude_candidates() -> None:
    t = np.linspace(1.0, 90.0, 1600)
    s = 0.5 + 1j * t
    lp, _ = finite_cell_log_and_derivative(s)
    lc, _ = arch_log_and_derivative(s)
    log_xi = lp + lc
    amp = np.log1p(np.abs(np.exp(log_xi - np.max(np.real(log_xi)))))
    phase = np.unwrap(np.imag(log_xi))
    phase_mod = np.mod(phase + np.pi / 2, np.pi) - np.pi / 2

    # Simple local minima of the scaled amplitude.
    idx = np.where((amp[1:-1] < amp[:-2]) & (amp[1:-1] < amp[2:]))[0] + 1
    idx = idx[amp[idx] < np.quantile(amp, 0.20)]
    idx = idx[:: max(1, len(idx) // 18)]

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.7), sharex=True)
    axes[0].plot(t, amp, color="#2563a5", lw=1.2)
    if len(idx):
        axes[0].scatter(t[idx], amp[idx], s=24, color="#f4a261", edgecolor="#162026", zorder=5, label="finite-atlas low-amplitude candidates")
    axes[0].set_ylabel("scaled completed amplitude")
    axes[0].grid(True, color="#d8d6ca", alpha=0.65)
    axes[0].legend(frameon=True, facecolor="#fbfbf7", edgecolor="#b9b5a7", loc="upper right")

    axes[1].plot(t, phase_mod, color="#c14953", lw=0.95)
    axes[1].axhline(0, color="#162026", lw=0.8)
    axes[1].set_xlabel(r"height $t$")
    axes[1].set_ylabel("phase mod pi")
    axes[1].grid(True, color="#d8d6ca", alpha=0.65)
    fig.suptitle("Finite completed atlas: zero candidates require amplitude collapse and phase alignment", y=0.98, fontsize=16, weight="bold")
    save(fig, "026_phase_amplitude_candidates")


def main() -> None:
    setup()
    fig_completed_atlas_architecture()
    fig_phase_velocity()
    fig_phase_amplitude_candidates()
    for path in sorted(OUT.glob("026_*")):
        print(path)


if __name__ == "__main__":
    main()
